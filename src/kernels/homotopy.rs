//! # Homotopy continuation with projective coordinates and a Cauchy endgame
//!
//! A robust solver for square polynomial systems built on the total-degree
//! homotopy of [`crate::kernels::rootfind::homotopy_system`]:
//!
//! * **projective tracking**: the system is homogenised with an extra
//!   variable `x_0` and the paths are tracked on the random affine patch
//!   `a . X = 1`, so that diverging paths stay bounded and solutions at
//!   infinity (`x_0 -> 0`) are detected and reported instead of being lost;
//! * **Cauchy endgame**: at `t = 1 - s_0` the path is continued around the
//!   circles `|1 - t| = s` in the complex `t` plane; the number of loops
//!   needed to return to the start is the winding number (cycle number) `c`
//!   of the endpoint and the mean of the values along the closed `c`-fold
//!   loop is the Cauchy-integral estimate of the endpoint. The radius is
//!   halved until two successive estimates agree and the target residual is
//!   tiny (the latter rejects the stable centroid of a cluster of close
//!   roots), which gives accurate
//!   endpoints even where the Jacobian is singular (multiple roots and
//!   solutions at infinity), where Newton's method fails;
//! * **path-jumping protection**: every step needs a small first Newton
//!   correction and a contracting corrector; afterwards non-singular finite
//!   endpoints are checked for duplicates (an isolated solution is reached
//!   by exactly one path) and the offending paths are re-tracked with
//!   successively smaller maximal steps.
#![allow(
    clippy::missing_const_for_fn,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::similar_names,
    clippy::float_cmp,
    clippy::suboptimal_flops,
    clippy::type_complexity,
    clippy::cognitive_complexity
)]

use num_complex::Complex64 as Cx;

use crate::kernels::random::Rng;
use crate::kernels::rootfind::{PolyTerm, RootError};

/// Options of [`homotopy_solve`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomotopyOptions {
    /// Seed of the random `gamma` and patch.
    pub seed: u64,
    /// `1 - t` at which the endgame starts.
    pub endgame_start: f64,
    /// Relative tolerance on the agreement of successive endgame estimates.
    pub endgame_tol: f64,
    /// Maximal step (in `t`) of the main tracking.
    pub max_step: f64,
    /// Number of re-tracking rounds against path jumping.
    pub max_retries: usize,
    /// An endpoint with `|x_0| <= infinity_tol max |X|` is at infinity.
    pub infinity_tol: f64,
}

impl Default for HomotopyOptions {
    fn default() -> Self {
        Self {
            seed: 0x9e37,
            endgame_start: 0.1,
            endgame_tol: 1e-10,
            max_step: 0.05,
            max_retries: 3,
            infinity_tol: 1e-6,
        }
    }
}

/// A classified path endpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct PathEnd {
    /// Affine coordinates `(re, im)` (NaN for points at infinity).
    pub x: Vec<(f64, f64)>,
    /// Homogeneous coordinates `(x_0, x_1, ..)` on the patch `a . X = 1`.
    pub homogeneous: Vec<(f64, f64)>,
    /// The endpoint lies at infinity.
    pub at_infinity: bool,
    /// Winding (cycle) number found by the endgame (1 for a regular path).
    pub winding: usize,
    /// The endpoint is singular (winding number > 1 or ill-conditioned
    /// Jacobian), i.e. a multiple root or a root at infinity.
    pub singular: bool,
    /// Scaled max-norm residual of the homogenised system at the endpoint.
    pub residual: f64,
    /// The tracking and the endgame completed.
    pub converged: bool,
    /// Number of re-tracking rounds the path needed.
    pub retries: usize,
    /// The endpoint still coincides with another regular endpoint after all
    /// retries (possible path jumping).
    pub suspect: bool,
}

struct Problem {
    n: usize,
    homog: Vec<Vec<PolyTerm>>,
    deg: Vec<u32>,
    gamma: Cx,
    patch: Vec<Cx>,
}

fn cpow(z: Cx, e: u32) -> Cx {
    let mut r = Cx::new(1.0, 0.0);
    for _ in 0..e {
        r *= z;
    }
    r
}

fn eval_poly(sys: &[Vec<PolyTerm>], x: &[Cx]) -> (Vec<Cx>, Vec<Vec<Cx>>) {
    let n = x.len();
    let mut f = vec![Cx::new(0.0, 0.0); sys.len()];
    let mut j = vec![vec![Cx::new(0.0, 0.0); n]; sys.len()];
    for (i, eq) in sys.iter().enumerate() {
        for term in eq {
            let pw: Vec<Cx> = (0..n).map(|k| cpow(x[k], term.exps[k])).collect();
            f[i] += pw.iter().fold(Cx::new(term.coef, 0.0), |a, &b| a * b);
            for k in 0..n {
                if term.exps[k] > 0 {
                    let mut d = Cx::new(term.coef * f64::from(term.exps[k]), 0.0);
                    for m in 0..n {
                        d *= if m == k { cpow(x[m], term.exps[m] - 1) } else { pw[m] };
                    }
                    j[i][k] += d;
                }
            }
        }
    }
    (f, j)
}

fn csolve(mut a: Vec<Vec<Cx>>, mut b: Vec<Cx>) -> Option<Vec<Cx>> {
    let n = b.len();
    for k in 0..n {
        let mut p = k;
        for i in k + 1..n {
            if a[i][k].norm() > a[p][k].norm() {
                p = i;
            }
        }
        if a[p][k].norm() < 1e-300 {
            return None;
        }
        a.swap(k, p);
        b.swap(k, p);
        for i in k + 1..n {
            let m = a[i][k] / a[k][k];
            if m.norm() != 0.0 {
                for c in k..n {
                    let v = a[k][c];
                    a[i][c] -= m * v;
                }
                let bk = b[k];
                b[i] -= m * bk;
            }
        }
    }
    let mut x = vec![Cx::new(0.0, 0.0); n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    if x.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) { None } else { Some(x) }
}

fn max_norm(v: &[Cx]) -> f64 {
    v.iter().map(|c| c.norm()).fold(0.0, f64::max)
}

impl Problem {
    /// `(H, H_x, H_t)` at `(X, t)` for complex `t`.
    fn eval(&self, x: &[Cx], t: Cx) -> (Vec<Cx>, Vec<Vec<Cx>>, Vec<Cx>) {
        let n = self.n;
        let (fv, fj) = eval_poly(&self.homog, x);
        let mut h = vec![Cx::new(0.0, 0.0); n + 1];
        let mut hx = vec![vec![Cx::new(0.0, 0.0); n + 1]; n + 1];
        let mut ht = vec![Cx::new(0.0, 0.0); n + 1];
        for i in 0..n {
            let d = self.deg[i];
            let g = cpow(x[i + 1], d) - cpow(x[0], d);
            let dgi = cpow(x[i + 1], d - 1) * f64::from(d);
            let dg0 = -cpow(x[0], d - 1) * f64::from(d);
            h[i] = self.gamma * g * (Cx::new(1.0, 0.0) - t) + fv[i] * t;
            ht[i] = fv[i] - self.gamma * g;
            for k in 0..=n {
                hx[i][k] = fj[i][k] * t;
            }
            hx[i][i + 1] += self.gamma * dgi * (Cx::new(1.0, 0.0) - t);
            hx[i][0] += self.gamma * dg0 * (Cx::new(1.0, 0.0) - t);
        }
        let mut p = Cx::new(-1.0, 0.0);
        for k in 0..=n {
            p += self.patch[k] * x[k];
            hx[n][k] = self.patch[k];
        }
        h[n] = p;
        (h, hx, ht)
    }

    fn tangent(&self, x: &[Cx], t: Cx) -> Option<Vec<Cx>> {
        let (_, hx, ht) = self.eval(x, t);
        csolve(hx, ht.iter().map(|v| -*v).collect())
    }

    /// Tracks along the straight segment `ta -> tb` in the `t` plane.
    fn track(&self, x0: &[Cx], ta: Cx, tb: Cx, max_step: f64) -> Option<Vec<Cx>> {
        let n1 = self.n + 1;
        let dt = tb - ta;
        let mut x = x0.to_vec();
        let mut tau = 0.0_f64;
        let mut step = (0.02_f64).min(max_step);
        let mut guard = 0usize;
        while tau < 1.0 {
            guard += 1;
            if guard > 400_000 {
                return None;
            }
            let s = step.min(1.0 - tau);
            let t0 = ta + dt * tau;
            let rk = || -> Option<Vec<Cx>> {
                let k1 = self.tangent(&x, t0)?;
                let xa: Vec<Cx> = (0..n1).map(|i| x[i] + k1[i] * dt * (0.5 * s)).collect();
                let k2 = self.tangent(&xa, ta + dt * (tau + 0.5 * s))?;
                let xb: Vec<Cx> = (0..n1).map(|i| x[i] + k2[i] * dt * (0.5 * s)).collect();
                let k3 = self.tangent(&xb, ta + dt * (tau + 0.5 * s))?;
                let xc: Vec<Cx> = (0..n1).map(|i| x[i] + k3[i] * dt * s).collect();
                let k4 = self.tangent(&xc, ta + dt * (tau + s))?;
                Some(
                    (0..n1)
                        .map(|i| x[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * dt * (s / 6.0))
                        .collect(),
                )
            };
            let tn = if s >= 1.0 - tau { tb } else { ta + dt * (tau + s) };
            let mut ok = false;
            let mut xn = x.clone();
            if let Some(pred) = rk() {
                xn = pred;
                let mut prev = f64::INFINITY;
                for it in 0..4 {
                    let (h, hx, _) = self.eval(&xn, tn);
                    let Some(dx) = csolve(hx, h) else { break };
                    for i in 0..n1 {
                        xn[i] -= dx[i];
                    }
                    let dn = max_norm(&dx);
                    if it == 0 && dn > 2e-3 * (1.0 + max_norm(&xn)) {
                        break; // predictor too inaccurate: possible jump, shrink the step
                    }
                    if dn < 1e-13 * (1.0 + max_norm(&xn)) {
                        ok = true;
                        break;
                    }
                    if it > 0 && dn > 0.5 * prev {
                        // stagnation at the rounding-noise floor of an ill-conditioned system
                        ok = dn < 1e-9 * (1.0 + max_norm(&xn)) && dn < 2.0 * prev;
                        break;
                    }
                    prev = dn;
                }
            }
            if ok && xn.iter().all(|c| c.re.is_finite() && c.im.is_finite()) {
                x = xn;
                tau += s;
                if s >= 1.0 - (tau - s) {
                    tau = 1.0;
                }
                step = (step * 1.7).min(max_step);
            } else {
                step *= 0.5;
                if step < 1e-12 {
                    return None;
                }
            }
        }
        Some(x)
    }

    fn cond_estimate(&self, x: &[Cx]) -> f64 {
        let (_, hx, _) = self.eval(x, Cx::new(1.0, 0.0));
        let n1 = self.n + 1;
        let nrm = |m: &[Vec<Cx>]| m.iter().map(|r| r.iter().map(|c| c.norm()).sum::<f64>()).fold(0.0, f64::max);
        let a_norm = nrm(&hx);
        let mut inv = vec![vec![Cx::new(0.0, 0.0); n1]; n1];
        for c in 0..n1 {
            let mut e = vec![Cx::new(0.0, 0.0); n1];
            e[c] = Cx::new(1.0, 0.0);
            match csolve(hx.clone(), e) {
                Some(col) => {
                    for r in 0..n1 {
                        inv[r][c] = col[r];
                    }
                }
                None => return f64::INFINITY,
            }
        }
        a_norm * nrm(&inv)
    }

    fn residual(&self, x: &[Cx]) -> f64 {
        let (fv, _) = eval_poly(&self.homog, x);
        let scale = 1.0 + max_norm(x).powi(self.deg.iter().copied().max().unwrap_or(1) as i32);
        max_norm(&fv) / scale
    }
}

struct Ended {
    x: Vec<Cx>,
    winding: usize,
    converged: bool,
    cond_ok: bool,
}

/// Tracks from the start point to `t = 1` with a Cauchy endgame.
fn run_path(pb: &Problem, start: &[Cx], opts: &HomotopyOptions, max_step: f64) -> Ended {
    let failed = |x: &[Cx]| Ended { x: x.to_vec(), winding: 1, converged: false, cond_ok: false };
    let s0 = opts.endgame_start.clamp(1e-3, 0.5);
    let mut x = match pb.track(start, Cx::new(0.0, 0.0), Cx::new(1.0 - s0, 0.0), max_step) {
        Some(v) => v,
        None => return failed(start),
    };
    let lambda = 0.5;
    let ver = 24usize;
    let mut s = s0;
    let mut prev: Option<(Vec<Cx>, usize)> = None;
    let mut best: Option<(Vec<Cx>, usize)> = None;
    for _ in 0..40 {
        // circle |1 - t| = s, counter-clockwise, vertices t_j = 1 - s e^{2 pi i j / ver}
        let base = x.clone();
        let mut pts: Vec<Vec<Cx>> = vec![x.clone()];
        let mut cur = x.clone();
        let mut winding = 0usize;
        let mut closed = false;
        let mut failed_loop = false;
        'loops: for _ in 0..16 {
            winding += 1;
            for j in 0..ver {
                let a = 1.0 - s * Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * j as f64 / ver as f64);
                let b = 1.0 - s * Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * (j + 1) as f64 / ver as f64);
                match pb.track(&cur, a, b, 0.5) {
                    Some(v) => cur = v,
                    None => {
                        failed_loop = true;
                        break 'loops;
                    }
                }
                pts.push(cur.clone());
            }
            let d = (0..cur.len()).map(|i| (cur[i] - base[i]).norm()).fold(0.0, f64::max);
            if d <= 1e-6 * (1.0 + max_norm(&base)) {
                closed = true;
                break;
            }
        }
        if failed_loop || !closed {
            // the circle meets a branch point of another path: use a smaller one
            let snew = s * 0.7;
            match pb.track(&base, Cx::new(1.0 - s, 0.0), Cx::new(1.0 - snew, 0.0), max_step) {
                Some(v) => {
                    x = v;
                    s = snew;
                    if s < 1e-9 {
                        break;
                    }
                    continue;
                }
                None => break,
            }
        }
        pts.pop(); // the last vertex repeats the first
        let cnt = pts.len() as f64;
        let n1 = base.len();
        let mut est = vec![Cx::new(0.0, 0.0); n1];
        for p in &pts {
            for i in 0..n1 {
                est[i] += p[i] / cnt;
            }
        }
        best = Some((est.clone(), winding));
        if let Some((pe, pw)) = &prev {
            let d = (0..n1).map(|i| (est[i] - pe[i]).norm()).fold(0.0, f64::max);
            if *pw == winding && d <= opts.endgame_tol * (1.0 + max_norm(&est)) && pb.residual(&est) < 1e-8 {
                return Ended { x: est, winding, converged: true, cond_ok: false };
            }
        }
        prev = Some((est, winding));
        // next, smaller radius along the real axis
        let snew = s * lambda;
        match pb.track(&base, Cx::new(1.0 - s, 0.0), Cx::new(1.0 - snew, 0.0), max_step) {
            Some(v) => x = v,
            None => break,
        }
        s = snew;
        if s < 1e-9 {
            break;
        }
    }
    match best {
        Some((est, w)) => Ended { x: est, winding: w, converged: false, cond_ok: false },
        None => failed(&x),
    }
}

fn classify(pb: &Problem, mut ended: Ended, opts: &HomotopyOptions, retries: usize) -> PathEnd {
    let n1 = pb.n + 1;
    // Newton polish of regular endpoints on the target (t = 1) system
    if ended.winding == 1 {
        for _ in 0..6 {
            let (h, hx, _) = pb.eval(&ended.x, Cx::new(1.0, 0.0));
            let Some(dx) = csolve(hx, h) else { break };
            let nrm = max_norm(&dx);
            if nrm > 1e-4 * (1.0 + max_norm(&ended.x)) {
                break;
            }
            for i in 0..n1 {
                ended.x[i] -= dx[i];
            }
            if nrm < 1e-15 * (1.0 + max_norm(&ended.x)) {
                break;
            }
        }
        ended.cond_ok = pb.cond_estimate(&ended.x) < 1e9;
    }
    let big = max_norm(&ended.x).max(1e-300);
    let at_inf = ended.x[0].norm() <= opts.infinity_tol * big;
    let x: Vec<(f64, f64)> = if at_inf {
        vec![(f64::NAN, f64::NAN); pb.n]
    } else {
        (1..n1)
            .map(|i| {
                let v = ended.x[i] / ended.x[0];
                (v.re, v.im)
            })
            .collect()
    };
    PathEnd {
        x,
        homogeneous: ended.x.iter().map(|c| (c.re, c.im)).collect(),
        at_infinity: at_inf,
        winding: ended.winding,
        singular: ended.winding > 1 || !ended.cond_ok,
        residual: pb.residual(&ended.x),
        converged: ended.converged || (ended.winding == 1 && ended.cond_ok && pb.residual(&ended.x) < 1e-8),
        retries,
        suspect: false,
    }
}

/// Solves the square polynomial system `system[i] = 0` (terms of equation
/// `i`, exponent vectors with one entry per variable) by total-degree
/// homotopy continuation in projective coordinates with a Cauchy endgame
/// and path-jumping protection; see the module documentation.
///
/// One [`PathEnd`] is returned per start solution (the product of the
/// degrees). Regular finite endpoints are Newton-polished; singular finite
/// endpoints and endpoints at infinity are classified with their winding
/// number.
///
/// # Errors
/// [`RootError::Invalid`] for a non-square or malformed system or more than
/// 100 000 start solutions.
pub fn homotopy_solve(system: &[Vec<PolyTerm>], opts: &HomotopyOptions) -> Result<Vec<PathEnd>, RootError> {
    let n = system.len();
    if n == 0 || system.iter().any(|eq| eq.is_empty() || eq.iter().any(|t| t.exps.len() != n)) {
        return Err(RootError::Invalid);
    }
    let deg: Vec<u32> = system
        .iter()
        .map(|eq| eq.iter().map(|t| t.exps.iter().sum::<u32>()).max().unwrap_or(0))
        .collect();
    if deg.contains(&0) {
        return Err(RootError::Invalid);
    }
    let total: usize = deg.iter().map(|&d| d as usize).product();
    if total > 100_000 {
        return Err(RootError::Invalid);
    }
    let homog: Vec<Vec<PolyTerm>> = system
        .iter()
        .zip(&deg)
        .map(|(eq, &d)| {
            eq.iter()
                .map(|t| {
                    let mut e = vec![d - t.exps.iter().sum::<u32>()];
                    e.extend_from_slice(&t.exps);
                    PolyTerm { coef: t.coef, exps: e }
                })
                .collect()
        })
        .collect();
    let mut rng = Rng::new(opts.seed);
    let gamma = Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * rng.uniform());
    let patch: Vec<Cx> = (0..=n)
        .map(|_| Cx::new(rng.uniform() * 2.0 - 1.0, rng.uniform() * 2.0 - 1.0))
        .collect();
    let pb = Problem { n, homog, deg: deg.clone(), gamma, patch };
    let starts: Vec<Vec<Cx>> = (0..total)
        .map(|mut idx| {
            let mut x = vec![Cx::new(1.0, 0.0)];
            for &d in &deg {
                let k = idx % d as usize;
                idx /= d as usize;
                x.push(Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * k as f64 / f64::from(d)));
            }
            let s: Cx = (0..=n).map(|i| pb.patch[i] * x[i]).sum();
            x.iter().map(|v| v / s).collect()
        })
        .collect();
    let mut out: Vec<PathEnd> = starts
        .iter()
        .map(|s| classify(&pb, run_path(&pb, s, opts, opts.max_step), opts, 0))
        .collect();
    // path-jumping protection: regular finite endpoints must be distinct
    let regular = |p: &PathEnd| p.converged && !p.singular && !p.at_infinity && p.winding == 1;
    let dist = |a: &PathEnd, b: &PathEnd| {
        a.x.iter().zip(&b.x).map(|(u, v)| (u.0 - v.0).hypot(u.1 - v.1)).fold(0.0, f64::max)
    };
    let mut step = opts.max_step;
    for round in 1..=opts.max_retries {
        let mut dup: Vec<usize> = Vec::new();
        for i in 0..out.len() {
            for j in i + 1..out.len() {
                if regular(&out[i]) && regular(&out[j]) {
                    let sc = 1.0 + out[i].x.iter().map(|c| c.0.hypot(c.1)).fold(0.0, f64::max);
                    if dist(&out[i], &out[j]) < 1e-7 * sc {
                        dup.push(i);
                        dup.push(j);
                    }
                }
            }
        }
        dup.sort_unstable();
        dup.dedup();
        if dup.is_empty() {
            break;
        }
        step *= 0.2;
        for &i in &dup {
            out[i] = classify(&pb, run_path(&pb, &starts[i], opts, step), opts, round);
        }
    }
    for i in 0..out.len() {
        for j in 0..out.len() {
            if i != j && regular(&out[i]) && regular(&out[j]) {
                let sc = 1.0 + out[i].x.iter().map(|c| c.0.hypot(c.1)).fold(0.0, f64::max);
                if dist(&out[i], &out[j]) < 1e-7 * sc {
                    out[i].suspect = true;
                }
            }
        }
    }
    Ok(out)
}
