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
//!   successively smaller maximal steps;
//! * **extended precision**: the whole tracker is generic over the number
//!   type. With [`Precision::DoubleDouble`] (or [`Precision::Adaptive`], which
//!   re-runs only the paths that failed, ended ill-conditioned or collided
//!   with another path) every path is tracked, run through the Cauchy endgame
//!   and Newton-polished in double-double arithmetic ([`crate::kernels::dd`],
//!   about 31 digits), which resolves tight clusters of roots whose spacing
//!   makes the `f64` Jacobian numerically singular. The extra digits of the
//!   endpoints are returned in [`PathEnd::x_lo`];
//! * **certification**: [`crate::kernels::homotopy_cert`] applies Smale's
//!   alpha-theory to the endpoints (approximate-zero certificates, disjoint
//!   zero balls, completeness against the Bezout number).

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::Complex64 as Cx;

use crate::kernels::dd::Cdd;
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
    /// Arithmetic of the path tracking and endgame.
    pub precision: Precision,
}

/// Number type used for tracking, endgame and polishing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Precision {
    /// IEEE double precision throughout.
    Double,
    /// Double-double (about 31 digits) throughout.
    DoubleDouble,
    /// Double precision first; paths that did not converge, ended with an
    /// ill-conditioned Jacobian at winding number one, or coincide with
    /// another regular endpoint are re-run in double-double.
    Adaptive,
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
            precision: Precision::Double,
        }
    }
}

/// A classified path endpoint.
#[derive(Debug, Clone, PartialEq)]
#[allow(clippy::struct_excessive_bools)] // independent public status flags (at_infinity, singular, converged, suspect)
pub struct PathEnd {
    /// Affine coordinates `(re, im)` (NaN for points at infinity).
    pub x: Vec<(f64, f64)>,
    /// Homogeneous coordinates `(x_0, x_1, ..)` on the patch `a . X = 1`.
    pub homogeneous: Vec<(f64, f64)>,
    /// Low-order parts of the affine coordinates: `x[i] + x_lo[i]` is the
    /// double-double value (zeros for double precision runs).
    pub x_lo: Vec<(f64, f64)>,
    /// The path was (re-)run in double-double arithmetic.
    pub extended: bool,
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

/// Complex number type the tracker is generic over.
trait Field:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    /// Relative tolerance of a converged Newton correction.
    const NEWTON_TOL: f64;
    /// Tolerance of a stagnating Newton iteration at the rounding floor.
    const STAG_TOL: f64;
    /// Newton iterations per step.
    const NEWTON_ITERS: usize;
    /// Condition number above which an endpoint is called singular.
    const COND_LIMIT: f64;
    /// Convergence tolerance of the final polish.
    const POLISH_TOL: f64;
    /// Residual below which an endgame estimate is accepted.
    const RESIDUAL_OK: f64;
    /// Smallest endgame radius `1 - t`.
    const MIN_RADIUS: f64;
    /// Ratio of successive endgame radii.
    const SHRINK: f64;
    /// Maximal number of endgame radii.
    const ENDGAME_ROUNDS: usize;
    fn make(re: f64, im: f64) -> Self;
    fn from_cx(z: Cx) -> Self;
    fn norm(self) -> f64;
    fn scale(self, s: f64) -> Self;
    fn finite(self) -> bool;
    /// `(re_hi, im_hi, re_lo, im_lo)`.
    fn parts(self) -> (f64, f64, f64, f64);
}

impl Field for Cx {
    const NEWTON_TOL: f64 = 1e-13;
    const STAG_TOL: f64 = 1e-9;
    const NEWTON_ITERS: usize = 4;
    const COND_LIMIT: f64 = 1e9;
    const POLISH_TOL: f64 = 1e-15;
    const RESIDUAL_OK: f64 = 1e-8;
    const MIN_RADIUS: f64 = 1e-9;
    const SHRINK: f64 = 0.5;
    const ENDGAME_ROUNDS: usize = 40;
    fn make(re: f64, im: f64) -> Self {
        Self::new(re, im)
    }
    fn from_cx(z: Cx) -> Self {
        z
    }
    fn norm(self) -> f64 {
        Self::norm(self)
    }
    fn scale(self, s: f64) -> Self {
        self * s
    }
    fn finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    fn parts(self) -> (f64, f64, f64, f64) {
        (self.re, self.im, 0.0, 0.0)
    }
}

impl Field for Cdd {
    const NEWTON_TOL: f64 = 1e-27;
    const STAG_TOL: f64 = 1e-20;
    const NEWTON_ITERS: usize = 7;
    const COND_LIMIT: f64 = 1e18;
    const POLISH_TOL: f64 = 1e-30;
    const RESIDUAL_OK: f64 = 1e-27;
    const MIN_RADIUS: f64 = 1e-28;
    const SHRINK: f64 = 0.1;
    const ENDGAME_ROUNDS: usize = 80;
    fn make(re: f64, im: f64) -> Self {
        Self::new(re, im)
    }
    fn from_cx(z: Cx) -> Self {
        Self::new(z.re, z.im)
    }
    fn norm(self) -> f64 {
        Self::norm(self)
    }
    fn scale(self, s: f64) -> Self {
        self * Self::new(s, 0.0)
    }
    fn finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    fn parts(self) -> (f64, f64, f64, f64) {
        (self.re.hi, self.im.hi, self.re.lo, self.im.lo)
    }
}

struct Problem<F: Field> {
    n: usize,
    homog: Vec<Vec<PolyTerm>>,
    deg: Vec<u32>,
    gamma: F,
    patch: Vec<F>,
}

fn one<F: Field>() -> F {
    F::make(1.0, 0.0)
}

fn zero<F: Field>() -> F {
    F::make(0.0, 0.0)
}

fn cpow<F: Field>(z: F, e: u32) -> F {
    let mut r = one::<F>();
    for _ in 0..e {
        r *= z;
    }
    r
}

fn eval_poly<F: Field>(sys: &[Vec<PolyTerm>], x: &[F]) -> (Vec<F>, Vec<Vec<F>>) {
    let n = x.len();
    let mut f = vec![zero::<F>(); sys.len()];
    let mut j = vec![vec![zero::<F>(); n]; sys.len()];
    for (i, eq) in sys.iter().enumerate() {
        for term in eq {
            let pw: Vec<F> = (0..n).map(|k| cpow(x[k], term.exps[k])).collect();
            f[i] += pw.iter().fold(F::make(term.coef, 0.0), |a, &b| a * b);
            for k in 0..n {
                if term.exps[k] > 0 {
                    let mut d = F::make(term.coef * f64::from(term.exps[k]), 0.0);
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

fn csolve<F: Field>(mut a: Vec<Vec<F>>, mut b: Vec<F>) -> Option<Vec<F>> {
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
    let mut x = vec![zero::<F>(); n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    if x.iter().all(|v| v.finite()) { Some(x) } else { None }
}

fn max_norm<F: Field>(v: &[F]) -> f64 {
    v.iter().map(|c| c.norm()).fold(0.0, f64::max)
}

impl<F: Field> Problem<F> {
    /// `(H, H_x, H_t)` at `(X, t)` for complex `t`.
    fn eval(&self, x: &[F], t: F) -> (Vec<F>, Vec<Vec<F>>, Vec<F>) {
        let n = self.n;
        let (fv, fj) = eval_poly(&self.homog, x);
        let mut h = vec![zero::<F>(); n + 1];
        let mut hx = vec![vec![zero::<F>(); n + 1]; n + 1];
        let mut ht = vec![zero::<F>(); n + 1];
        let omt = one::<F>() - t;
        for i in 0..n {
            let d = self.deg[i];
            let g = cpow(x[i + 1], d) - cpow(x[0], d);
            let dgi = cpow(x[i + 1], d - 1).scale(f64::from(d));
            let dg0 = -cpow(x[0], d - 1).scale(f64::from(d));
            h[i] = self.gamma * g * omt + fv[i] * t;
            ht[i] = fv[i] - self.gamma * g;
            for k in 0..=n {
                hx[i][k] = fj[i][k] * t;
            }
            hx[i][i + 1] += self.gamma * dgi * omt;
            hx[i][0] += self.gamma * dg0 * omt;
        }
        let mut p = F::make(-1.0, 0.0);
        for k in 0..=n {
            p += self.patch[k] * x[k];
            hx[n][k] = self.patch[k];
        }
        h[n] = p;
        (h, hx, ht)
    }

    fn tangent(&self, x: &[F], t: F) -> Option<Vec<F>> {
        let (_, hx, ht) = self.eval(x, t);
        csolve(hx, ht.iter().map(|v| -*v).collect())
    }

    /// Tracks along the straight segment `ta -> tb` in the `t` plane.
    fn track(&self, x0: &[F], ta: F, tb: F, max_step: f64) -> Option<Vec<F>> {
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
            let t0 = ta + dt.scale(tau);
            let rk = || -> Option<Vec<F>> {
                let k1 = self.tangent(&x, t0)?;
                let xa: Vec<F> = (0..n1).map(|i| x[i] + (k1[i] * dt).scale(0.5 * s)).collect();
                let k2 = self.tangent(&xa, ta + dt.scale(tau + 0.5 * s))?;
                let xb: Vec<F> = (0..n1).map(|i| x[i] + (k2[i] * dt).scale(0.5 * s)).collect();
                let k3 = self.tangent(&xb, ta + dt.scale(tau + 0.5 * s))?;
                let xc: Vec<F> = (0..n1).map(|i| x[i] + (k3[i] * dt).scale(s)).collect();
                let k4 = self.tangent(&xc, ta + dt.scale(tau + s))?;
                Some(
                    (0..n1)
                        .map(|i| {
                            x[i] + ((k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]) * dt).scale(s / 6.0)
                        })
                        .collect(),
                )
            };
            let tn = if s >= 1.0 - tau { tb } else { ta + dt.scale(tau + s) };
            let mut ok = false;
            let mut xn = x.clone();
            if let Some(pred) = rk() {
                xn = pred;
                let mut prev = f64::INFINITY;
                for it in 0..F::NEWTON_ITERS {
                    let (h, hx, _) = self.eval(&xn, tn);
                    let Some(dx) = csolve(hx, h) else { break };
                    for i in 0..n1 {
                        xn[i] -= dx[i];
                    }
                    let dn = max_norm(&dx);
                    if it == 0 && dn > 2e-3 * (1.0 + max_norm(&xn)) {
                        break; // predictor too inaccurate: possible jump, shrink the step
                    }
                    if dn < F::NEWTON_TOL * (1.0 + max_norm(&xn)) {
                        ok = true;
                        break;
                    }
                    if it > 0 && dn > 0.5 * prev {
                        // stagnation at the rounding-noise floor of an ill-conditioned system
                        ok = dn < F::STAG_TOL * (1.0 + max_norm(&xn)) && dn < 2.0 * prev;
                        break;
                    }
                    prev = dn;
                }
            }
            if ok && xn.iter().all(|c| c.finite()) {
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

    fn cond_estimate(&self, x: &[F]) -> f64 {
        let (_, hx, _) = self.eval(x, one::<F>());
        let n1 = self.n + 1;
        let nrm = |m: &[Vec<F>]| m.iter().map(|r| r.iter().map(|c| c.norm()).sum::<f64>()).fold(0.0, f64::max);
        let a_norm = nrm(&hx);
        let mut inv = vec![vec![zero::<F>(); n1]; n1];
        for c in 0..n1 {
            let mut e = vec![zero::<F>(); n1];
            e[c] = one::<F>();
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

    fn residual(&self, x: &[F]) -> f64 {
        let (fv, _) = eval_poly(&self.homog, x);
        let scale = 1.0 + max_norm(x).powi(self.deg.iter().copied().max().unwrap_or(1) as i32);
        max_norm(&fv) / scale
    }
}

struct Ended<F: Field> {
    x: Vec<F>,
    winding: usize,
    converged: bool,
    cond_ok: bool,
}

/// Tracks from the start point to `t = 1` with a Cauchy endgame.
fn run_path<F: Field>(pb: &Problem<F>, start: &[F], opts: &HomotopyOptions, max_step: f64) -> Ended<F> {
    let failed = |x: &[F]| Ended { x: x.to_vec(), winding: 1, converged: false, cond_ok: false };
    let s0 = opts.endgame_start.clamp(1e-3, 0.5);
    let Some(mut x) = pb.track(start, zero::<F>(), one::<F>() - F::make(s0, 0.0), max_step) else {
        return failed(start);
    };
    let lambda = F::SHRINK;
    let ver = 24usize;
    let mut s = s0;
    let mut prev: Option<(Vec<F>, usize)> = None;
    let mut best: Option<(Vec<F>, usize)> = None;
    let endgame_tol = if F::NEWTON_ITERS > 4 { opts.endgame_tol.min(1e-20) } else { opts.endgame_tol };
    // t = 1 - s w, formed in the working precision (s may be far below 1e-16)
    let one_minus = |s: f64, w: Cx| -> F { one::<F>() - F::make(s * w.re, s * w.im) };
    let circle = |s: f64, j: usize| -> F {
        one_minus(s, Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * j as f64 / ver as f64))
    };
    let on_axis = |s: f64| -> F { one_minus(s, Cx::new(1.0, 0.0)) };
    for _ in 0..F::ENDGAME_ROUNDS {
        // circle |1 - t| = s, counter-clockwise, vertices t_j = 1 - s e^{2 pi i j / ver}
        let base = x.clone();
        let mut pts: Vec<Vec<F>> = vec![x.clone()];
        let mut cur = x.clone();
        let mut winding = 0usize;
        let mut closed = false;
        let mut failed_loop = false;
        'loops: for _ in 0..16 {
            winding += 1;
            for j in 0..ver {
                let a = circle(s, j);
                let b = circle(s, j + 1);
                if let Some(v) = pb.track(&cur, a, b, 0.5) {
                    cur = v;
                } else {
                    failed_loop = true;
                    break 'loops;
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
            match pb.track(&base, on_axis(s), on_axis(snew), max_step) {
                Some(v) => {
                    x = v;
                    s = snew;
                    if s < F::MIN_RADIUS {
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
        let mut est = vec![zero::<F>(); n1];
        for p in &pts {
            for i in 0..n1 {
                est[i] += p[i].scale(1.0 / cnt);
            }
        }
        best = Some((est.clone(), winding));
        if let Some((pe, pw)) = &prev {
            let d = (0..n1).map(|i| (est[i] - pe[i]).norm()).fold(0.0, f64::max);
            if *pw == winding && d <= endgame_tol * (1.0 + max_norm(&est)) && pb.residual(&est) < F::RESIDUAL_OK {
                return Ended { x: est, winding, converged: true, cond_ok: false };
            }
        }
        prev = Some((est, winding));
        // next, smaller radius along the real axis
        let snew = s * lambda;
        match pb.track(&base, on_axis(s), on_axis(snew), max_step) {
            Some(v) => x = v,
            None => break,
        }
        s = snew;
        if s < F::MIN_RADIUS {
            break;
        }
    }
    match best {
        Some((est, w)) => Ended { x: est, winding: w, converged: false, cond_ok: false },
        None => failed(&x),
    }
}

fn classify<F: Field>(
    pb: &Problem<F>,
    mut ended: Ended<F>,
    opts: &HomotopyOptions,
    retries: usize,
    extended: bool,
) -> PathEnd {
    let n1 = pb.n + 1;
    // Newton polish of regular endpoints on the target (t = 1) system
    if ended.winding == 1 {
        for _ in 0..8 {
            let (h, hx, _) = pb.eval(&ended.x, one::<F>());
            let Some(dx) = csolve(hx, h) else { break };
            let nrm = max_norm(&dx);
            if nrm > 1e-4 * (1.0 + max_norm(&ended.x)) {
                break;
            }
            for i in 0..n1 {
                ended.x[i] -= dx[i];
            }
            if nrm < F::POLISH_TOL * (1.0 + max_norm(&ended.x)) {
                break;
            }
        }
        ended.cond_ok = pb.cond_estimate(&ended.x) < F::COND_LIMIT;
    }
    let big = max_norm(&ended.x).max(1e-300);
    let at_inf = ended.x[0].norm() <= opts.infinity_tol * big;
    let ratios: Vec<(f64, f64, f64, f64)> = if at_inf {
        vec![(f64::NAN, f64::NAN, 0.0, 0.0); pb.n]
    } else {
        (1..n1).map(|i| (ended.x[i] / ended.x[0]).parts()).collect()
    };
    let resid = pb.residual(&ended.x);
    PathEnd {
        x: ratios.iter().map(|r| (r.0, r.1)).collect(),
        x_lo: ratios.iter().map(|r| (r.2, r.3)).collect(),
        extended,
        homogeneous: ended.x.iter().map(|c| {
            let p = c.parts();
            (p.0, p.1)
        }).collect(),
        at_infinity: at_inf,
        winding: ended.winding,
        singular: ended.winding > 1 || !ended.cond_ok,
        residual: resid,
        converged: ended.converged || (ended.winding == 1 && ended.cond_ok && resid < F::RESIDUAL_OK),
        retries,
        suspect: false,
    }
}

fn build_problem<F: Field>(n: usize, homog: &[Vec<PolyTerm>], deg: &[u32], gamma: Cx, patch: &[Cx]) -> Problem<F> {
    Problem {
        n,
        homog: homog.to_vec(),
        deg: deg.to_vec(),
        gamma: F::from_cx(gamma),
        patch: patch.iter().map(|&c| F::from_cx(c)).collect(),
    }
}

fn run_and_classify<F: Field>(
    pb: &Problem<F>,
    start: &[Cx],
    opts: &HomotopyOptions,
    step: f64,
    retries: usize,
    extended: bool,
) -> PathEnd {
    let s: Vec<F> = start.iter().map(|&c| F::from_cx(c)).collect();
    classify(pb, run_path(pb, &s, opts, step), opts, retries, extended)
}

/// Solves a square polynomial system by total-degree homotopy continuation.
///
/// Solves the square polynomial system `system[i] = 0` (terms of equation
/// `i`, exponent vectors with one entry per variable) by total-degree
/// homotopy continuation in projective coordinates with a Cauchy endgame
/// and path-jumping protection; see the module documentation.
///
/// One [`PathEnd`] is returned per start solution (the product of the
/// degrees). Regular finite endpoints are Newton-polished; singular finite
/// endpoints and endpoints at infinity are classified with their winding
/// number. [`HomotopyOptions::precision`] selects double, double-double or
/// adaptive extended-precision tracking.
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
    let pb: Problem<Cx> = build_problem(n, &homog, &deg, gamma, &patch);
    let pb_dd: Problem<Cdd> = build_problem(n, &homog, &deg, gamma, &patch);
    let starts: Vec<Vec<Cx>> = (0..total)
        .map(|mut idx| {
            let mut x = vec![Cx::new(1.0, 0.0)];
            for &d in &deg {
                let k = idx % d as usize;
                idx /= d as usize;
                x.push(Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * k as f64 / f64::from(d)));
            }
            let s: Cx = (0..=n).map(|i| patch[i] * x[i]).sum();
            x.iter().map(|v| v / s).collect()
        })
        .collect();
    let track_one = |i: usize, step: f64, retries: usize, force_dd: bool| -> PathEnd {
        match opts.precision {
            Precision::DoubleDouble => run_and_classify(&pb_dd, &starts[i], opts, step, retries, true),
            Precision::Double => run_and_classify(&pb, &starts[i], opts, step, retries, false),
            Precision::Adaptive => {
                if force_dd {
                    return run_and_classify(&pb_dd, &starts[i], opts, step, retries, true);
                }
                let r = run_and_classify(&pb, &starts[i], opts, step, retries, false);
                let ill = !r.converged || (r.singular && !r.at_infinity);
                if ill { run_and_classify(&pb_dd, &starts[i], opts, step, retries, true) } else { r }
            }
        }
    };
    let mut out: Vec<PathEnd> = (0..total).map(|i| track_one(i, opts.max_step, 0, false)).collect();
    if opts.precision == Precision::Adaptive {
        // finite endpoints that are close to another one may belong to a tight
        // cluster that double precision merged or mis-tracked: redo them
        let close: Vec<usize> = (0..out.len())
            .filter(|&i| {
                let fin = |p: &PathEnd| !p.at_infinity && !p.extended;
                fin(&out[i])
                    && (0..out.len()).any(|j| {
                        j != i && fin(&out[j]) && {
                            let sc = 1.0 + out[i].x.iter().map(|c| c.0.hypot(c.1)).fold(0.0, f64::max);
                            out[i].x.iter().zip(&out[j].x).map(|(u, v)| (u.0 - v.0).hypot(u.1 - v.1)).fold(0.0, f64::max)
                                < 1e-3 * sc
                        }
                    })
            })
            .collect();
        for i in close {
            out[i] = track_one(i, opts.max_step, 0, true);
        }
    }
    // path-jumping protection: regular finite endpoints must be distinct
    let regular = |p: &PathEnd| p.converged && !p.singular && !p.at_infinity && p.winding == 1;
    // distance including the low-order parts; two extended endpoints are
    // duplicates only if they agree to about 20 digits
    let dist = |a: &PathEnd, b: &PathEnd| {
        a.x.iter()
            .zip(&b.x)
            .zip(a.x_lo.iter().zip(&b.x_lo))
            .map(|((u, v), (ul, vl))| ((u.0 - v.0) + (ul.0 - vl.0)).hypot((u.1 - v.1) + (ul.1 - vl.1)))
            .fold(0.0, f64::max)
    };
    let dup_tol = |a: &PathEnd, b: &PathEnd| if a.extended && b.extended { 1e-20 } else { 1e-7 };
    let mut step = opts.max_step;
    for round in 1..=opts.max_retries {
        let mut dup: Vec<usize> = Vec::new();
        for i in 0..out.len() {
            for j in i + 1..out.len() {
                if regular(&out[i]) && regular(&out[j]) {
                    let sc = 1.0 + out[i].x.iter().map(|c| c.0.hypot(c.1)).fold(0.0, f64::max);
                    if dist(&out[i], &out[j]) < dup_tol(&out[i], &out[j]) * sc {
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
            out[i] = track_one(i, step, round, round > 1);
        }
    }
    for i in 0..out.len() {
        for j in 0..out.len() {
            if i != j && regular(&out[i]) && regular(&out[j]) {
                let sc = 1.0 + out[i].x.iter().map(|c| c.0.hypot(c.1)).fold(0.0, f64::max);
                if dist(&out[i], &out[j]) < dup_tol(&out[i], &out[j]) * sc {
                    out[i].suspect = true;
                }
            }
        }
    }
    Ok(out)
}
