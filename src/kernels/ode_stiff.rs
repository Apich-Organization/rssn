//! # Stiff ODE integrators: Radau IIA(5/9/13) and variable-order BDF
//!
//! * [`radau5`] / [`radau5_jac`]: the three-stage, order-5, L-stable Radau IIA
//!   collocation method with simplified Newton iteration (one LU
//!   factorisation of the stage matrix reused across Newton iterations and
//!   steps), an embedded error estimator and adaptive step control;
//! * [`radau`]: the general driver with 3, 5 and 7 stages (orders 5, 9, 13)
//!   selected adaptively from the Newton convergence rate, dense output (the
//!   collocation polynomial), direction-aware/terminal event detection and
//!   mass matrices `M y' = f` (index-1 DAEs, singular `M` allowed);
//!   [`radau_fixed`] / [`radau5_fixed`] are constant-step variants used to
//!   verify the convergence orders;
//! * [`bdf`] / [`bdf_jac`] / [`bdf_ex`]: variable-order (1 to 5) backward
//!   differentiation formulas (fixed-leading-coefficient form, NDF
//!   modification, adaptive order/step, Jacobian/LU reuse); [`bdf_ex`] adds
//!   dense output (backward-difference interpolant) and event detection.
//!
//! Limitations: the stage systems are solved with one dense `s n x s n` LU
//! (no sparse or banded Jacobians); DAEs are limited to index 1; the error
//! estimate of the DAE case does not exclude algebraic components.

use crate::kernels::dense::{self, Lu, Mat};
use crate::kernels::ode_adaptive::{
    DenseKind, DenseStep, Event, OdeError, OdeOptions, OdeSolution, StepWindow, Tol, initial_step,
    process_events,
};

/// User-supplied Jacobian `J(t, y)` of the right-hand side.
pub type JacobianFn<'a> = dyn Fn(f64, &[f64]) -> Mat + 'a;

/// Integration interval `[t0, t1]` and initial state `y0` of an initial-value problem.
#[derive(Clone, Copy)]
struct IvpSpan<'a> {
    t0: f64,
    t1: f64,
    y0: &'a [f64],
}

/// Fixed inputs of one simplified-Newton solve for the Radau stage increments.
#[derive(Clone, Copy)]
struct NewtonCtx<'a> {
    tab: &'a RadauTab,
    mass: Option<&'a Mat>,
    lu: &'a Lu,
    scale: &'a [f64],
    tol: f64,
    maxit: usize,
}

fn rms(v: &[f64], scale: &[f64]) -> f64 {
    let n = scale.len().max(1);
    let s: f64 = v.iter().enumerate().map(|(i, x)| (x / scale[i % n]).powi(2)).sum();
    (s / v.len().max(1) as f64).sqrt()
}

fn fd_jacobian<F: Fn(f64, &[f64], &mut [f64])>(f: &F, t: f64, y: &[f64], f0: &[f64]) -> Mat {
    let n = y.len();
    let mut jac = Mat::zeros(n, n);
    let mut yp = y.to_vec();
    let mut fp = vec![0.0; n];
    for j in 0..n {
        let dj = f64::EPSILON.sqrt() * y[j].abs().max(1e-5);
        yp[j] = y[j] + dj;
        f(t, &yp, &mut fp);
        yp[j] = y[j];
        for i in 0..n {
            jac.set(i, j, (fp[i] - f0[i]) / dj);
        }
    }
    jac
}

/// Jacobian from the user callback or, without one, by finite differences.
fn jacobian_of<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    jac_fn: Option<&JacobianFn<'_>>,
    t: f64,
    y: &[f64],
    f0: &[f64],
    nfev: &mut usize,
) -> Mat {
    jac_fn.map_or_else(
        || {
            *nfev += y.len();
            fd_jacobian(f, t, y, f0)
        },
        |j| j(t, y),
    )
}

// ------------------------------------------------------------ Radau IIA

/// Legendre polynomial `P_k(x)` by the three-term recurrence.
fn legendre(k: usize, x: f64) -> f64 {
    let (mut p0, mut p1) = (1.0, x);
    if k == 0 {
        return 1.0;
    }
    for m in 1..k {
        let mf = m as f64;
        let p2 = ((2.0 * mf + 1.0) * x * p1 - mf * p0) / (mf + 1.0);
        p0 = p1;
        p1 = p2;
    }
    p1
}

/// Tableau of the `s`-stage Radau IIA collocation method (order `2s - 1`)
/// together with the data of the embedded error estimator.
#[derive(Debug, Clone)]
struct RadauTab {
    s: usize,
    c: Vec<f64>,
    a: Vec<Vec<f64>>,
    gamma0: f64,
    /// Weights of the stage values `Z_k` in the error estimator.
    w: Vec<f64>,
}

impl RadauTab {
    fn new(s: usize) -> Option<Self> {
        // Nodes: the roots of P_{s-1} - P_s other than x = 1, mapped to [0, 1].
        let g = |x: f64| legendre(s - 1, x) - legendre(s, x);
        let grid = 8000usize;
        let mut c = Vec::with_capacity(s);
        let mut xa = -1.0;
        let mut ga = g(xa);
        for i in 1..=grid {
            let xb = -1.0 + 2.0 * (1.0 - 1e-7) * i as f64 / grid as f64;
            let gb = g(xb);
            if ga == 0.0 {
                c.push(f64::midpoint(xa, 1.0));
            } else if ga * gb < 0.0 {
                let (mut lo, mut hi, mut glo) = (xa, xb, ga);
                for _ in 0..200 {
                    let mid = f64::midpoint(lo, hi);
                    let gm = g(mid);
                    if gm * glo <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        glo = gm;
                    }
                }
                c.push(f64::midpoint(f64::midpoint(lo, hi), 1.0));
            }
            xa = xb;
            ga = gb;
        }
        c.push(1.0);
        if c.len() != s {
            return None;
        }
        // A_ij = int_0^{c_i} l_j(tau) d tau, with l_j the Lagrange basis.
        let mut a = vec![vec![0.0; s]; s];
        for j in 0..s {
            let mut poly = vec![1.0];
            for k in 0..s {
                if k == j {
                    continue;
                }
                let den = c[j] - c[k];
                let mut np = vec![0.0; poly.len() + 1];
                for (m, pm) in poly.iter().enumerate() {
                    np[m + 1] += pm / den;
                    np[m] -= pm * c[k] / den;
                }
                poly = np;
            }
            for i in 0..s {
                a[i][j] = poly
                    .iter()
                    .enumerate()
                    .map(|(m, pm)| pm * c[i].powi(m as i32 + 1) / (m as f64 + 1.0))
                    .sum();
            }
        }
        let am = Mat::from_rows(&a).ok()?;
        let ainv = dense::lu_factor(&am).ok()?.inverse();
        let eig = dense::eigenvalues(&ainv).ok()?;
        let mu = eig.iter().min_by(|x, y| x.1.abs().total_cmp(&y.1.abs()))?.0;
        if mu <= 0.0 {
            return None;
        }
        let gamma0 = 1.0 / mu;
        // Embedded weights: gamma0 f(y0) + sum bh_j f_j integrates degree < s exactly.
        let mut v = Mat::zeros(s, s);
        let mut rhs = vec![0.0; s];
        for k in 0..s {
            for j in 0..s {
                v.set(k, j, c[j].powi(k as i32));
            }
            rhs[k] = if k == 0 { 1.0 - gamma0 } else { 1.0 / (k as f64 + 1.0) };
        }
        let bh = dense::solve(&v, &rhs).ok()?;
        let w: Vec<f64> = (0..s)
            .map(|k| (0..s).map(|j| (bh[j] - a[s - 1][j]) * ainv.at(j, k)).sum())
            .collect();
        Some(Self { s, c, a, gamma0, w })
    }
}

/// Builds `I (x) M - h (A (x) J)` (the `s n x s n` simplified-Newton matrix);
/// `mass = None` stands for the identity.
fn radau_matrix(tab: &RadauTab, mass: Option<&Mat>, jac: &Mat, h: f64) -> Mat {
    let n = jac.rows;
    let s = tab.s;
    let mut m = Mat::zeros(s * n, s * n);
    for bi in 0..s {
        for bj in 0..s {
            for p in 0..n {
                for q in 0..n {
                    let mut v = -h * tab.a[bi][bj] * jac.at(p, q);
                    if bi == bj {
                        v += mass.map_or_else(|| f64::from(u8::from(p == q)), |mm| mm.at(p, q));
                    }
                    m.set(bi * n + p, bj * n + q, v);
                }
            }
        }
    }
    m
}

/// Simplified Newton iteration for the stage increments `Z`
/// (`M Z_i = h sum_j a_ij f(t + c_j h, y + Z_j)`), started at `z`.
/// Returns the increments, the iteration count and the observed rate.
fn radau_newton<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
    ctx: &NewtonCtx<'_>,
    mut z: Vec<f64>,
    nfev: &mut usize,
) -> Option<(Vec<f64>, usize, f64)> {
    let NewtonCtx { tab, mass, lu, scale, tol, maxit } = *ctx;
    let n = y.len();
    let s = tab.s;
    let mut fv = vec![0.0; s * n];
    let mut yj = vec![0.0; n];
    let mut dn_old: Option<f64> = None;
    let mut theta: f64 = 0.0;
    for k in 0..maxit {
        for j in 0..s {
            for q in 0..n {
                yj[q] = y[q] + z[j * n + q];
            }
            f(t + tab.c[j] * h, &yj, &mut fv[j * n..(j + 1) * n]);
        }
        *nfev += s;
        if fv.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let mut rhs = vec![0.0; s * n];
        for i in 0..s {
            let mz = mass.map(|mm| mm.matvec(&z[i * n..(i + 1) * n]));
            for q in 0..n {
                let acc: f64 = (0..s).map(|j| tab.a[i][j] * fv[j * n + q]).sum();
                let lhs = mz.as_ref().map_or(z[i * n + q], |v| v[q]);
                rhs[i * n + q] = -(lhs - h * acc);
            }
        }
        let dz = lu.solve(&rhs);
        let dn = rms(&dz, scale);
        let rate = dn_old.map(|o| dn / o);
        if let Some(r) = rate {
            theta = theta.max(r);
            if r >= 1.0 || r.powi((maxit - k) as i32) / (1.0 - r) * dn > tol {
                return None;
            }
        }
        for (zi, d) in z.iter_mut().zip(&dz) {
            *zi += d;
        }
        if dn == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * dn < tol) {
            return Some((z, k + 1, theta));
        }
        dn_old = Some(dn);
    }
    None
}

/// Configuration of [`radau`]: number of stages (3, 5 or 7, i.e. method
/// orders 5, 9 and 13) and an optional mass matrix.
#[derive(Debug, Clone, Copy)]
pub struct RadauConfig<'a> {
    /// Mass matrix `M` of `M y' = f(t, y)` (`None` is the identity).
    /// `M` may be singular (index-1 differential-algebraic systems with
    /// consistent initial values).
    pub mass: Option<&'a Mat>,
    /// Lowest number of stages used (3, 5 or 7).
    pub min_stages: usize,
    /// Highest number of stages used (3, 5 or 7).
    pub max_stages: usize,
}

impl Default for RadauConfig<'_> {
    fn default() -> Self {
        Self { mass: None, min_stages: 3, max_stages: 7 }
    }
}

/// Radau IIA(5) integration of `y' = f(t, y)` with a finite-difference
/// Jacobian; see [`radau5_jac`].
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn radau5<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    let cfg = RadauConfig { mass: None, min_stages: 3, max_stages: 3 };
    radau_core(&f, None, IvpSpan { t0, t1, y0 }, opts, &[], &cfg)
}

/// Radau IIA(5) integration with a user-supplied Jacobian `jac(t, y)`.
///
/// The method is L-stable, stiffly accurate and of order 5. The stage
/// system is solved by simplified Newton iteration with a single LU
/// factorisation reused across iterations and across steps whose step size
/// is (nearly) unchanged; the Jacobian is refreshed only when the Newton
/// iteration fails to converge. The error estimate is the embedded
/// third-order formula of Hairer and Wanner.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn radau5_jac<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    let cfg = RadauConfig { mass: None, min_stages: 3, max_stages: 3 };
    radau_core(&f, Some(&jac), IvpSpan { t0, t1, y0 }, opts, &[], &cfg)
}

/// Variable-order Radau IIA integration (orders 5, 9 and 13) of
/// `M y' = f(t, y)` with dense output (the collocation polynomial) and
/// event detection.
///
/// The number of stages is raised while the simplified Newton iteration
/// converges very fast (smooth, tightly toleranced problems profit from
/// the high orders) and lowered when it converges slowly; the accepted
/// orders are reported in [`OdeSolution::orders`]. `jac` is the optional
/// Jacobian `df/dy`. Events are handled exactly as in
/// [`crate::kernels::ode_adaptive::dopri5`]: direction-aware, located on the
/// collocation polynomial, terminal events end the run.
///
/// With a singular mass matrix the system is an index-1 DAE; the initial
/// values must then be consistent, and unless `opts.h0` is given the first
/// step is `1e-5 |t1 - t0|`.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`] (also for unsupported stage counts or a mass
/// matrix of the wrong size).
#[allow(clippy::too_many_arguments)] // public IVP-solver signature (f, jac, t0, t1, y0, opts, events, cfg) is kept for API compatibility
pub fn radau<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    jac: Option<&JacobianFn<'_>>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
    cfg: &RadauConfig<'_>,
) -> Result<OdeSolution, OdeError> {
    radau_core(&f, jac, IvpSpan { t0, t1, y0 }, opts, events, cfg)
}

fn radau_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    jac_fn: Option<&JacobianFn<'_>>,
    span: IvpSpan<'_>,
    opts: &OdeOptions,
    events: &[Event<'_>],
    cfg: &RadauConfig<'_>,
) -> Result<OdeSolution, OdeError> {
    let IvpSpan { t0, t1, y0 } = span;
    let n = y0.len();
    if n == 0 || crate::kernels::exact_eq(t0, t1) {
        return Err(OdeError::Invalid);
    }
    let valid = |s: usize| matches!(s, 3 | 5 | 7);
    if !valid(cfg.min_stages) || !valid(cfg.max_stages) || cfg.min_stages > cfg.max_stages {
        return Err(OdeError::Invalid);
    }
    if let Some(m) = cfg.mass
        && (m.rows != n || m.cols != n)
    {
        return Err(OdeError::Invalid);
    }
    let mass = cfg.mass;
    let mut tabs: Vec<RadauTab> = Vec::new();
    for s in [3usize, 5, 7] {
        tabs.push(RadauTab::new(s).ok_or(OdeError::SolveFailed)?);
    }
    let tab_index = |s: usize| (s - 3) / 2;
    let maxit_of = |s: usize| match s {
        3 => 7,
        5 => 10,
        _ => 12,
    };
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let rtol = opts.rtol.max(100.0 * f64::EPSILON);
    let atol = opts.atol;
    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mut stages = cfg.min_stages;
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h_abs = match (opts.h0, mass) {
        (Some(h), _) => h.abs(),
        (None, None) => initial_step(f, t0, y0, &f0, dir, 5.0, Tol { rtol, atol }),
        (None, Some(_)) => 1e-5 * (t1 - t0).abs(),
    }
    .min(opts.hmax);
    let mut jac = jacobian_of(f, jac_fn, t, &y, &f0, &mut nfev);
    let mut current_jac = true;
    let mut lus: Option<(Lu, Lu)> = None;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev: Vec<f64> = events.iter().map(|e| (e.g)(t0, y0)).collect();
    let mut prev_dense: Option<crate::kernels::ode_adaptive::DenseStep> = None;
    let mut since_change = 0usize;
    let mut steps = 0;
    let mut rejected_last = false;
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        let min_step = 10.0 * f64::EPSILON * t.abs().max(1e-300);
        h_abs = h_abs.min(opts.hmax);
        let mut accepted = false;
        while !accepted {
            if h_abs < min_step || h_abs < 1e-300 {
                return Err(OdeError::StepLimit);
            }
            let tab = &tabs[tab_index(stages)];
            let s = tab.s;
            let mut h = h_abs * dir;
            let mut t_new = t + h;
            if (t_new - t1) * dir > 0.0 {
                t_new = t1;
                h = t_new - t;
                h_abs = h.abs();
                lus = None;
            }
            if lus.is_none() {
                let full = dense::lu_factor(&radau_matrix(tab, mass, &jac, h));
                let mut real = Mat::zeros(n, n);
                for p in 0..n {
                    for q in 0..n {
                        let m = mass.map_or_else(|| f64::from(u8::from(p == q)), |mm| mm.at(p, q));
                        real.set(p, q, m - h * tab.gamma0 * jac.at(p, q));
                    }
                }
                if let (Ok(l1), Ok(l2)) = (full, dense::lu_factor(&real)) {
                    lus = Some((l1, l2));
                } else {
                    h_abs *= 0.5;
                    sol.rejected += 1;
                    continue;
                }
            }
            let Some((lu_full, lu_real)) = lus.as_ref() else {
                return Err(OdeError::SolveFailed);
            };
            let scale0: Vec<f64> = y.iter().map(|v| atol + rtol * v.abs()).collect();
            let default_guess = || -> Vec<f64> {
                let mut z = vec![0.0; s * n];
                if mass.is_none() {
                    for i in 0..s {
                        for k in 0..n {
                            z[i * n + k] = tab.c[i] * h * f0[k];
                        }
                    }
                }
                z
            };
            let extrap_guess = || -> Option<Vec<f64>> {
                let pd = prev_dense.as_ref()?;
                let mut z = vec![0.0; s * n];
                for i in 0..s {
                    let v = crate::kernels::ode_adaptive::dense_eval(pd, t + tab.c[i] * h);
                    for k in 0..n {
                        z[i * n + k] = v[k] - y[k];
                    }
                }
                z.iter().all(|v| v.is_finite()).then_some(z)
            };
            let maxit = maxit_of(s);
            let ctx = NewtonCtx { tab, mass, lu: lu_full, scale: &scale0, tol: newton_tol, maxit };
            let mut newton = extrap_guess().and_then(|z| radau_newton(f, t, &y, h, &ctx, z, &mut nfev));
            if newton.is_none() {
                newton = radau_newton(f, t, &y, h, &ctx, default_guess(), &mut nfev);
            }
            let Some((z, iters, theta)) = newton else {
                if current_jac {
                    h_abs *= 0.5;
                } else {
                    jac = jacobian_of(f, jac_fn, t, &y, &f0, &mut nfev);
                    current_jac = true;
                }
                lus = None;
                sol.rejected += 1;
                rejected_last = true;
                continue;
            };
            let y_new: Vec<f64> = (0..n).map(|q| y[q] + z[(s - 1) * n + q]).collect();
            if y_new.iter().any(|v| !v.is_finite()) {
                return Err(OdeError::NonFinite);
            }
            // error estimate: (M - h g0 J)^{-1} (h g0 f0 + M sum_k w_k Z_k)
            let mut wz = vec![0.0; n];
            for k in 0..s {
                for q in 0..n {
                    wz[q] += tab.w[k] * z[k * n + q];
                }
            }
            let mwz = mass.map_or_else(|| wz.clone(), |mm| mm.matvec(&wz));
            let rhs: Vec<f64> = (0..n).map(|q| h * tab.gamma0 * f0[q] + mwz[q]).collect();
            let mut err = lu_real.solve(&rhs);
            let scale: Vec<f64> =
                (0..n).map(|q| atol + rtol * y[q].abs().max(y_new[q].abs())).collect();
            let mut en = rms(&err, &scale);
            if en > 1.0 && mass.is_none() && (rejected_last || sol.t.len() == 1) {
                let yp: Vec<f64> = (0..n).map(|q| y[q] + err[q]).collect();
                let mut fp = vec![0.0; n];
                f(t, &yp, &mut fp);
                nfev += 1;
                let rhs2: Vec<f64> = (0..n).map(|q| h * tab.gamma0 * fp[q] + wz[q]).collect();
                err = lu_real.solve(&rhs2);
                en = rms(&err, &scale);
            }
            let safety = 0.9 * (2.0 * maxit as f64 + 1.0) / (2.0 * maxit as f64 + iters as f64);
            let expo = 1.0 / (s as f64 + 1.0);
            let factor = if en == 0.0 {
                10.0
            } else {
                (safety * en.powf(-expo)).clamp(0.2, 10.0)
            };
            if en > 1.0 {
                h_abs *= factor.min(0.9);
                lus = None;
                sol.rejected += 1;
                rejected_last = true;
                continue;
            }
            // accept: collocation polynomial through (0, y) and (c_j, y + Z_j)
            let mut nodes = vec![0.0];
            nodes.extend_from_slice(&tab.c);
            let mut vals = vec![y.clone()];
            for j in 0..s {
                vals.push((0..n).map(|q| y[q] + z[j * n + q]).collect());
            }
            let step = crate::kernels::ode_adaptive::DenseStep {
                t0: t,
                h,
                y0: y.clone(),
                kind: crate::kernels::ode_adaptive::DenseKind::Lagrange { nodes, vals },
            };
            prev_dense = Some(step.clone());
            let t_old = t;
            let window = StepWindow { t_old, tn: t_new, ynew: &y_new, dir };
            if process_events(events, &mut gprev, &mut sol, step, &window) {
                sol.orders.push(2 * s - 1);
                sol.nfev = nfev;
                return Ok(sol);
            }
            accepted = true;
            t = t_new;
            y = y_new;
            sol.t.push(t);
            sol.y.push(y.clone());
            sol.orders.push(2 * s - 1);
            f(t, &y, &mut f0);
            nfev += 1;
            current_jac = false;
            rejected_last = false;
            since_change += 1;
            let h_changed = !(1.0..1.2).contains(&factor);
            if h_changed {
                h_abs *= factor;
                lus = None;
            }
            // order selection by the Newton convergence rate
            if since_change >= 5 {
                if theta < 2e-3 && stages < cfg.max_stages {
                    stages += 2;
                    since_change = 0;
                    lus = None;
                } else if theta > 0.2 && stages > cfg.min_stages {
                    stages -= 2;
                    since_change = 0;
                    lus = None;
                }
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Constant-step Radau IIA(5): `n_steps` equal steps from `t0` to `t1`.
///
/// The Jacobian is a finite-difference one refreshed every step and the
/// Newton iteration is fully converged. Returns the final state; used to
/// measure the order of convergence (5).
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`] or [`OdeError::NonFinite`].
pub fn radau5_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n_steps: usize,
) -> Result<Vec<f64>, OdeError> {
    radau_fixed(f, t0, t1, y0, n_steps, 3)
}

/// Constant-step Radau IIA with `stages` (3, 5 or 7) stages, i.e. order
/// `2 stages - 1` (5, 9, 13).
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`] or [`OdeError::NonFinite`].
pub fn radau_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n_steps: usize,
    stages: usize,
) -> Result<Vec<f64>, OdeError> {
    let n = y0.len();
    if n == 0 || n_steps == 0 || crate::kernels::exact_eq(t0, t1) || !matches!(stages, 3 | 5 | 7) {
        return Err(OdeError::Invalid);
    }
    let tab = RadauTab::new(stages).ok_or(OdeError::SolveFailed)?;
    let h = (t1 - t0) / n_steps as f64;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    let mut nfev = 0;
    let scale = vec![1.0; n];
    for s in 0..n_steps {
        let t = t0 + s as f64 * h;
        f(t, &y, &mut f0);
        let jac = fd_jacobian(&f, t, &y, &f0);
        let lu = dense::lu_factor(&radau_matrix(&tab, None, &jac, h))
            .map_err(|_| OdeError::SolveFailed)?;
        let mut z0 = vec![0.0; stages * n];
        for i in 0..stages {
            for k in 0..n {
                z0[i * n + k] = tab.c[i] * h * f0[k];
            }
        }
        let (z, _, _) = radau_newton(
            &f,
            t,
            &y,
            h,
            &NewtonCtx { tab: &tab, mass: None, lu: &lu, scale: &scale, tol: 1e-13, maxit: 60 },
            z0,
            &mut nfev,
        )
            .ok_or(OdeError::SolveFailed)?;
        for q in 0..n {
            y[q] += z[(stages - 1) * n + q];
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
    }
    Ok(y)
}

// ------------------------------------------------------------------ BDF

const BDF_MAX_ORDER: usize = 5;
const MAXIT: usize = 4;
const BDF_KAPPA: [f64; 6] = [0.0, -0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0];

fn bdf_consts() -> ([f64; 6], [f64; 6], [f64; 7]) {
    let mut gamma = [0.0; 6];
    for k in 1..=BDF_MAX_ORDER {
        gamma[k] = gamma[k - 1] + 1.0 / k as f64;
    }
    let mut alpha = [0.0; 6];
    let mut errc = [0.0; 7];
    for k in 0..6 {
        alpha[k] = (1.0 - BDF_KAPPA[k]) * gamma[k];
    }
    for k in 0..7 {
        let kap = if k < 6 { BDF_KAPPA[k] } else { 0.0 };
        let gam = if k < 6 { gamma[k] } else { 0.0 };
        errc[k] = kap * gam + 1.0 / (k as f64 + 1.0);
    }
    (gamma, alpha, errc)
}

fn bdf_r(order: usize, factor: f64) -> Vec<Vec<f64>> {
    let m = order + 1;
    let mut r = vec![vec![0.0; m]; m];
    for j in 0..m {
        r[0][j] = 1.0;
    }
    for i in 1..m {
        for j in 1..m {
            r[i][j] = (i as f64 - 1.0 - factor * j as f64) / i as f64;
        }
    }
    for i in 1..m {
        for j in 0..m {
            r[i][j] *= r[i - 1][j];
        }
    }
    r
}

fn bdf_change_d(d: &mut [Vec<f64>], order: usize, factor: f64) {
    let m = order + 1;
    let r = bdf_r(order, factor);
    let u = bdf_r(order, 1.0);
    let mut ru = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in 0..m {
            ru[i][j] = (0..m).map(|k| r[i][k] * u[k][j]).sum();
        }
    }
    let n = d[0].len();
    let old: Vec<Vec<f64>> = d[..m].to_vec();
    for i in 0..m {
        for q in 0..n {
            d[i][q] = (0..m).map(|k| ru[k][i] * old[k][q]).sum();
        }
    }
}

/// Variable-order BDF (orders 1 to 5, NDF coefficients) with a
/// finite-difference Jacobian; see [`bdf_jac`].
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn bdf<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
) -> Result<OdeSolution, OdeError> {
    bdf_core(&f, None, IvpSpan { t0, t1, y0 }, opts, max_order, &[])
}

/// Variable-order BDF with a user-supplied Jacobian `jac(t, y)`.
///
/// The solution is advanced in the fixed-leading-coefficient form on a
/// backward-difference array; the order (capped at `max_order`, `1..=5`)
/// and the step size are chosen after every `order + 1` equal steps by
/// comparing the error estimates of the neighbouring orders. The
/// Jacobian and its LU factorisation are reused while the Newton
/// iteration converges.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn bdf_jac<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    bdf_core(&f, Some(&jac), IvpSpan { t0, t1, y0 }, opts, max_order, &[])
}

/// Variable-order BDF with dense output and event detection.
///
/// Identical to [`bdf_jac`] (the Jacobian is optional) but additionally
/// locates `events` on the backward-difference interpolating polynomial of
/// the last step (direction-aware; a terminal event ends the run), and
/// [`OdeSolution::eval`] interpolates anywhere in the integrated range.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
#[allow(clippy::too_many_arguments)] // public IVP-solver signature (f, jac, t0, t1, y0, opts, max_order, events) is kept for API compatibility
pub fn bdf_ex<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    jac: Option<&JacobianFn<'_>>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    bdf_core(&f, jac, IvpSpan { t0, t1, y0 }, opts, max_order, events)
}

fn bdf_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    jac_fn: Option<&JacobianFn<'_>>,
    span: IvpSpan<'_>,
    opts: &OdeOptions,
    max_order: usize,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    let IvpSpan { t0, t1, y0 } = span;
    let n = y0.len();
    if n == 0 || crate::kernels::exact_eq(t0, t1) || max_order == 0 {
        return Err(OdeError::Invalid);
    }
    let max_order = max_order.min(BDF_MAX_ORDER);
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let (gamma, alpha, errc) = bdf_consts();
    let rtol = opts.rtol.max(100.0 * f64::EPSILON);
    let atol = opts.atol;
    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h_abs = opts
        .h0
        .unwrap_or_else(|| initial_step(f, t0, y0, &f0, dir, 1.0, Tol { rtol, atol }))
        .abs()
        .min(opts.hmax);
    let mut d = vec![vec![0.0; n]; BDF_MAX_ORDER + 3];
    d[0].clone_from(&y);
    for q in 0..n {
        d[1][q] = f0[q] * h_abs * dir;
    }
    let mut jac = jacobian_of(f, jac_fn, t, &y, &f0, &mut nfev);
    let mut order = 1usize;
    let mut n_equal = 0usize;
    let mut lu: Option<Lu> = None;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev: Vec<f64> = events.iter().map(|e| (e.g)(t0, y0)).collect();
    let mut steps = 0;
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        let min_step = 10.0 * f64::EPSILON * t.abs().max(1e-300);
        if h_abs > opts.hmax {
            bdf_change_d(&mut d, order, opts.hmax / h_abs);
            n_equal = 0;
            h_abs = opts.hmax;
        } else if h_abs < min_step {
            bdf_change_d(&mut d, order, min_step / h_abs);
            n_equal = 0;
            h_abs = min_step;
        }
        let mut current_jac = false;
        let mut accepted = false;
        let mut d_inc = vec![0.0; n];
        let mut y_new = vec![0.0; n];
        let mut t_new = t;
        let mut err_norm_v = 0.0;
        let mut safety = 0.9;
        let mut scale = vec![1.0; n];
        while !accepted {
            if h_abs < min_step {
                return Err(OdeError::StepLimit);
            }
            let mut h = h_abs * dir;
            t_new = t + h;
            if (t_new - t1) * dir > 0.0 {
                t_new = t1;
                bdf_change_d(&mut d, order, (t_new - t).abs() / h_abs);
                n_equal = 0;
                lu = None;
            }
            h = t_new - t;
            h_abs = h.abs();
            let y_pred: Vec<f64> = (0..n).map(|q| (0..=order).map(|k| d[k][q]).sum()).collect();
            let sc: Vec<f64> = y_pred.iter().map(|v| atol + rtol * v.abs()).collect();
            let psi: Vec<f64> = (0..n)
                .map(|q| (1..=order).map(|k| d[k][q] * gamma[k]).sum::<f64>() / alpha[order])
                .collect();
            let cc = h / alpha[order];
            let mut converged = false;
            let mut n_iter = 0;
            while !converged {
                if lu.is_none() {
                    let mut m = Mat::identity(n);
                    for p in 0..n {
                        for q in 0..n {
                            let v = m.at(p, q) - cc * jac.at(p, q);
                            m.set(p, q, v);
                        }
                    }
                    match dense::lu_factor(&m) {
                        Ok(l) => lu = Some(l),
                        Err(_) => break,
                    }
                }
                let Some(l) = lu.as_ref() else { break };
                // Newton iteration
                let mut yy = y_pred.clone();
                let mut dd = vec![0.0; n];
                let mut dn_old: Option<f64> = None;
                let mut fv = vec![0.0; n];
                converged = false;
                for k in 0..MAXIT {
                    f(t_new, &yy, &mut fv);
                    nfev += 1;
                    n_iter = k + 1;
                    if fv.iter().any(|v| !v.is_finite()) {
                        break;
                    }
                    let rhs: Vec<f64> = (0..n).map(|q| cc * fv[q] - psi[q] - dd[q]).collect();
                    let dy = l.solve(&rhs);
                    let dn = rms(&dy, &sc);
                    let rate = dn_old.map(|o| dn / o);
                    if let Some(r) = rate
                        && (r >= 1.0 || r.powi((MAXIT - k) as i32) / (1.0 - r) * dn > newton_tol) {
                            break;
                        }
                    for q in 0..n {
                        yy[q] += dy[q];
                        dd[q] += dy[q];
                    }
                    if dn == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * dn < newton_tol) {
                        converged = true;
                        break;
                    }
                    dn_old = Some(dn);
                }
                if converged {
                    y_new = yy;
                    d_inc = dd;
                } else {
                    if current_jac || lu.is_none() {
                        break;
                    }
                    jac = jac_fn.map_or_else(
                        || {
                            let mut fp = vec![0.0; n];
                            f(t_new, &y_pred, &mut fp);
                            nfev += n + 1;
                            fd_jacobian(f, t_new, &y_pred, &fp)
                        },
                        |j| j(t_new, &y_pred),
                    );
                    lu = None;
                    current_jac = true;
                }
            }
            if !converged {
                h_abs *= 0.5;
                bdf_change_d(&mut d, order, 0.5);
                n_equal = 0;
                lu = None;
                sol.rejected += 1;
                continue;
            }
            safety = 0.9 * (2.0 * MAXIT as f64 + 1.0) / (2.0 * MAXIT as f64 + n_iter as f64);
            scale = y_new.iter().map(|v| atol + rtol * v.abs()).collect();
            let err: Vec<f64> = d_inc.iter().map(|v| errc[order] * v).collect();
            err_norm_v = rms(&err, &scale);
            if err_norm_v > 1.0 {
                let factor = (safety * err_norm_v.powf(-1.0 / (order as f64 + 1.0))).max(0.2);
                h_abs *= factor;
                bdf_change_d(&mut d, order, factor);
                n_equal = 0;
                sol.rejected += 1;
            } else {
                accepted = true;
            }
        }
        n_equal += 1;
        let t_old = t;
        let y_old = y.clone();
        let h_step = t_new - t_old;
        t = t_new;
        y = y_new;
        // jac was refreshed during this step only if Newton failed; keep it
        for q in 0..n {
            d[order + 2][q] = d_inc[q] - d[order + 1][q];
            d[order + 1][q] = d_inc[q];
        }
        for i in (0..=order).rev() {
            for q in 0..n {
                let v = d[i + 1][q];
                d[i][q] += v;
            }
        }
        let step = DenseStep {
            t0: t_old,
            h: h_step,
            y0: y_old,
            kind: DenseKind::Backward(d[..=order].to_vec()),
        };
        if process_events(events, &mut gprev, &mut sol, step, &StepWindow { t_old, tn: t, ynew: &y, dir }) {
            sol.nfev = nfev;
            return Ok(sol);
        }
        sol.t.push(t);
        sol.y.push(y.clone());
        sol.orders.push(order);
        if n_equal < order + 1 {
            continue;
        }
        let err_m = if order > 1 {
            let e: Vec<f64> = d[order].iter().map(|v| errc[order - 1] * v).collect();
            rms(&e, &scale)
        } else {
            f64::INFINITY
        };
        let err_p = if order < max_order {
            let e: Vec<f64> = d[order + 2].iter().map(|v| errc[order + 1] * v).collect();
            rms(&e, &scale)
        } else {
            f64::INFINITY
        };
        let norms = [err_m, err_norm_v, err_p];
        let mut factors = [0.0; 3];
        for (k, nv) in norms.iter().enumerate() {
            factors[k] = nv.powf(-1.0 / (order as f64 + k as f64));
        }
        let mut best = 0;
        for k in 1..3 {
            if factors[k] > factors[best] {
                best = k;
            }
        }
        order = (order + best).saturating_sub(1).clamp(1, max_order);
        let factor = (safety * factors[best]).min(10.0);
        h_abs *= factor;
        bdf_change_d(&mut d, order, factor);
        n_equal = 0;
        lu = None;
    }
    sol.nfev = nfev;
    Ok(sol)
}
