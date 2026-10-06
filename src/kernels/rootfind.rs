//! # Root finding and nonlinear systems
//!
//! * scalar: Brent's method, Halley's method, damped Newton with
//!   bracketing safeguard;
//! * polynomials: Aberth-Ehrlich simultaneous iteration (complex roots)
//!   and companion-matrix eigenvalues;
//! * systems: Newton with backtracking line search, Broyden's method,
//!   Levenberg-Marquardt, and homotopy continuation for polynomial
//!   systems of one complex variable per equation (total-degree start
//!   system) via [`homotopy_univariate`], and square multivariate
//!   polynomial systems via [`homotopy_system`] (total-degree start
//!   system, gamma trick, RK4 predictor / Newton corrector tracking);
//! * polynomial roots by the three-stage Jenkins-Traub algorithm
//!   ([`polynomial_roots_jenkins_traub`]).
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use crate::kernels::dense::{self, Mat};

/// Failure modes of the iterative solvers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootError {
    /// The interval does not bracket a sign change.
    NotBracketed,
    /// The iteration limit was reached.
    NoConvergence,
    /// A singular Jacobian / derivative was met.
    Singular,
    /// Invalid input (empty, wrong length, ...).
    Invalid,
}

/// Brent's method on a bracketing interval `[a, b]`.
///
/// # Errors
/// [`RootError::NotBracketed`] if `f(a)` and `f(b)` have the same sign.
pub fn brent<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, RootError> {
    let (mut a, mut b) = (a, b);
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if fa * fb > 0.0 {
        return Err(RootError::NotBracketed);
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let (mut c, mut fc) = (a, fa);
    let mut d = c;
    let mut mflag = true;
    for _ in 0..max_iter {
        if fb == 0.0 || (b - a).abs() < tol {
            return Ok(b);
        }
        let mut s = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };
        let lo = (3.0 * a + b) / 4.0;
        let (lo, hi) = if lo < b { (lo, b) } else { (b, lo) };
        let cond = !(s > lo && s < hi)
            || (mflag && (s - b).abs() >= (b - c).abs() / 2.0)
            || (!mflag && (s - b).abs() >= (c - d).abs() / 2.0)
            || (mflag && (b - c).abs() < tol)
            || (!mflag && (c - d).abs() < tol);
        if cond {
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }
        let fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Err(RootError::NoConvergence)
}

/// Halley's method (cubic convergence) using `f`, `f'` and `f''`.
///
/// # Errors
/// [`RootError::NoConvergence`] / [`RootError::Singular`].
pub fn halley<F, D, E>(
    f: F,
    df: D,
    d2f: E,
    x0: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F: Fn(f64) -> f64,
    D: Fn(f64) -> f64,
    E: Fn(f64) -> f64,
{
    let mut x = x0;
    for _ in 0..max_iter {
        let (fx, dfx, d2) = (f(x), df(x), d2f(x));
        if fx == 0.0 {
            return Ok(x);
        }
        let denom = 2.0 * dfx * dfx - fx * d2;
        if denom == 0.0 {
            return Err(RootError::Singular);
        }
        let dx = 2.0 * fx * dfx / denom;
        x -= dx;
        if dx.abs() <= tol * (1.0 + x.abs()) {
            return Ok(x);
        }
    }
    Err(RootError::NoConvergence)
}

/// Safeguarded Newton: Newton steps are kept inside the bracket `[a, b]`
/// (bisection is used whenever a step would leave it or fails to shrink
/// the residual enough).
///
/// # Errors
/// [`RootError::NotBracketed`] / [`RootError::NoConvergence`].
pub fn newton_bracketed<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(
    f: F,
    df: D,
    a: f64,
    b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, RootError> {
    let (mut lo, mut hi) = (a, b);
    let (flo, fhi) = (f(lo), f(hi));
    if flo * fhi > 0.0 {
        return Err(RootError::NotBracketed);
    }
    if flo > 0.0 {
        std::mem::swap(&mut lo, &mut hi);
    }
    let mut x = 0.5 * (lo + hi);
    for _ in 0..max_iter {
        let fx = f(x);
        if fx == 0.0 {
            return Ok(x);
        }
        if fx < 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let d = df(x);
        let mut nx = if d == 0.0 { f64::NAN } else { x - fx / d };
        let (mn, mx) = if lo < hi { (lo, hi) } else { (hi, lo) };
        if !nx.is_finite() || nx <= mn || nx >= mx {
            nx = 0.5 * (lo + hi);
        }
        if (nx - x).abs() <= tol * (1.0 + x.abs()) {
            return Ok(nx);
        }
        x = nx;
    }
    Err(RootError::NoConvergence)
}

/// Complex number as `(re, im)` helper operations used by the polynomial solvers.
type C = (f64, f64);

fn cmul(a: C, b: C) -> C {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

fn cdiv(a: C, b: C) -> C {
    let d = b.0 * b.0 + b.1 * b.1;
    ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
}

/// Evaluates a polynomial (coefficients in ascending order) and its
/// derivative at a complex point with Horner's scheme.
#[must_use]
pub fn poly_eval_complex(coeffs: &[f64], z: (f64, f64)) -> ((f64, f64), (f64, f64)) {
    let mut p = (0.0, 0.0);
    let mut dp = (0.0, 0.0);
    for &c in coeffs.iter().rev() {
        dp = cmul(dp, z);
        dp = (dp.0 + p.0, dp.1 + p.1);
        p = cmul(p, z);
        p.0 += c;
    }
    (p, dp)
}

/// All complex roots of a polynomial (coefficients ascending) by the
/// Aberth-Ehrlich method with Cauchy-bound initial guesses.
///
/// Leading zero coefficients are dropped; roots at zero are handled by
/// the iteration itself.
///
/// # Errors
/// [`RootError::Invalid`] for constant polynomials,
/// [`RootError::NoConvergence`] if 500 iterations do not suffice.
pub fn polynomial_roots_aberth(coeffs: &[f64]) -> Result<Vec<(f64, f64)>, RootError> {
    let mut c = coeffs.to_vec();
    while c.last().is_some_and(|&v| v == 0.0) {
        c.pop();
    }
    if c.len() < 2 {
        return Err(RootError::Invalid);
    }
    let n = c.len() - 1;
    let lead = c[n];
    let cn: Vec<f64> = c.iter().map(|v| v / lead).collect();
    let radius = 1.0 + cn[..n].iter().map(|v| v.abs()).fold(0.0, f64::max);
    let mut z: Vec<C> = (0..n)
        .map(|k| {
            let ang = 2.0 * std::f64::consts::PI * k as f64 / n as f64 + 0.4;
            (radius * 0.5 * ang.cos(), radius * 0.5 * ang.sin())
        })
        .collect();
    for _ in 0..500 {
        let mut maxstep: f64 = 0.0;
        for i in 0..n {
            let (p, dp) = poly_eval_complex(&cn, z[i]);
            if p.0 == 0.0 && p.1 == 0.0 {
                continue;
            }
            let ratio = cdiv(p, dp);
            let mut s = (0.0, 0.0);
            for j in 0..n {
                if j != i {
                    let d = (z[i].0 - z[j].0, z[i].1 - z[j].1);
                    let inv = cdiv((1.0, 0.0), d);
                    s = (s.0 + inv.0, s.1 + inv.1);
                }
            }
            let denom = (1.0 - (ratio.0 * s.0 - ratio.1 * s.1), -(ratio.0 * s.1 + ratio.1 * s.0));
            let w = cdiv(ratio, denom);
            if !w.0.is_finite() || !w.1.is_finite() {
                continue;
            }
            z[i] = (z[i].0 - w.0, z[i].1 - w.1);
            maxstep = maxstep.max(w.0.hypot(w.1) / (1.0 + z[i].0.hypot(z[i].1)));
        }
        if maxstep < 1e-15 {
            return Ok(z);
        }
    }
    // Accept if residuals are small.
    let ok = z.iter().all(|&zi| {
        let (p, _) = poly_eval_complex(&cn, zi);
        p.0.hypot(p.1) < 1e-8 * (1.0 + zi.0.hypot(zi.1)).powi(n as i32)
    });
    if ok { Ok(z) } else { Err(RootError::NoConvergence) }
}

/// Roots of a polynomial as eigenvalues of its companion matrix.
///
/// # Errors
/// [`RootError::Invalid`] for constant polynomials,
/// [`RootError::NoConvergence`] if the QR iteration fails.
pub fn polynomial_roots_companion(coeffs: &[f64]) -> Result<Vec<(f64, f64)>, RootError> {
    let mut c = coeffs.to_vec();
    while c.last().is_some_and(|&v| v == 0.0) {
        c.pop();
    }
    if c.len() < 2 {
        return Err(RootError::Invalid);
    }
    let n = c.len() - 1;
    let lead = c[n];
    let mut m = Mat::zeros(n, n);
    for i in 1..n {
        m.set(i, i - 1, 1.0);
    }
    for i in 0..n {
        m.set(i, n - 1, -c[i] / lead);
    }
    dense::eigenvalues(&m).map_err(|_| RootError::NoConvergence)
}

/// Result of a nonlinear-system solve.
#[derive(Debug, Clone, PartialEq)]
pub struct SystemSolution {
    /// The approximate root.
    pub x: Vec<f64>,
    /// Infinity norm of the residual.
    pub residual: f64,
    /// Iterations used.
    pub iterations: usize,
}

fn norm_inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |a, &b| a.max(b.abs()))
}

/// Forward-difference Jacobian of `f: R^n -> R^m` at `x`.
pub fn numerical_jacobian<F: Fn(&[f64]) -> Vec<f64>>(f: &F, x: &[f64], fx: &[f64]) -> Mat {
    let n = x.len();
    let m = fx.len();
    let mut j = Mat::zeros(m, n);
    let mut xp = x.to_vec();
    for k in 0..n {
        let h = 1.5e-8 * x[k].abs().max(1.0);
        xp[k] = x[k] + h;
        let fp = f(&xp);
        xp[k] = x[k];
        for i in 0..m {
            j.set(i, k, (fp[i] - fx[i]) / h);
        }
    }
    j
}

/// Newton's method for `f(x) = 0` with a backtracking (Armijo) line search
/// on `|f|^2`. Pass `jac = None` to use finite differences.
///
/// # Errors
/// [`RootError::Singular`] on a singular Jacobian,
/// [`RootError::NoConvergence`] when the iteration limit is hit.
pub fn newton_system<F, J>(
    f: F,
    jac: Option<J>,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SystemSolution, RootError>
where
    F: Fn(&[f64]) -> Vec<f64>,
    J: Fn(&[f64]) -> Mat,
{
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    for it in 0..max_iter {
        let r = norm_inf(&fx);
        if r <= tol {
            return Ok(SystemSolution { x, residual: r, iterations: it });
        }
        let j = jac.as_ref().map_or_else(|| numerical_jacobian(&f, &x, &fx), |g| g(&x));
        let lu = dense::lu_factor(&j).map_err(|_| RootError::Singular)?;
        let neg: Vec<f64> = fx.iter().map(|v| -v).collect();
        let step = lu.solve(&neg);
        let f0: f64 = fx.iter().map(|v| v * v).sum();
        let mut t = 1.0;
        let mut accepted = false;
        for _ in 0..40 {
            let xn: Vec<f64> = x.iter().zip(&step).map(|(a, b)| a + t * b).collect();
            let fnw = f(&xn);
            let f1: f64 = fnw.iter().map(|v| v * v).sum();
            if f1.is_finite() && f1 <= (1.0 - 2e-4 * t) * f0 {
                x = xn;
                fx = fnw;
                accepted = true;
                break;
            }
            t *= 0.5;
        }
        if !accepted {
            return Err(RootError::NoConvergence);
        }
    }
    let r = norm_inf(&fx);
    if r <= tol {
        Ok(SystemSolution { x, residual: r, iterations: max_iter })
    } else {
        Err(RootError::NoConvergence)
    }
}

/// Broyden's "good" method: finite-difference Jacobian at the start, then
/// rank-one updates of its inverse-applied form (Sherman-Morrison on `B`).
///
/// # Errors
/// As [`newton_system`].
pub fn broyden<F: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SystemSolution, RootError> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut b = dense::lu_factor(&numerical_jacobian(&f, &x, &fx))
        .map_err(|_| RootError::Singular)?
        .inverse();
    for it in 0..max_iter {
        let r = norm_inf(&fx);
        if r <= tol {
            return Ok(SystemSolution { x, residual: r, iterations: it });
        }
        let step: Vec<f64> = b.matvec(&fx).iter().map(|v| -v).collect();
        let xn: Vec<f64> = x.iter().zip(&step).map(|(a, s)| a + s).collect();
        let fnw = f(&xn);
        if !norm_inf(&fnw).is_finite() {
            return Err(RootError::NoConvergence);
        }
        let y: Vec<f64> = fnw.iter().zip(&fx).map(|(a, c)| a - c).collect();
        let by = b.matvec(&y);
        let sty: f64 = step.iter().zip(&by).map(|(s, v)| s * v).sum();
        if sty.abs() > 1e-300 {
            // B += (s - B y)(s^T B) / (s^T B y)
            let mut stb = vec![0.0; n];
            for j in 0..n {
                for i in 0..n {
                    stb[j] += step[i] * b.at(i, j);
                }
            }
            for i in 0..n {
                let u = (step[i] - by[i]) / sty;
                for j in 0..n {
                    let v = b.at(i, j) + u * stb[j];
                    b.set(i, j, v);
                }
            }
        }
        x = xn;
        fx = fnw;
    }
    let r = norm_inf(&fx);
    if r <= tol {
        Ok(SystemSolution { x, residual: r, iterations: max_iter })
    } else {
        Err(RootError::NoConvergence)
    }
}

/// Levenberg-Marquardt for nonlinear least squares `min |f(x)|^2` with a
/// finite-difference Jacobian. Works for `m >= n` residuals.
///
/// # Errors
/// [`RootError::NoConvergence`] when the iteration limit is reached.
pub fn levenberg_marquardt<F: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SystemSolution, RootError> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut cost: f64 = fx.iter().map(|v| v * v).sum();
    let mut lambda = 1e-3;
    for it in 0..max_iter {
        let j = numerical_jacobian(&f, &x, &fx);
        let jt = j.transpose();
        let g = jt.matvec(&fx);
        if norm_inf(&g) <= tol * 1e-3 || cost.sqrt() <= tol {
            return Ok(SystemSolution { x, residual: cost.sqrt(), iterations: it });
        }
        let jtj = jt.matmul(&j).map_err(|_| RootError::Invalid)?;
        let mut improved = false;
        for _ in 0..30 {
            let mut a = jtj.clone();
            for i in 0..n {
                let d = jtj.at(i, i).max(1e-12);
                a.set(i, i, jtj.at(i, i) + lambda * d);
            }
            let neg: Vec<f64> = g.iter().map(|v| -v).collect();
            let Ok(step) = dense::solve(&a, &neg) else {
                lambda *= 10.0;
                continue;
            };
            let xn: Vec<f64> = x.iter().zip(&step).map(|(a, b)| a + b).collect();
            let fnw = f(&xn);
            let c1: f64 = fnw.iter().map(|v| v * v).sum();
            if c1.is_finite() && c1 < cost {
                let small = (cost - c1) <= 1e-16 * cost.max(1e-300);
                x = xn;
                fx = fnw;
                cost = c1;
                lambda = (lambda * 0.3).max(1e-15);
                improved = true;
                if small {
                    return Ok(SystemSolution { x, residual: cost.sqrt(), iterations: it });
                }
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            return Ok(SystemSolution { x, residual: cost.sqrt(), iterations: it });
        }
    }
    Err(RootError::NoConvergence)
}

/// Homotopy continuation for a polynomial with real coefficients: tracks
/// the roots of `H(z, t) = (1-t) * gamma * (z^n - 1) + t * p(z)` from the
/// roots of unity at `t = 0` to those of `p` at `t = 1` with predictor
/// (Euler) and Newton corrector steps and adaptive step size.
///
/// This is the one-variable instance of total-degree homotopy
/// continuation; the random complex `gamma` is fixed (deterministic) so
/// generic paths do not cross.
///
/// # Errors
/// [`RootError::Invalid`] for constant polynomials,
/// [`RootError::NoConvergence`] if a path fails to be tracked.
pub fn homotopy_univariate(coeffs: &[f64]) -> Result<Vec<(f64, f64)>, RootError> {
    let mut c = coeffs.to_vec();
    while c.last().is_some_and(|&v| v == 0.0) {
        c.pop();
    }
    if c.len() < 2 {
        return Err(RootError::Invalid);
    }
    let n = c.len() - 1;
    let lead = c[n];
    let p: Vec<f64> = c.iter().map(|v| v / lead).collect();
    let gamma: C = (0.6049_f64.cos() * 1.3, 0.6049_f64.sin() * 1.3);
    let h = |z: C, t: f64| -> (C, C, C) {
        // value, dH/dz, dH/dt
        let (pv, dpv) = poly_eval_complex(&p, z);
        let mut zn = (1.0, 0.0);
        let mut zn1 = (1.0, 0.0);
        for k in 0..n {
            if k == n - 1 {
                zn1 = zn;
            }
            zn = cmul(zn, z);
        }
        let g = ((zn.0 - 1.0), zn.1);
        let dg = (n as f64 * zn1.0, n as f64 * zn1.1);
        let gg = cmul(gamma, g);
        let gd = cmul(gamma, dg);
        let val = ((1.0 - t) * gg.0 + t * pv.0, (1.0 - t) * gg.1 + t * pv.1);
        let dz = ((1.0 - t) * gd.0 + t * dpv.0, (1.0 - t) * gd.1 + t * dpv.1);
        let dt = (pv.0 - gg.0, pv.1 - gg.1);
        (val, dz, dt)
    };
    let mut roots = Vec::with_capacity(n);
    for k in 0..n {
        let ang = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
        let mut z: C = (ang.cos(), ang.sin());
        let mut t = 0.0;
        let mut dt: f64 = 0.02;
        let mut guard = 0;
        while t < 1.0 {
            guard += 1;
            if guard > 200_000 {
                return Err(RootError::NoConvergence);
            }
            let step = dt.min(1.0 - t);
            let (_, dz, dtv) = h(z, t);
            let tangent = cdiv(dtv, dz);
            let mut zn = (z.0 - step * tangent.0, z.1 - step * tangent.1);
            let tn = t + step;
            let mut ok = false;
            for _ in 0..6 {
                let (v, d, _) = h(zn, tn);
                let corr = cdiv(v, d);
                zn = (zn.0 - corr.0, zn.1 - corr.1);
                if corr.0.hypot(corr.1) < 1e-13 * (1.0 + zn.0.hypot(zn.1)) {
                    ok = true;
                    break;
                }
            }
            let jump = (zn.0 - z.0).hypot(zn.1 - z.1);
            if ok && jump < 0.3 * (1.0 + z.0.hypot(z.1)) && zn.0.is_finite() {
                z = zn;
                t = tn;
                dt = (dt * 1.5).min(0.05);
            } else {
                dt *= 0.5;
                if dt < 1e-14 {
                    return Err(RootError::NoConvergence);
                }
            }
        }
        // Final polish on the target polynomial.
        for _ in 0..5 {
            let (v, d) = poly_eval_complex(&p, z);
            if d.0 == 0.0 && d.1 == 0.0 {
                break;
            }
            let corr = cdiv(v, d);
            z = (z.0 - corr.0, z.1 - corr.1);
        }
        roots.push(z);
    }
    Ok(roots)
}

// ------------------------------------------------------------------
// Jenkins-Traub three-stage polynomial root finder
// ------------------------------------------------------------------

type Cx = num_complex::Complex64;

/// Horner division of a descending-order polynomial by `(z - s)`:
/// returns the quotient and the remainder (the value at `s`).
fn jt_quot(p: &[Cx], s: Cx) -> (Vec<Cx>, Cx) {
    let mut q = Vec::with_capacity(p.len().saturating_sub(1));
    let mut acc = Cx::new(0.0, 0.0);
    for (i, &c) in p.iter().enumerate() {
        acc = acc * s + c;
        if i + 1 < p.len() {
            q.push(acc);
        }
    }
    (q, acc)
}

fn jt_eval(p: &[Cx], s: Cx) -> Cx {
    p.iter().fold(Cx::new(0.0, 0.0), |a, &c| a * s + c)
}

/// Lower bound on the moduli of the roots (Cauchy), by bisection on
/// `sum_{i>=1} |a_i| x^i = |a_0|`; `p` is in descending order with a
/// non-zero constant term.
fn jt_lower_bound(p: &[Cx]) -> f64 {
    let n = p.len() - 1;
    let a0 = p[n].norm();
    let lead = p[0].norm();
    let hi0 = 1.0 + p[1..].iter().map(|c| c.norm() / lead).fold(0.0, f64::max);
    let g = |x: f64| {
        let mut acc = 0.0;
        for c in &p[..n] {
            acc = (acc + c.norm()) * x;
        }
        acc - a0
    };
    let (mut lo, mut hi) = (0.0, hi0.max(1e-300));
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if g(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

/// `H_{new}` from `H` at shift `s` (normalised form `qp + t * qh`).
/// Returns `(H_new, flag)`; `flag` is true when `H(s)` vanishes.
fn jt_next_h(p_q: &[Cx], pv: Cx, h: &[Cx], s: Cx) -> (Vec<Cx>, bool) {
    let eta = f64::EPSILON;
    let (qh, hv) = jt_quot(h, s);
    let n = p_q.len(); // = deg P = number of coefficients of H
    let tiny = hv.norm() <= eta * 10.0 * h[h.len() - 1].norm();
    let mut out = vec![Cx::new(0.0, 0.0); n];
    if tiny {
        for j in 1..n {
            out[j] = qh[j - 1];
        }
    } else {
        let t = -pv / hv;
        out[0] = p_q[0];
        for j in 1..n {
            out[j] = t * qh[j - 1] + p_q[j];
        }
    }
    (out, tiny)
}

/// Newton-like correction `t = -P(s)/H(s)`; returns `(t, vanished)`.
fn jt_calct(h: &[Cx], s: Cx, pv: Cx) -> (Cx, bool) {
    let hv = jt_eval(h, s);
    let tiny = hv.norm() <= f64::EPSILON * 10.0 * h[h.len() - 1].norm();
    if tiny { (Cx::new(0.0, 0.0), true) } else { (-pv / hv, false) }
}

fn jt_converged(p: &[Cx], s: Cx) -> (bool, f64) {
    let eta = f64::EPSILON;
    let mre = 2.0 * 2.0_f64.sqrt() * eta;
    let (_, pv) = jt_quot(p, s);
    let mut e = p[0].norm() * mre / (eta + mre);
    let mut acc = Cx::new(0.0, 0.0);
    for &c in p {
        acc = acc * s + c;
        e = e * s.norm() + acc.norm();
    }
    let errev = e * (eta + mre) - mre * pv.norm();
    (pv.norm() <= 20.0 * errev.max(0.0), pv.norm())
}

fn jt_variable(p: &[Cx], h: &mut Vec<Cx>, s0: Cx) -> Option<Cx> {
    let mut s = s0;
    let mut b = false;
    let mut omp = 0.0;
    let mut relstp = 0.0_f64;
    for i in 1..=10 {
        let (ok, mp) = jt_converged(p, s);
        if ok {
            return Some(s);
        }
        if i != 1 {
            if !b && mp >= omp && relstp < 0.05 {
                // iteration is stalling on a cluster: perturb the shift
                b = true;
                let tp = relstp.max(f64::EPSILON);
                let r1 = tp.sqrt();
                let sp = s + s * Cx::new(r1, r1);
                let (qp, pv) = jt_quot(p, sp);
                for _ in 0..5 {
                    let (hn, _) = jt_next_h(&qp, pv, h, sp);
                    *h = hn;
                }
                omp = mp;
                continue;
            } else if mp * 0.1 > omp {
                return None;
            }
        }
        omp = mp;
        let (qp, pv) = jt_quot(p, s);
        let (hn, _) = jt_next_h(&qp, pv, h, s);
        *h = hn;
        let (t, tiny) = jt_calct(h, s, pv);
        if !tiny {
            relstp = t.norm() / s.norm().max(1e-300);
            s += t;
        }
        if !s.re.is_finite() || !s.im.is_finite() {
            return None;
        }
    }
    let (ok, _) = jt_converged(p, s);
    if ok { Some(s) } else { None }
}

fn jt_fixed(p: &[Cx], h: &mut Vec<Cx>, s: Cx, l2: usize) -> Option<Cx> {
    let (qp, pv) = jt_quot(p, s);
    let (mut t, _) = jt_calct(h, s, pv);
    let mut pasd = false;
    for j in 1..=l2 {
        let ot = t;
        let (hn, _) = jt_next_h(&qp, pv, h, s);
        *h = hn;
        let (tn, tiny) = jt_calct(h, s, pv);
        t = tn;
        let z = s + t;
        if !tiny && j != l2 {
            if (t - ot).norm() >= 0.5 * z.norm() {
                pasd = false;
            } else if !pasd {
                pasd = true;
            } else {
                return jt_variable(p, h, z);
            }
        }
    }
    jt_variable(p, h, s + t)
}

/// All complex roots of a polynomial with complex coefficients
/// (`(re, im)` pairs, ascending order) by the three-stage Jenkins-Traub
/// algorithm (no-shift, fixed-shift and variable-shift stages with
/// shifts on a circle just inside the smallest root modulus), polished
/// by Newton iteration against the original polynomial.
///
/// Roots are found in order of increasing modulus and deflated forward,
/// which is the numerically stable direction. Roots at zero are
/// detected from vanishing low-order coefficients.
///
/// # Errors
/// [`RootError::Invalid`] for constant polynomials,
/// [`RootError::NoConvergence`] if the shift iteration fails.
pub fn polynomial_roots_jenkins_traub_complex(coeffs: &[(f64, f64)]) -> Result<Vec<(f64, f64)>, RootError> {
    let mut asc: Vec<Cx> = coeffs.iter().map(|&(a, b)| Cx::new(a, b)).collect();
    while asc.last().is_some_and(|c| c.norm() == 0.0) {
        asc.pop();
    }
    if asc.len() < 2 || asc.iter().any(|c| !c.re.is_finite() || !c.im.is_finite()) {
        return Err(RootError::Invalid);
    }
    let orig: Vec<Cx> = asc.iter().rev().copied().collect(); // descending
    let mut roots: Vec<Cx> = Vec::new();
    // roots at the origin
    let mut lo = 0;
    while asc[lo].norm() == 0.0 {
        roots.push(Cx::new(0.0, 0.0));
        lo += 1;
    }
    let mut p: Vec<Cx> = asc[lo..].iter().rev().copied().collect();
    let cosr = (94.0_f64).to_radians().cos();
    let sinr = (94.0_f64).to_radians().sin();
    while p.len() > 2 {
        let n = p.len() - 1;
        let bnd = jt_lower_bound(&p);
        // H^(0) = P'
        let h0: Vec<Cx> =
            (0..n).map(|i| p[i] * Cx::new((n - i) as f64 / n as f64, 0.0)).collect();
        let mut h = h0;
        // stage 1: no-shift
        for _ in 0..5 {
            let c0 = h[h.len() - 1];
            let pc = p[p.len() - 1];
            let mut hn = vec![Cx::new(0.0, 0.0); n];
            if c0.norm() > f64::EPSILON * 10.0 * pc.norm() {
                // (H - H(0)/P(0) P) / z
                let t = -c0 / pc;
                hn[0] = t * p[0];
                for j in 1..n {
                    hn[j] = h[j - 1] + t * p[j];
                }
            } else {
                for j in 1..n {
                    hn[j] = h[j - 1];
                }
            }
            h = hn;
        }
        let saved = h.clone();
        let (mut xx, mut yy) = (0.5_f64.sqrt(), -(0.5_f64.sqrt()));
        let mut found: Option<Cx> = None;
        'outer: for _pass in 0..2 {
            for cnt2 in 1..=9 {
                let nx = cosr * xx - sinr * yy;
                yy = sinr * xx + cosr * yy;
                xx = nx;
                let s = Cx::new(bnd * xx, bnd * yy);
                let mut hh = saved.clone();
                if let Some(z) = jt_fixed(&p, &mut hh, s, 10 * cnt2) {
                    found = Some(z);
                    break 'outer;
                }
            }
        }
        let Some(z) = found else {
            return Err(RootError::NoConvergence);
        };
        roots.push(z);
        let (q, _) = jt_quot(&p, z);
        p = q;
    }
    if p.len() == 2 {
        roots.push(-p[1] / p[0]);
    }
    // polish against the original polynomial
    let dorig: Vec<Cx> = (0..orig.len() - 1).map(|i| orig[i] * Cx::new((orig.len() - 1 - i) as f64, 0.0)).collect();
    for r in &mut roots {
        for _ in 0..3 {
            let v = jt_eval(&orig, *r);
            let d = jt_eval(&dorig, *r);
            if d.norm() == 0.0 {
                break;
            }
            let step = v / d;
            if step.norm() > 1e-6 * (1.0 + r.norm()) {
                break;
            }
            let cand = *r - step;
            if jt_eval(&orig, cand).norm() < v.norm() {
                *r = cand;
            } else {
                break;
            }
        }
    }
    Ok(roots.into_iter().map(|c| (c.re, c.im)).collect())
}

/// All complex roots of a real-coefficient polynomial (ascending order)
/// by the Jenkins-Traub algorithm; see
/// [`polynomial_roots_jenkins_traub_complex`].
///
/// # Errors
/// [`RootError::Invalid`] for constant polynomials,
/// [`RootError::NoConvergence`] if the shift iteration fails.
pub fn polynomial_roots_jenkins_traub(coeffs: &[f64]) -> Result<Vec<(f64, f64)>, RootError> {
    let c: Vec<(f64, f64)> = coeffs.iter().map(|&v| (v, 0.0)).collect();
    polynomial_roots_jenkins_traub_complex(&c)
}

// ------------------------------------------------------------------
// Multivariate homotopy continuation
// ------------------------------------------------------------------

/// A monomial `coef * x_0^e_0 * x_1^e_1 * ...` of a polynomial system.
#[derive(Debug, Clone, PartialEq)]
pub struct PolyTerm {
    /// Real coefficient.
    pub coef: f64,
    /// Exponent of each variable (one entry per variable).
    pub exps: Vec<u32>,
}

/// A tracked homotopy path.
#[derive(Debug, Clone, PartialEq)]
pub struct PathResult {
    /// The endpoint `(re, im)` per variable (meaningful when `finite`).
    pub x: Vec<(f64, f64)>,
    /// Max-norm of the residual of the target system at `x`.
    pub residual: f64,
    /// `true` when the path converged to a finite solution; `false` for
    /// paths diverging to infinity or failing to be tracked.
    pub finite: bool,
}

fn cpow(z: Cx, e: u32) -> Cx {
    let mut r = Cx::new(1.0, 0.0);
    for _ in 0..e {
        r *= z;
    }
    r
}

/// Evaluates a polynomial system and its Jacobian at a complex point.
fn eval_system(sys: &[Vec<PolyTerm>], x: &[Cx]) -> (Vec<Cx>, Vec<Vec<Cx>>) {
    let n = x.len();
    let mut f = vec![Cx::new(0.0, 0.0); sys.len()];
    let mut j = vec![vec![Cx::new(0.0, 0.0); n]; sys.len()];
    for (i, eq) in sys.iter().enumerate() {
        for term in eq {
            let pw: Vec<Cx> = (0..n).map(|k| cpow(x[k], term.exps[k])).collect();
            let full: Cx = pw.iter().fold(Cx::new(1.0, 0.0), |a, &b| a * b);
            f[i] += full * term.coef;
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

/// Solves the complex linear system `A x = b` by Gaussian elimination
/// with partial pivoting; `None` if (numerically) singular.
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

/// `H(x, t) = (1-t) gamma g(x) + t f(x)` with `g_i = x_i^{d_i} - 1`:
/// returns `(H, dH/dx, dH/dt)`.
fn homotopy_h(
    sys: &[Vec<PolyTerm>],
    deg: &[u32],
    gamma: Cx,
    x: &[Cx],
    t: f64,
) -> (Vec<Cx>, Vec<Vec<Cx>>, Vec<Cx>) {
    let n = x.len();
    let (fv, fj) = eval_system(sys, x);
    let mut h = vec![Cx::new(0.0, 0.0); n];
    let mut hx = vec![vec![Cx::new(0.0, 0.0); n]; n];
    let mut ht = vec![Cx::new(0.0, 0.0); n];
    for i in 0..n {
        let g = cpow(x[i], deg[i]) - 1.0;
        let dg = if deg[i] == 0 { Cx::new(0.0, 0.0) } else { cpow(x[i], deg[i] - 1) * f64::from(deg[i]) };
        h[i] = gamma * g * (1.0 - t) + fv[i] * t;
        ht[i] = fv[i] - gamma * g;
        for k in 0..n {
            hx[i][k] = fj[i][k] * t;
        }
        hx[i][i] += gamma * dg * (1.0 - t);
    }
    (h, hx, ht)
}

fn path_tangent(
    sys: &[Vec<PolyTerm>],
    deg: &[u32],
    gamma: Cx,
    x: &[Cx],
    t: f64,
) -> Option<Vec<Cx>> {
    let (_, hx, ht) = homotopy_h(sys, deg, gamma, x, t);
    let rhs: Vec<Cx> = ht.iter().map(|v| -*v).collect();
    csolve(hx, rhs)
}

fn max_norm(v: &[Cx]) -> f64 {
    v.iter().map(|c| c.norm()).fold(0.0, f64::max)
}

fn track_path(
    sys: &[Vec<PolyTerm>],
    deg: &[u32],
    gamma: Cx,
    start: &[Cx],
) -> PathResult {
    let n = start.len();
    let mut x = start.to_vec();
    let mut t = 0.0_f64;
    let mut dt = 0.01_f64;
    let fail = |x: &[Cx]| PathResult {
        x: x.iter().map(|c| (c.re, c.im)).collect(),
        residual: f64::INFINITY,
        finite: false,
    };
    let mut guard = 0usize;
    while t < 1.0 {
        guard += 1;
        if guard > 200_000 {
            return fail(&x);
        }
        let step = dt.min(1.0 - t);
        // RK4 predictor on dx/dt = -H_x^{-1} H_t
        let rk = || -> Option<Vec<Cx>> {
            let k1 = path_tangent(sys, deg, gamma, &x, t)?;
            let xa: Vec<Cx> = (0..n).map(|i| x[i] + k1[i] * (0.5 * step)).collect();
            let k2 = path_tangent(sys, deg, gamma, &xa, t + 0.5 * step)?;
            let xb: Vec<Cx> = (0..n).map(|i| x[i] + k2[i] * (0.5 * step)).collect();
            let k3 = path_tangent(sys, deg, gamma, &xb, t + 0.5 * step)?;
            let xc: Vec<Cx> = (0..n).map(|i| x[i] + k3[i] * step).collect();
            let k4 = path_tangent(sys, deg, gamma, &xc, t + step)?;
            Some((0..n).map(|i| x[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * (step / 6.0)).collect())
        };
        let tn = if step >= 1.0 - t { 1.0 } else { t + step };
        let mut ok = false;
        let mut xn = x.clone();
        if let Some(pred) = rk() {
            xn = pred;
            let mut prev = f64::INFINITY;
            for it in 0..4 {
                let (h, hx, _) = homotopy_h(sys, deg, gamma, &xn, tn);
                let Some(dx) = csolve(hx, h) else { break };
                for i in 0..n {
                    xn[i] -= dx[i];
                }
                let dn = max_norm(&dx);
                if dn < 1e-12 * (1.0 + max_norm(&xn)) {
                    ok = true;
                    break;
                }
                if it > 0 && dn > 0.5 * prev {
                    break;
                }
                prev = dn;
            }
        }
        let jump = (0..n).map(|i| (xn[i] - x[i]).norm()).fold(0.0, f64::max);
        if ok && jump <= 0.25 * (1.0 + max_norm(&x)) && xn.iter().all(|c| c.re.is_finite() && c.im.is_finite()) {
            x = xn;
            t = tn;
            dt = (dt * 1.8).min(0.1);
            if max_norm(&x) > 1e9 {
                return fail(&x);
            }
        } else {
            dt *= 0.5;
            if dt < 1e-13 {
                // Accept a near-singular endpoint if the path is already close to t = 1.
                if t > 1.0 - 1e-4 && max_norm(&x) < 1e7 {
                    break;
                }
                return fail(&x);
            }
        }
    }
    // Newton polish on the target system.
    for _ in 0..8 {
        let (fv, fj) = eval_system(sys, &x);
        let Some(dx) = csolve(fj, fv) else { break };
        let nrm = max_norm(&dx);
        if nrm > 1e-3 * (1.0 + max_norm(&x)) {
            break;
        }
        for i in 0..n {
            x[i] -= dx[i];
        }
        if nrm < 1e-15 * (1.0 + max_norm(&x)) {
            break;
        }
    }
    let (fv, _) = eval_system(sys, &x);
    let res = max_norm(&fv);
    let scale = 1.0 + max_norm(&x).powi(deg.iter().copied().max().unwrap_or(1) as i32);
    PathResult {
        x: x.iter().map(|c| (c.re, c.im)).collect(),
        residual: res,
        finite: res.is_finite() && res < 1e-7 * scale && max_norm(&x) < 1e7,
    }
}

/// Solves a square polynomial system by total-degree homotopy
/// continuation with the default (fixed, generic) `gamma`; see
/// [`homotopy_system_with_gamma`].
///
/// # Errors
/// [`RootError::Invalid`] for a non-square or malformed system.
pub fn homotopy_system(system: &[Vec<PolyTerm>]) -> Result<Vec<PathResult>, RootError> {
    let gamma = Cx::from_polar(1.0, 0.9132);
    homotopy_system_with_gamma(system, (gamma.re, gamma.im))
}

/// Total-degree homotopy continuation for a square polynomial system
/// `f_i(x_1..x_n) = 0` (`system[i]` lists the terms of `f_i`).
///
/// The start system is `g_i = x_i^{d_i} - 1` with `d_i` the total degree
/// of `f_i`; its `prod d_i` solutions are tracked along
/// `H(x, t) = (1-t) gamma g(x) + t f(x)` from `t = 0` to `t = 1` with an
/// RK4 predictor, a Newton corrector and adaptive step control. The
/// complex constant `gamma` (the "gamma trick") makes the paths
/// non-crossing with probability one. One [`PathResult`] is returned per
/// start solution; paths diverging to infinity are marked
/// `finite == false`.
///
/// # Errors
/// [`RootError::Invalid`] for a non-square or malformed system.
pub fn homotopy_system_with_gamma(
    system: &[Vec<PolyTerm>],
    gamma: (f64, f64),
) -> Result<Vec<PathResult>, RootError> {
    let n = system.len();
    if n == 0 || system.iter().any(|eq| eq.is_empty() || eq.iter().any(|t| t.exps.len() != n)) {
        return Err(RootError::Invalid);
    }
    let deg: Vec<u32> = system
        .iter()
        .map(|eq| eq.iter().map(|t| t.exps.iter().sum::<u32>()).max().unwrap_or(0))
        .collect();
    if deg.iter().any(|&d| d == 0) {
        return Err(RootError::Invalid);
    }
    let total: usize = deg.iter().map(|&d| d as usize).product();
    if total > 200_000 {
        return Err(RootError::Invalid);
    }
    let g = Cx::new(gamma.0, gamma.1);
    let mut out = Vec::with_capacity(total);
    for mut idx in 0..total {
        let mut start = Vec::with_capacity(n);
        for &d in &deg {
            let k = idx % d as usize;
            idx /= d as usize;
            start.push(Cx::from_polar(1.0, 2.0 * std::f64::consts::PI * k as f64 / f64::from(d)));
        }
        out.push(track_path(system, &deg, g, &start));
    }
    Ok(out)
}
