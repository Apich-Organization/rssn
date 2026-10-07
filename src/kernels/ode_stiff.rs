//! # Stiff ODE integrators: Radau IIA(5) and variable-order BDF
//!
//! * [`radau5`] / [`radau5_jac`]: the three-stage, order-5, L-stable Radau IIA
//!   collocation method with a simplified Newton iteration (one LU
//!   factorisation of the `3n x 3n` stage matrix reused across Newton
//!   iterations and across steps with unchanged step size), an embedded
//!   error estimator and adaptive step control; [`radau5_fixed`] is the
//!   constant-step variant used to verify the order of convergence;
//! * [`bdf`] / [`bdf_jac`]: variable-order (1 to 5) backward differentiation
//!   formulas in the fixed-leading-coefficient form on a difference array,
//!   with the NDF modification of Klopfenstein and Shampine, adaptive
//!   order and step selection and Jacobian/LU reuse.
#![allow(
    clippy::missing_const_for_fn,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::too_many_arguments,
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::suboptimal_flops,
    clippy::float_cmp,
    clippy::similar_names,
    clippy::type_complexity,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::while_float,
    clippy::needless_pass_by_value,
    clippy::option_if_let_else
)]

use crate::kernels::dense::{self, Lu, Mat};
use crate::kernels::ode_adaptive::{OdeError, OdeOptions, OdeSolution, initial_step};

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
    jac_fn: Option<&dyn Fn(f64, &[f64]) -> Mat>,
    t: f64,
    y: &[f64],
    f0: &[f64],
    nfev: &mut usize,
) -> Mat {
    if let Some(j) = jac_fn {
        j(t, y)
    } else {
        *nfev += y.len();
        fd_jacobian(f, t, y, f0)
    }
}

const RADAU_C_SQ6: f64 = 2.449489742783178; // sqrt(6)

fn radau_tableau() -> ([[f64; 3]; 3], [f64; 3]) {
    let s6 = RADAU_C_SQ6;
    let a = [
        [(88.0 - 7.0 * s6) / 360.0, (296.0 - 169.0 * s6) / 1800.0, (-2.0 + 3.0 * s6) / 225.0],
        [(296.0 + 169.0 * s6) / 1800.0, (88.0 + 7.0 * s6) / 360.0, (-2.0 - 3.0 * s6) / 225.0],
        [(16.0 - s6) / 36.0, (16.0 + s6) / 36.0, 1.0 / 9.0],
    ];
    let c = [(4.0 - s6) / 10.0, (4.0 + s6) / 10.0, 1.0];
    (a, c)
}

/// Builds `I - h (A (x) J)` (the `3n x 3n` simplified-Newton matrix).
fn radau_matrix(a: &[[f64; 3]; 3], jac: &Mat, h: f64) -> Mat {
    let n = jac.rows;
    let mut m = Mat::zeros(3 * n, 3 * n);
    for bi in 0..3 {
        for bj in 0..3 {
            for p in 0..n {
                for q in 0..n {
                    let mut v = -h * a[bi][bj] * jac.at(p, q);
                    if bi == bj && p == q {
                        v += 1.0;
                    }
                    m.set(bi * n + p, bj * n + q, v);
                }
            }
        }
    }
    m
}

/// Simplified Newton iteration for the stage increments `Z`
/// (`Z_i = h sum_j a_ij f(t + c_j h, y + Z_j)`).
fn radau_newton<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
    a: &[[f64; 3]; 3],
    c: &[f64; 3],
    f0: &[f64],
    lu: &Lu,
    scale: &[f64],
    tol: f64,
    maxit: usize,
    nfev: &mut usize,
) -> Option<(Vec<f64>, usize)> {
    let n = y.len();
    let mut z = vec![0.0; 3 * n];
    for i in 0..3 {
        for k in 0..n {
            z[i * n + k] = c[i] * h * f0[k];
        }
    }
    let mut fv = vec![0.0; 3 * n];
    let mut yj = vec![0.0; n];
    let mut dn_old: Option<f64> = None;
    for k in 0..maxit {
        for j in 0..3 {
            for q in 0..n {
                yj[q] = y[q] + z[j * n + q];
            }
            f(t + c[j] * h, &yj, &mut fv[j * n..(j + 1) * n]);
        }
        *nfev += 3;
        if fv.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let mut rhs = vec![0.0; 3 * n];
        for i in 0..3 {
            for q in 0..n {
                let s: f64 = (0..3).map(|j| a[i][j] * fv[j * n + q]).sum();
                rhs[i * n + q] = -(z[i * n + q] - h * s);
            }
        }
        let dz = lu.solve(&rhs);
        let dn = rms(&dz, scale);
        let rate = dn_old.map(|o| dn / o);
        if let Some(r) = rate
            && (r >= 1.0 || r.powi((maxit - k) as i32) / (1.0 - r) * dn > tol) {
                return None;
            }
        for (zi, d) in z.iter_mut().zip(&dz) {
            *zi += d;
        }
        if dn == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * dn < tol) {
            return Some((z, k + 1));
        }
        dn_old = Some(dn);
    }
    None
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
    radau_core(&f, None, t0, t1, y0, opts)
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
    radau_core(&f, Some(&jac), t0, t1, y0, opts)
}

fn radau_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    jac_fn: Option<&dyn Fn(f64, &[f64]) -> Mat>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    let n = y0.len();
    if n == 0 || t0 == t1 {
        return Err(OdeError::Invalid);
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let (a, c) = radau_tableau();
    let rtol = opts.rtol.max(100.0 * f64::EPSILON);
    let atol = opts.atol;
    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mu_real = 3.0 + 3.0_f64.powf(2.0 / 3.0) - 3.0_f64.cbrt();
    let e_w = [(-13.0 - 7.0 * RADAU_C_SQ6) / 3.0, (-13.0 + 7.0 * RADAU_C_SQ6) / 3.0, -1.0 / 3.0];
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h_abs = opts
        .h0
        .unwrap_or_else(|| initial_step(f, t0, y0, &f0, dir, 5.0, rtol, atol))
        .abs()
        .min(opts.hmax);
    let mut jac = jacobian_of(f, jac_fn, t, &y, &f0, &mut nfev);
    let mut current_jac = true;
    let mut lus: Option<(Lu, Lu)> = None;
    let mut sol = OdeSolution::start(t0, y0);
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
            let mut h = h_abs * dir;
            let mut t_new = t + h;
            if (t_new - t1) * dir > 0.0 {
                t_new = t1;
                h = t_new - t;
                h_abs = h.abs();
                lus = None;
            }
            if lus.is_none() {
                let full = dense::lu_factor(&radau_matrix(&a, &jac, h));
                let mut real = Mat::identity(n);
                for p in 0..n {
                    for q in 0..n {
                        let v = if p == q { mu_real / h } else { 0.0 } - jac.at(p, q);
                        real.set(p, q, v);
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
            let Some((z, iters)) =
                radau_newton(f, t, &y, h, &a, &c, &f0, lu_full, &scale0, newton_tol, 7, &mut nfev)
            else {
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
            let y_new: Vec<f64> = (0..n).map(|q| y[q] + z[2 * n + q]).collect();
            if y_new.iter().any(|v| !v.is_finite()) {
                return Err(OdeError::NonFinite);
            }
            let ze: Vec<f64> = (0..n)
                .map(|q| (e_w[0] * z[q] + e_w[1] * z[n + q] + e_w[2] * z[2 * n + q]) / h)
                .collect();
            let rhs: Vec<f64> = (0..n).map(|q| f0[q] + ze[q]).collect();
            let mut err = lu_real.solve(&rhs);
            let scale: Vec<f64> =
                (0..n).map(|q| atol + rtol * y[q].abs().max(y_new[q].abs())).collect();
            let mut en = rms(&err, &scale);
            if en > 1.0 && (rejected_last || sol.t.len() == 1) {
                let yp: Vec<f64> = (0..n).map(|q| y[q] + err[q]).collect();
                let mut fp = vec![0.0; n];
                f(t, &yp, &mut fp);
                nfev += 1;
                let rhs2: Vec<f64> = (0..n).map(|q| fp[q] + ze[q]).collect();
                err = lu_real.solve(&rhs2);
                en = rms(&err, &scale);
            }
            let safety = 0.9 * 15.0 / (14.0 + iters as f64);
            let factor = if en == 0.0 {
                10.0
            } else {
                (safety * en.powf(-0.25)).clamp(0.2, 10.0)
            };
            if en > 1.0 {
                h_abs *= factor.min(0.9);
                lus = None;
                sol.rejected += 1;
                rejected_last = true;
                continue;
            }
            // accept
            accepted = true;
            t = t_new;
            y = y_new;
            sol.t.push(t);
            sol.y.push(y.clone());
            f(t, &y, &mut f0);
            nfev += 1;
            current_jac = false;
            rejected_last = false;
            if !(1.0..1.2).contains(&factor) {
                h_abs *= factor;
                lus = None;
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
    let n = y0.len();
    if n == 0 || n_steps == 0 || t0 == t1 {
        return Err(OdeError::Invalid);
    }
    let (a, c) = radau_tableau();
    let h = (t1 - t0) / n_steps as f64;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    let mut nfev = 0;
    let scale = vec![1e-14; n];
    for s in 0..n_steps {
        let t = t0 + s as f64 * h;
        f(t, &y, &mut f0);
        let jac = fd_jacobian(&f, t, &y, &f0);
        let lu = dense::lu_factor(&radau_matrix(&a, &jac, h)).map_err(|_| OdeError::SolveFailed)?;
        let (z, _) = radau_newton(&f, t, &y, h, &a, &c, &f0, &lu, &scale, 1e-3, 30, &mut nfev)
            .ok_or(OdeError::SolveFailed)?;
        for q in 0..n {
            y[q] += z[2 * n + q];
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
    bdf_core(&f, None, t0, t1, y0, opts, max_order)
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
    bdf_core(&f, Some(&jac), t0, t1, y0, opts, max_order)
}

fn bdf_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    jac_fn: Option<&dyn Fn(f64, &[f64]) -> Mat>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
) -> Result<OdeSolution, OdeError> {
    let n = y0.len();
    if n == 0 || t0 == t1 || max_order == 0 {
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
        .unwrap_or_else(|| initial_step(f, t0, y0, &f0, dir, 1.0, rtol, atol))
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
                    jac = if let Some(j) = jac_fn {
                        j(t_new, &y_pred)
                    } else {
                        let mut fp = vec![0.0; n];
                        f(t_new, &y_pred, &mut fp);
                        nfev += n + 1;
                        fd_jacobian(f, t_new, &y_pred, &fp)
                    };
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
        t = t_new;
        y = y_new;
        sol.t.push(t);
        sol.y.push(y.clone());
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
