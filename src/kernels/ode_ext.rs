//! # Exponential, IMEX, symplectic and variational integrators
//!
//! * [`etdrk4`] / [`etdrk4_diag`]: the exponential time differencing
//!   Runge-Kutta scheme of Cox and Matthews (fourth order) for
//!   `u' = L u + N(t, u)` with a stiff linear part; dense `L` uses the
//!   `phi` functions of [`crate::kernels::matfun`], diagonal `L` evaluates
//!   them by the contour-integral trick of Kassam and Trefethen;
//! * [`exprb32`]: adaptive exponential Rosenbrock method of order 3(2)
//!   (Hochbruck and Ostermann) for `y' = f(t, y)` with Jacobian
//!   linearisation;
//! * [`imex_rk`]: additive IMEX Runge-Kutta schemes for
//!   `y' = f_ex(t, y) + f_im(t, y)` (stiff part implicit): IMEX Euler,
//!   ARS(2,2,2) and ARS(3,4,3) of Ascher, Ruuth and Spiteri;
//! * [`symplectic_split`]: composition integrators of arbitrary even order
//!   for `H = T(p) + V(q)` built from Stormer-Verlet by the triple-jump
//!   construction of Yoshida and Suzuki, plus the optimised Yoshida sixth
//!   and eighth order sequences;
//! * [`variational_lobatto`]: Galerkin variational integrators of order
//!   `2 s - 2` (`s >= 2` Gauss-Lobatto nodes) for mechanical Lagrangians
//!   `L = q'^T M q' / 2 - V(q)`; they preserve the symplectic form and
//!   momentum maps exactly.
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
    clippy::similar_names,
    clippy::float_cmp,
    clippy::suboptimal_flops,
    clippy::type_complexity,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::while_float
)]

use num_complex::Complex64;

use crate::kernels::dense::{self, Mat};
use crate::kernels::matfun::{expm, phi_matrices};
use crate::kernels::ode_adaptive::{OdeError, OdeOptions, OdeSolution, err_norm};

// --------------------------------------------------------------- ETDRK4

fn mat_add_scaled(acc: &mut [f64], m: &Mat, v: &[f64], c: f64) {
    let r = m.matvec(v);
    for (a, x) in acc.iter_mut().zip(&r) {
        *a += c * x;
    }
}

/// ETDRK4 for `u' = L u + N(t, u)` with a dense matrix `L`.
///
/// Takes `steps` equal steps from `t0` to `t1` and returns the states at
/// the `steps + 1` grid points. Fourth-order accurate for any `L`,
/// independent of its stiffness.
///
/// # Errors
/// [`OdeError::Invalid`] for mismatched dimensions or `steps == 0`;
/// [`OdeError::SolveFailed`] if the matrix functions cannot be formed.
pub fn etdrk4<N: Fn(f64, &[f64], &mut [f64])>(
    l: &Mat,
    nl: N,
    t0: f64,
    t1: f64,
    u0: &[f64],
    steps: usize,
) -> Result<Vec<Vec<f64>>, OdeError> {
    let n = u0.len();
    if l.rows != n || l.cols != n || steps == 0 || n == 0 {
        return Err(OdeError::Invalid);
    }
    let h = (t1 - t0) / steps as f64;
    let mut hl = l.clone();
    hl.data.iter_mut().for_each(|x| *x *= h);
    let mut half = l.clone();
    half.data.iter_mut().for_each(|x| *x *= 0.5 * h);
    let p = phi_matrices(&hl, 3).map_err(|_| OdeError::SolveFailed)?;
    let ph = phi_matrices(&half, 1).map_err(|_| OdeError::SolveFailed)?;
    let e_full = &p[0];
    let e_half = &ph[0];
    let q = &ph[1]; // phi_1(hL/2); coefficient h/2
    let combo = |a: f64, b: f64, c: f64| -> Mat {
        let mut m = Mat::zeros(n, n);
        for i in 0..n * n {
            m.data[i] = a * p[1].data[i] + b * p[2].data[i] + c * p[3].data[i];
        }
        m
    };
    let f1 = combo(1.0, -3.0, 4.0);
    let f2 = combo(0.0, 1.0, -2.0);
    let f3 = combo(0.0, -1.0, 4.0);
    let mut u = u0.to_vec();
    let mut out = vec![u.clone()];
    let (mut nu, mut na, mut nb, mut nc) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for s in 0..steps {
        let t = t0 + s as f64 * h;
        nl(t, &u, &mut nu);
        let mut a = e_half.matvec(&u);
        mat_add_scaled(&mut a, q, &nu, 0.5 * h);
        nl(t + 0.5 * h, &a, &mut na);
        let mut b = e_half.matvec(&u);
        mat_add_scaled(&mut b, q, &na, 0.5 * h);
        nl(t + 0.5 * h, &b, &mut nb);
        let mut c = e_half.matvec(&a);
        let comb: Vec<f64> = (0..n).map(|i| 2.0 * nb[i] - nu[i]).collect();
        mat_add_scaled(&mut c, q, &comb, 0.5 * h);
        nl(t + h, &c, &mut nc);
        let mut un = e_full.matvec(&u);
        mat_add_scaled(&mut un, &f1, &nu, h);
        let sum_ab: Vec<f64> = (0..n).map(|i| na[i] + nb[i]).collect();
        mat_add_scaled(&mut un, &f2, &sum_ab, 2.0 * h);
        mat_add_scaled(&mut un, &f3, &nc, h);
        if un.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        u = un;
        out.push(u.clone());
    }
    Ok(out)
}

/// ETDRK4 for a diagonal linear part `u_i' = lambda_i u_i + N_i(t, u)`,
/// with the `phi` functions evaluated by Cauchy contour integrals
/// (Kassam and Trefethen), which avoids the cancellation of the closed
/// forms for `|h lambda| << 1`. Suited to spectral discretisations.
///
/// # Errors
/// [`OdeError::Invalid`] for mismatched dimensions or `steps == 0`.
pub fn etdrk4_diag<N: Fn(f64, &[f64], &mut [f64])>(
    lambda: &[f64],
    nl: N,
    t0: f64,
    t1: f64,
    u0: &[f64],
    steps: usize,
) -> Result<Vec<Vec<f64>>, OdeError> {
    let n = u0.len();
    if lambda.len() != n || steps == 0 || n == 0 {
        return Err(OdeError::Invalid);
    }
    let h = (t1 - t0) / steps as f64;
    let m = 64usize;
    let (mut e, mut e2, mut q, mut f1, mut f2, mut f3) =
        (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        let z = h * lambda[i];
        e[i] = z.exp();
        e2[i] = (z / 2.0).exp();
        let (mut sq, mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0, 0.0);
        for j in 0..m {
            let th = std::f64::consts::PI * (j as f64 + 0.5) / m as f64;
            let r = Complex64::new(z, 0.0) + Complex64::new(th.cos(), th.sin());
            let er = r.exp();
            let r3 = r * r * r;
            sq += (((r / 2.0).exp() - 1.0) / r).re;
            s1 += ((-4.0 - r + er * (4.0 - 3.0 * r + r * r)) / r3).re;
            s2 += ((2.0 + r + er * (-2.0 + r)) / r3).re;
            s3 += ((-4.0 - 3.0 * r - r * r + er * (4.0 - r)) / r3).re;
        }
        let mf = m as f64;
        q[i] = h * sq / mf;
        f1[i] = h * s1 / mf;
        f2[i] = h * s2 / mf;
        f3[i] = h * s3 / mf;
    }
    let mut u = u0.to_vec();
    let mut out = vec![u.clone()];
    let (mut nu, mut na, mut nb, mut nc) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for s in 0..steps {
        let t = t0 + s as f64 * h;
        nl(t, &u, &mut nu);
        let a: Vec<f64> = (0..n).map(|i| e2[i] * u[i] + q[i] * nu[i]).collect();
        nl(t + 0.5 * h, &a, &mut na);
        let b: Vec<f64> = (0..n).map(|i| e2[i] * u[i] + q[i] * na[i]).collect();
        nl(t + 0.5 * h, &b, &mut nb);
        let c: Vec<f64> = (0..n).map(|i| e2[i] * a[i] + q[i] * (2.0 * nb[i] - nu[i])).collect();
        nl(t + h, &c, &mut nc);
        let un: Vec<f64> = (0..n)
            .map(|i| e[i] * u[i] + f1[i] * nu[i] + 2.0 * f2[i] * (na[i] + nb[i]) + f3[i] * nc[i])
            .collect();
        if un.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        u = un;
        out.push(u.clone());
    }
    Ok(out)
}

// ------------------------------------------------------ Exp. Rosenbrock

/// Adaptive exponential Rosenbrock method `exprb32` (order 3 with an
/// embedded order-2 estimate) for `y' = f(t, y)`.
///
/// Each step linearises `f` at the current point (Jacobian `jac(t, y)` or a
/// finite-difference one; time dependence is handled by augmenting the state
/// with `t`), evaluates `phi_1` and `phi_3` of `h J` and applies
/// `y+ = y + h phi_1(hJ) f + 2 h phi_3(hJ) D` with
/// `D = N(U) - N(y)`, `U = y + h phi_1(hJ) f`, `N(y) = f(y) - J y`.
/// Linear problems are integrated exactly.
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::StepLimit`], [`OdeError::NonFinite`]
/// or [`OdeError::SolveFailed`].
pub fn exprb32<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    jac: Option<&dyn Fn(f64, &[f64]) -> Mat>,
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
    let m = n + 1;
    // augmented autonomous system Y = (y, t)
    let fa = |y: &[f64], out: &mut [f64]| {
        f(y[n], &y[..n], &mut out[..n]);
        out[n] = 1.0;
    };
    let mut sol = OdeSolution::start(t0, y0);
    let mut y = y0.to_vec();
    let mut t = t0;
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1usize;
    let mut h = opts
        .h0
        .unwrap_or_else(|| crate::kernels::ode_adaptive::initial_step(&f, t0, y0, &f0, dir, 3.0, opts.rtol, opts.atol))
        .abs()
        .min(opts.hmax)
        * dir;
    let mut steps = 0;
    while (t1 - t) * dir > 0.0 {
        steps += 1;
        if steps > opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        // augmented Jacobian
        let mut ya = y.clone();
        ya.push(t);
        let mut fa0 = vec![0.0; m];
        fa(&ya, &mut fa0);
        let mut ja = Mat::zeros(m, m);
        if let Some(jf) = jac {
            let jm = jf(t, &y);
            for i in 0..n {
                for j in 0..n {
                    ja.set(i, j, jm.at(i, j));
                }
            }
            // time derivative by differences
            let dt = 1.5e-8 * t.abs().max(1e-3);
            let mut fp = vec![0.0; n];
            f(t + dt, &y, &mut fp);
            for i in 0..n {
                ja.set(i, n, (fp[i] - f0[i]) / dt);
            }
            nfev += 1;
        } else {
            let mut yp = ya.clone();
            let mut fp = vec![0.0; m];
            for j in 0..m {
                let d = 1.5e-8 * ya[j].abs().max(1e-3);
                yp[j] = ya[j] + d;
                fa(&yp, &mut fp);
                yp[j] = ya[j];
                for i in 0..n {
                    ja.set(i, j, (fp[i] - fa0[i]) / d);
                }
            }
            nfev += m;
        }
        let mut hj = ja.clone();
        hj.data.iter_mut().for_each(|x| *x *= h);
        let ph = phi_matrices(&hj, 3).map_err(|_| OdeError::SolveFailed)?;
        let nterm = |yv: &[f64]| -> Vec<f64> {
            // N(Y) = F(Y) - J Y on the augmented state
            let mut fv = vec![0.0; m];
            fa(yv, &mut fv);
            let jy = ja.matvec(yv);
            (0..m).map(|i| fv[i] - jy[i]).collect()
        };
        let u_inc = ph[1].matvec(&fa0);
        let ua: Vec<f64> = (0..m).map(|i| ya[i] + h * u_inc[i]).collect();
        let n0 = nterm(&ya);
        let n1 = nterm(&ua);
        nfev += 2;
        let d: Vec<f64> = (0..m).map(|i| n1[i] - n0[i]).collect();
        let p3d = ph[3].matvec(&d);
        let ynew: Vec<f64> = (0..n).map(|i| ua[i] + 2.0 * h * p3d[i]).collect();
        let err: Vec<f64> = (0..n).map(|i| 2.0 * h * p3d[i]).collect();
        if ynew.iter().any(|v| !v.is_finite()) {
            h *= 0.25;
            sol.rejected += 1;
            if h.abs() < 1e-300 {
                return Err(OdeError::NonFinite);
            }
            continue;
        }
        let en = err_norm(&err, &y, &ynew, opts.rtol, opts.atol);
        if en <= 1.0 {
            t += h;
            y = ynew;
            f(t, &y, &mut f0);
            nfev += 1;
            sol.t.push(t);
            sol.y.push(y.clone());
            let fac = if en == 0.0 { 5.0 } else { (0.9 * en.powf(-1.0 / 3.0)).clamp(0.2, 5.0) };
            h *= fac;
            if h.abs() > opts.hmax {
                h = opts.hmax * dir;
            }
        } else {
            sol.rejected += 1;
            h *= (0.9 * en.powf(-1.0 / 3.0)).clamp(0.1, 0.9);
            if h.abs() < 1e-14 * (1.0 + t.abs()) {
                return Err(OdeError::StepLimit);
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

// ----------------------------------------------------------------- IMEX

/// IMEX Runge-Kutta scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImexScheme {
    /// IMEX (forward/backward) Euler, order 1.
    Euler,
    /// ARS(2,2,2), order 2, L-stable implicit part.
    Ars222,
    /// ARS(3,4,3), order 3.
    Ars343,
}

struct ImexTableau {
    ae: Vec<Vec<f64>>,
    ai: Vec<Vec<f64>>,
    be: Vec<f64>,
    bi: Vec<f64>,
    c: Vec<f64>,
}

fn imex_tableau(s: ImexScheme) -> ImexTableau {
    match s {
        ImexScheme::Euler => ImexTableau {
            ae: vec![vec![0.0, 0.0], vec![1.0, 0.0]],
            ai: vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            be: vec![1.0, 0.0],
            bi: vec![0.0, 1.0],
            c: vec![0.0, 1.0],
        },
        ImexScheme::Ars222 => {
            let g = 1.0 - 1.0 / 2f64.sqrt();
            let d = 1.0 - 1.0 / (2.0 * g);
            ImexTableau {
                ae: vec![vec![0.0, 0.0, 0.0], vec![g, 0.0, 0.0], vec![d, 1.0 - d, 0.0]],
                ai: vec![vec![0.0, 0.0, 0.0], vec![0.0, g, 0.0], vec![0.0, 1.0 - g, g]],
                be: vec![d, 1.0 - d, 0.0],
                bi: vec![0.0, 1.0 - g, g],
                c: vec![0.0, g, 1.0],
            }
        }
        ImexScheme::Ars343 => {
            let g = 0.435_866_521_508_459_0;
            let b1 = -1.5 * g * g + 4.0 * g - 0.25;
            let b2 = 1.5 * g * g - 5.0 * g + 1.25;
            ImexTableau {
                ae: vec![
                    vec![0.0, 0.0, 0.0, 0.0],
                    vec![g, 0.0, 0.0, 0.0],
                    vec![0.321_278_886_028_627_8, 0.396_654_374_725_601_7, 0.0, 0.0],
                    vec![-0.105_858_296_071_880, 0.552_929_148_035_939_8, 0.552_929_148_035_939_8, 0.0],
                ],
                ai: vec![
                    vec![0.0, 0.0, 0.0, 0.0],
                    vec![0.0, g, 0.0, 0.0],
                    vec![0.0, 0.5 * (1.0 - g), g, 0.0],
                    vec![0.0, b1, b2, g],
                ],
                be: vec![0.0, b1, b2, g],
                bi: vec![0.0, b1, b2, g],
                c: vec![0.0, g, 0.5 * (1.0 + g), 1.0],
            }
        }
    }
}

/// Additive IMEX Runge-Kutta integration of
/// `y' = f_ex(t, y) + f_im(t, y)` with `steps` equal steps.
///
/// `f_ex` (non-stiff) is treated explicitly, `f_im` (stiff) implicitly: each
/// stage solves `Y - h gamma f_im(t, Y) = rhs` by Newton iteration with the
/// Jacobian `jac_im(t, y)` of `f_im` (finite differences when `None`),
/// factored once per step. Returns the states on the `steps + 1` grid
/// points.
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`] or
/// [`OdeError::NonFinite`].
pub fn imex_rk<E, I>(
    f_ex: E,
    f_im: I,
    jac_im: Option<&dyn Fn(f64, &[f64]) -> Mat>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
    scheme: ImexScheme,
) -> Result<Vec<Vec<f64>>, OdeError>
where
    E: Fn(f64, &[f64], &mut [f64]),
    I: Fn(f64, &[f64], &mut [f64]),
{
    let n = y0.len();
    if n == 0 || steps == 0 {
        return Err(OdeError::Invalid);
    }
    let tab = imex_tableau(scheme);
    let s = tab.bi.len();
    let h = (t1 - t0) / steps as f64;
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    for st in 0..steps {
        let t = t0 + st as f64 * h;
        let jm = if let Some(jf) = jac_im {
            jf(t, &y)
        } else {
            let mut f0 = vec![0.0; n];
            f_im(t, &y, &mut f0);
            let mut jm = Mat::zeros(n, n);
            let mut yp = y.clone();
            let mut fp = vec![0.0; n];
            for j in 0..n {
                let d = 1.5e-8 * y[j].abs().max(1e-3);
                yp[j] = y[j] + d;
                f_im(t, &yp, &mut fp);
                yp[j] = y[j];
                for i in 0..n {
                    jm.set(i, j, (fp[i] - f0[i]) / d);
                }
            }
            jm
        };
        // factor I - h gamma J for the (common) diagonal coefficient
        let gamma = tab.ai.iter().enumerate().map(|(i, r)| r[i]).fold(0.0, f64::max);
        let mut w = Mat::identity(n);
        for i in 0..n {
            for j in 0..n {
                w.set(i, j, w.at(i, j) - h * gamma * jm.at(i, j));
            }
        }
        let lu = dense::lu_factor(&w).map_err(|_| OdeError::SolveFailed)?;
        let mut ke: Vec<Vec<f64>> = vec![vec![0.0; n]; s];
        let mut ki: Vec<Vec<f64>> = vec![vec![0.0; n]; s];
        for i in 0..s {
            let ti = t + tab.c[i] * h;
            // explicit part of the stage value
            let mut rhs = y.clone();
            for j in 0..i {
                for q in 0..n {
                    rhs[q] += h * (tab.ae[i][j] * ke[j][q] + tab.ai[i][j] * ki[j][q]);
                }
            }
            let aii = tab.ai[i][i];
            let yi = if aii == 0.0 {
                rhs
            } else {
                // Newton: Y = rhs + h aii f_im(ti, Y)
                let mut yv = rhs.clone();
                let mut fv = vec![0.0; n];
                let mut conv = false;
                for _ in 0..30 {
                    f_im(ti, &yv, &mut fv);
                    let res: Vec<f64> = (0..n).map(|q| -(yv[q] - rhs[q] - h * aii * fv[q])).collect();
                    let dy = lu.solve(&res);
                    // the factorisation uses the largest diagonal coefficient: scale the step
                    let scale = 1.0;
                    for q in 0..n {
                        yv[q] += scale * dy[q];
                    }
                    let dn = dy.iter().map(|v| v.abs()).fold(0.0, f64::max);
                    let yn = yv.iter().map(|v| v.abs()).fold(1.0, f64::max);
                    if dn <= 1e-13 * yn {
                        conv = true;
                        break;
                    }
                }
                if !conv {
                    return Err(OdeError::SolveFailed);
                }
                yv
            };
            f_ex(ti, &yi, &mut ke[i]);
            f_im(ti, &yi, &mut ki[i]);
            if ke[i].iter().chain(&ki[i]).any(|v| !v.is_finite()) {
                return Err(OdeError::NonFinite);
            }
        }
        for q in 0..n {
            let mut acc = 0.0;
            for i in 0..s {
                acc += tab.be[i] * ke[i][q] + tab.bi[i] * ki[i][q];
            }
            y[q] += h * acc;
        }
        out.push(y.clone());
    }
    Ok(out)
}

// ------------------------------------------------------------ Symplectic

/// Composition method for separable Hamiltonians.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitMethod {
    /// Stormer-Verlet (order 2).
    Verlet,
    /// Triple-jump composition of Verlet to the given even order (4, 6, 8, 10).
    TripleJump(usize),
    /// Yoshida's optimised 7-stage sixth-order sequence.
    Yoshida6,
    /// Yoshida's optimised 15-stage eighth-order sequence.
    Yoshida8,
}

fn composition_weights(m: SplitMethod) -> Vec<f64> {
    match m {
        SplitMethod::Verlet => vec![1.0],
        SplitMethod::TripleJump(order) => {
            // recursive: S_{k+2}(h) = S_k(g h) S_k((1 - 2g) h) S_k(g h), g = 1 / (2 - 2^{1/(k+1)})
            let mut w = vec![1.0];
            let mut k = 2usize;
            while k < order {
                let g = 1.0 / (2.0 - 2f64.powf(1.0 / (k as f64 + 1.0)));
                let mut nw = Vec::with_capacity(3 * w.len());
                nw.extend(w.iter().map(|x| x * g));
                nw.extend(w.iter().map(|x| x * (1.0 - 2.0 * g)));
                nw.extend(w.iter().map(|x| x * g));
                w = nw;
                k += 2;
            }
            w
        }
        SplitMethod::Yoshida6 => {
            let (w1, w2, w3) = (-1.177_679_984_178_87, 0.235_573_213_359_357, 0.784_513_610_477_560);
            let w0 = 1.0 - 2.0 * (w1 + w2 + w3);
            vec![w3, w2, w1, w0, w1, w2, w3]
        }
        SplitMethod::Yoshida8 => {
            let ws = [
                0.102_799_849_391_985,
                -1.960_610_232_975_49,
                1.938_139_137_622_76,
                -0.158_240_635_368_243,
                -1.444_852_236_860_48,
                0.253_693_336_566_229,
                0.914_844_246_229_740,
            ];
            let w0 = 1.0 - 2.0 * ws.iter().sum::<f64>();
            let mut v: Vec<f64> = ws.iter().rev().copied().collect();
            v.push(w0);
            v.extend(ws.iter().copied());
            v
        }
    }
}

/// Integrates `q' = p`, `p' = accel(q)` (unit masses; for `H = |p|^2/2 +
/// V(q)` the acceleration is `-grad V`) with a symplectic composition
/// method and `steps` equal steps from `t0` to `t1`.
///
/// Each Verlet sub-step of the composition is time-reversible and the whole
/// map is symplectic. Returns the final `(q, p)`.
///
/// # Errors
/// [`OdeError::Invalid`] for mismatched lengths, `steps == 0` or an odd
/// or too large order in [`SplitMethod::TripleJump`].
pub fn symplectic_split<A: Fn(&[f64]) -> Vec<f64>>(
    accel: A,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    steps: usize,
    method: SplitMethod,
) -> Result<(Vec<f64>, Vec<f64>), OdeError> {
    if q0.len() != p0.len() || steps == 0 {
        return Err(OdeError::Invalid);
    }
    if let SplitMethod::TripleJump(o) = method
        && (o < 2 || o % 2 != 0 || o > 12)
    {
        return Err(OdeError::Invalid);
    }
    let w = composition_weights(method);
    let h = (t1 - t0) / steps as f64;
    let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
    for _ in 0..steps {
        for &wi in &w {
            let dt = wi * h;
            let a = accel(&q);
            for i in 0..q.len() {
                p[i] += 0.5 * dt * a[i];
            }
            for i in 0..q.len() {
                q[i] += dt * p[i];
            }
            let a = accel(&q);
            for i in 0..q.len() {
                p[i] += 0.5 * dt * a[i];
            }
        }
    }
    Ok((q, p))
}

// ------------------------------------------------------------ Variational

fn lobatto_nodes(s: usize) -> (Vec<f64>, Vec<f64>) {
    // Gauss-Lobatto nodes on [0, 1]: endpoints and the roots of P'_{s-1}
    let mut x = vec![-1.0];
    if s > 2 {
        let m = s - 1;
        // roots of P'_m by Newton from Chebyshev-Gauss-Lobatto guesses
        for k in 1..s - 1 {
            let mut z = -(std::f64::consts::PI * k as f64 / m as f64).cos();
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, z);
                for j in 1..m {
                    let jf = j as f64;
                    let p2 = ((2.0 * jf + 1.0) * z * p1 - jf * p0) / (jf + 1.0);
                    p0 = p1;
                    p1 = p2;
                }
                // P_m = p1, P_{m-1} = p0; P'_m = m (p0 - z p1) / (1 - z^2)
                let dp = m as f64 * (p0 - z * p1) / (1.0 - z * z);
                let ddp = (2.0 * z * dp - m as f64 * (m as f64 + 1.0) * p1) / (1.0 - z * z);
                let dz = dp / ddp;
                z -= dz;
                if dz.abs() < 1e-16 {
                    break;
                }
            }
            x.push(z);
        }
    }
    x.push(1.0);
    x.sort_by(f64::total_cmp);
    let m = s - 1;
    let w: Vec<f64> = x
        .iter()
        .map(|&z| {
            let (mut p0, mut p1) = (1.0, z);
            for j in 1..m {
                let jf = j as f64;
                let p2 = ((2.0 * jf + 1.0) * z * p1 - jf * p0) / (jf + 1.0);
                p0 = p1;
                p1 = p2;
            }
            let pm = if m == 0 { 1.0 } else { p1 };
            2.0 / (s as f64 * (s as f64 - 1.0) * pm * pm)
        })
        .collect();
    (x.iter().map(|z| 0.5 * (z + 1.0)).collect(), w.iter().map(|v| 0.5 * v).collect())
}

/// Galerkin variational integrator with `stages >= 2` Gauss-Lobatto nodes
/// for `L(q, q') = q'^T M q' / 2 - V(q)`.
///
/// The path inside each step is the degree-`stages - 1` polynomial through
/// the node values; the discrete action is the Lobatto quadrature of the
/// Lagrangian. The discrete Euler-Lagrange equations are solved by Newton's
/// method (finite-difference Jacobian) for the interior values and the next
/// position; momenta follow from the discrete Legendre transforms. Order
/// `2 stages - 2`; symplectic and momentum-preserving by construction.
/// Returns the final `(q, p)`.
///
/// # Errors
/// [`OdeError::Invalid`] for bad dimensions or `stages < 2`;
/// [`OdeError::SolveFailed`] if Newton's iteration does not converge.
pub fn variational_lobatto<G: Fn(&[f64]) -> Vec<f64>>(
    mass: &Mat,
    grad_v: G,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    steps: usize,
    stages: usize,
) -> Result<(Vec<f64>, Vec<f64>), OdeError> {
    let d = q0.len();
    if d == 0 || p0.len() != d || mass.rows != d || mass.cols != d || steps == 0 || stages < 2 {
        return Err(OdeError::Invalid);
    }
    let s = stages;
    let h = (t1 - t0) / steps as f64;
    let (c, w) = lobatto_nodes(s);
    // differentiation matrix D[i][j] = l_j'(c_i)
    let mut dm = vec![vec![0.0; s]; s];
    for i in 0..s {
        for j in 0..s {
            if i != j {
                let mut num = 1.0;
                for k in 0..s {
                    if k != i && k != j {
                        num *= (c[i] - c[k]) / (c[j] - c[k]);
                    }
                }
                dm[i][j] = num / (c[j] - c[i]);
            }
        }
    }
    for i in 0..s {
        dm[i][i] = -(0..s).filter(|&j| j != i).map(|j| dm[i][j]).sum::<f64>();
    }
    // G_j(Q) = dL_d/dQ_j for nodal values Q (s x d)
    let gradient = |qn: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let mut g = vec![vec![0.0; d]; s];
        for i in 0..s {
            let mut qd = vec![0.0; d];
            for j in 0..s {
                for a in 0..d {
                    qd[a] += dm[i][j] * qn[j][a] / h;
                }
            }
            let mq = mass.matvec(&qd);
            for j in 0..s {
                for a in 0..d {
                    g[j][a] += h * w[i] * mq[a] * dm[i][j] / h;
                }
            }
            let gv = grad_v(&qn[i]);
            for a in 0..d {
                g[i][a] -= h * w[i] * gv[a];
            }
        }
        g
    };
    let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
    let dim_u = (s - 1) * d;
    for _ in 0..steps {
        // unknowns: nodal values 1..s-1 (interior and the end point)
        let mut x: Vec<f64> = Vec::with_capacity(dim_u);
        for j in 1..s {
            for a in 0..d {
                // linear initial guess using the momentum
                let mp = dense::solve(mass, &p).map_err(|_| OdeError::SolveFailed)?;
                x.push(q[a] + c[j] * h * mp[a]);
            }
        }
        let residual = |x: &[f64]| -> Vec<f64> {
            let mut qn: Vec<Vec<f64>> = vec![q.clone()];
            for j in 1..s {
                qn.push(x[(j - 1) * d..j * d].to_vec());
            }
            let g = gradient(&qn);
            let mut r = Vec::with_capacity(dim_u);
            // p_k + dL_d/dQ_0 = 0 would define p; instead interior: G_j = 0, node 0: p = -G_0
            for j in 1..s - 1 {
                r.extend_from_slice(&g[j]);
            }
            for a in 0..d {
                r.push(p[a] + g[0][a]);
            }
            r
        };
        let mut converged = false;
        for _ in 0..50 {
            let r = residual(&x);
            let rn = r.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if !rn.is_finite() {
                return Err(OdeError::NonFinite);
            }
            if rn <= 1e-13 * (1.0 + p.iter().map(|v| v.abs()).fold(0.0, f64::max)) {
                converged = true;
                break;
            }
            let mut jm = Mat::zeros(dim_u, dim_u);
            for col in 0..dim_u {
                let dx = 1e-7 * x[col].abs().max(1.0);
                let mut xp = x.clone();
                xp[col] += dx;
                let rp = residual(&xp);
                let mut xm = x.clone();
                xm[col] -= dx;
                let rm = residual(&xm);
                for row in 0..dim_u {
                    jm.set(row, col, (rp[row] - rm[row]) / (2.0 * dx));
                }
            }
            let lu = dense::lu_factor(&jm).map_err(|_| OdeError::SolveFailed)?;
            let neg: Vec<f64> = r.iter().map(|v| -v).collect();
            let dxv = lu.solve(&neg);
            for (xi, di) in x.iter_mut().zip(&dxv) {
                *xi += di;
            }
        }
        if !converged {
            return Err(OdeError::SolveFailed);
        }
        let mut qn: Vec<Vec<f64>> = vec![q.clone()];
        for j in 1..s {
            qn.push(x[(j - 1) * d..j * d].to_vec());
        }
        let g = gradient(&qn);
        p = g[s - 1].clone();
        q = qn[s - 1].clone();
    }
    Ok((q, p))
}

/// Matrix exponential helper re-exported for the exponential integrators.
pub fn exp_matrix(a: &Mat) -> Result<Mat, OdeError> {
    expm(a).map_err(|_| OdeError::SolveFailed)
}
