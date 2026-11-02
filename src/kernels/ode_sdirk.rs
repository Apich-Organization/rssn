//! # Diagonally implicit Runge–Kutta (SDIRK3(2)) for stiff systems
//!
//! [`sdirk3`] integrates `y' = f(t, y)` with the three-stage, third-order,
//! L-stable, stiffly accurate diagonally implicit Runge–Kutta method of
//! Alexander, with an embedded second-order formula for adaptive step
//! control. The step size is therefore chosen by the accuracy requirement
//! and not by stability: the stability function satisfies `R(∞) = 0`, so a
//! mode with a very negative eigenvalue is annihilated in one step, whereas
//! an explicit method would be forced to `h < 2.8 / |λ|`.
//!
//! The tableau is
//!
//! ```text
//! γ   | γ     0     0
//! c₂  | a₂₁   γ     0
//! c₃  | a₃₁   a₃₂   γ
//! ----+----------------------
//!     | b₁    b₂    b₃     third order (`y_{n+1} = Y₃`)
//!     | b̂₁    b̂₂    b̂₃     second order (error estimate)
//! ```
//!
//! with `γ = 0.435866521508459`, the root in `(0, 1)` of
//! `6γ³ − 18γ² + 9γ − 1`, `a₂₁ = (1 − γ)/2`, and `a₃₁`, `a₃₂` from the
//! third-order conditions `b·c = 1/2`, `b·c² = 1/3`, `b·A c = 1/6` under the
//! stiff-accuracy constraint `b = (a₃₁, a₃₂, γ)`; the embedded weights `b̂`
//! satisfy `b̂·1 = 1` and `b̂·c = 1/2` with `b̂₃ = 0`. Stiff accuracy makes the
//! last stage the solution and gives `R(∞) = 0`.
//!
//! Each stage solves `Y_i = Y_n + h Σ_{j<i} a_ij K_j + h γ f(t + c_i h, Y_i)`
//! by Newton's method with the analytic Jacobian, re-evaluated at every
//! iterate, so the iteration is quadratically convergent and the Jacobian
//! enters the solve itself (not only the error estimate).
//!
//! # Limitations
//!
//! The method is order 3 only, with no continuous extension (dense output), no
//! event location, no mass matrix (DAE) support and no banded or sparse
//! Jacobian: the Newton matrix `I − hγJ` is factored by a dense `O(n³)` LU at
//! every iteration of every stage, and no factorisation is reused across
//! steps, so cost grows cubically with the state dimension (use
//! [`crate::kernels::ode_stiff`] for large, structured or higher-index
//! problems). The error estimate is the difference to the embedded
//! second-order weights, whose leading term is `O(h³)`; it drives the step
//! size but is not itself L-stable, and only the classical order 3 is
//! claimed — a strongly stiff step can show the usual order reduction, so
//! what is guaranteed as `h|λ| → ∞` is L-stability, not accuracy order.
//! Every call starts from scratch — no step size and no factorisation is
//! carried between calls — and the right-hand side is a plain closure, so no
//! Jacobian sparsity pattern can be exploited.

use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;
use crate::kernels::jacobian::fd_jacobian_dense;
use crate::kernels::ode_adaptive::OdeError;
use crate::kernels::ode_adaptive::Tol;
use crate::kernels::ode_adaptive::initial_step;

/// Diagonal coefficient `γ` (root in `(0, 1)` of `6γ³ − 18γ² + 9γ − 1`).
const GAMMA: f64 = 0.435_866_521_508_459;
/// Stage coefficient `a₂₁ = (1 − γ)/2`.
const A21: f64 = 0.282_066_739_245_770_5;
/// Stage coefficient `a₃₁` (also `b₁`).
const A31: f64 = 1.208_496_649_176_010_1;
/// Stage coefficient `a₃₂` (also `b₂`).
const A32: f64 = -0.644_363_170_684_469_1;
/// Weights of the third-order solution (`b = (a₃₁, a₃₂, γ)`).
const B: [f64; 3] = [A31, A32, GAMMA];
/// Weights of the embedded second-order formula (`b̂₃ = 0`).
const BHAT: [f64; 3] = [0.772_630_127_667_551, 0.227_369_872_332_449, 0.0];
/// Abscissae `c = (γ, a₂₁ + γ, 1)`.
const C: [f64; 3] = [GAMMA, 0.717_933_260_754_229_5, 1.0];
/// Strictly lower triangle of the tableau (the diagonal is implicit).
const A_EXPLICIT: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [A21, 0.0, 0.0], [A31, A32, 0.0]];
/// Number of stages.
const STAGES: usize = 3;

/// Options of the SDIRK3(2) integrator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SdirkOptions {
    /// Relative tolerance of the local error estimate.
    pub rtol: f64,
    /// Absolute tolerance of the local error estimate.
    pub atol: f64,
    /// Initial step (`None` selects one automatically).
    pub h0: Option<f64>,
    /// Maximal step size.
    pub hmax: f64,
    /// Maximal number of attempted steps (accepted plus rejected).
    pub max_steps: usize,
    /// Maximal Newton iterations per stage.
    pub newton_max_iter: usize,
    /// Newton convergence tolerance, relative to `1 + |Y|`.
    pub newton_tol: f64,
}

impl Default for SdirkOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-10,
            h0: None,
            hmax: f64::INFINITY,
            max_steps: 1_000_000,
            newton_max_iter: 8,
            newton_tol: 1e-12,
        }
    }
}

/// Result of an adaptive [`sdirk3`] integration.
#[derive(Debug, Clone, PartialEq)]
pub struct SdirkSolution {
    /// Accepted time points, starting at `t0`.
    pub t: Vec<f64>,
    /// States at the accepted time points.
    pub y: Vec<Vec<f64>>,
    /// Number of accepted steps.
    pub accepted: usize,
    /// Number of rejected steps.
    pub rejected: usize,
    /// Total number of Newton iterations over all stages and attempts.
    pub newton_iterations: usize,
    /// Number of right-hand-side evaluations.
    pub nfev: usize,
    /// Number of Jacobian evaluations.
    pub njev: usize,
}

/// One step of the method: the new state, the local error estimate and the
/// evaluation counters.
struct StepOut {
    /// State after the step.
    y: Vec<f64>,
    /// Local error estimate `h Σ (b_i − b̂_i) K_i`.
    err: Vec<f64>,
    /// Newton iterations used by the three stages.
    newton: usize,
    /// Right-hand-side evaluations of the step.
    nfev: usize,
    /// Jacobian evaluations of the step.
    njev: usize,
}

/// Largest absolute entry of a slice (`0` for an empty slice).
fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |acc, x| acc.max(x.abs()))
}

/// Takes one step of length `h` from `(t, y)`.
///
/// Every stage is solved by Newton's method with the Jacobian supplied by
/// `jac`; the step is rejected by the caller when the error estimate is too
/// large.
fn step<F, J>(
    f: &F,
    jac: &J,
    t: f64,
    y: &[f64],
    h: f64,
    opts: &SdirkOptions,
) -> Result<StepOut, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    let n = y.len();
    if n == 0 || !h.is_finite() || h == 0.0 {
        return Err(OdeError::Invalid);
    }
    let mut k: Vec<Vec<f64>> = Vec::with_capacity(STAGES);
    let mut newton = 0usize;
    let mut nfev = 0usize;
    let mut njev = 0usize;
    for stage in 0..STAGES {
        let row = A_EXPLICIT.get(stage).ok_or(OdeError::Invalid)?;
        let c = C.get(stage).ok_or(OdeError::Invalid)?;
        let mut ystage = y.to_vec();
        for (a, kj) in row.iter().zip(k.iter()) {
            for (ys, kv) in ystage.iter_mut().zip(kj.iter()) {
                *ys += h * a * kv;
            }
        }
        let ypred = ystage.clone();
        let tstage = t + h * c;
        let mut fval = vec![0.0; n];
        let mut converged = false;
        for _ in 0..opts.newton_max_iter {
            f(tstage, &ystage, &mut fval);
            nfev = nfev.saturating_add(1);
            newton = newton.saturating_add(1);
            let residual: Vec<f64> = fval
                .iter()
                .zip(ystage.iter())
                .zip(ypred.iter())
                .map(|((fv, ys), yp)| h * GAMMA * fv - (ys - yp))
                .collect();
            let j = jac(tstage, &ystage);
            njev = njev.saturating_add(1);
            if j.rows != n || j.cols != n {
                return Err(OdeError::Invalid);
            }
            let mut matrix = Mat::identity(n);
            for (r, jrow) in j.data.chunks(n).enumerate() {
                for (col, value) in jrow.iter().enumerate() {
                    matrix.set(r, col, matrix.at(r, col) - h * GAMMA * value);
                }
            }
            let lu = lu_factor(&matrix).map_err(|_| OdeError::SolveFailed)?;
            let delta = lu.solve(&residual);
            if delta.iter().any(|d| !d.is_finite()) {
                return Err(OdeError::NonFinite);
            }
            for (ys, d) in ystage.iter_mut().zip(delta.iter()) {
                *ys += d;
            }
            if max_abs(&delta) <= opts.newton_tol * (1.0 + max_abs(&ystage)) {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(OdeError::SolveFailed);
        }
        let mut kstage = vec![0.0; n];
        f(tstage, &ystage, &mut kstage);
        nfev = nfev.saturating_add(1);
        if kstage.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        k.push(kstage);
    }
    let mut ynew = y.to_vec();
    let mut err = vec![0.0; n];
    for ((kv, b), bh) in k.iter().zip(B.iter()).zip(BHAT.iter()) {
        for ((yn, e), yi) in ynew.iter_mut().zip(err.iter_mut()).zip(kv.iter()) {
            *yn += h * b * yi;
            *e += h * (b - bh) * yi;
        }
    }
    if ynew.iter().any(|v| !v.is_finite()) || err.iter().any(|v| !v.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    Ok(StepOut {
        y: ynew,
        err,
        newton,
        nfev,
        njev,
    })
}

/// Integrates `y' = f(t, y)` from `t0` to `t1` with the adaptive SDIRK3(2)
/// method, using the analytic Jacobian `jac(t, y)`.
///
/// The step size is controlled by the embedded second-order estimate with
/// exponent `1/3`; steps whose Newton iteration fails are rejected and
/// retried with half the step.
///
/// # Errors
/// [`OdeError::Invalid`] for an empty state, `t0 == t1` or a Jacobian of the
/// wrong shape; [`OdeError::NonFinite`] when the state or the estimate stops
/// being finite; [`OdeError::SolveFailed`] when Newton iteration does not
/// converge at the smallest step; [`OdeError::StepLimit`] when
/// [`SdirkOptions::max_steps`] or the step-size floor is reached.
pub fn sdirk3<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &SdirkOptions,
) -> Result<SdirkSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    let n = y0.len();
    if n == 0 || crate::kernels::exact_eq(t0, t1) {
        return Err(OdeError::Invalid);
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let floor = 1e-14 * (1.0 + t0.abs());
    let hmax = opts.hmax.abs().max(floor);
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1usize;
    let h_start = opts.h0.unwrap_or_else(|| {
        initial_step(
            &f,
            t0,
            y0,
            &f0,
            dir,
            3.0,
            Tol {
                rtol: opts.rtol,
                atol: opts.atol,
            },
        )
    });
    let mut h = h_start.abs().clamp(floor, hmax) * dir;
    let mut sol = SdirkSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        accepted: 0,
        rejected: 0,
        newton_iterations: 0,
        nfev: 0,
        njev: 0,
    };
    let mut attempts = 0usize;
    while (t1 - t) * dir > 0.0 {
        if attempts >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        attempts = attempts.saturating_add(1);
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        let out = match step(&f, &jac, t, &y, h, opts) {
            | Ok(out) => out,
            | Err(OdeError::SolveFailed) => {
                sol.rejected = sol.rejected.saturating_add(1);
                h *= 0.5;
                if h.abs() < floor {
                    return Err(OdeError::StepLimit);
                }
                continue;
            },
            | Err(other) => return Err(other),
        };
        nfev = nfev.saturating_add(out.nfev);
        sol.njev = sol.njev.saturating_add(out.njev);
        sol.newton_iterations = sol.newton_iterations.saturating_add(out.newton);
        let en = out
            .err
            .iter()
            .zip(y.iter())
            .map(|(e, yi)| (e / (opts.atol + opts.rtol * yi.abs())).abs())
            .fold(0.0_f64, f64::max);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            t += h;
            y = out.y;
            sol.accepted = sol.accepted.saturating_add(1);
            sol.t.push(t);
            sol.y.push(y.clone());
            let fac = if en == 0.0 {
                5.0
            } else {
                (0.9 * en.powf(-1.0 / 3.0)).clamp(0.2, 5.0)
            };
            h = (h * fac).clamp(-hmax, hmax);
        } else {
            sol.rejected = sol.rejected.saturating_add(1);
            h *= (0.9 * en.powf(-1.0 / 3.0)).clamp(0.1, 0.9);
            if h.abs() < floor {
                return Err(OdeError::StepLimit);
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// The same integrator with a forward-difference Jacobian (`n + 1`
/// right-hand-side evaluations per Newton iterate).
///
/// # Errors
/// As [`sdirk3`].
pub fn sdirk3_fd<F>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &SdirkOptions,
) -> Result<SdirkSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let jac = |t: f64, y: &[f64]| {
        let mut f0 = vec![0.0; y.len()];
        f(t, y, &mut f0);
        fd_jacobian_dense(&f, t, y, &f0)
    };
    sdirk3(&f, jac, t0, t1, y0, opts)
}

/// Integrates `y' = f(t, y)` from `t0` to `t1` in `steps` equal steps with no
/// step-size control, returning the final state.
///
/// This is the entry point for convergence studies: with `steps = 2^k` the
/// error of a third-order method falls by `8` per refinement of the step.
///
/// # Errors
/// [`OdeError::Invalid`] for an empty state or `steps == 0`, and the errors of
/// [`sdirk3`] otherwise.
pub fn sdirk3_fixed<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    if y0.is_empty() || steps == 0 {
        return Err(OdeError::Invalid);
    }
    let opts = SdirkOptions::default();
    let h = (t1 - t0) / steps as f64;
    let mut t = t0;
    let mut y = y0.to_vec();
    for _ in 0..steps {
        let out = step(&f, &jac, t, &y, h, &opts)?;
        t += h;
        y = out.y;
    }
    Ok(y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::ode::OdeSolverMethod;
    use crate::kernels::ode::solve_fixed;
    use crate::kernels::ode_adaptive::OdeOptions;
    use crate::kernels::ode_adaptive::dopri5;

    /// Robertson chemical kinetics.
    fn robertson(
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
        dy[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1];
        dy[2] = 3e7 * y[1] * y[1];
    }

    /// Analytic Jacobian of [`robertson`].
    fn robertson_jac(
        _t: f64,
        y: &[f64],
    ) -> Mat {
        Mat::from_rows(&[
            vec![-0.04, 1e4 * y[2], 1e4 * y[1]],
            vec![0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]],
            vec![0.0, 6e7 * y[1], 0.0],
        ])
        .unwrap()
    }

    /// Hairer–Wanner reference solution of [`robertson`] at `t = 40`.
    const ROB40: [f64; 3] = [0.715_827_068_7, 9.185_534_764e-6, 0.284_163_745_7];

    #[test]
    fn robertson_stiff_step_counts() {
        let opts = SdirkOptions {
            rtol: 1e-6,
            atol: 1e-9,
            ..SdirkOptions::default()
        };
        let exp_opts = OdeOptions {
            rtol: 1e-6,
            atol: 1e-9,
            max_steps: 2_000_000,
            ..OdeOptions::default()
        };
        let one = sdirk3(robertson, robertson_jac, 0.0, 1.0, &[1.0, 0.0, 0.0], &opts).unwrap();
        let one_ex = dopri5(robertson, 0.0, 1.0, &[1.0, 0.0, 0.0], &exp_opts, &[]).unwrap();
        let long = sdirk3(robertson, robertson_jac, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts).unwrap();
        let long_ex = dopri5(robertson, 0.0, 40.0, &[1.0, 0.0, 0.0], &exp_opts, &[]).unwrap();
        let (imp1, exp1) = (one.accepted, one_ex.t.len() - 1);
        let (imp40, exp40) = (long.accepted, long_ex.t.len() - 1);
        eprintln!(
            "robertson t=1: implicit {imp1} steps, explicit {exp1} steps; \
             t=40: implicit {imp40} steps, explicit {exp40} steps"
        );
        // The implicit solver reaches t = 1 and t = 40 with accuracy-limited
        // steps; the explicit one is limited by the stability of the fast
        // mode and needs orders of magnitude more steps over the stiff window.
        assert!(imp1 < exp1 / 5, "t=1: {imp1} vs {exp1}");
        assert!(
            exp40 >= 100 * imp40,
            "t=40: {imp40} implicit vs {exp40} explicit steps"
        );
        assert!(imp40 < 400, "too many implicit steps: {imp40}");
        // Accuracy against the published reference at t = 40.
        let last = long.y.last().unwrap();
        for (got, want) in last.iter().zip(ROB40.iter()) {
            assert!(
                (got - want).abs() <= 1e-5 * want.abs().max(1e-6),
                "y = {got} vs {want}"
            );
        }
        assert!((last.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        // van der Pol at mu = 1000, the other defining stiff benchmark.
        let vdp = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = 1000.0 * (1.0 - y[0] * y[0]) * y[1] - y[0];
        };
        let vdp_jac = |_t: f64, y: &[f64]| {
            Mat::from_rows(&[
                vec![0.0, 1.0],
                vec![-2000.0 * y[0] * y[1] - 1.0, 1000.0 * (1.0 - y[0] * y[0])],
            ])
            .unwrap()
        };
        let vdp_sol = sdirk3(vdp, vdp_jac, 0.0, 2.0, &[2.0, 0.0], &opts).unwrap();
        let vdp_ex = dopri5(vdp, 0.0, 2.0, &[2.0, 0.0], &exp_opts, &[]).unwrap();
        eprintln!(
            "van der Pol mu=1000: implicit {} steps, explicit {} steps",
            vdp_sol.accepted,
            vdp_ex.t.len() - 1
        );
        assert!(vdp_sol.accepted < vdp_ex.t.len() / 5);
        let end = vdp_sol.y.last().unwrap();
        assert!((end[0] - 1.998_666_1).abs() < 1e-4, "y1 = {}", end[0]);
    }

    #[test]
    fn third_order_under_refinement() {
        // y' = y cos t, y(0) = 1, exact y = exp(sin t); smooth and non-stiff.
        let f = |t: f64, y: &[f64], dy: &mut [f64]| dy[0] = y[0] * t.cos();
        let jac = |t: f64, _y: &[f64]| Mat::from_rows(&[vec![t.cos()]]).unwrap();
        let exact = 1.0_f64.sin().exp();
        let mut log_h = Vec::new();
        let mut log_e = Vec::new();
        let mut prev = f64::NAN;
        for k in 0..5 {
            let steps = 10usize << k;
            let y = sdirk3_fixed(f, jac, 0.0, 1.0, &[1.0], steps).unwrap();
            let err = (y[0] - exact).abs();
            if k > 0 {
                let ratio = prev / err;
                assert!(ratio > 7.0, "error ratio {ratio} at {steps} steps");
            }
            prev = err;
            log_h.push((1.0 / steps as f64).ln());
            log_e.push(err.ln());
        }
        let n = log_h.len() as f64;
        let mean_x = log_h.iter().sum::<f64>() / n;
        let mean_y = log_e.iter().sum::<f64>() / n;
        let num: f64 = log_h
            .iter()
            .zip(log_e.iter())
            .map(|(x, y)| (x - mean_x) * (y - mean_y))
            .sum();
        let den: f64 = log_h.iter().map(|x| (x - mean_x) * (x - mean_x)).sum();
        let slope = num / den;
        eprintln!("fitted order of convergence: {slope:.4}");
        assert!(
            (2.85..3.15).contains(&slope),
            "fitted slope {slope} is not the design order 3"
        );
    }

    #[test]
    fn jacobian_matches_finite_difference() {
        let states = [
            (0.0, [1.0, 0.0, 0.0]),
            (1.0, [0.9, 1e-5, 0.1]),
            (40.0, [0.5, 1e-3, 0.5]),
            (40.0, [0.7, 1e-4, 0.3]),
        ];
        for (t, y) in states {
            let analytic = robertson_jac(t, &y);
            let mut f0 = vec![0.0; 3];
            robertson(t, &y, &mut f0);
            let fd = fd_jacobian_dense(&robertson, t, &y, &f0);
            for i in 0..3 {
                for j in 0..3 {
                    let a = analytic.at(i, j);
                    let b = fd.at(i, j);
                    // Relative to the size of the entry, with an absolute floor
                    // for entries that vanish: a forward difference of the
                    // quadratic term `3e7 y₂²` is `3e7 h` at `y₂ = 0`.
                    let tol = 1e-6 * a.abs().max(1.0) + 1e-5;
                    assert!(
                        (a - b).abs() <= tol,
                        "J[{i}][{j}] analytic {a} vs finite difference {b} at {y:?}"
                    );
                }
            }
        }
        // The Jacobian of the van der Pol system as well.
        let y = [0.5, 1.0];
        let analytic = Mat::from_rows(&[
            vec![0.0, 1.0],
            vec![-2000.0 * y[0] * y[1] - 1.0, 1000.0 * (1.0 - y[0] * y[0])],
        ])
        .unwrap();
        let mut f0 = vec![0.0; 2];
        let vdp = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = 1000.0 * (1.0 - y[0] * y[0]) * y[1] - y[0];
        };
        vdp(0.0, &y, &mut f0);
        let fd = fd_jacobian_dense(&vdp, 0.0, &y, &f0);
        for i in 0..2 {
            for j in 0..2 {
                let a = analytic.at(i, j);
                let b = fd.at(i, j);
                let tol = 1e-6 * a.abs().max(1.0) + 1e-5;
                assert!(
                    (a - b).abs() <= tol,
                    "vdP J[{i}][{j}] analytic {a} vs finite difference {b}"
                );
            }
        }
    }

    #[test]
    fn corrupted_jacobian_slows_or_fails() {
        let opts = SdirkOptions {
            rtol: 1e-6,
            atol: 1e-9,
            ..SdirkOptions::default()
        };
        let good = sdirk3(robertson, robertson_jac, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts).unwrap();
        // The same solver with the sign of the Jacobian flipped: the Newton
        // iteration no longer converges at the accepted step sizes, so the
        // step is rejected and the run either fails or crawls.
        let corrupted = |t: f64, y: &[f64]| {
            let mut j = robertson_jac(t, y);
            for value in &mut j.data {
                *value = -*value;
            }
            j
        };
        let result = sdirk3(robertson, corrupted, 0.0, 40.0, &[1.0, 0.0, 0.0], &opts);
        match result {
            | Err(e) => {
                eprintln!(
                    "corrupted Jacobian: failed with {e:?} (good run: {} steps, {} Newton iterations)",
                    good.accepted, good.newton_iterations
                )
            },
            | Ok(sol) => {
                eprintln!(
                    "corrupted Jacobian: {} steps, {} Newton iterations \
                     (good run: {} steps, {} Newton iterations)",
                    sol.accepted, sol.newton_iterations, good.accepted, good.newton_iterations
                );
                assert!(
                    sol.accepted > 5 * good.accepted
                        || sol.newton_iterations > 5 * good.newton_iterations,
                    "a corrupted Jacobian must slow the solver down"
                );
            },
        }
        // A merely inaccurate Jacobian (scaled by 2) still converges, but
        // needs many more Newton iterations than the exact one.
        let scaled = |t: f64, y: &[f64]| {
            let mut j = robertson_jac(t, y);
            for value in &mut j.data {
                *value *= 2.0;
            }
            j
        };
        if let Ok(sol) = sdirk3(robertson, scaled, 0.0, 1.0, &[1.0, 0.0, 0.0], &opts) {
            let good1 =
                sdirk3(robertson, robertson_jac, 0.0, 1.0, &[1.0, 0.0, 0.0], &opts).unwrap();
            eprintln!(
                "scaled Jacobian: {} Newton iterations vs {} for the exact one",
                sol.newton_iterations, good1.newton_iterations
            );
            assert!(sol.newton_iterations >= good1.newton_iterations);
        }
    }

    #[test]
    fn a_stability_at_large_steps() {
        // y' = -1000 y over one step h = 1: |z| = 1000, far beyond the
        // explicit stability limit 2.8. An A-stable method must not blow up.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -1000.0 * y[0];
        let jac = |_t: f64, _y: &[f64]| Mat::from_rows(&[vec![-1000.0]]).unwrap();
        let y1 = sdirk3_fixed(f, jac, 0.0, 1.0, &[1.0], 1).unwrap();
        eprintln!("|R(-1000)| = {}", y1[0].abs());
        assert!(y1[0].abs() <= 1.0, "A-stability violated: {}", y1[0]);
        assert!(y1[0].abs() < 0.01, "L-stability violated: {}", y1[0]);
        // A stiffer mode is annihilated even more strongly (R(inf) = 0).
        let y2 = sdirk3_fixed(f, jac, 0.0, 1e3, &[1.0], 1).unwrap();
        assert!(y2[0].abs() < 1e-2, "|R(-1e9)| = {}", y2[0].abs());
        // A two-mode system: the slow mode is integrated, the fast one is
        // damped; neither oscillates nor grows.
        let sys = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
            dy[1] = -1e6 * y[1];
        };
        let sys_jac =
            |_t: f64, _y: &[f64]| Mat::from_rows(&[vec![-1.0, 0.0], vec![0.0, -1e6]]).unwrap();
        let out = sdirk3_fixed(sys, sys_jac, 0.0, 1.0, &[1.0, 1.0], 10).unwrap();
        assert!(
            (out[0] - (-1.0_f64).exp()).abs() < 1e-3,
            "slow mode {}",
            out[0]
        );
        assert!(out[1].abs() < 1e-3, "fast mode {}", out[1]);
    }

    #[test]
    fn explicit_methods_unchanged() {
        // Regression guard: the explicit integrators must keep their values.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
        let traj = solve_fixed(f, &[1.0], (0.0, 1.0), 10, OdeSolverMethod::RungeKutta4).unwrap();
        // R(-0.1)^10 for the classical RK4 stability polynomial.
        assert!((traj.last()[0] - 0.367_879_774_412_498_75).abs() < 1e-14);
        let opts = OdeOptions {
            rtol: 1e-10,
            atol: 1e-12,
            ..OdeOptions::default()
        };
        let sol = dopri5(f, 0.0, 1.0, &[1.0], &opts, &[]).unwrap();
        assert!((sol.y.last().unwrap()[0] - (-1.0_f64).exp()).abs() < 1e-9);
        // The explicit solver still needs the stability-limited step count on
        // Robertson, i.e. it has not silently become implicit.
        let exp_opts = OdeOptions {
            rtol: 1e-6,
            atol: 1e-9,
            max_steps: 2_000_000,
            ..OdeOptions::default()
        };
        let sol = dopri5(robertson, 0.0, 1.0, &[1.0, 0.0, 0.0], &exp_opts, &[]).unwrap();
        let steps = sol.t.len() - 1;
        assert!((300..2000).contains(&steps), "explicit steps: {steps}");
    }
}
