//! # Stochastic differential equations
//!
//! Integrators for `dX = a(t, X) dt + g(t, X) dW` with an `n`-vector state
//! and an `m`-dimensional Wiener process; the diffusion callback writes the
//! `n x m` matrix `g` row-major.
//!
//! * [`euler_maruyama`]: strong order 1/2, weak order 1 (Ito);
//! * [`milstein`]: strong order 1 for *commutative* noise (scalar,
//!   diagonal, additive), with the Levy-area-free Milstein term evaluated
//!   by directional finite differences of `g`;
//! * [`srk_platen`]: the derivative-free strong order 1 stochastic
//!   Runge-Kutta scheme of Platen (same commutativity assumption);
//! * [`heun_stratonovich`]: the Stratonovich predictor-corrector scheme;
//! * [`sde_adaptive`]: embedded Euler-Maruyama/Milstein pair with
//!   automatic step size selection; rejected steps are re-drawn from the
//!   same Brownian path through Brownian-bridge splitting of the rejected
//!   increment (rejection sampling with memory), so the step controller
//!   does not bias the distribution of the path.
//!
//! The fixed-step schemes take the Brownian increments explicitly
//! ([`brownian_increments`] generates them) so that strong errors can be
//! measured against exact pathwise solutions.
//!
//! Limitation: the order-1 schemes need commutative noise (otherwise the
//! Levy areas would be required); for non-commutative noise use
//! [`euler_maruyama`].

use crate::kernels::random::Rng;

/// Errors of the SDE integrators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdeError {
    /// Inconsistent dimensions or parameters.
    Invalid,
    /// The state became non-finite.
    NonFinite,
    /// Step size underflow or too many steps in the adaptive solver.
    StepLimit,
}

/// Generates `steps` Brownian increments (`m` components each, variance
/// `h`) from `rng`.
#[must_use]
pub fn brownian_increments(rng: &mut Rng, m: usize, steps: usize, h: f64) -> Vec<Vec<f64>> {
    let sh = h.abs().sqrt();
    (0..steps).map(|_| (0..m).map(|_| sh * rng.normal()).collect()).collect()
}

fn check(y0: &[f64], m: usize, dws: &[Vec<f64>]) -> Result<(), SdeError> {
    if y0.is_empty() || m == 0 || dws.iter().any(|d| d.len() != m) {
        return Err(SdeError::Invalid);
    }
    Ok(())
}

fn finite(y: &[f64]) -> Result<(), SdeError> {
    if y.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(SdeError::NonFinite)
    }
}

/// Euler-Maruyama path with explicit increments; `dws[k]` is the increment
/// over `[t0 + k h, t0 + (k + 1) h]`. Returns the `dws.len() + 1` states.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn euler_maruyama<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        for i in 0..n {
            let mut s = a[i] * h;
            for j in 0..m {
                s += g[i * m + j] * dw[j];
            }
            y[i] += s;
        }
        finite(&y)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// One Milstein correction `(1/2) sum_{j,k} (D_j g)_{.k} (dW_j dW_k - delta_jk h)`
/// for commutative noise, by central differences along the columns of `g`.
fn milstein_term<G: Fn(f64, &[f64], &mut [f64])>(
    diff: &G,
    t: f64,
    y: &[f64],
    g: &[f64],
    m: usize,
    dw: &[f64],
    h: f64,
) -> Vec<f64> {
    let n = y.len();
    let mut corr = vec![0.0; n];
    let (mut gp, mut gm) = (vec![0.0; n * m], vec![0.0; n * m]);
    for j in 0..m {
        let col: Vec<f64> = (0..n).map(|i| g[i * m + j]).collect();
        let cn = col.iter().map(|v| v * v).sum::<f64>().sqrt();
        if cn == 0.0 {
            continue;
        }
        let scale = y.iter().map(|v| v.abs()).fold(1.0, f64::max);
        let eps = 1e-6 * scale / cn;
        let yp: Vec<f64> = (0..n).map(|i| y[i] + eps * col[i]).collect();
        let ym: Vec<f64> = (0..n).map(|i| y[i] - eps * col[i]).collect();
        diff(t, &yp, &mut gp);
        diff(t, &ym, &mut gm);
        for k in 0..m {
            let ito = dw[j] * dw[k] - if j == k { h } else { 0.0 };
            for i in 0..n {
                corr[i] += 0.5 * (gp[i * m + k] - gm[i * m + k]) / (2.0 * eps) * ito;
            }
        }
    }
    corr
}

/// Milstein scheme (strong order 1) for commutative noise; arguments as in
/// [`euler_maruyama`].
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn milstein<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = milstein_term(&diff, t, &y, &g, m, dw, h);
        for i in 0..n {
            let mut s = a[i] * h + corr[i];
            for j in 0..m {
                s += g[i * m + j] * dw[j];
            }
            y[i] += s;
        }
        finite(&y)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Platen's derivative-free stochastic Runge-Kutta scheme.
///
/// Derivative-free stochastic Runge-Kutta scheme of Platen (strong order 1,
/// commutative noise): the Milstein term is replaced by the difference of
/// the diffusion at the supporting values `x + a h + g_j sqrt(h)`.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn srk_platen<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let sh = h.abs().sqrt();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g, mut gs) = (vec![0.0; n], vec![0.0; n * m], vec![0.0; n * m]);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let mut corr = vec![0.0; n];
        for j in 0..m {
            let ups: Vec<f64> = (0..n).map(|i| y[i] + a[i] * h + g[i * m + j] * sh).collect();
            diff(t, &ups, &mut gs);
            for kk in 0..m {
                let ito = dw[j] * dw[kk] - if j == kk { h } else { 0.0 };
                for i in 0..n {
                    corr[i] += (gs[i * m + kk] - g[i * m + kk]) * ito / (2.0 * sh);
                }
            }
        }
        for i in 0..n {
            let mut s = a[i] * h + corr[i];
            for j in 0..m {
                s += g[i * m + j] * dw[j];
            }
            y[i] += s;
        }
        finite(&y)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Heun predictor-corrector for the *Stratonovich* interpretation
/// `dX = a dt + g o dW` (strong order 1 for commutative noise, converges to
/// the Stratonovich solution).
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn heun_stratonovich<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let (mut a2, mut g2) = (vec![0.0; n], vec![0.0; n * m]);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let pred: Vec<f64> = (0..n)
            .map(|i| y[i] + a[i] * h + (0..m).map(|j| g[i * m + j] * dw[j]).sum::<f64>())
            .collect();
        drift(t + h, &pred, &mut a2);
        diff(t + h, &pred, &mut g2);
        for i in 0..n {
            let mut s = f64::midpoint(a[i], a2[i]) * h;
            for j in 0..m {
                s += f64::midpoint(g[i * m + j], g2[i * m + j]) * dw[j];
            }
            y[i] += s;
        }
        finite(&y)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Result of [`sde_adaptive`].
#[derive(Debug, Clone, PartialEq)]
pub struct SdePath {
    /// Accepted time points.
    pub t: Vec<f64>,
    /// States at the accepted time points.
    pub y: Vec<Vec<f64>>,
    /// Number of rejected steps.
    pub rejected: usize,
}

/// Pending Brownian increments (a stack of `(length, increment)` pieces
/// that have been drawn but not yet consumed).
struct BrownianStack {
    pending: Vec<(f64, Vec<f64>)>,
}

impl BrownianStack {
    /// Increment over the next `want` units of time, assembled from pending
    /// pieces (split by Brownian bridge when too long) and fresh draws.
    fn take(&mut self, want: f64, m: usize, rng: &mut Rng) -> Vec<f64> {
        let mut acc = vec![0.0; m];
        let mut got = 0.0;
        while got < want * (1.0 - 1e-12) {
            let need = want - got;
            if let Some((len, inc)) = self.pending.pop() {
                if len <= need * (1.0 + 1e-12) {
                    for j in 0..m {
                        acc[j] += inc[j];
                    }
                    got += len;
                } else {
                    // bridge split: first part of length `need`
                    let frac = need / len;
                    let sd = (need * (len - need) / len).sqrt();
                    let first: Vec<f64> = (0..m).map(|j| frac * inc[j] + sd * rng.normal()).collect();
                    let rest: Vec<f64> = (0..m).map(|j| inc[j] - first[j]).collect();
                    self.pending.push((len - need, rest));
                    for j in 0..m {
                        acc[j] += first[j];
                    }
                    got += need;
                }
            } else {
                let sd = need.sqrt();
                for j in 0..m {
                    acc[j] += sd * rng.normal();
                }
                got += need;
            }
        }
        acc
    }

    /// Returns an unused increment of length `len` to the stack.
    fn give_back(&mut self, len: f64, inc: Vec<f64>) {
        self.pending.push((len, inc));
    }
}

/// Tolerances and limits of [`sde_adaptive_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SdeAdaptiveOptions {
    /// Relative tolerance of the local error.
    pub rtol: f64,
    /// Absolute tolerance of the local error.
    pub atol: f64,
    /// Seed of the Brownian path.
    pub seed: u64,
    /// Maximal number of steps.
    pub max_steps: usize,
}

/// Adaptive-step integration of an Ito SDE with commutative noise on
/// `[t0, t1]`.
///
/// Each step computes the Euler-Maruyama and the Milstein update from the
/// same Brownian increment; the Milstein correction plus the trapezoidal
/// drift correction `h |a(y_EM) - a(y)| / 2` estimate the local error of
/// the former (so additive noise and strongly varying drifts both control
/// the step). A step is accepted when the scaled error
/// `max_i |delta_i| / (atol + rtol |y_i|)` is at most one (the tolerance
/// bounds the local, not the global, pathwise error; the global error
/// accumulates like a random walk). Rejected increments are split by
/// Brownian bridges and re-used, so the realised path is a genuine
/// Brownian sample independent of the step sequence. The returned solution
/// is the Milstein one.
///
/// # Errors
/// [`SdeError::Invalid`], [`SdeError::NonFinite`] or [`SdeError::StepLimit`].
pub fn sde_adaptive_with<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &SdeAdaptiveOptions,
) -> Result<SdePath, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    let SdeAdaptiveOptions { rtol, atol, seed, max_steps } = *opts;
    if y0.is_empty() || m == 0 || t1 <= t0 || rtol <= 0.0 || atol < 0.0 {
        return Err(SdeError::Invalid);
    }
    let n = y0.len();
    let mut rng = Rng::new(seed);
    let mut stack = BrownianStack { pending: Vec::new() };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut path = SdePath { t: vec![t0], y: vec![y0.to_vec()], rejected: 0 };
    let mut h = ((t1 - t0) * 1e-3).max(1e-9);
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut a2 = vec![0.0; n];
    let mut count = 0;
    while t1 - t > 1e-14 * (1.0 + t1.abs()) {
        count += 1;
        if count > max_steps {
            return Err(SdeError::StepLimit);
        }
        h = h.min(t1 - t);
        let dw = stack.take(h, m, &mut rng);
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = milstein_term(&diff, t, &y, &g, m, &dw, h);
        let mut ynew = vec![0.0; n];
        let mut em_pred = vec![0.0; n];
        let mut err = 0.0_f64;
        for i in 0..n {
            let mut em = a[i] * h;
            for j in 0..m {
                em += g[i * m + j] * dw[j];
            }
            ynew[i] = y[i] + em + corr[i];
            em_pred[i] = y[i] + em;
        }
        // drift-curvature estimate: trapezoidal minus Euler drift increment
        drift(t + h, &em_pred, &mut a2);
        for i in 0..n {
            let sc = atol + rtol * y[i].abs().max(ynew[i].abs());
            let drift_err = 0.5 * h * (a2[i] - a[i]).abs();
            err = err.max((corr[i].abs() + drift_err) / sc.max(1e-300));
        }
        if !err.is_finite() || ynew.iter().any(|v| !v.is_finite()) {
            return Err(SdeError::NonFinite);
        }
        if err <= 1.0 {
            t += h;
            y = ynew;
            path.t.push(t);
            path.y.push(y.clone());
            let fac = if err == 0.0 { 2.0 } else { (0.9 / err.sqrt()).clamp(0.5, 2.0) };
            h *= fac;
        } else {
            path.rejected += 1;
            stack.give_back(h, dw);
            let fac = (0.9 / err.sqrt()).clamp(0.1, 0.9);
            h *= fac;
            if h < 1e-14 * (1.0 + t.abs()) {
                return Err(SdeError::StepLimit);
            }
        }
    }
    Ok(path)
}

/// [`sde_adaptive_with`] with the options passed as separate arguments.
///
/// # Errors
/// See [`sde_adaptive_with`].
#[deprecated(note = "use `sde_adaptive_with` with `SdeAdaptiveOptions`")]
#[allow(clippy::too_many_arguments)] // legacy signature kept as a thin wrapper
pub fn sde_adaptive<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    t1: f64,
    y0: &[f64],
    rtol: f64,
    atol: f64,
    seed: u64,
    max_steps: usize,
) -> Result<SdePath, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    sde_adaptive_with(drift, diff, m, t0, t1, y0, &SdeAdaptiveOptions { rtol, atol, seed, max_steps })
}
