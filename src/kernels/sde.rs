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
//! * [`milstein_noncommutative`]: strong order 1 for **non-commutative**
//!   noise, with the iterated Ito integrals `I_jl` (Levy areas) simulated
//!   conditionally on the increments by [`levy_areas`] (Fourier expansion of
//!   the Brownian bridge with the Gaussian tail correction of Wiktorsson);
//!   [`sample_noise`] draws increments and areas, [`aggregate_noise`] merges
//!   steps exactly so that coarse and fine schemes see the same path.
//!
//! * [`sde_adaptive_noncommutative_with`]: the adaptive Milstein solver for
//!   **non-commutative** noise. The Brownian path lives on a fixed fine grid of
//!   cells, each with its Levy areas; a step is a union of cells combined by
//!   Chen's relation, so rejecting or splitting a step keeps the path (and its
//!   areas) consistent. [`noncommutative_noise_grid`] returns the same path for
//!   a fixed-step reference.
//!
//! * [`sde_ensemble`]: many independent sample paths of a fixed-step scheme, run on the rayon
//!   pool; member `i` is seeded by [`ensemble_seed`]`(seed, i)`, so the ensemble does not depend on
//!   the thread count.
//!
//! Limitation: [`milstein`], [`srk_platen`] and [`sde_adaptive_with`] need
//! commutative noise; use [`milstein_noncommutative`] or
//! [`sde_adaptive_noncommutative_with`] otherwise. The non-commutative adaptive
//! solver cannot take a step shorter than its grid cell (it reports
//! [`SdeError::StepLimit`]; use more cells) and, being a pathwise Milstein
//! scheme with the Levy areas of [`levy_areas`], inherits their finite-mode
//! accuracy.

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
pub fn brownian_increments(
    rng: &mut Rng,
    m: usize,
    steps: usize,
    h: f64,
) -> Vec<Vec<f64>> {
    let sh = h.abs().sqrt();
    (0..steps)
        .map(|_| (0..m).map(|_| sh * rng.normal()).collect())
        .collect()
}

fn check(
    y0: &[f64],
    m: usize,
    dws: &[Vec<f64>],
) -> Result<(), SdeError> {
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

/// Milstein correction `sum_{j,k} (D_j g)_{.k} I_jk` for the iterated Ito
/// integrals `areas[j * m + k] = I_jk` (row-major), with the directional
/// derivatives `D_j g` taken by central differences along the columns of `g`.
fn milstein_term_areas<G: Fn(f64, &[f64], &mut [f64])>(
    diff: &G,
    t: f64,
    y: &[f64],
    g: &[f64],
    m: usize,
    areas: &[f64],
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
            for i in 0..n {
                corr[i] += (gp[i * m + k] - gm[i * m + k]) / (2.0 * eps) * areas[j * m + k];
            }
        }
    }
    corr
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
    let mut areas = vec![0.0; m * m];
    for j in 0..m {
        for k in 0..m {
            areas[j * m + k] = 0.5 * (dw[j] * dw[k] - if j == k { h } else { 0.0 });
        }
    }
    milstein_term_areas(diff, t, y, g, m, &areas)
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

/// Brownian increments together with the iterated Ito integrals over each
/// step, as needed by the Milstein scheme for non-commutative noise.
#[derive(Debug, Clone, PartialEq)]
pub struct NoiseSample {
    /// `dws[k][j]`: increment of the `j`-th Wiener process over step `k`.
    pub dws: Vec<Vec<f64>>,
    /// `areas[k][j * m + l] = I_jl = int int dW_j(u) dW_l(s)` over step `k`
    /// (inner integral `j`), so `I_jl + I_lj = dW_j dW_l` and
    /// `I_jj = (dW_j^2 - h) / 2`.
    pub areas: Vec<Vec<f64>>,
}

/// `sum_{r>p} 1/r^2`.
fn zeta2_tail(p: usize) -> f64 {
    if p >= 8 {
        let x = p as f64;
        1.0 / x - 0.5 / (x * x) + 1.0 / (6.0 * x.powi(3)) - 1.0 / (30.0 * x.powi(5))
    } else {
        let s: f64 = (1..=p).map(|r| 1.0 / (r as f64 * r as f64)).sum();
        std::f64::consts::PI.powi(2) / 6.0 - s
    }
}

/// Samples the iterated Ito integrals `I_jl` of one step of length `h`.
///
/// The integrals are drawn conditionally on the increments `dw` by the Fourier expansion of the
/// Brownian bridge (Kloeden, Platen and Wright) with the Gaussian
/// approximation of the truncated tail that reproduces its full conditional
/// covariance (Wiktorsson's improvement): with `xi = dw / sqrt(h)`,
///
/// `I_jl / h = xi_j xi_l / 2 + A_jl + sqrt(rho) (mu_j xi_l - mu_l xi_j) + sqrt(rho) nu_jl`
///
/// for `j < l`, where `A` is the Fourier sum over `terms` modes, `rho` the
/// tail sum divided by `2 pi^2`, and `mu`, `nu` standard normal. The
/// remainder error decays like `1 / terms^{3/2}` in distribution; a few
/// terms suffice for strong order one.
#[must_use]
pub fn levy_areas(
    rng: &mut Rng,
    m: usize,
    h: f64,
    dw: &[f64],
    terms: usize,
) -> Vec<f64> {
    let p = terms.max(1);
    let xi: Vec<f64> = dw.iter().map(|d| d / h.sqrt()).collect();
    let tail = zeta2_tail(p);
    let rho = tail / (2.0 * std::f64::consts::PI.powi(2));
    let srho = rho.sqrt();
    // Fourier coefficients
    let zeta: Vec<Vec<f64>> = (0..m)
        .map(|_| (0..p).map(|_| rng.normal()).collect())
        .collect();
    let eta: Vec<Vec<f64>> = (0..m)
        .map(|_| (0..p).map(|_| rng.normal()).collect())
        .collect();
    let mu: Vec<f64> = (0..m).map(|_| rng.normal()).collect();
    let mut areas = vec![0.0; m * m];
    for j in 0..m {
        areas[j * m + j] = 0.5 * (dw[j] * dw[j] - h);
        for l in j + 1..m {
            let mut a = 0.0;
            for r in 0..p {
                let rr = (r + 1) as f64;
                a += (zeta[j][r] * (std::f64::consts::SQRT_2 * xi[l] + eta[l][r])
                    - zeta[l][r] * (std::f64::consts::SQRT_2 * xi[j] + eta[j][r]))
                    / rr;
            }
            a /= 2.0 * std::f64::consts::PI;
            let nu = rng.normal();
            let ijl =
                h * (0.5 * xi[j] * xi[l] + a + srho * (mu[j] * xi[l] - mu[l] * xi[j]) + srho * nu);
            areas[j * m + l] = ijl;
            areas[l * m + j] = dw[j] * dw[l] - ijl;
        }
    }
    areas
}

/// Samples `steps` Brownian increments and their iterated integrals
/// ([`levy_areas`]) for `m` Wiener processes and step size `h`.
#[must_use]
pub fn sample_noise(
    rng: &mut Rng,
    m: usize,
    steps: usize,
    h: f64,
    terms: usize,
) -> NoiseSample {
    let dws = brownian_increments(rng, m, steps, h);
    let areas = dws
        .iter()
        .map(|dw| levy_areas(rng, m, h, dw, terms))
        .collect();
    NoiseSample { dws, areas }
}

/// Merges groups of `group` consecutive steps.
///
/// The increments add up and the iterated integrals combine exactly as
/// `I_jl = sum_k I_jl^k + sum_{k<k'} dW_j^k dW_l^{k'}`. Used to compare a
/// coarse scheme with a fine reference on the same Brownian path.
#[must_use]
pub fn aggregate_noise(
    sample: &NoiseSample,
    group: usize,
    m: usize,
) -> NoiseSample {
    let g = group.max(1);
    let mut dws = Vec::new();
    let mut areas = Vec::new();
    for (dc, ac) in sample.dws.chunks(g).zip(sample.areas.chunks(g)) {
        let mut dw = vec![0.0; m];
        let mut ar = vec![0.0; m * m];
        for (d, a) in dc.iter().zip(ac) {
            for j in 0..m {
                for l in 0..m {
                    ar[j * m + l] += a[j * m + l] + dw[j] * d[l];
                }
            }
            for j in 0..m {
                dw[j] += d[j];
            }
        }
        dws.push(dw);
        areas.push(ar);
    }
    NoiseSample { dws, areas }
}

/// Milstein scheme (strong order 1) for **non-commutative** noise.
///
/// The correction uses the sampled iterated integrals of `noise` instead of the
/// symmetric products `(dW_j dW_k - delta_jk h) / 2`. For commutative noise it
/// is equivalent to [`milstein`] in law.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn milstein_noncommutative<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    noise: &NoiseSample,
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, &noise.dws)?;
    if noise.areas.len() != noise.dws.len() || noise.areas.iter().any(|a| a.len() != m * m) {
        return Err(SdeError::Invalid);
    }
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    for (k, (dw, ar)) in noise.dws.iter().zip(&noise.areas).enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = milstein_term_areas(&diff, t, &y, &g, m, ar);
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
            let ups: Vec<f64> = (0..n)
                .map(|i| y[i] + a[i] * h + g[i * m + j] * sh)
                .collect();
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
    fn take(
        &mut self,
        want: f64,
        m: usize,
        rng: &mut Rng,
    ) -> Vec<f64> {
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
                    let first: Vec<f64> =
                        (0..m).map(|j| frac * inc[j] + sd * rng.normal()).collect();
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
    fn give_back(
        &mut self,
        len: f64,
        inc: Vec<f64>,
    ) {
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
    let SdeAdaptiveOptions {
        rtol,
        atol,
        seed,
        max_steps,
    } = *opts;
    if y0.is_empty() || m == 0 || t1 <= t0 || rtol <= 0.0 || atol < 0.0 {
        return Err(SdeError::Invalid);
    }
    let n = y0.len();
    let mut rng = Rng::new(seed);
    let mut stack = BrownianStack { pending: Vec::new() };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut path = SdePath {
        t: vec![t0],
        y: vec![y0.to_vec()],
        rejected: 0,
    };
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
            let fac = if err == 0.0 {
                2.0
            } else {
                (0.9 / err.sqrt()).clamp(0.5, 2.0)
            };
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

/// Options of [`sde_adaptive_noncommutative_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SdeNoncommutativeOptions {
    /// Tolerances, seed and step limit (as for [`sde_adaptive_with`]).
    pub base: SdeAdaptiveOptions,
    /// Number of cells of the Brownian grid on `[t0, t1]`. Steps are unions of
    /// whole cells, so this is the finest step `(t1 - t0) / cells`.
    pub cells: usize,
    /// Fourier modes of the Levy-area simulation ([`levy_areas`]).
    pub terms: usize,
}

/// The Brownian path of [`sde_adaptive_noncommutative_with`]: increments and
/// iterated integrals of the fine cells, drawn lazily in cell order, so cell `k`
/// depends only on the seed and `k` and never on the step sequence.
struct NoiseGrid {
    rng: Rng,
    m: usize,
    hc: f64,
    terms: usize,
    cells: NoiseSample,
}

impl NoiseGrid {
    #[allow(clippy::missing_const_for_fn)] // `Rng::new` is not `const`
    fn new(
        seed: u64,
        m: usize,
        hc: f64,
        terms: usize,
    ) -> Self {
        Self {
            rng: Rng::new(seed),
            m,
            hc,
            terms,
            cells: NoiseSample {
                dws: Vec::new(),
                areas: Vec::new(),
            },
        }
    }

    /// Draws cells until `upto` of them exist.
    fn ensure(
        &mut self,
        upto: usize,
    ) {
        while self.cells.dws.len() < upto {
            let dw = brownian_increments(&mut self.rng, self.m, 1, self.hc)
                .pop()
                .unwrap_or_default();
            let areas = levy_areas(&mut self.rng, self.m, self.hc, &dw, self.terms);
            self.cells.dws.push(dw);
            self.cells.areas.push(areas);
        }
    }

    /// Increment and iterated integrals over the `count` cells from `first`,
    /// combined exactly by Chen's relation
    /// `I_jl = I_jl^1 + I_jl^2 + dW_j^1 dW_l^2` (as [`aggregate_noise`]).
    fn span(
        &mut self,
        first: usize,
        count: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        self.ensure(first + count);
        let m = self.m;
        let mut dw = vec![0.0; m];
        let mut ar = vec![0.0; m * m];
        for k in first..first + count {
            let (d, a) = (&self.cells.dws[k], &self.cells.areas[k]);
            for j in 0..m {
                for l in 0..m {
                    ar[j * m + l] += a[j * m + l] + dw[j] * d[l];
                }
            }
            for j in 0..m {
                dw[j] += d[j];
            }
        }
        (dw, ar)
    }
}

/// The fine-cell Brownian path of [`sde_adaptive_noncommutative_with`] as a [`NoiseSample`].
///
/// Step `h = (t1 - t0) / cells` on `[t0, t1]`; use it to compare the adaptive solution with a
/// fine fixed-step [`milstein_noncommutative`] reference on the *same* path.
#[must_use]
pub fn noncommutative_noise_grid(
    m: usize,
    t0: f64,
    t1: f64,
    opts: &SdeNoncommutativeOptions,
) -> NoiseSample {
    let hc = (t1 - t0) / opts.cells.max(1) as f64;
    let mut grid = NoiseGrid::new(opts.base.seed, m, hc, opts.terms);
    grid.ensure(opts.cells);
    grid.cells
}

/// Adaptive-step Milstein integration of an Ito SDE with **non-commutative**
/// noise on `[t0, t1]`.
///
/// The Brownian path is fixed in advance on a grid of `opts.cells` cells, each
/// with its increment and its Levy areas ([`levy_areas`]); a step is a union of
/// cells, whose increment and areas are combined *exactly* by Chen's relation
/// `A_{0,2h} = A_1 + A_2 + (W_1 x W_2 - W_2 x W_1) / 2`. Splitting a rejected
/// step is therefore merging fewer of the same cells: the realised path is the
/// same for every step sequence, the errors of the controller cannot bias the
/// noise, and a fixed-step reference on the cells sees exactly the path the
/// adaptive solver used ([`noncommutative_noise_grid`]).
///
/// The error control is that of [`sde_adaptive_with`], with the Milstein
/// correction built from the areas (for commutative noise the correction is the
/// same function of the increments, so the estimate and the controller agree
/// with it in law). A step of a single cell that is still rejected ends the run
/// with [`SdeError::StepLimit`]: use more cells.
///
/// # Errors
/// [`SdeError::Invalid`], [`SdeError::NonFinite`] or [`SdeError::StepLimit`].
pub fn sde_adaptive_noncommutative_with<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &SdeNoncommutativeOptions,
) -> Result<SdePath, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    let SdeAdaptiveOptions {
        rtol,
        atol,
        seed,
        max_steps,
    } = opts.base;
    if y0.is_empty() || m == 0 || t1 <= t0 || rtol <= 0.0 || atol < 0.0 || opts.cells == 0 {
        return Err(SdeError::Invalid);
    }
    let n = y0.len();
    let hc = (t1 - t0) / opts.cells as f64;
    let mut grid = NoiseGrid::new(seed, m, hc, opts.terms);
    let mut cell = 0_usize;
    let mut y = y0.to_vec();
    let mut path = SdePath {
        t: vec![t0],
        y: vec![y0.to_vec()],
        rejected: 0,
    };
    let mut size = (opts.cells / 1000).max(1);
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut a2 = vec![0.0; n];
    let mut count = 0;
    while cell < opts.cells {
        count += 1;
        if count > max_steps {
            return Err(SdeError::StepLimit);
        }
        let k = size.clamp(1, opts.cells - cell);
        let h = k as f64 * hc;
        let t = t0 + cell as f64 * hc;
        let (dw, ar) = grid.span(cell, k);
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = milstein_term_areas(&diff, t, &y, &g, m, &ar);
        let mut ynew = vec![0.0; n];
        let mut em_pred = vec![0.0; n];
        for i in 0..n {
            let mut em = a[i] * h;
            for j in 0..m {
                em += g[i * m + j] * dw[j];
            }
            ynew[i] = y[i] + em + corr[i];
            em_pred[i] = y[i] + em;
        }
        drift(t + h, &em_pred, &mut a2);
        let mut err = 0.0_f64;
        for i in 0..n {
            let sc = atol + rtol * y[i].abs().max(ynew[i].abs());
            let drift_err = 0.5 * h * (a2[i] - a[i]).abs();
            err = err.max((corr[i].abs() + drift_err) / sc.max(1e-300));
        }
        if !err.is_finite() || ynew.iter().any(|v| !v.is_finite()) {
            return Err(SdeError::NonFinite);
        }
        if err <= 1.0 {
            cell += k;
            y = ynew;
            path.t.push(t0 + cell as f64 * hc);
            path.y.push(y.clone());
            let fac = if err == 0.0 {
                2.0
            } else {
                (0.9 / err.sqrt()).clamp(0.5, 2.0)
            };
            size = ((k as f64 * fac).round() as usize).max(1);
        } else {
            if k == 1 {
                return Err(SdeError::StepLimit);
            }
            path.rejected += 1;
            let fac = (0.9 / err.sqrt()).clamp(0.1, 0.9);
            size = ((k as f64 * fac).floor() as usize).clamp(1, k - 1);
        }
    }
    Ok(path)
}

/// The fixed-step scheme used by [`sde_ensemble`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnsembleScheme {
    /// [`euler_maruyama`].
    EulerMaruyama,
    /// [`milstein`] (commutative noise).
    Milstein,
    /// [`srk_platen`] (commutative noise).
    SrkPlaten,
    /// [`heun_stratonovich`].
    HeunStratonovich,
    /// [`milstein_noncommutative`] with [`sample_noise`] (`terms` Fourier modes).
    MilsteinNoncommutative {
        /// Fourier modes of the Levy areas.
        terms: usize,
    },
}

/// Seed of member `index` of an ensemble with base `seed`: a function of the pair only, so the
/// members do not depend on how many threads run them or in which order.
#[must_use]
pub const fn ensemble_seed(
    seed: u64,
    index: usize,
) -> u64 {
    let mut z = seed
        ^ (index as u64)
            .wrapping_add(1)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Paths from which [`sde_ensemble`] runs its members concurrently.
const SDE_ENSEMBLE_PAR_PATHS: usize = 4;

/// Simulates `paths` independent sample paths of `dX = a dt + g dW` with `steps` fixed steps of
/// size `h` and returns the end state of each, in path order.
///
/// Path `i` draws its Brownian increments from `Rng::new(ensemble_seed(seed, i))`, so the result
/// is independent of the thread count and equals running the paths one after the other. The
/// first error (in path order) is returned.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`] from a member.
#[allow(clippy::too_many_arguments)]
pub fn sde_ensemble<A, G>(
    drift: A,
    diff: G,
    m: usize,
    (t0, h, steps): (f64, f64, usize),
    y0: &[f64],
    paths: usize,
    seed: u64,
    scheme: EnsembleScheme,
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]) + Sync,
    G: Fn(f64, &[f64], &mut [f64]) + Sync,
{
    let member = |i: usize| -> Result<Vec<f64>, SdeError> {
        let mut rng = Rng::new(ensemble_seed(seed, i));
        let path = match scheme {
            | EnsembleScheme::MilsteinNoncommutative { terms } => {
                let noise = sample_noise(&mut rng, m, steps, h, terms);
                milstein_noncommutative(&drift, &diff, m, t0, h, y0, &noise)?
            },
            | _ => {
                let dws = brownian_increments(&mut rng, m, steps, h);
                match scheme {
                    | EnsembleScheme::EulerMaruyama => {
                        euler_maruyama(&drift, &diff, m, t0, h, y0, &dws)?
                    },
                    | EnsembleScheme::Milstein => milstein(&drift, &diff, m, t0, h, y0, &dws)?,
                    | EnsembleScheme::SrkPlaten => srk_platen(&drift, &diff, m, t0, h, y0, &dws)?,
                    | _ => heun_stratonovich(&drift, &diff, m, t0, h, y0, &dws)?,
                }
            },
        };
        path.last().cloned().ok_or(SdeError::Invalid)
    };
    crate::kernels::par::map_indexed(paths, SDE_ENSEMBLE_PAR_PATHS, member)
        .into_iter()
        .collect()
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
    sde_adaptive_with(
        drift,
        diff,
        m,
        t0,
        t1,
        y0,
        &SdeAdaptiveOptions {
            rtol,
            atol,
            seed,
            max_steps,
        },
    )
}
