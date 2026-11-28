//! Stochastic differential equations: Euler-Maruyama, Milstein, SRK, Heun and adaptive steps.
//!
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
//! * [`milstein_stratonovich`]: the same scheme for the *Stratonovich*
//!   interpretation, where the iterated integral loses the `-h` term
//!   (`dW_j dW_k / 2` instead of `(dW_j dW_k - delta_jk h) / 2`);
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
//!   pool. Each member uses [`ensemble_seed`] with the common seed and its index,
//!   so the ensemble does not depend on the thread count.
//!
//! * [`milstein_noncommutative_stratonovich`]: the Stratonovich Milstein scheme for
//!   non-commutative noise, with `J_jl = I_jl + delta_jl h / 2 = dW_j dW_l / 2 + A_jl`.
//!
//! * [`srk_noncommutative`] / [`srk_noncommutative_stratonovich`]: the derivative-free
//!   strong order 1 stochastic Runge-Kutta schemes of Kloeden-Platen (11.1.7 / 11.1.10)
//!   for non-commutative noise: supporting values `U^j = Y + a h + g^j sqrt(h)`,
//!   `Y' = Y + a h + sum_j g^j dW_j + h^{-1/2} sum_{j,l} (g^l(U^j) - g^l(Y)) I_jl`
//!   (`m` extra diffusion evaluations per step); the Stratonovich variant uses the
//!   `J_jl` and central supporting values `Y +- g^j sqrt(h)` (`2m` evaluations).
//!
//! * [`sde_adaptive_stratonovich_with`] / [`sde_adaptive_noncommutative_stratonovich_with`]:
//!   adaptive Stratonovich integration **by exact conversion to Ito form**. The drift is
//!   replaced by `a~_i = a_i + (1/2) sum_{j,k} g_{kj} d_k g_{ij}` ([`stratonovich_to_ito_drift`],
//!   directional central differences along the columns of `g`) and the Ito controllers run
//!   on `(a~, g)`. The returned path is the *controlled Ito-Milstein solution of the
//!   equivalent Ito SDE*; it converges to the Stratonovich solution of `(a, g)`.
//!
//! All schemes but the Stratonovich ones ([`heun_stratonovich`],
//! [`milstein_stratonovich`], [`milstein_noncommutative_stratonovich`],
//! [`srk_noncommutative_stratonovich`], and the two adaptive Stratonovich solvers)
//! integrate the Ito SDE; the Stratonovich ones integrate `dX = a dt + g o dW` and
//! converge to the Stratonovich solution, which for `dX = sigma X dW` is
//! `X_0 exp(sigma W_T)` rather than the Ito `X_0 exp(sigma W_T - sigma^2 T / 2)`.
//!
//! Limitation: [`milstein`], [`milstein_stratonovich`], [`srk_platen`],
//! [`heun_stratonovich`] (strong order 1/2 only off the commutative case) and
//! [`sde_adaptive_with`] are strong order one only for commutative noise and do not
//! simulate the Levy areas non-commutative noise needs; every SDE is commutative when
//! `m = 1`, so scalar noise is always covered, and otherwise one needs the
//! `*_noncommutative*` schemes and their Levy areas. The
//! non-commutative adaptive solvers cannot take a step shorter than their grid cell
//! (they report [`SdeError::StepBelowCell`]; use more cells) and, being pathwise
//! Milstein schemes with the Levy areas of [`levy_areas`], inherit their
//! finite-mode accuracy. The adaptive Stratonovich solvers are as accurate as the
//! Ito controller on the converted drift (their drift carries a finite-difference
//! error of order `1e-10`). No scheme here is better than strong order one, and
//! the Milstein variants are weak order one as well (they
//! are not weak order two), so ensemble means converge no faster than `O(h)` in
//! general.
//!
//! # Conventions: Ito or Stratonovich
//!
//! The interpretation is encoded in the scheme's *name*, and the difference is
//! a drift term, not an increment term. For `dX = sigma X dW`:
//!
//! | functions | interpretation | solution |
//! |---|---|---|
//! | [`euler_maruyama`], [`milstein`], [`srk_platen`], [`milstein_noncommutative`], [`srk_noncommutative`], [`sde_adaptive`], [`sde_ensemble`] | Ito | `X_0 exp(sigma W_T - sigma^2 T / 2)` |
//! | [`milstein_stratonovich`], [`heun_stratonovich`], [`milstein_noncommutative_stratonovich`], [`srk_noncommutative_stratonovich`], [`sde_adaptive_stratonovich_with`], [`sde_adaptive_noncommutative_stratonovich_with`] | Stratonovich | `X_0 exp(sigma W_T)` |
//!
//! The increments are convention-free — both interpretations sample the same
//! Brownian path. Rust callers supply coefficient callbacks, and the native
//! ABI uses retained coefficient tapes and noise samples. A caller must not
//! pair a Stratonovich scheme with coefficients
//! meant for Ito (or the reverse) without the `g g' h / 2` conversion; the two
//! then differ at strong order one, which is the same order as the schemes
//! themselves.
//!
//! The noise operators `sde_wiener` and `sde_levy_areas` in
//! [`crate::rules::sde`](mod@crate::rules::sde) draw exactly the increments and iterated integrals the
//! fixed-step solvers accept, so a host can sample once and integrate several
//! schemes on the same path.

use std::fmt;

use crate::kernels::numeric_common::scale_product;
use crate::kernels::random::Rng;

/// Errors of the SDE integrators.
///
/// Each carries the data a caller needs to act on it: which invariant failed,
/// which step became non-finite, which bound was reached. The variants used to
/// be empty, which left a host matching on the variant and on nothing else.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SdeError {
    /// Inconsistent dimensions or parameters; `what` names the input class.
    Invalid {
        /// The invariant that failed, for the diagnostic.
        what: &'static str,
    },
    /// The state became non-finite at the given step index.
    NonFinite {
        /// Index of the step (0-based) at which the state left the reals.
        step: usize,
    },
    /// The adaptive solver reached its step budget.
    StepLimit {
        /// The bound that was hit (`SdeAdaptiveOptions::max_steps`).
        limit: usize,
    },
    /// The adaptive step size shrank below the floating-point floor.
    StepUnderflow {
        /// The time reached when the step could shrink no further.
        t: f64,
        /// The step size that underflowed.
        h: f64,
    },
    /// The non-commutative adaptive solver cannot step below one grid cell.
    StepBelowCell {
        /// Number of cells of its Brownian grid; use more.
        cells: usize,
    },
}

impl fmt::Display for SdeError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Invalid { what } => write!(f, "invalid SDE input: {what}"),
            | Self::NonFinite { step } => {
                write!(f, "the SDE state became non-finite at step {step}")
            },
            | Self::StepLimit { limit } => {
                write!(f, "the SDE step limit of {limit} steps was reached")
            },
            | Self::StepUnderflow { t, h } => {
                write!(f, "the SDE step size underflowed at t = {t} (h = {h:e})")
            },
            | Self::StepBelowCell { cells } => {
                write!(
                    f,
                    "the adaptive non-commutative SDE solver cannot step below one of its \
                     {cells} grid cells; use more cells"
                )
            },
        }
    }
}

impl std::error::Error for SdeError {}

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
    if y0.is_empty() || m == 0 {
        return Err(SdeError::Invalid {
            what: "the state is empty or the noise dimension is zero",
        });
    }
    if dws.iter().any(|d| d.len() != m) {
        return Err(SdeError::Invalid {
            what: "a Brownian increment does not have the noise dimension",
        });
    }
    Ok(())
}

fn finite(
    y: &[f64],
    step: usize,
) -> Result<(), SdeError> {
    if y.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(SdeError::NonFinite { step })
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
    euler_maruyama_core::<_, _, true>(drift, diff, m, t0, h, y0, dws).map(|(_, states)| states)
}

// Keep the path retention choice inside the kernel interface.
#[allow(clippy::single_call_fn)]
pub(crate) fn euler_maruyama_terminal<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<f64>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    euler_maruyama_core::<_, _, false>(drift, diff, m, t0, h, y0, dws).map(|(state, _)| state)
}

fn euler_maruyama_core<A, G, const RETAIN_PATH: bool>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<(Vec<f64>, Vec<Vec<f64>>), SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = if RETAIN_PATH {
        vec![y.clone()]
    } else {
        Vec::new()
    };
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
        finite(&y, k)?;
        if RETAIN_PATH {
            out.push(y.clone());
        }
    }
    Ok((y, out))
}

struct ProbeWidth {
    step: f64,
    sides: f64,
    column_scale: Option<f64>,
}

impl ProbeWidth {
    fn weighted_difference(
        &self,
        plus: f64,
        minus: f64,
        weight: f64,
    ) -> f64 {
        if !(plus.is_finite() && minus.is_finite() && weight.is_finite()) {
            return f64::NAN;
        }
        let difference = plus - minus;
        if difference == 0.0 || weight == 0.0 {
            return 0.0;
        }
        if self.column_scale.is_none() && difference.is_finite() {
            let slope = difference / (self.sides * self.step);
            let term = slope * weight;
            if slope.is_normal() && term.is_normal() {
                return term;
            }
        }
        // Keep the column magnitude and integration weight separate from the
        // quotient so an intermediate overflow or underflow cannot erase it.
        let (difference, sides) = if difference.is_finite() {
            (difference, self.sides)
        } else {
            (0.5 * plus - 0.5 * minus, 0.5 * self.sides)
        };
        let signed_difference = if weight.is_sign_negative() {
            -difference
        } else {
            difference
        };
        scale_product(
            signed_difference,
            &[self.column_scale.unwrap_or(1.0), weight.abs()],
            &[self.step, sides],
        )
        .unwrap_or(f64::NAN)
    }
}

struct DiffusionProbes {
    plus_state: Vec<f64>,
    minus_state: Vec<f64>,
    plus_diffusion: Vec<f64>,
    minus_diffusion: Vec<f64>,
}

impl DiffusionProbes {
    fn new(
        n: usize,
        m: usize,
    ) -> Self {
        Self {
            plus_state: vec![0.0; n],
            minus_state: vec![0.0; n],
            plus_diffusion: vec![0.0; n * m],
            minus_diffusion: vec![0.0; n * m],
        }
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the diffusion call and selected column are explicit"
    )]
    fn measure<G: Fn(f64, &[f64], &mut [f64])>(
        &mut self,
        diff: &G,
        t: f64,
        y: &[f64],
        g: &[f64],
        m: usize,
        j: usize,
        radius: f64,
    ) -> Result<Option<ProbeWidth>, ()> {
        let column = g[j..].iter().step_by(m).copied();
        let square_sum = column.clone().map(|value| value * value).sum::<f64>();
        let direct_step = radius / square_sum.sqrt();
        let (step, column_scale) =
            if square_sum.is_normal() && direct_step.is_normal() && direct_step <= 0.5 * f64::MAX {
                (direct_step, None)
            } else {
                if square_sum.is_nan() {
                    return Err(());
                }
                let peak = column.clone().map(f64::abs).fold(0.0, f64::max);
                if peak == 0.0 {
                    return Ok(None);
                }
                if !peak.is_finite() || !radius.is_finite() {
                    return Err(());
                }
                let norm = column
                    .clone()
                    .map(|value| (value / peak).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (radius / norm, Some(peak))
            };

        let mut plus_finite = true;
        let mut minus_finite = true;
        for (((&value, component), plus), minus) in y
            .iter()
            .zip(column)
            .zip(&mut self.plus_state)
            .zip(&mut self.minus_state)
        {
            let offset = column_scale.map_or_else(
                || step * component,
                |peak| {
                    let ratio = component / peak;
                    if ratio.is_normal() || component == 0.0 {
                        step * ratio
                    } else {
                        scale_product(component, &[step], &[peak]).unwrap_or(f64::NAN)
                    }
                },
            );
            *plus = value + offset;
            *minus = value - offset;
            plus_finite &= plus.is_finite();
            minus_finite &= minus.is_finite();
        }
        if !(plus_finite || minus_finite) {
            return Err(());
        }
        // At a finite range boundary a one-sided probe uses the already
        // evaluated diffusion at the unperturbed state for its other value.
        if plus_finite {
            diff(t, &self.plus_state, &mut self.plus_diffusion);
        } else {
            self.plus_diffusion.copy_from_slice(g);
        }
        if minus_finite {
            diff(t, &self.minus_state, &mut self.minus_diffusion);
        } else {
            self.minus_diffusion.copy_from_slice(g);
        }
        Ok(Some(ProbeWidth {
            step,
            sides: if plus_finite && minus_finite {
                2.0
            } else {
                1.0
            },
            column_scale,
        }))
    }
}

fn probe_radius(y: &[f64]) -> f64 {
    1e-6 * y.iter().map(|v| v.abs()).fold(1.0, f64::max)
}

struct MilsteinWork {
    probes: DiffusionProbes,
    correction: Vec<f64>,
}

impl MilsteinWork {
    fn new(
        n: usize,
        m: usize,
    ) -> Self {
        Self {
            probes: DiffusionProbes::new(n, m),
            correction: vec![0.0; n],
        }
    }

    fn correction<G: Fn(f64, &[f64], &mut [f64])>(
        &mut self,
        diff: &G,
        t: f64,
        y: &[f64],
        g: &[f64],
        m: usize,
        integral: impl Fn(usize, usize) -> f64,
    ) -> &[f64] {
        self.correction.fill(0.0);
        self.probes.plus_diffusion.fill(0.0);
        self.probes.minus_diffusion.fill(0.0);
        let radius = probe_radius(y);
        for j in 0..m {
            let width = match self.probes.measure(diff, t, y, g, m, j, radius) {
                | Ok(Some(width)) => width,
                | Ok(None) => continue,
                | Err(()) => {
                    self.correction.fill(f64::NAN);
                    break;
                },
            };
            // Bounds on the operands keep both the quotient and the weighted
            // term away from the exponent limits. Check the whole direction
            // once so the common accumulation loop has no per-entry branch.
            let direct = width.column_scale.is_none()
                && (1e-50..=1e50).contains(&width.step)
                && self
                    .probes
                    .plus_diffusion
                    .iter()
                    .zip(&self.probes.minus_diffusion)
                    .fold(true, |valid, (plus, minus)| {
                        let magnitude = (plus - minus).abs();
                        valid & (magnitude == 0.0 || (1e-100..=1e100).contains(&magnitude))
                    });
            let denominator = width.sides * width.step;
            for k in 0..m {
                let weight = integral(j, k);
                if direct && (weight == 0.0 || (1e-50..=1e50).contains(&weight.abs())) {
                    for (i, correction) in self.correction.iter_mut().enumerate() {
                        let at = i * m + k;
                        let difference =
                            self.probes.plus_diffusion[at] - self.probes.minus_diffusion[at];
                        *correction += difference / denominator * weight;
                    }
                } else {
                    for (i, correction) in self.correction.iter_mut().enumerate() {
                        let at = i * m + k;
                        *correction += width.weighted_difference(
                            self.probes.plus_diffusion[at],
                            self.probes.minus_diffusion[at],
                            weight,
                        );
                    }
                }
            }
        }
        &self.correction
    }
}

/// Milstein correction `sum_{j,k} (D_j g)_{.k} I_jk` for the iterated Ito
/// integrals `areas[j * m + k] = I_jk` (row-major).
pub(crate) fn milstein_term_areas<G: Fn(f64, &[f64], &mut [f64])>(
    diff: &G,
    t: f64,
    y: &[f64],
    g: &[f64],
    m: usize,
    areas: &[f64],
) -> Vec<f64> {
    let mut work = MilsteinWork::new(y.len(), m);
    work.correction(diff, t, y, g, m, |j, k| areas[j * m + k]);
    work.correction
}

fn ito_product(
    left: f64,
    right: f64,
    diagonal: f64,
) -> f64 {
    let difference = left * right - diagonal;
    if difference.is_finite() {
        0.5 * difference
    } else {
        (0.5 * left).mul_add(right, -0.5 * diagonal)
    }
}

fn stratonovich_product(
    left: f64,
    right: f64,
) -> f64 {
    let half = 0.5 * left;
    if half.is_normal() || left == 0.0 {
        half * right
    } else {
        0.5 * (left * right)
    }
}

/// One Milstein correction `(1/2) sum_{j,k} (D_j g)_{.k} (dW_j dW_k - delta_jk h)`
/// for commutative noise, by central differences along the columns of `g`.
pub(crate) fn milstein_term<G: Fn(f64, &[f64], &mut [f64])>(
    diff: &G,
    t: f64,
    y: &[f64],
    g: &[f64],
    m: usize,
    dw: &[f64],
    h: f64,
) -> Vec<f64> {
    let mut work = MilsteinWork::new(y.len(), m);
    work.correction(diff, t, y, g, m, |j, k| {
        ito_product(dw[j], dw[k], if j == k { h } else { 0.0 })
    });
    work.correction
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
    milstein_core::<_, _, true>(drift, diff, m, t0, h, y0, dws).map(|(_, states)| states)
}

// Keep the path retention choice inside the kernel interface.
#[allow(clippy::single_call_fn)]
pub(crate) fn milstein_terminal<A, G>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<Vec<f64>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    milstein_core::<_, _, false>(drift, diff, m, t0, h, y0, dws).map(|(state, _)| state)
}

fn milstein_core<A, G, const RETAIN_PATH: bool>(
    drift: A,
    diff: G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    dws: &[Vec<f64>],
) -> Result<(Vec<f64>, Vec<Vec<f64>>), SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check(y0, m, dws)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = if RETAIN_PATH {
        vec![y.clone()]
    } else {
        Vec::new()
    };
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut work = MilsteinWork::new(n, m);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, k| {
            ito_product(dw[j], dw[k], if j == k { h } else { 0.0 })
        });
        for i in 0..n {
            let mut s = a[i] * h + corr[i];
            for j in 0..m {
                s += g[i * m + j] * dw[j];
            }
            y[i] += s;
        }
        finite(&y, k)?;
        if RETAIN_PATH {
            out.push(y.clone());
        }
    }
    Ok((y, out))
}

/// Milstein scheme (strong order 1) for the *Stratonovich* interpretation
/// `dX = a dt + g o dW` of the same coefficients; arguments as in
/// [`euler_maruyama`].
///
/// It differs from the Ito [`milstein`] only in the iterated integral of the
/// correction, `I_jk = dW_j dW_k / 2` instead of
/// `(dW_j dW_k - delta_jk h) / 2`: for scalar noise the correction is
/// `g g' dW^2 / 2` rather than `g g' (dW^2 - h) / 2`, so the two differ by the
/// Ito-Stratonovich drift correction `g g' h / 2` and converge to different
/// solutions. For `dX = sigma X dW` they are `X_0 exp(sigma W_T)` (here) and
/// `X_0 exp(sigma W_T - sigma^2 T / 2)` ([`milstein`]).
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn milstein_stratonovich<A, G>(
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
    let mut work = MilsteinWork::new(n, m);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, k| {
            stratonovich_product(dw[j], dw[k])
        });
        for (((yi, ai), ci), grow) in y.iter_mut().zip(&a).zip(corr).zip(g.chunks(m)) {
            *yi += ai * h + ci + grow.iter().zip(dw).map(|(gij, dwj)| gij * dwj).sum::<f64>();
        }
        finite(&y, k)?;
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
///
/// `dw` contains the `m` supplied increments. The returned matrix is flat
/// and row-major. Pair evaluation can run in parallel, with the same RNG
/// draws and floating-point summation order for every worker count.
#[must_use]
pub fn levy_areas(
    rng: &mut Rng,
    m: usize,
    h: f64,
    dw: &[f64],
    terms: usize,
) -> Vec<f64> {
    const PAIRS_PER_CHUNK: usize = 128;
    let p = terms.max(1);
    if m == 0 {
        return Vec::new();
    }
    if m == 1 {
        // A scalar integral has no skew part. Keep the original coefficient
        // draws so subsequent steps see the same normal stream.
        for _ in 0..p {
            let _ = rng.normal();
            let _ = rng.normal();
        }
        let _ = rng.normal();
        return vec![ito_product(dw[0], dw[0], h)];
    }
    let sh = h.sqrt();
    let xi: Vec<f64> = dw.iter().map(|d| d / sh).collect();
    let srho = (zeta2_tail(p) / (2.0 * std::f64::consts::PI.powi(2))).sqrt();
    // Draw all zeta coefficients, then all eta coefficients, in row order.
    // Cache the shifted eta used in every pair incident to a given row.
    let count = m * p;
    let mut coefficients = crate::kernels::reuse::f64_scratch(3 * count);
    coefficients.extend((0..2 * count).map(|_| rng.normal()));
    coefficients.resize(3 * count, 0.0);
    let (zeta, other) = coefficients.split_at_mut(count);
    let (eta, shifted) = other.split_at_mut(count);
    for ((shifted, eta), &increment) in shifted
        .chunks_exact_mut(p)
        .zip(eta.chunks_exact(p))
        .zip(&xi)
    {
        let offset = std::f64::consts::SQRT_2 * increment;
        for (value, &entry) in shifted.iter_mut().zip(eta) {
            *value = offset + entry;
        }
    }
    let mu: Vec<f64> = (0..m).map(|_| rng.normal()).collect();
    let bridge = LevyBridge {
        h,
        sh,
        srho,
        p,
        dw,
        xi: &xi,
        mu: &mu,
        zeta,
        eta,
        shifted,
    };
    let mut areas = vec![0.0; m * m];
    let pairs = m * (m - 1) / 2;
    if pairs.saturating_mul(p) >= 65_536
        && crate::execution::should_split(pairs.div_ceil(PAIRS_PER_CHUNK), 4)
    {
        // Draw the independent tail normals in the original row-major order.
        // Each disjoint chunk then replaces its draws with both integrals.
        let mut packed = crate::kernels::reuse::f64_scratch(2 * pairs);
        for _ in 0..pairs {
            packed.extend([rng.normal(), 0.0]);
        }
        crate::kernels::par::for_chunks_mut(
            &mut packed,
            2 * PAIRS_PER_CHUNK,
            4,
            |chunk, values| {
                let start = chunk * PAIRS_PER_CHUNK;
                let (mut first, mut end) = (0, m - 1);
                while first + 1 < end {
                    let middle = usize::midpoint(first, end);
                    if levy_row_start(m, middle) <= start {
                        first = middle;
                    } else {
                        end = middle;
                    }
                }
                let mut j = first;
                let mut l = j + 1 + start - levy_row_start(m, j);
                for pair in values.as_chunks_mut::<2>().0 {
                    *pair = bridge.pair(j, l, pair[0]);
                    l += 1;
                    if l == m {
                        j += 1;
                        l = j + 1;
                    }
                }
            },
        );
        let mut pairs = packed.as_chunks::<2>().0.iter();
        for j in 0..m {
            areas[j * m + j] = ito_product(dw[j], dw[j], h);
            for (l, pair) in (j + 1..m).zip(&mut pairs) {
                areas[j * m + l] = pair[0];
                areas[l * m + j] = pair[1];
            }
        }
        return areas;
    }
    for j in 0..m {
        areas[j * m + j] = ito_product(dw[j], dw[j], h);
        for l in j + 1..m {
            let [forward, reverse] = bridge.pair(j, l, rng.normal());
            areas[j * m + l] = forward;
            areas[l * m + j] = reverse;
        }
    }
    areas
}

const fn levy_row_start(
    m: usize,
    row: usize,
) -> usize {
    row * (2 * m - row - 1) / 2
}

struct LevyBridge<'a> {
    h: f64,
    sh: f64,
    srho: f64,
    p: usize,
    dw: &'a [f64],
    xi: &'a [f64],
    mu: &'a [f64],
    zeta: &'a [f64],
    eta: &'a [f64],
    shifted: &'a [f64],
}

impl LevyBridge<'_> {
    fn pair(
        &self,
        j: usize,
        l: usize,
        tail_normal: f64,
    ) -> [f64; 2] {
        let Self {
            h,
            sh,
            srho,
            p,
            dw,
            xi,
            mu,
            zeta,
            eta,
            shifted,
        } = *self;
        let row_j = j * p..(j + 1) * p;
        let row_l = l * p..(l + 1) * p;
        let (zj, zl) = (&zeta[row_j.clone()], &zeta[row_l.clone()]);
        let (ej, el) = (&shifted[row_j.clone()], &shifted[row_l.clone()]);
        let mut a = 0.0;
        for (r, ((&zj, &zl), (&ej, &el))) in zj.iter().zip(zl).zip(ej.iter().zip(el)).enumerate() {
            a += (zj * el - zl * ej) / (r + 1) as f64;
        }
        a /= 2.0 * std::f64::consts::PI;
        let ijl = h
            * (0.5 * xi[j] * xi[l]
                + a
                + srho * (mu[j] * xi[l] - mu[l] * xi[j])
                + srho * tail_normal);
        let ilj = dw[j] * dw[l] - ijl;
        if ijl.is_finite() && ilj.is_finite() || !h.is_finite() || h <= 0.0 {
            return [ijl, ilj];
        }
        // The normalized increments or their product can overflow
        // even when the physical integrals are both representable.
        let symmetric = stratonovich_product(dw[j], dw[l]);
        let (ej, el) = (&eta[row_j], &eta[row_l]);
        let mut skew = 0.0;
        for (r, ((&zj, &zl), (&ej, &el))) in zj.iter().zip(zl).zip(ej.iter().zip(el)).enumerate() {
            let denominator = (r + 1) as f64 * (2.0 * std::f64::consts::PI);
            let increment_scale = std::f64::consts::SQRT_2 * (sh / denominator);
            let mode_scale = h / denominator;
            let left = dw[j] * increment_scale + ej * mode_scale;
            let right = dw[l] * increment_scale + el * mode_scale;
            skew += zj * right - zl * left;
        }
        let tail_scale = srho * sh;
        skew += mu[j] * (dw[l] * tail_scale) - mu[l] * (dw[j] * tail_scale);
        skew += (h * srho) * tail_normal;
        [symmetric + skew, symmetric - skew]
    }
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
    check_sample(y0, m, noise)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut work = MilsteinWork::new(n, m);
    for (k, (dw, ar)) in noise.dws.iter().zip(&noise.areas).enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, l| ar[j * m + l]);
        for i in 0..n {
            let mut s = a[i] * h + corr[i];
            for j in 0..m {
                s += g[i * m + j] * dw[j];
            }
            y[i] += s;
        }
        finite(&y, k)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Validates a [`NoiseSample`] against the state and the noise dimension.
fn check_sample(
    y0: &[f64],
    m: usize,
    noise: &NoiseSample,
) -> Result<(), SdeError> {
    check(y0, m, &noise.dws)?;
    if noise.areas.len() != noise.dws.len() || noise.areas.iter().any(|a| a.len() != m * m) {
        return Err(SdeError::Invalid {
            what: "the Levy areas do not match the increments (one m*m block per step)",
        });
    }
    Ok(())
}

/// Stratonovich iterated integrals `J_jl = I_jl + delta_jl h / 2` from the Ito ones.
fn stratonovich_areas(
    ito: &[f64],
    m: usize,
    h: f64,
) -> Vec<f64> {
    let mut j = ito.to_vec();
    for k in 0..m {
        j[k * m + k] += 0.5 * h;
    }
    j
}

/// Milstein scheme (strong order 1) for **non-commutative** noise in the
/// *Stratonovich* interpretation `dX = a dt + g o dW`.
///
/// The correction is `sum_{j,l} (D_j g)_{.l} J_jl` with the Stratonovich iterated
/// integrals `J_jl = I_jl + delta_jl h / 2 = dW_j dW_l / 2 + A_jl` (`A` the Levy
/// area), i.e. the same sampled areas as [`milstein_noncommutative`] with the
/// diagonal shifted by `h / 2`. The drift is the Stratonovich drift, unchanged.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn milstein_noncommutative_stratonovich<A, G>(
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
    check_sample(y0, m, noise)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut work = MilsteinWork::new(n, m);
    for (k, (dw, ar)) in noise.dws.iter().zip(&noise.areas).enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, l| {
            let value = ar[j * m + l];
            if j == l {
                value + 0.5 * h
            } else {
                value
            }
        });
        for i in 0..n {
            let mut s = a[i] * h + corr[i];
            for l in 0..m {
                s += g[i * m + l] * dw[l];
            }
            y[i] += s;
        }
        finite(&y, k)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Shared body of [`srk_noncommutative`] and [`srk_noncommutative_stratonovich`].
#[allow(clippy::too_many_arguments)] // the kernel's common argument list plus the calculus flag
fn srk_noncommutative_run<A, G>(
    drift: &A,
    diff: &G,
    m: usize,
    t0: f64,
    h: f64,
    y0: &[f64],
    noise: &NoiseSample,
    stratonovich: bool,
) -> Result<Vec<Vec<f64>>, SdeError>
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    check_sample(y0, m, noise)?;
    let n = y0.len();
    let sh = h.abs().sqrt();
    let mut y = y0.to_vec();
    let mut out = vec![y.clone()];
    let (mut a, mut g, mut gs, mut gm) = (
        vec![0.0; n],
        vec![0.0; n * m],
        vec![0.0; n * m],
        vec![0.0; n * m],
    );
    for (k, (dw, ar)) in noise.dws.iter().zip(&noise.areas).enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let areas = if stratonovich {
            stratonovich_areas(ar, m, h)
        } else {
            ar.clone()
        };
        let mut corr = vec![0.0; n];
        for j in 0..m {
            if stratonovich {
                // Central supporting values: the mean of J_jj is h / 2, so the
                // difference quotient must be accurate to O(h) (a forward one leaves
                // a systematic O(sqrt(h)) error and only strong order 1/2).
                let up: Vec<f64> = (0..n).map(|i| y[i] + g[i * m + j] * sh).collect();
                let um: Vec<f64> = (0..n).map(|i| y[i] - g[i * m + j] * sh).collect();
                diff(t, &up, &mut gs);
                diff(t, &um, &mut gm);
                for l in 0..m {
                    let w = areas[j * m + l] / (2.0 * sh);
                    for i in 0..n {
                        corr[i] += (gs[i * m + l] - gm[i * m + l]) * w;
                    }
                }
            } else {
                let ups: Vec<f64> = (0..n)
                    .map(|i| y[i] + a[i] * h + g[i * m + j] * sh)
                    .collect();
                diff(t, &ups, &mut gs);
                for l in 0..m {
                    let w = areas[j * m + l] / sh;
                    for i in 0..n {
                        corr[i] += (gs[i * m + l] - g[i * m + l]) * w;
                    }
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
        finite(&y, k)?;
        out.push(y.clone());
    }
    Ok(out)
}

/// Derivative-free stochastic Runge-Kutta scheme of strong order 1 for
/// **non-commutative** Ito noise (Kloeden-Platen, eq. 11.1.7).
///
/// With the supporting values `U^j = Y + a h + g^j sqrt(h)` (`g^j` the `j`-th
/// column of `g`, one per noise channel),
///
/// `Y' = Y + a h + sum_j g^j dW_j + h^{-1/2} sum_{j,l} (g^l(U^j) - g^l(Y)) I_jl`,
///
/// where `I_jl` are the sampled iterated Ito integrals of `noise`
/// ([`sample_noise`]). It needs `m` extra diffusion evaluations per step and no
/// derivative of `g`. On commutative noise the Levy-area part cancels and it
/// agrees with [`srk_platen`] up to the `O(h)` of the supporting values.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn srk_noncommutative<A, G>(
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
    srk_noncommutative_run(&drift, &diff, m, t0, h, y0, noise, false)
}

/// The Stratonovich counterpart of [`srk_noncommutative`].
///
/// After Kloeden-Platen, eq. 11.1.10: with the Stratonovich iterated integrals
/// `J_jl = I_jl + delta_jl h / 2` and the Stratonovich drift it converges to the
/// solution of `dX = a dt + g o dW` with strong order 1.
///
/// Because `J_jj` has mean `h / 2` (unlike `I_jj`), the difference quotient must
/// be accurate to `O(h)`: the scheme uses the *central* supporting values
/// `U^j_{+-} = Y +- g^j sqrt(h)` and
/// `Y' = Y + a h + sum_j g^j dW_j + (2 sqrt(h))^{-1} sum_{j,l} (g^l(U^j_+) - g^l(U^j_-)) J_jl`
/// (a forward quotient would leave a systematic `O(sqrt(h))` error and only
/// strong order 1/2). It costs `2m` extra diffusion evaluations per step.
///
/// # Errors
/// [`SdeError::Invalid`] or [`SdeError::NonFinite`].
pub fn srk_noncommutative_stratonovich<A, G>(
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
    srk_noncommutative_run(&drift, &diff, m, t0, h, y0, noise, true)
}

/// The Ito-equivalent drift of a Stratonovich SDE `dX = a dt + g o dW`:
/// `a~_i = a_i + (1/2) sum_j sum_k g_{kj} d_k g_{ij}`.
///
/// The correction is `(1/2) sum_j (D_j g)_{.j}` with the directional derivative of
/// column `j` along itself, with the same scaled probes as the Milstein term.
/// The ordinary step is `eps = 1e-6 max(1, |y|_inf) / |g^j|`. A probe uses one
/// side at the floating-point range boundary when only that side is finite. The closure writes
/// `a~` into its output, so any Ito integrator applied to `(a~, g)` converges to
/// the Stratonovich solution of `(a, g)`.
pub fn stratonovich_to_ito_drift<'a, A, G>(
    drift: &'a A,
    diff: &'a G,
    m: usize,
) -> impl Fn(f64, &[f64], &mut [f64]) + 'a
where
    A: Fn(f64, &[f64], &mut [f64]),
    G: Fn(f64, &[f64], &mut [f64]),
{
    move |t: f64, y: &[f64], out: &mut [f64]| {
        drift(t, y, out);
        let n = y.len();
        let mut g = vec![0.0; n * m];
        diff(t, y, &mut g);
        let mut probes = DiffusionProbes::new(n, m);
        let radius = probe_radius(y);
        for j in 0..m {
            let width = match probes.measure(diff, t, y, &g, m, j, radius) {
                | Ok(Some(width)) => width,
                | Ok(None) => continue,
                | Err(()) => {
                    out.fill(f64::NAN);
                    return;
                },
            };
            for (i, value) in out.iter_mut().enumerate() {
                let at = i * m + j;
                *value += width.weighted_difference(
                    probes.plus_diffusion[at],
                    probes.minus_diffusion[at],
                    0.5,
                );
            }
        }
    }
}

// Keep both the diffusion difference and the iterated integral in range
// until the final correction is restored to its original scale.
#[cold]
fn platen_scaled_correction(
    support: f64,
    base: f64,
    (weight, scales): (f64, [f64; 2]),
    denominator: f64,
) -> f64 {
    let difference = support - base;
    let (difference, scale) = if difference.is_infinite() {
        (0.5 * support - 0.5 * base, 2.0)
    } else {
        (difference, 1.0)
    };
    if weight == 0.0 {
        return difference * weight / denominator;
    }
    let signed = if weight.is_sign_negative() {
        -difference
    } else {
        difference
    };
    scale_product(
        signed,
        &[weight.abs(), scales[0], scales[1], scale],
        &[denominator],
    )
    .unwrap_or(f64::NAN)
}

#[cold]
fn platen_scaled_integral(
    left: f64,
    right: f64,
    diagonal: f64,
    ito: f64,
) -> (f64, [f64; 2]) {
    const UP: f64 = f64::from_bits((1023 + 600) << 52);
    const DOWN: f64 = f64::from_bits((1023 - 600) << 52);

    if diagonal == 0.0 {
        let signed = if right.is_sign_negative() {
            -left
        } else {
            left
        };
        return (signed, [right.abs(), 1.0]);
    }
    // A diagonal term has left == right. Scaling both factors by the same
    // power of two also scales the subtraction, including tiny cancellations.
    let (scale, inverse) = if ito.is_finite() {
        (UP, DOWN)
    } else {
        (DOWN, UP)
    };
    (
        (scale * left).mul_add(scale * right, -(scale * diagonal) * scale),
        [inverse; 2],
    )
}

fn platen_bounded_diffusion(coefficients: &[f64]) -> bool {
    // For nonzero magnitudes in these exponent bins (2^-250 through 2^251),
    // every nonzero difference needs at most exponents -302 through 251.
    coefficients.iter().fold(true, |bounded, value| {
        let magnitude = value.abs().to_bits();
        let exponent = magnitude >> 52;
        bounded & (magnitude == 0 || ((1023 - 250)..=(1023 + 250)).contains(&exponent))
    })
}

/// Platen's derivative-free stochastic Runge-Kutta scheme for commutative
/// Itô noise, with strong order 1.
///
/// Each supplied Brownian increment has `m` components. The scheme replaces
/// the Milstein term with differences of the diffusion at supporting values
/// `x + a h + g_j sqrt(abs(h))`.
/// Corrections at extreme scales retain an out-of-range intermediate product
/// in scaled form when the weighted result is representable. An empty
/// increment sequence returns just `y0` without evaluating either callback.
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
    if dws.is_empty() {
        return Ok(vec![y0.to_vec()]);
    }
    let n = y0.len();
    let sh = h.abs().sqrt();
    let denominator = 2.0 * sh;
    let mut y = y0.to_vec();
    let mut out = Vec::with_capacity(dws.len() + 1);
    out.push(y.clone());
    let (mut a, mut g, mut gs) = (vec![0.0; n], vec![0.0; n * m], vec![0.0; n * m]);
    let mut scratch = super::reuse::f64_scratch(2 * n);
    scratch.resize(2 * n, 0.0);
    let (corr, ups) = scratch.split_at_mut(n);
    for (k, dw) in dws.iter().enumerate() {
        let t = t0 + k as f64 * h;
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        corr.fill(0.0);
        let bounded_base = n >= 32 && m >= 2 && platen_bounded_diffusion(&g);
        for j in 0..m {
            for (i, value) in ups.iter_mut().enumerate() {
                *value = y[i] + a[i] * h + g[i * m + j] * sh;
            }
            diff(t, ups, &mut gs);
            let bounded_direction = bounded_base && platen_bounded_diffusion(&gs);
            for kk in 0..m {
                let diagonal = if j == kk { h } else { 0.0 };
                let noise_product = dw[j] * dw[kk];
                let ito = noise_product - diagonal;
                let direct = ito.is_normal()
                    || (ito == 0.0 && (noise_product.is_normal() || dw[j] == 0.0 || dw[kk] == 0.0));
                // A bounded direction and weight cannot overflow or underflow
                // their product. Keep this state loop free of per-entry checks.
                if bounded_direction && (1e-100..=1e100).contains(&ito.abs()) {
                    for (i, value) in corr.iter_mut().enumerate() {
                        let at = i * m + kk;
                        *value += (gs[at] - g[at]) * ito / denominator;
                    }
                } else if direct {
                    for (i, value) in corr.iter_mut().enumerate() {
                        let at = i * m + kk;
                        let difference = gs[at] - g[at];
                        let product = difference * ito;
                        *value += if product.is_normal() || difference == 0.0 {
                            product / denominator
                        } else {
                            platen_scaled_correction(gs[at], g[at], (ito, [1.0; 2]), denominator)
                        };
                    }
                } else {
                    let weight = platen_scaled_integral(dw[j], dw[kk], diagonal, ito);
                    for (i, value) in corr.iter_mut().enumerate() {
                        let at = i * m + kk;
                        *value += platen_scaled_correction(gs[at], g[at], weight, denominator);
                    }
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
        finite(&y, k)?;
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
        finite(&y, k)?;
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
    /// Number of accepted steps (`t.len() - 1`).
    pub accepted: usize,
    /// Step size of the last accepted step; `0.0` when none was accepted.
    pub last_h: f64,
    /// Seed the path's Brownian sample was drawn from.
    ///
    /// A host that records it can reproduce the *increments* with
    /// `sde_wiener`-style sampling if it also knows the step sequence; the
    /// adaptive solver's own sequence depends on the error control, so this is
    /// provenance, not a recipe for an independent re-run.
    pub seed: u64,
}

impl SdePath {
    /// A path that has accepted no step yet: origin `(t0, y0)`, drawn from
    /// `seed`.
    fn new(
        seed: u64,
        t0: f64,
        y0: &[f64],
    ) -> Self {
        Self {
            t: vec![t0],
            y: vec![y0.to_vec()],
            rejected: 0,
            accepted: 0,
            last_h: 0.0,
            seed,
        }
    }
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
        return Err(SdeError::Invalid {
            what: "empty state, zero noise dimension, inverted interval or non-positive tolerance",
        });
    }
    let n = y0.len();
    let mut rng = Rng::new(seed);
    let mut stack = BrownianStack { pending: Vec::new() };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut path = SdePath::new(seed, t0, y0);
    let mut h = ((t1 - t0) * 1e-3).max(1e-9);
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut work = MilsteinWork::new(n, m);
    let mut a2 = vec![0.0; n];
    let mut count = 0;
    while t1 - t > 1e-14 * (1.0 + t1.abs()) {
        count += 1;
        if count > max_steps {
            return Err(SdeError::StepLimit { limit: max_steps });
        }
        h = h.min(t1 - t);
        let dw = stack.take(h, m, &mut rng);
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, k| {
            ito_product(dw[j], dw[k], if j == k { h } else { 0.0 })
        });
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
            return Err(SdeError::NonFinite { step: count - 1 });
        }
        if err <= 1.0 {
            t += h;
            y = ynew;
            path.t.push(t);
            path.y.push(y.clone());
            path.accepted += 1;
            path.last_h = h;
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
                return Err(SdeError::StepUnderflow { t, h });
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
        return Err(SdeError::Invalid {
            what: "empty state, zero noise dimension, inverted interval, zero cells or \
                   non-positive tolerance",
        });
    }
    let n = y0.len();
    let hc = (t1 - t0) / opts.cells as f64;
    let mut grid = NoiseGrid::new(seed, m, hc, opts.terms);
    let mut cell = 0_usize;
    let mut y = y0.to_vec();
    let mut path = SdePath::new(seed, t0, y0);
    let mut size = (opts.cells / 1000).max(1);
    let (mut a, mut g) = (vec![0.0; n], vec![0.0; n * m]);
    let mut work = MilsteinWork::new(n, m);
    let mut a2 = vec![0.0; n];
    let mut count = 0;
    while cell < opts.cells {
        count += 1;
        if count > max_steps {
            return Err(SdeError::StepLimit { limit: max_steps });
        }
        let k = size.clamp(1, opts.cells - cell);
        let h = k as f64 * hc;
        let t = t0 + cell as f64 * hc;
        let (dw, ar) = grid.span(cell, k);
        drift(t, &y, &mut a);
        diff(t, &y, &mut g);
        let corr = work.correction(&diff, t, &y, &g, m, |j, l| ar[j * m + l]);
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
            return Err(SdeError::NonFinite { step: count - 1 });
        }
        if err <= 1.0 {
            cell += k;
            y = ynew;
            path.t.push(t0 + cell as f64 * hc);
            path.y.push(y.clone());
            path.accepted += 1;
            path.last_h = h;
            let fac = if err == 0.0 {
                2.0
            } else {
                (0.9 / err.sqrt()).clamp(0.5, 2.0)
            };
            size = ((k as f64 * fac).round() as usize).max(1);
        } else {
            if k == 1 {
                return Err(SdeError::StepBelowCell { cells: opts.cells });
            }
            path.rejected += 1;
            let fac = (0.9 / err.sqrt()).clamp(0.1, 0.9);
            size = ((k as f64 * fac).floor() as usize).clamp(1, k - 1);
        }
    }
    Ok(path)
}

/// Adaptive-step integration of a **Stratonovich** SDE `dX = a dt + g o dW` with
/// commutative noise on `[t0, t1]`, by exact conversion to Ito form.
///
/// The drift is replaced by the Ito-equivalent drift
/// ([`stratonovich_to_ito_drift`]) and [`sde_adaptive_with`] runs on the Ito SDE
/// `dX = a~ dt + g dW`, which has the same solution process. The returned path is
/// therefore the *controlled Ito-Milstein solution of the equivalent Ito SDE*; it
/// converges to the Stratonovich solution of the original coefficients (for
/// `dX = sigma X o dW` to `X_0 exp(sigma W_T)`), with the error control of
/// [`sde_adaptive_with`] applied to the converted drift. The conversion uses
/// finite-difference derivatives of `g`, which costs `2m + 1` extra diffusion
/// evaluations per drift evaluation.
///
/// # Errors
/// As [`sde_adaptive_with`].
pub fn sde_adaptive_stratonovich_with<A, G>(
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
    let ito = stratonovich_to_ito_drift(&drift, &diff, m);
    sde_adaptive_with(ito, &diff, m, t0, t1, y0, opts)
}

/// Adaptive-step integration of a **Stratonovich** SDE with **non-commutative**
/// noise on `[t0, t1]`, by exact conversion to Ito form.
///
/// As [`sde_adaptive_stratonovich_with`], with the controller of
/// [`sde_adaptive_noncommutative_with`] (fine cell grid, Chen's relation); the
/// returned path is the controlled Ito-Milstein solution of the equivalent Ito
/// SDE, which converges to the Stratonovich solution. The Brownian path is the
/// same as for the Ito solver ([`noncommutative_noise_grid`]).
///
/// # Errors
/// As [`sde_adaptive_noncommutative_with`].
pub fn sde_adaptive_noncommutative_stratonovich_with<A, G>(
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
    let ito = stratonovich_to_ito_drift(&drift, &diff, m);
    sde_adaptive_noncommutative_with(ito, &diff, m, t0, t1, y0, opts)
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
    /// [`milstein_noncommutative_stratonovich`] with [`sample_noise`].
    MilsteinNoncommutativeStratonovich {
        /// Fourier modes of the Levy areas.
        terms: usize,
    },
    /// [`srk_noncommutative`] with [`sample_noise`].
    SrkNoncommutative {
        /// Fourier modes of the Levy areas.
        terms: usize,
    },
    /// [`srk_noncommutative_stratonovich`] with [`sample_noise`].
    SrkNoncommutativeStratonovich {
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
            | EnsembleScheme::MilsteinNoncommutativeStratonovich { terms } => {
                let noise = sample_noise(&mut rng, m, steps, h, terms);
                milstein_noncommutative_stratonovich(&drift, &diff, m, t0, h, y0, &noise)?
            },
            | EnsembleScheme::SrkNoncommutative { terms } => {
                let noise = sample_noise(&mut rng, m, steps, h, terms);
                srk_noncommutative(&drift, &diff, m, t0, h, y0, &noise)?
            },
            | EnsembleScheme::SrkNoncommutativeStratonovich { terms } => {
                let noise = sample_noise(&mut rng, m, steps, h, terms);
                srk_noncommutative_stratonovich(&drift, &diff, m, t0, h, y0, &noise)?
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
        path.last().cloned().ok_or(SdeError::Invalid {
            what: "an ensemble path came back empty",
        })
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
