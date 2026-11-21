//! Extended SDE kernels: drift-implicit schemes for stiff SDEs,
//! Poisson jump sampling, and Chambers-Mallows-Stuck alpha-stable sampling.

use std::num::NonZeroUsize;

use crate::kernels::numeric_common::trajectory_sample_count;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::f64_scratch;

/// Error returned by extended SDE kernels.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SdeExtError {
    /// Invalid parameter or dimension.
    Invalid(String),
    /// Non-finite value encountered.
    NonFinite(String),
    /// Step or iteration limit reached.
    StepLimit(String),
}

impl SdeExtError {
    /// Convenience constructor for [`Self::Invalid`].
    #[must_use]
    pub fn invalid(msg: impl Into<String>) -> Self {
        Self::Invalid(msg.into())
    }
}

impl std::fmt::Display for SdeExtError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Invalid(s) => write!(f, "invalid sde parameter: {s}"),
            | Self::NonFinite(s) => write!(f, "non-finite value: {s}"),
            | Self::StepLimit(s) => write!(f, "step limit: {s}"),
        }
    }
}

impl std::error::Error for SdeExtError {}

/// Drift-implicit schemes for stiff SDEs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DriftImplicitScheme {
    /// Drift-implicit backward Euler (`theta = 1.0`):
    /// `Y_{k+1} = Y_k + a(t_{k+1}, Y_{k+1}) dt + G(t_k, Y_k) dW_k`.
    /// L-stable in the drift; strong order 0.5, weak order 1.0.
    BackwardEuler,
    /// Symmetric drift-implicit trapezoidal / Crank–Nicolson (`theta = 0.5`):
    /// `Y_{k+1} = Y_k + 0.5 * (a(t_k, Y_k) + a(t_{k+1}, Y_{k+1})) dt + G(t_k, Y_k) dW_k`.
    /// A-stable in the drift; strong order 0.5, weak order 1.0.
    Trapezoidal,
    /// General `theta`-Maruyama scheme (`theta in [0.0, 1.0]`):
    /// `Y_{k+1} = Y_k + ((1 - theta) a(t_k, Y_k) + theta a(t_{k+1}, Y_{k+1})) dt + G(t_k, Y_k) dW_k`.
    Theta(f64),
    /// Drift-implicit Milstein (`theta = 1.0` on drift + Itô Milstein diffusion correction):
    /// L-stable in the drift and strong order 1.0 on commutative noise.
    DriftImplicitMilstein,
    /// Higham–Mao–Stuart Split-Step Backward Euler (SSBE):
    /// first solves `Y^* = Y_k + a(t_{k+1}, Y^*) dt` implicitly, then applies
    /// `Y_{k+1} = Y^* + G(t_k, Y^*) dW_k`. Strongly convergent for stiff nonlinear
    /// one-sided Lipschitz drifts (e.g. stochastic Ginzburg–Landau).
    SplitStepBackwardEuler,
}

fn solve_linear_sde(
    a: &mut [f64],
    b: &mut [f64],
) -> Result<(), SdeExtError> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))
            .unwrap_or(col);
        let value = a[pivot * n + col];
        if value == 0.0 || !value.is_finite() {
            return Err(SdeExtError::StepLimit(
                "singular Jacobian in drift-implicit Newton solver".into(),
            ));
        }
        if pivot != col {
            let (before, after) = a.split_at_mut(pivot * n);
            before[col * n..(col + 1) * n].swap_with_slice(&mut after[..n]);
            b.swap(col, pivot);
        }
        let (top, below) = a.split_at_mut((col + 1) * n);
        let pivot_row = &top[col * n..];
        let pivot_rhs = b[col];
        for (row, rhs) in below.chunks_exact_mut(n).zip(&mut b[col + 1..]) {
            let factor = row[col] / pivot_row[col];
            for (value, &pivot_value) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *value -= factor * pivot_value;
            }
            *rhs -= factor * pivot_rhs;
        }
    }
    for row in (0..n).rev() {
        let start = row * n;
        let tail: f64 = a[start + row + 1..start + n]
            .iter()
            .zip(&b[row + 1..])
            .map(|(a, b)| a * b)
            .sum();
        b[row] = (b[row] - tail) / a[start + row];
    }
    if b.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::NonFinite(
            "non-finite update in drift-implicit Newton solver".into(),
        ));
    }
    Ok(())
}

#[derive(Default)]
struct NewtonWork {
    rate: Vec<f64>,
    perturbed_rate: Vec<f64>,
    residual: Vec<f64>,
    trial: Vec<f64>,
    jacobian: Option<Scratch<'static, f64>>,
}

fn residuals<'a>(
    state: &'a [f64],
    base: &'a [f64],
    rate: &'a [f64],
    scale: f64,
) -> impl Iterator<Item = f64> + 'a {
    state
        .iter()
        .zip(base)
        .zip(rate)
        .map(move |((&state, &base), &rate)| state - base - scale * rate)
}

fn residual_norm(
    state: &[f64],
    base: &[f64],
    rate: &[f64],
    scale: f64,
    residual: &mut [f64],
) -> Result<(f64, bool), SdeExtError> {
    let mut max = 0.0_f64;
    let mut converged = true;
    for (((slot, value), state), base) in residual
        .iter_mut()
        .zip(residuals(state, base, rate, scale))
        .zip(state)
        .zip(base)
    {
        if !value.is_finite() {
            return Err(SdeExtError::NonFinite(
                "non-finite residual in implicit drift solver".into(),
            ));
        }
        *slot = value;
        max = max.max(value.abs());
        let tolerance = (64.0 * f64::EPSILON * state.abs().max(base.abs())).max(f64::from_bits(1));
        converged &= value.abs() <= tolerance;
    }
    Ok((max, converged))
}

/// Solves `Y - base - scale * a(t_imp, Y) = 0` for `Y` via damped Newton-Raphson.
fn solve_implicit_drift<D>(
    drift: &D,
    t_imp: f64,
    base: &[f64],
    y: &mut [f64],
    scale: f64,
    work: &mut NewtonWork,
) -> Result<(), SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
{
    let n = base.len();
    if scale == 0.0 {
        return Ok(());
    }
    let NewtonWork {
        rate,
        perturbed_rate,
        residual,
        trial,
        jacobian,
    } = work;
    rate.clear();
    rate.resize(n, 0.0);
    residual.resize(n, 0.0);

    for _ in 0..25 {
        drift(t_imp, y, rate);
        let (norm, converged) = residual_norm(y, base, rate, scale, residual)?;
        if converged {
            return Ok(());
        }
        let jac = jacobian.get_or_insert_with(|| {
            let mut jac = f64_scratch(n * n);
            jac.resize(n * n, 0.0);
            jac
        });
        perturbed_rate.clear();
        perturbed_rate.resize(n, 0.0);
        trial.clear();
        trial.extend_from_slice(y);
        for (col, &value) in y.iter().enumerate() {
            let mut eps = f64::EPSILON.sqrt() * value.abs().max(1.0);
            if !(value + eps).is_finite() {
                eps = -eps;
            }
            trial[col] = value + eps;
            drift(t_imp, trial, perturbed_rate);
            trial[col] = value;
            if !perturbed_rate.iter().all(|v| v.is_finite()) {
                return Err(SdeExtError::NonFinite(
                    "drift Jacobian became non-finite".into(),
                ));
            }
            for (row, (perturbed, reference)) in perturbed_rate.iter().zip(rate.iter()).enumerate()
            {
                let da = (perturbed - reference) / eps;
                jac[row * n + col] = if row == col {
                    1.0 - scale * da
                } else {
                    -scale * da
                };
            }
        }
        for value in residual.iter_mut() {
            *value = -*value;
        }
        solve_linear_sde(jac, residual)?;
        let mut alpha = 1.0;
        let mut accepted = false;
        for _ in 0..10 {
            for ((slot, &value), &delta) in trial.iter_mut().zip(y.iter()).zip(residual.iter()) {
                *slot = value + alpha * delta;
            }
            if trial.iter().all(|v| v.is_finite()) {
                drift(t_imp, trial, perturbed_rate);
                let trial_norm = residuals(trial, base, perturbed_rate, scale)
                    .try_fold(0.0_f64, |max, value| {
                        value.is_finite().then_some(max.max(value.abs()))
                    });
                if trial_norm.is_some_and(|trial_norm| trial_norm < norm) {
                    y.copy_from_slice(trial);
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if !accepted {
            for (value, delta) in y.iter_mut().zip(residual.iter()) {
                *value += delta;
            }
            if y.iter().any(|v| !v.is_finite()) {
                return Err(SdeExtError::NonFinite(
                    "non-finite state in implicit drift solver".into(),
                ));
            }
        }
    }
    drift(t_imp, y, rate);
    if residual_norm(y, base, rate, scale, residual)?.1 {
        Ok(())
    } else {
        Err(SdeExtError::StepLimit(
            "drift-implicit Newton iteration did not converge".into(),
        ))
    }
}

/// One drift-implicit step with prescribed Brownian increments and iterated Itô integrals.
///
/// Noncommutative Milstein uses the supplied full iterated-integral matrix;
/// Euler/theta schemes do not need those integrals. Invalid shapes and
/// nonconvergent Newton iterations return errors.
pub fn drift_implicit_step_with_areas<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t: f64,
    h: f64,
    y: &[f64],
    dw: &[f64],
    areas: &[f64],
    scheme: DriftImplicitScheme,
) -> Result<Vec<f64>, SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    validate_workspace(y.len(), noise, scheme)?;
    let theta = step_theta(scheme)?;
    if noise.checked_mul(noise) != Some(areas.len()) || areas.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::invalid(
            "invalid implicit Milstein iterated-integral matrix",
        ));
    }
    if !h.is_finite()
        || h <= 0.0
        || !step_advances(t, h)
        || dw.len() != noise
        || y.iter().chain(dw).any(|value| !value.is_finite())
    {
        return Err(SdeExtError::invalid(
            "invalid drift-implicit step initial data, grid, or increments",
        ));
    }
    let mut next = y.to_vec();
    ImplicitWork::new(y.len(), noise, scheme, theta).step(
        drift,
        diffusion,
        t,
        h,
        y,
        dw,
        Some(areas),
        &mut next,
    )?;
    Ok(next)
}

fn validate_workspace(
    n: usize,
    noise: usize,
    scheme: DriftImplicitScheme,
) -> Result<(), SdeExtError> {
    const MAX_WORKSPACE: usize = 8_388_608;
    if n == 0
        || noise == 0
        || n.checked_mul(noise).is_none_or(|size| size > MAX_WORKSPACE)
        || n.checked_mul(n).is_none_or(|size| size > MAX_WORKSPACE)
        || (scheme == DriftImplicitScheme::DriftImplicitMilstein
            && noise
                .checked_mul(noise)
                .is_none_or(|size| size > MAX_WORKSPACE))
    {
        return Err(SdeExtError::invalid(
            "invalid or oversized drift-implicit step workspace",
        ));
    }
    Ok(())
}

fn step_theta(scheme: DriftImplicitScheme) -> Result<f64, SdeExtError> {
    let theta = match scheme {
        | DriftImplicitScheme::Theta(theta) => theta,
        | DriftImplicitScheme::Trapezoidal => 0.5,
        | _ => 1.0,
    };
    if !(0.0..=1.0).contains(&theta) {
        return Err(SdeExtError::invalid("theta must be in [0,1]"));
    }
    Ok(theta)
}

fn step_advances(
    t: f64,
    h: f64,
) -> bool {
    let next = t + h;
    next.is_finite() && next > t
}

struct ImplicitWork {
    base: Vec<f64>,
    diffusion: Vec<f64>,
    newton: NewtonWork,
    noise: usize,
    scheme: DriftImplicitScheme,
    theta: f64,
}

impl ImplicitWork {
    fn new(
        n: usize,
        noise: usize,
        scheme: DriftImplicitScheme,
        theta: f64,
    ) -> Self {
        Self {
            base: vec![0.0; n],
            diffusion: Vec::new(),
            newton: NewtonWork::default(),
            noise,
            scheme,
            theta,
        }
    }

    fn step<D, F>(
        &mut self,
        drift: &D,
        diffusion: &F,
        t: f64,
        h: f64,
        state: &[f64],
        dw: &[f64],
        areas: Option<&[f64]>,
        next: &mut [f64],
    ) -> Result<(), SdeExtError>
    where
        D: Fn(f64, &[f64], &mut [f64]),
        F: Fn(f64, &[f64], &mut [f64]),
    {
        let base = &mut self.base;
        let g = &mut self.diffusion;
        let work = &mut self.newton;
        base.copy_from_slice(state);
        let split = self.scheme == DriftImplicitScheme::SplitStepBackwardEuler;
        if split {
            solve_implicit_drift(drift, t + h, state, base, h, work)?;
        } else if self.theta < 1.0 {
            work.rate.clear();
            work.rate.resize(state.len(), 0.0);
            drift(t, state, &mut work.rate);
            if work.rate.iter().any(|v| !v.is_finite()) {
                return Err(SdeExtError::NonFinite("non-finite explicit drift".into()));
            }
            let scale = (1.0 - self.theta) * h;
            for (value, rate) in base.iter_mut().zip(&work.rate) {
                *value += scale * rate;
            }
        }
        let noise = self.noise;
        g.clear();
        g.resize(state.len() * noise, 0.0);
        diffusion(t, if split { base } else { state }, g);
        if g.iter().any(|v| !v.is_finite()) {
            return Err(SdeExtError::NonFinite(
                "non-finite implicit diffusion".into(),
            ));
        }
        for (value, row) in base.iter_mut().zip(g.chunks_exact(noise)) {
            for (coefficient, increment) in row.iter().zip(dw) {
                *value += coefficient * increment;
            }
        }
        if self.scheme == DriftImplicitScheme::DriftImplicitMilstein {
            let correction = if let Some(areas) = areas {
                crate::kernels::sde::milstein_term_areas(diffusion, t, state, g, noise, areas)
            } else {
                crate::kernels::sde::milstein_term(diffusion, t, state, g, noise, dw, h)
            };
            for (value, correction) in base.iter_mut().zip(correction) {
                *value += correction;
            }
        }
        if base.iter().any(|v| !v.is_finite()) {
            return Err(SdeExtError::NonFinite("non-finite implicit base".into()));
        }
        next.copy_from_slice(base);
        if !split {
            solve_implicit_drift(drift, t + h, base, next, self.theta * h, work)?;
        }
        Ok(())
    }
}

/// Fixed-step drift-implicit integration over prescribed Wiener increments.
///
/// Rejects inconsistent shapes, invalid timesteps, non-finite coefficients and
/// nonconvergent implicit equations. Milstein includes all commutative cross terms.
pub fn drift_implicit_core<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t0: f64,
    dt: f64,
    y0: &[f64],
    increments: &[Vec<f64>],
    scheme: DriftImplicitScheme,
) -> Result<Vec<Vec<f64>>, SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    drift_implicit_sampled(drift, diffusion, noise, t0, dt, y0, increments, scheme, 1)
}

/// Runs every fixed drift-implicit step while retaining selected states.
///
/// The result starts at `y0`, includes the state after every `save_every`
/// increments, and ends at the final state without duplicating it. A path
/// with no increments returns only `y0`. Saved states use the same arithmetic
/// and callback evaluations as [`drift_implicit_core`]. No randomness is drawn
/// internally. `t0` and `dt` are times, `increments` are Wiener increments,
/// and the states retain the coordinate units of `y0`.
///
/// # Errors
/// Rejects zero retention stride, inconsistent shapes, a non-advancing grid,
/// non-finite coefficients or states, and nonconvergent implicit equations.
pub fn drift_implicit_sampled<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t0: f64,
    dt: f64,
    y0: &[f64],
    increments: &[Vec<f64>],
    scheme: DriftImplicitScheme,
    save_every: usize,
) -> Result<Vec<Vec<f64>>, SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    let theta = step_theta(scheme)?;
    let stride = NonZeroUsize::new(save_every)
        .ok_or_else(|| SdeExtError::invalid("the snapshot stride must be positive"))?;
    if y0.is_empty()
        || noise == 0
        || !t0.is_finite()
        || !dt.is_finite()
        || dt <= 0.0
        || y0.iter().any(|v| !v.is_finite())
    {
        return Err(SdeExtError::invalid(
            "invalid drift-implicit initial data, grid, or save stride",
        ));
    }
    if increments.is_empty() {
        return Ok(vec![y0.to_vec()]);
    }
    validate_workspace(y0.len(), noise, scheme)?;
    let mut previous = f64::NEG_INFINITY;
    for (index, dw) in increments.iter().enumerate() {
        let t = t0 + index as f64 * dt;
        if t <= previous
            || !step_advances(t, dt)
            || dw.len() != noise
            || dw.iter().any(|v| !v.is_finite())
        {
            return Err(SdeExtError::invalid(
                "invalid or non-advancing drift-implicit grid or increments",
            ));
        }
        previous = t;
    }
    let total = increments.len();
    let retained = trajectory_sample_count(total, stride)
        .ok_or_else(|| SdeExtError::invalid("too many drift-implicit snapshots"))?;
    let mut states = Vec::new();
    states.try_reserve_exact(retained).map_err(|_| {
        SdeExtError::invalid("drift-implicit snapshot allocation exceeds available capacity")
    })?;
    states.push(y0.to_vec());
    let mut work = ImplicitWork::new(y0.len(), noise, scheme, theta);
    let mut state = y0.to_vec();
    let mut next = vec![0.0; y0.len()];
    for (index, dw) in increments.iter().enumerate() {
        let t = t0 + index as f64 * dt;
        work.step(drift, diffusion, t, dt, &state, dw, None, &mut next)?;
        std::mem::swap(&mut state, &mut next);
        let completed = index + 1;
        if completed == total || completed.is_multiple_of(save_every) {
            states.push(state.clone());
        }
    }
    Ok(states)
}

/// Exact Poisson sampling: Knuth multiplication for small means and
/// Hoermann transformed rejection (PTRS) for large means.
/// The normal callback is retained for source compatibility and is unused.
pub fn sample_poisson_with(
    next_unit: &mut dyn FnMut() -> f64,
    _next_normal: &mut dyn FnMut() -> f64,
    mean: f64,
) -> usize {
    if mean <= 0.0 || !mean.is_finite() {
        return 0;
    }
    if mean <= 30.0 {
        let threshold = (-mean).exp();
        let mut k = 0usize;
        let mut p = 1.0f64;
        loop {
            k += 1;
            p *= next_unit().max(f64::MIN_POSITIVE);
            if p <= threshold {
                return k - 1;
            }
        }
    } else {
        let b = 0.931 + 2.53 * mean.sqrt();
        let a = -0.059 + 0.02483 * b;
        let inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = next_unit() - 0.5;
            let v = next_unit().max(f64::MIN_POSITIVE);
            let us = 0.5 - u.abs();
            if us <= 0.0 {
                continue;
            }
            let k = ((2.0 * a / us + b) * u + mean + 0.43).floor();
            if k < 0.0 {
                continue;
            }
            if us >= 0.07 && v <= vr {
                return k as usize;
            }
            if us < 0.013 && v > us {
                continue;
            }
            if (v * inv_alpha / (a / (us * us) + b)).ln()
                <= -mean + k * mean.ln() - crate::kernels::special::ln_gamma_numerical(k + 1.0)
            {
                return k as usize;
            }
        }
    }
}

/// Samples an alpha-stable random variable `S_alpha(scale, beta, location)` via
/// Chambers-Mallows-Stuck (1976).
///
/// At `alpha = 2` this uses one standard normal draw. Otherwise it consumes
/// two uniform draws from `[0, 1)` and uses the S1 scale and location
/// parametrization, including its logarithmic correction at `alpha = 1`.
///
/// # Errors
/// Invalid parameters, invalid uniform draws, or a non-finite sample.
pub fn sample_alpha_stable_with(
    next_unit: &mut dyn FnMut() -> f64,
    next_normal: &mut dyn FnMut() -> f64,
    alpha: f64,
    beta: f64,
    scale: f64,
    location: f64,
) -> Result<f64, SdeExtError> {
    if !(0.0 < alpha
        && alpha <= 2.0
        && (-1.0..=1.0).contains(&beta)
        && scale > 0.0
        && scale.is_finite()
        && location.is_finite())
    {
        return Err(SdeExtError::invalid(
            "alpha-stable parameters require alpha in (0, 2], beta in [-1, 1], and scale > 0",
        ));
    }
    if crate::kernels::exact_eq(alpha, 2.0) {
        let normal = next_normal();
        let gaussian_scale = std::f64::consts::SQRT_2 * scale;
        let centered = if gaussian_scale.is_finite() {
            gaussian_scale * normal
        } else {
            (std::f64::consts::SQRT_2 * normal) * scale
        };
        return finite_sample(location + centered);
    }
    let uniform = next_unit();
    let exponential = next_unit();
    if !(0.0..1.0).contains(&uniform) || !(0.0..1.0).contains(&exponential) {
        return Err(SdeExtError::invalid(
            "alpha-stable uniform draws must be in [0,1)",
        ));
    }
    let u = uniform.clamp(1e-15, 1.0 - 1e-15);
    let v = std::f64::consts::PI * (u - 0.5);
    let w = -exponential.max(f64::MIN_POSITIVE).ln();

    let half_pi = std::f64::consts::FRAC_PI_2;
    if crate::kernels::exact_eq(alpha, 1.0) {
        let term1 = (half_pi + beta * v) * v.tan();
        let term2 = beta * ((half_pi * w * v.cos()) / (half_pi + beta * v)).ln();
        let x = (term1 - term2) / half_pi;
        finite_sample(
            scale * x + location + (2.0 / std::f64::consts::PI) * beta * scale * scale.ln(),
        )
    } else {
        let b_tan = beta * (alpha * half_pi).tan();
        let theta0 = b_tan.atan() / alpha;
        let factor = (1.0 + b_tan * b_tan).powf(1.0 / (2.0 * alpha));
        let num = (alpha * (v + theta0)).sin();
        let denom = v.cos().powf(1.0 / alpha);
        let cos_part = ((v - alpha * (v + theta0)).cos() / w).powf((1.0 - alpha) / alpha);
        let x = factor * (num / denom) * cos_part;
        finite_sample(scale * x + location)
    }
}

fn finite_sample(value: f64) -> Result<f64, SdeExtError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SdeExtError::NonFinite(
            "alpha-stable sample is not representable".into(),
        ))
    }
}
