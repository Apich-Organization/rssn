//! Extended SDE kernels: drift-implicit schemes for stiff SDEs,
//! Poisson jump sampling, and Chambers-Mallows-Stuck alpha-stable sampling.

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
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Result<Vec<f64>, SdeExtError> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap_or(col);
        if a[pivot][col].abs() < 1e-300 || !a[pivot][col].is_finite() {
            return Err(SdeExtError::StepLimit(
                "singular Jacobian in drift-implicit Newton solver".into(),
            ));
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            for k in col..n {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::NonFinite(
            "non-finite update in drift-implicit Newton solver".into(),
        ));
    }
    Ok(x)
}

/// Solves `Y - base - scale * a(t_imp, Y) = 0` for `Y` via damped Newton-Raphson.
fn solve_implicit_drift<D>(
    drift: &D,
    t_imp: f64,
    base: &[f64],
    y_init: &[f64],
    scale: f64,
) -> Result<Vec<f64>, SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
{
    let n = base.len();
    if scale == 0.0 {
        return Ok(base.to_vec());
    }
    let mut y = y_init.to_vec();
    let mut a_val = vec![0.0; n];
    let tol = 64.0 * f64::EPSILON * (1.0 + base.iter().fold(0.0f64, |acc, x| acc.max(x.abs())));

    for _ in 0..25 {
        drift(t_imp, &y, &mut a_val);
        if !a_val.iter().all(|v| v.is_finite()) {
            return Err(SdeExtError::NonFinite(
                "drift became non-finite in implicit solver".into(),
            ));
        }
        let res: Vec<f64> = (0..n).map(|i| y[i] - base[i] - scale * a_val[i]).collect();
        let norm = res.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        if norm <= tol {
            return Ok(y);
        }
        let mut jac = vec![vec![0.0; n]; n];
        let mut a_pert = vec![0.0; n];
        for col in 0..n {
            let eps = f64::EPSILON.sqrt() * y[col].abs().max(1.0);
            let mut y_p = y.clone();
            y_p[col] += eps;
            drift(t_imp, &y_p, &mut a_pert);
            if !a_pert.iter().all(|v| v.is_finite()) {
                return Err(SdeExtError::NonFinite(
                    "drift Jacobian became non-finite".into(),
                ));
            }
            for row in 0..n {
                let da = (a_pert[row] - a_val[row]) / eps;
                jac[row][col] = if row == col {
                    1.0 - scale * da
                } else {
                    -scale * da
                };
            }
        }
        let neg_res: Vec<f64> = res.iter().map(|v| -v).collect();
        let delta = solve_linear_sde(jac, neg_res)?;
        let mut alpha = 1.0;
        let mut accepted = false;
        for _ in 0..10 {
            let y_trial: Vec<f64> = (0..n).map(|i| y[i] + alpha * delta[i]).collect();
            drift(t_imp, &y_trial, &mut a_pert);
            if a_pert.iter().all(|v| v.is_finite()) {
                let trial_norm = (0..n)
                    .map(|i| (y_trial[i] - base[i] - scale * a_pert[i]).abs())
                    .fold(0.0f64, f64::max);
                if trial_norm < norm {
                    y = y_trial;
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if !accepted {
            for i in 0..n {
                y[i] += delta[i];
            }
        }
    }
    drift(t_imp, &y, &mut a_val);
    let final_norm = (0..n)
        .map(|i| (y[i] - base[i] - scale * a_val[i]).abs())
        .fold(0.0f64, f64::max);
    if final_norm.is_finite() && final_norm <= tol {
        Ok(y)
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
    if noise.checked_mul(noise) != Some(areas.len()) || areas.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::invalid(
            "invalid implicit Milstein iterated-integral matrix",
        ));
    }
    implicit_step(drift, diffusion, noise, t, h, y, dw, Some(areas), scheme)
}

fn implicit_step<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t: f64,
    h: f64,
    y: &[f64],
    dw: &[f64],
    areas: Option<&[f64]>,
    scheme: DriftImplicitScheme,
) -> Result<Vec<f64>, SdeExtError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    let n = y.len();
    if n == 0
        || noise == 0
        || dw.len() != noise
        || !h.is_finite()
        || h <= 0.0
        || !t.is_finite()
        || !(t + h).is_finite()
        || t + h <= t
        || y.iter().chain(dw).any(|v| !v.is_finite())
        || n.checked_mul(noise).is_none_or(|size| size > 8_388_608)
        || n.checked_mul(n).is_none_or(|size| size > 8_388_608)
    {
        return Err(SdeExtError::invalid(
            "invalid or oversized drift-implicit step",
        ));
    }
    let theta = match scheme {
        | DriftImplicitScheme::Theta(theta) => theta,
        | DriftImplicitScheme::Trapezoidal => 0.5,
        | _ => 1.0,
    };
    if !(0.0..=1.0).contains(&theta) {
        return Err(SdeExtError::invalid("theta must be in [0,1]"));
    }
    let mut base = y.to_vec();
    let mut a = vec![0.0; n];
    let mut g = vec![0.0; n * noise];
    let split = scheme == DriftImplicitScheme::SplitStepBackwardEuler;
    if split {
        base = solve_implicit_drift(drift, t + h, y, y, h)?;
    } else if theta < 1.0 {
        drift(t, y, &mut a);
        if a.iter().any(|v| !v.is_finite()) {
            return Err(SdeExtError::NonFinite("non-finite explicit drift".into()));
        }
        for i in 0..n {
            base[i] += (1.0 - theta) * h * a[i];
        }
    }
    diffusion(t, if split { &base } else { y }, &mut g);
    if g.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::NonFinite(
            "non-finite implicit diffusion".into(),
        ));
    }
    for i in 0..n {
        for j in 0..noise {
            base[i] += g[i * noise + j] * dw[j];
        }
    }
    if scheme == DriftImplicitScheme::DriftImplicitMilstein {
        let correction = if let Some(areas) = areas {
            crate::kernels::sde::milstein_term_areas(diffusion, t, y, &g, noise, areas)
        } else {
            crate::kernels::sde::milstein_term(diffusion, t, y, &g, noise, dw, h)
        };
        for i in 0..n {
            base[i] += correction[i];
        }
    }
    if base.iter().any(|v| !v.is_finite()) {
        return Err(SdeExtError::NonFinite("non-finite implicit base".into()));
    }
    if split {
        Ok(base)
    } else {
        solve_implicit_drift(drift, t + h, &base, &base, theta * h)
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
    if y0.is_empty()
        || noise == 0
        || !t0.is_finite()
        || !dt.is_finite()
        || dt <= 0.0
        || y0.iter().any(|v| !v.is_finite())
        || increments
            .iter()
            .any(|dw| dw.len() != noise || dw.iter().any(|v| !v.is_finite()))
    {
        return Err(SdeExtError::invalid(
            "invalid drift-implicit initial data, grid, or increments",
        ));
    }
    let mut states = vec![y0.to_vec()];
    for (index, dw) in increments.iter().enumerate() {
        let t = t0 + index as f64 * dt;
        let state = implicit_step(
            drift,
            diffusion,
            noise,
            t,
            dt,
            states
                .last()
                .ok_or_else(|| SdeExtError::invalid("empty implicit trajectory"))?,
            dw,
            None,
            scheme,
        )?;
        states.push(state);
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
        return Ok(location + 2.0f64.sqrt() * scale * next_normal());
    }
    let u = next_unit().clamp(1e-15, 1.0 - 1e-15);
    let v = std::f64::consts::PI * (u - 0.5);
    let w = -next_unit().max(f64::MIN_POSITIVE).ln();

    let half_pi = std::f64::consts::FRAC_PI_2;
    if crate::kernels::exact_eq(alpha, 1.0) {
        let term1 = (half_pi + beta * v) * v.tan();
        let term2 = beta * ((half_pi * w * v.cos()) / (half_pi + beta * v)).ln();
        let x = (term1 - term2) / half_pi;
        Ok(scale * x + location + (2.0 / std::f64::consts::PI) * beta * scale * scale.ln())
    } else {
        let b_tan = beta * (alpha * half_pi).tan();
        let theta0 = b_tan.atan() / alpha;
        let factor = (1.0 + b_tan * b_tan).powf(1.0 / (2.0 * alpha));
        let num = (alpha * (v + theta0)).sin();
        let denom = v.cos().powf(1.0 / alpha);
        let cos_part = ((v - alpha * (v + theta0)).cos() / w).powf((1.0 - alpha) / alpha);
        let x = factor * (num / denom) * cos_part;
        Ok(scale * x + location)
    }
}
