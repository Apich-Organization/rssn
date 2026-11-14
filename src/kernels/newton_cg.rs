//! Matrix-free trust-region Newton-CG (Steihaug truncation).
//!
//! Only gradient and Hessian-vector actions are required; no dense Hessian,
//! inverse or factorization is formed. Truncated conjugate gradients follow
//! negative curvature to the trust boundary. The outer ratio test rejects
//! poor models and non-finite trial objectives, shrinking the region.
//!
//! # Limitations
//!
//! Smooth, unconstrained objectives only. Hessian actions must be linear and
//! symmetric and match the supplied gradient. Convergence is first-order
//! stationarity, so a zero gradient at a saddle is not a minimum certificate.
//! Floating-point intermediate products must be representable. There is no
//! preconditioner or constrained step. Opaque tape calls without second-order
//! rules cannot supply a Hessian action.

use crate::backend::tape::Tape;

/// Invalid configuration or unusable derivative information.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonCgError {
    /// Invalid point, radius, tolerance, or iteration configuration.
    InvalidInput,
    /// Initial objective or gradient is non-finite.
    NonFiniteInitial,
    /// Hessian action or conjugate-gradient arithmetic is non-finite.
    NonFiniteModel,
}

impl std::fmt::Display for NewtonCgError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::InvalidInput => "invalid Newton-CG input or options",
            | Self::NonFiniteInitial => "initial objective or gradient is not finite",
            | Self::NonFiniteModel => "Hessian action or Newton-CG arithmetic is not finite",
        })
    }
}
impl std::error::Error for NewtonCgError {}

/// Termination of matrix-free minimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonCgStatus {
    /// Gradient norm reached tolerance; this is not a second-order certificate.
    FirstOrderConverged,
    /// Outer iteration budget exhausted.
    IterationLimit,
    /// Rejected steps or an ineffective model shrank the radius below its limit.
    RadiusTooSmall,
}

/// Configuration of [`minimize`].
#[derive(Debug, Clone, Copy)]
pub struct NewtonCgOptions {
    /// Absolute stopping tolerance for the Euclidean gradient norm.
    pub gradient_tolerance: f64,
    /// Maximum outer iterations (zero returns the initial point).
    pub max_iterations: usize,
    /// Maximum Hessian actions in each truncated CG solve, excluding model check.
    pub max_cg_iterations: usize,
    /// Relative CG residual tolerance, tightened near stationarity.
    pub cg_relative_tolerance: f64,
    /// Initial trust-region radius.
    pub initial_radius: f64,
    /// Upper bound for the trust-region radius.
    pub max_radius: f64,
    /// Stop if the region must be shrunk below this positive radius.
    pub min_radius: f64,
    /// Accept actual/predicted decrease ratios above this threshold in `[0, .25)`.
    pub acceptance_threshold: f64,
}

impl Default for NewtonCgOptions {
    fn default() -> Self {
        Self {
            gradient_tolerance: 1e-8,
            max_iterations: 200,
            max_cg_iterations: 100,
            cg_relative_tolerance: 0.1,
            initial_radius: 1.0,
            max_radius: 1000.0,
            min_radius: 1e-14,
            acceptance_threshold: 0.1,
        }
    }
}

/// The last accepted iterate and explicit convergence/work counters.
#[derive(Debug, Clone)]
pub struct NewtonCgResult {
    /// Last accepted point; rejected trial points never replace it.
    pub x: Vec<f64>,
    /// Objective at `x`.
    pub value: f64,
    /// Euclidean gradient norm at `x`.
    pub gradient_norm: f64,
    /// Reason for stopping.
    pub status: NewtonCgStatus,
    /// Completed outer iterations.
    pub iterations: usize,
    /// Objective evaluations, including rejected and non-finite trials.
    pub objective_evaluations: usize,
    /// Gradient evaluations, including rejected trials with invalid gradients.
    pub gradient_evaluations: usize,
    /// Hessian-vector evaluations, including explicit predicted-decrease checks.
    pub hessian_vector_evaluations: usize,
    /// CG solves stopped on non-positive curvature.
    pub negative_curvature_steps: usize,
    /// CG solves truncated at the trust boundary.
    pub boundary_steps: usize,
    /// Number of accepted outer steps.
    pub accepted_steps: usize,
}

/// Minimize a smooth objective using only vector actions.
///
/// `gradient(x, out)` and `hessian_vector(x, direction, out)` must fill every
/// output element. Buffers are reset to NaN before each callback, so missing
/// writes are detected. Non-finite trial objectives/gradients are rejected;
/// the last accepted point is retained. No parallel callback calls occur.
///
/// # Errors
/// Returns an error for invalid options/input, a non-finite initial value or
/// gradient, or a non-finite Hessian action/CG intermediate.
pub fn minimize<F, G, H>(
    mut objective: F,
    mut gradient: G,
    mut hessian_vector: H,
    initial: &[f64],
    options: &NewtonCgOptions,
) -> Result<NewtonCgResult, NewtonCgError>
where
    F: FnMut(&[f64]) -> f64,
    G: FnMut(&[f64], &mut [f64]),
    H: FnMut(&[f64], &[f64], &mut [f64]),
{
    validate_options(initial, options)?;
    let mut x = initial.to_vec();
    let mut value = objective(&x);
    let mut g = vec![f64::NAN; x.len()];
    gradient(&x, &mut g);
    if !value.is_finite() || !finite(&g) || !length(&g).is_finite() {
        return Err(NewtonCgError::NonFiniteInitial);
    }
    let mut result = NewtonCgResult {
        x: Vec::new(),
        value,
        gradient_norm: length(&g),
        status: NewtonCgStatus::IterationLimit,
        iterations: 0,
        objective_evaluations: 1,
        gradient_evaluations: 1,
        hessian_vector_evaluations: 0,
        negative_curvature_steps: 0,
        boundary_steps: 0,
        accepted_steps: 0,
    };
    let mut radius = options.initial_radius;
    let mut model_action = vec![f64::NAN; x.len()];
    let mut trial_gradient = vec![f64::NAN; x.len()];
    for iteration in 0..options.max_iterations {
        if length(&g) <= options.gradient_tolerance {
            break;
        }
        let step = truncated_cg(&x, &g, radius, options, &mut hessian_vector)?;
        result.iterations = iteration + 1;
        result.hessian_vector_evaluations += step.actions;
        result.negative_curvature_steps += usize::from(step.negative_curvature);
        result.boundary_steps += usize::from(step.boundary);
        model_action.fill(f64::NAN);
        hessian_vector(&x, &step.vector, &mut model_action);
        result.hessian_vector_evaluations += 1;
        if !finite(&model_action) {
            return Err(NewtonCgError::NonFiniteModel);
        }
        let predicted = -dot(&g, &step.vector) - 0.5 * dot(&step.vector, &model_action);
        let trial: Vec<_> = x.iter().zip(&step.vector).map(|(x, p)| x + p).collect();
        let mut ratio = f64::NEG_INFINITY;
        let mut trial_value = f64::NAN;
        if predicted.is_finite() && predicted > 0.0 && finite(&trial) {
            trial_value = objective(&trial);
            result.objective_evaluations += 1;
            if trial_value.is_finite() {
                ratio = (value - trial_value) / predicted;
            }
        }
        if ratio > options.acceptance_threshold {
            trial_gradient.fill(f64::NAN);
            gradient(&trial, &mut trial_gradient);
            result.gradient_evaluations += 1;
            if finite(&trial_gradient) && length(&trial_gradient).is_finite() {
                x = trial;
                value = trial_value;
                g.clone_from(&trial_gradient);
                result.accepted_steps += 1;
            } else {
                ratio = f64::NEG_INFINITY;
            }
        }
        if ratio < 0.25 || ratio.is_nan() {
            radius *= 0.25;
        } else if ratio > 0.75 && step.boundary {
            radius = (radius * 2.0).min(options.max_radius);
        }
        if radius < options.min_radius {
            result.status = NewtonCgStatus::RadiusTooSmall;
            break;
        }
    }
    result.gradient_norm = length(&g);
    if result.gradient_norm <= options.gradient_tolerance {
        result.status = NewtonCgStatus::FirstOrderConverged;
    }
    result.x = x;
    result.value = value;
    Ok(result)
}

/// Minimize one tape output using analytic reverse gradients and HVPs.
///
/// Parameters are constant. The objective tape must have matching input arity.
/// No finite-difference Hessian is substituted when a second-order rule is absent.
///
/// # Errors
/// Returns the same errors as [`minimize`], or `InvalidInput` for an invalid
/// output index/arity. Opaque second-order rules report `NonFiniteModel`.
pub fn minimize_tape(
    tape: &Tape,
    output: usize,
    initial: &[f64],
    params: &[f64],
    options: &NewtonCgOptions,
) -> Result<NewtonCgResult, NewtonCgError> {
    if tape.inputs != initial.len() || output >= tape.outputs.len() || tape.params != params.len() {
        return Err(NewtonCgError::InvalidInput);
    }
    let mut values = vec![f64::NAN; tape.outputs.len()];
    minimize(
        |x| {
            tape.eval(x, params, &mut values);
            values[output]
        },
        |x, out| out.copy_from_slice(&tape.gradient_reverse(output, x, params)),
        |x, direction, out| {
            out.copy_from_slice(&tape.hessian_vector_product(output, x, params, direction));
        },
        initial,
        options,
    )
}

fn validate_options(
    initial: &[f64],
    o: &NewtonCgOptions,
) -> Result<(), NewtonCgError> {
    if !finite(initial)
        || !o.gradient_tolerance.is_finite()
        || o.gradient_tolerance < 0.0
        || o.max_cg_iterations == 0
        || !o.cg_relative_tolerance.is_finite()
        || o.cg_relative_tolerance <= 0.0
        || o.cg_relative_tolerance >= 1.0
        || !o.min_radius.is_finite()
        || o.min_radius <= 0.0
        || !o.initial_radius.is_finite()
        || o.initial_radius < o.min_radius
        || !o.max_radius.is_finite()
        || o.max_radius < o.initial_radius
        || !o.acceptance_threshold.is_finite()
        || !(0.0..0.25).contains(&o.acceptance_threshold)
    {
        Err(NewtonCgError::InvalidInput)
    } else {
        Ok(())
    }
}

struct Step {
    vector: Vec<f64>,
    actions: usize,
    boundary: bool,
    negative_curvature: bool,
}

fn truncated_cg<H>(
    x: &[f64],
    gradient: &[f64],
    radius: f64,
    o: &NewtonCgOptions,
    action: &mut H,
) -> Result<Step, NewtonCgError>
where
    H: FnMut(&[f64], &[f64], &mut [f64]),
{
    let mut z = vec![0.0; x.len()];
    let mut r = gradient.to_vec();
    let mut d: Vec<_> = r.iter().map(|r| -r).collect();
    let mut hd = vec![f64::NAN; x.len()];
    let tolerance = o.cg_relative_tolerance.min(length(gradient).sqrt()) * length(gradient);
    let mut actions = 0;
    let mut boundary = false;
    let mut negative_curvature = false;
    for _ in 0..o.max_cg_iterations {
        let rnorm = length(&r);
        if rnorm <= tolerance {
            break;
        }
        hd.fill(f64::NAN);
        action(x, &d, &mut hd);
        actions += 1;
        if !finite(&hd) {
            return Err(NewtonCgError::NonFiniteModel);
        }
        let curvature = dot(&d, &hd);
        if !curvature.is_finite() {
            return Err(NewtonCgError::NonFiniteModel);
        }
        if curvature <= 0.0 {
            z = boundary_point(&z, &d, radius)?;
            boundary = true;
            negative_curvature = true;
            break;
        }
        let alpha = rnorm * (rnorm / curvature);
        let next: Vec<_> = z.iter().zip(&d).map(|(z, d)| z + alpha * d).collect();
        if !finite(&next) || length(&next) >= radius {
            z = boundary_point(&z, &d, radius)?;
            boundary = true;
            break;
        }
        z = next;
        for (r, h) in r.iter_mut().zip(&hd) {
            *r += alpha * h;
        }
        let next_norm = length(&r);
        if !next_norm.is_finite() {
            return Err(NewtonCgError::NonFiniteModel);
        }
        if next_norm <= tolerance {
            break;
        }
        let beta = (next_norm / rnorm).powi(2);
        for (d, r) in d.iter_mut().zip(&r) {
            *d = -r + beta * *d;
        }
        if !finite(&d) {
            return Err(NewtonCgError::NonFiniteModel);
        }
    }
    Ok(Step {
        vector: z,
        actions,
        boundary,
        negative_curvature,
    })
}

fn boundary_point(
    z: &[f64],
    d: &[f64],
    radius: f64,
) -> Result<Vec<f64>, NewtonCgError> {
    let dnorm = length(d);
    if !dnorm.is_finite() || dnorm == 0.0 {
        return Err(NewtonCgError::NonFiniteModel);
    }
    let u: Vec<_> = z.iter().map(|z| z / radius).collect();
    let w: Vec<_> = d.iter().map(|d| d / dnorm).collect();
    let along = dot(&u, &w);
    let slack = (1.0 - dot(&u, &u)).max(0.0);
    let root = along.hypot(slack.sqrt());
    let distance = if root + along == 0.0 {
        0.0
    } else if along >= 0.0 {
        slack / (root + along)
    } else {
        root - along
    };
    let result: Vec<_> = z
        .iter()
        .zip(w)
        .map(|(z, w)| z + (radius * distance) * w)
        .collect();
    if finite(&result) {
        Ok(result)
    } else {
        Err(NewtonCgError::NonFiniteModel)
    }
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn length(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |norm, x| norm.hypot(*x))
}
fn finite(x: &[f64]) -> bool {
    x.iter().all(|x| x.is_finite())
}
