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
//! Hessian actions must fit in floating point, as must the accepted gradient
//! and the model decrease used by the trust-region test. There is no
//! preconditioner or constrained step. Opaque tape calls without second-order
//! rules cannot supply a Hessian action.

use super::numeric_common::scale_product;
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
    objective: F,
    gradient: G,
    hessian_vector: H,
    initial: &[f64],
    options: &NewtonCgOptions,
) -> Result<NewtonCgResult, NewtonCgError>
where
    F: FnMut(&[f64]) -> f64,
    G: FnMut(&[f64], &mut [f64]),
    H: FnMut(&[f64], &[f64], &mut [f64]),
{
    NewtonCgWorkspace::default().minimize(objective, gradient, hessian_vector, initial, options)
}

/// Reusable scratch storage for repeated matrix-free Newton–CG solves.
///
/// Reuse one workspace for related minimizations, including problems with
/// different dimensions. Each call starts from its own initial point. Values
/// returned by one solve are owned and independent of subsequent calls.
#[derive(Debug, Default)]
pub struct NewtonCgWorkspace {
    cg: CgWorkspace,
    g: Vec<f64>,
    trial: Vec<f64>,
    trial_gradient: Vec<f64>,
    model_action: Vec<f64>,
}

impl NewtonCgWorkspace {
    /// Minimize with the callback contract of [`minimize`], reusing buffers.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`minimize`].
    pub fn minimize<F, G, H>(
        &mut self,
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
        let Self {
            cg,
            g,
            trial,
            trial_gradient,
            model_action,
        } = self;
        g.resize(x.len(), f64::NAN);
        g.fill(f64::NAN);
        gradient(&x, g);
        let mut gradient_norm = length(g);
        if !value.is_finite() || !gradient_norm.is_finite() {
            return Err(NewtonCgError::NonFiniteInitial);
        }
        let mut result = NewtonCgResult {
            x: Vec::new(),
            value,
            gradient_norm,
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
        if options.max_iterations > 0 && gradient_norm > options.gradient_tolerance {
            trial.resize(x.len(), 0.0);
            model_action.resize(x.len(), f64::NAN);
            trial_gradient.resize(x.len(), f64::NAN);
        }
        for iteration in 0..options.max_iterations {
            if gradient_norm <= options.gradient_tolerance {
                break;
            }
            let step = cg.solve(&x, g, gradient_norm, radius, options, &mut hessian_vector)?;
            result.iterations = iteration + 1;
            result.hessian_vector_evaluations += step.actions;
            result.negative_curvature_steps += usize::from(step.negative_curvature);
            result.boundary_steps += usize::from(step.boundary);
            model_action.fill(f64::NAN);
            hessian_vector(&x, &cg.z, model_action);
            result.hessian_vector_evaluations += 1;
            if !finite(model_action) {
                return Err(NewtonCgError::NonFiniteModel);
            }
            let predicted = -dot(g, &cg.z) - 0.5 * dot(&cg.z, model_action);
            for ((trial, x), p) in trial.iter_mut().zip(&x).zip(&cg.z) {
                *trial = x + p;
            }
            let mut ratio = f64::NEG_INFINITY;
            let mut trial_value = f64::NAN;
            if predicted.is_finite() && predicted > 0.0 && finite(trial) {
                trial_value = objective(trial);
                result.objective_evaluations += 1;
                if trial_value.is_finite() {
                    ratio = (value - trial_value) / predicted;
                }
            }
            if ratio > options.acceptance_threshold {
                trial_gradient.fill(f64::NAN);
                gradient(trial, trial_gradient);
                result.gradient_evaluations += 1;
                let trial_norm = length(trial_gradient);
                if trial_norm.is_finite() {
                    std::mem::swap(&mut x, trial);
                    std::mem::swap(g, trial_gradient);
                    value = trial_value;
                    gradient_norm = trial_norm;
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
        result.gradient_norm = gradient_norm;
        if result.gradient_norm <= options.gradient_tolerance {
            result.status = NewtonCgStatus::FirstOrderConverged;
        }
        result.x = x;
        result.value = value;
        Ok(result)
    }
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
    NewtonCgWorkspace::default().minimize_tape(tape, output, initial, params, options)
}

impl NewtonCgWorkspace {
    /// Minimize a tape output with reusable Newton and CG buffers.
    ///
    /// # Errors
    /// Returns the same errors as [`minimize_tape`].
    pub fn minimize_tape(
        &mut self,
        tape: &Tape,
        output: usize,
        initial: &[f64],
        params: &[f64],
        options: &NewtonCgOptions,
    ) -> Result<NewtonCgResult, NewtonCgError> {
        if tape.inputs != initial.len()
            || output >= tape.outputs.len()
            || tape.params != params.len()
        {
            return Err(NewtonCgError::InvalidInput);
        }
        let mut values = vec![f64::NAN; tape.outputs.len()];
        self.minimize(
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
    actions: usize,
    boundary: bool,
    negative_curvature: bool,
}

#[derive(Debug, Default)]
struct CgWorkspace {
    z: Vec<f64>,
    r: Vec<f64>,
    d: Vec<f64>,
    hd: Vec<f64>,
    next: Vec<f64>,
}

impl CgWorkspace {
    fn solve<H>(
        &mut self,
        x: &[f64],
        gradient: &[f64],
        gradient_norm: f64,
        radius: f64,
        o: &NewtonCgOptions,
        action: &mut H,
    ) -> Result<Step, NewtonCgError>
    where
        H: FnMut(&[f64], &[f64], &mut [f64]),
    {
        let Self { z, r, d, hd, next } = self;
        z.resize(x.len(), 0.0);
        z.fill(0.0);
        r.resize(x.len(), 0.0);
        r.copy_from_slice(gradient);
        d.resize(x.len(), 0.0);
        for (d, r) in d.iter_mut().zip(&*r) {
            *d = -r;
        }
        hd.resize(x.len(), f64::NAN);
        next.resize(x.len(), 0.0);
        let tolerance = o.cg_relative_tolerance.min(gradient_norm.sqrt()) * gradient_norm;
        let mut rnorm = gradient_norm;
        let mut step = Step {
            actions: 0,
            boundary: false,
            negative_curvature: false,
        };
        for _ in 0..o.max_cg_iterations {
            if rnorm <= tolerance {
                break;
            }
            hd.fill(f64::NAN);
            action(x, d, hd);
            step.actions += 1;
            if !finite(hd) {
                return Err(NewtonCgError::NonFiniteModel);
            }
            let Some(scale) = step_scale(rnorm, d, hd) else {
                boundary_point(z, d, radius, hd, next)?;
                step.boundary = true;
                step.negative_curvature = true;
                break;
            };
            scale.step.assign_to(next, z, d);
            if !finite(next) || length(next) >= radius {
                boundary_point(z, d, radius, hd, next)?;
                step.boundary = true;
                break;
            }
            std::mem::swap(z, next);
            scale.residual.add_to(r, hd);
            let next_norm = length(r);
            if !next_norm.is_finite() {
                return Err(NewtonCgError::NonFiniteModel);
            }
            if next_norm <= tolerance {
                break;
            }
            let beta = (next_norm / rnorm).powi(2);
            rnorm = next_norm;
            for (d, r) in d.iter_mut().zip(&*r) {
                *d = -r + beta * *d;
            }
            if !finite(d) {
                return Err(NewtonCgError::NonFiniteModel);
            }
        }
        Ok(step)
    }
}

#[derive(Clone, Copy)]
struct Update {
    factor: f64,
    divisor: f64,
}

impl Update {
    fn assign_to(
        self,
        target: &mut [f64],
        base: &[f64],
        source: &[f64],
    ) {
        if super::exact_eq(self.divisor, 1.0) {
            for ((target, base), source) in target.iter_mut().zip(base).zip(source) {
                *target = base + self.factor * source;
            }
        } else {
            for ((target, base), source) in target.iter_mut().zip(base).zip(source) {
                *target = base + self.factor * (source / self.divisor);
            }
        }
    }

    fn add_to(
        self,
        target: &mut [f64],
        source: &[f64],
    ) {
        if super::exact_eq(self.divisor, 1.0) {
            for (target, source) in target.iter_mut().zip(source) {
                *target += self.factor * source;
            }
        } else {
            for (target, source) in target.iter_mut().zip(source) {
                *target += self.factor * (source / self.divisor);
            }
        }
    }
}

struct StepScale {
    step: Update,
    residual: Update,
}

// Fall back before underflow can be mistaken for non-positive curvature.
// Apply separate scale factors to the direction and Hessian action: alpha
// need not fit in a float when the actual vector increments do.
fn step_scale(
    rnorm: f64,
    d: &[f64],
    hd: &[f64],
) -> Option<StepScale> {
    let curvature = dot(d, hd);
    if curvature.is_normal() {
        if curvature < 0.0 {
            return None;
        }
        let alpha = rnorm * (rnorm / curvature);
        if alpha.is_normal() {
            let update = Update {
                factor: alpha,
                divisor: 1.0,
            };
            return Some(StepScale {
                step: update,
                residual: update,
            });
        }
    }

    let d_scale = d.iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
    let h_scale = hd.iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
    if d_scale == 0.0 || h_scale == 0.0 {
        return None;
    }
    let scaled: f64 = d
        .iter()
        .zip(hd)
        .map(|(d, h)| (d / d_scale) * (h / h_scale))
        .sum();
    let step = |divisor, other_scale| {
        let factor = if scaled > 0.0 {
            scale_product(rnorm, &[rnorm], &[other_scale, scaled])
        } else {
            scale_product(rnorm, &[rnorm, divisor], &[curvature])
        };
        Update {
            factor: factor.unwrap_or(f64::INFINITY),
            divisor,
        }
    };
    // A representable raw product can resolve the sign in a mixed-scale
    // direction whose doubly normalized terms have all rounded to zero.
    if scaled > 0.0 || (scaled == 0.0 && curvature.is_finite() && curvature > 0.0) {
        Some(StepScale {
            step: step(d_scale, h_scale),
            residual: step(h_scale, d_scale),
        })
    } else {
        None
    }
}

fn boundary_point(
    z: &mut [f64],
    d: &[f64],
    radius: f64,
    w: &mut [f64],
    u: &mut [f64],
) -> Result<(), NewtonCgError> {
    let mut dnorm = length(d);
    let mut scale = 1.0;
    if !dnorm.is_normal() {
        scale = d.iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
        if scale == 0.0 {
            return Err(NewtonCgError::NonFiniteModel);
        }
        dnorm = d.iter().fold(0.0_f64, |norm, x| norm.hypot(x / scale));
    }
    if !dnorm.is_finite() || dnorm == 0.0 {
        return Err(NewtonCgError::NonFiniteModel);
    }
    for ((u, w), (z, d)) in u.iter_mut().zip(&mut *w).zip(z.iter().zip(d)) {
        *u = z / radius;
        *w = (d / scale) / dnorm;
    }
    let along = dot(u, w);
    let slack = (1.0 - dot(u, u)).max(0.0);
    let root = along.hypot(slack.sqrt());
    let distance = if along < 0.0 {
        root - along
    } else if root + along > 0.0 {
        slack / (root + along)
    } else {
        0.0
    };
    let extent = radius * distance;
    if extent.is_finite() {
        crate::kernels::simd::axpy_f64(extent, w, z);
    } else {
        for ((z, u), w) in z.iter_mut().zip(u).zip(w) {
            *z = radius * (*u + distance * *w);
        }
    }
    if finite(z) {
        Ok(())
    } else {
        Err(NewtonCgError::NonFiniteModel)
    }
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    crate::kernels::simd::dot_f64(a, b)
}
fn length(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |norm, x| norm.hypot(*x))
}
fn finite(x: &[f64]) -> bool {
    x.iter().all(|x| x.is_finite())
}

#[cfg(test)]
mod tests {
    use super::boundary_point;

    #[test]
    fn inward_ray_from_rounded_boundary_reaches_the_far_root() {
        let radius = 0.75 * f64::MAX;
        let mut z = [-radius];
        boundary_point(&mut z, &[1.0], radius, &mut [0.0], &mut [0.0]).unwrap();
        assert_eq!(z[0].to_bits(), radius.to_bits());
    }

    #[test]
    fn subnormal_boundary_direction_is_normalized_in_scaled_coordinates() {
        let mut z = [0.0; 2];
        boundary_point(
            &mut z,
            &[f64::from_bits(1); 2],
            1.0,
            &mut [0.0; 2],
            &mut [0.0; 2],
        )
        .unwrap();
        assert!(z.iter().all(|coordinate| {
            (coordinate - std::f64::consts::FRAC_1_SQRT_2).abs() <= f64::EPSILON
        }));
    }
}
