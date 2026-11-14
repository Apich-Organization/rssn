//! Scaled, damped nonlinear least squares without normal equations.
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

/// Termination reason for a nonlinear least-squares solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LeastSquaresStatus {
    /// The residual norm met the requested absolute tolerance.
    ResidualConverged,
    /// The scaled first-order gradient met the requested tolerance.
    GradientConverged,
    /// The iteration limit was reached.
    IterationLimit,
    /// The residual-evaluation budget was reached.
    EvaluationLimit,
    /// No useful representable step could be accepted.
    Stagnated,
}
impl LeastSquaresStatus {
    /// Whether an explicit residual or first-order convergence test passed.
    #[must_use]
    pub const fn converged(self) -> bool {
        matches!(self, Self::ResidualConverged | Self::GradientConverged)
    }
}
/// Controls for the damped least-squares algorithm.
#[derive(Clone, Copy, Debug)]
pub struct LeastSquaresOptions {
    /// Maximum accepted or rejected outer iterations.
    pub max_iterations: usize,
    /// Maximum calls to the residual function, including finite differences.
    pub max_evaluations: usize,
    /// Absolute Euclidean residual tolerance.
    pub residual_tolerance: f64,
    /// Dimensionless tolerance on the column-scaled gradient.
    pub gradient_tolerance: f64,
    /// Small relative parameter steps trigger stagnation unless convergence passed.
    pub step_tolerance: f64,
    /// Initial damping parameter, restricted to `1e-300..=1e16`.
    pub initial_damping: f64,
}
impl Default for LeastSquaresOptions {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            max_evaluations: 10_000,
            residual_tolerance: 1e-10,
            gradient_tolerance: 1e-10,
            step_tolerance: 1e-12,
            initial_damping: 1e-3,
        }
    }
}
/// Invalid input or a malformed callback result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LeastSquaresError {
    /// Options or initial parameters are invalid.
    InvalidInput,
    /// Residual vector is empty, changes length, or is non-finite initially.
    InvalidResidual,
    /// A Jacobian has the wrong shape or contains non-finite entries.
    InvalidJacobian,
    /// A stable SVD could not be obtained.
    Decomposition,
}
/// Result and work accounting. A nonconverged result retains the best point.
#[derive(Clone, Debug, PartialEq)]
pub struct LeastSquaresResult {
    /// Best accepted parameters.
    pub x: Vec<f64>,
    /// Euclidean residual norm (infinity if the norm is not representable).
    pub residual_norm: f64,
    /// Number of attempted outer iterations.
    pub iterations: usize,
    /// Number of residual function evaluations.
    pub evaluations: usize,
    /// Number of Jacobian evaluations, including finite difference sweeps.
    pub jacobian_evaluations: usize,
    /// Number of accepted trial steps.
    pub accepted_steps: usize,
    /// Number of rejected trial steps.
    pub rejected_steps: usize,
    /// Reason for termination.
    pub status: LeastSquaresStatus,
}
impl LeastSquaresResult {
    /// Whether the reported termination passed a convergence test.
    #[must_use]
    pub const fn converged(&self) -> bool {
        self.status.converged()
    }
}
fn norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |acc, &x| acc.hypot(x))
}
fn finite(v: &[f64]) -> bool {
    v.iter().all(|v| v.is_finite())
}
fn report(
    x: Vec<f64>,
    residual: &[f64],
    iterations: usize,
    work: (usize, usize, usize, usize),
    status: LeastSquaresStatus,
) -> LeastSquaresResult {
    let (evaluations, jacobians, accepted, rejected) = work;
    LeastSquaresResult {
        x,
        residual_norm: norm(residual),
        iterations,
        evaluations,
        jacobian_evaluations: jacobians,
        accepted_steps: accepted,
        rejected_steps: rejected,
        status,
    }
}
// Compare a non-underflowing normalized magnitude against a finite tolerance.
// A positive exact gradient smaller than the minimum subnormal still fails a
// zero tolerance rather than being rounded into a fictitious stationary point.
fn magnitude_within(
    fraction: f64,
    exponent: i32,
    tolerance: f64,
) -> bool {
    if fraction == 0.0 {
        return true;
    }
    if tolerance == 0.0 {
        return false;
    }
    let (fraction, exponent) = if fraction.abs() >= 2.0 {
        (fraction * 0.5, exponent + 1)
    } else {
        (fraction, exponent)
    };
    let (normal, adjustment) = if tolerance < f64::MIN_POSITIVE {
        (tolerance * 18014398509481984.0, -54)
    } else {
        (tolerance, 0)
    };
    let bits = normal.to_bits();
    let tol_exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + adjustment;
    let tol_fraction = f64::from_bits((bits & ((1u64 << 52) - 1)) | (1023u64 << 52));
    exponent < tol_exponent || (exponent == tol_exponent && fraction.abs() <= tol_fraction)
}
// Divide an unbounded-exponent product sum by a positive finite singular
// value before materializing it as f64. This keeps a tiny projection that can
// become representable again after division by a small singular value.
fn scaled_quotient(
    fraction: f64,
    exponent: i32,
    denominator: f64,
) -> Option<f64> {
    if fraction == 0.0 {
        return Some(0.0);
    }
    let (normal, adjust) = if denominator < f64::MIN_POSITIVE {
        (denominator * 18014398509481984.0, -54)
    } else {
        (denominator, 0)
    };
    let bits = normal.to_bits();
    let divisor_exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + adjust;
    let divisor_fraction = f64::from_bits((bits & ((1u64 << 52) - 1)) | (1023u64 << 52));
    let mut value = fraction / divisor_fraction;
    let mut shift = exponent - divisor_exponent;
    while shift > 512 {
        value *= 2.0_f64.powi(512);
        shift -= 512;
        if !value.is_finite() {
            return None;
        }
    }
    while shift < -512 {
        value *= 2.0_f64.powi(-512);
        shift += 512;
        if value == 0.0 {
            return Some(0.0);
        }
    }
    let value = value * 2.0_f64.powi(shift);
    value.is_finite().then_some(value)
}
// Multiply/divide in a safe order; if no finite ordering exists the candidate
// is rejected rather than claiming convergence from an invalid step.
fn mul_ratio(
    v: f64,
    numerator: f64,
    denominator: f64,
) -> Option<f64> {
    let choices = [
        (v * numerator) / denominator,
        (v / denominator) * numerator,
        v * (numerator / denominator),
    ];
    choices
        .into_iter()
        .find(|x| x.is_finite() && (*x != 0.0 || v == 0.0))
        .or_else(|| (v == 0.0).then_some(0.0))
}
/// Analytic Jacobian callback for nonlinear least squares.
pub type LeastSquaresJacobian<'a> = dyn Fn(&[f64]) -> Mat + 'a;

/// Minimise one half of the squared residual norm using a damped, column-scaled
/// augmented SVD. The Jacobian may be analytic or safely differenced.
///
/// The reported scaled gradient is formed from `J[:,j]/max(abs(J[:,j]))` and
/// `f/max(abs(f))`; it is a dimensionless first-order test, not a root proof.
/// Trial steps are accepted only on positive actual and predicted decrease.
/// Residual length must remain constant on every evaluation.
///
/// # Errors
/// Invalid input, malformed callback results, or a failed decomposition return
/// a [`LeastSquaresError`]. Nonconvergence is reported as a result status.
pub fn least_squares<F>(
    f: F,
    jac: Option<&LeastSquaresJacobian<'_>>,
    x0: &[f64],
    options: LeastSquaresOptions,
) -> Result<LeastSquaresResult, LeastSquaresError>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let LeastSquaresOptions {
        max_iterations,
        max_evaluations,
        residual_tolerance: ftol,
        gradient_tolerance: gtol,
        step_tolerance: xtol,
        initial_damping,
    } = options;
    if x0.is_empty()
        || !finite(x0)
        || max_evaluations == 0
        || !ftol.is_finite()
        || ftol < 0.0
        || !gtol.is_finite()
        || gtol < 0.0
        || !xtol.is_finite()
        || xtol < 0.0
        || !initial_damping.is_finite()
        || !(1e-300..=1e16).contains(&initial_damping)
    {
        return Err(LeastSquaresError::InvalidInput);
    }
    let mut x = x0.to_vec();
    let mut residual = f(&x);
    let mut evaluations = 1;
    let mut jacobians = 0;
    let mut accepted = 0;
    let mut rejected = 0;
    if residual.is_empty() || !finite(&residual) {
        return Err(LeastSquaresError::InvalidResidual);
    }
    let m = residual.len();
    let n = x.len();
    if m.checked_add(n)
        .and_then(|rows| rows.checked_mul(n))
        .is_none()
    {
        return Err(LeastSquaresError::InvalidInput);
    }
    let mut lambda = initial_damping;
    let mut pending_small = false;
    macro_rules! finish {
        ($it:expr, $status:expr) => {
            return Ok(report(
                x,
                &residual,
                $it,
                (evaluations, jacobians, accepted, rejected),
                $status,
            ))
        };
    }
    for iteration in 0..=max_iterations {
        let scale = residual.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 || norm(&residual) <= ftol {
            finish!(iteration, LeastSquaresStatus::ResidualConverged);
        }
        if iteration == max_iterations {
            finish!(iteration, LeastSquaresStatus::IterationLimit);
        }
        jacobians += 1;
        let j = if let Some(jac) = jac {
            jac(&x)
        } else {
            let mut matrix = Mat::zeros(m, n);
            for col in 0..n {
                let h = f64::EPSILON.sqrt() * x[col].abs().max(1.0);
                let mut trial = x.clone();
                let mut sample = None;
                for sign in [1.0_f64, -1.0] {
                    trial[col] = sign.mul_add(h, x[col]);
                    let delta = trial[col] - x[col];
                    if !trial[col].is_finite() || delta == 0.0 || !delta.is_finite() {
                        continue;
                    }
                    if evaluations == max_evaluations {
                        finish!(iteration, LeastSquaresStatus::EvaluationLimit);
                    }
                    let value = f(&trial);
                    evaluations += 1;
                    if value.len() != m {
                        return Err(LeastSquaresError::InvalidResidual);
                    }
                    if finite(&value) {
                        sample = Some((value, delta));
                        break;
                    }
                }
                let Some((sample, delta)) = sample else {
                    return Err(LeastSquaresError::InvalidJacobian);
                };
                for row in 0..m {
                    matrix.set(row, col, (sample[row] - residual[row]) / delta);
                }
            }
            matrix
        };
        if j.rows != m || j.cols != n || j.data.len() != m * n || !finite(&j.data) {
            return Err(LeastSquaresError::InvalidJacobian);
        }
        let mut col_scale = vec![0.0_f64; n];
        for row in 0..m {
            for col in 0..n {
                col_scale[col] = col_scale[col].max(j.at(row, col).abs());
            }
        }
        let mut b = Mat::zeros(m, n);
        for col in 0..n {
            if col_scale[col] == 0.0 {
                col_scale[col] = 1.0;
            }
            for row in 0..m {
                b.set(row, col, j.at(row, col) / col_scale[col]);
            }
        }
        let r: Vec<_> = residual.iter().map(|v| v / scale).collect();
        let mut gradient_converged = true;
        for col in 0..n {
            let (fraction, exponent) = crate::kernels::krylov::exact_product_sum_scaled(
                (0..m).map(|row| (b.at(row, col), r[row])),
            );
            if !fraction.is_finite() {
                return Err(LeastSquaresError::InvalidJacobian);
            }
            if !magnitude_within(fraction, exponent, gtol) {
                gradient_converged = false;
            }
        }
        if gradient_converged {
            finish!(iteration, LeastSquaresStatus::GradientConverged);
        }
        if pending_small {
            finish!(iteration, LeastSquaresStatus::Stagnated);
        }
        let mut success = false;
        let mut last_step_small = false;
        let initial_trial_damping = lambda;
        for attempt in 0..64 {
            if attempt == 32 {
                lambda = (initial_trial_damping * 0.1).max(1e-300);
            }
            let next_damping = |value: f64| {
                if attempt < 32 {
                    (value * 4.0).min(1e16)
                } else {
                    (value * 0.1).max(1e-300)
                }
            };
            let mut augmented = Mat::zeros(m + n, n);
            for row in 0..m {
                for col in 0..n {
                    augmented.set(row, col, b.at(row, col));
                }
            }
            for col in 0..n {
                augmented.set(m + col, col, lambda.sqrt());
            }
            let svd = dense::svd(&augmented).map_err(|_| LeastSquaresError::Decomposition)?;
            let mut z = vec![0.0; n];
            for k in 0..svd.s.len() {
                if svd.s[k] <= 0.0 || !svd.s[k].is_finite() {
                    return Err(LeastSquaresError::Decomposition);
                }
                let (projection, exponent) = crate::kernels::krylov::exact_product_sum_scaled(
                    (0..m).map(|row| (svd.u.at(row, k), r[row])),
                );
                let Some(coefficient) = scaled_quotient(-projection, exponent, svd.s[k]) else {
                    return Err(LeastSquaresError::Decomposition);
                };
                for col in 0..n {
                    z[col] += svd.v.at(col, k) * coefficient;
                }
            }
            let mut trial = x.clone();
            let mut step = vec![0.0; n];
            let mut valid = true;
            for col in 0..n {
                if let Some(value) = mul_ratio(z[col], scale, col_scale[col]) {
                    step[col] = value;
                    trial[col] = x[col] + value;
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid || !finite(&trial) || !finite(&step) {
                lambda = next_damping(lambda);
                rejected += 1;
                continue;
            }
            let step_scale = step.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            let parameter_scale = x.iter().map(|v| v.abs()).fold(1.0_f64, f64::max);
            last_step_small = step_scale / parameter_scale <= xtol;
            let predicted: f64 = (0..m)
                .map(|row| {
                    let model = r[row] + (0..n).map(|col| b.at(row, col) * z[col]).sum::<f64>();
                    (r[row] - model) * (r[row] + model)
                })
                .sum();
            if !predicted.is_finite() || predicted <= 0.0 {
                lambda = next_damping(lambda);
                rejected += 1;
                continue;
            }
            if evaluations == max_evaluations {
                finish!(iteration, LeastSquaresStatus::EvaluationLimit);
            }
            let trial_r = f(&trial);
            evaluations += 1;
            if trial_r.len() != m {
                return Err(LeastSquaresError::InvalidResidual);
            }
            if !finite(&trial_r) {
                lambda = next_damping(lambda);
                rejected += 1;
                continue;
            }
            let actual: f64 = residual
                .iter()
                .zip(&trial_r)
                .map(|(old, new)| {
                    let a = old / scale;
                    let b = new / scale;
                    (a - b) * (a + b)
                })
                .sum();
            if actual.is_finite() && actual > 0.0 {
                let rho = actual / predicted;
                if rho > 1e-4 {
                    x = trial;
                    residual = trial_r;
                    accepted += 1;
                    success = true;
                    if rho > 0.75 {
                        lambda = (lambda * 0.3).max(1e-300);
                    } else if rho < 0.25 {
                        lambda = (lambda * 4.0).min(1e16);
                    }
                    break;
                }
            }
            lambda = next_damping(lambda);
            rejected += 1;
        }
        if !success {
            finish!(iteration + 1, LeastSquaresStatus::Stagnated);
        }
        pending_small = last_step_small;
    }
    unreachable!("bounded iteration loop returns")
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn subnormal_projection_survives_division_before_materialization() {
        let divisor = 2.0_f64.powi(-1000);
        assert_eq!(
            scaled_quotient(1.0, -1500, divisor),
            Some(2.0_f64.powi(-500))
        );
        assert!(!magnitude_within(1.0, -2000, 0.0));
        assert!(magnitude_within(1.0, -2000, f64::from_bits(1)));
    }
}
