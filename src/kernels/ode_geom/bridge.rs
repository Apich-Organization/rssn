//! Fixed output scheduling, Hermite segments, and generic shooting iteration.
use super::GeometricMethod;
use super::OdeGeomError;
use super::Rate;
use super::adaptive::GeometricStepper;

/// A shooting residual failure retains the caller's original error type.
#[derive(Debug)]
pub enum ShootingError<E> {
    /// The integration or boundary callback failed.
    Residual(E),
    /// The numerical iteration failed independently of the caller.
    Numerical(OdeGeomError),
}
fn residual_value<E, F: FnMut(&[f64]) -> Result<Vec<f64>, E>>(
    residual: &mut F,
    point: &[f64],
) -> Result<Vec<f64>, ShootingError<E>> {
    let values = residual(point).map_err(ShootingError::Residual)?;
    if values.len() != point.len() {
        return Err(ShootingError::Numerical(OdeGeomError::ShapeMismatch(
            format!(
                "boundary condition returned {} residuals for {} components",
                values.len(),
                point.len()
            ),
        )));
    }
    if values.iter().any(|x| !x.is_finite()) {
        return Err(ShootingError::Numerical(OdeGeomError::Domain(
            "boundary condition returned nonfinite residuals".into(),
        )));
    }
    Ok(values)
}
/// Initial state solving a boundary residual by damped finite-difference Newton.
///
/// The caller supplies an integration-backed residual. At most 60 Newton trials,
/// each with at most 30 backtracking residuals, preserve the legacy work bound.
///
/// # Errors
/// Original residual errors or a malformed, singular, nonconverged iteration.
pub fn shoot_initial<E, F: FnMut(&[f64]) -> Result<Vec<f64>, E>>(
    residual: &mut F,
    guess: &[f64],
) -> Result<Vec<f64>, ShootingError<E>> {
    if guess.is_empty() {
        return Err(ShootingError::Numerical(OdeGeomError::Empty));
    }
    if guess.iter().any(|x| !x.is_finite()) {
        return Err(ShootingError::Numerical(OdeGeomError::Domain(
            "shooting guess must be finite".into(),
        )));
    }
    let norm = |v: &[f64]| v.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    let n = guess.len();
    let mut p = guess.to_vec();
    let mut r = residual_value(residual, &p)?;
    for _ in 0..60 {
        if norm(&r) <= 1e-9 * (1.0 + norm(&p)) {
            return Ok(p);
        }
        let mut jacobian = vec![vec![0.0; n]; n];
        for column in 0..n {
            let step = f64::EPSILON.sqrt() * p[column].abs().max(1.0);
            let mut shifted = p.clone();
            shifted[column] += step;
            let shifted_r = residual_value(residual, &shifted)?;
            for row in 0..n {
                jacobian[row][column] = (shifted_r[row] - r[row]) / step;
            }
        }
        let delta = super::solve_dense(jacobian, r.iter().map(|x| -x).collect())
            .map_err(ShootingError::Numerical)?;
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..30 {
            let trial: Vec<_> = p.iter().zip(&delta).map(|(x, d)| x + scale * d).collect();
            if let Ok(trial_r) = residual_value(residual, &trial)
                && norm(&trial_r) < norm(&r)
            {
                p = trial;
                r = trial_r;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            return Err(ShootingError::Numerical(OdeGeomError::NotConverged));
        }
    }
    if norm(&r) <= 1e-9 * (1.0 + norm(&p)) {
        Ok(p)
    } else {
        Err(ShootingError::Numerical(OdeGeomError::NotConverged))
    }
}
/// Sample a requested fixed method with `steps` substeps per output interval.
///
/// # Errors
/// Invalid grids, callback values, or geometric-stage convergence failures.
pub fn sample_fixed_geometric(
    rhs: &mut Rate<'_>,
    y0: &[f64],
    ts: &[f64],
    method: GeometricMethod,
    steps: usize,
) -> Result<Vec<Vec<f64>>, OdeGeomError> {
    if ts.is_empty() {
        return Err(OdeGeomError::Empty);
    }
    if ts.iter().chain(y0).any(|x| !x.is_finite())
        || ts
            .windows(2)
            .any(|w| w[1] < w[0] || !(w[1] - w[0]).is_finite())
    {
        return Err(OdeGeomError::Domain(
            "sample times must be finite and increasing".into(),
        ));
    }
    let mut y = y0.to_vec();
    let mut t = ts[0];
    let mut out = vec![y.clone()];
    let mut stepper = GeometricStepper::new(method);
    let steps = steps.max(1);
    for target in &ts[1..] {
        let h = (target - t) / steps as f64;
        let mut time = t;
        for _ in 0..steps {
            y = stepper.step(rhs, time, &y, h)?.0;
            time += h;
        }
        t = *target;
        out.push(y.clone());
    }
    Ok(out)
}
/// Evaluate one checked cubic Hermite segment in linear state-vector work.
///
/// # Errors
/// Empty or mismatched vectors, nonfinite inputs, a query outside the segment,
/// or arithmetic overflow during the evaluation.
pub fn hermite_segment(
    t0: f64,
    t1: f64,
    y0: &[f64],
    y1: &[f64],
    f0: &[f64],
    f1: &[f64],
    t: f64,
) -> Result<Vec<f64>, OdeGeomError> {
    let n = y0.len();
    if n == 0 || [y1.len(), f0.len(), f1.len()].iter().any(|len| *len != n) {
        return Err(OdeGeomError::ShapeMismatch(
            "Hermite endpoint vectors must have equal nonzero lengths".into(),
        ));
    }
    if !t0.is_finite() || !t1.is_finite() || !t.is_finite() || t1 <= t0 || !(t0..=t1).contains(&t) {
        return Err(OdeGeomError::Domain(
            "invalid Hermite segment or query".into(),
        ));
    }
    let width = t1 - t0;
    let (h, scale, s) = if width.is_finite() {
        (width, 1.0, (t - t0) / width)
    } else {
        // Halving each endpoint keeps the interval width representable.
        // Restore that factor of two in the dimensionless slope weights.
        let half_start = 0.5 * t0;
        let half_width = 0.5 * t1 - half_start;
        (half_width, 2.0, (0.5 * t - half_start) / half_width)
    };
    let (h00, h10, h01, h11) = (
        2.0 * s * s * s - 3.0 * s * s + 1.0,
        s * s * s - 2.0 * s * s + s,
        -2.0 * s * s * s + 3.0 * s * s,
        s * s * s - s * s,
    );
    let left_rate = (h10 * scale) * h;
    let right_rate = (h11 * scale) * h;
    let output: Vec<_> = (0..n)
        .map(|i| h00 * y0[i] + left_rate * f0[i] + h01 * y1[i] + right_rate * f1[i])
        .collect();
    // All four inputs contribute even at an endpoint. A nonfinite input
    // propagates through multiplication and addition into a nonfinite output.
    if output.iter().any(|v| !v.is_finite()) {
        return Err(OdeGeomError::Domain(
            "nonfinite Hermite input or evaluation".into(),
        ));
    }
    Ok(output)
}
