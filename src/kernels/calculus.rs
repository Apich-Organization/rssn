//! # Numerical Calculus
//!
//! Finite-difference derivatives of plain closures: partial derivatives,
//! gradients, Jacobians and Hessians. The graph engine uses these when a
//! derivative is wanted numerically and no symbolic rule applies.

use super::numeric_common::NumericError;
use super::numeric_common::scale_product;
use super::numeric_common::scaled_difference;

/// Step for a central difference at `x`: scaled to the magnitude of `x`,
/// near the optimum `eps^(1/3)` for a second-order formula.
fn step(x: f64) -> f64 {
    6.0e-6 * x.abs().max(1.0)
}

/// Central-difference derivative of `f` at `x`.
pub fn derivative(
    f: impl Fn(f64) -> f64,
    x: f64,
) -> f64 {
    let h = step(x);
    (f(x + h) - f(x - h)) / (2.0 * h)
}

/// Central-difference partial derivative of `f` with respect to coordinate
/// `index` at `point`. Returns NaN when `index` is out of range.
pub fn partial_derivative(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
    index: usize,
) -> f64 {
    let Some(&x) = point.get(index) else {
        return f64::NAN;
    };
    let h = step(x);
    let mut shifted = point.to_vec();
    let mut at = |value: f64| {
        if let Some(slot) = shifted.get_mut(index) {
            *slot = value;
        }
        f(&shifted)
    };
    (at(x + h) - at(x - h)) / (2.0 * h)
}

/// Gradient of `f` at `point`.
pub fn gradient(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
) -> Vec<f64> {
    (0..point.len())
        .map(|i| partial_derivative(&f, point, i))
        .collect()
}

/// Jacobian of `f: R^n -> R^m` at `point`, in row-major order (`m` rows of
/// `n` entries). `f` writes its `m` outputs into its second argument.
pub fn jacobian(
    f: impl Fn(&[f64], &mut [f64]),
    point: &[f64],
    m: usize,
) -> Vec<f64> {
    let n = point.len();
    let mut out = vec![0.0; m.saturating_mul(n)];
    let mut shifted = point.to_vec();
    let (mut plus, mut minus) = (vec![0.0; m], vec![0.0; m]);
    for (j, &x) in point.iter().enumerate() {
        let h = step(x);
        if let Some(slot) = shifted.get_mut(j) {
            *slot = x + h;
        }
        f(&shifted, &mut plus);
        if let Some(slot) = shifted.get_mut(j) {
            *slot = x - h;
        }
        f(&shifted, &mut minus);
        if let Some(slot) = shifted.get_mut(j) {
            *slot = x;
        }
        for (i, (p, q)) in plus.iter().zip(&minus).enumerate() {
            if let Some(slot) = out.get_mut(i.saturating_mul(n).saturating_add(j)) {
                *slot = (p - q) / (2.0 * h);
            }
        }
    }
    out
}

/// Hessian of `f` at `point`, in row-major order.
///
/// Uses central differences, with an inward stencil for diagonal entries
/// at a finite floating-point boundary. The result is symmetric.
pub fn hessian(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
) -> Vec<f64> {
    let n = point.len();
    let mut out = vec![0.0; n.saturating_mul(n)];
    let mut p = point.to_vec();
    let centre = f(point);
    // A second difference needs a larger step than a first one.
    let steps: Vec<f64> = point.iter().copied().map(hessian_step).collect();
    let eval = |p: &mut [f64], moves: &[(usize, f64)]| {
        for &(i, delta) in moves {
            p[i] = point[i] + delta;
        }
        let value = f(p);
        for &(i, _) in moves {
            p[i] = point[i];
        }
        value
    };
    for (i, &hi) in steps.iter().enumerate() {
        out[i * n + i] = hessian_diagonal_at(&f, &mut p, i, centre, hi);
        for (j, &hj) in steps.iter().enumerate().skip(i + 1) {
            let value = (eval(&mut p, &[(i, hi), (j, hj)])
                - eval(&mut p, &[(i, hi), (j, -hj)])
                - eval(&mut p, &[(i, -hi), (j, hj)])
                + eval(&mut p, &[(i, -hi), (j, -hj)]))
                / (4.0 * hi * hj);
            for (r, c) in [(i, j), (j, i)] {
                if let Some(slot) = out.get_mut(r.saturating_mul(n).saturating_add(c)) {
                    *slot = value;
                }
            }
        }
    }
    out
}

fn hessian_step(x: f64) -> f64 {
    1.0e-4 * x.abs().max(1.0)
}

fn distinct_offset(
    x: f64,
    delta: f64,
) -> f64 {
    let candidate = x + delta;
    if super::exact_eq(candidate, x) {
        if delta.is_sign_negative() {
            x.next_down()
        } else {
            x.next_up()
        }
    } else {
        candidate
    }
}

/// Restores the scale after a dimensionless second difference. Computing
/// raw differences first preserves cancellation between close ordinates.
fn scale_curvature<const N: usize>(
    samples: [f64; N],
    span: f64,
    stencil: impl Fn(&[f64; N], f64) -> f64,
) -> f64 {
    let numerator = stencil(&samples, 1.0);
    if numerator.is_normal()
        || (numerator == 0.0 && !samples.iter().any(|value| value.is_subnormal()))
        || samples.iter().any(|value| !value.is_finite())
    {
        return (numerator / span) / span;
    }
    let scale = samples
        .iter()
        .fold(0.0_f64, |max, value| max.max(value.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let normalized = stencil(&samples, scale);
    scale_product(normalized, &[scale], &[span, span])
        .unwrap_or_else(|_| f64::INFINITY.copysign(normalized))
}

fn hessian_diagonal_at(
    f: &impl Fn(&[f64]) -> f64,
    shifted: &mut [f64],
    index: usize,
    center: f64,
    step: f64,
) -> f64 {
    let x = shifted[index];
    if !x.is_finite() {
        return f64::NAN;
    }
    let lower = distinct_offset(x, -step);
    let upper = distinct_offset(x, step);
    let mut sample = |position| {
        shifted[index] = position;
        f(shifted)
    };
    let result = if lower.is_finite() && upper.is_finite() {
        let before = sample(lower);
        let after = sample(upper);
        let left = x - lower;
        let right = upper - x;
        let span = left.max(right);
        if super::exact_eq(left, right) {
            scale_curvature([before, center, after], span, |[a, b, c], k| {
                if super::exact_eq(k, 1.0) {
                    (c - b) - (b - a)
                } else {
                    scaled_difference(*c, *b, k) - scaled_difference(*b, *a, k)
                }
            })
        } else {
            let p = left / span;
            let q = right / span;
            scale_curvature([before, center, after], span, |[a, b, c], k| {
                let r = scaled_difference(*c, *b, k) / q;
                let l = scaled_difference(*b, *a, k) / p;
                (r - l) / f64::midpoint(p, q)
            })
        }
    } else {
        // Three inward points keep second-order accuracy at the finite
        // limits. Restrict their total range even for a very large step.
        let delta = (0.25 * step).copysign(-x);
        let x1 = distinct_offset(x, delta);
        let x2 = distinct_offset(x1, delta);
        let x3 = distinct_offset(x2, delta);
        let span = (x3 - x).abs();
        let a = (x1 - x).abs() / span;
        let b = (x2 - x).abs() / span;
        scale_curvature(
            [center, sample(x1), sample(x2), sample(x3)],
            span,
            |[f0, f1, f2, f3], k| {
                let first = scaled_difference(*f1, *f0, k) / a;
                let second = scaled_difference(*f2, *f1, k) / (b - a);
                let third = scaled_difference(*f3, *f2, k) / (1.0 - b);
                let low = (second - first) / b;
                let high = (third - second) / (1.0 - a);
                2.0 * (low - (a + b) * (high - low))
            },
        )
    };
    shifted[index] = x;
    result
}

fn visit_diagonal(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
    steps: Option<&[f64]>,
    mut visit: impl FnMut(f64),
) {
    if point.is_empty() {
        return;
    }
    let center = f(point);
    let mut shifted = point.to_vec();
    for (index, &x) in point.iter().enumerate() {
        let step = steps.map_or_else(|| hessian_step(x), |steps| steps[index]);
        visit(hessian_diagonal_at(&f, &mut shifted, index, center, step));
    }
}

/// Computes the second partial derivatives along the coordinate axes.
///
/// Uses one central evaluation and two evaluations per axis on an interior
/// stencil. At a finite floating-point boundary it uses three inward samples.
/// An empty point yields an empty diagonal.
pub fn hessian_diagonal(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
) -> Vec<f64> {
    let mut diagonal = Vec::with_capacity(point.len());
    visit_diagonal(f, point, None, |value| diagonal.push(value));
    diagonal
}

/// Computes the Hessian diagonal with a separate nominal step for each axis.
///
/// If a step would leave a coordinate unchanged, it uses the adjacent
/// representable coordinate and the actual distance in the difference.
///
/// # Errors
/// Returns [`NumericError::ShapeMismatch`] for a mismatched step count and
/// [`NumericError::Domain`] for nonfinite coordinates or unusable steps.
pub fn hessian_diagonal_with_steps(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
    steps: &[f64],
) -> Result<Vec<f64>, NumericError> {
    if point.len() != steps.len() {
        return Err(NumericError::ShapeMismatch(
            "expected one step per coordinate".into(),
        ));
    }
    if point.iter().any(|x| !x.is_finite())
        || steps.iter().any(|step| !step.is_finite() || *step <= 0.0)
    {
        return Err(NumericError::Domain(
            "diagonal Hessian requires finite coordinates and positive finite steps".into(),
        ));
    }
    let mut diagonal = Vec::with_capacity(point.len());
    visit_diagonal(f, point, Some(steps), |value| diagonal.push(value));
    Ok(diagonal)
}

/// The step a *second* difference uses.
///
/// A second difference subtracts three values that nearly cancel, so the
/// rounding noise is `eps/h^2` while the truncation error is `h^2`: the
/// balance is at `eps^(1/4)`, which is a hundred times larger than the first
/// difference's step. Using `eps^(1/3)` here loses the answer in the noise.
fn second_difference_step(x: f64) -> f64 {
    let scale = x.abs().max(1.0);
    f64::EPSILON.sqrt().sqrt() * scale
}

/// The first derivative of `f` at `x`, by a central difference with one
/// Richardson extrapolation step.
pub fn richardson_derivative<F>(
    f: &mut F,
    x: f64,
) -> f64
where
    F: FnMut(f64) -> f64,
{
    let h = f64::EPSILON.cbrt() * x.abs().max(1.0);
    let coarse = (f(x + h) - f(x - h)) / (2.0 * h);
    let fine = (f(x + h / 2.0) - f(x - h / 2.0)) / h;
    // 4/3 * fine - 1/3 * coarse cancels the h^2 term.
    (4.0 * fine - coarse) / 3.0
}

/// The second derivative of `f` at `x`.
pub fn second_derivative<F>(
    f: &mut F,
    x: f64,
) -> f64
where
    F: FnMut(f64) -> f64,
{
    let h = second_difference_step(x);
    let coarse = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    let fine = (f(x + h / 2.0) - 2.0 * f(x) + f(x - h / 2.0)) / (h * h / 4.0);
    (4.0 * fine - coarse) / 3.0
}

/// Divergence of a vector field, using central differences.
pub fn divergence(
    f: impl Fn(&[f64]) -> Vec<f64>,
    x: &[f64],
) -> f64 {
    let fx = f(x);
    let j = crate::kernels::rootfind::numerical_jacobian(&f, x, &fx);
    (0..j.rows.min(j.cols)).map(|i| j.at(i, i)).sum()
}

/// Curl of a three-dimensional vector field. Returns `None` for wrong dimensions.
pub fn curl(
    f: impl Fn(&[f64]) -> Vec<f64>,
    x: &[f64],
) -> Option<Vec<f64>> {
    if x.len() != 3 {
        return None;
    }
    let fx = f(x);
    let j = crate::kernels::rootfind::numerical_jacobian(&f, x, &fx);
    if j.rows != 3 {
        return None;
    }
    Some(vec![
        j.at(2, 1) - j.at(1, 2),
        j.at(0, 2) - j.at(2, 0),
        j.at(1, 0) - j.at(0, 1),
    ])
}

/// Laplacian of a scalar field, using only the diagonal Hessian stencil.
pub fn laplacian(
    f: impl Fn(&[f64]) -> f64,
    x: &[f64],
) -> f64 {
    let mut trace = 0.0;
    visit_diagonal(f, x, None, |value| trace += value);
    trace
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_derivatives() {
        assert!((derivative(f64::sin, 0.3) - 0.3_f64.cos()).abs() < 1e-9);
        let f = |p: &[f64]| p[0] * p[0] * p[1] + p[1].sin();
        let g = gradient(f, &[2.0, 0.5]);
        assert!((g[0] - 2.0).abs() < 1e-8);
        assert!((g[1] - (4.0 + 0.5_f64.cos())).abs() < 1e-8);
        assert!(partial_derivative(f, &[2.0, 0.5], 7).is_nan());
    }

    #[test]
    fn jacobian_of_a_map() {
        let f = |p: &[f64], out: &mut [f64]| {
            out[0] = p[0] * p[1];
            out[1] = p[0] + p[1];
            out[2] = p[0].exp();
        };
        let j = jacobian(f, &[1.0, 2.0], 3);
        let want = [2.0, 1.0, 1.0, 1.0, std::f64::consts::E, 0.0];
        for (got, want) in j.iter().zip(want) {
            assert!((got - want).abs() < 1e-8, "{got} vs {want}");
        }
    }

    #[test]
    fn hessian_is_symmetric_and_accurate() {
        let f = |p: &[f64]| p[0] * p[0] * p[1] + p[1] * p[1] * p[1];
        let h = hessian(f, &[1.0, 2.0]);
        let want = [4.0, 2.0, 2.0, 12.0];
        for (got, want) in h.iter().zip(want) {
            assert!((got - want).abs() < 1e-5, "{got} vs {want}");
        }
        assert_eq!(h[1].to_bits(), h[2].to_bits());
    }
}
