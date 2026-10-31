//! # Numerical Calculus
//!
//! Finite-difference derivatives of plain closures: partial derivatives,
//! gradients, Jacobians and Hessians. The graph engine uses these when a
//! derivative is wanted numerically and no symbolic rule applies.

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

/// Hessian of `f` at `point`, in row-major order, by second-order central
/// differences. The result is symmetric by construction.
pub fn hessian(
    f: impl Fn(&[f64]) -> f64,
    point: &[f64],
) -> Vec<f64> {
    let n = point.len();
    let mut out = vec![0.0; n.saturating_mul(n)];
    let mut p = point.to_vec();
    let centre = f(point);
    // A second difference needs a larger step than a first one.
    let steps: Vec<f64> = point.iter().map(|x| 1.0e-4 * x.abs().max(1.0)).collect();
    let eval = |p: &mut Vec<f64>, moves: &[(usize, f64)]| {
        for &(i, delta) in moves {
            if let Some(slot) = p.get_mut(i) {
                *slot += delta;
            }
        }
        let value = f(p);
        for &(i, delta) in moves {
            if let Some(slot) = p.get_mut(i) {
                *slot -= delta;
            }
        }
        value
    };
    for i in 0..n {
        let hi = steps.get(i).copied().unwrap_or(1.0e-4);
        for j in i..n {
            let hj = steps.get(j).copied().unwrap_or(1.0e-4);
            let value = if i == j {
                (eval(&mut p, &[(i, hi)]) - 2.0 * centre + eval(&mut p, &[(i, -hi)])) / (hi * hi)
            } else {
                (eval(&mut p, &[(i, hi), (j, hj)])
                    - eval(&mut p, &[(i, hi), (j, -hj)])
                    - eval(&mut p, &[(i, -hi), (j, hj)])
                    + eval(&mut p, &[(i, -hi), (j, -hj)]))
                    / (4.0 * hi * hj)
            };
            for (r, c) in [(i, j), (j, i)] {
                if let Some(slot) = out.get_mut(r.saturating_mul(n).saturating_add(c)) {
                    *slot = value;
                }
            }
        }
    }
    out
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
