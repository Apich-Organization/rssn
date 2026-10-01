//! # Numeric summation
//!
//! Partial sums, and infinite sums with convergence acceleration.

use crate::kernels::convergence::richardson_extrapolation_with;
use crate::kernels::convergence::wynn_epsilon;

/// Result of summing an infinite series numerically.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeriesSum {
    /// Estimate of the sum.
    pub value: f64,
    /// Estimated absolute error.
    pub error: f64,
    /// Number of terms evaluated.
    pub terms: usize,
}

/// Evaluates the power series `sum_n coefficients[n] * (x - center)^n` by
/// Horner's scheme. An empty coefficient list gives 0.
#[must_use]
pub fn evaluate_power_series(
    coefficients: &[f64],
    center: f64,
    x: f64,
) -> f64 {
    let dx = x - center;
    coefficients.iter().rev().fold(0.0, |acc, &c| acc.mul_add(dx, c))
}

/// Sum of `term(k)` for `k` from `from` to `to` inclusive, with Kahan
/// compensation. Empty when `to < from`.
pub fn sum_range(
    term: impl Fn(f64) -> f64,
    from: i64,
    to: i64,
) -> f64 {
    let (mut sum, mut carry) = (0.0_f64, 0.0_f64);
    let mut k = from;
    while k <= to {
        #[allow(clippy::cast_precision_loss)]
        let y = term(k as f64) - carry;
        let t = sum + y;
        carry = (t - sum) - y;
        sum = t;
        k += 1;
    }
    sum
}

/// Sum of `term(k)` for `k = from, from + 1, ...` to infinity.
///
/// Partial sums are accumulated up to doubling checkpoints. At each one
/// two extrapolations are tried — Wynn's epsilon algorithm on the recent
/// partial sums (geometric and alternating convergence) and Richardson
/// extrapolation over the checkpoints (algebraic convergence such as
/// `1/k^2`) — and the process stops when the better of the two changes by
/// less than `tolerance`, or `max_terms` have been used. The error
/// estimate is that last change, so a series that does not converge shows
/// up as a large error rather than as a confident wrong number.
pub fn sum_to_infinity(
    term: impl Fn(f64) -> f64,
    from: i64,
    tolerance: f64,
    max_terms: usize,
) -> SeriesSum {
    let mut recent: Vec<f64> = Vec::new();
    let mut checkpoints: Vec<f64> = Vec::new();
    let mut sum = 0.0_f64;
    let mut k = from;
    let mut used = 0_usize;
    let mut best = (f64::NAN, f64::INFINITY);
    let mut previous = (f64::NAN, f64::NAN);
    let mut target = 16_usize;
    while used < max_terms {
        while used < target {
            #[allow(clippy::cast_precision_loss)]
            let t = term(k as f64);
            if !t.is_finite() {
                return SeriesSum { value: f64::NAN, error: f64::INFINITY, terms: used };
            }
            sum += t;
            k += 1;
            used += 1;
            if target - used < 24 {
                recent.push(sum);
            }
        }
        checkpoints.push(sum);
        let wynn = wynn_epsilon(&recent).last().copied().unwrap_or(0.0);
        let rich = richardson_extrapolation_with(&checkpoints, 2.0)
            .last()
            .copied()
            .unwrap_or(0.0);
        // Each extrapolation is judged by how much it moved since the
        // previous checkpoint.
        let candidates = [(wynn, (wynn - previous.0).abs()), (rich, (rich - previous.1).abs())];
        previous = (wynn, rich);
        for (value, change) in candidates {
            if value.is_finite() && change.is_finite() && change < best.1 {
                best = (value, change);
            }
        }
        if best.1 <= tolerance.max(f64::EPSILON * best.0.abs()) {
            break;
        }
        recent.clear();
        target *= 2;
    }
    SeriesSum { value: best.0, error: best.1, terms: used }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_series_evaluation() {
        // 1 + x + x^2/2 at x = 1 around 0.
        assert!((evaluate_power_series(&[1.0, 1.0, 0.5], 0.0, 1.0) - 2.5).abs() < 1e-15);
        // (x - 2)^2 expanded around 2: 0 + 0 x + 1 x^2 at x = 5.
        assert!((evaluate_power_series(&[0.0, 0.0, 1.0], 2.0, 5.0) - 9.0).abs() < 1e-9);
        assert!((evaluate_power_series(&[], 1.0, 3.0) - 0.0).abs() < 1e-9);
        assert!((evaluate_power_series(&[4.0], 1.0, 3.0) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn finite_sums() {
        assert!((sum_range(|k| k, 1, 100) - 5050.0).abs() < 1e-9);
        assert!((sum_range(|k| k, 5, 4) - 0.0).abs() < 1e-9);
        assert!((sum_range(|k| 0.1 * k.powi(0), 1, 10) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn geometric_and_alternating_series_accelerate() {
        let s = sum_to_infinity(|k| 0.5_f64.powf(k), 0, 1e-13, 10_000);
        assert!((s.value - 2.0).abs() < 1e-12, "{s:?}");
        // ln 2 = 1 - 1/2 + 1/3 - ... converges very slowly unaccelerated.
        let s = sum_to_infinity(|k| (-1.0_f64).powf(k + 1.0) / k, 1, 1e-12, 10_000);
        assert!((s.value - std::f64::consts::LN_2).abs() < 1e-10, "{s:?}");
        assert!(s.terms < 200, "acceleration should need few terms: {s:?}");
    }

    #[test]
    fn slowly_converging_positive_series() {
        // zeta(2): the error estimate must not claim more than is true.
        let s = sum_to_infinity(|k| 1.0 / (k * k), 1, 1e-10, 200_000);
        let exact = std::f64::consts::PI.powi(2) / 6.0;
        assert!((s.value - exact).abs() < 1e-9, "{s:?}");
        assert!(s.terms <= 70_000, "Richardson should not need many terms: {s:?}");
    }

    #[test]
    fn divergence_is_not_reported_as_convergence() {
        let s = sum_to_infinity(|k| 1.0 / k, 1, 1e-10, 5_000);
        assert!(s.error > 1e-6 || !s.value.is_finite(), "{s:?}");
        let s = sum_to_infinity(|k| 1.0 / (k - 3.0), 1, 1e-10, 100);
        assert!(s.value.is_nan());
    }
}
