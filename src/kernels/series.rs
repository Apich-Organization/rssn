//! # Numeric summation
//!
//! Partial sums, and infinite sums with convergence acceleration.

use crate::kernels::convergence::LevinVariant;
use crate::kernels::convergence::levin_transform_terms;
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
    coefficients
        .iter()
        .rev()
        .fold(0.0, |acc, &c| acc.mul_add(dx, c))
}

/// Sum of `term(k)` for `k` from `from` to `to` inclusive, with Kahan compensation.
///
/// Calls `term` once per integer in increasing order. A reversed interval
/// returns zero without calling `term`. Distinct integer indices with
/// magnitude above `2^53` may round to the same `f64` callback value.
#[allow(
    clippy::cast_precision_loss,
    reason = "The public callback takes an f64 index."
)]
pub fn sum_range(
    term: impl Fn(f64) -> f64,
    from: i64,
    to: i64,
) -> f64 {
    if from > to {
        return 0.0;
    }
    let (mut sum, mut carry) = (0.0_f64, 0.0_f64);
    // Handle `to` separately so the loop never advances past i64::MAX.
    for k in from..to {
        let y = term(k as f64) - carry;
        let t = sum + y;
        carry = (t - sum) - y;
        sum = t;
    }
    sum + (term(to as f64) - carry)
}

/// Leading terms kept for the Levin transform.
const LEVIN_TERMS: usize = 40;
const WYNN_TERMS: usize = 24;

/// Sum of `term(k)` for `k = from, from + 1, ...` to infinity.
///
/// The integer index advances through `i64::MAX`. See [`sum_terms_to_infinity`]
/// for the stopping rules. Use that entry point for recurrences or to keep
/// indices exact beyond the integer precision of the `f64` callback.
pub fn sum_to_infinity(
    term: impl Fn(f64) -> f64,
    from: i64,
    tolerance: f64,
    max_terms: usize,
) -> SeriesSum {
    #[allow(clippy::cast_precision_loss)]
    let terms = (from..=i64::MAX).map(|k| term(k as f64));
    sum_terms_to_infinity(terms, tolerance, max_terms)
}

/// Estimates an infinite sum from terms supplied in order.
///
/// This accepts an unbounded generator or a finite prefix. At each doubling
/// checkpoint it compares Wynn's epsilon algorithm on the most recent partial
/// sums, Richardson extrapolation on the doubling checkpoints, and the
/// Levin `u` transform of the first terms. A final partial block may still
/// contribute Wynn or Levin estimates. The calculation stops on the requested
/// tolerance, input exhaustion, or after evaluating `max_terms` terms.
///
/// The reported error is a change between estimates and is heuristic. When
/// extrapolation has insufficient data, the result is the available partial
/// sum with infinite error. A non-finite term or partial sum yields `NaN`
/// and infinite error, counting the triggering term in `terms`.
pub fn sum_terms_to_infinity(
    terms: impl IntoIterator<Item = f64>,
    tolerance: f64,
    max_terms: usize,
) -> SeriesSum {
    let mut terms = terms.into_iter();
    let mut recent = Vec::with_capacity(WYNN_TERMS.min(max_terms));
    let mut richardson = Vec::new();
    let mut head = Vec::with_capacity(LEVIN_TERMS.min(max_terms));
    let mut sum = 0.0_f64;
    let mut used = 0_usize;
    let mut best = (f64::NAN, f64::INFINITY);
    let (mut previous_wynn, mut previous_richardson) = (f64::NAN, f64::NAN);
    let mut levin_head = 0;
    let mut target = 16_usize;
    while used < max_terms {
        let end = target.min(max_terms);
        let start = used;
        while used < end {
            let Some(term) = terms.next() else {
                break;
            };
            used += 1;
            sum += term;
            if !sum.is_finite() {
                return SeriesSum {
                    value: f64::NAN,
                    error: f64::INFINITY,
                    terms: used,
                };
            }
            if end - used < WYNN_TERMS {
                recent.push(sum);
            }
            if head.len() < LEVIN_TERMS {
                head.push(term);
            }
        }
        if used == start {
            break;
        }
        if let Some(&wynn) = wynn_epsilon(&recent).last() {
            consider_estimate(&mut best, wynn, previous_wynn);
            previous_wynn = wynn;
        }
        // Richardson assumes an exact factor of two between checkpoints.
        if used == target && target.is_power_of_two() {
            let rich = richardson_step(&mut richardson, sum);
            consider_estimate(&mut best, rich, previous_richardson);
            previous_richardson = rich;
        }
        if head.len() >= 8 && head.len() != levin_head {
            let transforms = levin_transform_terms(&head, LevinVariant::U, 1.0);
            if let [.., before, last] = transforms.as_slice()
                && transforms.len() >= 8
            {
                consider_estimate(&mut best, *last, *before);
            }
            levin_head = head.len();
        }
        if used < end || best.1 <= tolerance.max(f64::EPSILON * best.0.abs()) {
            break;
        }
        recent.clear();
        target = target.saturating_mul(2);
    }
    SeriesSum {
        value: if best.0.is_finite() || used == 0 {
            best.0
        } else {
            sum
        },
        error: best.1,
        terms: used,
    }
}

fn consider_estimate(
    best: &mut (f64, f64),
    value: f64,
    previous: f64,
) {
    let change = (value - previous).abs();
    if value.is_finite() && change.is_finite() && change < best.1 {
        *best = (value, change);
    }
}

// The newest row of the Richardson table suffices: each entry uses the
// preceding column from this row and from the previous row.
#[allow(
    clippy::single_call_fn,
    reason = "Advance only the new Richardson row at a doubling checkpoint."
)]
fn richardson_step(
    row: &mut Vec<f64>,
    mut estimate: f64,
) -> f64 {
    let mut factor = 2.0_f64;
    for previous in row.iter_mut() {
        let previous = std::mem::replace(previous, estimate);
        estimate = factor.mul_add(estimate, -previous) / (factor - 1.0);
        factor *= 2.0;
    }
    row.push(estimate);
    estimate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_richardson_row_matches_the_full_dyadic_table() {
        let partial = [15.0, 15.5, 16.125, 16.3125, 16.413, 16.4417, 16.455_119];
        let reference = crate::kernels::convergence::richardson_extrapolation_with(&partial, 2.0);
        let mut row = Vec::new();
        for (sample, expected) in partial.into_iter().zip(reference) {
            assert_eq!(
                richardson_step(&mut row, sample).to_bits(),
                expected.to_bits()
            );
        }
    }

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
        // zeta(2) converges logarithmically.
        let s = sum_to_infinity(|k| 1.0 / (k * k), 1, 1e-13, 10_000);
        assert!(
            (s.value - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-11,
            "{s:?}"
        );
    }

    #[test]
    fn slowly_converging_positive_series() {
        // zeta(2): the error estimate must not claim more than is true.
        let s = sum_to_infinity(|k| 1.0 / (k * k), 1, 1e-10, 200_000);
        let exact = std::f64::consts::PI.powi(2) / 6.0;
        assert!((s.value - exact).abs() < 1e-9, "{s:?}");
        assert!(
            s.terms <= 70_000,
            "Richardson should not need many terms: {s:?}"
        );
    }

    #[test]
    fn divergence_is_not_reported_as_convergence() {
        let s = sum_to_infinity(|k| 1.0 / k, 1, 1e-10, 5_000);
        assert!(s.error > 1e-6 || !s.value.is_finite(), "{s:?}");
        let s = sum_to_infinity(|k| 1.0 / (k - 3.0), 1, 1e-10, 100);
        assert!(s.value.is_nan());
    }
}
