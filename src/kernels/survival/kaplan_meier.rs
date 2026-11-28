//! Kaplan-Meier product-limit estimator of the survival function.

use crate::kernels::distributions::Distribution;
use crate::kernels::survival::common::check_inputs;
use crate::kernels::survival::common::risk_table;

/// Product-limit estimate of a survival function.
///
/// Every vector has one entry per distinct time at which an event occurred, in
/// increasing order of time.
#[derive(Clone, Debug, PartialEq)]
pub struct KmEstimate {
    /// Distinct event times, increasing.
    pub times: Vec<f64>,
    /// Subjects at risk at each time.
    pub at_risk: Vec<usize>,
    /// Events observed at each time.
    pub events: Vec<usize>,
    /// Estimated survival probability `S(t)` at each time.
    pub survival: Vec<f64>,
    /// Greenwood variance of `S(t)`; `0.0` once the risk set is exhausted by
    /// events (where `S(t) = 0`).
    pub variance: Vec<f64>,
    /// Lower confidence bound at the requested level.
    pub lower: Vec<f64>,
    /// Upper confidence bound at the requested level.
    pub upper: Vec<f64>,
}

/// Complementary log-log confidence interval for one survival probability.
///
/// The transform `log(-log S)` is used because it keeps the interval inside
/// `[0, 1]` and behaves better than the plain scale near the tails.
fn loglog_interval(
    survival: f64,
    variance: f64,
    z: f64,
) -> (f64, f64) {
    let point = survival.clamp(0.0, 1.0);
    if !survival.is_finite()
        || survival <= 0.0
        || survival >= 1.0
        || !variance.is_finite()
        || variance < 0.0
    {
        return (point, point);
    }
    let ln_s = survival.ln();
    let se = variance.sqrt() / (survival * ln_s.abs());
    if !se.is_finite() {
        return (point, point);
    }
    let centre = (-ln_s).ln();
    // `S = exp(-exp(L))` decreases in `L`, so the upper end of the interval
    // for `L` gives the lower end for `S`.
    let lower = (-(centre + z * se).exp()).exp().clamp(0.0, 1.0);
    let upper = (-(centre - z * se).exp()).exp().clamp(0.0, 1.0);
    (lower, upper)
}

/// Kaplan-Meier product-limit estimate of the survival function.
///
/// `time` holds the observed times and `event` is true when the observation
/// ended in the event of interest and false when it was censored. The inputs
/// are not modified; the estimate is computed on an internal sorted copy.
/// `conf` is the two-sided confidence level and must lie strictly between 0
/// and 1.
///
/// # Errors
/// Returns an error for empty or mismatched inputs, non-finite times, a
/// confidence level outside `(0, 1)`, or a sample in which every observation
/// is censored.
pub fn kaplan_meier(
    time: &[f64],
    event: &[bool],
    conf: f64,
) -> Result<KmEstimate, String> {
    check_inputs(time, event)?;
    if !(conf > 0.0 && conf < 1.0) {
        return Err("the confidence level must lie strictly between 0 and 1".to_string());
    }
    let rows = risk_table(time, event);
    if rows.is_empty() {
        return Err("no events observed".to_string());
    }
    let z = Distribution::Normal { mu: 0.0, sigma: 1.0 }
        .quantile(f64::midpoint(1.0, conf))
        .map_err(|_| "the confidence level is not usable".to_string())?;

    let mut estimate = KmEstimate {
        times: Vec::with_capacity(rows.len()),
        at_risk: Vec::with_capacity(rows.len()),
        events: Vec::with_capacity(rows.len()),
        survival: Vec::with_capacity(rows.len()),
        variance: Vec::with_capacity(rows.len()),
        lower: Vec::with_capacity(rows.len()),
        upper: Vec::with_capacity(rows.len()),
    };
    let mut survival = 1.0;
    let mut greenwood = 0.0;
    let mut exhausted = false;
    for row in &rows {
        let n = row.at_risk as f64;
        let d = row.events as f64;
        survival *= 1.0 - d / n;
        if row.at_risk > row.events {
            greenwood += d / (n * (n - d));
        } else {
            // The risk set is emptied by events, so S(t) = 0 and the
            // limiting Greenwood variance is 0.
            survival = 0.0;
            exhausted = true;
        }
        let variance = if exhausted {
            0.0
        } else {
            survival * survival * greenwood
        };
        let (lower, upper) = loglog_interval(survival, variance, z);
        estimate.times.push(row.time);
        estimate.at_risk.push(row.at_risk);
        estimate.events.push(row.events);
        estimate.survival.push(survival);
        estimate.variance.push(variance);
        estimate.lower.push(lower);
        estimate.upper.push(upper);
    }
    Ok(estimate)
}

/// Smallest observed time at which the estimated survival drops to at most
/// one half, or `None` when the curve never reaches one half.
#[must_use]
pub fn median_survival(km: &KmEstimate) -> Option<f64> {
    km.times
        .iter()
        .zip(&km.survival)
        .find(|&(_, &s)| s <= 0.5)
        .map(|(&t, _)| t)
}

#[cfg(test)]
mod tests {
    use super::kaplan_meier;
    use super::median_survival;

    /// Textbook set: times 2, 3, 5, 6, 8, 9 with events at 2, 5 and 6.
    fn sample() -> ([f64; 6], [bool; 6]) {
        (
            [2.0, 3.0, 5.0, 6.0, 8.0, 9.0],
            [true, false, true, true, false, false],
        )
    }

    #[test]
    fn survival_matches_hand_computation() {
        let (time, event) = sample();
        let km = kaplan_meier(&time, &event, 0.95).expect("estimate");
        assert_eq!(km.times, vec![2.0, 5.0, 6.0]);
        assert_eq!(km.at_risk, vec![6, 4, 3]);
        assert_eq!(km.events, vec![1, 1, 1]);
        // S = 5/6, then 5/6 * 3/4 = 0.625, then 0.625 * 2/3.
        let expected = [5.0 / 6.0, 0.625, 0.625 * 2.0 / 3.0];
        for (got, want) in km.survival.iter().zip(expected) {
            assert!(
                (got - want).abs() < 1e-15,
                "survival {got} differs from {want}"
            );
        }
        // The inputs are left untouched.
        assert_eq!(time, [2.0, 3.0, 5.0, 6.0, 8.0, 9.0]);
    }

    #[test]
    fn greenwood_variance_matches_hand_computation() {
        let (time, event) = sample();
        let km = kaplan_meier(&time, &event, 0.95).expect("estimate");
        // Increments 1/30, 1/12 and 1/6 of the Greenwood sum.
        let sums = [
            1.0 / 30.0,
            1.0 / 30.0 + 1.0 / 12.0,
            1.0 / 30.0 + 1.0 / 12.0 + 1.0 / 6.0,
        ];
        let survival = [5.0 / 6.0, 0.625, 0.625 * 2.0 / 3.0];
        for ((got, &s), &sum) in km.variance.iter().zip(&survival).zip(&sums) {
            let want = s * s * sum;
            assert!(
                (got - want).abs() < 1e-14,
                "variance {got} differs from {want}"
            );
        }
        // The log-log interval brackets the point estimate.
        for ((&lower, &upper), &s) in km.lower.iter().zip(&km.upper).zip(&km.survival) {
            assert!(
                lower <= s && s <= upper,
                "interval [{lower}, {upper}] misses {s}"
            );
        }
    }

    #[test]
    fn median_and_error_paths() {
        let (time, event) = sample();
        let km = kaplan_meier(&time, &event, 0.95).expect("estimate");
        // S drops to 0.4167 at t = 6, which is the first value at most 0.5.
        assert_eq!(median_survival(&km), Some(6.0));
        assert!(kaplan_meier(&time, &event, 1.0).is_err());
        assert!(kaplan_meier(&time, &event, 0.0).is_err());
        assert!(kaplan_meier(&time, &[false; 6], 0.95).is_err());
        assert!(kaplan_meier(&[1.0], &[true, false], 0.95).is_err());
        // S stays at 2/3, so the curve never reaches one half.
        let never = kaplan_meier(&[1.0, 2.0, 3.0], &[true, false, false], 0.95).expect("estimate");
        assert_eq!(median_survival(&never), None);
    }

    #[test]
    fn exhausted_risk_set_has_zero_variance_and_zero_interval() {
        let km = kaplan_meier(&[1.0, 2.0, 3.0], &[true, true, true], 0.95).expect("estimate");
        assert_eq!(km.times, vec![1.0, 2.0, 3.0]);
        assert!(km.variance.iter().all(|v| v.is_finite()));
        assert!((km.survival[2] - 0.0).abs() < 1e-15);
        assert_eq!(km.variance[2], 0.0);
        assert_eq!(km.lower[2], 0.0);
        assert_eq!(km.upper[2], 0.0);

        // Tied events exhausting the remaining risk set.
        let tied = kaplan_meier(&[1.0, 2.0, 2.0], &[true, true, true], 0.95).expect("estimate");
        assert_eq!(tied.times, vec![1.0, 2.0]);
        assert!((tied.survival[1] - 0.0).abs() < 1e-15);
        assert!(tied.variance[1].is_finite() && tied.variance[1] == 0.0);
        assert_eq!(tied.lower[1], 0.0);
        assert_eq!(tied.upper[1], 0.0);
    }
}
