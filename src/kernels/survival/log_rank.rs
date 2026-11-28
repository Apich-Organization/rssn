//! The log-rank test for equality of two survival curves.

use crate::kernels::distributions::Distribution;
use crate::kernels::survival::common::check_inputs;

/// Result of the log-rank test comparing two groups.
#[derive(Clone, Debug, PartialEq)]
pub struct LogRank {
    /// Chi-square statistic with one degree of freedom.
    pub chi2: f64,
    /// Upper-tail p-value of the statistic.
    pub p_value: f64,
    /// Observed number of events in each group.
    pub observed: [f64; 2],
    /// Number of events expected under the null hypothesis in each group.
    pub expected: [f64; 2],
}

/// Sorts paired times and event flags by time.
fn sorted_pairs(
    time: &[f64],
    event: &[bool],
) -> Vec<(f64, bool)> {
    let mut pairs: Vec<(f64, bool)> = time.iter().copied().zip(event.iter().copied()).collect();
    pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    pairs
}

/// Log-rank test comparing the survival curves of two independent groups.
///
/// At every distinct event time the observed number of events in the first
/// group is compared with its hypergeometric expectation over the pooled risk
/// set; the squared standardised difference is referred to a chi-square
/// distribution with one degree of freedom.
///
/// # Errors
/// Returns an error for empty or mismatched inputs, non-finite times, or when
/// the pooled risk sets carry no variance (for example when a single subject
/// is at risk at every event time).
#[allow(clippy::float_cmp)] // exact comparison groups tied observation times
pub fn log_rank(
    time_a: &[f64],
    event_a: &[bool],
    time_b: &[f64],
    event_b: &[bool],
) -> Result<LogRank, String> {
    check_inputs(time_a, event_a)?;
    check_inputs(time_b, event_b)?;
    let a = sorted_pairs(time_a, event_a);
    let b = sorted_pairs(time_b, event_b);

    let mut candidates: Vec<f64> = a
        .iter()
        .filter(|pair| pair.1)
        .map(|pair| pair.0)
        .chain(b.iter().filter(|pair| pair.1).map(|pair| pair.0))
        .collect();
    candidates.sort_unstable_by(f64::total_cmp);
    candidates.dedup();

    let mut ia = 0usize;
    let mut ib = 0usize;
    let mut observed_a = 0.0;
    let mut expected_a = 0.0;
    let mut variance = 0.0;
    let mut events_total = 0.0;
    for &t in &candidates {
        while ia < a.len() && a.get(ia).is_some_and(|pair| pair.0 < t) {
            ia += 1;
        }
        while ib < b.len() && b.get(ib).is_some_and(|pair| pair.0 < t) {
            ib += 1;
        }
        let at_risk_a = a.len() - ia;
        let at_risk_b = b.len() - ib;
        let mut events_a = 0usize;
        let mut k = ia;
        while k < a.len() && a.get(k).is_some_and(|pair| pair.0 == t) {
            if a.get(k).is_some_and(|pair| pair.1) {
                events_a += 1;
            }
            k += 1;
        }
        let mut events_b = 0usize;
        let mut k = ib;
        while k < b.len() && b.get(k).is_some_and(|pair| pair.0 == t) {
            if b.get(k).is_some_and(|pair| pair.1) {
                events_b += 1;
            }
            k += 1;
        }

        let n = (at_risk_a + at_risk_b) as f64;
        let d = (events_a + events_b) as f64;
        observed_a += events_a as f64;
        events_total += d;
        if d > 0.0 && n > 1.0 {
            let share = at_risk_a as f64 / n;
            expected_a += d * share;
            variance += d * share * (1.0 - share) * (n - d) / (n - 1.0);
        }
    }

    if !variance.is_finite() || variance <= 0.0 {
        return Err("the pooled risk sets carry no variance".to_string());
    }
    let chi2 = (observed_a - expected_a).powi(2) / variance;
    let upper_tail = 1.0_f64 - Distribution::ChiSquared { df: 1.0 }.cdf(chi2);
    let p_value = upper_tail.clamp(0.0, 1.0);
    Ok(LogRank {
        chi2,
        p_value,
        observed: [observed_a, events_total - observed_a],
        expected: [expected_a, events_total - expected_a],
    })
}

#[cfg(test)]
mod tests {
    use super::log_rank;

    #[test]
    fn statistic_matches_hand_computation() {
        // Group A: two events at 1 and 2. Group B: two censored observations.
        // At t = 1: n1 = 2, n2 = 2, d = 1 -> E1 = 0.5, V = 0.25.
        // At t = 2: n1 = 1, n2 = 1, d = 1 -> E1 = 0.5, V = 0.25.
        // O1 = 2, E1 = 1, V = 0.5, so chi2 = 1 / 0.5 = 2.
        let result =
            log_rank(&[1.0, 2.0], &[true, true], &[1.0, 2.0], &[false, false]).expect("test");
        assert!((result.chi2 - 2.0).abs() < 1e-12, "chi2 = {}", result.chi2);
        assert!(
            (result.p_value - 0.157_299_207_1).abs() < 1e-8,
            "p = {}",
            result.p_value
        );
        assert!((result.observed[0] - 2.0).abs() < 1e-12);
        assert!((result.observed[1] - 0.0).abs() < 1e-12);
        assert!((result.expected[0] - 1.0).abs() < 1e-12);
        assert!((result.expected[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn identical_groups_give_a_zero_statistic() {
        let time = [1.0, 2.0, 3.0, 4.0];
        let event = [true, true, false, true];
        let result = log_rank(&time, &event, &time, &event).expect("test");
        assert!(result.chi2.abs() < 1e-12);
        assert!((result.p_value - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_inputs_are_rejected() {
        // Group B is no longer at risk at the only event time, so the pooled
        // risk set has a single member and no variance.
        assert!(log_rank(&[1.0], &[true], &[0.5], &[false]).is_err());
        assert!(log_rank(&[], &[], &[1.0], &[true]).is_err());
        assert!(log_rank(&[1.0, 2.0], &[true], &[1.0], &[true]).is_err());
    }
}
