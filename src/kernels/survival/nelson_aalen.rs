//! Nelson-Aalen estimator of the cumulative hazard.

use crate::kernels::survival::common::check_inputs;
use crate::kernels::survival::common::risk_table;

/// Nelson-Aalen estimate of the cumulative hazard function.
///
/// Every vector has one entry per distinct time at which an event occurred, in
/// increasing order of time.
#[derive(Clone, Debug, PartialEq)]
pub struct NaEstimate {
    /// Distinct event times, increasing.
    pub times: Vec<f64>,
    /// Hazard increment `d_i / n_i` at each time.
    pub hazard: Vec<f64>,
    /// Cumulative hazard, the running sum of [`NaEstimate::hazard`].
    pub cumulative_hazard: Vec<f64>,
}

/// Nelson-Aalen estimate of the cumulative hazard function.
///
/// The increment at a time with `d` events among `n` subjects at risk is
/// `d / n`. The inputs are not modified.
///
/// # Errors
/// Returns an error for empty or mismatched inputs, non-finite times, or a
/// sample in which every observation is censored.
pub fn nelson_aalen(
    time: &[f64],
    event: &[bool],
) -> Result<NaEstimate, String> {
    check_inputs(time, event)?;
    let rows = risk_table(time, event);
    if rows.is_empty() {
        return Err("no events observed".to_string());
    }
    let mut estimate = NaEstimate {
        times: Vec::with_capacity(rows.len()),
        hazard: Vec::with_capacity(rows.len()),
        cumulative_hazard: Vec::with_capacity(rows.len()),
    };
    let mut cumulative = 0.0;
    for row in &rows {
        let increment = row.events as f64 / row.at_risk as f64;
        cumulative += increment;
        estimate.times.push(row.time);
        estimate.hazard.push(increment);
        estimate.cumulative_hazard.push(cumulative);
    }
    Ok(estimate)
}

#[cfg(test)]
mod tests {
    use super::nelson_aalen;

    #[test]
    fn cumulative_hazard_is_the_running_sum() {
        // Times 2, 3, 5, 6, 8, 9 with events at 2, 5 and 6.
        let time = [2.0, 3.0, 5.0, 6.0, 8.0, 9.0];
        let event = [true, false, true, true, false, false];
        let na = nelson_aalen(&time, &event).expect("estimate");
        assert_eq!(na.times, vec![2.0, 5.0, 6.0]);
        let increments = [1.0 / 6.0, 1.0 / 4.0, 1.0 / 3.0];
        for (got, want) in na.hazard.iter().zip(increments) {
            assert!(
                (got - want).abs() < 1e-15,
                "hazard {got} differs from {want}"
            );
        }
        let mut running = 0.0;
        for (got, &inc) in na.cumulative_hazard.iter().zip(&increments) {
            running += inc;
            assert!(
                (got - running).abs() < 1e-15,
                "cumulative hazard {got} differs from {running}"
            );
        }
    }

    #[test]
    fn tied_events_share_one_row() {
        let time = [1.0, 1.0, 1.0, 2.0];
        let event = [true, true, false, true];
        let na = nelson_aalen(&time, &event).expect("estimate");
        assert_eq!(na.times, vec![1.0, 2.0]);
        assert!((na.hazard[0] - 0.5).abs() < 1e-15);
        assert!((na.hazard[1] - 1.0).abs() < 1e-15);
        assert!((na.cumulative_hazard[1] - 1.5).abs() < 1e-15);
    }

    #[test]
    fn censored_only_sample_is_rejected() {
        assert!(nelson_aalen(&[1.0, 2.0], &[false, false]).is_err());
        assert!(nelson_aalen(&[1.0], &[true, true]).is_err());
    }
}
