//! Shared validation and risk-table construction for the survival kernels.

/// One row of a risk table: a distinct observed time together with the number
/// of subjects still at risk and the number of events at that time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct RiskRow {
    /// The observation time of the row.
    pub(super) time: f64,
    /// Subjects with an observation time of at least [`RiskRow::time`].
    pub(super) at_risk: usize,
    /// Subjects whose event happened exactly at [`RiskRow::time`].
    pub(super) events: usize,
}

/// Checks that paired time and event slices are usable and returns their
/// common length.
pub(super) fn check_inputs(
    time: &[f64],
    event: &[bool],
) -> Result<usize, String> {
    if time.len() != event.len() {
        return Err("time and event must have the same length".to_string());
    }
    if time.is_empty() {
        return Err("no observations".to_string());
    }
    if time.iter().any(|t| !t.is_finite()) {
        return Err("observation times must be finite".to_string());
    }
    Ok(time.len())
}

/// Indices ordered by observed time. Original indices break exact ties.
pub(super) fn time_order(
    time: &[f64],
    descending: bool,
) -> Vec<usize> {
    let mut order: Vec<_> = (0..time.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        let comparison = time[a].total_cmp(&time[b]);
        let comparison = if descending {
            comparison.reverse()
        } else {
            comparison
        };
        comparison.then(a.cmp(&b))
    });
    order
}

/// Finds the exclusive end of the time group beginning at `start`.
#[allow(clippy::float_cmp)] // exact comparison groups tied observation times
pub(super) fn group_end(
    time: &[f64],
    order: &[usize],
    start: usize,
) -> usize {
    let t = time[order[start]];
    start
        + order[start..]
            .iter()
            .take_while(|&&index| time[index] == t)
            .count()
}

/// Builds the risk table over the distinct times at which at least one event
/// occurred.
///
/// Subjects censored at a time `t` still count towards the risk set of `t`,
/// which is the usual convention of the Kaplan-Meier and Nelson-Aalen
/// estimators.
pub(super) fn risk_table(
    time: &[f64],
    event: &[bool],
) -> Vec<RiskRow> {
    let n = time.len();
    let order = time_order(time, false);
    let mut rows = Vec::new();
    let mut start = 0usize;
    while start < n {
        let t = time[order[start]];
        let end = group_end(time, &order, start);
        let mut events = 0usize;
        for &k in &order[start..end] {
            if event.get(k).copied().unwrap_or(false) {
                events += 1;
            }
        }
        if events > 0 {
            rows.push(RiskRow {
                time: t,
                at_risk: n - start,
                events,
            });
        }
        start = end;
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::check_inputs;
    use super::risk_table;

    #[test]
    fn risk_table_counts_censored_subjects_at_risk() {
        // 6 subjects; the censored observations at 3 and 9 stay at risk at
        // their own time but do not produce rows of their own.
        let time = [2.0, 3.0, 5.0, 6.0, 8.0, 9.0];
        let event = [true, false, true, true, false, false];
        let rows = risk_table(&time, &event);
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].at_risk, 6);
        assert_eq!(rows[1].at_risk, 4);
        assert_eq!(rows[2].at_risk, 3);
        assert_eq!(rows[0].events, 1);
        assert_eq!(rows[2].events, 1);
    }

    #[test]
    fn mismatched_lengths_are_rejected() {
        assert!(check_inputs(&[1.0, 2.0], &[true]).is_err());
        assert!(check_inputs(&[], &[]).is_err());
        assert!(check_inputs(&[f64::NAN], &[true]).is_err());
        assert!(check_inputs(&[1.0], &[true]).is_ok());
    }
}
