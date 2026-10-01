//! # Convergence acceleration
//!
//! Sequence transformations that speed up slowly converging sequences and
//! series: Aitken's delta-squared process, Richardson extrapolation and
//! Wynn's epsilon algorithm, plus drivers that find the limit of a sequence
//! and sum a series given as closures over the index.
//! [`crate::kernels::series::sum_to_infinity`] is built on these.

/// Sums `term(n)` for `n = start, start + 1, ...` (at most `max_terms` terms)
/// and stops, before adding it, at the first term whose magnitude is below
/// `tolerance`.
pub fn sum_series_numerical(
    term: impl Fn(f64) -> f64,
    start: usize,
    max_terms: usize,
    tolerance: f64,
) -> f64 {
    let mut sum = 0.0;
    for n in start..start.saturating_add(max_terms) {
        #[allow(clippy::cast_precision_loss)]
        let value = term(n as f64);
        if value.abs() < tolerance {
            break;
        }
        sum += value;
    }
    sum
}

/// One pass of Aitken's delta-squared process over `sequence`: entry `i` is
/// `s_i - (s_{i+1} - s_i)^2 / (s_{i+2} - 2 s_{i+1} + s_i)`. Entries whose
/// second difference is below `1e-9` in magnitude (already converged, or
/// numerically meaningless) are skipped, so the result can be shorter than
/// `len - 2`; fewer than three inputs give an empty vector.
#[must_use]
pub fn aitken_acceleration(sequence: &[f64]) -> Vec<f64> {
    sequence
        .windows(3)
        .filter_map(|w| {
            let denominator = 2.0f64.mul_add(-w[1], w[2]) + w[0];
            (denominator.abs() > 1e-9).then(|| w[0] - (w[1] - w[0]).powi(2) / denominator)
        })
        .collect()
}

/// Finds the limit of the sequence `term(0), term(1), ..., term(max_terms-1)`
/// by applying Aitken's process repeatedly until the last two accelerated
/// values differ by less than `tolerance`.
///
/// # Errors
///
/// Returns an error if a term is not finite, or if the tolerance is not
/// reached with the given number of terms.
pub fn find_sequence_limit(
    term: impl Fn(f64) -> f64,
    max_terms: usize,
    tolerance: f64,
) -> Result<f64, String> {
    #[allow(clippy::cast_precision_loss)]
    let sequence: Vec<f64> = (0..max_terms).map(|n| term(n as f64)).collect();
    if sequence.iter().any(|v| !v.is_finite()) {
        return Err("sequence has a non-finite term".to_string());
    }
    let mut accelerated = aitken_acceleration(&sequence);
    while accelerated.len() > 1 {
        let last = accelerated[accelerated.len() - 1];
        let second_last = accelerated[accelerated.len() - 2];
        if (last - second_last).abs() < tolerance {
            return Ok(last);
        }
        accelerated = aitken_acceleration(&accelerated);
    }
    Err("convergence not found".to_string())
}

/// Richardson extrapolation of approximations `A(h), A(h/2), A(h/4), ...`
/// whose error has an expansion in powers of `h^2` (Romberg integration,
/// central differences): the standard factors `4^j`. The result holds the
/// diagonal of the extrapolation table, the last entry being the highest
/// order estimate.
#[must_use]
pub fn richardson_extrapolation(sequence: &[f64]) -> Vec<f64> {
    richardson_extrapolation_with(sequence, 4.0)
}

/// Richardson extrapolation with column `j` of the table eliminating one
/// error term by the factor `growth^j`: `growth = 4` for errors in `h^2,
/// h^4, ...` as `h` halves, `growth = 2` for errors in `1/n, 1/n^2, ...` as
/// the number of terms `n` doubles. Returns the diagonal of the table.
#[must_use]
pub fn richardson_extrapolation_with(
    sequence: &[f64],
    growth: f64,
) -> Vec<f64> {
    let n = sequence.len();
    let mut table = vec![vec![0.0; n]; n];
    for (i, &value) in sequence.iter().enumerate() {
        table[i][0] = value;
    }
    let mut factor = 1.0;
    for j in 1..n {
        factor *= growth;
        for i in j..n {
            table[i][j] = factor.mul_add(table[i][j - 1], -table[i - 1][j - 1]) / (factor - 1.0);
        }
    }
    (0..n).map(|i| table[i][i]).collect()
}

/// Wynn's epsilon algorithm (Shanks transformation) on `sequence`.
///
/// Returns the accelerated estimates: the first entry is the last element of
/// the sequence itself, and each further entry is the last element of the
/// next even epsilon column, so later entries are higher-order estimates of
/// the limit. Estimates that are not finite are dropped, and the process
/// stops when two consecutive values coincide (exact convergence).
#[must_use]
pub fn wynn_epsilon(sequence: &[f64]) -> Vec<f64> {
    let Some(&first) = sequence.last() else {
        return Vec::new();
    };
    let mut estimates = vec![first];
    let mut previous = vec![0.0_f64; sequence.len() + 1];
    let mut current: Vec<f64> = sequence.to_vec();
    // `even` tells whether `current` is an even epsilon column, which holds
    // estimates of the limit; odd columns are auxiliary.
    let mut even = true;
    while current.len() > 1 {
        let mut next = Vec::with_capacity(current.len() - 1);
        for i in 0..current.len() - 1 {
            let delta = current[i + 1] - current[i];
            if delta == 0.0 {
                if even {
                    estimates.push(current[i + 1]);
                }
                return estimates;
            }
            next.push(previous[i + 1] + 1.0 / delta);
        }
        previous = current;
        current = next;
        even = !even;
        if even {
            if let Some(&last) = current.last() {
                if last.is_finite() {
                    estimates.push(last);
                }
            }
        }
    }
    estimates
}
