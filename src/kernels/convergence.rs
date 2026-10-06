//! # Convergence acceleration
//!
//! Sequence transformations that speed up slowly converging sequences and
//! series: Aitken's delta-squared process, Richardson extrapolation and
//! Wynn's epsilon algorithm, plus drivers that find the limit of a sequence
//! and sum a series given as closures over the index.
//! [`crate::kernels::series::sum_to_infinity`] is built on these.
//!
//! The Levin family ([`levin_transform`] with the `u`, `t` and `v`
//! remainder estimates) accelerates both alternating and logarithmically
//! convergent series, where Aitken and Wynn are weak or unstable.

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

/// One pass of Aitken's delta-squared process over `sequence`.
///
/// Entry `i` is
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

/// Richardson extrapolation of approximations `A(h), A(h/2), A(h/4), ...`.
///
/// The error of the approximations is assumed to have an expansion in powers of `h^2`; the
/// error has an expansion in powers of `h^2` (Romberg integration,
/// central differences): the standard factors `4^j`. The result holds the
/// diagonal of the extrapolation table, the last entry being the highest
/// order estimate.
#[must_use]
pub fn richardson_extrapolation(sequence: &[f64]) -> Vec<f64> {
    richardson_extrapolation_with(sequence, 4.0)
}

/// Richardson extrapolation with a configurable error-term growth factor.
///
/// Column `j` of the table eliminates one
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

/// Remainder estimates `omega_n` of the Levin-type transformations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LevinVariant {
    /// Levin `u`: `omega_n = (beta + n) a_n`; best for logarithmic and
    /// general series.
    U,
    /// Levin `t`: `omega_n = a_n`; best for alternating series.
    T,
    /// Levin `v` (Smith-Ford): `omega_n = a_n a_{n+1} / (a_n - a_{n+1})`.
    V,
}

/// Levin-type sequence transformation of the partial sums `s_0, s_1, ...`.
///
/// With the term `a_n = s_n - s_{n-1}` (`a_0 = s_0`) and the remainder
/// estimate `omega_n` of the chosen [`LevinVariant`], the entry `k` of the
/// result is the transform `L_k^{(0)}` computed with the stable
/// Fessler-Ford-Smith recurrence for numerator and denominator,
/// `L_k^{(n)} = N_k^{(n)} / D_k^{(n)}` with
/// `N_{k+1}^{(n)} = N_k^{(n+1)} - (beta+n)(beta+n+k)^{k-1} / (beta+n+k+1)^k N_k^{(n)}`.
/// Entry `0` is the (weighted) first partial sum; later entries use more
/// terms. `beta` is the shift (`1.0` is the usual choice). Entries that
/// are not finite are dropped.
#[must_use]
pub fn levin_transform(partial_sums: &[f64], variant: LevinVariant, beta: f64) -> Vec<f64> {
    let m = partial_sums.len();
    let usable = if variant == LevinVariant::V { m.saturating_sub(1) } else { m };
    if usable == 0 {
        return Vec::new();
    }
    let a: Vec<f64> = (0..m)
        .map(|n| if n == 0 { partial_sums[0] } else { partial_sums[n] - partial_sums[n - 1] })
        .collect();
    let mut num = Vec::with_capacity(usable);
    let mut den = Vec::with_capacity(usable);
    for n in 0..usable {
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        let omega = match variant {
            LevinVariant::U => (beta + nf) * a[n],
            LevinVariant::T => a[n],
            LevinVariant::V => a[n] * a[n + 1] / (a[n] - a[n + 1]),
        };
        if omega == 0.0 || !omega.is_finite() {
            break;
        }
        num.push(partial_sums[n] / omega);
        den.push(1.0 / omega);
    }
    let mut out = Vec::new();
    let mut k = 0usize;
    loop {
        if let (Some(&nu), Some(&de)) = (num.first(), den.first()) {
            let v = nu / de;
            if v.is_finite() {
                out.push(v);
            }
        }
        if num.len() < 2 {
            break;
        }
        #[allow(clippy::cast_precision_loss)]
        let kf = k as f64;
        let mut nn = Vec::with_capacity(num.len() - 1);
        let mut dd = Vec::with_capacity(num.len() - 1);
        for n in 0..num.len() - 1 {
            #[allow(clippy::cast_precision_loss)]
            let b = beta + n as f64;
            let c = b * (b + kf).powf(kf - 1.0) / (b + kf + 1.0).powf(kf);
            nn.push(num[n + 1] - c * num[n]);
            dd.push(den[n + 1] - c * den[n]);
        }
        num = nn;
        den = dd;
        k += 1;
    }
    out
}

/// Sums `term(0) + term(1) + ...` with a Levin transformation of the
/// first `n_terms` partial sums.
///
/// Returns `(estimate, error_estimate)` where the error estimate is the
/// difference of the two highest-order transforms; the best (smallest
/// successive difference) transform is reported.
///
/// # Errors
/// Returns an error if fewer than four terms are requested or no finite
/// transform could be formed.
pub fn levin_sum(
    term: impl Fn(usize) -> f64,
    n_terms: usize,
    variant: LevinVariant,
) -> Result<(f64, f64), String> {
    if n_terms < 4 {
        return Err("levin_sum needs at least four terms".to_string());
    }
    let mut s = 0.0;
    let partial: Vec<f64> = (0..n_terms)
        .map(|i| {
            s += term(i);
            s
        })
        .collect();
    let est = levin_transform(&partial, variant, 1.0);
    if est.len() < 3 {
        return Err("levin transform did not produce enough finite estimates".to_string());
    }
    let mut best = (est[est.len() - 1], f64::INFINITY);
    for i in 2..est.len() {
        let d = (est[i] - est[i - 1]).abs();
        if d < best.1 {
            best = (est[i], d);
        }
    }
    Ok(best)
}
