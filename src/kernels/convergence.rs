//! Sequence and series convergence acceleration (Aitken, Richardson, Wynn,
//! Levin on partial sums and on terms).
//!
//! # Convergence acceleration
//!
//! Sequence transformations that speed up slowly converging sequences and
//! series: Aitken's delta-squared process, Richardson extrapolation, Wynn's
//! epsilon algorithm and the Levin family, plus drivers that find the limit
//! of a sequence and sum a series given as closures, term slices or lazy
//! iterators over the terms. [`crate::kernels::series::sum_to_infinity`] is
//! built on these.
//!
//! The Levin family (`u`, `t` and `v` remainder estimates) accelerates both
//! alternating and logarithmically convergent series, where Aitken and Wynn
//! are weak or unstable. It is available on partial sums
//! ([`levin_transform`]) and, without the cancellation of recovering
//! `a_n = s_n - s_{n-1}` from rounded sums, directly on the terms
//! ([`levin_transform_terms`], [`levin_sum_lazy`]).
//! [`levin_condition_numbers`] estimates roundoff amplification at each
//! order using the Fessler–Ford–Smith recurrence in one working column.
//! [Algorithm 602: HURRY](https://doi.org/10.1145/356044.356052)
//!
//! # Limitations
//!
//! * **Real scalars only.** Every transform takes `f64` partial sums or
//!   terms; there is no complex, vector, matrix or ball/interval variant and
//!   no multi-precision path, so the accuracy floor is `kappa * eps`, not a
//!   proof.
//! * **Convergence is assumed, not detected.** [`levin_sum`] and
//!   [`levin_sum_lazy`] return the successive difference of the two highest
//!   orders as the error estimate. It is a heuristic, not a bound: a series
//!   that does not converge shows up as a large error estimate, not as an
//!   error return. [`find_sequence_limit`] is stricter and fails when its
//!   tolerance is not reached.
//! * **Divergent series only inside the model.** Levin sums divergent series
//!   whose terms alternate and fit the remainder model: `sum_{n>=1} (-1)^n n`
//!   gives `-1/4` and `sum (-1)^n n!` its Borel sum `0.59635` to `1e-11`.
//!   Monotone divergent series (`sum n`, `sum 1/n`) are not summed: the
//!   transform returns a large error estimate, a wrong value or nothing.
//!   There is no rigorous resummation theory behind this, only the
//!   transform's empirical behaviour on a model of the remainder.
//! * **Order is limited by conditioning.** The transform is a weighted mean
//!   `L = sum_j w_j s_j` with `sum_j w_j = 1`; the weights solve the
//!   Vandermonde-type system that the recurrence inverts, and
//!   `kappa = sum_j |w_j|` ([`levin_condition_number`]) is its condition
//!   number. For positive logarithmically convergent terms `kappa` grows by
//!   a factor `~3.5` per order (`kappa(11) ~ 4.6e5`, `kappa(16) ~ 2.3e8` for
//!   `sum 1/n^2`), so `kappa * eps` is `1e-10` at order 11 and `5e-8` at
//!   order 16. The measured optimum for `sum 1/n^2` is order 11 (`1.8e-11`);
//!   past order 12 the error grows again (`1e-10` at 14, `1e-8` at 17,
//!   `3.5e-7` at 19) and `kappa * eps` passes 1 near order 30, where the
//!   solve carries no digits at all. Alternating terms give a convex
//!   combination (`kappa = 1`), so conditioning is not their limit —
//!   truncation is, and they reach `1e-16` at order 13.
//! * **Fixed transforms.** Aitken needs a geometric error, Richardson needs
//!   the error expansion (powers of `h` or `h^2`), and Wynn is poor on
//!   logarithmic (`1/n`) tails. There is no automatic choice of transform,
//!   variant, order or `beta` beyond the successive-difference heuristic in
//!   [`levin_sum`], and no acceleration of infinite products or of sequences
//!   available only through a recurrence.
//! * **The early terms are needed.** The transforms consume a prefix of the
//!   series; they cannot skip a slow start or extrapolate from a tail alone.

use super::numeric_common::scale_product;
use super::numeric_common::scaled_difference;

const RICHARDSON_POWER_BLOCK: f64 = f64::from_bits(1535_u64 << 52); // 2^512

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
/// For errors that expand in powers of `h^2` (Romberg integration,
/// central differences), the standard factors are `4^j`. The result is the
/// diagonal of the extrapolation table, with the highest order estimate last.
#[must_use]
pub fn richardson_extrapolation(sequence: &[f64]) -> Vec<f64> {
    richardson_extrapolation_with(sequence, 4.0)
}

/// Richardson extrapolation with a configurable error-term growth factor.
///
/// Column `j` of the table eliminates one
/// error term by the factor `growth^j`: `growth = 4` for errors in `h^2,
/// h^4, ...` as `h` halves, `growth = 2` for errors in `1/n, 1/n^2, ...` as
/// the number of terms `n` doubles. Returns the diagonal and retains only
/// one working column. For finite `growth > 1`, a change of scale avoids
/// overflow in intermediate powers or products when the estimate is finite.
#[must_use]
pub fn richardson_extrapolation_with(
    sequence: &[f64],
    growth: f64,
) -> Vec<f64> {
    let Some(&first) = sequence.first() else {
        return Vec::new();
    };
    let n = sequence.len();
    let mut diagonal = Vec::with_capacity(n);
    diagonal.push(first);
    if n == 1 {
        return diagonal;
    }
    let mut current = sequence.to_vec();
    let stable_growth = growth.is_finite() && growth > 1.0;
    let mut factor = 1.0;
    let mut scale_blocks = 0;
    for remaining in (1..n).rev() {
        if scale_blocks < 6 {
            let larger = factor * growth;
            if stable_growth && !larger.is_finite() {
                factor = (factor / RICHARDSON_POWER_BLOCK) * (growth / RICHARDSON_POWER_BLOCK);
                scale_blocks += 2;
            } else {
                factor = larger;
            }
        }
        if scale_blocks != 0
            && richardson_correction_rounds_to_zero(&current[..=remaining], factor, scale_blocks)
        {
            diagonal.extend_from_slice(&current[1..=remaining]);
            break;
        }
        let active = &mut current[..=remaining];
        for offset in 0..remaining {
            active[offset] = richardson_entry(
                active[offset],
                active[offset + 1],
                factor,
                scale_blocks,
                stable_growth,
            );
        }
        diagonal.push(active[0]);
    }
    diagonal
}

fn richardson_entry(
    previous: f64,
    current: f64,
    factor: f64,
    scale_blocks: i32,
    stable_growth: bool,
) -> f64 {
    if scale_blocks == 0 {
        let ordinary = factor.mul_add(current, -previous) / (factor - 1.0);
        if ordinary.is_finite() || !stable_growth || !previous.is_finite() || !current.is_finite() {
            return ordinary;
        }
        return current + scaled_difference(current, previous, factor - 1.0);
    }
    if !previous.is_finite() || !current.is_finite() {
        return f64::NAN;
    }
    if scale_blocks >= 6 {
        return current;
    }
    // One block prevents overflow in the difference. Once the factor is
    // outside the float range, subtracting one has negligible relative effect.
    let difference = scaled_difference(current, previous, RICHARDSON_POWER_BLOCK);
    let factors = [
        factor,
        RICHARDSON_POWER_BLOCK,
        RICHARDSON_POWER_BLOCK,
        RICHARDSON_POWER_BLOCK,
    ];
    let denominator = &factors[..scale_blocks as usize];
    // Finite operands and a positive denominator greater than one cannot overflow.
    let correction = scale_product(difference, &[], denominator).unwrap_or(f64::NAN);
    if correction == 0.0 {
        current
    } else {
        current + correction
    }
}

fn richardson_correction_rounds_to_zero(
    current: &[f64],
    factor: f64,
    scale_blocks: i32,
) -> bool {
    let max = current.iter().try_fold(0.0_f64, |max, value| {
        value.is_finite().then(|| max.max(value.abs()))
    });
    let exponent = |value: f64| ((value.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    // The upper bound for a difference is 2^(exponent(max) + 2).
    // Leave a spare bit below half of the smallest subnormal.
    max.is_some_and(|max| {
        max == 0.0 || exponent(max) + 2 - exponent(factor) - 512 * scale_blocks <= -1076
    })
}

/// Wynn's epsilon algorithm (Shanks transformation) on `sequence`.
///
/// Returns the accelerated estimates: the first entry is the last element of
/// the sequence itself, and each further entry is the last element of the
/// next even epsilon column. Non-finite estimates are dropped. A repeated
/// even-column value supplies a final estimate, and a repeated auxiliary
/// value ends the recurrence without adding an estimate.
#[must_use]
pub fn wynn_epsilon(sequence: &[f64]) -> Vec<f64> {
    let Some(&last) = sequence.last() else {
        return Vec::new();
    };
    let mut estimates = vec![last];
    if sequence.len() < 2 {
        return estimates;
    }
    let mut work = super::reuse::f64_scratch(sequence.len().saturating_mul(2));
    work.resize(sequence.len(), 0.0);
    work.extend_from_slice(sequence);
    let (previous, current) = work.split_at_mut(sequence.len());
    wynn_epsilon_columns(
        previous,
        current,
        |delta| (delta != 0.0).then(|| 1.0 / delta),
        |value| {
            if value.is_finite() {
                estimates.push(value);
            }
        },
    );
    estimates
}

/// Visits the last estimate in each even column of Wynn's epsilon table.
/// `previous` starts as the zero auxiliary column and `current` as the input.
pub(super) fn wynn_epsilon_columns<'a, T>(
    mut previous: &'a mut [T],
    mut current: &'a mut [T],
    reciprocal: impl Fn(T) -> Option<T>,
    mut estimate: impl FnMut(T),
) where
    T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T>,
{
    let mut even = true;
    while current.len() > 1 {
        let mut slots = previous.iter_mut();
        let Some(mut slot) = slots.next() else {
            return;
        };
        // The next column only reads the previous column from one position
        // to the right, so writing from the left reuses its storage safely.
        for (&[left, right], before) in current.array_windows::<2>().zip(slots) {
            let Some(inverse) = reciprocal(right - left) else {
                if even {
                    estimate(right);
                }
                return;
            };
            *slot = *before + inverse;
            slot = before;
        }
        previous = previous.split_at_mut(current.len().saturating_sub(1)).0;
        std::mem::swap(&mut previous, &mut current);
        even = !even;
        if even && let Some(&last) = current.last() {
            estimate(last);
        }
    }
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

const fn levin_usable_len(
    len: usize,
    variant: LevinVariant,
) -> usize {
    match variant {
        | LevinVariant::V => len.saturating_sub(1),
        | _ => len,
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "the Smith–Ford remainder has its own numerical range"
)]
const fn levin_v_remainder(
    term: f64,
    next: f64,
) -> f64 {
    if !term.is_finite() || !next.is_finite() {
        return f64::NAN;
    }
    let product = term * next;
    let difference = term - next;
    if product.is_normal() && difference.is_finite() {
        return product / difference;
    }
    // Divide through by the larger term before the multiplication can
    // overflow or underflow. The quotient inside the denominator is bounded.
    if term.abs() <= next.abs() {
        term / (term / next - 1.0)
    } else {
        next / (1.0 - next / term)
    }
}

fn levin_rows(
    terms: &[f64],
    variant: LevinVariant,
    beta: f64,
    mut row: impl FnMut(usize, f64) -> [f64; 2],
) -> Vec<[f64; 2]> {
    let usable = levin_usable_len(terms.len(), variant);
    let mut rows = Vec::with_capacity(usable);
    for (n, &a) in terms.iter().take(usable).enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        let omega = match variant {
            | LevinVariant::U => (beta + nf) * a,
            | LevinVariant::T => a,
            | LevinVariant::V => {
                let Some(&next) = terms.get(n.saturating_add(1)) else {
                    break;
                };
                levin_v_remainder(a, next)
            },
        };
        if omega == 0.0 || !omega.is_finite() {
            break;
        }
        rows.push(row(n, omega));
    }
    rows
}

/// The first row of each column in the Fessler–Ford–Smith recurrence.
/// Each step reads the row to its right before that row is overwritten.
fn levin_columns(
    mut rows: Vec<[f64; 2]>,
    beta: f64,
) -> impl Iterator<Item = [f64; 2]> {
    let mut order = 0usize;
    std::iter::from_fn(move || {
        let first = *rows.first()?;
        #[allow(clippy::cast_precision_loss)]
        let k = order as f64;
        let next_len = rows.len().saturating_sub(1);
        let mut tail = rows.as_mut_slice();
        for n in 0..next_len {
            let Some((left, right)) = tail.split_first_mut() else {
                break;
            };
            let Some(&[nhi, dhi]) = right.first() else {
                break;
            };
            #[allow(clippy::cast_precision_loss)]
            let b = beta + n as f64;
            let x = b + k;
            let numerator = b * x.powf(k - 1.0);
            let denominator = (x + 1.0).powf(k);
            // Preserve the usual evaluation order. Divide the powers first
            // when their separate evaluation overflows for a positive shift.
            let factor =
                if b > 0.0 && b.is_finite() && (!numerator.is_finite() || !denominator.is_finite())
                {
                    (b / x) * (x / (x + 1.0)).powf(k)
                } else {
                    numerator / denominator
                };
            let [nlo, dlo] = *left;
            *left = [nhi - factor * nlo, dhi - factor * dlo];
            tail = right;
        }
        rows.truncate(next_len);
        order = order.saturating_add(1);
        Some(first)
    })
}

#[allow(
    clippy::indexing_slicing,
    reason = "the seed range is bounded by both slices"
)]
fn levin_core(
    terms: &[f64],
    sums: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> Vec<f64> {
    let terms = &terms[..terms.len().min(sums.len())];
    let rows = levin_rows(terms, variant, beta, |n, omega| {
        [sums[n] / omega, 1.0 / omega]
    });
    let mut out = Vec::with_capacity(rows.len());
    out.extend(
        levin_columns(rows, beta).filter_map(|[numerator, denominator]| {
            let estimate = numerator / denominator;
            estimate.is_finite().then_some(estimate)
        }),
    );
    out
}

/// Levin-type sequence transformation of the partial sums `s_0, s_1, ...`.
///
/// With the term `a_n = s_n - s_{n-1}` (`a_0 = s_0`) and the remainder
/// estimate `omega_n` of the chosen [`LevinVariant`], the transform
/// `L_k^{(0)}` is computed with the stable
/// Fessler-Ford-Smith recurrence for numerator and denominator,
/// `L_k^{(n)} = N_k^{(n)} / D_k^{(n)}` with
/// `N_{k+1}^{(n)} = N_k^{(n+1)} - (beta+n)(beta+n+k)^{k-1} / (beta+n+k+1)^k N_k^{(n)}`.
/// Order zero uses the first partial sum and later orders use more
/// terms. `beta` is the shift (`1.0` is the usual choice). The result lists
/// finite estimates in order and omits any non-finite estimates.
///
/// The terms are recovered by differencing, which cancels the leading
/// digits that consecutive partial sums share; when the terms are available
/// directly, [`levin_transform_terms`] avoids that loss.
#[must_use]
pub fn levin_transform(
    partial_sums: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> Vec<f64> {
    let mut previous = 0.0;
    let terms: Vec<f64> = partial_sums
        .iter()
        .enumerate()
        .map(|(n, &s)| {
            let term = if n == 0 { s } else { s - previous };
            previous = s;
            term
        })
        .collect();
    levin_core(&terms, partial_sums, variant, beta)
}

/// Levin-type transformation of the terms `a_0, a_1, ...` themselves.
///
/// Same transform as [`levin_transform`], but the terms are used as given
/// instead of being recovered from rounded partial sums by
/// `a_n = s_n - s_{n-1}`. That difference throws away the leading digits
/// `s_n` and `s_{n-1}` share — most of them late in a slowly converging
/// series — so this entry point keeps the full precision of the terms. The
/// partial sums are accumulated in the order the terms are given.
///
/// This is the entry point for a series whose terms are known in closed form
/// or produced one at a time (see [`levin_sum_lazy`]).
#[must_use]
pub fn levin_transform_terms(
    terms: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> Vec<f64> {
    let mut sum = 0.0;
    let sums: Vec<f64> = terms
        .iter()
        .map(|&t| {
            sum += t;
            sum
        })
        .collect();
    levin_core(terms, &sums, variant, beta)
}

/// Picks the estimate with the smallest successive difference: the best
/// agreement between two consecutive orders is the most reliable sign that
/// the transform has converged.
fn best_estimate(estimates: &[f64]) -> Result<(f64, f64), String> {
    if estimates.len() < 3 {
        return Err("levin transform did not produce enough finite estimates".to_string());
    }
    let mut best = (estimates.last().copied().unwrap_or(f64::NAN), f64::INFINITY);
    for pair in estimates.windows(2).skip(1) {
        let (Some(&previous), Some(&value)) = (pair.first(), pair.last()) else {
            continue;
        };
        let difference = (value - previous).abs();
        if difference < best.1 {
            best = (value, difference);
        }
    }
    Ok(best)
}

/// Sums `term(0) + term(1) + ...` with a Levin transformation of the
/// first `n_terms` partial sums.
///
/// Returns `(estimate, error_estimate)` where the error estimate is the
/// difference of the two highest-order transforms; the best (smallest
/// successive difference) transform is reported. Prefer
/// [`levin_sum_lazy`] when the terms are available directly, since this
/// path differences the partial sums to recover them.
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
    best_estimate(&levin_transform(&partial, variant, 1.0))
}

/// Sums a series whose terms arrive lazily, one at a time.
///
/// Pulls at most `max_terms` items from `terms` — an unbounded generator is
/// fine — and transforms the terms directly with [`levin_transform_terms`],
/// so the acceleration never depends on differencing rounded partial sums.
/// `beta` is the shift of the remainder estimate (`1.0` is usual).
///
/// Returns `(estimate, error_estimate)`, chosen by the same
/// smallest-successive-difference rule as [`levin_sum`].
///
/// # Errors
/// Returns an error if `max_terms` is below four, if the iterator yields
/// fewer than four terms, if a term is not finite, or if no finite transform
/// can be formed.
pub fn levin_sum_lazy(
    terms: impl Iterator<Item = f64>,
    max_terms: usize,
    variant: LevinVariant,
    beta: f64,
) -> Result<(f64, f64), String> {
    if max_terms < 4 {
        return Err("levin_sum_lazy needs at least four terms".to_string());
    }
    let collected: Vec<f64> = terms.take(max_terms).collect();
    if collected.len() < 4 {
        return Err(format!(
            "levin_sum_lazy got {} terms, needs four",
            collected.len()
        ));
    }
    if collected.iter().any(|t| !t.is_finite()) {
        return Err("series has a non-finite term".to_string());
    }
    best_estimate(&levin_transform_terms(&collected, variant, beta))
}

fn levin_condition_rows(
    terms: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> Vec<[f64; 2]> {
    if !beta.is_finite() || beta <= 0.0 {
        return Vec::new();
    }
    levin_rows(terms, variant, beta, |n, omega| {
        let reciprocal = 1.0 / omega;
        // The coefficient of a term has a fixed sign in each column. With
        // positive beta the row recurrence subtracts a nonnegative multiple
        // of its left neighbor, so these alternating seeds add magnitudes.
        let norm = if n.is_multiple_of(2) {
            reciprocal.abs()
        } else {
            -reciprocal.abs()
        };
        [norm, reciprocal]
    })
}

const fn levin_condition([norm, denominator]: [f64; 2]) -> f64 {
    let kappa = (norm / denominator).abs();
    if kappa.is_finite() && denominator.is_finite() {
        kappa.max(1.0)
    } else {
        f64::INFINITY
    }
}

/// Condition number `kappa = sum_j |w_j|` at the highest order the term
/// count supports. [`LevinVariant::V`] consumes one extra look-ahead term.
///
/// The transform is a weighted mean of the partial sums,
/// `L = sum_j w_j s_j` with `sum_j w_j = 1`. Those weights are the solution
/// of the Vandermonde-type system in the nodes `beta + n + j` that the
/// Fessler-Ford-Smith recurrence inverts implicitly, so `kappa` — the 1-norm
/// of that solution — is the factor by which rounding in `s_j` is amplified
/// into the estimate, and a lower bound for the 1-norm condition number of
/// the system. The accuracy floor of the transform is `kappa * eps`.
///
/// `kappa` depends on the series, which is why `terms` is an argument:
/// alternating terms give a convex combination with `kappa = 1`, while
/// positive logarithmically convergent terms give weights of alternating
/// sign whose absolute sum grows by a factor `~3.5` per order (for
/// `1/(n+1)^2`: `4.6e5` at order 11, `2.3e8` at order 16).
/// Returns infinity for a nonpositive or non-finite `beta`, an invalid
/// remainder in the required prefix, or a singular or unrepresentable
/// condition number. Fewer than two terms return `1.0` by convention.
#[must_use]
pub fn levin_condition_number(
    terms: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> f64 {
    if terms.len() < 2 {
        return 1.0;
    }
    let rows = levin_condition_rows(terms, variant, beta);
    if rows.len() != levin_usable_len(terms.len(), variant) {
        return f64::INFINITY;
    }
    levin_columns(rows, beta)
        .last()
        .map_or(f64::INFINITY, levin_condition)
}

/// Estimates the Levin weight condition number for each order starting at zero.
///
/// Each entry is the absolute sum of the normalized weights for its order.
/// This measures sensitivity to the partial sums with fixed remainder
/// estimates.
/// A singular or unrepresentable order contributes infinity at its original
/// position. The profile ends at the first invalid remainder estimate.
/// [`LevinVariant::V`] uses one additional term beyond each order's prefix.
/// A nonpositive or non-finite `beta` gives an empty profile. The recurrence
/// uses quadratic time in the valid prefix length and linear workspace in the
/// number of terms.
///
/// ```
/// use rssn::kernels::convergence::LevinVariant;
/// use rssn::kernels::convergence::levin_condition_numbers;
/// let profile = levin_condition_numbers(&[1.0, -0.5, 0.25], LevinVariant::T, 1.0);
/// assert_eq!(profile.len(), 3);
/// assert!(profile.iter().all(|kappa| (kappa - 1.0).abs() < 1e-14));
/// ```
#[must_use]
pub fn levin_condition_numbers(
    terms: &[f64],
    variant: LevinVariant,
    beta: f64,
) -> Vec<f64> {
    let rows = levin_condition_rows(terms, variant, beta);
    let mut profile = Vec::with_capacity(rows.len());
    profile.extend(levin_columns(rows, beta).map(levin_condition));
    profile
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `levin_transform` verbatim as it stood before the term-based entry
    /// points were added. The refactor moved the recurrence into
    /// [`levin_core`]; this frozen copy proves the partial-sum path still
    /// returns the same bits.
    fn legacy_levin_transform(
        partial_sums: &[f64],
        variant: LevinVariant,
        beta: f64,
    ) -> Vec<f64> {
        let m = partial_sums.len();
        let usable = if variant == LevinVariant::V {
            m.saturating_sub(1)
        } else {
            m
        };
        if usable == 0 {
            return Vec::new();
        }
        let a: Vec<f64> = (0..m)
            .map(|n| {
                if n == 0 {
                    partial_sums[0]
                } else {
                    partial_sums[n] - partial_sums[n - 1]
                }
            })
            .collect();
        let mut num = Vec::with_capacity(usable);
        let mut den = Vec::with_capacity(usable);
        for n in 0..usable {
            let nf = n as f64;
            let omega = match variant {
                | LevinVariant::U => (beta + nf) * a[n],
                | LevinVariant::T => a[n],
                | LevinVariant::V => a[n] * a[n + 1] / (a[n] - a[n + 1]),
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
            let kf = k as f64;
            let mut nn = Vec::with_capacity(num.len() - 1);
            let mut dd = Vec::with_capacity(num.len() - 1);
            for n in 0..num.len() - 1 {
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

    fn partial_sums_of(terms: &[f64]) -> Vec<f64> {
        let mut sum = 0.0;
        terms
            .iter()
            .map(|&t| {
                sum += t;
                sum
            })
            .collect()
    }

    fn zeta2_terms(count: usize) -> Vec<f64> {
        (0..count).map(|n| 1.0 / ((n + 1) as f64).powi(2)).collect()
    }

    fn alternating_harmonic(count: usize) -> Vec<f64> {
        (0..count)
            .map(|n| {
                let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                sign / ((n + 1) as f64)
            })
            .collect()
    }

    #[test]
    fn levin_series_refactor_is_bit_identical() {
        let cases = [
            zeta2_terms(40),
            alternating_harmonic(40),
            (0..24).map(|n| 0.5_f64.powi(n)).collect(),
            (1..40).map(|n| 1.0 / f64::from(n).powi(3)).collect(),
        ];
        for terms in cases {
            let sums = partial_sums_of(&terms);
            for variant in [LevinVariant::U, LevinVariant::T, LevinVariant::V] {
                let old = legacy_levin_transform(&sums, variant, 1.0);
                let new = levin_transform(&sums, variant, 1.0);
                assert_eq!(old.len(), new.len(), "{variant:?}");
                for (a, b) in old.iter().zip(&new) {
                    assert_eq!(a.to_bits(), b.to_bits(), "{variant:?}: {a} vs {b}");
                }
            }
        }
    }

    #[test]
    fn levin_series_accuracy_table() {
        let pi2 = std::f64::consts::PI * std::f64::consts::PI / 6.0;
        let ln2 = std::f64::consts::LN_2;
        let zeta3 = 1.202_056_903_159_594_2;
        let eta_half = 0.604_898_643_421_630_4;

        let cases: [(Vec<f64>, f64, &str); 4] = [
            (zeta2_terms(32), pi2, "sum 1/n^2"),
            (alternating_harmonic(32), ln2, "sum (-1)^(n+1)/n"),
            (
                (0..32).map(|n| 1.0 / f64::from(n + 1).powi(3)).collect(),
                zeta3,
                "sum 1/n^3",
            ),
            (
                (0..32)
                    .map(|n| {
                        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                        sign / f64::from(n + 1).sqrt()
                    })
                    .collect(),
                eta_half,
                "sum (-1)^n/sqrt(n+1)",
            ),
        ];
        for (terms, exact, name) in cases {
            let raw = (partial_sums_of(&terms).last().copied().unwrap_or(f64::NAN) - exact).abs();
            let est = levin_transform_terms(&terms, LevinVariant::U, 1.0);
            let best = est
                .iter()
                .map(|e| (e - exact).abs())
                .fold(f64::INFINITY, f64::min);
            let (driver, _) =
                levin_sum_lazy(terms.iter().copied(), terms.len(), LevinVariant::U, 1.0)
                    .unwrap_or_else(|e| panic!("{e}"));
            println!(
                "{name:24} raw {raw:.3e}  levin_u best {best:.3e}  driver {:.3e}",
                (driver - exact).abs()
            );
            assert!(best < raw * 1e-5, "{name}: {best} vs {raw}");
            assert!(
                (driver - exact).abs() < raw * 1e-5,
                "{name}: {driver} vs {raw}"
            );
        }
    }

    #[test]
    fn levin_series_abel_divergent_sum() {
        // sum_{n>=1} (-1)^n n = -1/4 in the Abel sense.
        let terms = (1..25_usize).map(|n| {
            let sign = if n % 2 == 1 { -1.0 } else { 1.0 };
            sign * (n as f64)
        });
        let (value, _) =
            levin_sum_lazy(terms, 24, LevinVariant::U, 1.0).unwrap_or_else(|e| panic!("{e}"));
        assert!((value + 0.25).abs() < 1e-12, "{value}");

        // sum (-1)^n n! is Borel summable to 0.5963473623231941.
        let factorial = (0..20).map(|n| {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            sign * (1..=n).map(f64::from).product::<f64>()
        });
        let (value, _) =
            levin_sum_lazy(factorial, 20, LevinVariant::T, 1.0).unwrap_or_else(|e| panic!("{e}"));
        assert!((value - 0.596_347_362_323_194_1).abs() < 1e-9, "{value}");
    }

    #[test]
    fn levin_series_lazy_iterator_stops_at_max_terms() {
        let mut produced = 0_usize;
        let series = (0..).map(|n: usize| {
            produced += 1;
            1.0 / ((n + 1) as f64).powi(2)
        });
        let (value, error) =
            levin_sum_lazy(series, 24, LevinVariant::U, 1.0).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(produced, 24);
        let pi2 = std::f64::consts::PI * std::f64::consts::PI / 6.0;
        assert!((value - pi2).abs() < 1e-9, "{value} {error}");
        // Too few terms and non-finite terms are rejected.
        assert!(levin_sum_lazy((0..3).map(f64::from), 3, LevinVariant::U, 1.0).is_err());
        assert!(
            levin_sum_lazy(
                (0..8).map(|n| if n == 3 { f64::NAN } else { 1.0 }),
                8,
                LevinVariant::U,
                1.0
            )
            .is_err()
        );
    }

    #[test]
    fn levin_series_condition_number_grows_with_order() {
        let terms = zeta2_terms(17);
        let mut previous = 0.0;
        for order in [8_usize, 11, 14, 16] {
            let kappa = levin_condition_number(&terms[..=order], LevinVariant::U, 1.0);
            println!(
                "order {order:2} kappa {kappa:.3e} kappa*eps {:.2e}",
                kappa * f64::EPSILON
            );
            assert!(kappa > previous);
            previous = kappa;
        }
        // Alternating terms give a convex combination: kappa = 1.
        let alt = alternating_harmonic(17);
        let kappa = levin_condition_number(&alt, LevinVariant::U, 1.0);
        assert!((kappa - 1.0).abs() < 1e-12, "{kappa}");
    }

    #[test]
    fn levin_condition_profile_matches_explicit_weights_in_all_variants() {
        let terms = [1.0, -0.41, 0.31, 0.17, -0.13, 0.11, -0.07, 0.04];
        for variant in [LevinVariant::U, LevinVariant::T, LevinVariant::V] {
            for beta in [0.25, 1.0, 4.5] {
                let profile = levin_condition_numbers(&terms, variant, beta);
                let look_ahead = usize::from(variant == LevinVariant::V);
                assert_eq!(profile.len(), terms.len().saturating_sub(look_ahead));
                for (order, &actual) in profile.iter().enumerate() {
                    let count = order.saturating_add(1).saturating_add(look_ahead);
                    let prefix = terms.get(..count).unwrap_or_default();
                    let expected: f64 = (0..=order)
                        .map(|column| {
                            let mut basis = vec![0.0; count];
                            if let Some(value) = basis.get_mut(column) {
                                *value = 1.0;
                            }
                            levin_core(prefix, &basis, variant, beta)
                                .last()
                                .map_or(f64::NAN, |value| value.abs())
                        })
                        .sum();
                    assert!(
                        (actual - expected).abs() < 1e-12 * expected.max(1.0),
                        "{variant:?} beta={beta} order={order}: {actual} vs {expected}"
                    );
                }
                assert_eq!(
                    profile.last().copied(),
                    Some(levin_condition_number(&terms, variant, beta))
                );
            }
        }
    }

    #[test]
    fn levin_condition_profile_identifies_singular_and_missing_orders() {
        use LevinVariant::T;
        use LevinVariant::V;

        let singular = [1.0, 1.0, 0.5];
        let profile = levin_condition_numbers(&singular, T, 1.0);
        assert_eq!(profile.len(), 3);
        assert_eq!(profile.first().copied(), Some(1.0));
        assert_eq!(profile.get(1).copied(), Some(f64::INFINITY));
        assert!(profile.last().is_some_and(|v| v.is_finite()));
        assert_eq!(levin_condition_numbers(&singular, V, 1.0), [] as [f64; 0]);
        assert!(levin_condition_number(&singular, V, 1.0).is_infinite());
        assert!(levin_condition_number(&singular, T, 0.0).is_infinite());
        assert_eq!(
            levin_condition_numbers(&singular, T, f64::NAN),
            [] as [f64; 0]
        );

        let incomplete = [1.0, 0.0, 0.5];
        assert_eq!(levin_condition_numbers(&incomplete, T, 1.0), [1.0]);
        assert!(levin_condition_number(&incomplete, T, 1.0).is_infinite());
    }

    #[test]
    fn levin_v_works_when_the_product_of_terms_is_out_of_range() {
        for scale in [1e-200, 1e200] {
            let terms = [scale, scale * 0.5, scale * 0.25, scale * 0.125];
            let estimate = levin_transform_terms(&terms, LevinVariant::V, 1.0);
            assert_eq!(estimate.len(), 3);
            assert!((estimate.last().copied().unwrap_or(f64::NAN) / scale - 2.0).abs() < 1e-14);
        }
    }

    #[test]
    fn levin_retains_high_orders_when_separate_coefficient_powers_overflow() {
        let terms = alternating_harmonic(192);
        let estimate = levin_transform_terms(&terms, LevinVariant::T, 1.0);
        assert_eq!(estimate.len(), terms.len());
        let estimate = estimate.last().copied().unwrap_or(f64::NAN);
        assert!((estimate - std::f64::consts::LN_2).abs() < 1e-13);
        assert!((levin_condition_number(&terms, LevinVariant::T, 1.0) - 1.0).abs() < 1e-13);
    }

    #[test]
    fn levin_series_weights_form_the_weighted_mean() {
        let terms = zeta2_terms(9);
        let n = terms.len();
        let mut weights = Vec::with_capacity(n);
        for j in 0..n {
            let mut basis = vec![0.0_f64; n];
            basis[j] = 1.0;
            weights.push(
                *levin_core(&terms, &basis, LevinVariant::U, 1.0)
                    .last()
                    .unwrap(),
            );
        }
        // The weights are a partition of unity; at order 8 (kappa = 1.2e4)
        // they satisfy it to kappa * eps.
        let total: f64 = weights.iter().sum();
        assert!((total - 1.0).abs() < 1e-11, "{total}");
        let sums = partial_sums_of(&terms);
        let estimate = *levin_transform_terms(&terms, LevinVariant::U, 1.0)
            .last()
            .unwrap();
        let weighted: f64 = weights.iter().zip(&sums).map(|(w, s)| w * s).sum();
        assert!(
            (weighted - estimate).abs() < 1e-11 * estimate.abs(),
            "{weighted} vs {estimate}"
        );
        // The reported condition number is the 1-norm of exactly these weights.
        let kappa: f64 = weights.iter().map(|w| w.abs()).sum();
        assert!(
            (levin_condition_number(&terms, LevinVariant::U, 1.0) - kappa).abs() < 1e-9 * kappa
        );
    }

    #[test]
    fn levin_series_terms_beat_differenced_sums_at_high_order() {
        let terms = zeta2_terms(17);
        let sums = partial_sums_of(&terms);
        let pi2 = std::f64::consts::PI * std::f64::consts::PI / 6.0;
        let by_terms = *levin_transform_terms(&terms, LevinVariant::U, 1.0)
            .last()
            .unwrap();
        let by_sums = *levin_transform(&sums, LevinVariant::U, 1.0).last().unwrap();
        println!(
            "order 16: terms err {:.3e}  differenced sums err {:.3e}",
            (by_terms - pi2).abs(),
            (by_sums - pi2).abs()
        );
        assert!((by_terms - pi2).abs() < (by_sums - pi2).abs());
    }
}
