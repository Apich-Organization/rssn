//! Descriptive sample statistics over floating-point and exact rational data.
//!
//! These kernels preserve the language-independent type-7 quantile, sample
//! covariance and population moment conventions. Adapters choose units and
//! numeric representations without owning these numerical algorithms.
// Slice indices follow explicit nonempty/length/range checks below. Rational
// operations cannot overflow, and all divisors are checked positive counts.
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)]
use super::exact_eq;
use super::numeric_common::NumericError;
use super::numeric_common::scale_product;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

fn finite_range(values: &[f64]) -> Result<(f64, f64), NumericError> {
    let (&first, rest) = values.split_first().ok_or(NumericError::Empty)?;
    if !first.is_finite() {
        return Err(NumericError::Domain(
            "the sample must have finite values".into(),
        ));
    }
    let (mut low, mut high) = (first, first);
    for &value in rest {
        if !value.is_finite() {
            return Err(NumericError::Domain(
                "the sample must have finite values".into(),
            ));
        }
        low = low.min(value);
        high = high.max(value);
    }
    Ok((low, high))
}

pub(super) fn interpolate(
    low: f64,
    high: f64,
    fraction: f64,
) -> f64 {
    if fraction <= 0.0 {
        return low;
    }
    if fraction >= 1.0 {
        return high;
    }
    let range = high - low;
    if range.is_finite() {
        range.mul_add(fraction, low).clamp(low, high)
    } else {
        low * (1.0 - fraction) + high * fraction
    }
}

#[derive(Default)]
struct CompensatedSum {
    total: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(
        &mut self,
        value: f64,
    ) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }

    fn finish(self) -> f64 {
        if self.total.is_finite() {
            self.total + self.correction
        } else {
            self.total
        }
    }
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = CompensatedSum::default();
    for value in values {
        sum.add(value);
    }
    sum.finish()
}

struct ScaledSample<'a> {
    values: &'a [f64],
    center: f64,
    scale: f64,
    scaled_mean: f64,
}

impl<'a> ScaledSample<'a> {
    fn from_range(
        values: &'a [f64],
        low: f64,
        high: f64,
    ) -> Self {
        let center = f64::midpoint(low, high);
        let scale = (low - center).abs().max((high - center).abs());
        let scaled_mean = if scale > 0.0 {
            compensated_sum(values.iter().map(|value| (value - center) / scale))
                / values.len() as f64
        } else {
            0.0
        };
        Self {
            values,
            center,
            scale,
            scaled_mean,
        }
    }

    fn new(values: &'a [f64]) -> Result<Self, NumericError> {
        let (low, high) = finite_range(values)?;
        Ok(Self::from_range(values, low, high))
    }

    fn deviations(&self) -> impl Iterator<Item = f64> + '_ {
        let scale = if self.scale > 0.0 {
            self.scale
        } else {
            1.0
        };
        self.values
            .iter()
            .map(move |value| (value - self.center) / scale - self.scaled_mean)
    }

    fn variance(
        &self,
        population: bool,
    ) -> f64 {
        let total = compensated_sum(self.deviations().map(|deviation| deviation * deviation));
        total / divisor(self.values.len(), population)
    }
}

const fn divisor(
    count: usize,
    population: bool,
) -> f64 {
    if population || count == 1 {
        count as f64
    } else {
        (count - 1) as f64
    }
}

fn check_paired_lengths(
    left: usize,
    right: usize,
) -> Result<(), NumericError> {
    if left == right {
        Ok(())
    } else {
        Err(NumericError::ShapeMismatch(format!(
            "cov needs two samples of one length, found {left} and {right}"
        )))
    }
}

fn sort_exact(values: &mut [BigRational]) {
    values.sort();
}

/// The median of an exact sample.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample.
pub fn median_exact(values: &mut [BigRational]) -> Result<BigRational, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    sort_exact(values);
    let len = values.len();
    if len % 2 == 1 {
        return Ok(values[len / 2].clone());
    }
    Ok((values[len / 2 - 1].clone() + values[len / 2].clone())
        / BigRational::from_integer(BigInt::from(2)))
}

/// The median of a float sample.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] for non-finite values.
pub fn median(values: &mut [f64]) -> Result<f64, NumericError> {
    finite_range(values)?;
    values.sort_by(f64::total_cmp);
    let len = values.len();
    if len % 2 == 1 {
        return Ok(values[len / 2]);
    }
    Ok(f64::midpoint(values[len / 2 - 1], values[len / 2]))
}

/// The `q`-quantile of a sorted exact sample, `0 ≤ q ≤ 1`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] when `q` is outside `[0, 1]`.
pub fn quantile_exact(
    sorted: &[BigRational],
    q: &BigRational,
) -> Result<BigRational, NumericError> {
    if sorted.is_empty() {
        return Err(NumericError::Empty);
    }
    if q.is_negative() || *q > BigRational::one() {
        return Err(NumericError::Domain(
            "a quantile needs `0 <= q <= 1`".to_string(),
        ));
    }
    let position = q.clone() * BigRational::from_integer(BigInt::from(sorted.len() - 1));
    // Do not round an exact position through f64: a value just below an
    // integer can otherwise select the following interval.
    let lower = position.to_integer().to_usize().ok_or_else(|| {
        NumericError::Domain("quantile position is outside the sample index range".into())
    })?;
    if lower + 1 >= sorted.len() {
        return Ok(sorted[lower].clone());
    }
    let fraction = position - BigRational::from_integer(BigInt::from(lower));
    let below = sorted[lower].clone();
    let above = sorted[lower + 1].clone();
    Ok(below.clone() + fraction * (above - below))
}

/// The `q`-quantile of a float sample, `0 ≤ q ≤ 1`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] for non-finite values or a `q` outside `[0, 1]`.
pub fn quantile(
    values: &[f64],
    q: f64,
) -> Result<f64, NumericError> {
    let (low, high) = finite_range(values)?;
    check_probabilities(&[q])?;
    // Preserve the total-order treatment of signed zeros at the endpoints.
    if q == 0.0 && low != 0.0 {
        return Ok(low);
    }
    if q >= 1.0 && high != 0.0 {
        return Ok(high);
    }
    let (rank, fraction) = type7_position(values.len(), q);
    let mut sample = values.to_vec();
    let (_, lower, upper) = sample.select_nth_unstable_by(rank, f64::total_cmp);
    if fraction == 0.0 {
        return Ok(*lower);
    }
    let next = upper
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .unwrap_or(*lower);
    Ok(interpolate(*lower, next, fraction))
}

fn check_probabilities(probabilities: &[f64]) -> Result<(), NumericError> {
    if probabilities.iter().any(|p| !(0.0..=1.0).contains(p)) {
        return Err(NumericError::Domain(
            "a quantile needs `0 <= q <= 1`".into(),
        ));
    }
    Ok(())
}

pub(super) fn type7_position(
    len: usize,
    q: f64,
) -> (usize, f64) {
    let last = len - 1;
    let position = q * last as f64;
    let lower = position.floor().to_usize().unwrap_or(last).min(last);
    (lower, position - lower as f64)
}

/// Type-7 quantiles of a finite sample, in the supplied probability order.
///
/// A single query selects the two adjacent ranks in linear time. Multiple
/// queries share one ordered copy of the sample.
///
/// # Errors
/// Returns [`NumericError::Empty`] for an empty sample and [`NumericError::Domain`]
/// for non-finite observations or a probability outside `[0, 1]`.
pub fn quantiles(
    values: &[f64],
    probabilities: &[f64],
) -> Result<Vec<f64>, NumericError> {
    if let [probability] = probabilities {
        return Ok(vec![quantile(values, *probability)?]);
    }
    finite_range(values)?;
    check_probabilities(probabilities)?;
    if probabilities.is_empty() {
        return Ok(Vec::new());
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    Ok(probabilities
        .iter()
        .map(|&q| {
            let (rank, fraction) = type7_position(sorted.len(), q);
            interpolate(
                sorted[rank],
                sorted[(rank + 1).min(sorted.len() - 1)],
                fraction,
            )
        })
        .collect())
}

/// The arithmetic mean of a float sample.
///
/// Non-finite inputs follow IEEE addition in sample order.
///
/// # Errors
///
/// Returns [`NumericError::Empty`] for an empty sample.
pub fn mean(values: &[f64]) -> Result<f64, NumericError> {
    match values {
        | [] => return Err(NumericError::Empty),
        | [value] => return Ok(*value),
        | _ => {},
    }
    let count = values.len() as f64;
    let total = compensated_sum(values.iter().copied());
    if total.is_finite() {
        return Ok(total / count);
    }
    let Ok((low, high)) = finite_range(values) else {
        return Ok(total / count);
    };
    // Divide terms first when the finite sample has an overflowing sum.
    let mean = compensated_sum(values.iter().map(|value| value / count));
    if mean.is_finite() {
        return Ok(mean.clamp(low, high));
    }
    // Rounding the individual quotients can still overflow near f64::MAX.
    let sample = ScaledSample::from_range(values, low, high);
    Ok(sample
        .scaled_mean
        .mul_add(sample.scale, sample.center)
        .clamp(low, high))
}

/// The variance of a float sample.
///
/// The population divisor is the sample length. The sample divisor is one
/// less than the length for samples with two or more observations. A singleton
/// has zero variance.
///
/// # Errors
///
/// Returns [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] for non-finite values or a result outside the
/// finite range.
pub fn var(
    values: &[f64],
    population: bool,
) -> Result<f64, NumericError> {
    let sample = ScaledSample::new(values)?;
    if sample.scale > 0.0 {
        scale_product(
            sample.variance(population),
            &[sample.scale, sample.scale],
            &[],
        )
    } else {
        Ok(0.0)
    }
}

/// The standard deviation of a float sample.
///
/// Uses the same divisor as [`var`]. The standard deviation can be finite
/// even when the variance exceeds the float range.
///
/// # Errors
///
/// Returns [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] for non-finite values or a result outside the
/// finite range.
pub fn std_dev(
    values: &[f64],
    population: bool,
) -> Result<f64, NumericError> {
    let (normalized, scale) = std_dev_parts(values, population)?;
    if scale > 0.0 {
        scale_product(normalized, &[scale], &[])
    } else {
        Ok(0.0)
    }
}

pub(super) fn std_dev_parts(
    values: &[f64],
    population: bool,
) -> Result<(f64, f64), NumericError> {
    let sample = ScaledSample::new(values)?;
    let normalized = if sample.scale > 0.0 {
        sample.variance(population).sqrt()
    } else {
        0.0
    };
    Ok((normalized, sample.scale))
}

/// The sample covariance of two equally long samples.
///
/// Use `population` to divide by the sample length. Otherwise, divide by one
/// less than the length for samples with two or more observations. The
/// covariance of a pair of singleton samples is zero.
///
/// # Errors
///
/// Returns [`NumericError::ShapeMismatch`] for different lengths,
/// [`NumericError::Empty`] for empty samples, and [`NumericError::Domain`]
/// for non-finite values or a result outside the finite range.
pub fn cov(
    left: &[f64],
    right: &[f64],
    population: bool,
) -> Result<f64, NumericError> {
    check_paired_lengths(left.len(), right.len())?;
    let (left, right) = (ScaledSample::new(left)?, ScaledSample::new(right)?);
    if left.scale <= 0.0 || right.scale <= 0.0 {
        return Ok(0.0);
    }
    let total = compensated_sum(
        left.deviations()
            .zip(right.deviations())
            .map(|(a, b)| a * b),
    );
    let count = divisor(left.values.len(), population);
    if total.abs() < f64::MIN_POSITIVE && left.scale.max(right.scale) > 1.0 {
        // Separate scaling can underflow a tiny observation before it is
        // multiplied by a large partner. Recover such products in the
        // original units when the raw dot product is representable.
        let (left_mean, right_mean) = (mean(left.values)?, mean(right.values)?);
        let raw_total = compensated_sum(
            left.values
                .iter()
                .zip(right.values)
                .map(|(a, b)| (a - left_mean) * (b - right_mean)),
        );
        if raw_total.is_finite() && !exact_eq(raw_total, 0.0) {
            return Ok(raw_total / count);
        }
    }
    scale_product(total / count, &[left.scale, right.scale], &[])
}

/// The Pearson correlation coefficient of two equally long samples.
///
/// Population and sample normalization give the same coefficient.
///
/// # Errors
///
/// Returns [`NumericError::ShapeMismatch`] for different lengths,
/// [`NumericError::Empty`] for empty samples, and [`NumericError::Domain`]
/// for non-finite values or a sample with no spread.
pub fn corr(
    left: &[f64],
    right: &[f64],
    _population: bool,
) -> Result<f64, NumericError> {
    check_paired_lengths(left.len(), right.len())?;
    let (left, right) = (ScaledSample::new(left)?, ScaledSample::new(right)?);
    if left.scale <= 0.0 || right.scale <= 0.0 {
        return Err(NumericError::Domain(
            "a correlation needs two samples with a spread".into(),
        ));
    }
    let mut covariance = CompensatedSum::default();
    let (mut left_squared, mut right_squared) = (0.0, 0.0);
    for (a, b) in left.deviations().zip(right.deviations()) {
        covariance.add(a * b);
        left_squared = a.mul_add(a, left_squared);
        right_squared = b.mul_add(b, right_squared);
    }
    Ok((covariance.finish() / (left_squared.sqrt() * right_squared.sqrt())).clamp(-1.0, 1.0))
}

fn moments(values: &[f64]) -> Result<(f64, f64, f64), NumericError> {
    let sample = ScaledSample::new(values)?;
    let count = values.len() as f64;
    let (mut m2, mut m3, mut m4) = (0.0, CompensatedSum::default(), 0.0);
    for deviation in sample.deviations() {
        let squared = deviation * deviation;
        m2 += squared;
        m3.add(squared * deviation);
        m4 = squared.mul_add(squared, m4);
    }
    Ok((m2 / count, m3.finish() / count, m4 / count))
}

/// The skewness `m3 / m2^{3/2}`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample, and [`NumericError::Domain`]
/// for a constant sample or non-finite values.
pub fn skewness(values: &[f64]) -> Result<f64, NumericError> {
    let (m2, m3, _) = moments(values)?;
    if exact_eq(m2, 0.0) {
        return Err(NumericError::Domain(
            "skewness of a constant sample is undefined".to_string(),
        ));
    }
    Ok(m3 / (m2 * m2.sqrt()))
}

/// The excess kurtosis `m4 / m2^2 - 3`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample, and [`NumericError::Domain`]
/// for a constant sample or non-finite values.
pub fn kurtosis(values: &[f64]) -> Result<f64, NumericError> {
    let (m2, _, m4) = moments(values)?;
    if exact_eq(m2, 0.0) {
        return Err(NumericError::Domain(
            "kurtosis of a constant sample is undefined".to_string(),
        ));
    }
    Ok(m4 / (m2 * m2) - 3.0)
}

/// A uniform histogram with `bins` buckets.
///
/// Returns `(counts, edges)`. For a nonconstant sample, the `bins + 1`
/// edges span the sample range. An interior edge belongs to the bucket above
/// it, and the highest edge belongs to the last bucket. Edges can coincide
/// when a bucket is narrower than the representable spacing. Constant samples
/// use unit-width buckets starting at the value minus 0.5, with all observations
/// in the first bucket.
///
/// Large samples with a bounded number of buckets can be counted in parallel.
/// The counts and edges are independent of the active Rayon pool size.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and [`NumericError::Domain`]
/// for non-finite values, a zero bucket count, or buckets that cannot be allocated.
pub fn histogram(
    values: &[f64],
    bins: usize,
) -> Result<(Vec<i64>, Vec<f64>), NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    if bins == 0 {
        return Err(NumericError::Domain(
            "a histogram needs a bucket".to_string(),
        ));
    }
    let (low, high) = finite_range(values)?;
    let edge_capacity = bins
        .checked_add(1)
        .ok_or_else(|| NumericError::Domain("histogram has too many buckets".into()))?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(edge_capacity)
        .map_err(|_| NumericError::Domain("histogram has too many buckets".into()))?;
    let mut counts = empty_histogram_counts(bins)?;
    if exact_eq(low, high) {
        edges.extend((0..edge_capacity).map(|index| (low - 0.5) + index as f64));
        counts[0] = i64::try_from(values.len())
            .map_err(|_| NumericError::Domain("histogram count exceeds Int range".into()))?;
        return Ok((counts, edges));
    }

    let bins_as_float = bins as f64;
    edges.extend(
        (0..edge_capacity).map(|index| interpolate(low, high, index as f64 / bins_as_float)),
    );
    let scale = if (high - low).is_finite() {
        1.0
    } else {
        0.5
    };
    let start = low * scale;
    let range = high * scale - start;
    let last = bins - 1;
    count_histogram(values, &mut counts, |value| {
        // Converting the nonnegative estimate to an index already truncates it.
        let estimate = (value * scale - start) / range * bins_as_float;
        let index = estimate.to_usize().unwrap_or(last).min(last);
        // The direct estimate is generally enough. Consult the reported edges
        // if rounding moves it across a boundary or yields duplicate edges.
        if value >= edges[index] && (index == last || value < edges[index + 1]) {
            index
        } else {
            edges[1..bins].partition_point(|&edge| edge <= value)
        }
    });
    Ok((counts, edges))
}

/// Histogram counts for caller-supplied, possibly nonuniform bucket edges.
///
/// Returns one count per adjacent pair of edges. Every bucket includes its
/// lower edge, and the last bucket also includes its upper edge. Every sample
/// appears in one bucket, including observations exactly on the outer edges.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample. [`NumericError::Domain`] if
/// the edges contain fewer than two values, fail to increase strictly, or
/// would require unallocatable counter storage. Non-finite observations or
/// edges and observations outside the supplied range also give a domain error.
pub fn histogram_edges(
    values: &[f64],
    edges: &[f64],
) -> Result<Vec<i64>, NumericError> {
    let (minimum, maximum) = finite_range(values)?;
    let [first, interior @ .., last] = edges else {
        return Err(NumericError::Domain(
            "a histogram needs at least two edges".into(),
        ));
    };
    if !first.is_finite() || !last.is_finite() || !edges.windows(2).all(|pair| pair[0] < pair[1]) {
        return Err(NumericError::Domain(
            "histogram edges must be finite and strictly increasing".into(),
        ));
    }
    if minimum < *first || maximum > *last {
        return Err(NumericError::Domain(
            "histogram edges must enclose every observation".into(),
        ));
    }
    let mut counts = empty_histogram_counts(edges.len() - 1)?;
    count_histogram(values, &mut counts, |value| {
        interior.partition_point(|&edge| edge <= value)
    });
    Ok(counts)
}

fn empty_histogram_counts(bins: usize) -> Result<Vec<i64>, NumericError> {
    let mut counts = Vec::new();
    counts
        .try_reserve_exact(bins)
        .map_err(|_| NumericError::Domain("histogram has too many buckets".into()))?;
    counts.resize(bins, 0);
    Ok(counts)
}

fn count_histogram(
    values: &[f64],
    counts: &mut [i64],
    bucket: impl Fn(f64) -> usize + Sync,
) {
    const MIN_SHARD_VALUES: usize = 64 * 1024;
    // Limit private counter storage and final merge traffic for wide histograms.
    const MAX_PRIVATE_BUCKETS: usize = 32 * 1024;
    let bins = counts.len();
    let shard_limit = (values.len() / MIN_SHARD_VALUES).min(MAX_PRIVATE_BUCKETS / bins);
    // A slice of f64 has fewer than i64::MAX elements on supported platforms.
    let count = |input: &[f64], output: &mut [i64]| {
        for &value in input {
            output[bucket(value)] += 1;
        }
    };
    if shard_limit < 2 {
        count(values, counts);
        return;
    }
    let shards = shard_limit.min(crate::execution::cpu_threads());
    if shards < 2 {
        count(values, counts);
        return;
    }
    let shard_len = values.len().div_ceil(shards);
    let Ok(mut private_counts) = empty_histogram_counts(values.len().div_ceil(shard_len) * bins)
    else {
        count(values, counts);
        return;
    };
    super::par::for_chunks_mut(&mut private_counts, bins, 2, |index, row| {
        let start = index * shard_len;
        count(
            &values[start..start.saturating_add(shard_len).min(values.len())],
            row,
        );
    });
    for row in private_counts.chunks_exact(bins) {
        for (total, subtotal) in counts.iter_mut().zip(row) {
            *total += subtotal;
        }
    }
}

/// The most frequent value of an exact sample.
///
/// Ties go to the smallest value, so the answer does not depend on the order of
/// the input.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample.
pub fn mode_exact(values: &[BigRational]) -> Result<BigRational, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    let mut sorted = values.to_vec();
    sort_exact(&mut sorted);
    let mut best = sorted[0].clone();
    let mut best_count = 1usize;
    let mut current = sorted[0].clone();
    let mut current_count = 1usize;
    for value in &sorted[1..] {
        if *value == current {
            current_count += 1;
        } else {
            current = value.clone();
            current_count = 1;
        }
        if current_count > best_count {
            best_count = current_count;
            best = current.clone();
        }
    }
    Ok(best)
}

/// The most frequent value of a float sample; ties go to the smallest.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and
/// [`NumericError::Domain`] for non-finite values.
pub fn mode(values: &[f64]) -> Result<f64, NumericError> {
    finite_range(values)?;
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut best = sorted[0];
    let mut best_count = 1usize;
    let mut current = sorted[0];
    let mut current_count = 1usize;
    for value in &sorted[1..] {
        if exact_eq(*value, current) {
            current_count += 1;
        } else {
            current = *value;
            current_count = 1;
        }
        if current_count > best_count {
            best_count = current_count;
            best = current;
        }
    }
    Ok(best)
}

/// The arithmetic mean of an exact sample.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample.
pub fn mean_exact(values: &[BigRational]) -> Result<BigRational, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    let total = values
        .iter()
        .fold(BigRational::zero(), |acc, value| acc + value.clone());
    Ok(total / BigRational::from_integer(BigInt::from(values.len())))
}

/// The variance of an exact sample.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample.
pub fn var_exact(
    values: &[BigRational],
    population: bool,
) -> Result<BigRational, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    let mean = mean_exact(values)?;
    let total = values.iter().fold(BigRational::zero(), |acc, value| {
        let deviation = value.clone() - mean.clone();
        acc + deviation.clone() * deviation
    });
    let count = if population || values.len() == 1 {
        values.len()
    } else {
        values.len() - 1
    };
    Ok(total / BigRational::from_integer(BigInt::from(count)))
}

/// The covariance of two exact samples.
///
/// # Errors
///
/// [`NumericError::ShapeMismatch`] when the lengths differ, and
/// [`NumericError::Empty`] when there are no elements.
pub fn cov_exact(
    left: &[BigRational],
    right: &[BigRational],
    population: bool,
) -> Result<BigRational, NumericError> {
    if left.len() != right.len() {
        return Err(NumericError::ShapeMismatch(format!(
            "cov needs two samples of one length, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    if left.is_empty() {
        return Err(NumericError::Empty);
    }
    let mean_left = mean_exact(left)?;
    let mean_right = mean_exact(right)?;
    let total = left
        .iter()
        .zip(right)
        .fold(BigRational::zero(), |acc, (a, b)| {
            acc + (a.clone() - mean_left.clone()) * (b.clone() - mean_right.clone())
        });
    let count = if population || left.len() == 1 {
        left.len()
    } else {
        left.len() - 1
    };
    Ok(total / BigRational::from_integer(BigInt::from(count)))
}
