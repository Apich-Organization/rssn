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
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

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
/// [`NumericError::Empty`] for an empty sample.
pub fn median(values: &mut [f64]) -> Result<f64, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
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
/// [`NumericError::Domain`] when `q` is outside `[0, 1]`.
pub fn quantile(
    values: &[f64],
    q: f64,
) -> Result<f64, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
    if !(0.0..=1.0).contains(&q) {
        return Err(NumericError::Domain(
            "a quantile needs `0 <= q <= 1`".to_string(),
        ));
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let position = q * (sorted.len() - 1) as f64;
    let lower = position.floor().max(0.0) as usize;
    if lower + 1 >= sorted.len() {
        return Ok(sorted[lower]);
    }
    let fraction = position - lower as f64;
    Ok(sorted[lower] + fraction * (sorted[lower + 1] - sorted[lower]))
}

/// The sample covariance of two equally long samples.
///
/// # Errors
///
/// [`NumericError::ShapeMismatch`] when the lengths differ, and
/// [`NumericError::Empty`] when there are no elements.
pub fn cov(
    left: &[f64],
    right: &[f64],
    population: bool,
) -> Result<f64, NumericError> {
    if left.len() != right.len() {
        return Err(NumericError::ShapeMismatch(format!(
            "cov needs two samples of one length, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    let count = left.len();
    if count == 0 {
        return Err(NumericError::Empty);
    }
    if count == 1 {
        return Ok(0.0);
    }
    let mean_left = left.iter().sum::<f64>() / count as f64;
    let mean_right = right.iter().sum::<f64>() / count as f64;
    let total: f64 = left
        .iter()
        .zip(right)
        .map(|(a, b)| (a - mean_left) * (b - mean_right))
        .sum();
    let divisor = if population {
        count as f64
    } else {
        (count - 1) as f64
    };
    Ok(total / divisor)
}

/// The Pearson correlation coefficient of two equally long samples.
///
/// # Errors
///
/// [`NumericError::ShapeMismatch`] when the lengths differ, and
/// [`NumericError::Domain`] when one sample has no spread.
pub fn corr(
    left: &[f64],
    right: &[f64],
    population: bool,
) -> Result<f64, NumericError> {
    let covariance = cov(left, right, population)?;
    let variance_left = cov(left, left, population)?;
    let variance_right = cov(right, right, population)?;
    let denominator = (variance_left * variance_right).sqrt();
    if exact_eq(denominator, 0.0) || !denominator.is_finite() {
        return Err(NumericError::Domain(
            "a correlation needs two samples with a spread".to_string(),
        ));
    }
    Ok((covariance / denominator).clamp(-1.0, 1.0))
}

fn moments(values: &[f64]) -> Option<(f64, f64, f64, f64)> {
    if values.is_empty() {
        return None;
    }
    let count = values.len() as f64;
    let mean = values.iter().sum::<f64>() / count;
    let mut m2 = 0.0;
    let mut m3 = 0.0;
    let mut m4 = 0.0;
    for value in values {
        let deviation = value - mean;
        m2 += deviation * deviation;
        m3 += deviation * deviation * deviation;
        m4 += deviation * deviation * deviation * deviation;
    }
    Some((mean, m2 / count, m3 / count, m4 / count))
}

/// The skewness `m3 / m2^{3/2}`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample, and [`NumericError::Domain`]
/// when every element is the same (the ratio is undefined).
pub fn skewness(values: &[f64]) -> Result<f64, NumericError> {
    let (_, m2, m3, _) = moments(values).ok_or(NumericError::Empty)?;
    if exact_eq(m2, 0.0) {
        return Err(NumericError::Domain(
            "skewness of a constant sample is undefined".to_string(),
        ));
    }
    Ok(m3 / m2.powf(1.5))
}

/// The excess kurtosis `m4 / m2^2 - 3`.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample, and [`NumericError::Domain`]
/// when every element is the same.
pub fn kurtosis(values: &[f64]) -> Result<f64, NumericError> {
    let (_, m2, _, m4) = moments(values).ok_or(NumericError::Empty)?;
    if exact_eq(m2, 0.0) {
        return Err(NumericError::Domain(
            "kurtosis of a constant sample is undefined".to_string(),
        ));
    }
    Ok(m4 / (m2 * m2) - 3.0)
}

/// A uniform histogram with `bins` buckets.
///
/// Returns `(counts, edges)`; `edges` has `bins + 1` entries spanning the
/// sample's range, and a value on an interior edge belongs to the bucket above
/// it, so the last bucket is closed on both ends.
///
/// # Errors
///
/// [`NumericError::Empty`] for an empty sample and [`NumericError::Domain`]
/// when `bins` is zero.
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
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for value in values {
        low = low.min(*value);
        high = high.max(*value);
    }
    if !low.is_finite() || !high.is_finite() {
        return Err(NumericError::Domain(
            "a histogram needs finite values".to_string(),
        ));
    }
    let width = if high > low {
        (high - low) / bins as f64
    } else {
        1.0
    };
    let start = if high > low { low } else { low - 0.5 };
    let edges: Vec<f64> = (0..=bins)
        .map(|index| start + width * index as f64)
        .collect();
    let mut counts = vec![0i64; bins];
    for value in values {
        let mut index = ((value - start) / width).floor() as isize;
        if index < 0 {
            index = 0;
        }
        if index >= bins as isize {
            index = bins as isize - 1;
        }
        counts[index as usize] += 1;
    }
    Ok((counts, edges))
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
/// [`NumericError::Empty`] for an empty sample.
pub fn mode(values: &[f64]) -> Result<f64, NumericError> {
    if values.is_empty() {
        return Err(NumericError::Empty);
    }
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
