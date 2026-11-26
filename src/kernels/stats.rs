//! # Numerical Statistics
//!
//! This module provides numerical statistical functions, leveraging the `statrs` crate
//! for robust implementations. It includes descriptive statistics (mean, variance, median,
//! percentiles, covariance, correlation), probability distributions (Normal, Uniform, Binomial,
//! Poisson, Exponential, Gamma), and statistical inference methods (ANOVA, t-tests).
//!
//! # Limitations and coverage
//!
//! * Only the six wrapped distributions exist (`Normal`, `Uniform`,
//!   `Binomial`, `Poisson`, `Exponential`, `Gamma`); there is no beta,
//!   log-normal, chi-square, Student-t or multivariate distribution, no
//!   parameter fitting and no sampling (the samplers live in
//!   [`crate::kernels::random`]).
//! * The inference helpers (`one_way_anova`, `two_sample_t_test`,
//!   `welch_t_test`, `chi_squared_test`) report the statistic and the
//!   asymptotic p-value only: no confidence intervals,
//!   effect sizes, multiple-comparison correction, paired or repeated-measures
//!   designs, and no non-parametric alternative.
//! * `mean` returns `0` for an empty slice and `std_dev` follows the sample
//!   convention; NaN and infinite inputs propagate through the descriptive
//!   functions unless stated otherwise, and the in-place helpers (`median`,
//!   `min`, `max`, `mode`, `skewness`, `kurtosis`, `iqr`) sort or partially
//!   sort their input.
//! * The distribution objects wrap `statrs`; the valid parameter ranges and
//!   the numerical accuracy are inherited from it. `inverse_cdf` returns `NaN`
//!   outside `[0, 1]` instead of panicking. The estimators in this module take
//!   the full sample; streaming means and covariances are available through
//!   [`crate::kernels::online_stats::OnlineCovariance`].

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use statrs::distribution::Binomial;
use statrs::distribution::Continuous;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Discrete;
use statrs::distribution::DiscreteCDF;
use statrs::distribution::Normal;
use statrs::distribution::Uniform;
use statrs::statistics::Data;
use statrs::statistics::Max;
use statrs::statistics::Min;
use statrs::statistics::OrderStatistics;
use statrs::statistics::Statistics;

use crate::kernels::krylov::exact_product_sum_integer;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::numeric_common::scaled_difference;
pub use crate::kernels::online_stats::OnlineCovariance;
pub use crate::kernels::online_stats::OnlineStatsError;

/// Computes the mean of a slice of data.
///
/// # Arguments
/// * `data` - A unmutable slice of `f64` data points.
///
/// # Returns
/// The mean of the data as an `f64`.
#[must_use]
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }

    let n = data.len() as f64;
    let direct = crate::kernels::simd::sum_f64(data) / n;
    if direct.is_finite() {
        return direct;
    }
    // The plain sum overflowed (or the data hold a non-finite value). Redo it
    // with scaled terms; a genuine NaN/inf in the data still propagates.
    if data.iter().all(|x| x.is_finite()) {
        data.iter().map(|x| x / n).sum()
    } else {
        direct
    }
}

/// Represents the type of variance to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarianceType {
    /// Population variance (divided by n).
    Population,
    /// Sample variance (divided by n - 1).
    Sample,
}

/// Computes the variance of a slice of data with a specified variance type (Population or Sample).
///
/// This implementation uses Welford's algorithm to compute the variance in a single pass,
/// which is highly numerically stable even for very large numbers or shifted datasets.
///
/// # Arguments
/// * `data` - A slice of `f64` data points.
/// * `var_type` - The type of variance to compute (Population or Sample).
///
/// # Returns
/// `Some(variance)` if `data.len() >= 2`, otherwise `None`.
#[must_use]
pub fn variance_with_type(
    data: &[f64],
    var_type: VarianceType,
) -> Option<f64> {
    if data.is_empty() || (var_type == VarianceType::Sample && data.len() < 2) {
        return None;
    }

    let mut n = 0.0;
    let mut mean = 0.0;
    let mut m2 = 0.0;

    for &x in data {
        n += 1.0;
        let delta = x - mean;
        let delta_n = delta / n;
        m2 += delta * delta_n * (n - 1.0);
        mean += delta_n;
    }

    let divisor = match var_type {
        | VarianceType::Population => n,
        | VarianceType::Sample => n - 1.0,
    };

    Some(m2 / divisor)
}

/// Computes the variance of a slice of f64 values (defaults to population variance).
///
/// # Arguments
/// * `data` - An immutable slice of `f64` data points.
///
/// # Returns
/// The variance of the data as an `f64`.
#[must_use]
pub fn variance(data: &[f64]) -> f64 {
    variance_with_type(data, VarianceType::Population).unwrap_or(0.0)
}

/// Computes the standard deviation of a slice of data.
///
/// # Arguments
/// * `data` - A slice of `f64` data points.
///
/// # Returns
/// The standard deviation of the data as an `f64`.
#[must_use]
pub fn std_dev(data: &[f64]) -> f64 {
    Statistics::std_dev(data.iter())
}

/// Computes the median of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
///
/// # Returns
/// The median of the data as an `f64`.
pub fn median(data: &mut [f64]) -> f64 {
    if data.is_empty() || data.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }

    let n = data.len();
    let mid = n / 2;
    let (_, upper, _) = data.select_nth_unstable_by(mid, f64::total_cmp);
    let upper_val = *upper;
    if n % 2 == 1 {
        upper_val
    } else {
        let lower_val = data[..mid]
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        f64::midpoint(lower_val, upper_val)
    }
}

/// Computes the p-th percentile of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
/// * `p` - The desired percentile (e.g., 50.0 for median).
///
/// # Returns
/// The p-th percentile of the data (fractional `p` allowed), or NaN for empty
/// data, data containing NaN, or `p` outside `0..=100`.
pub fn percentile(
    data: &mut [f64],
    p: f64,
) -> f64 {
    if data.is_empty() || !(0.0..=100.0).contains(&p) || data.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }

    let mut data_container = Data::new(data);

    data_container.quantile(p / 100.0)
}

/// Computes the sample covariance of two slices of data.
///
/// # Arguments
/// * `data1` - The first slice of `f64` data points.
/// * `data2` - The second slice of `f64` data points.
///
/// # Returns
/// The covariance of the two data sets, or NaN when the slices differ in
/// length or hold fewer than two points.
#[must_use]
pub fn covariance(
    data1: &[f64],
    data2: &[f64],
) -> f64 {
    if data1.len() != data2.len() || data1.len() < 2 {
        return f64::NAN;
    }

    let mean1 = mean(data1);

    let mean2 = mean(data2);

    let n = data1.len() as f64;

    crate::kernels::simd::centered_stats_f64(data1, data2, mean1, mean2).0 / (n - 1.0)
}

/// Computes the Pearson correlation coefficient of two slices of data.
///
/// # Arguments
/// * `data1` - The first slice of `f64` data points.
/// * `data2` - The second slice of `f64` data points.
///
/// # Returns
/// The Pearson correlation coefficient as an `f64`.
#[must_use]
pub fn correlation(
    data1: &[f64],
    data2: &[f64],
) -> f64 {
    if data1.len() != data2.len() || data1.len() < 2 {
        return f64::NAN;
    }

    let mean1 = mean(data1);
    let mean2 = mean(data2);
    let (cov, var1, var2) = crate::kernels::simd::centered_stats_f64(data1, data2, mean1, mean2);
    cov / (var1 * var2).sqrt()
}

/// Represents a Normal (Gaussian) distribution.
pub struct NormalDist(Normal);

impl NormalDist {
    /// Creates a new `NormalDist` instance.
    ///
    /// # Errors
    /// Returns an error unless both parameters are finite and `std_dev > 0`.
    pub fn new(
        mean: f64,
        std_dev: f64,
    ) -> Result<Self, String> {
        Normal::new(mean, std_dev)
            .map(NormalDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability density function (PDF) value at `x`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// The PDF value as an `f64`.
    #[must_use]
    pub fn pdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.pdf(x)
    }

    /// Returns the cumulative distribution function (CDF) value at `x`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the CDF.
    ///
    /// # Returns
    /// The CDF value as an `f64`.
    #[must_use]
    pub fn cdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.cdf(x)
    }

    /// Returns the quantile (inverse CDF) at probability `p`.
    ///
    /// Returns `NaN` when `p` is outside `[0, 1]`; the boundary values are the
    /// (possibly infinite) support limits rather than errors.
    #[must_use]
    pub fn inverse_cdf(
        &self,
        p: f64,
    ) -> f64 {
        if !(0.0..=1.0).contains(&p) {
            return f64::NAN;
        }
        self.0.inverse_cdf(p)
    }
}

/// Represents a Uniform distribution.
pub struct UniformDist(Uniform);

impl UniformDist {
    /// Creates a new `UniformDist` instance.
    ///
    /// # Arguments
    /// * `min` - The minimum value of the distribution.
    /// * `max` - The maximum value of the distribution.
    ///
    /// # Returns
    /// A `Result` containing the `UniformDist` instance, or an error string if parameters are invalid.
    ///
    /// # Errors
    /// Returns an error if `min >= max`.
    pub fn new(
        min: f64,
        max: f64,
    ) -> Result<Self, String> {
        Uniform::new(min, max)
            .map(UniformDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability density function (PDF) value at `x`.
    #[must_use]
    pub fn pdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.pdf(x)
    }

    /// Returns the cumulative distribution function (CDF) value at `x`.
    #[must_use]
    pub fn cdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.cdf(x)
    }

    /// Returns the quantile (inverse CDF) at probability `p`.
    ///
    /// Returns `NaN` when `p` is outside `[0, 1]`.
    #[must_use]
    pub fn inverse_cdf(
        &self,
        p: f64,
    ) -> f64 {
        if !(0.0..=1.0).contains(&p) {
            return f64::NAN;
        }
        self.0.inverse_cdf(p)
    }
}

/// Represents a Binomial distribution.
pub struct BinomialDist(Binomial);

impl BinomialDist {
    /// Creates a new `BinomialDist` instance.
    ///
    /// # Errors
    /// Returns an error if `p < 0.0` or `p > 1.0`.
    pub fn new(
        n: u64,
        p: f64,
    ) -> Result<Self, String> {
        Binomial::new(p, n)
            .map(BinomialDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability mass function (PMF) value at `k`.
    #[must_use]
    pub fn pmf(
        &self,
        k: u64,
    ) -> f64 {
        self.0.pmf(k)
    }

    /// Returns the cumulative distribution function (CDF) value at `k`.
    #[must_use]
    pub fn cdf(
        &self,
        k: u64,
    ) -> f64 {
        self.0.cdf(k)
    }
}

/// Fits the least-squares line `y = intercept + slope * x` to paired samples.
///
/// Returns `(slope, intercept)`. Empty and nonfinite samples return two
/// NaNs. A constant predictor returns a NaN slope and the response mean.
/// The calculation borrows the observations and rescales centered moments
/// when the unscaled moments exceed the Float range.
#[must_use]
pub fn simple_linear_regression(data: &[(f64, f64)]) -> (f64, f64) {
    linear_fit(data.iter().copied())
}

pub(crate) fn simple_linear_regression_columns(
    x: &[f64],
    y: &[f64],
) -> (f64, f64) {
    debug_assert_eq!(x.len(), y.len());
    linear_fit(x.iter().copied().zip(y.iter().copied()))
}

fn linear_fit(data: impl ExactSizeIterator<Item = (f64, f64)> + Clone) -> (f64, f64) {
    let n = data.len() as f64;
    if data.len() == 0 {
        return (f64::NAN, f64::NAN);
    }

    let (sum_x, sum_y) = data
        .clone()
        .fold((-0.0, -0.0), |(sx, sy), (x, y)| (sx + x, sy + y));
    let (mean_x, mean_y) = (sum_x / n, sum_y / n);
    let (xx, xy) = fit_moments(data.clone().map(|(x, y)| (x - mean_x, y - mean_y)));
    let slope = xy / xx;
    let intercept = slope.mul_add(-mean_x, mean_y);
    // A rounded mean can distort the slope of nearly coincident predictors.
    let origin_resolution = mean_x * (1.0 / 16_777_216.0);
    if xx.is_normal()
        && xy.is_normal()
        && slope.is_normal()
        && intercept.is_finite()
        && xx >= origin_resolution * origin_resolution * n
    {
        return (slope, intercept);
    }
    scaled_linear_fit(data, n)
}

fn fit_moments(centered: impl Iterator<Item = (f64, f64)>) -> (f64, f64) {
    centered.fold((-0.0, -0.0), |(xx, xy), (x, y)| (xx + x * x, xy + x * y))
}

#[cold]
fn scaled_linear_fit(
    data: impl ExactSizeIterator<Item = (f64, f64)> + Clone,
    n: f64,
) -> (f64, f64) {
    let Some((first_x, first_y)) = data.clone().next() else {
        return (f64::NAN, f64::NAN);
    };
    let mut max_x = 0.0_f64;
    let mut max_y = 0.0_f64;
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    for (x, y) in data.clone() {
        if !x.is_finite() || !y.is_finite() {
            return (f64::NAN, f64::NAN);
        }
        max_x = max_x.max(x.abs());
        max_y = max_y.max(y.abs());
        if x != 0.0 {
            min_x = min_x.min(x.abs());
        }
        if y != 0.0 {
            min_y = min_y.min(y.abs());
        }
    }
    // Scaling to a single Float per axis can discard small observations
    // before their products are formed when an axis spans many exponents.
    let wide =
        |smallest: f64, largest: f64| (largest.to_bits() >> 52) > (smallest.to_bits() >> 52) + 400;
    if wide(min_x, max_x) || wide(min_y, max_y) {
        return exact_linear_fit(data);
    }

    let binary_scale = |largest: f64| {
        let exponent = largest.to_bits() & (0x7ff0_u64 << 48);
        if largest == 0.0 {
            1.0
        } else {
            f64::from_bits(exponent.max(f64::MIN_POSITIVE.to_bits()))
        }
    };
    let sx = binary_scale(max_x);
    let sy = binary_scale(max_y);
    let shifted = data.map(|(x, y)| {
        (
            scaled_difference(x, first_x, sx),
            scaled_difference(y, first_y, sy),
        )
    });
    // Center near the first sampled point before taking the means. Closely
    // spaced observations remain distinct even if their mean rounds to a sample.
    let (sum_dx, sum_dy) = shifted
        .clone()
        .fold((0.0, 0.0), |(dx, dy), (x, y)| (dx + x, dy + y));
    let (dx, dy) = (sum_dx / n, sum_dy / n);
    let (xx, xy) = fit_moments(shifted.map(|(x, y)| (x - dx, y - dy)));
    let normalized_y = first_y / sy + dy;
    if xx == 0.0 {
        return (f64::NAN, normalized_y * sy);
    }

    let normalized_slope = xy / xx;
    let slope = scale_product(normalized_slope, &[sy], &[sx])
        .unwrap_or_else(|_| f64::INFINITY.copysign(normalized_slope));
    let fitted_at_anchor = normalized_slope.mul_add(-dx, dy) + first_y / sy;
    let normalized_intercept = normalized_slope.mul_add(-(first_x / sx), fitted_at_anchor);
    (slope, normalized_intercept * sy)
}

fn exact_linear_fit(data: impl ExactSizeIterator<Item = (f64, f64)> + Clone) -> (f64, f64) {
    const SAMPLE_SHIFT: usize = 1074;
    let n = BigInt::from(data.len());
    // Product bins use 2^-2148. For sums of individual samples the
    // coefficient is divisible by 2^1074, so x and y use 2^-1074.
    let x = exact_product_sum_integer(data.clone().map(|(x, _)| (x, 1.0))) >> SAMPLE_SHIFT;
    let y = exact_product_sum_integer(data.clone().map(|(_, y)| (y, 1.0))) >> SAMPLE_SHIFT;
    let xx = &n * exact_product_sum_integer(data.clone().map(|(x, _)| (x, x))) - &x * &x;
    let xy = &n * exact_product_sum_integer(data) - &x * &y;
    let response_scale = n << SAMPLE_SHIFT;
    if xx.sign() == num_bigint::Sign::NoSign {
        let mean_y = BigRational::new_raw(y, response_scale)
            .to_f64()
            .unwrap_or(f64::NAN);
        return (f64::NAN, mean_y);
    }
    let numerator = y * &xx - &xy * x;
    let denominator = response_scale * &xx;
    let intercept = BigRational::new_raw(numerator, denominator)
        .to_f64()
        .unwrap_or(f64::NAN);
    let slope = BigRational::new_raw(xy, xx).to_f64().unwrap_or(f64::NAN);
    (slope, intercept)
}

/// Computes the minimum value of a slice of data.
pub fn min(data: &mut [f64]) -> f64 {
    if data.is_empty() || data.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }

    let data_container = Data::new(data);

    data_container.min()
}

/// Computes the maximum value of a slice of data.
pub fn max(data: &mut [f64]) -> f64 {
    if data.is_empty() || data.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }

    let data_container = Data::new(data);

    data_container.max()
}

/// Computes the adjusted Fisher-Pearson sample skewness coefficient `G1`
/// (the estimator used by Excel's `SKEW` and by `scipy.stats.skew(bias=False)`).
///
/// With `g1 = m3 / m2^(3/2)` (central moments divided by `n`),
/// `G1 = sqrt(n (n - 1)) / (n - 2) * g1`.
///
/// Returns `NaN` for fewer than three points and `0.0` for constant data.
pub fn skewness(data: &mut [f64]) -> f64 {
    let n = data.len() as f64;

    if n < 3.0 {
        return f64::NAN;
    }

    let mean = mean(data);

    let m2 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;

    let m3 = data.iter().map(|&x| (x - mean).powi(3)).sum::<f64>() / n;

    if m2 == 0.0 {
        return 0.0;
    }

    let g1 = m3 / m2.powf(1.5);

    (n * (n - 1.0)).sqrt() / (n - 2.0) * g1
}

/// Computes the unbiased sample excess kurtosis `G2` (the estimator used by
/// Excel's `KURT` and `scipy.stats.kurtosis(fisher=True, bias=False)`); a
/// normal distribution has expected value 0.
///
/// With `g2 = m4 / m2^2 - 3` (central moments divided by `n`),
/// `G2 = (n - 1) / ((n - 2)(n - 3)) * ((n + 1) g2 + 6)`.
///
/// Returns `NaN` for fewer than four points and `0.0` for constant data.
pub fn kurtosis(data: &mut [f64]) -> f64 {
    let n = data.len() as f64;

    if n < 4.0 {
        return f64::NAN;
    }

    let mean = mean(data);

    let m2 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;

    let m4 = data.iter().map(|&x| (x - mean).powi(4)).sum::<f64>() / n;

    if m2 == 0.0 {
        return 0.0;
    }

    let g2 = m4 / m2.powi(2) - 3.0;

    (n - 1.0) / ((n - 2.0) * (n - 3.0)) * (n + 1.0).mul_add(g2, 6.0)
}

/// Represents a Poisson distribution.
pub struct PoissonDist(statrs::distribution::Poisson);

impl PoissonDist {
    /// Creates a new `PoissonDist` instance.
    ///
    /// # Errors
    /// Returns an error if `rate <= 0.0`.
    pub fn new(rate: f64) -> Result<Self, String> {
        statrs::distribution::Poisson::new(rate)
            .map(PoissonDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability mass function (PMF) value at `k`.
    #[must_use]
    pub fn pmf(
        &self,
        k: u64,
    ) -> f64 {
        self.0.pmf(k)
    }

    /// Returns the cumulative distribution function (CDF) value at `k`.
    #[must_use]
    pub fn cdf(
        &self,
        k: u64,
    ) -> f64 {
        self.0.cdf(k)
    }
}

/// Represents an Exponential distribution.
pub struct ExponentialDist(statrs::distribution::Exp);

impl ExponentialDist {
    /// Creates a new `ExponentialDist` instance.
    ///
    /// # Errors
    /// Returns an error if `rate <= 0.0`.
    pub fn new(rate: f64) -> Result<Self, String> {
        statrs::distribution::Exp::new(rate)
            .map(ExponentialDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability density function (PDF) value at `x`.
    #[must_use]
    pub fn pdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.pdf(x)
    }

    /// Returns the cumulative distribution function (CDF) value at `x`.
    #[must_use]
    pub fn cdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.cdf(x)
    }

    /// Returns the quantile (inverse CDF) at probability `p`.
    ///
    /// Returns `NaN` when `p` is outside `[0, 1]`.
    #[must_use]
    pub fn inverse_cdf(
        &self,
        p: f64,
    ) -> f64 {
        if !(0.0..=1.0).contains(&p) {
            return f64::NAN;
        }
        self.0.inverse_cdf(p)
    }
}

/// Represents a Gamma distribution.
pub struct GammaDist(statrs::distribution::Gamma);

impl GammaDist {
    /// Creates a new `GammaDist` instance.
    ///
    /// # Errors
    /// Returns an error if `shape <= 0.0` or `rate <= 0.0`.
    pub fn new(
        shape: f64,
        rate: f64,
    ) -> Result<Self, String> {
        statrs::distribution::Gamma::new(shape, rate)
            .map(GammaDist)
            .map_err(|e| e.to_string())
    }

    /// Returns the probability density function (PDF) value at `x`.
    #[must_use]
    pub fn pdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.pdf(x)
    }

    /// Returns the cumulative distribution function (CDF) value at `x`.
    #[must_use]
    pub fn cdf(
        &self,
        x: f64,
    ) -> f64 {
        self.0.cdf(x)
    }

    /// Returns the quantile (inverse CDF) at probability `p`.
    ///
    /// Returns `NaN` when `p` is outside `[0, 1]`.
    #[must_use]
    pub fn inverse_cdf(
        &self,
        p: f64,
    ) -> f64 {
        if !(0.0..=1.0).contains(&p) {
            return f64::NAN;
        }
        self.0.inverse_cdf(p)
    }
}

/// Performs a One-Way Analysis of Variance (ANOVA) test.
/// Tests the null hypothesis that the means of two or more groups are equal.
///
/// # Arguments
/// * `groups` - A slice of slices, where each inner slice is a group of data.
///
/// # Returns
/// A tuple containing the F-statistic and the p-value.
pub fn one_way_anova(groups: &mut [&mut [f64]]) -> (f64, f64) {
    let k = groups.len() as f64;

    if k < 2.0 {
        return (f64::NAN, f64::NAN);
    }

    let all_data: Vec<f64> = groups.iter().flat_map(|g| g.iter()).copied().collect();

    let n_total = all_data.len() as f64;

    let grand_mean = mean(&all_data);

    let mut ss_between = 0.0;

    for group in groups.iter_mut() {
        let n_group = group.len() as f64;

        let mean_group = mean(group);

        ss_between += n_group * (mean_group - grand_mean).powi(2);
    }

    let df_between = k - 1.0;

    let ms_between = ss_between / df_between;

    let mut ss_within = 0.0;

    for group in groups.iter_mut() {
        let mean_group = mean(group);

        ss_within += group.iter().map(|&x| (x - mean_group).powi(2)).sum::<f64>();
    }

    let df_within = n_total - k;

    let ms_within = ss_within / df_within;

    if ms_within == 0.0 {
        return (f64::INFINITY, 0.0);
    }

    let f_stat = ms_between / ms_within;

    let Ok(f_dist) = statrs::distribution::FisherSnedecor::new(df_between, df_within) else {
        return (f64::NAN, f64::NAN);
    };

    if f_stat.is_nan() {
        return (f_stat, f64::NAN);
    }
    let p_value = f_dist.sf(f_stat);

    (f_stat, p_value)
}

/// Performs an independent two-sample t-test to determine if two samples have different means.
///
/// # Returns
/// A tuple containing the t-statistic and the p-value.
#[must_use]
pub fn two_sample_t_test(
    sample1: &[f64],
    sample2: &[f64],
) -> (f64, f64) {
    let n1 = sample1.len() as f64;

    let n2 = sample2.len() as f64;

    let sample1_vec = sample1.to_vec();

    let sample2_vec = sample2.to_vec();

    let mean1 = mean(&sample1_vec);

    let mean2 = mean(&sample2_vec);

    // The pooled estimator weights *sample* (n - 1) variances by (n - 1).
    let (Some(var1), Some(var2)) = (
        variance_with_type(&sample1_vec, VarianceType::Sample),
        variance_with_type(&sample2_vec, VarianceType::Sample),
    ) else {
        return (f64::NAN, f64::NAN);
    };

    let s_p_sq = (n1 - 1.0).mul_add(var1, (n2 - 1.0) * var2) / (n1 + n2 - 2.0);

    let t_stat = (mean1 - mean2) / (s_p_sq * (1.0 / n1 + 1.0 / n2)).sqrt();

    let df = n1 + n2 - 2.0;

    let Ok(t_dist) = statrs::distribution::StudentsT::new(0.0, 1.0, df) else {
        return (f64::NAN, f64::NAN);
    };

    if t_stat.is_nan() {
        return (t_stat, f64::NAN);
    }
    let p_value = 2.0 * t_dist.sf(t_stat.abs());

    (t_stat, p_value)
}

/// Computes the Shannon entropy of a discrete probability distribution.
/// H(X) = -Σ p(x) * log2(p(x))
#[must_use]
pub fn shannon_entropy(probabilities: &[f64]) -> f64 {
    if probabilities.iter().any(|p| p.is_nan()) {
        return f64::NAN;
    }

    probabilities
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.log2())
        .sum()
}

/// Computes the geometric mean of a slice of positive data.
/// GM = (x1 * x2 * ... * xn)^(1/n)
#[must_use]
pub fn geometric_mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }

    // Use log to avoid overflow
    let log_sum: f64 = data.iter().map(|&x| x.ln()).sum();

    (log_sum / data.len() as f64).exp()
}

/// Computes the harmonic mean of a slice of positive data.
/// HM = n / (1/x1 + 1/x2 + ... + 1/xn)
#[must_use]
pub fn harmonic_mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }

    let reciprocal_sum: f64 = data.iter().map(|&x| 1.0 / x).sum();

    data.len() as f64 / reciprocal_sum
}

/// Computes the range (max - min) of a slice of data.
#[must_use]
pub fn range(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }

    let max_val = data.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let min_val = data.iter().copied().fold(f64::INFINITY, f64::min);

    max_val - min_val
}

/// Computes the Interquartile Range (IQR) = Q3 - Q1.
pub fn iqr(data: &mut [f64]) -> f64 {
    if data.len() < 4 || data.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }

    let mut data_container = Data::new(data);

    let q1 = data_container.percentile(25);

    let q3 = data_container.percentile(75);

    q3 - q1
}

/// Computes the z-scores (standard scores) for each data point.
/// z = (x - mean) / `std_dev`
#[must_use]
pub fn z_scores(data: &[f64]) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }

    let m = mean(data);

    let s = std_dev(data);

    if s == 0.0 || s.is_nan() {
        return vec![0.0; data.len()];
    }

    data.iter().map(|&x| (x - m) / s).collect()
}

/// Finds the mode (most frequent value) of a slice of data.
/// For continuous data, rounds to a specified number of decimal places.
/// Returns None if no mode exists (all values unique or empty).
#[must_use]
pub fn mode(
    data: &[f64],
    decimal_places: u32,
) -> Option<f64> {
    if data.is_empty() {
        return None;
    }

    let factor = 10_f64.powi(decimal_places as i32);

    let mut counts = std::collections::HashMap::new();

    for &val in data {
        let rounded = (val * factor).round() as i64;

        *counts.entry(rounded).or_insert(0) += 1;
    }

    let max_count = *counts.values().max()?;

    if max_count == 1 {
        return None; // No mode if all values are unique
    }

    let mode_key = counts.into_iter().find(|(_, v)| *v == max_count)?.0;

    Some(mode_key as f64 / factor)
}

/// Performs Welch's t-test for two samples with potentially unequal variances.
/// Returns (t-statistic, p-value).
#[must_use]
pub fn welch_t_test(
    sample1: &[f64],
    sample2: &[f64],
) -> (f64, f64) {
    let n1 = sample1.len() as f64;

    let n2 = sample2.len() as f64;

    if n1 < 2.0 || n2 < 2.0 {
        return (f64::NAN, f64::NAN);
    }

    let mean1 = mean(sample1);

    let mean2 = mean(sample2);

    let var1 = variance(sample1);

    let var2 = variance(sample2);

    // Correct for sample variance (using n-1)
    let s1_sq = var1 * n1 / (n1 - 1.0);

    let s2_sq = var2 * n2 / (n2 - 1.0);

    let se = (s1_sq / n1 + s2_sq / n2).sqrt();

    if se == 0.0 {
        return (f64::NAN, f64::NAN);
    }

    let t_stat = (mean1 - mean2) / se;

    // Welch-Satterthwaite degrees of freedom
    let num = (s1_sq / n1 + s2_sq / n2).powi(2);

    let denom = (s1_sq / n1).powi(2) / (n1 - 1.0) + (s2_sq / n2).powi(2) / (n2 - 1.0);

    let df = num / denom;

    let Ok(t_dist) = statrs::distribution::StudentsT::new(0.0, 1.0, df) else {
        return (f64::NAN, f64::NAN);
    };

    if t_stat.is_nan() {
        return (t_stat, f64::NAN);
    }
    let p_value = 2.0 * t_dist.sf(t_stat.abs());

    (t_stat, p_value)
}

/// Performs a chi-squared goodness-of-fit test.
/// Tests if observed frequencies match expected frequencies.
/// Returns (chi-squared statistic, p-value).
#[must_use]
pub fn chi_squared_test(
    observed: &[f64],
    expected: &[f64],
) -> (f64, f64) {
    if observed.len() != expected.len() || observed.is_empty() {
        return (f64::NAN, f64::NAN);
    }

    let chi_sq: f64 = observed
        .iter()
        .zip(expected.iter())
        .map(|(&o, &e)| {
            if e == 0.0 {
                0.0
            } else {
                (o - e).powi(2) / e
            }
        })
        .sum();

    let df = (observed.len() - 1) as f64;

    if df <= 0.0 {
        return (chi_sq, f64::NAN);
    }

    let Ok(chi_dist) = statrs::distribution::ChiSquared::new(df) else {
        return (chi_sq, f64::NAN);
    };

    if chi_sq.is_nan() {
        return (chi_sq, f64::NAN);
    }
    let p_value = chi_dist.sf(chi_sq);

    (chi_sq, p_value)
}

/// Computes the coefficient of variation (CV = `std_dev` / mean).
/// Useful for comparing variability across datasets with different means.
#[must_use]
pub fn coefficient_of_variation(data: &[f64]) -> f64 {
    let m = mean(data);

    if m == 0.0 {
        return f64::NAN;
    }

    std_dev(data) / m
}

/// Computes the sample standard error of the mean (SEM = `std_dev` / sqrt(n)).
#[must_use]
pub fn standard_error(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }

    std_dev(data) / (data.len() as f64).sqrt()
}
