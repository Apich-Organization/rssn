//! # Numerical Statistics
//!
//! This module provides numerical statistical functions, leveraging the `statrs` crate
//! for robust implementations. It includes descriptive statistics (mean, variance, median,
//! percentiles, covariance, correlation), probability distributions (Normal, Uniform, Binomial,
//! Poisson, Exponential, Gamma), and statistical inference methods (ANOVA, t-tests).

use statrs::distribution::Binomial;
use statrs::distribution::Continuous;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Discrete;
use statrs::distribution::DiscreteCDF;
use statrs::distribution::Normal;
use statrs::distribution::Uniform;
use statrs::statistics::Data;
use statrs::statistics::Distribution;
use statrs::statistics::Max;
use statrs::statistics::Min;
use statrs::statistics::OrderStatistics;

/// Computes the mean of a slice of data.
///
/// # Arguments
/// * `data` - A unmutable slice of `f64` data points.
///
/// # Returns
/// The mean of the data as an `f64`.
#[must_use]

pub fn mean(data : &[f64]) -> f64 {

    if data.is_empty() {

        return 0.0;
    }

    data.iter()
        .sum::<f64>()
        / (data.len() as f64)
}

/// Computes the variance of a slice of f64 values.
///
/// # Arguments
/// * `data` - A unmutable slice of `f64` data points.
///
/// # Returns
/// The variance of the data as an `f64`.
#[must_use]

pub fn variance(data : &[f64]) -> f64 {

    if data.is_empty() {

        return 0.0;
    }

    let mean_val = mean(data);

    data.iter()
        .map(|&val| (val - mean_val).powi(2))
        .sum::<f64>()
        / (data.len() as f64)
}

/// Computes the standard deviation of a slice of data.
///
/// # Arguments
/// * `data` - A slice of `f64` data points.
///
/// # Returns
/// The standard deviation of the data as an `f64`.
#[must_use]

pub fn std_dev(data : &[f64]) -> f64 {

    let data_vec : Vec<f64> = data.to_vec();

    let _data_container = Data::new(data_vec);

    let data_vec : Vec<f64> = data.to_vec();

    let data_container = Data::new(data_vec);

    data_container
        .std_dev()
        .unwrap_or(f64::NAN)
}

/// Computes the median of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
///
/// # Returns
/// The median of the data as an `f64`.

pub fn median(data : &mut [f64]) -> f64 {

    let mut data_container = Data::new(data);

    data_container.median()
}

/// Computes the p-th percentile of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
/// * `p` - The desired percentile (e.g., 50.0 for median).
///
/// # Returns
/// The p-th percentile of the data as an `f64`.

pub fn percentile(
    data : &mut [f64],
    p : f64,
) -> f64 {

    let mut data_container = Data::new(data);

    data_container.percentile(p as usize)
}

/// Computes the covariance of two slices of data.
///
/// # Arguments
/// * `data1` - The first slice of `f64` data points.
/// * `data2` - The second slice of `f64` data points.
///
/// # Returns
/// The covariance of the two data sets as an `f64`.
#[must_use]

pub fn covariance(
    data1 : &[f64],
    data2 : &[f64],
) -> f64 {

    let data1_vec = data1.to_vec();

    let data2_vec = data2.to_vec();

    let mean1 = mean(&data1_vec);

    let mean2 = mean(&data2_vec);

    let n = data1.len() as f64;

    data1
        .iter()
        .zip(data2.iter())
        .map(|(&x, &y)| (x - mean1) * (y - mean2))
        .sum::<f64>()
        / (n - 1.0)
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
    data1 : &[f64],
    data2 : &[f64],
) -> f64 {

    let cov = covariance(data1, data2);

    let std_dev1 = std_dev(data1);

    let std_dev2 = std_dev(data2);

    cov / (std_dev1 * std_dev2)
}

/// Represents a Normal (Gaussian) distribution.

pub struct NormalDist(Normal);

impl NormalDist {
    /// Creates a new `NormalDist` instance.
    ///
    /// # Arguments
    /// * `mean` - The mean `μ` of the distribution.
    /// * `std_dev` - The standard deviation `σ` of the distribution.
    ///
    /// # Returns
    /// A `Result` containing the `NormalDist` instance, or an error string if parameters are invalid.
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
        x : f64,
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
        x : f64,
    ) -> f64 {

        self.0.cdf(x)
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

    pub fn new(
        min : f64,
        max : f64,
    ) -> Result<Self, String> {

        Uniform::new(min, max)
            .map(UniformDist)
            .map_err(|e| e.to_string())
    }

    #[must_use]

    pub fn pdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.pdf(x)
    }

    #[must_use]

    pub fn cdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.cdf(x)
    }
}

/// Represents a Binomial distribution.

pub struct BinomialDist(Binomial);

impl BinomialDist {
    pub fn new(
        n : u64,
        p : f64,
    ) -> Result<Self, String> {

        Binomial::new(p, n)
            .map(BinomialDist)
            .map_err(|e| e.to_string())
    }

    #[must_use]

    pub fn pmf(
        &self,
        k : u64,
    ) -> f64 {

        self.0.pmf(k)
    }

    #[must_use]

    pub fn cdf(
        &self,
        k : u64,
    ) -> f64 {

        self.0.cdf(k)
    }
}

/// Performs a simple linear regression on a set of 2D points.
/// Returns the slope (b1) and intercept (b0) of the best-fit line y = b0 + b1*x.
#[must_use]

pub fn simple_linear_regression(data : &[(f64, f64)]) -> (f64, f64) {

    if data.is_empty() {

        return (f64::NAN, f64::NAN);
    }

    let (xs, ys) : (Vec<_>, Vec<_>) = data
        .iter()
        .copied()
        .unzip();

    let mean_x = mean(&xs);

    let mean_y = mean(&ys);

    // Compute slope using consistent formula: b1 = Σ(xi - x̄)(yi - ȳ) / Σ(xi - x̄)²
    let numerator : f64 = xs
        .iter()
        .zip(ys.iter())
        .map(|(&x, &y)| (x - mean_x) * (y - mean_y))
        .sum();

    let denominator : f64 = xs
        .iter()
        .map(|&x| (x - mean_x).powi(2))
        .sum();

    if denominator == 0.0 {

        return (f64::NAN, mean_y);
    }

    let b1 = numerator / denominator;

    let b0 = mean_y - b1 * mean_x;

    (b1, b0)
}

/// Computes the minimum value of a slice of data.

pub fn min(data : &mut [f64]) -> f64 {

    let data_container = Data::new(data);

    data_container.min()
}

/// Computes the maximum value of a slice of data.

pub fn max(data : &mut [f64]) -> f64 {

    let data_container = Data::new(data);

    data_container.max()
}

/// Computes the skewness of a slice of data.

pub fn skewness(data : &mut [f64]) -> f64 {

    let data_container = Data::new(data);

    data_container
        .skewness()
        .unwrap_or(f64::NAN)
}

/// Computes the sample kurtosis (Fisher's g2) of a slice of data.

pub fn kurtosis(data : &mut [f64]) -> f64 {

    let n = data.len() as f64;

    if n < 4.0 {

        return f64::NAN;
    }

    let mean = mean(data);

    let m2 = data
        .iter()
        .map(|&x| (x - mean).powi(2))
        .sum::<f64>()
        / n;

    let m4 = data
        .iter()
        .map(|&x| (x - mean).powi(4))
        .sum::<f64>()
        / n;

    if m2 == 0.0 {

        return 0.0;
    }

    let g2 = m4 / m2.powi(2) - 3.0;

    let term1 = n.mul_add(n, -1.0) / ((n - 2.0) * (n - 3.0));

    let term2 = (g2 + 3.0) - 3.0 * (n - 1.0).powi(2) / ((n - 2.0) * (n - 3.0));

    term1 * term2
}

/// Represents a Poisson distribution.

pub struct PoissonDist(statrs::distribution::Poisson);

impl PoissonDist {
    pub fn new(rate : f64) -> Result<Self, String> {

        statrs::distribution::Poisson::new(rate)
            .map(PoissonDist)
            .map_err(|e| e.to_string())
    }

    #[must_use]

    pub fn pmf(
        &self,
        k : u64,
    ) -> f64 {

        self.0.pmf(k)
    }

    #[must_use]

    pub fn cdf(
        &self,
        k : u64,
    ) -> f64 {

        self.0.cdf(k)
    }
}

/// Represents an Exponential distribution.

pub struct ExponentialDist(statrs::distribution::Exp);

impl ExponentialDist {
    pub fn new(rate : f64) -> Result<Self, String> {

        statrs::distribution::Exp::new(rate)
            .map(ExponentialDist)
            .map_err(|e| e.to_string())
    }

    #[must_use]

    pub fn pdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.pdf(x)
    }

    #[must_use]

    pub fn cdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.cdf(x)
    }
}

/// Represents a Gamma distribution.

pub struct GammaDist(statrs::distribution::Gamma);

impl GammaDist {
    pub fn new(
        shape : f64,
        rate : f64,
    ) -> Result<Self, String> {

        statrs::distribution::Gamma::new(shape, rate)
            .map(GammaDist)
            .map_err(|e| e.to_string())
    }

    #[must_use]

    pub fn pdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.pdf(x)
    }

    #[must_use]

    pub fn cdf(
        &self,
        x : f64,
    ) -> f64 {

        self.0.cdf(x)
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

pub fn one_way_anova(groups : &mut [&mut [f64]]) -> (f64, f64) {

    let k = groups.len() as f64;

    if k < 2.0 {

        return (f64::NAN, f64::NAN);
    }

    let all_data : Vec<f64> = groups
        .iter()
        .flat_map(|g| g.iter())
        .copied()
        .collect();

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

        ss_within += group
            .iter()
            .map(|&x| (x - mean_group).powi(2))
            .sum::<f64>();
    }

    let df_within = n_total - k;

    let ms_within = ss_within / df_within;

    if ms_within == 0.0 {

        return (f64::INFINITY, 0.0);
    }

    let f_stat = ms_between / ms_within;

    let f_dist = match statrs::distribution::FisherSnedecor::new(
        df_between,
        df_within,
    ) {
        | Ok(dist) => dist,
        | Err(_) => return (f64::NAN, f64::NAN),
    };

    let p_value = 1.0 - f_dist.cdf(f_stat);

    (f_stat, p_value)
}

/// Performs an independent two-sample t-test to determine if two samples have different means.
///
/// # Returns
/// A tuple containing the t-statistic and the p-value.
#[must_use]

pub fn two_sample_t_test(
    sample1 : &[f64],
    sample2 : &[f64],
) -> (f64, f64) {

    let n1 = sample1.len() as f64;

    let n2 = sample2.len() as f64;

    let sample1_vec = sample1.to_vec();

    let sample2_vec = sample2.to_vec();

    let mean1 = mean(&sample1_vec);

    let mean2 = mean(&sample2_vec);

    let var1 = variance(&sample1_vec);

    let var2 = variance(&sample2_vec);

    let s_p_sq = (n1 - 1.0).mul_add(
        var1,
        (n2 - 1.0) * var2,
    ) / (n1 + n2 - 2.0);

    let t_stat = (mean1 - mean2) / (s_p_sq * (1.0 / n1 + 1.0 / n2)).sqrt();

    let df = n1 + n2 - 2.0;

    let t_dist = match statrs::distribution::StudentsT::new(0.0, 1.0, df) {
        | Ok(dist) => dist,
        | Err(_) => return (f64::NAN, f64::NAN),
    };

    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    (t_stat, p_value)
}

/// Computes the Shannon entropy of a discrete probability distribution.
/// H(X) = -Σ p(x) * log2(p(x))
#[must_use]

pub fn shannon_entropy(probabilities : &[f64]) -> f64 {

    probabilities
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.log2())
        .sum()
}

/// Computes the geometric mean of a slice of positive data.
/// GM = (x1 * x2 * ... * xn)^(1/n)
#[must_use]

pub fn geometric_mean(data : &[f64]) -> f64 {

    if data.is_empty() {

        return f64::NAN;
    }

    // Use log to avoid overflow
    let log_sum : f64 = data
        .iter()
        .map(|&x| x.ln())
        .sum();

    (log_sum / data.len() as f64).exp()
}

/// Computes the harmonic mean of a slice of positive data.
/// HM = n / (1/x1 + 1/x2 + ... + 1/xn)
#[must_use]

pub fn harmonic_mean(data : &[f64]) -> f64 {

    if data.is_empty() {

        return f64::NAN;
    }

    let reciprocal_sum : f64 = data
        .iter()
        .map(|&x| 1.0 / x)
        .sum();

    data.len() as f64 / reciprocal_sum
}

/// Computes the range (max - min) of a slice of data.
#[must_use]

pub fn range(data : &[f64]) -> f64 {

    if data.is_empty() {

        return f64::NAN;
    }

    let max_val = data
        .iter()
        .cloned()
        .fold(
            f64::NEG_INFINITY,
            f64::max,
        );

    let min_val = data
        .iter()
        .cloned()
        .fold(
            f64::INFINITY,
            f64::min,
        );

    max_val - min_val
}

/// Computes the Interquartile Range (IQR) = Q3 - Q1.

pub fn iqr(data : &mut [f64]) -> f64 {

    if data.len() < 4 {

        return f64::NAN;
    }

    let mut data_container = Data::new(data);

    let q1 = data_container.percentile(25);

    let q3 = data_container.percentile(75);

    q3 - q1
}

/// Computes the z-scores (standard scores) for each data point.
/// z = (x - mean) / std_dev
#[must_use]

pub fn z_scores(data : &[f64]) -> Vec<f64> {

    if data.is_empty() {

        return vec![];
    }

    let m = mean(data);

    let s = std_dev(data);

    if s == 0.0 || s.is_nan() {

        return vec![0.0; data.len()];
    }

    data.iter()
        .map(|&x| (x - m) / s)
        .collect()
}

/// Finds the mode (most frequent value) of a slice of data.
/// For continuous data, rounds to a specified number of decimal places.
/// Returns None if no mode exists (all values unique or empty).
#[must_use]

pub fn mode(
    data : &[f64],
    decimal_places : u32,
) -> Option<f64> {

    if data.is_empty() {

        return None;
    }

    let factor = 10_f64.powi(decimal_places as i32);

    let mut counts = std::collections::HashMap::new();

    for &val in data {

        let rounded = (val * factor).round() as i64;

        *counts
            .entry(rounded)
            .or_insert(0) += 1;
    }

    let max_count = *counts
        .values()
        .max()?;

    if max_count == 1 {

        return None; // No mode if all values are unique
    }

    let mode_key = counts
        .into_iter()
        .find(|(_, v)| *v == max_count)?
        .0;

    Some(mode_key as f64 / factor)
}

/// Performs Welch's t-test for two samples with potentially unequal variances.
/// Returns (t-statistic, p-value).
#[must_use]

pub fn welch_t_test(
    sample1 : &[f64],
    sample2 : &[f64],
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

    let t_dist = match statrs::distribution::StudentsT::new(0.0, 1.0, df) {
        | Ok(dist) => dist,
        | Err(_) => return (f64::NAN, f64::NAN),
    };

    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    (t_stat, p_value)
}

/// Performs a chi-squared goodness-of-fit test.
/// Tests if observed frequencies match expected frequencies.
/// Returns (chi-squared statistic, p-value).
#[must_use]

pub fn chi_squared_test(
    observed : &[f64],
    expected : &[f64],
) -> (f64, f64) {

    if observed.len() != expected.len() || observed.is_empty() {

        return (f64::NAN, f64::NAN);
    }

    let chi_sq : f64 = observed
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

    let chi_dist = match statrs::distribution::ChiSquared::new(df) {
        | Ok(dist) => dist,
        | Err(_) => return (chi_sq, f64::NAN),
    };

    let p_value = 1.0 - chi_dist.cdf(chi_sq);

    (chi_sq, p_value)
}

/// Computes the coefficient of variation (CV = std_dev / mean).
/// Useful for comparing variability across datasets with different means.
#[must_use]

pub fn coefficient_of_variation(data : &[f64]) -> f64 {

    let m = mean(data);

    if m == 0.0 {

        return f64::NAN;
    }

    std_dev(data) / m
}

/// Computes the sample standard error of the mean (SEM = std_dev / sqrt(n)).
#[must_use]

pub fn standard_error(data : &[f64]) -> f64 {

    if data.is_empty() {

        return f64::NAN;
    }

    std_dev(data) / (data.len() as f64).sqrt()
}
