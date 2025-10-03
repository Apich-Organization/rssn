//! # Numerical Statistics
//!
//! This module provides numerical statistical functions, leveraging the `statrs` crate
//! for robust implementations. It includes descriptive statistics (mean, variance, median,
//! percentiles, covariance, correlation), probability distributions (Normal, Uniform, Binomial,
//! Poisson, Exponential, Gamma), and statistical inference methods (ANOVA, t-tests).

use statrs::distribution::{Binomial, Continuous, Discrete, Normal, Uniform};
//use statrs::statistics::{Data, Median, OrderStatistics, Statistics};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::DiscreteCDF;
use statrs::statistics::Distribution;
use statrs::statistics::{Data, Max, Min, OrderStatistics};

// =====================================================================================
// region: Descriptive Statistics
// =====================================================================================

/// Computes the mean of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
///
/// # Returns
/// The mean of the data as an `f64`.
pub fn mean(data: &mut [f64]) -> f64 {
    let data_container = Data::new(data);
    data_container.mean().unwrap_or(f64::NAN)
}

/// Computes the variance of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
///
/// # Returns
/// The variance of the data as an `f64`.
pub fn variance(data: &mut [f64]) -> f64 {
    let data_container = Data::new(data);
    data_container.variance().unwrap_or(f64::NAN)
}

/// Computes the standard deviation of a slice of data.
///
/// # Arguments
/// * `data` - A slice of `f64` data points.
///
/// # Returns
/// The standard deviation of the data as an `f64`.
pub fn std_dev(data: &[f64]) -> f64 {
    // Assuming 'data' is the immutable &[f64] being passed in:
    let data_vec: Vec<f64> = data.to_vec(); // Create an owned, mutable copy of the data

    // Pass the owned vector to Data::new. Vec<f64> satisfies the trait bound.
    let _data_container = Data::new(data_vec); // ✅ FIXES THE ERROR

    // Assuming 'data' is the input &[f64]
    // 1. Clone the immutable slice into an owned, mutable vector (Vec<f64>).
    let data_vec: Vec<f64> = data.to_vec();

    // 2. Create the Data container using the Vec<f64>.
    // The type is now Data<Vec<f64>>, which satisfies the AsMut bound.
    let data_container = Data::new(data_vec);

    // 3. Now std_dev() will work because the required trait bounds are satisfied.

    data_container.std_dev().unwrap_or(f64::NAN)
}

/// Computes the median of a slice of data.
///
/// # Arguments
/// * `data` - A mutable slice of `f64` data points.
///
/// # Returns
/// The median of the data as an `f64`.
pub fn median(data: &mut [f64]) -> f64 {
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
pub fn percentile(data: &mut [f64], p: f64) -> f64 {
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
pub fn covariance(data1: &[f64], data2: &[f64]) -> f64 {
    let mut data1_vec = data1.to_vec();
    let mut data2_vec = data2.to_vec();
    let mean1 = mean(&mut data1_vec);
    let mean2 = mean(&mut data2_vec);
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
pub fn correlation(data1: &[f64], data2: &[f64]) -> f64 {
    let cov = covariance(data1, data2);
    let std_dev1 = std_dev(data1);
    let std_dev2 = std_dev(data2);
    cov / (std_dev1 * std_dev2)
}

// =====================================================================================
// region: Probability Distributions
// =====================================================================================

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
    pub fn pdf(&self, x: f64) -> f64 {
        self.0.pdf(x)
    }
    /// Returns the cumulative distribution function (CDF) value at `x`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the CDF.
    ///
    /// # Returns
    /// The CDF value as an `f64`.
    pub fn cdf(&self, x: f64) -> f64 {
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
    pub fn new(min: f64, max: f64) -> Result<Self, String> {
        Uniform::new(min, max)
            .map(UniformDist)
            .map_err(|e| e.to_string())
    }
    pub fn pdf(&self, x: f64) -> f64 {
        self.0.pdf(x)
    }
    pub fn cdf(&self, x: f64) -> f64 {
        self.0.cdf(x)
    }
}

/// Represents a Binomial distribution.
pub struct BinomialDist(Binomial);

impl BinomialDist {
    pub fn new(n: u64, p: f64) -> Result<Self, String> {
        Binomial::new(p, n)
            .map(BinomialDist)
            .map_err(|e| e.to_string())
    }
    pub fn pmf(&self, k: u64) -> f64 {
        self.0.pmf(k)
    }
    pub fn cdf(&self, k: u64) -> f64 {
        self.0.cdf(k)
    }
}

// =====================================================================================
// region: Regression
// =====================================================================================

/// Performs a simple linear regression on a set of 2D points.
/// Returns the slope (b1) and intercept (b0) of the best-fit line y = b0 + b1*x.
pub fn simple_linear_regression(data: &[(f64, f64)]) -> (f64, f64) {
    let (mut xs, mut ys): (Vec<_>, Vec<_>) = data.iter().cloned().unzip();
    let mean_x = mean(&mut xs);
    let mean_y = mean(&mut ys);
    let cov_xy = covariance(&mut xs, &mut ys);
    let var_x = variance(&mut xs);

    let b1 = cov_xy / var_x;
    let b0 = mean_y - b1 * mean_x;

    (b1, b0)
}

// =====================================================================================
// region: More Descriptive Statistics
// =====================================================================================

/// Computes the minimum value of a slice of data.
pub fn min(data: &mut [f64]) -> f64 {
    let data_container = Data::new(data);
    data_container.min()
}

/// Computes the maximum value of a slice of data.
pub fn max(data: &mut [f64]) -> f64 {
    let data_container = Data::new(data);
    data_container.max()
}

/// Computes the skewness of a slice of data.
pub fn skewness(data: &mut [f64]) -> f64 {
    let data_container = Data::new(data);
    data_container.skewness().unwrap_or(f64::NAN)
}

/// Computes the sample kurtosis (Fisher's g2) of a slice of data.
pub fn kurtosis(data: &mut [f64]) -> f64 {
    let n = data.len() as f64;
    if n < 4.0 {
        return f64::NAN; // Kurtosis is not well-defined for small samples
    }
    let mean = mean(data);
    let m2 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let m4 = data.iter().map(|&x| (x - mean).powi(4)).sum::<f64>() / n;

    if m2 == 0.0 {
        return 0.0;
    }

    // Fisher's kurtosis
    let g2 = m4 / m2.powi(2) - 3.0;

    // Apply sample size correction for an unbiased estimator
    let _correction = (n + 1.0) * (n - 1.0) / ((n - 2.0) * (n - 3.0));
    let term1 = (n * n - 1.0) / ((n - 2.0) * (n - 3.0));
    let term2 = (g2 + 3.0) - 3.0 * (n - 1.0).powi(2) / ((n - 2.0) * (n - 3.0));
    term1 * term2
}

// =====================================================================================
// region: More Probability Distributions
// =====================================================================================

/// Represents a Poisson distribution.
pub struct PoissonDist(statrs::distribution::Poisson);

impl PoissonDist {
    pub fn new(rate: f64) -> Result<Self, String> {
        statrs::distribution::Poisson::new(rate)
            .map(PoissonDist)
            .map_err(|e| e.to_string())
    }
    pub fn pmf(&self, k: u64) -> f64 {
        self.0.pmf(k)
    }
    pub fn cdf(&self, k: u64) -> f64 {
        self.0.cdf(k)
    }
}

/// Represents an Exponential distribution.
pub struct ExponentialDist(statrs::distribution::Exp);

impl ExponentialDist {
    pub fn new(rate: f64) -> Result<Self, String> {
        statrs::distribution::Exp::new(rate)
            .map(ExponentialDist)
            .map_err(|e| e.to_string())
    }
    pub fn pdf(&self, x: f64) -> f64 {
        self.0.pdf(x)
    }
    pub fn cdf(&self, x: f64) -> f64 {
        self.0.cdf(x)
    }
}

/// Represents a Gamma distribution.
pub struct GammaDist(statrs::distribution::Gamma);

impl GammaDist {
    pub fn new(shape: f64, rate: f64) -> Result<Self, String> {
        statrs::distribution::Gamma::new(shape, rate)
            .map(GammaDist)
            .map_err(|e| e.to_string())
    }
    pub fn pdf(&self, x: f64) -> f64 {
        self.0.pdf(x)
    }
    pub fn cdf(&self, x: f64) -> f64 {
        self.0.cdf(x)
    }
}

// =====================================================================================
// region: Advanced Inference
// =====================================================================================

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

    //let all_data: Vec<f64> = groups.iter().flat_map(|g| *g).copied().collect();
    //let all_data: Vec<f64> = groups.iter().flat_map(|g| *g).copied().collect();
    // Original code (Error):
    // let all_data: Vec<f64> = groups.iter().flat_map(|g| *g).copied().collect();

    // Corrected Code:
    let all_data: Vec<f64> = groups
        .iter()
        // g is &&mut [f64]. *g is &mut [f64]. .iter() on the slice yields &mut f64.
        // The inner closure *x creates a copy of the f64 which is then referenced.
        .flat_map(|g| g.iter())
        .copied()
        .collect();
    let n_total = all_data.len() as f64;
    let grand_mean = mean(&mut all_data.clone());

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
    let f_dist = statrs::distribution::FisherSnedecor::new(df_between, df_within).unwrap();
    let p_value = 1.0 - f_dist.cdf(f_stat);

    (f_stat, p_value)
}

// =====================================================================================
// region: Statistical Inference
// =====================================================================================

/// Performs an independent two-sample t-test to determine if two samples have different means.
///
/// # Returns
/// A tuple containing the t-statistic and the p-value.
pub fn two_sample_t_test(sample1: &[f64], sample2: &[f64]) -> (f64, f64) {
    let n1 = sample1.len() as f64;
    let n2 = sample2.len() as f64;
    let mut sample1_vec = sample1.to_vec();
    let mut sample2_vec = sample2.to_vec();
    let mean1 = mean(&mut sample1_vec);
    let mean2 = mean(&mut sample2_vec);
    let var1 = variance(&mut sample1_vec);
    let var2 = variance(&mut sample2_vec);

    let s_p_sq = ((n1 - 1.0) * var1 + (n2 - 1.0) * var2) / (n1 + n2 - 2.0);
    let t_stat = (mean1 - mean2) / (s_p_sq * (1.0 / n1 + 1.0 / n2)).sqrt();

    let df = n1 + n2 - 2.0;
    let t_dist = statrs::distribution::StudentsT::new(0.0, 1.0, df).unwrap();

    // Two-tailed test
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    (t_stat, p_value)
}

// =====================================================================================
// region: Information Theory
// =====================================================================================

/// Computes the Shannon entropy of a discrete probability distribution.
/// H(X) = -Σ p(x) * log2(p(x))
pub fn shannon_entropy(probabilities: &[f64]) -> f64 {
    probabilities
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.log2())
        .sum()
}
