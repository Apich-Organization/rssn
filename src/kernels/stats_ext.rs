//! Hypothesis testing, p-value adjustment, regression (OLS, ridge, lasso, logistic), PCA,
//! rank correlations (Spearman, Kendall), KDE, and bootstrap / jackknife resampling.

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::descriptive;
use crate::kernels::distributions::Distribution;
use crate::kernels::numeric_common::scaled_difference;
use crate::kernels::par;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::robust;
use crate::kernels::special::ln_beta_numerical;
use crate::kernels::special::regularized_beta;
use std::borrow::Cow;
use std::cmp::Ordering;
use std::f64::consts::PI;

/// Error returned by extended statistics kernels.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StatsExtError {
    /// Input slice was empty.
    Empty,
    /// Mismatched dimensions or lengths.
    ShapeMismatch(String),
    /// Invalid domain or insufficient observations.
    Domain(String),
    /// Linear system was singular.
    Singular,
    /// Iterative solver did not converge.
    NotConverged,
}

impl std::fmt::Display for StatsExtError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("empty input"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
            | Self::Singular => f.write_str("singular matrix"),
            | Self::NotConverged => f.write_str("did not converge"),
        }
    }
}

impl std::error::Error for StatsExtError {}

fn map_fit_error(error: DenseError) -> StatsExtError {
    match error {
        | DenseError::Singular => StatsExtError::Singular,
        | DenseError::Cancelled | DenseError::NoConvergence => StatsExtError::NotConverged,
        | DenseError::Shape => {
            StatsExtError::ShapeMismatch("inconsistent regression matrix".into())
        },
        | DenseError::NonFinite => {
            StatsExtError::Domain("regression produced non-finite data".into())
        },
    }
}

/// The outcome of a statistical hypothesis test.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TestResult {
    /// The test statistic (e.g. `t`, `chi2`, `F`, `D`, or `U`).
    pub statistic: f64,
    /// The p-value under the null hypothesis.
    pub p_value: f64,
    /// Degrees of freedom, when applicable.
    pub df: Option<f64>,
}

fn sorted_p_values(probabilities: &[f64]) -> Result<Vec<(usize, f64)>, StatsExtError> {
    if probabilities.iter().any(|p| !(0.0..=1.0).contains(p)) {
        return Err(StatsExtError::Domain(
            "p-values must be finite and in [0, 1]".into(),
        ));
    }
    let mut order: Vec<_> = probabilities.iter().copied().enumerate().collect();
    order.sort_unstable_by(|(_, a), (_, b)| a.total_cmp(b));
    Ok(order)
}

fn adjust_p_step_up(
    probabilities: &[f64],
    dependence_factor: f64,
) -> Result<Vec<f64>, StatsExtError> {
    let order = sorted_p_values(probabilities)?;
    let mut adjusted = vec![0.0; order.len()];
    let tests = order.len() as f64;
    let mut upper: f64 = 1.0;
    for (rank, (index, probability)) in order.into_iter().enumerate().rev() {
        let multiplier = dependence_factor * (tests / (rank + 1) as f64);
        upper = upper.min(multiplier * probability);
        adjusted[index] = upper;
    }
    Ok(adjusted)
}

/// Benjamini–Hochberg adjusted p-values for false discovery rate control.
///
/// The guarantee holds under independence or positive regression dependence
/// of the null tests.
/// Results retain the original order. An empty input yields an empty vector.
///
/// # Errors
/// [`StatsExtError::Domain`] if any p-value is outside `[0, 1]` or non-finite.
pub fn adjust_p_bh(probabilities: &[f64]) -> Result<Vec<f64>, StatsExtError> {
    adjust_p_step_up(probabilities, 1.0)
}

/// Benjamini–Yekutieli adjusted p-values for false discovery rate control
/// under arbitrary dependence. Results retain the original order.
/// An empty input yields an empty vector.
///
/// # Errors
/// [`StatsExtError::Domain`] if any p-value is outside `[0, 1]` or non-finite.
pub fn adjust_p_by(probabilities: &[f64]) -> Result<Vec<f64>, StatsExtError> {
    let harmonic = (1..=probabilities.len())
        .map(|rank| 1.0 / rank as f64)
        .sum();
    adjust_p_step_up(probabilities, harmonic)
}

/// Holm step-down adjusted p-values for family-wise error rate control.
/// Results retain the original order. An empty input yields an empty vector.
///
/// # Errors
/// [`StatsExtError::Domain`] if any p-value is outside `[0, 1]` or non-finite.
pub fn adjust_p_holm(probabilities: &[f64]) -> Result<Vec<f64>, StatsExtError> {
    let order = sorted_p_values(probabilities)?;
    let mut adjusted = vec![0.0; order.len()];
    let tests = order.len();
    let mut lower: f64 = 0.0;
    for (rank, (index, probability)) in order.into_iter().enumerate() {
        lower = lower.max((tests - rank) as f64 * probability).min(1.0);
        adjusted[index] = lower;
    }
    Ok(adjusted)
}

/// The result of a linear regression (`ols`, `ridge`, or `lasso`).
#[derive(Clone, Debug, PartialEq)]
pub struct LinearRegressionResult {
    /// Fitted slope coefficients (length `cols`).
    pub coefficients: Vec<f64>,
    /// Fitted intercept (0.0 when `fit_intercept` is `false`).
    pub intercept: f64,
    /// Coefficient of determination `R²`.
    pub r_squared: f64,
    /// Residuals `y - y_hat` (length `rows`).
    pub residuals: Vec<f64>,
}

/// The result of a binary logistic regression.
#[derive(Clone, Debug, PartialEq)]
pub struct LogisticRegressionResult {
    /// Fitted feature weights (length `cols`).
    pub coefficients: Vec<f64>,
    /// Fitted intercept (0.0 when `fit_intercept` is `false`).
    pub intercept: f64,
    /// Final log-likelihood.
    pub log_likelihood: f64,
    /// Number of Newton-Raphson iterations performed.
    pub iterations: usize,
}

/// The result of Principal Component Analysis.
#[derive(Clone, Debug, PartialEq)]
pub struct PcaResult {
    /// Column means subtracted before decomposition (length `cols`).
    pub mean: Vec<f64>,
    /// Top `n_components` principal axes, row-major `n_components × cols`.
    pub components: Vec<f64>,
    /// Variance explained by each retained component.
    pub explained_variance: Vec<f64>,
    /// Fraction of total variance explained by each retained component.
    pub explained_variance_ratio: Vec<f64>,
    /// Projected data, row-major `rows × n_components`.
    pub transformed: Vec<f64>,
    /// Number of retained components.
    pub n_components: usize,
}

/// The result of bootstrap resampling.
#[derive(Clone, Debug, PartialEq)]
pub struct BootstrapResult {
    /// Statistic evaluated on the original sample.
    pub estimate: f64,
    /// Mean of the bootstrap replicates.
    pub mean: f64,
    /// Bootstrap standard error (sample std dev of replicates).
    pub std_error: f64,
    /// Lower 95% percentile confidence bound (2.5th percentile).
    pub ci_low: f64,
    /// Upper 95% percentile confidence bound (97.5th percentile).
    pub ci_high: f64,
    /// All bootstrap replicate values.
    pub replicates: Vec<f64>,
}

/// The result of leave-one-out jackknife resampling.
#[derive(Clone, Debug, PartialEq)]
pub struct JackknifeResult {
    /// Statistic evaluated on the original sample.
    pub estimate: f64,
    /// Jackknife bias estimate `(n - 1) * (mean_jack - estimate)`.
    pub bias: f64,
    /// Jackknife standard error.
    pub std_error: f64,
    /// Leave-one-out estimates `θ_{(-i)}`.
    pub jackknife_estimates: Vec<f64>,
}

fn finite(values: &[f64]) -> Result<(), StatsExtError> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err(StatsExtError::Domain(
            "observations must be finite".to_string(),
        ));
    }
    Ok(())
}

// Moments in unitless coordinates avoid overflowing the sum or the squares of
// finite measurements. Keep the scale separate when the final statistic is
// dimensionless or only needs the standard deviation.
struct ScaledMoments {
    mean: f64,
    variance: f64,
    scale: f64,
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|value| value.abs()).fold(0.0, f64::max)
}

fn scaled_moments(values: &[f64]) -> Result<ScaledMoments, StatsExtError> {
    finite(values)?;
    let n = values.len();
    let first = *values.first().ok_or(StatsExtError::Empty)?;
    let scale = max_abs(values).max(f64::MIN_POSITIVE);
    let anchor = first / scale;
    // An anchored sum keeps a constant sample exactly constant and avoids
    // losing a small spread to rounding of a large common offset.
    let mean_offset = values
        .iter()
        .map(|&v| scaled_difference(v, first, scale))
        .sum::<f64>()
        / n as f64;
    let mean = (anchor + mean_offset).clamp(-1.0, 1.0);
    let variance = if n > 1 {
        values
            .iter()
            .map(|&v| (scaled_difference(v, first, scale) - mean_offset).powi(2))
            .sum::<f64>()
            / n.saturating_sub(1) as f64
    } else {
        0.0
    };
    Ok(ScaledMoments {
        mean,
        variance,
        scale,
    })
}

/// Average 1-based ranks with tie averaging, plus `sum(t^3 - t)` for tie groups.
fn average_ranks_and_tie_sum(values: &[f64]) -> (Vec<f64>, f64) {
    let n = values.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&a, &b| values[a].total_cmp(&values[b]));

    let mut ranks = vec![0.0; n];
    let mut tie_sum = 0.0;
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && crate::kernels::exact_eq(values[order[j]], values[order[i]]) {
            j += 1;
        }
        let count = (j - i) as f64;
        if count > 1.0 {
            tie_sum += count * count * count - count;
        }
        let avg_rank = f64::midpoint((i + 1) as f64, j as f64);
        for k in i..j {
            ranks[order[k]] = avg_rank;
        }
        i = j;
    }
    (ranks, tie_sum)
}

fn two_sided_t_p_value(
    t_stat: f64,
    df: f64,
) -> f64 {
    let abs_t = t_stat.abs();
    let root_df = df.sqrt();
    let hypotenuse = root_df.hypot(abs_t);
    let ratio = root_df / hypotenuse;
    let a = df * 0.5;
    if abs_t < root_df * 0.5 {
        let complement = (abs_t / hypotenuse).powi(2);
        return (1.0 - regularized_beta(complement, 0.5, a)).clamp(0.0, 1.0);
    }
    if ratio < f64::MIN_POSITIVE.sqrt() {
        // In the far tail I_x(a, b) = x^a / (a B(a, b)) to floating-point
        // accuracy. Here x may underflow even though x^a is representable.
        let log_x = 2.0 * (root_df.ln() - abs_t.ln());
        return (a * log_x - a.ln() - ln_beta_numerical(a, 0.5)).exp();
    }
    regularized_beta(ratio * ratio, a, 0.5).clamp(0.0, 1.0)
}

/// One-sample Student's t-test against population mean `mu0`.
///
/// # Errors
/// [`StatsExtError::Domain`] for fewer than two observations, non-finite
/// observations or null mean, or a sample with zero variance.
pub fn t_test_one_sample(
    sample: &[f64],
    mu0: f64,
) -> Result<TestResult, StatsExtError> {
    if sample.len() < 2 || !mu0.is_finite() {
        return Err(StatsExtError::Domain(
            "t_test_one_sample requires at least 2 observations and a finite mean".into(),
        ));
    }
    let moments = scaled_moments(sample)?;
    if moments.variance == 0.0 {
        return Err(StatsExtError::Domain(
            "t_test_one_sample is undefined for zero-variance samples".into(),
        ));
    }
    let n = sample.len() as f64;
    let se = (moments.variance / n).sqrt();
    let statistic = (moments.mean - mu0 / moments.scale) / se;
    let df = n - 1.0;
    Ok(TestResult {
        statistic,
        p_value: two_sided_t_p_value(statistic, df),
        df: Some(df),
    })
}

/// Two-sample t-test with pooled variance or Welch's unequal-variance formula.
///
/// # Errors
/// [`StatsExtError::Domain`] for fewer than two observations per sample,
/// non-finite data, or zero combined standard error.
pub fn t_test_two_sample(
    left: &[f64],
    right: &[f64],
    equal_var: bool,
) -> Result<TestResult, StatsExtError> {
    if left.len() < 2 || right.len() < 2 {
        return Err(StatsExtError::Domain(
            "t_test_two_sample requires at least 2 observations per sample".into(),
        ));
    }
    let left_moments = scaled_moments(left)?;
    let right_moments = scaled_moments(right)?;
    let customary_scale = |sample: &ScaledMoments| {
        if sample.variance > 0.0 {
            sample.scale
        } else {
            0.0
        }
    };
    let scale = customary_scale(&left_moments).max(customary_scale(&right_moments));
    if scale == 0.0 {
        return Err(StatsExtError::Domain(
            "t_test_two_sample is undefined when both samples have zero variance".into(),
        ));
    }
    let std_at_scale =
        |sample: &ScaledMoments| sample.variance.sqrt() * (customary_scale(sample) / scale);
    let left_std = std_at_scale(&left_moments);
    let right_std = std_at_scale(&right_moments);
    let n1 = left.len() as f64;
    let n2 = right.len() as f64;
    let (se, df) = if equal_var {
        let df = n1 + n2 - 2.0;
        let pooled =
            (left_std * ((n1 - 1.0) / df).sqrt()).hypot(right_std * ((n2 - 1.0) / df).sqrt());
        (pooled * (1.0 / n1 + 1.0 / n2).sqrt(), df)
    } else {
        let s1 = left_std / n1.sqrt();
        let s2 = right_std / n2.sqrt();
        let se = s1.hypot(s2);
        let df = 1.0 / ((s1 / se).powi(4) / (n1 - 1.0) + (s2 / se).powi(4) / (n2 - 1.0));
        (se, df)
    };
    let left_mean = left_moments.mean * left_moments.scale;
    let right_mean = right_moments.mean * right_moments.scale;
    let statistic = scaled_difference(left_mean, right_mean, scale) / se;
    Ok(TestResult {
        statistic,
        p_value: two_sided_t_p_value(statistic, df),
        df: Some(df),
    })
}

/// Paired-samples t-test on `left - right`.
///
/// # Errors
///
/// [`StatsExtError::ShapeMismatch`] if lengths differ, or [`StatsExtError::Domain`]
/// if fewer than 2 pairs or differences have zero variance.
pub fn t_test_paired(
    left: &[f64],
    right: &[f64],
) -> Result<TestResult, StatsExtError> {
    if left.len() != right.len() {
        return Err(StatsExtError::ShapeMismatch(format!(
            "t_test_paired needs equal-length samples, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    let scale = max_abs(left).max(max_abs(right)).max(f64::MIN_POSITIVE);
    let diffs: Vec<f64> = left
        .iter()
        .zip(right)
        .map(|(&a, &b)| scaled_difference(a, b, scale))
        .collect();
    t_test_one_sample(&diffs, 0.0)
}

/// Pearson's chi-squared goodness-of-fit test.
///
/// # Errors
///
/// [`StatsExtError::ShapeMismatch`] if lengths differ, or [`StatsExtError::Domain`]
/// if fewer than 2 categories or any expected count is `<= 0`.
pub fn chi_squared_goodness_of_fit(
    observed: &[f64],
    expected: &[f64],
) -> Result<TestResult, StatsExtError> {
    finite(observed)?;
    finite(expected)?;
    if observed.len() != expected.len() {
        return Err(StatsExtError::ShapeMismatch(format!(
            "chi_squared_goodness_of_fit needs equal lengths, found {} and {}",
            observed.len(),
            expected.len()
        )));
    }
    if observed.len() < 2 {
        return Err(StatsExtError::Domain(
            "chi_squared_goodness_of_fit requires at least 2 categories".to_string(),
        ));
    }
    let mut stat = 0.0;
    for (&o, &e) in observed.iter().zip(expected) {
        if e <= 0.0 || o < 0.0 {
            return Err(StatsExtError::Domain(
                "observed counts must be >= 0 and expected counts must be > 0".to_string(),
            ));
        }
        let diff = o - e;
        stat += (diff / e.sqrt()).powi(2);
    }
    let df = (observed.len() - 1) as f64;
    let p_value =
        crate::kernels::special::regularized_upper_gamma(df * 0.5, stat * 0.5).clamp(0.0, 1.0);
    Ok(TestResult {
        statistic: stat,
        p_value,
        df: Some(df),
    })
}

/// Pearson's chi-squared test of independence on a row-major `rows × cols` contingency table.
///
/// # Errors
///
/// [`StatsExtError::ShapeMismatch`] if `observed.len() != rows * cols`, or
/// [`StatsExtError::Domain`] if `rows < 2`, `cols < 2`, any entry is negative,
/// or any row/column marginal total is zero.
pub fn chi_squared_independence(
    rows: usize,
    cols: usize,
    observed: &[f64],
) -> Result<TestResult, StatsExtError> {
    finite(observed)?;
    if rows.checked_mul(cols) != Some(observed.len()) {
        return Err(StatsExtError::ShapeMismatch(format!(
            "expected {rows}x{cols} contingency table, found {} entries",
            observed.len()
        )));
    }
    if rows < 2 || cols < 2 {
        return Err(StatsExtError::Domain(
            "chi_squared_independence requires at least a 2x2 table".to_string(),
        ));
    }
    let scale = max_abs(observed);
    if scale == 0.0 || observed.iter().any(|&value| value < 0.0) {
        return Err(StatsExtError::Domain(
            "contingency table needs nonnegative counts and positive marginals".into(),
        ));
    }
    let mut row_sums = vec![0.0; rows];
    let mut col_sums = vec![0.0; cols];
    for (row, row_sum) in observed.chunks_exact(cols).zip(&mut row_sums) {
        for (col_sum, &value) in col_sums.iter_mut().zip(row) {
            let relative = value / scale;
            *row_sum += relative;
            *col_sum += relative;
        }
    }
    if row_sums.contains(&0.0) || col_sums.contains(&0.0) {
        return Err(StatsExtError::Domain(
            "contingency table marginals must all be positive at the working scale".into(),
        ));
    }
    let total: f64 = row_sums.iter().sum();
    let mut scaled_statistic = 0.0;
    for (row, &row_sum) in observed.chunks_exact(cols).zip(&row_sums) {
        for (&value, &col_sum) in row.iter().zip(&col_sums) {
            let expected = (row_sum / total) * col_sum;
            let relative = value / scale;
            let standardized = if expected > 0.0 {
                (relative - expected) / expected.sqrt()
            } else {
                // Taking roots first retains expected counts that are too
                // small to represent in the scaled table itself.
                let root_expected = (row_sum.sqrt() / total.sqrt()) * col_sum.sqrt();
                relative / root_expected - root_expected
            };
            scaled_statistic += standardized.powi(2);
        }
    }
    let stat = scale * scaled_statistic;
    let df = ((rows - 1) * (cols - 1)) as f64;
    let p_value =
        crate::kernels::special::regularized_upper_gamma(df * 0.5, stat * 0.5).clamp(0.0, 1.0);
    Ok(TestResult {
        statistic: stat,
        p_value,
        df: Some(df),
    })
}

/// One-way Analysis of Variance (ANOVA) across `k >= 2` groups.
///
/// Returns `TestResult` with `statistic = F`, `df = Some(df_between)`, and upper-tail `p_value`.
///
/// # Errors
///
/// [`StatsExtError::Domain`] if fewer than 2 groups, any group is empty,
/// total observations `<= groups.len()`, or within-group variance is zero.
pub fn anova_one_way(groups: &[&[f64]]) -> Result<TestResult, StatsExtError> {
    let k = groups.len();
    if k < 2 {
        return Err(StatsExtError::Domain(
            "anova_one_way requires at least 2 groups".to_string(),
        ));
    }
    let mut n_total = 0usize;
    let mut samples = Vec::with_capacity(k);
    for group in groups {
        if group.is_empty() {
            return Err(StatsExtError::Domain(
                "anova_one_way groups must be non-empty".into(),
            ));
        }
        n_total = n_total.checked_add(group.len()).ok_or_else(|| {
            StatsExtError::Domain("anova_one_way observation count overflows".into())
        })?;
        samples.push((group.len(), scaled_moments(group)?));
    }
    if n_total <= k {
        return Err(StatsExtError::Domain(
            "anova_one_way requires more total observations than groups".into(),
        ));
    }
    let scale = samples
        .iter()
        .map(|(_, sample)| sample.scale)
        .fold(0.0, f64::max);
    let scaled_mean = |sample: &ScaledMoments| sample.mean * (sample.scale / scale);
    let grand_mean: f64 = samples
        .iter()
        .map(|(n, sample)| scaled_mean(sample) * (*n as f64 / n_total as f64))
        .sum();
    let mut ss_between = 0.0;
    let mut ss_within = 0.0;
    let mut variable = false;
    for (n, sample) in &samples {
        ss_between += *n as f64 * (scaled_mean(sample) - grand_mean).powi(2);
        let standard_deviation = sample.variance.sqrt() * (sample.scale / scale);
        ss_within += n.saturating_sub(1) as f64 * standard_deviation.powi(2);
        variable |= sample.variance > 0.0;
    }
    if !variable {
        return Err(StatsExtError::Domain(
            "anova_one_way is undefined when within-group variance is zero".into(),
        ));
    }
    let df1 = k.saturating_sub(1) as f64;
    let df2 = n_total.saturating_sub(k) as f64;
    let f_stat = if ss_within == 0.0 {
        f64::INFINITY
    } else {
        (ss_between / df1) / (ss_within / df2)
    };
    let ratio = df2 / df1;
    let p_value = regularized_beta(ratio / (ratio + f_stat), df2 * 0.5, df1 * 0.5).clamp(0.0, 1.0);
    Ok(TestResult {
        statistic: f_stat,
        p_value,
        df: Some(df1),
    })
}

fn ks_asymptotic_p_value(
    d: f64,
    n_eff: f64,
) -> f64 {
    if d <= 0.0 {
        return 1.0;
    }
    let en = n_eff.sqrt();
    let lambda = (en + 0.12 + 0.11 / en) * d;
    // The alternating survival series converges slowly as lambda -> 0.
    // Use its Jacobi-theta transform for the lower tail instead.
    if lambda < 1.18 {
        let mut lower_sum = 0.0;
        for k in 1..=100 {
            let odd = f64::from(2 * k - 1);
            let term = (-PI * PI * odd * odd / (8.0 * lambda * lambda)).exp();
            lower_sum += term;
            if term < 1e-16 {
                break;
            }
        }
        return (1.0 - (2.0 * PI).sqrt() * lower_sum / lambda).clamp(0.0, 1.0);
    }
    let mut sum = 0.0;
    for k in 1..=100 {
        let kf = f64::from(k);
        let term = (-2.0 * kf * kf * lambda * lambda).exp();
        if k % 2 == 1 {
            sum += term;
        } else {
            sum -= term;
        }
        if term < 1e-16 {
            break;
        }
    }
    (2.0 * sum).clamp(0.0, 1.0)
}

/// One-sample Kolmogorov–Smirnov test against a reference continuous CDF.
///
/// # Errors
///
/// [`StatsExtError::Empty`] when `sample` is empty.
pub fn ks_test_one_sample(
    sample: &[f64],
    cdf: &dyn Fn(f64) -> f64,
) -> Result<TestResult, StatsExtError> {
    finite(sample)?;
    if sample.is_empty() {
        return Err(StatsExtError::Empty);
    }
    let mut sorted = sample.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len() as f64;
    let mut d_max = 0.0f64;
    for (i, &x) in sorted.iter().enumerate() {
        let f_x = cdf(x);
        if !(0.0..=1.0).contains(&f_x) {
            return Err(StatsExtError::Domain(
                "reference CDF must return a value in [0, 1]".to_string(),
            ));
        }
        let d_plus = ((i + 1) as f64 / n - f_x).abs();
        let d_minus = (f_x - i as f64 / n).abs();
        d_max = d_max.max(d_plus).max(d_minus);
    }
    Ok(TestResult {
        statistic: d_max,
        p_value: ks_asymptotic_p_value(d_max, n),
        df: None,
    })
}

/// Two-sample Kolmogorov–Smirnov test comparing two empirical distributions.
///
/// # Errors
///
/// [`StatsExtError::Empty`] when either sample is empty.
pub fn ks_test_two_sample(
    left: &[f64],
    right: &[f64],
) -> Result<TestResult, StatsExtError> {
    finite(left)?;
    finite(right)?;
    if left.is_empty() || right.is_empty() {
        return Err(StatsExtError::Empty);
    }
    let mut a = left.to_vec();
    let mut b = right.to_vec();
    a.sort_by(f64::total_cmp);
    b.sort_by(f64::total_cmp);
    let n1 = a.len();
    let n2 = b.len();
    let mut i = 0usize;
    let mut j = 0usize;
    let mut d_max = 0.0f64;
    while i < n1 && j < n2 {
        let val = if a[i] <= b[j] { a[i] } else { b[j] };
        while i < n1 && crate::kernels::exact_eq(a[i], val) {
            i += 1;
        }
        while j < n2 && crate::kernels::exact_eq(b[j], val) {
            j += 1;
        }
        let cdf1 = i as f64 / n1 as f64;
        let cdf2 = j as f64 / n2 as f64;
        d_max = d_max.max((cdf1 - cdf2).abs());
    }
    let n_eff = (n1 as f64 * n2 as f64) / (n1 as f64 + n2 as f64);
    Ok(TestResult {
        statistic: d_max,
        p_value: ks_asymptotic_p_value(d_max, n_eff),
        df: None,
    })
}

/// Two-sided Mann–Whitney U test (Wilcoxon rank-sum test) with tie-corrected
/// normal approximation.
///
/// Returns `statistic = min(U1, U2)` and the two-sided asymptotic `p_value`.
///
/// # Errors
///
/// [`StatsExtError::Empty`] when either sample is empty, or [`StatsExtError::Domain`]
/// when all combined values are identical.
pub fn mann_whitney_u(
    left: &[f64],
    right: &[f64],
) -> Result<TestResult, StatsExtError> {
    finite(left)?;
    finite(right)?;
    if left.is_empty() || right.is_empty() {
        return Err(StatsExtError::Empty);
    }
    let n1 = left.len() as f64;
    let n2 = right.len() as f64;
    let n = n1 + n2;
    let mut combined = Vec::with_capacity(left.len() + right.len());
    combined.extend_from_slice(left);
    combined.extend_from_slice(right);

    let (ranks, tie_sum) = average_ranks_and_tie_sum(&combined);
    let r1: f64 = ranks[..left.len()].iter().sum();
    let u1 = r1 - n1 * (n1 + 1.0) / 2.0;
    let u2 = n1 * n2 - u1;
    let u_stat = u1.min(u2);

    let mean_u = n1 * n2 / 2.0;
    let var_u = (n1 * n2 / 12.0) * ((n + 1.0) - tie_sum / (n * (n - 1.0)));
    if var_u <= 0.0 {
        return Err(StatsExtError::Domain(
            "mann_whitney_u is undefined when all observations are tied".to_string(),
        ));
    }
    let z = (u1 - mean_u).abs() / var_u.sqrt();
    let p_value = (2.0 * Distribution::Normal { mu: 0.0, sigma: 1.0 }.cdf(-z)).clamp(0.0, 1.0);
    Ok(TestResult {
        statistic: u_stat,
        p_value,
        df: None,
    })
}

// ----- regression -----------------------------------------------------------

fn validate_regression_inputs(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
) -> Result<(), StatsExtError> {
    finite(x)?;
    finite(y)?;
    if rows == 0 || cols == 0 {
        return Err(StatsExtError::Empty);
    }
    if rows.checked_mul(cols) != Some(x.len()) {
        return Err(StatsExtError::ShapeMismatch(format!(
            "expected {rows}x{cols} design matrix, found {} entries",
            x.len()
        )));
    }
    if y.len() != rows {
        return Err(StatsExtError::ShapeMismatch(format!(
            "expected response vector of length {rows}, found {}",
            y.len()
        )));
    }
    Ok(())
}

fn center_design(
    rows: usize,
    cols: usize,
    x: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let mut x_mean = vec![0.0; cols];
    for row in x.chunks_exact(cols) {
        for (mean, value) in x_mean.iter_mut().zip(row) {
            *mean += *value;
        }
    }
    for (column, mean) in x_mean.iter_mut().enumerate() {
        *mean /= rows as f64;
        if !mean.is_finite() {
            // Fall back to a scaled sum when a finite column mean overflowed
            // on the ordinary summation path.
            let values = || x.iter().skip(column).step_by(cols);
            let scale = values().map(|value| value.abs()).fold(0.0, f64::max);
            let normalized =
                (values().map(|value| value / scale).sum::<f64>() / rows as f64).clamp(-1.0, 1.0);
            *mean = normalized * scale;
        }
    }
    let centered = x
        .chunks_exact(cols)
        .flat_map(|row| row.iter().zip(&x_mean).map(|(&value, &mean)| value - mean))
        .collect();
    (centered, x_mean)
}

type CenteredData<'a> = (Cow<'a, [f64]>, Cow<'a, [f64]>, Vec<f64>, f64);

fn center_xy<'a>(
    rows: usize,
    cols: usize,
    x: &'a [f64],
    y: &'a [f64],
    fit_intercept: bool,
) -> Result<CenteredData<'a>, StatsExtError> {
    if !fit_intercept {
        return Ok((Cow::Borrowed(x), Cow::Borrowed(y), Vec::new(), 0.0));
    }
    let (xc, x_mean) = center_design(rows, cols, x);
    let moments = scaled_moments(y)?;
    let y_mean = moments.mean * moments.scale;
    let yc: Vec<f64> = y.iter().map(|&v| v - y_mean).collect();
    finite(&xc)?;
    finite(&yc)?;
    Ok((Cow::Owned(xc), Cow::Owned(yc), x_mean, y_mean))
}

fn build_linear_result(
    cols: usize,
    x: &[f64],
    y: &[f64],
    coefficients: Vec<f64>,
    x_mean: &[f64],
    y_mean: f64,
    fit_intercept: bool,
) -> Result<LinearRegressionResult, StatsExtError> {
    finite(&coefficients)?;
    let intercept = if fit_intercept {
        y_mean
            - x_mean
                .iter()
                .zip(&coefficients)
                .map(|(m, b)| m * b)
                .sum::<f64>()
    } else {
        0.0
    };
    let residuals: Vec<_> = x
        .chunks_exact(cols)
        .zip(y)
        .map(|(row, &target)| {
            target
                - row
                    .iter()
                    .zip(&coefficients)
                    .fold(intercept, |pred, (value, weight)| pred + value * weight)
        })
        .collect();
    finite(&residuals)?;

    let response = scaled_moments(y)?;
    let sum_squares = response.variance * y.len().saturating_sub(1) as f64;
    let residual_scale = max_abs(&residuals);
    let r_squared = if sum_squares == 0.0 {
        // A constant response is explained only by an exact constant fit.
        if residual_scale == 0.0 {
            1.0
        } else {
            0.0
        }
    } else if residual_scale == 0.0 {
        1.0
    } else {
        let residual_norm = residuals
            .iter()
            .map(|&r| (r / residual_scale).powi(2))
            .sum::<f64>()
            .sqrt();
        let relative = (residual_scale / response.scale) * (residual_norm / sum_squares.sqrt());
        1.0 - relative.powi(2)
    };
    Ok(LinearRegressionResult {
        coefficients,
        intercept,
        r_squared,
        residuals,
    })
}

/// Ordinary least squares with a minimum-norm solution for deficient designs.
///
/// # Errors
///
/// Rejects invalid dimensions, non-finite data or results, or failed decomposition.
pub fn ols(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    fit_intercept: bool,
) -> Result<LinearRegressionResult, StatsExtError> {
    validate_regression_inputs(rows, cols, x, y)?;
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept)?;
    let design = Mat {
        rows,
        cols,
        data: xc.into_owned(),
    };
    let coefficients = dense::lstsq_svd(&design, &yc).map_err(map_fit_error)?;
    build_linear_result(cols, x, y, coefficients, &x_mean, y_mean, fit_intercept)
}

/// Ridge (L2-regularized) linear regression minimizing `||y - Xβ||₂² + alpha ||β||₂²`.
///
/// Uses the thin SVD of the design without allocating a square feature matrix.
///
/// # Errors
///
/// Rejects invalid dimensions, non-finite data or results, and failed decomposition.
/// A negative or non-finite `alpha` is invalid. For `alpha == 0`, a design without
/// full column rank returns [`StatsExtError::Singular`].
pub fn ridge(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    alpha: f64,
    fit_intercept: bool,
) -> Result<LinearRegressionResult, StatsExtError> {
    validate_regression_inputs(rows, cols, x, y)?;
    if alpha < 0.0 || !alpha.is_finite() {
        return Err(StatsExtError::Domain(
            "ridge regularization parameter alpha must be >= 0".to_string(),
        ));
    }
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept)?;
    let design = Mat {
        rows,
        cols,
        data: xc.into_owned(),
    };
    let factors = dense::svd(&design).map_err(map_fit_error)?;
    let k = factors.s.len();
    if alpha == 0.0 {
        let largest = factors.s.first().copied().unwrap_or(0.0);
        let cutoff = (rows.max(cols) as f64 * f64::EPSILON) * largest;
        if k < cols || factors.s.last().is_none_or(|&value| value <= cutoff) {
            return Err(StatsExtError::Singular);
        }
    }

    let root_alpha = alpha.sqrt();
    let mut projected = vec![0.0; k];
    for (row, &target) in factors.u.data.chunks_exact(k).zip(yc.iter()) {
        for (coordinate, value) in projected.iter_mut().zip(row) {
            *coordinate += target * value;
        }
    }
    for (coordinate, &sigma) in projected.iter_mut().zip(&factors.s) {
        // Evaluate sigma / (sigma² + alpha) without squaring sigma or alpha.
        if sigma > root_alpha {
            *coordinate = (*coordinate / sigma) / (1.0 + (root_alpha / sigma).powi(2));
        } else if alpha > 0.0 {
            *coordinate *= (sigma / alpha) / (1.0 + (sigma / root_alpha).powi(2));
        }
    }
    let coefficients = factors
        .v
        .data
        .chunks_exact(k)
        .map(|row| row.iter().zip(&projected).map(|(v, w)| v * w).sum())
        .collect();
    build_linear_result(cols, x, y, coefficients, &x_mean, y_mean, fit_intercept)
}

fn soft_threshold(
    rho: f64,
    lambda: f64,
) -> f64 {
    if rho > lambda {
        rho - lambda
    } else if rho < -lambda {
        rho + lambda
    } else {
        0.0
    }
}

/// Lasso (L1-regularized) linear regression via coordinate descent,
/// minimizing `(1 / (2n)) ||y - Xβ||₂² + alpha ||β||₁`.
///
/// # Errors
///
/// [`StatsExtError::Empty`], [`StatsExtError::ShapeMismatch`], [`StatsExtError::Domain`]
/// when `alpha < 0` or `tol <= 0`, or [`StatsExtError::NotConverged`].
pub fn lasso(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    alpha: f64,
    fit_intercept: bool,
    max_iter: usize,
    tol: f64,
) -> Result<LinearRegressionResult, StatsExtError> {
    validate_regression_inputs(rows, cols, x, y)?;
    if alpha < 0.0 || !alpha.is_finite() {
        return Err(StatsExtError::Domain(
            "lasso regularization parameter alpha must be >= 0".to_string(),
        ));
    }
    if !tol.is_finite() || tol <= 0.0 || max_iter == 0 {
        return Err(StatsExtError::Domain(
            "lasso requires tol > 0 and max_iter > 0".to_string(),
        ));
    }
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept)?;
    let n = rows as f64;
    let mut col_sq_norms = vec![0.0; cols];
    for r in 0..rows {
        for c in 0..cols {
            let v = xc[r * cols + c];
            col_sq_norms[c] += v * v;
        }
    }

    finite(&col_sq_norms)?;
    let mut beta = vec![0.0; cols];
    let mut residual = yc.into_owned();
    let mut converged = false;

    for _ in 0..max_iter {
        let mut max_delta = 0.0f64;
        for j in 0..cols {
            let norm_sq = col_sq_norms[j];
            if norm_sq == 0.0 {
                continue;
            }
            let old_b = beta[j];
            let mut rho = 0.0;
            for r in 0..rows {
                let xij = xc[r * cols + j];
                rho += xij * (residual[r] + old_b * xij);
            }
            let new_b = soft_threshold(rho / n, alpha) / (norm_sq / n);
            if !new_b.is_finite() {
                return Err(StatsExtError::NotConverged);
            }
            let diff = new_b - old_b;
            if diff != 0.0 {
                for r in 0..rows {
                    residual[r] -= diff * xc[r * cols + j];
                }
                beta[j] = new_b;
                max_delta = max_delta.max(diff.abs());
            }
        }
        if max_delta < tol {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(StatsExtError::NotConverged);
    }
    build_linear_result(cols, x, y, beta, &x_mean, y_mean, fit_intercept)
}

/// Binary logistic regression via Newton-Raphson (IRLS).
///
/// Expects `y[i]` in `[0, 1]` (typically `0.0` or `1.0`).
///
/// # Errors
///
/// [`StatsExtError::Empty`], [`StatsExtError::ShapeMismatch`], [`StatsExtError::Domain`]
/// for invalid `y` or tolerances, or [`StatsExtError::NotConverged`].
pub fn logistic_regression(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    fit_intercept: bool,
    max_iter: usize,
    tol: f64,
) -> Result<LogisticRegressionResult, StatsExtError> {
    validate_regression_inputs(rows, cols, x, y)?;
    if y.iter().any(|&v| !(0.0..=1.0).contains(&v)) {
        return Err(StatsExtError::Domain(
            "logistic_regression targets must be in [0, 1]".to_string(),
        ));
    }
    if !tol.is_finite() || tol <= 0.0 || max_iter == 0 {
        return Err(StatsExtError::Domain(
            "logistic_regression requires tol > 0 and max_iter > 0".to_string(),
        ));
    }

    let width_error =
        || StatsExtError::Domain("logistic design is too wide for its Hessian".into());
    let p = cols
        .checked_add(usize::from(fit_intercept))
        .ok_or_else(width_error)?;
    let hessian_len = p.checked_mul(p).ok_or_else(width_error)?;
    let mut data = Vec::new();
    data.try_reserve_exact(hessian_len)
        .map_err(|_| width_error())?;
    data.resize(hessian_len, 0.0);
    let mut hessian = Mat {
        data,
        rows: p,
        cols: p,
    };
    let mut grad = vec![0.0; p];
    let mut coefficients = vec![0.0; p];
    let mut features = vec![1.0; p];

    for iter in 1..=max_iter {
        grad.fill(0.0);
        hessian.data.fill(0.0);
        for (row, &target) in x.chunks_exact(cols).zip(y) {
            let values = if fit_intercept {
                features[1..].copy_from_slice(row);
                features.as_slice()
            } else {
                row
            };
            let eta = values
                .iter()
                .zip(&coefficients)
                .map(|(v, w)| v * w)
                .sum::<f64>();
            if !eta.is_finite() {
                return Err(StatsExtError::NotConverged);
            }
            let z = (-eta.abs()).exp();
            let small_probability = z / (1.0 + z);
            let weight = (z / (1.0 + z).powi(2)).max(1e-12);
            let error = if eta >= 0.0 {
                (1.0 - target) - small_probability
            } else {
                small_probability - target
            };
            for (j, (grow, (hrow, xj))) in grad
                .iter_mut()
                .zip(hessian.data.chunks_exact_mut(p).zip(values))
                .enumerate()
            {
                *grow += xj * error;
                for (entry, xk) in hrow.iter_mut().zip(values).take(j + 1) {
                    *entry += (weight * xj) * xk;
                }
            }
        }
        // The Hessian is symmetric. The lower triangle is enough work for each row.
        for j in 0..p {
            for k in 0..j {
                hessian.data[k * p + j] = hessian.data[j * p + k];
            }
            hessian.data[j * p + j] += 1e-8;
        }
        let step = dense::solve(&hessian, &grad).map_err(map_fit_error)?;
        let mut max_step = 0.0_f64;
        for (weight, update) in coefficients.iter_mut().zip(step) {
            *weight -= update;
            max_step = max_step.max(update.abs());
        }
        finite(&coefficients)?;
        if max_step < tol {
            return logistic_result(cols, x, y, coefficients, fit_intercept, iter);
        }
    }
    Err(StatsExtError::NotConverged)
}

fn logistic_result(
    cols: usize,
    x: &[f64],
    y: &[f64],
    mut coefficients: Vec<f64>,
    fit_intercept: bool,
    iterations: usize,
) -> Result<LogisticRegressionResult, StatsExtError> {
    let intercept = if fit_intercept {
        coefficients.remove(0)
    } else {
        0.0
    };
    let mut log_likelihood = 0.0;
    for (row, &target) in x.chunks_exact(cols).zip(y) {
        let eta = row
            .iter()
            .zip(&coefficients)
            .fold(intercept, |sum, (value, weight)| sum + value * weight);
        if !eta.is_finite() {
            return Err(StatsExtError::NotConverged);
        }
        let loss = if eta >= 0.0 {
            (1.0 - target) * eta
        } else {
            -target * eta
        };
        log_likelihood -= loss + (-eta.abs()).exp().ln_1p();
    }
    Ok(LogisticRegressionResult {
        coefficients,
        intercept,
        log_likelihood,
        iterations,
    })
}

// ----- PCA ------------------------------------------------------------------

/// Principal Component Analysis via SVD of the mean-centered data matrix.
///
/// `data` is row-major `rows × cols`. If `n_components` is `None`, retains
/// `min(rows, cols)` components.
///
/// # Errors
///
/// [`StatsExtError::Empty`], [`StatsExtError::ShapeMismatch`], [`StatsExtError::Domain`]
/// when `n_components` is 0 or exceeds `min(rows, cols)`, or [`StatsExtError::NotConverged`].
pub fn pca(
    rows: usize,
    cols: usize,
    data: &[f64],
    n_components: Option<usize>,
) -> Result<PcaResult, StatsExtError> {
    finite(data)?;
    if rows == 0 || cols == 0 {
        return Err(StatsExtError::Empty);
    }
    if rows.checked_mul(cols) != Some(data.len()) {
        return Err(StatsExtError::ShapeMismatch(format!(
            "expected {rows}x{cols} data matrix, found {} entries",
            data.len()
        )));
    }
    let max_comp = rows.min(cols);
    let k = n_components.unwrap_or(max_comp);
    if k == 0 || k > max_comp {
        return Err(StatsExtError::Domain(format!(
            "n_components must be in 1..={max_comp}, got {k}"
        )));
    }

    let (data, mean) = center_design(rows, cols, data);
    let centered = Mat { rows, cols, data };
    let factors = dense::svd(&centered).map_err(map_fit_error)?;
    let (u, s, vt) = (factors.u.data, factors.s, factors.v.transpose().data);

    let standard_deviation_divisor = (rows.saturating_sub(1).max(1) as f64).sqrt();
    let explained_variance = s
        .iter()
        .take(k)
        .map(|value| (value / standard_deviation_divisor).powi(2))
        .collect();
    let singular_scale = s.first().copied().unwrap_or(0.0);
    let explained_variance_ratio = if singular_scale > 0.0 {
        let scaled: Vec<_> = s
            .iter()
            .map(|value| (value / singular_scale).powi(2))
            .collect();
        let total: f64 = scaled.iter().sum();
        scaled
            .into_iter()
            .take(k)
            .map(|value| value / total)
            .collect()
    } else {
        vec![0.0; k]
    };
    let components = vt[..k * cols].to_vec();

    // U is rows × max_comp; scores = U[:, :k] * diag(s[:k]).
    let mut transformed = vec![0.0; rows * k];
    for r in 0..rows {
        for c in 0..k {
            transformed[r * k + c] = u[r * max_comp + c] * s[c];
        }
    }

    Ok(PcaResult {
        mean,
        components,
        explained_variance,
        explained_variance_ratio,
        transformed,
        n_components: k,
    })
}

// ----- rank correlations ----------------------------------------------------

/// Spearman's rank correlation coefficient (with average ranks for ties).
///
/// # Errors
///
/// [`StatsExtError::ShapeMismatch`] when lengths differ, or [`StatsExtError::Domain`]
/// when either sample has zero variance.
pub fn spearman(
    left: &[f64],
    right: &[f64],
) -> Result<f64, StatsExtError> {
    finite(left)?;
    finite(right)?;
    if left.len() != right.len() {
        return Err(StatsExtError::ShapeMismatch(format!(
            "spearman needs two samples of one length, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    if left.is_empty() {
        return Err(StatsExtError::Empty);
    }
    let rank = |sample: &&[f64]| average_ranks_and_tie_sum(sample).0;
    let ranks = par::map_slice(
        &[left, right],
        if left.len() >= 16_384 { 2 } else { 3 },
        rank,
    );
    descriptive::corr(&ranks[0], &ranks[1], false)
        .map_err(|error| StatsExtError::Domain(error.to_string()))
}

/// Kendall's tau-b rank correlation, adjusted for ties in either sample.
///
/// Lexicographic sorting and merge-based inversion counting take
/// `O(n log n)` time and `O(n)` scratch space.
///
/// # Errors
/// Returns an error for mismatched lengths, fewer than two observations,
/// non-finite observations, or a constant sample.
pub fn kendall(
    left: &[f64],
    right: &[f64],
) -> Result<f64, StatsExtError> {
    finite(left)?;
    finite(right)?;
    if left.len() != right.len() {
        return Err(StatsExtError::ShapeMismatch(format!(
            "kendall needs two samples of one length, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    if left.is_empty() {
        return Err(StatsExtError::Empty);
    }
    if left.len() < 2 {
        return Err(StatsExtError::Domain(
            "kendall requires at least 2 observations".into(),
        ));
    }

    let mut observations: Vec<(f64, f64)> =
        left.iter().copied().zip(right.iter().copied()).collect();
    observations.sort_unstable_by(|(ax, ay), (bx, by)| {
        ax.partial_cmp(bx)
            .unwrap_or(Ordering::Equal)
            .then_with(|| ay.partial_cmp(by).unwrap_or(Ordering::Equal))
    });
    let mut left_ties = 0;
    let mut joint_ties = 0;
    for left_group in observations.chunk_by(|a, b| crate::kernels::exact_eq(a.0, b.0)) {
        left_ties += pair_count(left_group.len());
        for joint_group in left_group.chunk_by(|a, b| crate::kernels::exact_eq(a.1, b.1)) {
            joint_ties += pair_count(joint_group.len());
        }
    }
    let pairs = pair_count(observations.len());
    if left_ties == pairs {
        return Err(StatsExtError::Domain(
            "kendall is undefined for a constant sample".into(),
        ));
    }

    // Within each left tie, the right values already ascend. Strict inversions
    // of the right values therefore count only genuinely discordant pairs.
    let mut sorted: Vec<f64> = observations.into_iter().map(|(_, y)| y).collect();
    let mut scratch = vec![0.0; sorted.len()];
    let discordant = sort_and_count_inversions(&mut sorted, &mut scratch);
    let right_ties: u128 = sorted
        .chunk_by(|&a, &b| crate::kernels::exact_eq(a, b))
        .map(|group| pair_count(group.len()))
        .sum();
    if right_ties == pairs {
        return Err(StatsExtError::Domain(
            "kendall is undefined for a constant sample".into(),
        ));
    }

    let comparable = pairs - left_ties - (right_ties - joint_ties);
    let concordant = comparable - discordant;
    let denom = ((pairs - left_ties) as f64).sqrt() * ((pairs - right_ties) as f64).sqrt();
    Ok((((concordant as i128 - discordant as i128) as f64) / denom).clamp(-1.0, 1.0))
}

const fn pair_count(count: usize) -> u128 {
    count as u128 * count.saturating_sub(1) as u128 / 2
}

fn sort_and_count_inversions(
    values: &mut [f64],
    scratch: &mut [f64],
) -> u128 {
    if values.len() < 2 {
        return 0;
    }
    let (left, right) = values.split_at_mut(values.len() / 2);
    let (left_scratch, right_scratch) = scratch.split_at_mut(left.len());
    let mut inversions = sort_and_count_inversions(left, left_scratch)
        + sort_and_count_inversions(right, right_scratch);
    if left
        .last()
        .zip(right.first())
        .is_some_and(|(lower, upper)| lower <= upper)
    {
        return inversions;
    }
    let mut left = &*left;
    let mut right = &*right;
    for slot in &mut *scratch {
        match (left.split_first(), right.split_first()) {
            | (Some((&a, rest)), Some((&b, _))) if a <= b => {
                *slot = a;
                left = rest;
            },
            | (_, Some((&b, rest))) => {
                *slot = b;
                right = rest;
                inversions += left.len() as u128;
            },
            | (Some((&a, rest)), None) => {
                *slot = a;
                left = rest;
            },
            | (None, None) => break,
        }
    }
    values.copy_from_slice(scratch);
    inversions
}

// ----- KDE ------------------------------------------------------------------

/// Gaussian Kernel Density Estimation evaluated at `points`.
///
/// When `bandwidth` is `None`, uses Silverman's rule of thumb.
///
/// # Errors
///
/// [`StatsExtError::Empty`] if `samples` is empty, or [`StatsExtError::Domain`]
/// if `bandwidth <= 0` or cannot be determined from a zero-spread sample.
pub fn kde(
    samples: &[f64],
    points: &[f64],
    bandwidth: Option<f64>,
) -> Result<Vec<f64>, StatsExtError> {
    finite(samples)?;
    finite(points)?;
    if samples.is_empty() {
        return Err(StatsExtError::Empty);
    }
    let h = match bandwidth {
        | Some(bw) => {
            if bw <= 0.0 || !bw.is_finite() {
                return Err(StatsExtError::Domain(
                    "kde bandwidth must be positive and finite".to_string(),
                ));
            }
            bw
        },
        | None => {
            let bw = robust::bandwidth_silverman(samples);
            if bw <= 0.0 || !bw.is_finite() {
                return Err(StatsExtError::Domain(
                    "cannot infer kde bandwidth for a zero-spread sample".to_string(),
                ));
            }
            bw
        },
    };

    Ok(robust::gaussian_kde_finite(samples, points, h))
}

// ----- resampling (bootstrap & jackknife) -----------------------------------

fn finite_statistic(
    statistic: &dyn Fn(&[f64]) -> Result<f64, StatsExtError>,
    sample: &[f64],
) -> Result<f64, StatsExtError> {
    let value = statistic(sample)?;
    if !value.is_finite() {
        return Err(StatsExtError::Domain(
            "resampling statistic must return a finite value".into(),
        ));
    }
    Ok(value)
}

/// Non-parametric bootstrap resampling with `n_resamples` draws.
///
/// # Errors
///
/// [`StatsExtError::Empty`] when `samples` is empty, [`StatsExtError::Domain`]
/// when the resample count is zero or cannot be allocated, a uniform draw is
/// outside `[0, 1)`, or a statistic is non-finite. Callback errors propagate.
pub fn bootstrap(
    samples: &[f64],
    n_resamples: usize,
    statistic: &dyn Fn(&[f64]) -> Result<f64, StatsExtError>,
    uniform: &mut dyn FnMut() -> f64,
) -> Result<BootstrapResult, StatsExtError> {
    if samples.is_empty() {
        return Err(StatsExtError::Empty);
    }
    if n_resamples == 0 {
        return Err(StatsExtError::Domain(
            "bootstrap requires n_resamples > 0".to_string(),
        ));
    }
    let mut replicates = Vec::new();
    replicates.try_reserve_exact(n_resamples).map_err(|_| {
        StatsExtError::Domain("bootstrap resample count exceeds available storage".into())
    })?;
    let estimate = finite_statistic(statistic, samples)?;
    let n = samples.len();
    let mut resample = f64_scratch(n.max(n_resamples));
    resample.resize(n, 0.0);

    for _ in 0..n_resamples {
        for slot in resample.iter_mut() {
            let draw = uniform();
            if !(0.0..1.0).contains(&draw) {
                return Err(StatsExtError::Domain(
                    "uniform source must return values in [0, 1)".to_string(),
                ));
            }
            let idx = (draw * n as f64).floor() as usize;
            *slot = samples[idx.min(n.saturating_sub(1))];
        }
        replicates.push(finite_statistic(statistic, &resample)?);
    }

    let moments = scaled_moments(&replicates)?;
    resample.clear();
    resample.extend_from_slice(&replicates);
    let ci_low = robust::quantile_select(&mut resample, 0.025, 7);
    let ci_high = robust::quantile_select(&mut resample, 0.975, 7);

    Ok(BootstrapResult {
        estimate,
        mean: moments.mean * moments.scale,
        std_error: moments.variance.sqrt() * moments.scale,
        ci_low,
        ci_high,
        replicates,
    })
}

/// Leave-one-out jackknife resampling for bias and standard error estimation.
///
/// # Errors
///
/// [`StatsExtError::Domain`] for fewer than two observations or a non-finite
/// statistic. Callback errors propagate.
pub fn jackknife(
    samples: &[f64],
    statistic: &dyn Fn(&[f64]) -> Result<f64, StatsExtError>,
) -> Result<JackknifeResult, StatsExtError> {
    let n = samples.len();
    if n < 2 {
        return Err(StatsExtError::Domain(
            "jackknife requires at least 2 observations".to_string(),
        ));
    }
    let estimate = finite_statistic(statistic, samples)?;
    let mut jackknife_estimates = Vec::with_capacity(n);
    let mut loo = f64_scratch(n - 1);
    loo.extend_from_slice(&samples[1..]);

    for (i, &sample) in samples.iter().enumerate() {
        jackknife_estimates.push(finite_statistic(statistic, &loo)?);
        // Moving the omitted index from i to i+1 changes just one element.
        if let Some(slot) = loo.get_mut(i) {
            *slot = sample;
        }
    }

    let moments = scaled_moments(&jackknife_estimates)?;
    let correction = n.saturating_sub(1) as f64;
    let bias = correction * (moments.mean * moments.scale - estimate);
    let std_error = moments.scale * (correction * (moments.variance / n as f64).sqrt());

    Ok(JackknifeResult {
        estimate,
        bias,
        std_error,
        jackknife_estimates,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    fn close(
        actual: f64,
        expected: f64,
    ) {
        assert!(
            (actual - expected).abs() <= expected.abs().max(1.0) * 5e-12,
            "{actual} != {expected}",
        );
    }

    #[test]
    fn t_statistics_are_scale_invariant_and_preserve_small_tails() -> Result<(), StatsExtError> {
        let left = [1.0, 2.0, 4.0];
        let right = [3.0, 5.0, 7.0];
        let big_left = left.map(|value| value * 1e250);
        let big_right = right.map(|value| value * 1e250);
        let small = t_test_one_sample(&left, 0.0)?;
        let big = t_test_one_sample(&big_left, 0.0)?;
        close(big.statistic, small.statistic);
        close(big.p_value, small.p_value);
        for equal_variance in [false, true] {
            let reference = t_test_two_sample(&left, &right, equal_variance)?;
            let result = t_test_two_sample(&big_left, &big_right, equal_variance)?;
            close(result.statistic, reference.statistic);
            close(result.p_value, reference.p_value);
        }
        assert!(t_test_one_sample(&left, f64::INFINITY).is_err());

        let tail = t_test_one_sample(&[-1.0, 1.0], -1e200)?;
        let exact_cauchy_tail = (2.0 / PI) * (1.0 / tail.statistic).atan();
        assert!(tail.p_value > 0.0);
        close(tail.p_value / exact_cauchy_tail, 1.0);
        Ok(())
    }

    #[test]
    fn hypothesis_tests_keep_large_counts_and_measurements_finite() -> Result<(), StatsExtError> {
        let biggest = f64::MAX;
        let reference = t_test_one_sample(&[-2.0, 1.5], 0.0)?;
        let paired = t_test_paired(&[-biggest, biggest], &[biggest, -biggest / 2.0])?;
        close(paired.statistic, reference.statistic);

        let balanced = chi_squared_independence(2, 2, &[1e200; 4])?;
        close(balanced.statistic, 0.0);
        close(balanced.p_value, 1.0);
        let fit = chi_squared_goodness_of_fit(&[1e200, 2e200], &[2e200, 1e200])?;
        close(fit.statistic, 1.5e200);

        let normal = anova_one_way(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]])?;
        let large = anova_one_way(&[&[1e200, 2e200, 3e200], &[4e200, 5e200, 6e200]])?;
        close(large.statistic, normal.statistic);
        close(large.p_value, normal.p_value);
        Ok(())
    }

    #[test]
    fn kendall_merge_count_agrees_with_pairwise_reference() -> Result<(), StatsExtError> {
        let left = [-0.0, 4.0, 3.0, 0.0, -2.0, 3.0, 7.0, 4.0];
        let right = [2.0, -4.0, 5.0, 2.0, 3.0, -0.0, 0.0, 8.0];
        let mut numerator = 0.0;
        let mut left_untied = 0.0_f64;
        let mut right_untied = 0.0_f64;
        for (i, (&x, &y)) in left.iter().zip(&right).enumerate() {
            for (&other_x, &other_y) in left.iter().zip(&right).skip(i + 1) {
                let sign = |a: f64, b: f64| {
                    if a < b {
                        -1.0
                    } else if a > b {
                        1.0
                    } else {
                        0.0
                    }
                };
                let a: f64 = sign(x, other_x);
                let b: f64 = sign(y, other_y);
                left_untied += a.abs();
                right_untied += b.abs();
                numerator += a * b;
            }
        }
        close(
            kendall(&left, &right)?,
            numerator / (left_untied * right_untied).sqrt(),
        );
        assert!(kendall(&[-0.0, 0.0], &[1.0, -1.0]).is_err());
        let sorted: Vec<_> = (0..10_000).map(f64::from).collect();
        let descending: Vec<_> = sorted.iter().rev().copied().collect();
        close(kendall(&sorted, &descending)?, -1.0);
        Ok(())
    }

    #[test]
    fn singular_and_wide_designs_have_well_defined_regularized_solutions()
    -> Result<(), StatsExtError> {
        let x = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
        let y = [2.0, 4.0, 6.0];
        let least_squares = ols(3, 2, &x, &y, false)?;
        for value in least_squares.coefficients {
            close(value, 1.0);
        }
        assert_eq!(
            ridge(3, 2, &x, &y, 0.0, false),
            Err(StatsExtError::Singular)
        );
        let wide = ridge(1, 2, &[1.0, 2.0], &[3.0], 1.0, false)?;
        for (&actual, expected) in wide.coefficients.iter().zip([0.5, 1.0]) {
            close(actual, expected);
        }
        let large = ridge(2, 1, &[1e200, -1e200], &[2e200, -2e200], 1.0, false)?;
        close(large.coefficients.first().copied().unwrap_or_default(), 2.0);
        close(large.r_squared, 1.0);
        Ok(())
    }

    #[test]
    fn large_scale_pca_kde_and_logistic_outputs_stay_meaningful() -> Result<(), StatsExtError> {
        let principal = pca(3, 1, &[1e154, -1e154, 0.0], None)?;
        close(
            principal
                .explained_variance
                .first()
                .copied()
                .unwrap_or_default(),
            1e308,
        );
        assert_eq!(principal.explained_variance_ratio, [1.0]);

        let densities = kde(&[-f64::MAX], &[f64::MAX], Some(f64::MAX))?;
        let density = densities.first().copied().unwrap_or_default();
        close(density * f64::MAX, (-2.0_f64).exp() / (2.0 * PI).sqrt());

        let fit = logistic_regression(
            6,
            1,
            &[-2.0, -1.0, -0.5, 0.5, 1.0, 2.0],
            &[0.0, 0.0, 1.0, 0.0, 1.0, 1.0],
            true,
            100,
            1e-8,
        )?;
        assert!(fit.coefficients.first().is_some_and(|&value| value > 0.5));
        assert!(fit.log_likelihood < 0.0 && fit.log_likelihood.is_finite());
        close(fit.intercept, 0.0);
        Ok(())
    }

    #[test]
    fn resampling_rejects_impossible_storage_and_nonfinite_statistics() {
        let invoked = Cell::new(false);
        let statistic = |_: &[f64]| {
            invoked.set(true);
            Ok(1.0)
        };
        assert!(bootstrap(&[1.0], usize::MAX, &statistic, &mut || 0.0).is_err());
        assert!(!invoked.get());
        let bad = |_: &[f64]| Ok(f64::NAN);
        assert!(bootstrap(&[1.0], 1, &bad, &mut || 0.0).is_err());
        assert!(jackknife(&[1.0, 2.0], &bad).is_err());
    }
}
