//! Hypothesis testing, regression (OLS, ridge, lasso, logistic), PCA,
//! rank correlations (Spearman, Kendall), KDE, and bootstrap / jackknife resampling.

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::distributions::Distribution;
use crate::kernels::robust;
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

fn solve(
    rows: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, StatsExtError> {
    let m = Mat {
        rows,
        cols,
        data: a.to_vec(),
    };
    dense::solve(&m, b).map_err(|e| {
        match e {
            | DenseError::Singular => StatsExtError::Singular,
            | _ => StatsExtError::Domain(format!("solve failed: {e:?}")),
        }
    })
}

fn lstsq(
    rows: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, StatsExtError> {
    let m = Mat {
        rows,
        cols,
        data: a.to_vec(),
    };
    dense::lstsq(&m, b).map_err(|e| StatsExtError::Domain(format!("lstsq failed: {e:?}")))
}

type SvdFactors = (Vec<f64>, Vec<f64>, Vec<f64>);

fn svd(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<SvdFactors, StatsExtError> {
    let m = Mat {
        rows,
        cols,
        data: a.to_vec(),
    };
    let res = dense::svd(&m).map_err(|e| StatsExtError::Domain(format!("svd failed: {e:?}")))?;
    Ok((res.u.data, res.s, res.v.transpose().data))
}

fn quantile(
    values: &[f64],
    p: f64,
) -> Result<f64, StatsExtError> {
    if values.is_empty() {
        return Err(StatsExtError::Empty);
    }
    if !(0.0..=1.0).contains(&p) || p.is_nan() {
        return Err(StatsExtError::Domain(
            "quantile p must be in [0, 1]".to_string(),
        ));
    }
    Ok(robust::quantile(values, p, 7))
}

fn corr(
    x: &[f64],
    y: &[f64],
    _population: bool,
) -> Result<f64, StatsExtError> {
    if x.len() != y.len() {
        return Err(StatsExtError::ShapeMismatch(
            "corr length mismatch".to_string(),
        ));
    }
    if x.len() < 2 {
        return Err(StatsExtError::Domain(
            "corr requires at least 2 observations".to_string(),
        ));
    }
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mx;
        let dy = yi - my;
        num += dx * dy;
        sx += dx * dx;
        sy += dy * dy;
    }
    if sx == 0.0 || sy == 0.0 {
        return Err(StatsExtError::Domain("zero variance in corr".to_string()));
    }
    Ok(num / (sx.sqrt() * sy.sqrt()))
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

fn mean_and_var(values: &[f64]) -> Result<(f64, f64), StatsExtError> {
    finite(values)?;
    let n = values.len();
    if n == 0 {
        return Err(StatsExtError::Empty);
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    if n == 1 {
        return Ok((mean, 0.0));
    }
    let ss: f64 = values.iter().map(|&v| (v - mean) * (v - mean)).sum();
    Ok((mean, ss / (n - 1) as f64))
}

/// Average 1-based ranks with tie averaging, plus `sum(t^3 - t)` for tie groups.
fn average_ranks_and_tie_sum(values: &[f64]) -> (Vec<f64>, f64) {
    let n = values.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));

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
    let dist = Distribution::StudentT { df };
    (2.0 * dist.cdf(-t_stat.abs())).clamp(0.0, 1.0)
}

/// One-sample Student's t-test against population mean `mu0`.
///
/// # Errors
///
/// [`StatsExtError::Domain`] when `sample.len() < 2` or the sample variance is zero.
pub fn t_test_one_sample(
    sample: &[f64],
    mu0: f64,
) -> Result<TestResult, StatsExtError> {
    if sample.len() < 2 {
        return Err(StatsExtError::Domain(
            "t_test_one_sample requires at least 2 observations".to_string(),
        ));
    }
    let (mean, var) = mean_and_var(sample)?;
    if var == 0.0 {
        return Err(StatsExtError::Domain(
            "t_test_one_sample is undefined for zero-variance samples".to_string(),
        ));
    }
    let n = sample.len() as f64;
    let se = (var / n).sqrt();
    let t_stat = (mean - mu0) / se;
    let df = n - 1.0;
    Ok(TestResult {
        statistic: t_stat,
        p_value: two_sided_t_p_value(t_stat, df),
        df: Some(df),
    })
}

/// Two-sample t-test (Student's pooled-variance when `equal_var` is `true`,
/// Welch's unequal-variance otherwise).
///
/// # Errors
///
/// [`StatsExtError::Domain`] when either sample has fewer than 2 elements or
/// the combined standard error is zero.
pub fn t_test_two_sample(
    left: &[f64],
    right: &[f64],
    equal_var: bool,
) -> Result<TestResult, StatsExtError> {
    if left.len() < 2 || right.len() < 2 {
        return Err(StatsExtError::Domain(
            "t_test_two_sample requires at least 2 observations per sample".to_string(),
        ));
    }
    let (m1, v1) = mean_and_var(left)?;
    let (m2, v2) = mean_and_var(right)?;
    let n1 = left.len() as f64;
    let n2 = right.len() as f64;

    let (se, df) = if equal_var {
        let df = n1 + n2 - 2.0;
        let pooled = ((n1 - 1.0) * v1 + (n2 - 1.0) * v2) / df;
        ((pooled * (1.0 / n1 + 1.0 / n2)).sqrt(), df)
    } else {
        let s1 = v1 / n1;
        let s2 = v2 / n2;
        let sum = s1 + s2;
        if sum == 0.0 {
            return Err(StatsExtError::Domain(
                "t_test_two_sample is undefined when both samples have zero variance".to_string(),
            ));
        }
        let denom = (s1 * s1) / (n1 - 1.0) + (s2 * s2) / (n2 - 1.0);
        let df = (sum * sum) / denom;
        (sum.sqrt(), df)
    };

    if se == 0.0 {
        return Err(StatsExtError::Domain(
            "t_test_two_sample is undefined when standard error is zero".to_string(),
        ));
    }
    let t_stat = (m1 - m2) / se;
    Ok(TestResult {
        statistic: t_stat,
        p_value: two_sided_t_p_value(t_stat, df),
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
    let diffs: Vec<f64> = left.iter().zip(right).map(|(a, b)| a - b).collect();
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
        stat += (diff * diff) / e;
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
    let mut row_sums = vec![0.0; rows];
    let mut col_sums = vec![0.0; cols];
    let mut total = 0.0;
    for r in 0..rows {
        for c in 0..cols {
            let v = observed[r * cols + c];
            if v < 0.0 {
                return Err(StatsExtError::Domain(
                    "contingency table counts cannot be negative".to_string(),
                ));
            }
            row_sums[r] += v;
            col_sums[c] += v;
            total += v;
        }
    }
    if total <= 0.0 || row_sums.iter().any(|&s| s <= 0.0) || col_sums.iter().any(|&s| s <= 0.0) {
        return Err(StatsExtError::Domain(
            "contingency table marginals must all be positive".to_string(),
        ));
    }
    let mut stat = 0.0;
    for r in 0..rows {
        for c in 0..cols {
            let expected = row_sums[r] * col_sums[c] / total;
            let diff = observed[r * cols + c] - expected;
            stat += (diff * diff) / expected;
        }
    }
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
    let mut grand_sum = 0.0;
    for g in groups {
        finite(g)?;
        if g.is_empty() {
            return Err(StatsExtError::Domain(
                "anova_one_way groups must be non-empty".to_string(),
            ));
        }
        n_total += g.len();
        grand_sum += g.iter().sum::<f64>();
    }
    if n_total <= k {
        return Err(StatsExtError::Domain(
            "anova_one_way requires more total observations than groups".to_string(),
        ));
    }
    let grand_mean = grand_sum / n_total as f64;
    let mut ss_between = 0.0;
    let mut ss_within = 0.0;
    for g in groups {
        let g_mean = g.iter().sum::<f64>() / g.len() as f64;
        ss_between += g.len() as f64 * (g_mean - grand_mean) * (g_mean - grand_mean);
        for &x in *g {
            ss_within += (x - g_mean) * (x - g_mean);
        }
    }
    let df1 = (k - 1) as f64;
    let df2 = (n_total - k) as f64;
    let ms_within = ss_within / df2;
    if ms_within == 0.0 {
        return Err(StatsExtError::Domain(
            "anova_one_way is undefined when within-group variance is zero".to_string(),
        ));
    }
    let f_stat = (ss_between / df1) / ms_within;
    let p_value =
        crate::kernels::special::regularized_beta(df2 / (df2 + df1 * f_stat), df2 * 0.5, df1 * 0.5)
            .clamp(0.0, 1.0);
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

fn center_xy(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    fit_intercept: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64) {
    if !fit_intercept {
        return (x.to_vec(), y.to_vec(), vec![0.0; cols], 0.0);
    }
    let mut x_mean = vec![0.0; cols];
    for r in 0..rows {
        for c in 0..cols {
            x_mean[c] += x[r * cols + c];
        }
    }
    for m in &mut x_mean {
        *m /= rows as f64;
    }
    let y_mean = y.iter().sum::<f64>() / rows as f64;
    let mut xc = vec![0.0; rows * cols];
    for r in 0..rows {
        for c in 0..cols {
            xc[r * cols + c] = x[r * cols + c] - x_mean[c];
        }
    }
    let yc: Vec<f64> = y.iter().map(|&v| v - y_mean).collect();
    (xc, yc, x_mean, y_mean)
}

fn build_linear_result(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    coefficients: Vec<f64>,
    x_mean: &[f64],
    y_mean: f64,
    fit_intercept: bool,
) -> LinearRegressionResult {
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
    let mut residuals = Vec::with_capacity(rows);
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    let raw_y_mean = y.iter().sum::<f64>() / rows as f64;
    for r in 0..rows {
        let mut pred = intercept;
        for c in 0..cols {
            pred += x[r * cols + c] * coefficients[c];
        }
        let res = y[r] - pred;
        residuals.push(res);
        ss_res += res * res;
        ss_tot += (y[r] - raw_y_mean) * (y[r] - raw_y_mean);
    }
    let r_squared = if ss_tot == 0.0 {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    };
    LinearRegressionResult {
        coefficients,
        intercept,
        r_squared,
        residuals,
    }
}

/// Ordinary Least Squares linear regression using SVD least-squares.
///
/// # Errors
///
/// [`StatsExtError::Empty`], [`StatsExtError::ShapeMismatch`], or [`StatsExtError::NotConverged`].
pub fn ols(
    rows: usize,
    cols: usize,
    x: &[f64],
    y: &[f64],
    fit_intercept: bool,
) -> Result<LinearRegressionResult, StatsExtError> {
    validate_regression_inputs(rows, cols, x, y)?;
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept);
    let coefficients = lstsq(rows, cols, &xc, &yc)?;
    Ok(build_linear_result(
        rows,
        cols,
        x,
        y,
        coefficients,
        &x_mean,
        y_mean,
        fit_intercept,
    ))
}

/// Ridge (L2-regularized) linear regression minimizing `||y - Xβ||₂² + alpha ||β||₂²`.
///
/// # Errors
///
/// [`StatsExtError::Empty`], [`StatsExtError::ShapeMismatch`], [`StatsExtError::Domain`]
/// when `alpha < 0`, or [`StatsExtError::Singular`].
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
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept);
    let mut xtx = vec![0.0; cols * cols];
    let mut xty = vec![0.0; cols];
    for r in 0..rows {
        for c1 in 0..cols {
            let v1 = xc[r * cols + c1];
            xty[c1] += v1 * yc[r];
            for c2 in 0..cols {
                xtx[c1 * cols + c2] += v1 * xc[r * cols + c2];
            }
        }
    }
    for c in 0..cols {
        xtx[c * cols + c] += alpha;
    }
    let coefficients = solve(cols, cols, &xtx, &xty)?;
    Ok(build_linear_result(
        rows,
        cols,
        x,
        y,
        coefficients,
        &x_mean,
        y_mean,
        fit_intercept,
    ))
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
    let (xc, yc, x_mean, y_mean) = center_xy(rows, cols, x, y, fit_intercept);
    let n = rows as f64;
    let mut col_sq_norms = vec![0.0; cols];
    for r in 0..rows {
        for c in 0..cols {
            let v = xc[r * cols + c];
            col_sq_norms[c] += v * v;
        }
    }

    let mut beta = vec![0.0; cols];
    let mut residual = yc;
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
    Ok(build_linear_result(
        rows,
        cols,
        x,
        y,
        beta,
        &x_mean,
        y_mean,
        fit_intercept,
    ))
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

    let p = if fit_intercept {
        cols + 1
    } else {
        cols
    };
    let mut w = vec![0.0; p];
    let mut converged = false;
    let mut iterations = 0usize;

    let feature = |r: usize, j: usize| -> f64 {
        if fit_intercept && j == 0 {
            1.0
        } else if fit_intercept {
            x[r * cols + (j - 1)]
        } else {
            x[r * cols + j]
        }
    };

    for iter in 1..=max_iter {
        iterations = iter;
        let mut grad = vec![0.0; p];
        let mut hess = vec![0.0; p * p];

        for r in 0..rows {
            let mut eta = 0.0;
            for j in 0..p {
                eta += w[j] * feature(r, j);
            }
            let prob = 1.0 / (1.0 + (-eta).exp());
            let weight = (prob * (1.0 - prob)).max(1e-12);
            let err = prob - y[r];
            for j in 0..p {
                let xj = feature(r, j);
                grad[j] += xj * err;
                for k in 0..p {
                    hess[j * p + k] += weight * xj * feature(r, k);
                }
            }
        }
        // Small ridge damping for numerical stability near separation.
        for j in 0..p {
            hess[j * p + j] += 1e-8;
        }
        let step = solve(p, p, &hess, &grad)?;
        let mut max_step = 0.0f64;
        for j in 0..p {
            w[j] -= step[j];
            max_step = max_step.max(step[j].abs());
        }
        if max_step < tol {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(StatsExtError::NotConverged);
    }

    let mut log_likelihood = 0.0;
    for r in 0..rows {
        let mut eta = 0.0;
        for j in 0..p {
            eta += w[j] * feature(r, j);
        }
        let prob = (1.0 / (1.0 + (-eta).exp())).clamp(1e-15, 1.0 - 1e-15);
        log_likelihood += y[r] * prob.ln() + (1.0 - y[r]) * (1.0 - prob).ln();
    }

    let (intercept, coefficients) = if fit_intercept {
        (w[0], w[1..].to_vec())
    } else {
        (0.0, w)
    };
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

    let (centered, _, mean, _) = center_xy(rows, cols, data, &vec![0.0; rows], true);
    let (u, s, vt) = svd(rows, cols, &centered)?;

    let denom = if rows > 1 {
        (rows - 1) as f64
    } else {
        1.0
    };
    let all_var: Vec<f64> = s.iter().map(|&sv| (sv * sv) / denom).collect();
    let total_var: f64 = all_var.iter().sum();

    let explained_variance = all_var[..k].to_vec();
    let explained_variance_ratio = explained_variance
        .iter()
        .map(|&v| {
            if total_var > 0.0 {
                v / total_var
            } else {
                0.0
            }
        })
        .collect();
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
    let (r_left, _) = average_ranks_and_tie_sum(left);
    let (r_right, _) = average_ranks_and_tie_sum(right);
    corr(&r_left, &r_right, false)
}

/// Kendall's tau-b rank correlation coefficient (adjusted for ties).
///
/// # Errors
///
/// [`StatsExtError::ShapeMismatch`] when lengths differ, [`StatsExtError::Empty`]
/// when empty, or [`StatsExtError::Domain`] when all values in either sample are tied.
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
    let n = left.len();
    if n == 0 {
        return Err(StatsExtError::Empty);
    }
    if n < 2 {
        return Err(StatsExtError::Domain(
            "kendall requires at least 2 observations".to_string(),
        ));
    }
    let mut concordant = 0i64;
    let mut discordant = 0i64;
    let mut ties_x = 0i64;
    let mut ties_y = 0i64;
    let mut ties_xy = 0i64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = left[i].total_cmp(&left[j]);
            let dy = right[i].total_cmp(&right[j]);
            match (dx, dy) {
                | (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => ties_xy += 1,
                | (std::cmp::Ordering::Equal, _) => ties_x += 1,
                | (_, std::cmp::Ordering::Equal) => ties_y += 1,
                | (a, b) if a == b => concordant += 1,
                | _ => discordant += 1,
            }
        }
    }
    let n0 = (n * (n - 1) / 2) as f64;
    let n1 = (ties_x + ties_xy) as f64;
    let n2 = (ties_y + ties_xy) as f64;
    let denom = ((n0 - n1) * (n0 - n2)).sqrt();
    if denom == 0.0 {
        return Err(StatsExtError::Domain(
            "kendall is undefined when all elements in a sample are tied".to_string(),
        ));
    }
    Ok(((concordant - discordant) as f64 / denom).clamp(-1.0, 1.0))
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
            let (_, var) = mean_and_var(samples)?;
            let std = var.sqrt();
            let q75 = quantile(samples, 0.75)?;
            let q25 = quantile(samples, 0.25)?;
            let iqr = q75 - q25;
            let scale = if iqr > 0.0 {
                std.min(iqr / 1.34)
            } else {
                std
            };
            let bw = 0.9 * scale * (samples.len() as f64).powf(-0.2);
            if bw <= 0.0 || !bw.is_finite() {
                return Err(StatsExtError::Domain(
                    "cannot infer kde bandwidth for a zero-spread sample".to_string(),
                ));
            }
            bw
        },
    };

    let n = samples.len() as f64;
    let norm = 1.0 / (n * h * (2.0 * PI).sqrt());
    let mut densities = Vec::with_capacity(points.len());
    for &x in points {
        let sum: f64 = samples
            .iter()
            .map(|&xi| {
                let u = (x - xi) / h;
                (-0.5 * u * u).exp()
            })
            .sum();
        densities.push(norm * sum);
    }
    Ok(densities)
}

// ----- resampling (bootstrap & jackknife) -----------------------------------

/// Non-parametric bootstrap resampling with `n_resamples` draws.
///
/// # Errors
///
/// [`StatsExtError::Empty`] when `samples` is empty, [`StatsExtError::Domain`]
/// when `n_resamples == 0`, or any error returned by `statistic`.
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
    let estimate = statistic(samples)?;
    let n = samples.len();
    let mut replicates = Vec::with_capacity(n_resamples);
    let mut resample = vec![0.0; n];

    for _ in 0..n_resamples {
        for slot in &mut resample {
            let draw = uniform();
            if !(0.0..1.0).contains(&draw) {
                return Err(StatsExtError::Domain(
                    "uniform source must return values in [0, 1)".to_string(),
                ));
            }
            let idx = (draw * n as f64).floor() as usize;
            *slot = samples[idx];
        }
        replicates.push(statistic(&resample)?);
    }

    let (mean, var) = mean_and_var(&replicates)?;
    let ci_low = quantile(&replicates, 0.025)?;
    let ci_high = quantile(&replicates, 0.975)?;

    Ok(BootstrapResult {
        estimate,
        mean,
        std_error: var.sqrt(),
        ci_low,
        ci_high,
        replicates,
    })
}

/// Leave-one-out jackknife resampling for bias and standard error estimation.
///
/// # Errors
///
/// [`StatsExtError::Domain`] when `samples.len() < 2`, or any error returned by `statistic`.
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
    let estimate = statistic(samples)?;
    let mut jackknife_estimates = Vec::with_capacity(n);
    let mut loo = Vec::with_capacity(n - 1);

    for i in 0..n {
        loo.clear();
        loo.extend_from_slice(&samples[..i]);
        loo.extend_from_slice(&samples[i + 1..]);
        jackknife_estimates.push(statistic(&loo)?);
    }

    let nf = n as f64;
    let mean_jack = jackknife_estimates.iter().sum::<f64>() / nf;
    let bias = (nf - 1.0) * (mean_jack - estimate);
    let ss: f64 = jackknife_estimates
        .iter()
        .map(|&v| (v - mean_jack) * (v - mean_jack))
        .sum();
    let std_error = (((nf - 1.0) / nf) * ss).sqrt();

    Ok(JackknifeResult {
        estimate,
        bias,
        std_error,
        jackknife_estimates,
    })
}
