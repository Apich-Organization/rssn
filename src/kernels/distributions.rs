//! Univariate and multivariate probability distributions.

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::special::erfc_numerical as erfc;
use crate::kernels::special::gamma_numerical as gamma;
use crate::kernels::special::ln_gamma_numerical as ln_gamma;
use crate::kernels::special::regularized_beta;
use crate::kernels::special::regularized_lower_gamma as regularized_gamma_p;
use crate::kernels::special::regularized_upper_gamma as regularized_gamma_q;
use statrs::distribution::DiscreteCDF;
use std::f64::consts::PI;

/// Error returned by distribution operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DistError {
    /// Empty input.
    Empty,
    /// Mismatched dimensions or vector lengths.
    ShapeMismatch(String),
    /// Invalid distribution parameter or argument domain.
    Domain(String),
    /// Matrix is not positive definite.
    NotPositiveDefinite,
}

impl std::fmt::Display for DistError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("empty input"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
            | Self::NotPositiveDefinite => f.write_str("matrix is not positive definite"),
        }
    }
}

impl std::error::Error for DistError {}

fn cholesky(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<Vec<f64>, DistError> {
    if rows != cols || data.len() != rows * cols {
        return Err(DistError::ShapeMismatch(
            "expected square matrix".to_string(),
        ));
    }
    for i in 0..rows {
        for j in 0..cols {
            if !data[i * cols + j].is_finite()
                || (data[i * cols + j] - data[j * cols + i]).abs()
                    > 1e-12
                        * data[i * cols + j]
                            .abs()
                            .max(data[j * cols + i].abs())
                            .max(1.0)
            {
                return Err(DistError::NotPositiveDefinite);
            }
        }
    }
    let m = Mat {
        rows,
        cols,
        data: data.to_vec(),
    };
    dense::cholesky(&m)
        .map(|l| l.data)
        .map_err(|_| DistError::NotPositiveDefinite)
}


// Bounded integer bisection avoids unbounded mass summation and overflow in
// a doubling search for quantiles outside the u64 count range.
fn discrete_quantile(
    cdf: impl Fn(u64) -> f64,
    u: f64,
    maximum: u64,
) -> Result<f64, DistError> {
    if u <= cdf(0) {
        return Ok(0.0);
    }
    let mut low = 0u64;
    let mut high = maximum.min(2);
    while cdf(high) < u {
        if high == maximum {
            return Err(DistError::Domain(
                "quantile exceeds the supported integer range".to_string(),
            ));
        }
        low = high;
        high = high.saturating_mul(2).min(maximum);
    }
    if !cdf(high).is_finite() {
        return Err(DistError::Domain(
            "distribution CDF did not produce a finite probability".to_string(),
        ));
    }
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if cdf(middle) >= u {
            high = middle;
        } else {
            low = middle;
        }
    }
    Ok(high as f64)
}

/// One of the univariate distributions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Distribution {
    /// `Normal<d> = { mu, sigma }`.
    Normal {
        /// The mean.
        mu: f64,
        /// The standard deviation.
        sigma: f64,
    },
    /// `Uniform<d> = { lo, hi }`.
    Uniform {
        /// The lower bound.
        lo: f64,
        /// The upper bound.
        hi: f64,
    },
    /// `Exponential<d> = { rate }`.
    Exponential {
        /// The rate, in events per unit.
        rate: f64,
    },
    /// `Gamma<d> = { shape, rate }`.
    Gamma {
        /// The shape `k`.
        shape: f64,
        /// The rate `λ`.
        rate: f64,
    },
    /// `Beta = { alpha, beta }`.
    Beta {
        /// The first shape.
        alpha: f64,
        /// The second shape.
        beta: f64,
    },
    /// `ChiSquared = { df }`.
    ChiSquared {
        /// The degrees of freedom.
        df: f64,
    },
    /// `StudentT = { df }`.
    StudentT {
        /// The degrees of freedom.
        df: f64,
    },
    /// `Poisson = { lambda }`.
    Poisson {
        /// The mean rate.
        lambda: f64,
    },
    /// `Binomial = { n, p }`.
    Binomial {
        /// The number of trials.
        n: u64,
        /// The success probability.
        p: f64,
    },
    /// `LogNormal<d> = { mu, sigma }`.
    LogNormal {
        /// The mean of the natural logarithm.
        mu: f64,
        /// The standard deviation of the natural logarithm.
        sigma: f64,
    },
    /// `Weibull<d> = { shape, scale }`.
    Weibull {
        /// The shape parameter `k > 0`.
        shape: f64,
        /// The scale parameter `λ > 0`.
        scale: f64,
    },
    /// `Cauchy<d> = { loc, scale }`.
    Cauchy {
        /// The location parameter `x₀`.
        loc: f64,
        /// The scale parameter `γ > 0`.
        scale: f64,
    },
    /// `F = { df1, df2 }` (Snedecor's F-distribution).
    F {
        /// The numerator degrees of freedom.
        df1: f64,
        /// The denominator degrees of freedom.
        df2: f64,
    },
    /// `NegativeBinomial = { r, p }` (number of failures before `r` successes).
    NegativeBinomial {
        /// The number of successes `r > 0`.
        r: f64,
        /// The success probability `0 < p <= 1`.
        p: f64,
    },
    /// `Geometric = { p }` (number of failures before the first success).
    Geometric {
        /// The success probability `0 < p <= 1`.
        p: f64,
    },
}

impl Distribution {
    fn validate(&self) -> Result<(), DistError> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        let valid = match *self {
            | Self::Normal { mu, sigma } | Self::LogNormal { mu, sigma } => {
                mu.is_finite() && positive(sigma)
            },
            | Self::Uniform { lo, hi } => lo.is_finite() && hi.is_finite() && lo < hi,
            | Self::Exponential { rate } => positive(rate),
            | Self::Gamma { shape, rate } => positive(shape) && positive(rate),
            | Self::Beta { alpha, beta } => positive(alpha) && positive(beta),
            | Self::ChiSquared { df } | Self::StudentT { df } => positive(df),
            | Self::Poisson { lambda } => lambda.is_finite() && lambda >= 0.0,
            | Self::Binomial { p, .. } => (0.0..=1.0).contains(&p),
            | Self::Weibull { shape, scale } => positive(shape) && positive(scale),
            | Self::Cauchy { loc, scale } => loc.is_finite() && positive(scale),
            | Self::F { df1, df2 } => positive(df1) && positive(df2),
            | Self::NegativeBinomial { r, p } => positive(r) && positive(p) && p <= 1.0,
            | Self::Geometric { p } => positive(p) && p <= 1.0,
        };
        if valid {
            Ok(())
        } else {
            Err(DistError::Domain(
                "invalid distribution parameters".to_string(),
            ))
        }
    }

    /// The probability density (or mass) at `x`.
    #[must_use]
    pub fn pdf(
        &self,
        x: f64,
    ) -> f64 {
        if self.validate().is_err() || x.is_nan() {
            return f64::NAN;
        }
        if x.is_infinite() {
            return 0.0;
        }
        match *self {
            | Self::Normal { mu, sigma } => {
                if sigma <= 0.0 {
                    return f64::NAN;
                }
                let z = (x - mu) / sigma;
                (-0.5 * z * z).exp() / (sigma * (2.0 * PI).sqrt())
            },
            | Self::Uniform { lo, hi } => {
                if hi <= lo {
                    f64::NAN
                } else if x >= lo && x <= hi {
                    1.0 / (hi - lo)
                } else {
                    0.0
                }
            },
            | Self::Exponential { rate } => {
                if rate <= 0.0 {
                    f64::NAN
                } else if x < 0.0 {
                    0.0
                } else {
                    rate * (-rate * x).exp()
                }
            },
            | Self::Gamma { shape, rate } => {
                if shape <= 0.0 || rate <= 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if x == 0.0 {
                    return if shape < 1.0 {
                        f64::INFINITY
                    } else if crate::kernels::exact_eq(shape, 1.0) {
                        rate
                    } else {
                        0.0
                    };
                }
                ((shape - 1.0) * x.ln() - rate * x + shape * rate.ln() - ln_gamma(shape)).exp()
            },
            | Self::Beta { alpha, beta } => {
                if alpha <= 0.0 || beta <= 0.0 {
                    return f64::NAN;
                }
                if !(0.0..=1.0).contains(&x) {
                    return 0.0;
                }
                let log =
                    (alpha - 1.0) * x.ln() + (beta - 1.0) * (1.0 - x).ln() + ln_gamma(alpha + beta)
                        - ln_gamma(alpha)
                        - ln_gamma(beta);
                log.exp()
            },
            | Self::ChiSquared { df } => {
                if df <= 0.0 {
                    return f64::NAN;
                }
                Self::Gamma {
                    shape: df / 2.0,
                    rate: 0.5,
                }
                .pdf(x)
            },
            | Self::StudentT { df } => {
                if df <= 0.0 {
                    return f64::NAN;
                }
                (ln_gamma(f64::midpoint(df, 1.0))
                    - ln_gamma(df / 2.0)
                    - 0.5 * (df * PI).ln()
                    - f64::midpoint(df, 1.0) * (x * x / df).ln_1p())
                .exp()
            },
            | Self::Poisson { lambda } => {
                if lambda < 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 || x.fract() != 0.0 {
                    return 0.0;
                }
                if lambda == 0.0 {
                    return if x == 0.0 { 1.0 } else { 0.0 };
                }
                (-lambda + x * lambda.ln() - ln_gamma(x + 1.0)).exp()
            },
            | Self::Binomial { n, p } => {
                if !(0.0..=1.0).contains(&p) {
                    return f64::NAN;
                }
                if x < 0.0 || x.fract() != 0.0 || x > n as f64 {
                    return 0.0;
                }
                let k = x;
                (ln_gamma(n as f64 + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n as f64 - k + 1.0)
                    + k * p.ln()
                    + (n as f64 - k) * (1.0 - p).ln())
                .exp()
            },
            | Self::LogNormal { mu, sigma } => {
                if sigma <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    return 0.0;
                }
                let z = (x.ln() - mu) / sigma;
                (-0.5 * z * z).exp() / (x * sigma * (2.0 * PI).sqrt())
            },
            | Self::Weibull { shape, scale } => {
                if shape <= 0.0 || scale <= 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if x == 0.0 {
                    return if shape < 1.0 {
                        f64::INFINITY
                    } else if crate::kernels::exact_eq(shape, 1.0) {
                        1.0 / scale
                    } else {
                        0.0
                    };
                }
                let ratio = x / scale;
                (shape / scale) * ratio.powf(shape - 1.0) * (-ratio.powf(shape)).exp()
            },
            | Self::Cauchy { loc, scale } => {
                if scale <= 0.0 {
                    return f64::NAN;
                }
                let z = (x - loc) / scale;
                1.0 / (PI * scale * (1.0 + z * z))
            },
            | Self::F { df1, df2 } => {
                if df1 <= 0.0 || df2 <= 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if x == 0.0 {
                    return if df1 < 2.0 {
                        f64::INFINITY
                    } else if crate::kernels::exact_eq(df1, 2.0) {
                        1.0
                    } else {
                        0.0
                    };
                }
                (f64::midpoint(df1 * df1.ln(), df2 * df2.ln()) + (0.5 * df1 - 1.0) * x.ln()
                    - f64::midpoint(df1, df2) * (df1 * x + df2).ln()
                    + ln_gamma(f64::midpoint(df1, df2))
                    - ln_gamma(0.5 * df1)
                    - ln_gamma(0.5 * df2))
                .exp()
            },
            | Self::NegativeBinomial { r, p } => {
                if r <= 0.0 || p <= 0.0 || p > 1.0 {
                    return f64::NAN;
                }
                if x < 0.0 || x.fract() != 0.0 {
                    return 0.0;
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    return if x == 0.0 { 1.0 } else { 0.0 };
                }
                (ln_gamma(x + r) - ln_gamma(x + 1.0) - ln_gamma(r)
                    + r * p.ln()
                    + x * (1.0 - p).ln())
                .exp()
            },
            | Self::Geometric { p } => {
                if p <= 0.0 || p > 1.0 {
                    return f64::NAN;
                }
                if x < 0.0 || x.fract() != 0.0 {
                    return 0.0;
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    return if x == 0.0 { 1.0 } else { 0.0 };
                }
                p * (1.0 - p).powf(x)
            },
        }
    }

    /// The cumulative distribution function at `x`.
    #[must_use]
    pub fn cdf(
        &self,
        x: f64,
    ) -> f64 {
        if self.validate().is_err() || x.is_nan() {
            return f64::NAN;
        }
        if x == f64::NEG_INFINITY {
            return 0.0;
        }
        if x == f64::INFINITY {
            return 1.0;
        }
        match *self {
            | Self::Normal { mu, sigma } => {
                if sigma <= 0.0 {
                    return f64::NAN;
                }
                let z = (x - mu) / sigma;
                0.5 * erfc(-z / std::f64::consts::SQRT_2)
            },
            | Self::Uniform { lo, hi } => {
                if hi <= lo {
                    return f64::NAN;
                }
                ((x - lo) / (hi - lo)).clamp(0.0, 1.0)
            },
            | Self::Exponential { rate } => {
                if rate <= 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    0.0
                } else {
                    -(-rate * x).exp_m1()
                }
            },
            | Self::Gamma { shape, rate } => {
                if shape <= 0.0 || rate <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    0.0
                } else {
                    regularized_gamma_p(shape, rate * x)
                }
            },
            | Self::Beta { alpha, beta } => regularized_beta(x.clamp(0.0, 1.0), alpha, beta),
            | Self::ChiSquared { df } => {
                if df <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    0.0
                } else {
                    regularized_gamma_p(df / 2.0, x / 2.0)
                }
            },
            | Self::StudentT { df } => {
                if df <= 0.0 {
                    return f64::NAN;
                }
                let tail = 0.5 * regularized_beta(df / (df + x * x), df / 2.0, 0.5);
                if x >= 0.0 { 1.0 - tail } else { tail }
            },
            | Self::Poisson { lambda } => {
                if lambda < 0.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                // P(X <= k) = Q(k+1, lambda).
                if x == f64::INFINITY {
                    1.0
                } else {
                    regularized_gamma_q(x.floor() + 1.0, lambda)
                }
            },
            | Self::Binomial { n, p } => {
                if !(0.0..=1.0).contains(&p) {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if x >= n as f64 {
                    return 1.0;
                }
                // P(X <= k) = I_{1-p}(n-k, k+1).
                regularized_beta(1.0 - p, n as f64 - x.floor(), x.floor() + 1.0)
            },
            | Self::LogNormal { mu, sigma } => {
                if sigma <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    return 0.0;
                }
                Self::Normal { mu, sigma }.cdf(x.ln())
            },
            | Self::Weibull { shape, scale } => {
                if shape <= 0.0 || scale <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    0.0
                } else {
                    -(-(x / scale).powf(shape)).exp_m1()
                }
            },
            | Self::Cauchy { loc, scale } => {
                if scale <= 0.0 {
                    return f64::NAN;
                }
                0.5 + ((x - loc) / scale).atan() / PI
            },
            | Self::F { df1, df2 } => {
                if df1 <= 0.0 || df2 <= 0.0 {
                    return f64::NAN;
                }
                if x <= 0.0 {
                    0.0
                } else {
                    let arg = (df1 * x) / (df1 * x + df2);
                    regularized_beta(arg, 0.5 * df1, 0.5 * df2)
                }
            },
            | Self::NegativeBinomial { r, p } => {
                if r <= 0.0 || p <= 0.0 || p > 1.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    return 1.0;
                }
                regularized_beta(p, r, x.floor() + 1.0)
            },
            | Self::Geometric { p } => {
                if p <= 0.0 || p > 1.0 {
                    return f64::NAN;
                }
                if x < 0.0 {
                    return 0.0;
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    return 1.0;
                }
                -((x.floor() + 1.0) * (-p).ln_1p()).exp_m1()
            },
        }
    }

    /// The mean, when it exists.
    #[must_use]
    pub fn mean(&self) -> f64 {
        if self.validate().is_err() {
            return f64::NAN;
        }
        match *self {
            | Self::Normal { mu, .. } => mu,
            | Self::Uniform { lo, hi } => f64::midpoint(lo, hi),
            | Self::Exponential { rate } => 1.0 / rate,
            | Self::Gamma { shape, rate } => shape / rate,
            | Self::Beta { alpha, beta } => alpha / (alpha + beta),
            | Self::ChiSquared { df } => df,
            | Self::StudentT { df } => {
                if df > 1.0 {
                    0.0
                } else {
                    f64::NAN
                }
            },
            | Self::Poisson { lambda } => lambda,
            | Self::Binomial { n, p } => n as f64 * p,
            | Self::LogNormal { mu, sigma } => (mu + 0.5 * sigma * sigma).exp(),
            | Self::Weibull { shape, scale } => scale * gamma(1.0 + 1.0 / shape),
            | Self::Cauchy { .. } => f64::NAN,
            | Self::F { df2, .. } => {
                if df2 > 2.0 {
                    df2 / (df2 - 2.0)
                } else {
                    f64::NAN
                }
            },
            | Self::NegativeBinomial { r, p } => r * (1.0 - p) / p,
            | Self::Geometric { p } => (1.0 - p) / p,
        }
    }

    /// The variance, when it exists.
    #[must_use]
    pub fn variance(&self) -> f64 {
        if self.validate().is_err() {
            return f64::NAN;
        }
        match *self {
            | Self::Normal { sigma, .. } => sigma * sigma,
            | Self::Uniform { lo, hi } => (hi - lo) * (hi - lo) / 12.0,
            | Self::Exponential { rate } => 1.0 / (rate * rate),
            | Self::Gamma { shape, rate } => shape / (rate * rate),
            | Self::Beta { alpha, beta } => {
                alpha * beta / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0))
            },
            | Self::ChiSquared { df } => 2.0 * df,
            | Self::StudentT { df } => {
                if df > 2.0 {
                    df / (df - 2.0)
                } else {
                    f64::NAN
                }
            },
            | Self::Poisson { lambda } => lambda,
            | Self::Binomial { n, p } => n as f64 * p * (1.0 - p),
            | Self::LogNormal { mu, sigma } => {
                let s2 = sigma * sigma;
                s2.exp_m1() * (2.0 * mu + s2).exp()
            },
            | Self::Weibull { shape, scale } => {
                let g1 = gamma(1.0 + 1.0 / shape);
                let g2 = gamma(1.0 + 2.0 / shape);
                scale * scale * (g2 - g1 * g1)
            },
            | Self::Cauchy { .. } => f64::NAN,
            | Self::F { df1, df2 } => {
                if df2 > 4.0 {
                    2.0 * df2 * df2 * (df1 + df2 - 2.0)
                        / (df1 * (df2 - 2.0) * (df2 - 2.0) * (df2 - 4.0))
                } else {
                    f64::NAN
                }
            },
            | Self::NegativeBinomial { r, p } => r * (1.0 - p) / (p * p),
            | Self::Geometric { p } => (1.0 - p) / (p * p),
        }
    }

    /// The `u`-quantile, `0 ≤ u ≤ 1`: the inverse of [`Distribution::cdf`].
    ///
    /// This is what `sample` uses at run time, because the interpreter's RNG
    /// hands out one uniform per action (§6.10) and an inverse-CDF sample needs
    /// exactly one. Continuous distributions invert by bisection on the CDF,
    /// which converges to machine precision and needs no
    /// distribution-specific algebra; discrete ones invert their CDF using bracketed integer searches.
    ///
    /// # Errors
    ///
    /// [`DistError::Domain`] when `u` is outside `[0, 1]` or a
    /// parameter is unusable.
    pub fn quantile(
        &self,
        u: f64,
    ) -> Result<f64, DistError> {
        self.validate()?;
        if !(0.0..=1.0).contains(&u) {
            return Err(DistError::Domain(
                "a quantile needs `0 <= u <= 1`".to_string(),
            ));
        }
        match *self {
            | Self::Normal { mu, sigma } => {
                Ok(mu
                    - sigma
                        * std::f64::consts::SQRT_2
                        * crate::kernels::special::inverse_erfc(2.0 * u))
            },
            | Self::LogNormal { mu, sigma } => Ok(Self::Normal { mu, sigma }.quantile(u)?.exp()),
            | Self::Uniform { lo, hi } => Ok(lo * (1.0 - u) + hi * u),
            | Self::Exponential { rate } => Ok(-(-u).ln_1p() / rate),

            | Self::Poisson { lambda } => {
                if lambda == 0.0 || u == 0.0 {
                    return Ok(0.0);
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(f64::INFINITY);
                }
                let distribution = statrs::distribution::Poisson::new(lambda)
                    .map_err(|e| DistError::Domain(e.to_string()))?;
                discrete_quantile(|k| distribution.cdf(k), u, u64::MAX)
            },
            | Self::Binomial { n, p } => {
                let distribution = statrs::distribution::Binomial::new(p, n)
                    .map_err(|e| DistError::Domain(e.to_string()))?;
                discrete_quantile(|k| distribution.cdf(k), u, n)
            },
            | Self::NegativeBinomial { r, p } => {
                if crate::kernels::exact_eq(p, 1.0) || u == 0.0 {
                    return Ok(0.0);
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(f64::INFINITY);
                }
                let distribution = statrs::distribution::NegativeBinomial::new(r, p)
                    .map_err(|e| DistError::Domain(e.to_string()))?;
                discrete_quantile(|k| distribution.cdf(k), u, u64::MAX)
            },
            | Self::Geometric { p } => {
                if crate::kernels::exact_eq(p, 1.0) || u == 0.0 {
                    return Ok(0.0);
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(f64::INFINITY);
                }
                // This API counts failures (support 0, 1, ...); statrs Geometric
                // counts trials, starting at 1. Use stable logarithms directly.
                let mut q = ((-u).ln_1p() / (-p).ln_1p()).ceil().max(1.0) - 1.0;
                if q < 9_007_199_254_740_992.0 {
                    if self.cdf(q) < u {
                        q += 1.0;
                    }
                    if q > 0.0 && self.cdf(q - 1.0) >= u {
                        q -= 1.0;
                    }
                }
                Ok(q)
            },
            | Self::Weibull { shape, scale } => {
                if shape <= 0.0 || scale <= 0.0 {
                    return Err(DistError::Domain(
                        "Weibull requires positive shape and scale".to_string(),
                    ));
                }
                if u == 0.0 {
                    return Ok(0.0);
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(f64::INFINITY);
                }
                Ok(scale * (-(1.0 - u).ln()).powf(1.0 / shape))
            },
            | Self::Cauchy { loc, scale } => {
                if scale <= 0.0 {
                    return Err(DistError::Domain(
                        "Cauchy requires positive scale".to_string(),
                    ));
                }
                if u == 0.0 {
                    return Ok(f64::NEG_INFINITY);
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(f64::INFINITY);
                }
                Ok(loc + scale * (PI * (u - 0.5)).tan())
            },
            | _ => {
                if u == 0.0 {
                    return Ok(self.support_low());
                }
                if crate::kernels::exact_eq(u, 1.0) {
                    return Ok(self.support_high());
                }
                let mut low = self.support_low();
                let mut high = self.support_high();
                if !low.is_finite() {
                    low = -1.0;
                    for _ in 0..1024 {
                        if self.cdf(low) <= u {
                            break;
                        }
                        low = (low * 2.0).max(-f64::MAX);
                    }
                }
                if !high.is_finite() {
                    high = 1.0;
                    for _ in 0..1024 {
                        if self.cdf(high) >= u {
                            break;
                        }
                        high = (high * 2.0).min(f64::MAX);
                    }
                }
                if self.cdf(high) < u {
                    return Ok(f64::INFINITY);
                }
                if self.cdf(low) > u {
                    return Ok(f64::NEG_INFINITY);
                }
                for _ in 0..1100 {
                    let middle = low * 0.5 + high * 0.5;
                    if crate::kernels::exact_eq(middle, low)
                        || crate::kernels::exact_eq(middle, high)
                    {
                        break;
                    }
                    let probability = self.cdf(middle);
                    if !probability.is_finite() {
                        return Err(DistError::Domain(
                            "distribution CDF did not produce a finite probability".to_string(),
                        ));
                    }
                    if probability < u {
                        low = middle;
                    } else {
                        high = middle;
                    }
                }
                Ok(low * 0.5 + high * 0.5)
            },
        }
    }

    /// The lower end of the support, or `-inf`.
    const fn support_low(&self) -> f64 {
        match *self {
            | Self::Uniform { lo, .. } => lo,
            | Self::Exponential { .. }
            | Self::Gamma { .. }
            | Self::ChiSquared { .. }
            | Self::LogNormal { .. }
            | Self::Weibull { .. }
            | Self::F { .. }
            | Self::Beta { .. } => 0.0,
            | _ => f64::NEG_INFINITY,
        }
    }

    /// The upper end of the support, or `+inf`.
    const fn support_high(&self) -> f64 {
        match *self {
            | Self::Uniform { hi, .. } => hi,
            | Self::Beta { .. } => 1.0,
            | _ => f64::INFINITY,
        }
    }

    /// One sample, from a uniform source on `[0, 1)`.
    ///
    /// # Errors
    ///
    /// [`DistError::Domain`] for a parameter the sampler cannot use.
    pub fn sample(
        &self,
        uniform: &mut dyn FnMut() -> f64,
    ) -> Result<f64, DistError> {
        self.validate()?;
        let one = |u: f64| {
            if u <= 0.0 {
                f64::MIN_POSITIVE
            } else {
                u
            }
        };
        Ok(match *self {
            | Self::Normal { mu, sigma } => {
                // Box-Muller; the second uniform is discarded.
                let u1 = one(uniform());
                let u2 = uniform();
                mu + sigma * (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
            },
            | Self::Uniform { lo, hi } => lo + (hi - lo) * uniform(),
            | Self::Exponential { rate } => -one(uniform()).ln() / rate,
            | Self::Gamma { shape, rate } => gamma_sample(shape, uniform)? / rate,
            | Self::Beta { alpha, beta } => {
                let x = gamma_sample(alpha, uniform)?;
                let y = gamma_sample(beta, uniform)?;
                if x + y == 0.0 {
                    0.5
                } else {
                    x / (x + y)
                }
            },
            | Self::ChiSquared { df } => 2.0 * gamma_sample(df / 2.0, uniform)?,
            | Self::StudentT { df } => {
                let z = Self::Normal { mu: 0.0, sigma: 1.0 }.sample(uniform)?;
                let v = Self::ChiSquared { df }.sample(uniform)?;
                z / (v / df).sqrt()
            },
            | Self::Poisson { lambda } => {
                if lambda == 0.0 {
                    0.0
                } else {
                    self.quantile(uniform())?
                }
            },
            | Self::Binomial { n, p } => {
                let mut successes = 0.0;
                for _ in 0..n {
                    if uniform() < p {
                        successes += 1.0;
                    }
                }
                successes
            },
            | Self::LogNormal { mu, sigma } => (Self::Normal { mu, sigma }.sample(uniform)?).exp(),
            | Self::Weibull { shape, scale } => {
                if shape <= 0.0 || scale <= 0.0 {
                    return Err(DistError::Domain(
                        "Weibull requires positive shape and scale".to_string(),
                    ));
                }
                scale * (-one(uniform()).ln()).powf(1.0 / shape)
            },
            | Self::Cauchy { loc, scale } => {
                if scale <= 0.0 {
                    return Err(DistError::Domain(
                        "Cauchy requires positive scale".to_string(),
                    ));
                }
                loc + scale * (PI * (uniform() - 0.5)).tan()
            },
            | Self::F { df1, df2 } => {
                if df1 <= 0.0 || df2 <= 0.0 {
                    return Err(DistError::Domain(
                        "F requires positive degrees of freedom".to_string(),
                    ));
                }
                let u1 = Self::ChiSquared { df: df1 }.sample(uniform)?;
                let u2 = Self::ChiSquared { df: df2 }.sample(uniform)?;
                (u1 / df1) / (u2 / df2)
            },
            | Self::NegativeBinomial { r, p } => {
                if r <= 0.0 || p <= 0.0 || p > 1.0 {
                    return Err(DistError::Domain(
                        "NegativeBinomial requires `r > 0` and `0 < p <= 1`".to_string(),
                    ));
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    0.0
                } else {
                    let lambda = Self::Gamma {
                        shape: r,
                        rate: p / (1.0 - p),
                    }
                    .sample(uniform)?;
                    Self::Poisson { lambda }.sample(uniform)?
                }
            },
            | Self::Geometric { p } => {
                if p <= 0.0 || p > 1.0 {
                    return Err(DistError::Domain(
                        "Geometric requires `0 < p <= 1`".to_string(),
                    ));
                }
                if crate::kernels::exact_eq(p, 1.0) {
                    0.0
                } else {
                    (one(uniform()).ln() / (-p).ln_1p()).floor()
                }
            },
        })
    }
}

/// Marsaglia–Tsang for `Gamma(shape, 1)`.
fn gamma_sample(
    shape: f64,
    uniform: &mut dyn FnMut() -> f64,
) -> Result<f64, DistError> {
    if shape <= 0.0 {
        return Err(DistError::Domain(
            "a Gamma shape must be positive".to_string(),
        ));
    }
    if shape < 1.0 {
        // Boost by one and scale by U^(1/shape).
        let boosted = gamma_sample(shape + 1.0, uniform)?;
        let u = uniform().max(f64::MIN_POSITIVE);
        return Ok(boosted * u.powf(1.0 / shape));
    }
    let d = shape - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let mut x;
        let mut v;
        loop {
            let normal = Distribution::Normal { mu: 0.0, sigma: 1.0 }.sample(uniform)?;
            x = normal;
            v = 1.0 + c * x;
            if v > 0.0 {
                break;
            }
        }
        v = v * v * v;
        let u = uniform().max(f64::MIN_POSITIVE);
        if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
            return Ok(d * v);
        }
    }
}

/// Multivariate probability distributions.
#[derive(Clone, Debug, PartialEq)]
pub enum MultivariateDistribution {
    /// `MultivariateNormal = { mu, cov }` where `cov` is a row-major `d × d`
    /// symmetric positive-definite covariance matrix.
    MultivariateNormal {
        /// The `d`-dimensional mean vector.
        mu: Vec<f64>,
        /// The `d × d` row-major covariance matrix.
        cov: Vec<f64>,
    },
    /// `Dirichlet = { alpha }` on the `(d - 1)`-simplex.
    Dirichlet {
        /// The concentration parameters `alpha_i > 0`.
        alpha: Vec<f64>,
    },
}

impl MultivariateDistribution {
    /// The dimension `d` of the distribution.
    #[must_use]
    pub const fn dim(&self) -> usize {
        match self {
            | Self::MultivariateNormal { mu, .. } => mu.len(),
            | Self::Dirichlet { alpha } => alpha.len(),
        }
    }

    /// The natural logarithm of the probability density at `x`.
    ///
    /// # Errors
    ///
    /// [`DistError::Empty`] if the distribution has dimension 0,
    /// [`DistError::ShapeMismatch`] if `x` or `cov` has the wrong length,
    /// [`DistError::Domain`] if `Dirichlet` parameters are non-positive, or
    /// [`DistError::NotPositiveDefinite`] if `cov` is not positive definite.
    pub fn log_pdf(
        &self,
        x: &[f64],
    ) -> Result<f64, DistError> {
        if x.iter().any(|v| !v.is_finite()) {
            return Err(DistError::Domain(
                "density coordinates must be finite".to_string(),
            ));
        }
        match self {
            | Self::MultivariateNormal { mu, cov } => {
                if mu.iter().any(|v| !v.is_finite()) {
                    return Err(DistError::Domain(
                        "mean coordinates must be finite".to_string(),
                    ));
                }
                let d = mu.len();
                if d == 0 {
                    return Err(DistError::Empty);
                }
                if x.len() != d {
                    return Err(DistError::ShapeMismatch(format!(
                        "expected vector of length {d}, found {}",
                        x.len()
                    )));
                }
                if cov.len() != d * d {
                    return Err(DistError::ShapeMismatch(format!(
                        "expected {d}x{d} covariance matrix, found {} entries",
                        cov.len()
                    )));
                }
                let l = cholesky(d, d, cov)?;
                // Forward substitution: L * y = x - mu
                let mut y = vec![0.0; d];
                for i in 0..d {
                    let mut sum = x[i] - mu[i];
                    for j in 0..i {
                        sum -= l[i * d + j] * y[j];
                    }
                    let diag = l[i * d + i];
                    if diag <= 0.0 || !diag.is_finite() {
                        return Err(DistError::NotPositiveDefinite);
                    }
                    y[i] = sum / diag;
                }
                let log_det: f64 = 2.0 * (0..d).map(|i| l[i * d + i].ln()).sum::<f64>();
                let quad: f64 = y.iter().map(|v| v * v).sum();
                Ok(-0.5 * (d as f64 * (2.0 * PI).ln() + log_det + quad))
            },
            | Self::Dirichlet { alpha } => {
                let d = alpha.len();
                if d == 0 {
                    return Err(DistError::Empty);
                }
                if alpha.iter().any(|&a| a <= 0.0 || !a.is_finite()) {
                    return Err(DistError::Domain(
                        "Dirichlet concentration parameters must be positive".to_string(),
                    ));
                }
                if x.len() != d {
                    return Err(DistError::ShapeMismatch(format!(
                        "expected vector of length {d}, found {}",
                        x.len()
                    )));
                }
                let sum_x: f64 = x.iter().sum();
                if x.iter().any(|&xi| xi < 0.0) || (sum_x - 1.0).abs() > 1e-8 {
                    return Ok(f64::NEG_INFINITY);
                }
                let alpha_sum: f64 = alpha.iter().sum();
                let log_b: f64 =
                    alpha.iter().map(|&a| ln_gamma(a)).sum::<f64>() - ln_gamma(alpha_sum);
                let mut log_kernel = 0.0;
                for (&a, &xi) in alpha.iter().zip(x.iter()) {
                    if !crate::kernels::exact_eq(a, 1.0) {
                        if xi <= 0.0 {
                            return Ok(if a > 1.0 {
                                f64::NEG_INFINITY
                            } else {
                                f64::INFINITY
                            });
                        }
                        log_kernel += (a - 1.0) * xi.ln();
                    }
                }
                Ok(log_kernel - log_b)
            },
        }
    }

    /// The probability density at `x`.
    ///
    /// # Errors
    ///
    /// Propagates validation errors from [`MultivariateDistribution::log_pdf`].
    pub fn pdf(
        &self,
        x: &[f64],
    ) -> Result<f64, DistError> {
        Ok(self.log_pdf(x)?.exp())
    }

    /// The mean vector.
    ///
    /// # Errors
    ///
    /// [`DistError::Empty`] or [`DistError::Domain`] on invalid parameters.
    pub fn mean(&self) -> Result<Vec<f64>, DistError> {
        match self {
            | Self::MultivariateNormal { mu, .. } => {
                if mu.is_empty() {
                    return Err(DistError::Empty);
                }
                Ok(mu.clone())
            },
            | Self::Dirichlet { alpha } => {
                if alpha.is_empty() {
                    return Err(DistError::Empty);
                }
                if alpha.iter().any(|&a| a <= 0.0 || !a.is_finite()) {
                    return Err(DistError::Domain(
                        "Dirichlet concentration parameters must be positive".to_string(),
                    ));
                }
                let sum: f64 = alpha.iter().sum();
                Ok(alpha.iter().map(|&a| a / sum).collect())
            },
        }
    }

    /// The row-major `d × d` covariance matrix.
    ///
    /// # Errors
    ///
    /// [`DistError::Empty`], [`DistError::ShapeMismatch`], or [`DistError::Domain`]
    /// on invalid parameters.
    pub fn cov(&self) -> Result<Vec<f64>, DistError> {
        match self {
            | Self::MultivariateNormal { mu, cov } => {
                let d = mu.len();
                if d == 0 {
                    return Err(DistError::Empty);
                }
                if cov.len() != d * d {
                    return Err(DistError::ShapeMismatch(format!(
                        "expected {d}x{d} covariance matrix, found {} entries",
                        cov.len()
                    )));
                }
                Ok(cov.clone())
            },
            | Self::Dirichlet { alpha } => {
                let d = alpha.len();
                if d == 0 {
                    return Err(DistError::Empty);
                }
                if alpha.iter().any(|&a| a <= 0.0 || !a.is_finite()) {
                    return Err(DistError::Domain(
                        "Dirichlet concentration parameters must be positive".to_string(),
                    ));
                }
                let a0: f64 = alpha.iter().sum();
                let denom = a0 * a0 * (a0 + 1.0);
                let mut out = vec![0.0; d * d];
                for i in 0..d {
                    for j in 0..d {
                        out[i * d + j] = if i == j {
                            alpha[i] * (a0 - alpha[i]) / denom
                        } else {
                            -alpha[i] * alpha[j] / denom
                        };
                    }
                }
                Ok(out)
            },
        }
    }

    /// Draw one `d`-dimensional sample using a uniform source on `[0, 1)`.
    ///
    /// # Errors
    ///
    /// Propagates domain, shape, or positive-definiteness errors.
    pub fn sample(
        &self,
        uniform: &mut dyn FnMut() -> f64,
    ) -> Result<Vec<f64>, DistError> {
        match self {
            | Self::MultivariateNormal { mu, cov } => {
                if mu.iter().any(|v| !v.is_finite()) {
                    return Err(DistError::Domain(
                        "mean coordinates must be finite".to_string(),
                    ));
                }
                let d = mu.len();
                if d == 0 {
                    return Err(DistError::Empty);
                }
                if cov.len() != d * d {
                    return Err(DistError::ShapeMismatch(format!(
                        "expected {d}x{d} covariance matrix, found {} entries",
                        cov.len()
                    )));
                }
                let l = cholesky(d, d, cov)?;
                let std_normal = Distribution::Normal { mu: 0.0, sigma: 1.0 };
                let mut z = Vec::with_capacity(d);
                for _ in 0..d {
                    z.push(std_normal.sample(uniform)?);
                }
                let mut out = mu.clone();
                for i in 0..d {
                    let mut sum = 0.0;
                    for j in 0..=i {
                        sum += l[i * d + j] * z[j];
                    }
                    out[i] += sum;
                }
                Ok(out)
            },
            | Self::Dirichlet { alpha } => {
                if alpha.is_empty() {
                    return Err(DistError::Empty);
                }
                if alpha.iter().any(|&a| a <= 0.0 || !a.is_finite()) {
                    return Err(DistError::Domain(
                        "Dirichlet concentration parameters must be positive".to_string(),
                    ));
                }
                // Normalize in log space. Scaling each log by the smallest
                // concentration also avoids overflow of ln(U)/alpha for tiny alpha.
                let log_scale = alpha.iter().copied().fold(1.0, f64::min);
                let mut logs = Vec::with_capacity(alpha.len());
                for &a in alpha {
                    let log = if a < 1.0 {
                        gamma_sample(a + 1.0, uniform)?.ln() * log_scale
                            + uniform().max(f64::MIN_POSITIVE).ln() * (log_scale / a)
                    } else {
                        gamma_sample(a, uniform)?.ln() * log_scale
                    };
                    logs.push(log);
                }
                let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let mut weights: Vec<f64> = logs
                    .iter()
                    .map(|v| ((v - max_log) / log_scale).exp())
                    .collect();
                let total: f64 = weights.iter().sum();
                if !total.is_finite() || total <= 0.0 {
                    return Err(DistError::Domain(
                        "Dirichlet sample exceeded floating-point range".to_string(),
                    ));
                }
                for weight in &mut weights {
                    *weight /= total;
                }
                Ok(weights)
            },
        }
    }
}

/// Dimensionless multivariate normal log density with a row-major covariance matrix.
///
/// Rejects non-finite data, mismatched dimensions and non-positive-definite covariance.
pub fn multivariate_normal_log_pdf(
    mu: &[f64],
    covariance: &Mat,
    x: &[f64],
) -> Result<f64, DistError> {
    if covariance.rows != mu.len() || covariance.cols != mu.len() {
        return Err(DistError::ShapeMismatch(
            "covariance dimensions must match mean length".to_string(),
        ));
    }
    MultivariateDistribution::MultivariateNormal {
        mu: mu.to_vec(),
        cov: covariance.data.clone(),
    }
    .log_pdf(x)
}

/// Dimensionless multivariate normal density; covariance dimensions must match the mean.
pub fn multivariate_normal_pdf(
    mu: &[f64],
    covariance: &Mat,
    x: &[f64],
) -> Result<f64, DistError> {
    Ok(multivariate_normal_log_pdf(mu, covariance, x)?.exp())
}

/// One dimensionless multivariate normal sample from a deterministic seed.
///
/// Rejects mismatched dimensions and invalid covariance matrices.
pub fn multivariate_normal_sample(
    mu: &[f64],
    covariance: &Mat,
    seed: u64,
) -> Result<Vec<f64>, DistError> {
    if covariance.rows != mu.len() || covariance.cols != mu.len() {
        return Err(DistError::ShapeMismatch(
            "covariance dimensions must match mean length".to_string(),
        ));
    }
    let mut rng = crate::kernels::random::Rng::new(seed);
    MultivariateDistribution::MultivariateNormal {
        mu: mu.to_vec(),
        cov: covariance.data.clone(),
    }
    .sample(&mut || rng.uniform())
}

/// Dirichlet log density on the unit simplex, with positive finite concentrations.
pub fn dirichlet_log_pdf(
    alpha: &[f64],
    x: &[f64],
) -> Result<f64, DistError> {
    MultivariateDistribution::Dirichlet {
        alpha: alpha.to_vec(),
    }
    .log_pdf(x)
}

/// Dirichlet density on the unit simplex, with positive finite concentrations.
pub fn dirichlet_pdf(
    alpha: &[f64],
    x: &[f64],
) -> Result<f64, DistError> {
    Ok(dirichlet_log_pdf(alpha, x)?.exp())
}

/// One Dirichlet sample from a deterministic seed and positive finite concentrations.
pub fn dirichlet_sample(
    alpha: &[f64],
    seed: u64,
) -> Result<Vec<f64>, DistError> {
    let mut rng = crate::kernels::random::Rng::new(seed);
    MultivariateDistribution::Dirichlet {
        alpha: alpha.to_vec(),
    }
    .sample(&mut || rng.uniform())
}
