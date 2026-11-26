//! The `statrs` distribution family behind one enum.
//!
//! [`crate::kernels::distributions::Distribution`] holds the hand-written
//! kernels the rule layer evaluates. This module is the wider surface: every
//! univariate distribution `statrs` provides, addressed by name and a parameter
//! vector, with a common set of queries.
//!
//! Each query returns `Result` because a distribution whose parameters are
//! invalid, or a point outside the support, has no answer to give.

use statrs::distribution as sd;
use statrs::distribution::Continuous;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::Discrete;
use statrs::distribution::DiscreteCDF;
use statrs::statistics::Distribution as StatsDistribution;

/// One distribution from the `statrs` family.
#[derive(Debug, Clone, PartialEq)]
pub enum StatrsDistribution {
    /// Normal `N(mu, sigma^2)`.
    Normal {
        /// Mean.
        mu: f64,
        /// Standard deviation.
        sigma: f64,
    },
    /// Uniform on `[lo, hi]`.
    Uniform {
        /// Lower bound.
        lo: f64,
        /// Upper bound.
        hi: f64,
    },
    /// Beta on `[0, 1]`.
    Beta {
        /// First shape.
        alpha: f64,
        /// Second shape.
        beta: f64,
    },
    /// Cauchy with location and scale.
    Cauchy {
        /// Location.
        x0: f64,
        /// Scale.
        gamma: f64,
    },
    /// Chi distribution with `k` degrees of freedom.
    Chi {
        /// Degrees of freedom.
        k: f64,
    },
    /// Chi-squared with `k` degrees of freedom.
    ChiSquared {
        /// Degrees of freedom.
        k: f64,
    },
    /// Erlang with integer shape and rate.
    Erlang {
        /// Shape.
        shape: f64,
        /// Rate.
        rate: f64,
    },
    /// Exponential with rate.
    Exponential {
        /// Rate.
        rate: f64,
    },
    /// Fisher–Snedecor `F(d1, d2)`.
    FisherSnedecor {
        /// Numerator degrees of freedom.
        d1: f64,
        /// Denominator degrees of freedom.
        d2: f64,
    },
    /// Gamma with shape and rate.
    Gamma {
        /// Shape.
        shape: f64,
        /// Rate.
        rate: f64,
    },
    /// Gumbel with location and scale.
    Gumbel {
        /// Location.
        location: f64,
        /// Scale.
        scale: f64,
    },
    /// Inverse gamma with shape and rate.
    InverseGamma {
        /// Shape.
        shape: f64,
        /// Rate.
        rate: f64,
    },
    /// Laplace with location and scale.
    Laplace {
        /// Location.
        mu: f64,
        /// Scale.
        b: f64,
    },
    /// Lévy with location and scale.
    Levy {
        /// Location.
        mu: f64,
        /// Scale.
        c: f64,
    },
    /// Log-normal with log-mean and log-standard deviation.
    LogNormal {
        /// Log-mean.
        mu: f64,
        /// Log-standard deviation.
        sigma: f64,
    },
    /// Pareto with scale and shape.
    Pareto {
        /// Scale.
        scale: f64,
        /// Shape.
        shape: f64,
    },
    /// Student's t with `df` degrees of freedom.
    StudentsT {
        /// Degrees of freedom.
        df: f64,
    },
    /// Triangular on `[min, max]` with a mode.
    Triangular {
        /// Lower bound.
        min: f64,
        /// Upper bound.
        max: f64,
        /// Mode.
        mode: f64,
    },
    /// Weibull with shape and scale.
    Weibull {
        /// Shape.
        shape: f64,
        /// Scale.
        scale: f64,
    },
    /// Bernoulli with success probability.
    Bernoulli {
        /// Success probability.
        p: f64,
    },
    /// Binomial with `n` trials and success probability.
    Binomial {
        /// Number of trials.
        n: u64,
        /// Success probability.
        p: f64,
    },
    /// Categorical over a probability vector.
    Categorical {
        /// Probability mass, indexed by outcome.
        weights: Vec<f64>,
    },
    /// Uniform over the integers `lo..=hi`.
    DiscreteUniform {
        /// Lower bound.
        lo: i64,
        /// Upper bound.
        hi: i64,
    },
    /// Geometric with success probability.
    Geometric {
        /// Success probability.
        p: f64,
    },
    /// Hypergeometric draw without replacement.
    Hypergeometric {
        /// Population size.
        population: u64,
        /// Successes in the population.
        successes: u64,
        /// Draws.
        draws: u64,
    },
    /// Negative binomial with `r` failures and success probability.
    NegativeBinomial {
        /// Target number of failures.
        r: f64,
        /// Success probability.
        p: f64,
    },
    /// Poisson with rate.
    Poisson {
        /// Rate.
        lambda: f64,
    },
}

/// The integer point a discrete query lands on.
///
/// A point below zero has no mass, and a point beyond `u64` is clamped to the
/// largest representable index.
fn discrete_point(x: f64) -> Result<u64, String> {
    if x.is_nan() {
        return Err("the query point is not a number".to_string());
    }
    if x < 0.0 {
        return Ok(0);
    }
    let floor = x.floor();
    Ok(if floor >= 18_446_744_073_709_551_615.0 {
        u64::MAX
    } else {
        floor as u64
    })
}

/// The signed integer point a discrete query lands on.
fn signed_discrete_point(x: f64) -> Result<i64, String> {
    if x.is_nan() {
        return Err("the query point is not a number".to_string());
    }
    let floor = x.floor();
    if floor >= i64::MAX as f64 {
        Ok(i64::MAX)
    } else if floor <= i64::MIN as f64 {
        Ok(i64::MIN)
    } else {
        Ok(floor as i64)
    }
}

fn error(what: &str) -> String {
    format!("invalid {what} distribution parameters")
}

/// The continuous family, dispatched with one macro so every query shares the
/// same constructor list.
macro_rules! continuous {
    ($value:expr, |$d:ident| $body:expr) => {
        match $value {
            | StatrsDistribution::Normal { mu, sigma } => {
                let $d = sd::Normal::new(*mu, *sigma).map_err(|_| error("normal"))?;
                $body
            },
            | StatrsDistribution::Uniform { lo, hi } => {
                let $d = sd::Uniform::new(*lo, *hi).map_err(|_| error("uniform"))?;
                $body
            },
            | StatrsDistribution::Beta { alpha, beta } => {
                let $d = sd::Beta::new(*alpha, *beta).map_err(|_| error("beta"))?;
                $body
            },
            | StatrsDistribution::Cauchy { x0, gamma } => {
                let $d = sd::Cauchy::new(*x0, *gamma).map_err(|_| error("cauchy"))?;
                $body
            },
            | StatrsDistribution::Chi { k } => {
                if *k < 1.0 || k.fract() != 0.0 {
                    return Err(error("chi"));
                }
                let $d = sd::Chi::new(*k as u64).map_err(|_| error("chi"))?;
                $body
            },
            | StatrsDistribution::ChiSquared { k } => {
                let $d = sd::ChiSquared::new(*k).map_err(|_| error("chi-squared"))?;
                $body
            },
            | StatrsDistribution::Erlang { shape, rate } => {
                if *shape < 0.0 || shape.fract() != 0.0 {
                    return Err(error("erlang"));
                }
                let $d = sd::Erlang::new(*shape as u64, *rate).map_err(|_| error("erlang"))?;
                $body
            },
            | StatrsDistribution::Exponential { rate } => {
                let $d = sd::Exp::new(*rate).map_err(|_| error("exponential"))?;
                $body
            },
            | StatrsDistribution::FisherSnedecor { d1, d2 } => {
                let $d = sd::FisherSnedecor::new(*d1, *d2).map_err(|_| error("fisher-snedecor"))?;
                $body
            },
            | StatrsDistribution::Gamma { shape, rate } => {
                let $d = sd::Gamma::new(*shape, *rate).map_err(|_| error("gamma"))?;
                $body
            },
            | StatrsDistribution::Gumbel { location, scale } => {
                let $d = sd::Gumbel::new(*location, *scale).map_err(|_| error("gumbel"))?;
                $body
            },
            | StatrsDistribution::InverseGamma { shape, rate } => {
                let $d =
                    sd::InverseGamma::new(*shape, *rate).map_err(|_| error("inverse-gamma"))?;
                $body
            },
            | StatrsDistribution::Laplace { mu, b } => {
                let $d = sd::Laplace::new(*mu, *b).map_err(|_| error("laplace"))?;
                $body
            },
            | StatrsDistribution::Levy { mu, c } => {
                let $d = sd::Levy::new(*mu, *c).map_err(|_| error("levy"))?;
                $body
            },
            | StatrsDistribution::LogNormal { mu, sigma } => {
                let $d = sd::LogNormal::new(*mu, *sigma).map_err(|_| error("log-normal"))?;
                $body
            },
            | StatrsDistribution::Pareto { scale, shape } => {
                let $d = sd::Pareto::new(*scale, *shape).map_err(|_| error("pareto"))?;
                $body
            },
            | StatrsDistribution::StudentsT { df } => {
                let $d = sd::StudentsT::new(0.0, 1.0, *df).map_err(|_| error("student's t"))?;
                $body
            },
            | StatrsDistribution::Triangular { min, max, mode } => {
                let $d = sd::Triangular::new(*min, *max, *mode).map_err(|_| error("triangular"))?;
                $body
            },
            | StatrsDistribution::Weibull { shape, scale } => {
                let $d = sd::Weibull::new(*shape, *scale).map_err(|_| error("weibull"))?;
                $body
            },
            | other => Err(format!("{other:?} is discrete")),
        }
    };
}

/// The discrete family, dispatched the same way.
macro_rules! discrete {
    ($value:expr, |$d:ident| $body:expr) => {
        match $value {
            | StatrsDistribution::Bernoulli { p } => {
                let $d = sd::Bernoulli::new(*p).map_err(|_| error("bernoulli"))?;
                $body
            },
            | StatrsDistribution::Binomial { n, p } => {
                let $d = sd::Binomial::new(*p, *n).map_err(|_| error("binomial"))?;
                $body
            },
            | StatrsDistribution::Categorical { weights } => {
                let $d = sd::Categorical::new(weights).map_err(|_| error("categorical"))?;
                $body
            },
            | StatrsDistribution::Geometric { p } => {
                let $d = sd::Geometric::new(*p).map_err(|_| error("geometric"))?;
                $body
            },
            | StatrsDistribution::Hypergeometric {
                population,
                successes,
                draws,
            } => {
                let $d = sd::Hypergeometric::new(*population, *successes, *draws)
                    .map_err(|_| error("hypergeometric"))?;
                $body
            },
            | StatrsDistribution::NegativeBinomial { r, p } => {
                let $d =
                    sd::NegativeBinomial::new(*r, *p).map_err(|_| error("negative-binomial"))?;
                $body
            },
            | StatrsDistribution::Poisson { lambda } => {
                let $d = sd::Poisson::new(*lambda).map_err(|_| error("poisson"))?;
                $body
            },
            | other => Err(format!("{other:?} is continuous")),
        }
    };
}

/// The one discrete family whose index type is signed.
macro_rules! discrete_signed {
    ($value:expr, |$d:ident| $body:expr) => {
        match $value {
            | StatrsDistribution::DiscreteUniform { lo, hi } => {
                let $d =
                    sd::DiscreteUniform::new(*lo, *hi).map_err(|_| error("discrete-uniform"))?;
                $body
            },
            | other => Err(format!("{other:?} does not use a signed index")),
        }
    };
}

/// The arms whose distribution implements `statrs::statistics::Distribution`.
///
/// `NegativeBinomial` is the one exception: `statrs` gives it mass and CDF but
/// no mean, variance or entropy, so those report an error.
macro_rules! moments {
    ($value:expr, |$d:ident| $body:expr) => {
        match $value {
            | StatrsDistribution::Normal { mu, sigma } => {
                let $d = sd::Normal::new(*mu, *sigma).map_err(|_| error("normal"))?;
                $body
            },
            | StatrsDistribution::Uniform { lo, hi } => {
                let $d = sd::Uniform::new(*lo, *hi).map_err(|_| error("uniform"))?;
                $body
            },
            | StatrsDistribution::Beta { alpha, beta } => {
                let $d = sd::Beta::new(*alpha, *beta).map_err(|_| error("beta"))?;
                $body
            },
            | StatrsDistribution::Cauchy { x0, gamma } => {
                let $d = sd::Cauchy::new(*x0, *gamma).map_err(|_| error("cauchy"))?;
                $body
            },
            | StatrsDistribution::Chi { k } => {
                if *k < 1.0 || k.fract() != 0.0 {
                    return Err(error("chi"));
                }
                let $d = sd::Chi::new(*k as u64).map_err(|_| error("chi"))?;
                $body
            },
            | StatrsDistribution::ChiSquared { k } => {
                let $d = sd::ChiSquared::new(*k).map_err(|_| error("chi-squared"))?;
                $body
            },
            | StatrsDistribution::Erlang { shape, rate } => {
                if *shape < 0.0 || shape.fract() != 0.0 {
                    return Err(error("erlang"));
                }
                let $d = sd::Erlang::new(*shape as u64, *rate).map_err(|_| error("erlang"))?;
                $body
            },
            | StatrsDistribution::Exponential { rate } => {
                let $d = sd::Exp::new(*rate).map_err(|_| error("exponential"))?;
                $body
            },
            | StatrsDistribution::FisherSnedecor { d1, d2 } => {
                let $d = sd::FisherSnedecor::new(*d1, *d2).map_err(|_| error("fisher-snedecor"))?;
                $body
            },
            | StatrsDistribution::Gamma { shape, rate } => {
                let $d = sd::Gamma::new(*shape, *rate).map_err(|_| error("gamma"))?;
                $body
            },
            | StatrsDistribution::Gumbel { location, scale } => {
                let $d = sd::Gumbel::new(*location, *scale).map_err(|_| error("gumbel"))?;
                $body
            },
            | StatrsDistribution::InverseGamma { shape, rate } => {
                let $d =
                    sd::InverseGamma::new(*shape, *rate).map_err(|_| error("inverse-gamma"))?;
                $body
            },
            | StatrsDistribution::Laplace { mu, b } => {
                let $d = sd::Laplace::new(*mu, *b).map_err(|_| error("laplace"))?;
                $body
            },
            | StatrsDistribution::Levy { mu, c } => {
                let $d = sd::Levy::new(*mu, *c).map_err(|_| error("levy"))?;
                $body
            },
            | StatrsDistribution::LogNormal { mu, sigma } => {
                let $d = sd::LogNormal::new(*mu, *sigma).map_err(|_| error("log-normal"))?;
                $body
            },
            | StatrsDistribution::Pareto { scale, shape } => {
                let $d = sd::Pareto::new(*scale, *shape).map_err(|_| error("pareto"))?;
                $body
            },
            | StatrsDistribution::StudentsT { df } => {
                let $d = sd::StudentsT::new(0.0, 1.0, *df).map_err(|_| error("student-t"))?;
                $body
            },
            | StatrsDistribution::Triangular { min, max, mode } => {
                let $d = sd::Triangular::new(*min, *max, *mode).map_err(|_| error("triangular"))?;
                $body
            },
            | StatrsDistribution::Weibull { shape, scale } => {
                let $d = sd::Weibull::new(*shape, *scale).map_err(|_| error("weibull"))?;
                $body
            },
            | StatrsDistribution::Bernoulli { p } => {
                let $d = sd::Bernoulli::new(*p).map_err(|_| error("bernoulli"))?;
                $body
            },
            | StatrsDistribution::Binomial { n, p } => {
                let $d = sd::Binomial::new(*p, *n).map_err(|_| error("binomial"))?;
                $body
            },
            | StatrsDistribution::Categorical { weights } => {
                let $d = sd::Categorical::new(weights).map_err(|_| error("categorical"))?;
                $body
            },
            | StatrsDistribution::DiscreteUniform { lo, hi } => {
                let $d =
                    sd::DiscreteUniform::new(*lo, *hi).map_err(|_| error("discrete-uniform"))?;
                $body
            },
            | StatrsDistribution::Geometric { p } => {
                let $d = sd::Geometric::new(*p).map_err(|_| error("geometric"))?;
                $body
            },
            | StatrsDistribution::Hypergeometric {
                population,
                successes,
                draws,
            } => {
                let $d = sd::Hypergeometric::new(*population, *successes, *draws)
                    .map_err(|_| error("hypergeometric"))?;
                $body
            },
            | StatrsDistribution::Poisson { lambda } => {
                let $d = sd::Poisson::new(*lambda).map_err(|_| error("poisson"))?;
                $body
            },
            | StatrsDistribution::NegativeBinomial { .. } => {
                Err("statrs defines no moments for the negative binomial".to_string())
            },
        }
    };
}

impl StatrsDistribution {
    /// Whether the distribution indexes with a signed integer.
    #[must_use]
    pub const fn is_signed_discrete(&self) -> bool {
        matches!(self, Self::DiscreteUniform { .. })
    }

    /// Whether the distribution is discrete.
    #[must_use]
    pub const fn is_discrete(&self) -> bool {
        matches!(
            self,
            Self::Bernoulli { .. }
                | Self::Binomial { .. }
                | Self::Categorical { .. }
                | Self::DiscreteUniform { .. }
                | Self::Geometric { .. }
                | Self::Hypergeometric { .. }
                | Self::NegativeBinomial { .. }
                | Self::Poisson { .. }
        )
    }

    /// Builds a distribution from a name and its parameters.
    ///
    /// The names are the `statrs` type names in kebab case: `normal`,
    /// `uniform`, `beta`, `cauchy`, `chi`, `chi-squared`, `erlang`,
    /// `exponential`, `fisher-snedecor`, `gamma`, `gumbel`, `inverse-gamma`,
    /// `laplace`, `levy`, `log-normal`, `pareto`, `student-t`, `triangular`,
    /// `weibull`, `bernoulli`, `binomial`, `categorical`,
    /// `discrete-uniform`, `geometric`, `hypergeometric`,
    /// `negative-binomial` and `poisson`.
    ///
    /// # Errors
    /// Returns an error for an unknown name or a wrong parameter count.
    pub fn from_name(
        name: &str,
        params: &[f64],
    ) -> Result<Self, String> {
        let p = |index: usize| -> Result<f64, String> {
            params
                .get(index)
                .copied()
                .ok_or_else(|| format!("{name} needs parameter {index}"))
        };
        let u = |index: usize| -> Result<u64, String> {
            let value = p(index)?;
            if value < 0.0 || value.fract() != 0.0 {
                return Err(format!(
                    "{name} parameter {index} must be a non-negative integer"
                ));
            }
            Ok(value as u64)
        };
        let i = |index: usize| -> Result<i64, String> {
            let value = p(index)?;
            if value.fract() != 0.0 {
                return Err(format!("{name} parameter {index} must be an integer"));
            }
            Ok(value as i64)
        };
        Ok(match name {
            | "normal" => {
                Self::Normal {
                    mu: p(0)?,
                    sigma: p(1)?,
                }
            },
            | "uniform" => Self::Uniform { lo: p(0)?, hi: p(1)? },
            | "beta" => {
                Self::Beta {
                    alpha: p(0)?,
                    beta: p(1)?,
                }
            },
            | "cauchy" => {
                Self::Cauchy {
                    x0: p(0)?,
                    gamma: p(1)?,
                }
            },
            | "chi" => Self::Chi { k: p(0)? },
            | "chi-squared" => Self::ChiSquared { k: p(0)? },
            | "erlang" => {
                Self::Erlang {
                    shape: p(0)?,
                    rate: p(1)?,
                }
            },
            | "exponential" => Self::Exponential { rate: p(0)? },
            | "fisher-snedecor" => Self::FisherSnedecor { d1: p(0)?, d2: p(1)? },
            | "gamma" => {
                Self::Gamma {
                    shape: p(0)?,
                    rate: p(1)?,
                }
            },
            | "gumbel" => {
                Self::Gumbel {
                    location: p(0)?,
                    scale: p(1)?,
                }
            },
            | "inverse-gamma" => {
                Self::InverseGamma {
                    shape: p(0)?,
                    rate: p(1)?,
                }
            },
            | "laplace" => Self::Laplace { mu: p(0)?, b: p(1)? },
            | "levy" => Self::Levy { mu: p(0)?, c: p(1)? },
            | "log-normal" => {
                Self::LogNormal {
                    mu: p(0)?,
                    sigma: p(1)?,
                }
            },
            | "pareto" => {
                Self::Pareto {
                    scale: p(0)?,
                    shape: p(1)?,
                }
            },
            | "student-t" => Self::StudentsT { df: p(0)? },
            | "triangular" => {
                Self::Triangular {
                    min: p(0)?,
                    max: p(1)?,
                    mode: p(2)?,
                }
            },
            | "weibull" => {
                Self::Weibull {
                    shape: p(0)?,
                    scale: p(1)?,
                }
            },
            | "bernoulli" => Self::Bernoulli { p: p(0)? },
            | "binomial" => Self::Binomial { n: u(0)?, p: p(1)? },
            | "categorical" => {
                Self::Categorical {
                    weights: params.to_vec(),
                }
            },
            | "discrete-uniform" => Self::DiscreteUniform { lo: i(0)?, hi: i(1)? },
            | "geometric" => Self::Geometric { p: p(0)? },
            | "hypergeometric" => {
                Self::Hypergeometric {
                    population: u(0)?,
                    successes: u(1)?,
                    draws: u(2)?,
                }
            },
            | "negative-binomial" => Self::NegativeBinomial { r: p(0)?, p: p(1)? },
            | "poisson" => Self::Poisson { lambda: p(0)? },
            | other => return Err(format!("unknown distribution `{other}`")),
        })
    }

    /// The density (continuous) or mass (discrete) at `x`.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid.
    pub fn pdf(
        &self,
        x: f64,
    ) -> Result<f64, String> {
        if self.is_signed_discrete() {
            let point = signed_discrete_point(x)?;
            discrete_signed!(self, |d| Ok(d.pmf(point)))
        } else if self.is_discrete() {
            let point = discrete_point(x)?;
            discrete!(self, |d| Ok(d.pmf(point)))
        } else {
            continuous!(self, |d| Ok(d.pdf(x)))
        }
    }

    /// The cumulative probability at `x`.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid.
    pub fn cdf(
        &self,
        x: f64,
    ) -> Result<f64, String> {
        if self.is_signed_discrete() {
            let point = signed_discrete_point(x)?;
            discrete_signed!(self, |d| Ok(d.cdf(point)))
        } else if self.is_discrete() {
            let point = discrete_point(x)?;
            discrete!(self, |d| Ok(d.cdf(point)))
        } else {
            continuous!(self, |d| Ok(d.cdf(x)))
        }
    }

    /// The survival probability `1 - cdf(x)`, computed without cancellation.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid.
    pub fn sf(
        &self,
        x: f64,
    ) -> Result<f64, String> {
        if self.is_signed_discrete() {
            let point = signed_discrete_point(x)?;
            discrete_signed!(self, |d| Ok(d.sf(point)))
        } else if self.is_discrete() {
            let point = discrete_point(x)?;
            discrete!(self, |d| Ok(d.sf(point)))
        } else {
            continuous!(self, |d| Ok(d.sf(x)))
        }
    }

    /// The mean, when it exists.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid or the mean is
    /// infinite.
    pub fn mean(&self) -> Result<f64, String> {
        let value = moments!(self, |d| Ok(d.mean().map_or(f64::NAN, |v| v)))?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err("the mean does not exist".to_string())
        }
    }

    /// The variance, when it exists.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid or the variance is
    /// infinite.
    pub fn variance(&self) -> Result<f64, String> {
        let value = moments!(self, |d| Ok(d.variance().map_or(f64::NAN, |v| v)))?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err("the variance does not exist".to_string())
        }
    }

    /// The differential (or discrete) entropy, when `statrs` defines it.
    ///
    /// # Errors
    /// Returns an error when the parameters are invalid.
    pub fn entropy(&self) -> Result<f64, String> {
        let value = moments!(self, |d| Ok(d.entropy().map_or(f64::NAN, |v| v)))?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err("the entropy is not defined for this distribution".to_string())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::StatrsDistribution;

    #[test]
    fn normal_matches_closed_forms() {
        let normal = StatrsDistribution::from_name("normal", &[0.0, 1.0]).expect("normal");
        let density = normal.pdf(0.0).expect("pdf");
        assert!((density - 1.0 / (2.0 * std::f64::consts::PI).sqrt()).abs() < 1e-12);
        assert!((normal.cdf(0.0).expect("cdf") - 0.5).abs() < 1e-12);
        assert!((normal.sf(0.0).expect("sf") - 0.5).abs() < 1e-12);
        assert!(normal.mean().expect("mean").abs() < 1e-12);
        assert!((normal.variance().expect("variance") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn discrete_mass_and_cdf_agree() {
        let binomial = StatrsDistribution::from_name("binomial", &[10.0, 0.5]).expect("binomial");
        let mut total = 0.0;
        for k in 0..=10 {
            total += binomial.pdf(f64::from(k)).expect("pmf");
        }
        assert!((total - 1.0).abs() < 1e-12, "total {total}");
        assert!((binomial.cdf(10.0).expect("cdf") - 1.0).abs() < 1e-12);
        assert!((binomial.mean().expect("mean") - 5.0).abs() < 1e-12);
    }

    #[test]
    fn the_wider_family_is_reachable_by_name() {
        for (name, params) in [
            ("cauchy", vec![0.0, 1.0]),
            ("gumbel", vec![0.0, 1.0]),
            ("laplace", vec![0.0, 1.0]),
            ("levy", vec![0.0, 1.0]),
            ("pareto", vec![1.0, 2.0]),
            ("triangular", vec![0.0, 2.0, 1.0]),
            ("weibull", vec![2.0, 1.0]),
            ("chi", vec![3.0]),
            ("erlang", vec![2.0, 1.0]),
            ("inverse-gamma", vec![3.0, 2.0]),
            ("geometric", vec![0.25]),
            ("hypergeometric", vec![50.0, 10.0, 5.0]),
            ("negative-binomial", vec![3.0, 0.4]),
            ("discrete-uniform", vec![1.0, 6.0]),
            ("categorical", vec![0.2, 0.3, 0.5]),
        ] {
            let distribution = StatrsDistribution::from_name(name, &params)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            let point = 1.0;
            let density = distribution.pdf(point).expect("pdf");
            let cdf = distribution.cdf(point).expect("cdf");
            assert!(density.is_finite(), "{name} density {density}");
            assert!((0.0..=1.0).contains(&cdf), "{name} cdf {cdf}");
        }
    }

    #[test]
    fn invalid_parameters_are_reported() {
        assert!(StatrsDistribution::from_name("normal", &[0.0]).is_err());
        assert!(StatrsDistribution::from_name("nope", &[]).is_err());
        let bad = StatrsDistribution::from_name("normal", &[0.0, -1.0]).expect("shape is valid");
        assert!(bad.pdf(0.0).is_err());
    }
}
