//! Link functions mapping the linear predictor `eta` to the mean `mu`.
//!
//! Every link provides three pieces of information that the iteratively
//! reweighted least-squares loop needs: the map `g`, its inverse `g⁻¹` and the
//! derivative `g'`. The derivative enters the working weights
//! `w = 1 / (V(mu) g'(mu)²)` and the working response
//! `z = eta + (y - mu) g'(mu)`.

use crate::kernels::special::erfc_numerical;

/// The link function connecting the linear predictor to the mean response.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Link {
    /// Identity link, `eta = mu`; canonical for the Gaussian family.
    Identity,
    /// Log link, `eta = ln(mu)`; canonical for the Poisson family.
    Log,
    /// Logit link, `eta = ln(mu / (1 - mu))`; canonical for the binomial family.
    Logit,
    /// Probit link, `eta = Phi⁻¹(mu)` with `Phi` the standard normal CDF.
    Probit,
    /// Inverse link, `eta = 1 / mu`; canonical for the gamma family.
    Inverse,
    /// Log-log link, `eta = -ln(-ln(mu))`.
    LogLog,
    /// Complementary log-log link, `eta = ln(-ln(1 - mu))`.
    ComplementaryLogLog,
}

impl Link {
    pub(crate) fn named(name: &str) -> Option<Self> {
        match name {
            | "identity" => Some(Self::Identity),
            | "log" => Some(Self::Log),
            | "logit" => Some(Self::Logit),
            | "probit" => Some(Self::Probit),
            | "inverse" => Some(Self::Inverse),
            | "loglog" | "log_log" => Some(Self::LogLog),
            | "cloglog" | "complementary_loglog" => Some(Self::ComplementaryLogLog),
            | _ => None,
        }
    }
}

/// Standard normal probability density at `z`.
#[must_use]
fn normal_density(z: f64) -> f64 {
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// Standard normal cumulative distribution at `z`.
#[must_use]
fn normal_cdf(z: f64) -> f64 {
    0.5 * erfc_numerical(-z / std::f64::consts::SQRT_2)
}

/// Inverse standard normal distribution `Phi^{-1}(p)`.
///
/// Delegates to [`crate::kernels::special::inverse_erfc`] for full double
/// precision. Returns `-inf` for `p <= 0` and `+inf` for `p >= 1`; `NaN` for `NaN`.
#[must_use]
pub(crate) fn normal_quantile(p: f64) -> f64 {
    if p.is_nan() {
        return f64::NAN;
    }
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    -std::f64::consts::SQRT_2 * crate::kernels::special::inverse_erfc(2.0 * p)
}

/// Evaluates the link function `g(mu)`.
///
/// Out-of-range arguments follow the limits of the function: the logit link
/// returns `±inf` at `mu = 0` or `mu = 1`, and the log link returns `-inf`
/// at `mu = 0`.
#[must_use]
pub fn link(
    kind: Link,
    mu: f64,
) -> f64 {
    match kind {
        | Link::Identity => mu,
        | Link::Log => mu.ln(),
        | Link::Logit => (mu / (1.0 - mu)).ln(),
        | Link::Probit => normal_quantile(mu),
        | Link::Inverse => 1.0 / mu,
        | Link::LogLog => -(-mu.ln()).ln(),
        | Link::ComplementaryLogLog => (-(-mu).ln_1p()).ln(),
    }
}

/// Evaluates the inverse link `g⁻¹(eta)`, i.e. the mean implied by a linear
/// predictor.
#[must_use]
pub fn link_inverse(
    kind: Link,
    eta: f64,
) -> f64 {
    match kind {
        | Link::Identity => eta,
        | Link::Log => eta.exp(),
        | Link::Logit if eta >= 0.0 => 1.0 / (1.0 + (-eta).exp()),
        | Link::Logit => {
            let positive = eta.exp();
            positive / (1.0 + positive)
        },
        | Link::Probit => normal_cdf(eta),
        | Link::Inverse => 1.0 / eta,
        | Link::LogLog => (-(-eta).exp()).exp(),
        | Link::ComplementaryLogLog => -(-eta.exp()).exp_m1(),
    }
}

/// Evaluates the derivative `g'(mu) = d eta / d mu`.
///
/// The reciprocal link has a negative derivative. IRLS uses the square of the
/// derivative in its working weights `1 / (V(mu) g'(mu)²)`.
#[must_use]
pub fn link_derivative(
    kind: Link,
    mu: f64,
) -> f64 {
    match kind {
        | Link::Identity => 1.0,
        | Link::Log => 1.0 / mu,
        | Link::Logit => 1.0 / (mu * (1.0 - mu)),
        | Link::Probit => {
            let eta = normal_quantile(mu);
            1.0 / normal_density(eta)
        },
        | Link::Inverse => -1.0 / (mu * mu),
        | Link::LogLog => -1.0 / (mu * mu.ln()),
        | Link::ComplementaryLogLog => -1.0 / ((1.0 - mu) * (-mu).ln_1p()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL: [Link; 7] = [
        Link::Identity,
        Link::Log,
        Link::Logit,
        Link::Probit,
        Link::Inverse,
        Link::LogLog,
        Link::ComplementaryLogLog,
    ];

    /// Points where every link is finite and its derivative is finite.
    fn interior(kind: Link) -> f64 {
        match kind {
            | Link::Logit | Link::Probit | Link::LogLog | Link::ComplementaryLogLog => 0.3,
            | Link::Log | Link::Inverse | Link::Identity => 1.7,
        }
    }

    #[test]
    fn inverse_round_trips() {
        for kind in ALL {
            let mu = interior(kind);
            let back = link_inverse(kind, link(kind, mu));
            assert!((back - mu).abs() < 1e-9, "{kind:?}: {back} vs {mu}");
        }
    }

    #[test]
    fn derivative_matches_finite_difference() {
        for kind in ALL {
            let mu = interior(kind);
            let h = 1e-6 * mu.max(1.0);
            let numeric = (link(kind, mu + h) - link(kind, mu - h)) / (2.0 * h);
            let analytic = link_derivative(kind, mu);
            assert!(
                (numeric - analytic).abs() < 1e-5 * analytic.abs().max(1.0),
                "{kind:?}: numeric {numeric} analytic {analytic}"
            );
        }
    }

    #[test]
    fn probability_links_keep_representable_left_tails() {
        let predictor: f64 = -710.0;
        let tiny = predictor.exp();
        assert!(tiny > 0.0);
        assert_eq!(link_inverse(Link::Logit, predictor), tiny);
        assert_eq!(link_inverse(Link::ComplementaryLogLog, predictor), tiny);
        assert!((link(Link::ComplementaryLogLog, tiny) - predictor).abs() < 1e-12);

        let mean: f64 = 1e-20;
        assert_eq!(link(Link::ComplementaryLogLog, mean), mean.ln());
        let derivative = link_derivative(Link::ComplementaryLogLog, mean);
        assert!((derivative * mean - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normal_quantile_known_values() {
        assert!(normal_quantile(0.5).abs() < 1e-12);
        assert!((normal_quantile(0.975) - 1.959_963_984_540_054).abs() < 1e-9);
        assert!((normal_quantile(0.025) + 1.959_963_984_540_054).abs() < 1e-9);
        assert!((normal_quantile(0.001) + 3.090_232_306_174_118).abs() < 1e-8);
    }

    #[test]
    fn probit_uses_normal_cdf() {
        assert!((link_inverse(Link::Probit, 0.0) - 0.5).abs() < 1e-12);
        assert!((link_inverse(Link::Probit, 1.0) - 0.841_344_746_068_543).abs() < 1e-9);
    }

    #[test]
    fn loglog_and_cloglog_limits() {
        assert!(link_inverse(Link::LogLog, 10.0) > 0.999);
        assert!(link_inverse(Link::ComplementaryLogLog, -10.0) < 1e-4);
    }
}
