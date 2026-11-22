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

/// Inverse standard normal distribution (Acklam's rational approximation
/// followed by one Halley refinement step).
///
/// Returns `-inf` for `p <= 0` and `+inf` for `p >= 1`; `NaN` for `NaN`.
#[must_use]
pub(crate) fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p.is_nan() {
        return f64::NAN;
    }
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }

    let mut z = if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };

    // Halley refinement: z <- z - (Phi(z) - p) / phi(z), corrected for the
    // density derivative so the result is accurate to a few ulps.
    let err = normal_cdf(z) - p;
    let density = normal_density(z);
    if density > 0.0 && err.is_finite() {
        let u = err / density;
        z -= u / (1.0 + 0.5 * z * u);
    }
    z
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
        | Link::ComplementaryLogLog => (-(1.0 - mu).ln()).ln(),
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
        | Link::Logit => 1.0 / (1.0 + (-eta).exp()),
        | Link::Probit => normal_cdf(eta),
        | Link::Inverse => 1.0 / eta,
        | Link::LogLog => (-(-eta).exp()).exp(),
        | Link::ComplementaryLogLog => 1.0 - (-eta.exp()).exp(),
    }
}

/// Evaluates the derivative `g'(mu) = d eta / d mu`.
///
/// The derivative is always positive on the interior of the mean space, so
/// the working weights `1 / (V(mu) g'(mu)²)` stay positive.
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
        | Link::ComplementaryLogLog => -1.0 / ((1.0 - mu) * (1.0 - mu).ln()),
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
