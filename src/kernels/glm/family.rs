//! Response families: variance functions, canonical links and deviance
//! contributions.
//!
//! A family fixes the mean–variance relationship `V(mu)`, the unit deviance
//! `d(y, mu)` and the canonical link. Together with the link these determine
//! the iteratively reweighted least-squares iteration and the model deviance.

/// The distributional family of the response.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Family {
    /// Normal response with constant variance; canonical link identity.
    Gaussian,
    /// Poisson counts; canonical link log.
    Poisson,
    /// Binomial proportions in `[0, 1]`; canonical link logit.
    Binomial,
    /// Positive continuous response with variance proportional to `mu²`;
    /// canonical link inverse.
    Gamma,
    /// Positive continuous response with variance proportional to `mu³`.
    InverseGaussian,
}

/// The variance function `V(mu)` of the family.
#[must_use]
pub fn variance(
    family: Family,
    mu: f64,
) -> f64 {
    match family {
        | Family::Gaussian => 1.0,
        | Family::Poisson => mu,
        | Family::Binomial => mu * (1.0 - mu),
        | Family::Gamma => mu * mu,
        | Family::InverseGaussian => mu * mu * mu,
    }
}

/// Whether the family has a known (unit) dispersion.
#[must_use]
pub(crate) fn has_unit_dispersion(family: Family) -> bool {
    matches!(family, |Family::Poisson| Family::Binomial)
}

/// Unit deviance `d(y, mu)`, i.e. twice the log-likelihood gap between the
/// saturated model and a fit with mean `mu`.
#[must_use]
pub(crate) fn unit_deviance(
    family: Family,
    y: f64,
    mu: f64,
) -> f64 {
    match family {
        | Family::Gaussian => {
            let r = y - mu;
            r * r
        },
        | Family::Poisson => {
            if y <= 0.0 {
                2.0 * mu
            } else {
                2.0 * (y * (y / mu).ln() - (y - mu))
            }
        },
        | Family::Binomial => {
            let first = if y <= 0.0 {
                0.0
            } else {
                y * (y / mu).ln()
            };
            let second = if y >= 1.0 {
                0.0
            } else {
                (1.0 - y) * ((1.0 - y) / (1.0 - mu)).ln()
            };
            2.0 * (first + second)
        },
        | Family::Gamma => 2.0 * ((y - mu) / mu - (y / mu).ln()),
        | Family::InverseGaussian => {
            let r = y - mu;
            r * r / (y * mu * mu)
        },
    }
}

/// Clamps a mean to the interior of the family's support so that the link,
/// the variance function and the deviance stay finite.
#[must_use]
pub(crate) fn clamp_mean(
    family: Family,
    mu: f64,
) -> f64 {
    /// Smallest mean magnitude used to keep logs and divisions finite.
    const FLOOR: f64 = 1e-12;
    match family {
        | Family::Gaussian => mu,
        | Family::Poisson | Family::Gamma | Family::InverseGaussian => mu.max(FLOOR),
        | Family::Binomial => mu.clamp(FLOOR, 1.0 - FLOOR),
    }
}

/// A sensible starting mean for each observation, used to seed the iteration.
#[must_use]
pub(crate) fn starting_mean(
    family: Family,
    y: f64,
) -> f64 {
    match family {
        | Family::Gaussian => y,
        | Family::Poisson => y + 0.1,
        | Family::Binomial => 0.5 * y + 0.25,
        | Family::Gamma | Family::InverseGaussian => y.max(1e-3),
    }
}

/// Validates the response vector against the family's support.
///
/// # Errors
/// Returns a lowercase message naming the offending constraint.
pub(crate) fn validate_response(
    family: Family,
    y: &[f64],
) -> Result<(), String> {
    match family {
        | Family::Gaussian => {
            if y.iter().any(|v| !v.is_finite()) {
                return Err("response contains a non-finite value".to_string());
            }
        },
        | Family::Poisson => {
            if y.iter().any(|v| !v.is_finite() || *v < 0.0) {
                return Err("poisson response must be finite and non-negative".to_string());
            }
        },
        | Family::Binomial => {
            if y.iter().any(|v| !v.is_finite() || *v < 0.0 || *v > 1.0) {
                return Err("binomial response must lie in [0, 1]".to_string());
            }
        },
        | Family::Gamma | Family::InverseGaussian => {
            if y.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err("response must be finite and strictly positive".to_string());
            }
        },
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variance_functions() {
        assert!((variance(Family::Gaussian, 3.0) - 1.0).abs() < 1e-12);
        assert!((variance(Family::Poisson, 3.0) - 3.0).abs() < 1e-12);
        assert!((variance(Family::Binomial, 0.25) - 0.1875).abs() < 1e-12);
        assert!((variance(Family::Gamma, 2.0) - 4.0).abs() < 1e-12);
        assert!((variance(Family::InverseGaussian, 2.0) - 8.0).abs() < 1e-12);
    }

    #[test]
    fn gaussian_deviance_is_squared_error() {
        assert!((unit_deviance(Family::Gaussian, 3.0, 1.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn poisson_deviance_at_saturation_is_zero() {
        for y in [0.0, 1.0, 4.0, 9.5] {
            assert!(unit_deviance(Family::Poisson, y, y).abs() < 1e-12);
        }
        assert!((unit_deviance(Family::Poisson, 0.0, 2.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn binomial_deviance_at_saturation_is_zero() {
        for y in [0.0, 0.2, 0.5, 1.0] {
            assert!(unit_deviance(Family::Binomial, y, y).abs() < 1e-12);
        }
    }

    #[test]
    fn gamma_and_inverse_gaussian_vanish_at_saturation() {
        for y in [0.5, 1.0, 3.0] {
            assert!(unit_deviance(Family::Gamma, y, y).abs() < 1e-12);
            assert!(unit_deviance(Family::InverseGaussian, y, y).abs() < 1e-12);
        }
    }

    #[test]
    fn invalid_responses_are_rejected() {
        assert!(validate_response(Family::Poisson, &[-1.0]).is_err());
        assert!(validate_response(Family::Binomial, &[1.5]).is_err());
        assert!(validate_response(Family::Gamma, &[0.0]).is_err());
        assert!(validate_response(Family::InverseGaussian, &[-2.0]).is_err());
        assert!(validate_response(Family::Gaussian, &[f64::NAN]).is_err());
    }

    #[test]
    fn unit_dispersion_flags() {
        assert!(has_unit_dispersion(Family::Poisson));
        assert!(has_unit_dispersion(Family::Binomial));
        assert!(!has_unit_dispersion(Family::Gaussian));
    }

    #[test]
    fn clamping_stays_inside_the_support() {
        assert!(clamp_mean(Family::Poisson, -1.0) > 0.0);
        let p = clamp_mean(Family::Binomial, 1.5);
        assert!(p < 1.0 && p > 0.0);
        assert!(clamp_mean(Family::Gamma, 0.0) > 0.0);
    }
}
