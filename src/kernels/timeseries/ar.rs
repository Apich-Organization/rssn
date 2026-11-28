//! Autoregressive models estimated by the Yule–Walker equations.

use crate::kernels::timeseries::acf::mean_of;

/// A fitted autoregressive model with an intercept.
#[derive(Debug, Clone, PartialEq)]
pub struct ArFit {
    /// Autoregressive coefficients `phi_1 … phi_p`, in lag order.
    pub phi: Vec<f64>,
    /// Innovation variance.
    pub sigma2: f64,
    /// Akaike information criterion of the Gaussian likelihood.
    pub aic: f64,
    /// Intercept of the mean-reverting representation,
    /// `mu (1 - Σ phi_j)`.
    pub intercept: f64,
}

/// Fits an autoregression of order `p` by the Yule–Walker equations, solved
/// with the Levinson–Durbin recursion.
///
/// The series is centred on its sample mean before the autocovariances are
/// formed; the returned [`ArFit::intercept`] maps the centred model back to
/// the original scale.
///
/// # Errors
/// Returns a lowercase message when `p` is zero, the series is too short for
/// the requested order, it has no variance, or the recursion meets a
/// non-positive partial variance (a sign that the sample is not stationary).
pub fn ar_yule_walker(
    x: &[f64],
    p: usize,
) -> Result<ArFit, String> {
    if p == 0 {
        return Err("autoregressive order must be positive".to_string());
    }
    let n = x.len();
    if n < p + 2 {
        return Err("series too short for the requested order".to_string());
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err("series contains a non-finite value".to_string());
    }

    let mean = mean_of(x);
    let mut autocovariance = vec![0.0; p + 1];
    for lag in 0..=p {
        let mut sum = 0.0;
        for t in lag..n {
            sum += (x[t] - mean) * (x[t - lag] - mean);
        }
        autocovariance[lag] = sum / n as f64;
    }
    if autocovariance[0] <= 0.0 || !autocovariance[0].is_finite() {
        return Err("series has no variance".to_string());
    }
    let rho: Vec<f64> = autocovariance
        .iter()
        .map(|g| g / autocovariance[0])
        .collect();

    let mut coefficients = vec![0.0; p + 1];
    let mut previous = vec![0.0; p + 1];
    let mut partial_variance = autocovariance[0];
    for order in 1..=p {
        let mut numerator = rho[order];
        let mut denominator = 1.0;
        for j in 1..order {
            numerator -= previous[j] * rho[order - j];
            denominator -= previous[j] * rho[j];
        }
        if denominator <= 1e-12 {
            return Err("series is not stationary at the requested order".to_string());
        }
        let reflection = numerator / denominator;
        coefficients[order] = reflection;
        for j in 1..order {
            coefficients[j] = previous[j] - reflection * previous[order - j];
        }
        partial_variance *= 1.0 - reflection * reflection;
        if partial_variance <= 0.0 {
            return Err("series is not stationary at the requested order".to_string());
        }
        previous.copy_from_slice(&coefficients);
    }

    let sigma2 = partial_variance;
    let sum_phi: f64 = coefficients[1..].iter().sum();
    let intercept = mean * (1.0 - sum_phi);
    let loglik = -0.5 * n as f64 * ((2.0 * std::f64::consts::PI * sigma2).ln() + 1.0);
    let aic = -2.0 * loglik + 2.0 * (p + 2) as f64;

    Ok(ArFit {
        phi: coefficients[1..].to_vec(),
        sigma2,
        aic,
        intercept,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic linear congruential generator mapped to `[-1, 1)`.
    fn lcg(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let unit = (state >> 11) as f64 / (1_u64 << 53) as f64;
            unit * 2.0 - 1.0
        }
    }

    #[test]
    fn ar_one_recovers_its_coefficient() {
        let mut noise = lcg(9_876_543);
        let mut x = vec![0.0; 400];
        for t in 1..400 {
            x[t] = 0.6 * x[t - 1] + noise();
        }
        let fit = ar_yule_walker(&x, 1).expect("ar(1) fit");
        assert!(
            (fit.phi[0] - 0.6).abs() < 0.05,
            "estimated phi {}",
            fit.phi[0]
        );
        assert!(fit.sigma2 > 0.0);
        assert!(fit.aic.is_finite());
    }

    #[test]
    fn intercept_recovers_a_non_zero_mean() {
        let mut noise = lcg(135_791);
        let mut x = vec![0.0; 500];
        // Stationary AR(1) with mean 4: x_t = 0.5 x_{t-1} + 2 + e_t.
        for t in 1..500 {
            x[t] = 0.5 * x[t - 1] + 2.0 + noise();
        }
        let fit = ar_yule_walker(&x, 1).expect("ar(1) fit");
        // The intercept is mu (1 - phi) = 4 * 0.5 = 2 for this process.
        assert!((fit.intercept - 2.0).abs() < 0.15, "{}", fit.intercept);
        let mu = fit.intercept / (1.0 - fit.phi[0]);
        assert!((mu - 4.0).abs() < 0.3, "{mu}");
    }

    #[test]
    fn ar_two_recovers_both_coefficients() {
        let mut noise = lcg(246_802);
        let mut x = vec![0.0; 2000];
        for t in 2..2000 {
            x[t] = 0.4 * x[t - 1] - 0.25 * x[t - 2] + noise();
        }
        let fit = ar_yule_walker(&x, 2).expect("ar(2) fit");
        assert!((fit.phi[0] - 0.4).abs() < 0.05, "{:?}", fit.phi);
        assert!((fit.phi[1] + 0.25).abs() < 0.05, "{:?}", fit.phi);
    }

    #[test]
    fn invalid_orders_are_rejected() {
        assert!(ar_yule_walker(&[1.0, 2.0, 3.0], 0).is_err());
        assert!(ar_yule_walker(&[1.0, 2.0], 5).is_err());
        assert!(ar_yule_walker(&[1.0; 50], 2).is_err());
    }
}
