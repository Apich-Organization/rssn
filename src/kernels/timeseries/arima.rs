//! ARIMA models: Hannan–Rissanen initialisation followed by conditional
//! sum-of-squares refinement.

use crate::kernels::dense::Mat;
use crate::kernels::dense::lstsq;
use crate::kernels::timeseries::acf::difference;
use crate::kernels::timeseries::acf::mean_of;

/// A fitted ARIMA model for the differenced series.
#[derive(Debug, Clone, PartialEq)]
pub struct ArimaFit {
    /// Autoregressive coefficients `phi_1 … phi_p`.
    pub ar: Vec<f64>,
    /// Moving-average coefficients `theta_1 … theta_q`.
    pub ma: Vec<f64>,
    /// Innovation variance of the differenced series.
    pub sigma2: f64,
    /// Akaike information criterion.
    pub aic: f64,
    /// Gaussian log-likelihood at the estimate.
    pub loglik: f64,
    /// Intercept of the differenced series.
    pub mean: f64,
}

/// Largest autoregressive order accepted by [`arima_fit`].
const MAX_AR: usize = 6;
/// Largest moving-average order accepted by [`arima_fit`].
const MAX_MA: usize = 4;

/// One-step-ahead residuals and their sum of squares for a parameter vector
/// `[c, phi_1 … phi_p, theta_1 … theta_q]`.
///
/// Pre-sample values are replaced by the sample mean of the series and
/// pre-sample shocks by zero, the usual conditional-sum-of-squares convention.
fn conditional_residuals(
    w: &[f64],
    parameters: &[f64],
    p: usize,
    q: usize,
) -> (Vec<f64>, f64) {
    let n = w.len();
    let intercept = parameters.first().copied().unwrap_or(0.0);
    let mean = mean_of(w);
    let mut residuals = vec![0.0; n];
    let mut sse = 0.0;
    for t in 0..n {
        let mut prediction = intercept;
        for i in 1..=p {
            let lagged = if t >= i { w[t - i] } else { mean };
            prediction += parameters.get(i).copied().unwrap_or(0.0) * lagged;
        }
        for j in 1..=q {
            let lagged = if t >= j {
                residuals[t - j]
            } else {
                0.0
            };
            prediction += parameters.get(p + j).copied().unwrap_or(0.0) * lagged;
        }
        let error = w[t] - prediction;
        residuals[t] = error;
        sse += error * error;
    }
    (residuals, sse)
}

/// Fits an ARIMA(`p`, `d`, `q`) model.
///
/// The series is differenced `d` times and an ARMA(`p`, `q`) is fitted to the
/// result by the Hannan–Rissanen procedure: a long autoregression provides
/// residuals, a regression on lagged values and lagged residuals gives initial
/// parameters, and a Gauss–Newton search on the conditional sum of squares
/// refines them. Orders are capped at `p ≤ 6` and `q ≤ 4`.
///
/// # Errors
/// Returns a lowercase message when an order exceeds the cap, the series is
/// too short for the requested differencing or order, or a least-squares stage
/// cannot be solved.
pub fn arima_fit(
    x: &[f64],
    p: usize,
    d: usize,
    q: usize,
) -> Result<ArimaFit, String> {
    if p > MAX_AR {
        return Err(format!(
            "autoregressive order above {MAX_AR} is not supported"
        ));
    }
    if q > MAX_MA {
        return Err(format!(
            "moving-average order above {MAX_MA} is not supported"
        ));
    }
    let n = x.len();
    if n == 0 {
        return Err("series is empty".to_string());
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err("series contains a non-finite value".to_string());
    }
    if d >= n {
        return Err("series too short for the requested differencing".to_string());
    }
    let w = difference(x, d);
    let len = w.len();
    let k = 1 + p + q;
    if len < k + 2 {
        return Err("series too short for the requested order".to_string());
    }

    // Stage 1: a long autoregression whose residuals proxy the innovations.
    let target = (p.max(q + 1) + 1).max((len as f64).ln().ceil() as usize + 1);
    let max_m = len.saturating_sub(k + 2).min((len - 1) / 2);
    if max_m == 0 {
        return Err("series too short for the requested order".to_string());
    }
    let m = target.min(max_m).max(1);

    let mut ar_design = Mat::zeros(len - m, m + 1);
    let mut ar_rhs = vec![0.0; len - m];
    for (row, t) in (m..len).enumerate() {
        ar_design.set(row, 0, 1.0);
        for j in 1..=m {
            ar_design.set(row, j, w[t - j]);
        }
        ar_rhs[row] = w[t];
    }
    let ar_coefficients =
        lstsq(&ar_design, &ar_rhs).map_err(|err| format!("long autoregression failed: {err}"))?;

    let mean = mean_of(&w);
    let mut innovations = vec![0.0; len];
    for t in m..len {
        let mut prediction = ar_coefficients.first().copied().unwrap_or(0.0);
        for j in 1..=m {
            let lagged = if t >= j { w[t - j] } else { mean };
            prediction += ar_coefficients.get(j).copied().unwrap_or(0.0) * lagged;
        }
        innovations[t] = w[t] - prediction;
    }

    // Stage 2: regress on lagged values and lagged residuals.
    let first_row = m + 1;
    let rows = len - first_row;
    let mut design = Mat::zeros(rows, k);
    let mut rhs = vec![0.0; rows];
    for (row, t) in (first_row..len).enumerate() {
        design.set(row, 0, 1.0);
        for i in 1..=p {
            design.set(row, i, w[t - i]);
        }
        for j in 1..=q {
            design.set(row, p + j, innovations[t - j]);
        }
        rhs[row] = w[t];
    }
    let mut parameters = match lstsq(&design, &rhs) {
        | Ok(values) => values,
        | Err(_) => {
            let mut fallback = vec![0.0; k];
            fallback[0] = mean;
            fallback
        },
    };

    // Stage 3: Gauss–Newton on the conditional sum of squares.
    let (mut residuals, mut sse) = conditional_residuals(&w, &parameters, p, q);
    for _ in 0..100 {
        let mut jacobian = Mat::zeros(len, k);
        for column in 0..k {
            let step = 1e-6 * (1.0 + parameters[column].abs());
            let mut up = parameters.clone();
            up[column] += step;
            let mut down = parameters.clone();
            down[column] -= step;
            let (up_residuals, _) = conditional_residuals(&w, &up, p, q);
            let (down_residuals, _) = conditional_residuals(&w, &down, p, q);
            for t in 0..len {
                let slope = (up_residuals[t] - down_residuals[t]) / (2.0 * step);
                jacobian.set(t, column, slope);
            }
        }
        let mut target_vector = vec![0.0; len];
        for t in 0..len {
            target_vector[t] = -residuals[t];
        }
        let delta = match lstsq(&jacobian, &target_vector) {
            | Ok(values) => values,
            | Err(_) => break,
        };

        let mut scale = 1.0;
        let mut improved = false;
        for _ in 0..40 {
            let candidate: Vec<f64> = parameters
                .iter()
                .zip(delta.iter())
                .map(|(value, step)| value + scale * step)
                .collect();
            let (candidate_residuals, candidate_sse) = conditional_residuals(&w, &candidate, p, q);
            if candidate_sse.is_finite() && candidate_sse <= sse {
                parameters = candidate;
                residuals = candidate_residuals;
                sse = candidate_sse;
                improved = true;
                break;
            }
            scale *= 0.5;
        }
        if !improved {
            break;
        }
        let largest = delta.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if largest < 1e-10 {
            break;
        }
    }

    let sigma2 = (sse / len as f64).max(f64::MIN_POSITIVE);
    let loglik = -0.5 * len as f64 * ((2.0 * std::f64::consts::PI * sigma2).ln() + 1.0);
    let aic = -2.0 * loglik + 2.0 * (k + 1) as f64;

    let mut ar = vec![0.0; p];
    for (slot, value) in ar.iter_mut().zip(parameters.iter().skip(1)) {
        *slot = *value;
    }
    let mut ma = vec![0.0; q];
    for (slot, value) in ma.iter_mut().zip(parameters.iter().skip(1 + p)) {
        *slot = *value;
    }
    let mean = parameters.first().copied().unwrap_or(0.0);

    Ok(ArimaFit {
        ar,
        ma,
        sigma2,
        aic,
        loglik,
        mean,
    })
}

/// Forecasts an ARIMA model forward from a series.
///
/// `x` must be the series the model was fitted to, i.e. already differenced
/// `d` times, so the returned values are forecasts of the same differenced
/// series. Future shocks are set to zero, so the forecasts follow the
/// conditional mean recursion.
///
/// # Errors
/// Returns a lowercase message when the horizon is zero and the series is
/// empty, or the series is shorter than the model order.
pub fn arima_forecast(
    fit: &ArimaFit,
    x: &[f64],
    horizon: usize,
) -> Result<Vec<f64>, String> {
    if horizon == 0 {
        return Ok(Vec::new());
    }
    if x.is_empty() {
        return Err("series is empty".to_string());
    }
    let p = fit.ar.len();
    let q = fit.ma.len();
    if x.len() < p.max(q) + 1 {
        return Err("series too short for the model order".to_string());
    }

    let mut parameters = Vec::with_capacity(1 + p + q);
    parameters.push(fit.mean);
    parameters.extend_from_slice(&fit.ar);
    parameters.extend_from_slice(&fit.ma);
    let (mut residuals, _) = conditional_residuals(x, &parameters, p, q);

    let mut history = x.to_vec();
    let mut out = Vec::with_capacity(horizon);
    for _ in 0..horizon {
        let t = history.len();
        let mut prediction = fit.mean;
        for i in 1..=p {
            prediction += fit.ar[i - 1] * history[t - i];
        }
        for j in 1..=q {
            prediction += fit.ma[j - 1] * residuals[t - j];
        }
        out.push(prediction);
        history.push(prediction);
        residuals.push(0.0);
    }
    Ok(out)
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
    fn arima_one_one_zero_continues_a_linear_trend() {
        // A random walk with drift: differencing leaves white noise around the
        // drift, which is exactly the stationary series ARIMA(1,1,0) models.
        let mut noise = lcg(3_141_592);
        let n = 300;
        let mut x = vec![0.0; n];
        x[0] = 2.0;
        for t in 1..n {
            x[t] = x[t - 1] + 0.5 + 0.05 * noise();
        }

        let fit = arima_fit(&x, 1, 1, 0).expect("arima(1,1,0)");
        assert!((fit.mean - 0.5).abs() < 0.1, "drift {}", fit.mean);
        assert!(fit.ar[0].abs() < 0.3, "ar {}", fit.ar[0]);

        let w = difference(&x, 1);
        let increments = arima_forecast(&fit, &w, 10).expect("forecast");
        assert_eq!(increments.len(), 10);
        for value in &increments {
            assert!((value - 0.5).abs() < 0.1, "increment {value}");
        }

        // Integrating the increments from the last observation continues the
        // trend: each step advances by about the fitted drift.
        let mut path = vec![x[n - 1]];
        for value in &increments {
            let last = *path.last().unwrap_or(&0.0);
            path.push(last + value);
        }
        for window in path.windows(2) {
            assert!((window[1] - window[0] - 0.5).abs() < 0.1);
        }
    }

    #[test]
    fn white_noise_gives_a_small_autoregressive_coefficient() {
        let mut noise = lcg(2_718_281);
        let x: Vec<f64> = (0..400).map(|_| noise()).collect();
        let fit = arima_fit(&x, 1, 0, 0).expect("arima(1,0,0)");
        assert!(fit.ar[0].abs() < 0.2, "ar {}", fit.ar[0]);
        assert!(fit.sigma2 > 0.0);
        assert!(fit.aic.is_finite());
        assert!(fit.loglik.is_finite());
    }

    #[test]
    fn moving_average_model_fits_and_forecasts() {
        let mut noise = lcg(1_618_033);
        let mut x = vec![0.0; 500];
        let mut previous = 0.0;
        for value in &mut x {
            let shock = noise();
            *value = shock + 0.5 * previous;
            previous = shock;
        }
        let fit = arima_fit(&x, 0, 0, 1).expect("arima(0,0,1)");
        assert!((fit.ma[0] - 0.5).abs() < 0.15, "ma {:?}", fit.ma);
        let forecast = arima_forecast(&fit, &x, 3).expect("forecast");
        assert_eq!(forecast.len(), 3);
        assert!(forecast.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn invalid_orders_are_rejected() {
        let x: Vec<f64> = (0..50).map(|t| t as f64).collect();
        assert!(arima_fit(&x, 7, 0, 0).is_err());
        assert!(arima_fit(&x, 0, 0, 5).is_err());
        assert!(arima_fit(&[], 1, 0, 0).is_err());
        assert!(arima_fit(&[1.0, 2.0, 3.0], 1, 0, 0).is_err());
        assert!(arima_fit(&[1.0, 2.0, 3.0], 0, 3, 0).is_err());
    }

    #[test]
    fn forecast_of_an_empty_series_is_rejected() {
        let fit = ArimaFit {
            ar: vec![0.5],
            ma: Vec::new(),
            sigma2: 1.0,
            aic: 0.0,
            loglik: 0.0,
            mean: 0.0,
        };
        assert!(arima_forecast(&fit, &[], 3).is_err());
        assert!(arima_forecast(&fit, &[1.0], 3).is_err());
        assert!(arima_forecast(&fit, &[1.0, 2.0], 0).is_ok());
    }
}
