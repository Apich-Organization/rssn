//! Holt–Winters exponential smoothing with an additive or multiplicative
//! seasonal component.

use crate::kernels::timeseries::acf::mean_of;

/// The output of a Holt–Winters smoothing pass.
#[derive(Debug, Clone, PartialEq)]
pub struct HoltWinters {
    /// Level estimate after each observation.
    pub level: Vec<f64>,
    /// Trend estimate after each observation.
    pub trend: Vec<f64>,
    /// Seasonal component: the first `season` entries are the initial indices
    /// and entry `season + t` is the index updated at time `t`.
    pub season: Vec<f64>,
    /// One-step-ahead fitted value for each observation.
    pub fitted: Vec<f64>,
    /// Sum of squared one-step-ahead errors.
    pub sse: f64,
}

/// Validates the smoothing parameters and the series length.
fn validate(
    x: &[f64],
    season: usize,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<(), String> {
    if season == 0 {
        return Err("seasonal period must be positive".to_string());
    }
    if x.len() < 2 * season {
        return Err("series must cover at least two full seasons".to_string());
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err("series contains a non-finite value".to_string());
    }
    for (name, value) in [("alpha", alpha), ("beta", beta), ("gamma", gamma)] {
        if !(0.0..=1.0).contains(&value) {
            return Err(format!("{name} must lie in [0, 1]"));
        }
    }
    Ok(())
}

/// Starting level, trend and seasonal indices from the first two seasons.
fn initial_state(
    x: &[f64],
    season: usize,
    multiplicative: bool,
) -> Result<(f64, f64, Vec<f64>), String> {
    let first = mean_of(&x[..season]);
    let second = mean_of(&x[season..2 * season]);
    if multiplicative && (first <= 0.0 || second <= 0.0) {
        return Err("multiplicative smoothing needs strictly positive data".to_string());
    }
    let trend = (second - first) / season as f64;
    let mut indices = Vec::with_capacity(season);
    for value in &x[..season] {
        indices.push(if multiplicative {
            value / first
        } else {
            value - first
        });
    }
    Ok((first, trend, indices))
}

/// Holt–Winters smoothing with an additive seasonal component.
///
/// The recursion is
/// `l_t = α(y_t - s_{t-m}) + (1-α)(l_{t-1} + b_{t-1})`,
/// `b_t = β(l_t - l_{t-1}) + (1-β)b_{t-1}` and
/// `s_t = γ(y_t - l_{t-1} - b_{t-1}) + (1-γ)s_{t-m}`, with the first two
/// seasons used for the starting values.
///
/// # Errors
/// Returns a lowercase message when the seasonal period is zero, the series
/// covers fewer than two seasons, a value is not finite, or a smoothing
/// parameter leaves `[0, 1]`.
pub fn holt_winters_additive(
    x: &[f64],
    season: usize,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<HoltWinters, String> {
    validate(x, season, alpha, beta, gamma)?;
    let (level0, trend0, indices) = initial_state(x, season, false)?;
    let n = x.len();

    let mut level = vec![0.0; n];
    let mut trend = vec![0.0; n];
    let mut seasonal = vec![0.0; n + season];
    let mut fitted = vec![0.0; n];
    let mut sse = 0.0;
    for (slot, value) in seasonal.iter_mut().zip(indices.iter()) {
        *slot = *value;
    }

    for t in 0..n {
        let previous_level = if t == 0 {
            level0
        } else {
            level[t - 1]
        };
        let previous_trend = if t == 0 {
            trend0
        } else {
            trend[t - 1]
        };
        let previous_season = seasonal[t];

        let forecast = previous_level + previous_trend + previous_season;
        fitted[t] = forecast;
        let error = x[t] - forecast;
        sse += error * error;

        let new_level =
            alpha * (x[t] - previous_season) + (1.0 - alpha) * (previous_level + previous_trend);
        let new_trend = beta * (new_level - previous_level) + (1.0 - beta) * previous_trend;
        let new_season =
            gamma * (x[t] - previous_level - previous_trend) + (1.0 - gamma) * previous_season;

        level[t] = new_level;
        trend[t] = new_trend;
        seasonal[t + season] = new_season;
    }

    Ok(HoltWinters {
        level,
        trend,
        season: seasonal,
        fitted,
        sse,
    })
}

/// Holt–Winters smoothing with a multiplicative seasonal component.
///
/// The recursion is
/// `l_t = α(y_t / s_{t-m}) + (1-α)(l_{t-1} + b_{t-1})`,
/// `b_t = β(l_t - l_{t-1}) + (1-β)b_{t-1}` and
/// `s_t = γ(y_t / (l_{t-1} + b_{t-1})) + (1-γ)s_{t-m}`.
///
/// # Errors
/// Returns a lowercase message when the seasonal period is zero, the series
/// covers fewer than two seasons, a value is not finite or not strictly
/// positive, a smoothing parameter leaves `[0, 1]`, or a level or seasonal
/// index is zero where the recursion divides by it.
pub fn holt_winters_multiplicative(
    x: &[f64],
    season: usize,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<HoltWinters, String> {
    validate(x, season, alpha, beta, gamma)?;
    if x.iter().any(|v| *v <= 0.0) {
        return Err("multiplicative smoothing needs strictly positive data".to_string());
    }
    let (level0, trend0, indices) = initial_state(x, season, true)?;
    let n = x.len();

    let mut level = vec![0.0; n];
    let mut trend = vec![0.0; n];
    let mut seasonal = vec![0.0; n + season];
    let mut fitted = vec![0.0; n];
    let mut sse = 0.0;
    for (slot, value) in seasonal.iter_mut().zip(indices.iter()) {
        *slot = *value;
    }

    for t in 0..n {
        let previous_level = if t == 0 {
            level0
        } else {
            level[t - 1]
        };
        let previous_trend = if t == 0 {
            trend0
        } else {
            trend[t - 1]
        };
        let previous_season = seasonal[t];

        let base = previous_level + previous_trend;
        if previous_season == 0.0 || base == 0.0 {
            return Err("multiplicative recursion met a zero level or index".to_string());
        }

        let forecast = base * previous_season;
        fitted[t] = forecast;
        let error = x[t] - forecast;
        sse += error * error;

        let new_level = alpha * (x[t] / previous_season) + (1.0 - alpha) * base;
        let new_trend = beta * (new_level - previous_level) + (1.0 - beta) * previous_trend;
        let new_season = gamma * (x[t] / base) + (1.0 - gamma) * previous_season;

        level[t] = new_level;
        trend[t] = new_trend;
        seasonal[t + season] = new_season;
    }

    Ok(HoltWinters {
        level,
        trend,
        season: seasonal,
        fitted,
        sse,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATTERN: [f64; 4] = [10.0, 14.0, 9.0, 13.0];

    fn repeating(n: usize) -> Vec<f64> {
        (0..n).map(|t| PATTERN[t % 4]).collect()
    }

    #[test]
    fn additive_smoothing_recovers_a_constant_seasonal_pattern() {
        let x = repeating(64);
        let fit = holt_winters_additive(&x, 4, 0.5, 0.1, 0.5).expect("additive fit");
        assert_eq!(fit.fitted.len(), x.len());
        assert_eq!(fit.season.len(), x.len() + 4);
        for t in 0..x.len() {
            assert!(
                (fit.fitted[t] - x[t]).abs() < 1e-3,
                "t {t}: {} vs {}",
                fit.fitted[t],
                x[t]
            );
        }
        assert!(fit.sse < 1e-6, "sse {}", fit.sse);
    }

    #[test]
    fn additive_smoothing_keeps_the_seasonal_shape() {
        let x = repeating(40);
        let fit = holt_winters_additive(&x, 4, 0.4, 0.1, 0.4).expect("additive fit");
        let level = fit.level[39];
        let season: Vec<f64> = (0..4).map(|i| fit.season[40 + i]).collect();
        for i in 0..4 {
            let reconstructed = level + season[i];
            assert!((reconstructed - PATTERN[i]).abs() < 1e-3, "{reconstructed}");
        }
    }

    #[test]
    fn multiplicative_smoothing_recovers_a_proportional_pattern() {
        let base = [10.0, 14.0, 9.0, 13.0];
        let x: Vec<f64> = (0..64).map(|t| base[t % 4]).collect();
        let fit = holt_winters_multiplicative(&x, 4, 0.5, 0.1, 0.5).expect("multiplicative fit");
        for t in 0..x.len() {
            assert!(
                (fit.fitted[t] - x[t]).abs() < 1e-3,
                "t {t}: {} vs {}",
                fit.fitted[t],
                x[t]
            );
        }
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let x = repeating(12);
        assert!(holt_winters_additive(&x, 0, 0.5, 0.1, 0.5).is_err());
        assert!(holt_winters_additive(&x, 4, 1.5, 0.1, 0.5).is_err());
        assert!(holt_winters_additive(&x, 8, 0.5, 0.1, 0.5).is_err());
        assert!(holt_winters_additive(&[f64::NAN; 12], 4, 0.5, 0.1, 0.5).is_err());
        assert!(holt_winters_multiplicative(&[0.0; 12], 4, 0.5, 0.1, 0.5).is_err());
    }
}
