//! Time-series models.
//!
//! | operator | value |
//! |---|---|
//! | `acf(x, max_lag)` | the sample autocorrelation of `x` from lag 0 to `max_lag` |
//! | `pacf(x, max_lag)` | the partial autocorrelation, by Durbin–Levinson |
//! | `difference(x, d)` | `d`-fold first difference |
//! | `ar_fit(x, p)` | `list(phi, sigma2, aic, intercept)` of an AR(`p`) by Yule–Walker |
//! | `arima_fit(x, p, d, q)` | `list(ar, ma, sigma2, aic, loglik, mean)` of an ARIMA(`p`, `d`, `q`) by Hannan–Rissanen with a conditional-sum-of-squares refinement |
//! | `arima_forecast(x, p, d, q, horizon)` | the next `horizon` values of the series |
//! | `holt_winters(x, season, alpha, beta, gamma)` | `list(fitted, level, trend, season, sse)` of additive Holt–Winters |
//! | `holt_winters_mul(x, season, alpha, beta, gamma)` | the multiplicative form |
//! | `kalman_filter(y, transition, observation, process_cov, noise_cov)` | `list(filtered, loglik)` of a linear Gaussian state-space model from a zero initial state and unit initial covariance |
//!
//! The orders are capped (`p <= 6`, `q <= 4`) and a series too short for the
//! requested order is a [`DiagKind::Domain`] finding rather than a fit of a
//! model the data cannot identify.

use super::bridge::Data;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::require;
use super::bridge::rows;
use super::bridge::scalar;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;

/// Longest series a model accepts.
const MAX_LEN: usize = 100_000;
/// Largest autoregressive order.
const MAX_P: usize = 6;
/// Largest moving-average order.
const MAX_Q: usize = 4;

/// Every operator of this family.
reviewed_entries! {
    (
        "acf", "vd", acf, timeseries::acf,
        "`acf(x, max_lag)`: the sample autocorrelation of `x` for lags `0..=max_lag`."
    ),
    (
        "pacf", "vd", pacf, timeseries::pacf,
        "`pacf(x, max_lag)`: the partial autocorrelation of `x` for lags `1..=max_lag`, by the Durbin–Levinson recursion."
    ),
    (
        "difference", "vd", difference, timeseries::difference,
        "`difference(x, d)`: the `d`-fold first difference of `x`."
    ),
    (
        "ar_fit", "vd", ar_fit, timeseries::ar_yule_walker,
        "`ar_fit(x, p)`: `list(phi, sigma2, aic, intercept)` of an autoregression of order `p` fitted by the Yule–Walker equations."
    ),
    (
        "arima_fit", "vddd", arima_fit, timeseries::arima_fit,
        "`arima_fit(x, p, d, q)`: `list(ar, ma, sigma2, aic, loglik, mean)` of an ARIMA(`p`, `d`, `q`) fitted by Hannan–Rissanen and refined by conditional sum of squares."
    ),
    (
        "arima_forecast", "vdddd", arima_forecast, timeseries::arima_forecast,
        "`arima_forecast(x, p, d, q, horizon)`: the next `horizon` values of `x` under an ARIMA(`p`, `d`, `q`) refitted to the whole series."
    ),
    (
        "holt_winters", "vdddd", holt_winters, timeseries::holt_winters_additive,
        "`holt_winters(x, season, alpha, beta, gamma)`: `list(fitted, level, trend, season, sse)` of additive Holt–Winters smoothing."
    ),
    (
        "holt_winters_mul", "vdddd", holt_winters_mul, timeseries::holt_winters_multiplicative,
        "`holt_winters_mul(x, season, alpha, beta, gamma)`: the multiplicative Holt–Winters form, `list(fitted, level, trend, season, sse)`."
    ),
    (
        "kalman_filter", "ttttt", kalman_filter, timeseries::kalman_filter,
        "`kalman_filter(y, transition, observation, process_cov, noise_cov)`: `list(filtered, loglik)` of a linear Gaussian state-space model from a zero initial state and the unit initial covariance; `y` is a list of observation vectors."
    ),
}

/// The sample at `at`.
fn sample(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    let values = vector(args, at)?;
    require(
        !values.is_empty() && values.len() <= MAX_LEN,
        Some(at),
        DiagKind::Domain,
        format_args!("expected 1 to {MAX_LEN} samples, got {}", values.len()),
    )?;
    Some(values)
}

/// A non-negative order argument.
fn order(
    args: &[Data],
    at: usize,
    max: usize,
) -> Option<usize> {
    count(args, at, max)
}

fn acf(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let lags = order(args, 1, 64)?;
    super::bridge::vec_out(crate::kernels::timeseries::acf(values, lags))
}

fn pacf(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let lags = order(args, 1, 64)?;
    super::bridge::vec_out(crate::kernels::timeseries::pacf(values, lags))
}

fn difference(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let order = order(args, 1, 8)?;
    super::bridge::vec_out(crate::kernels::timeseries::difference(values, order))
}

fn ar_fit(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let p = order(args, 1, MAX_P)?;
    require(
        p > 0,
        Some(1),
        DiagKind::Domain,
        format_args!("the order must be at least 1"),
    )?;
    let fit = match crate::kernels::timeseries::ar_yule_walker(values, p) {
        | Ok(fit) => fit,
        | Err(message) => return bad(0, DiagKind::Domain, format_args!("{message}")),
    };
    Some(Data::Tree(vec![
        Data::vector(fit.phi),
        Data::Scalar(fit.sigma2),
        Data::Scalar(fit.aic),
        Data::Scalar(fit.intercept),
    ]))
}

fn arima_fit(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let p = order(args, 1, MAX_P)?;
    let d = order(args, 2, 2)?;
    let q = order(args, 3, MAX_Q)?;
    let fit = match crate::kernels::timeseries::arima_fit(values, p, d, q) {
        | Ok(fit) => fit,
        | Err(message) => return bad(0, DiagKind::Domain, format_args!("{message}")),
    };
    Some(Data::Tree(vec![
        Data::vector(fit.ar),
        Data::vector(fit.ma),
        Data::Scalar(fit.sigma2),
        Data::Scalar(fit.aic),
        Data::Scalar(fit.loglik),
        Data::Scalar(fit.mean),
    ]))
}

fn arima_forecast(args: &[Data]) -> Option<Data> {
    let values = sample(args, 0)?;
    let p = order(args, 1, MAX_P)?;
    let d = order(args, 2, 2)?;
    let q = order(args, 3, MAX_Q)?;
    let horizon = order(args, 4, 4096)?;
    let fit = match crate::kernels::timeseries::arima_fit(values, p, d, q) {
        | Ok(fit) => fit,
        | Err(message) => return bad(0, DiagKind::Domain, format_args!("{message}")),
    };
    match crate::kernels::timeseries::arima_forecast(&fit, values, horizon) {
        | Ok(ahead) => super::bridge::vec_out(ahead),
        | Err(message) => bad(4, DiagKind::Domain, format_args!("{message}")),
    }
}

fn smoothing(
    args: &[Data],
    multiplicative: bool,
) -> Option<Data> {
    let values = sample(args, 0)?;
    let season = order(args, 1, 512)?;
    let alpha = scalar(args, 2)?;
    let beta = scalar(args, 3)?;
    let gamma = scalar(args, 4)?;
    for (at, value) in [(2, alpha), (3, beta), (4, gamma)] {
        require(
            (0.0..=1.0).contains(&value),
            Some(at),
            DiagKind::Domain,
            format_args!("a smoothing weight must lie in [0, 1], got {value}"),
        )?;
    }
    let fitted = if multiplicative {
        crate::kernels::timeseries::holt_winters_multiplicative(values, season, alpha, beta, gamma)
    } else {
        crate::kernels::timeseries::holt_winters_additive(values, season, alpha, beta, gamma)
    };
    let fit = match fitted {
        | Ok(fit) => fit,
        | Err(message) => return bad(0, DiagKind::Domain, format_args!("{message}")),
    };
    Some(Data::Tree(vec![
        Data::vector(fit.fitted),
        Data::vector(fit.level),
        Data::vector(fit.trend),
        Data::vector(fit.season),
        Data::Scalar(fit.sse),
    ]))
}

fn holt_winters(args: &[Data]) -> Option<Data> {
    smoothing(args, false)
}

fn holt_winters_mul(args: &[Data]) -> Option<Data> {
    smoothing(args, true)
}

/// A square matrix argument.
fn square(
    args: &[Data],
    at: usize,
    name: &str,
) -> Option<Vec<Vec<f64>>> {
    let Some(matrix) = args.get(at).and_then(rows) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected the {name} as a list of rows"),
        );
    };
    let size = matrix.len();
    if size == 0 || matrix.iter().any(|row| row.len() != size) {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("the {name} must be a non-empty square matrix"),
        );
    }
    Some(matrix)
}

fn kalman_filter(args: &[Data]) -> Option<Data> {
    let Some(observations) = args.first().and_then(rows) else {
        return bad(
            0,
            DiagKind::Shape,
            format_args!("expected a list of observation vectors"),
        );
    };
    require(
        !observations.is_empty() && observations.len() <= MAX_LEN,
        Some(0),
        DiagKind::Domain,
        format_args!(
            "expected 1 to {MAX_LEN} observations, got {}",
            observations.len()
        ),
    )?;
    let transition = square(args, 1, "transition")?;
    let Some(observation) = args.get(2).and_then(rows) else {
        return bad(
            2,
            DiagKind::Shape,
            format_args!("expected the observation matrix as a list of rows"),
        );
    };
    let process = square(args, 3, "process covariance")?;
    let noise = square(args, 4, "noise covariance")?;
    let states = transition.len();
    let start = vec![0.0; states];
    let mut initial = vec![vec![0.0; states]; states];
    for (index, row) in initial.iter_mut().enumerate() {
        if let Some(entry) = row.get_mut(index) {
            *entry = 1.0;
        }
    }
    match crate::kernels::timeseries::kalman_filter(
        &observations,
        &transition,
        &observation,
        &process,
        &noise,
        &start,
        &initial,
    ) {
        | Ok(result) => {
            let filtered = result
                .filtered
                .into_iter()
                .map(Data::vector)
                .collect::<Vec<Data>>();
            Some(Data::Tree(vec![
                Data::Tree(filtered),
                Data::Scalar(result.loglik),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
