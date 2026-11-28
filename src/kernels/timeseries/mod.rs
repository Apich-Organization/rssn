//! Time-series models: autocorrelation, ARIMA, exponential smoothing and state-space filtering.
//!
//! The module is organised by model class: [`acf`] holds the sample
//! autocorrelation, partial autocorrelation and differencing primitives that
//! the rest builds on, [`ar`] the Yule–Walker autoregression, [`arima`] the
//! Hannan–Rissanen / conditional-sum-of-squares ARIMA fit, [`holtwinters`]
//! exponential smoothing, and [`kalman`] linear Gaussian state-space filtering
//! and smoothing.
//!
//! [`acf`]: mod@crate::kernels::timeseries::acf
//! [`ar`]: crate::kernels::timeseries::ar
//! [`arima`]: crate::kernels::timeseries::arima
//! [`holtwinters`]: crate::kernels::timeseries::holtwinters
//! [`kalman`]: crate::kernels::timeseries::kalman

pub mod acf;
pub mod ar;
pub mod arima;
pub mod holtwinters;
pub mod kalman;
pub mod rolling;

pub use acf::acf;
pub use acf::difference;
pub use acf::integrate;
pub use acf::pacf;
pub use ar::ArFit;
pub use ar::ar_yule_walker;
pub use arima::ArimaFit;
pub use arima::arima_fit;
pub use arima::arima_forecast;
pub use holtwinters::HoltWinters;
pub use holtwinters::holt_winters_additive;
pub use holtwinters::holt_winters_multiplicative;
pub use kalman::KalmanResult;
pub use kalman::kalman_filter;
pub use kalman::kalman_smooth;
