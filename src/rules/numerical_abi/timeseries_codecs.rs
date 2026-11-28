// Dimensional contracts for time-series statistics and Gaussian state-space models.
use crate::kernels::timeseries;

struct KalmanCovariances(Vec<Vec<Vec<f64>>>);

impl Output for KalmanCovariances {
    fn ty() -> TypeExpr {
        tensor_type(3)
    }

    fn encode(self) -> Result<Value, OpError> {
        let state = self.0.first().map_or(0, Vec::len);
        output_tensor(
            vec![self.0.len(), state, state],
            self.0.into_iter().flatten().flatten().collect(),
        )
    }
}

struct KalmanHistory {
    filtered: Vec<Vec<f64>>,
    predicted: Vec<Vec<f64>>,
    covariance: KalmanCovariances,
    loglik: f64,
}
record!(KalmanHistory { filtered: Vec<Vec<f64>>, predicted: Vec<Vec<f64>>, covariance: KalmanCovariances, loglik: f64 });

fn kalman_error(message: String) -> OpError {
    match crate::cancellation::checkpoint() {
        | Err(cancelled) => cancelled.into(),
        | Ok(()) => OpError::Domain(message),
    }
}

fn kalman_history_budget(
    observations: usize,
    state: usize,
) -> Result<(), OpError> {
    let entries = observations.saturating_mul(state.saturating_mul(state.saturating_add(2)));
    if entries >= super::native_object_api::VALUE_LIMIT {
        return Err(OpError::Domain(
            "Kalman history exceeds the native tensor size limit".into(),
        ));
    }
    Ok(())
}

binding!(
    timeseries_integrate,
    "num_timeseries_integrate",
    "numerical_vector",
    "Undo repeated first differences using the original leading samples. Input vectors share a dimension, which the result preserves. The returned length is the sum of the input lengths.",
    "crate::kernels::timeseries::acf::integrate",
    DimensionPolicy::EqualPreserve,
    (differences: Vec<f64>, initial: Vec<f64>) -> Result<Vec<f64>, OpError>,
    |differences: Vec<f64>, initial: Vec<f64>| {
        if differences.len().saturating_add(initial.len()) >= crate::rules::numerical_abi::native_object_api::VALUE_LIMIT {
            return Err(OpError::Domain("integrated series exceeds the native tensor size limit".into()));
        }
        Ok(super::timeseries::integrate(&differences, &initial))
    }
);

binding!(
    timeseries_kalman_filter,
    "num_timeseries_kalman_filter",
    "numerical_vector",
    "Filter a homogeneous-dimensional linear Gaussian state-space model with explicit initial mean and covariance. Return filtered means, predicted means, filtered covariance matrices, and log likelihood relative to a positive observation scale.",
    "crate::kernels::timeseries::kalman::kalman_filter",
    DimensionPolicy::KalmanFilter,
    (y: Vec<Vec<f64>>, f: Vec<Vec<f64>>, h: Vec<Vec<f64>>, q: Vec<Vec<f64>>, r: Vec<Vec<f64>>, x0: Vec<f64>, p0: Vec<Vec<f64>>, observation_scale: f64) -> Result<super::KalmanHistory, OpError>,
    |y: Vec<Vec<f64>>, f: Vec<Vec<f64>>, h: Vec<Vec<f64>>, q: Vec<Vec<f64>>, r: Vec<Vec<f64>>, x0: Vec<f64>, p0: Vec<Vec<f64>>, observation_scale: f64| {
        super::kalman_history_budget(y.len(), x0.len())?;
        if observation_scale <= 0.0 {
            return Err(OpError::Domain("observation scale must be positive".into()));
        }
        let result = super::timeseries::kalman_filter(&y, &f, &h, &q, &r, &x0, &p0)
            .map_err(super::kalman_error)?;
        // Density of y/scale relative to the normalized observation axis.
        let scalars = y.len().saturating_mul(h.len()) as f64;
        Ok(super::KalmanHistory {
            filtered: result.filtered,
            predicted: result.predicted,
            covariance: super::KalmanCovariances(result.covariance),
            loglik: scalars.mul_add(observation_scale.ln(), result.loglik),
        })
    }
);

binding!(
    timeseries_kalman_smooth,
    "num_timeseries_kalman_smooth",
    "numerical_vector",
    "Rauch–Tung–Striebel smoothed state means with explicit initial mean and covariance. State components share a dimension and observations may have a different common dimension. Returns one state per observation.",
    "crate::kernels::timeseries::kalman::kalman_smooth",
    DimensionPolicy::KalmanSmooth,
    (y: Vec<Vec<f64>>, f: Vec<Vec<f64>>, h: Vec<Vec<f64>>, q: Vec<Vec<f64>>, r: Vec<Vec<f64>>, x0: Vec<f64>, p0: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, OpError>,
    |y: Vec<Vec<f64>>, f: Vec<Vec<f64>>, h: Vec<Vec<f64>>, q: Vec<Vec<f64>>, r: Vec<Vec<f64>>, x0: Vec<f64>, p0: Vec<Vec<f64>>| {
        super::kalman_history_budget(y.len(), x0.len())?;
        super::timeseries::kalman_smooth(&y, &f, &h, &q, &r, &x0, &p0).map_err(super::kalman_error)
    }
);

struct RollingMoments {
    mean: Vec<f64>,
    variance: Vec<f64>,
}
record!(RollingMoments { mean: Vec<f64>, variance: Vec<f64> });

binding!(
    timeseries_rolling_moments,
    "num_timeseries_rolling_moments",
    "numerical_vector",
    "Compute means and variances for complete windows of a finite real time series. Use ddof 0 for population variance or ddof 1 for sample variance. The result contains empty vectors when the positive window width exceeds the input length. Means retain sample units and variances square them.",
    "crate::kernels::timeseries::rolling::rolling_moments",
    DimensionPolicy::RollingMoments,
    (samples: Vec<f64>, window: usize, ddof: usize) -> Result<super::RollingMoments, OpError>,
    |samples: Vec<f64>, window: usize, ddof: usize| {
        if window > ddof && samples.len().checked_sub(window).is_some_and(|tail| {
            (tail + 1).saturating_mul(2) >= crate::rules::numerical_abi::native_object_api::VALUE_LIMIT
        }) {
            return Err(OpError::Domain("rolling moments exceed the native tensor size limit".into()));
        }
        let moments = super::timeseries::rolling::rolling_moments(&samples, window, ddof)
            .map_err(|error| match error {
                super::timeseries::rolling::RollingError::Cancelled => crate::cancellation::Cancelled.into(),
                _ => OpError::Domain(error.to_string()),
            })?;
        let (mean, variance) = moments
            .into_iter()
            .map(|moment| (moment.mean, moment.variance))
            .unzip();
        Ok(super::RollingMoments { mean, variance })
    }
);
