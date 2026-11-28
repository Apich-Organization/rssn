// Checked zero-phase filter with dimensionless coefficients and measured samples.
binding!(
    signal_filtfilt,
    "num_signal_filtfilt",
    "numerical_vector",
    "`num_signal_filtfilt(b, a, samples)` filters forward and backward with odd reflection. Coefficients are dimensionless. The result preserves the length and units of the samples.",
    "crate::kernels::signal::checked_filtfilt",
    DimensionPolicy::PreserveThird,
    (b: Vec<f64>, a: Vec<f64>, samples: Vec<f64>) -> Result<Vec<f64>, OpError>,
    |b: Vec<f64>, a: Vec<f64>, samples: Vec<f64>| {
        crate::kernels::signal::checked_filtfilt(&b, &a, &samples)
            .map_err(|error| OpError::Domain(error.to_string()))
    }
);
