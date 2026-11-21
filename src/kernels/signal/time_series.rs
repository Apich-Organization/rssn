use super::SignalError;
use crate::kernels::fftconv::autocorrelate_prefix;
use crate::kernels::numeric_common::scaled_difference;
use crate::kernels::par;
use crate::kernels::reuse::f64_scratch;

#[allow(
    clippy::arithmetic_side_effects,
    reason = "covariances are calculated on bounded samples before restoring physical units"
)]
pub(super) fn autocovariance(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    let (mut sums, scale) = scaled_lag_sums(samples, max_lag)?;
    for sum in &mut sums {
        let covariance = ((*sum / samples.len() as f64) * scale) * scale;
        if !covariance.is_finite() {
            return Err(SignalError::NonFiniteResult);
        }
        *sum = covariance;
    }
    Ok(sums)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "normalized covariances have bounded magnitudes and positive denominators"
)]
pub(super) fn autocorrelation(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    let (mut sums, _) = scaled_lag_sums(samples, max_lag)?;
    let Some((zero_lag, rest)) = sums.split_first_mut() else {
        return Err(SignalError::DegenerateVariance);
    };
    let variance = *zero_lag;
    if variance <= 0.0 {
        return Err(SignalError::DegenerateVariance);
    }
    for sum in rest {
        *sum = (*sum / variance).clamp(-1.0, 1.0);
    }
    *zero_lag = 1.0;
    Ok(sums)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the lag count is bounded by the sample size and rescaled samples are bounded"
)]
fn scaled_lag_sums(
    samples: &[f64],
    max_lag: usize,
) -> Result<(Vec<f64>, f64), SignalError> {
    let Some(&origin) = samples.first() else {
        return Err(SignalError::InsufficientSamples);
    };
    if max_lag >= samples.len() {
        return Err(SignalError::InvalidLag);
    }
    let scale = samples
        .iter()
        .try_fold(0.0_f64, |peak, &value| {
            value.is_finite().then(|| peak.max(value.abs()))
        })
        .ok_or(SignalError::NonFiniteInput)?;
    let count = max_lag + 1;
    if scale <= 0.0 {
        return Ok((vec![0.0; count], scale));
    }

    let mut centered = f64_scratch(samples.len());
    centered.extend(
        samples
            .iter()
            .map(|&value| scaled_difference(value, origin, scale)),
    );
    let mean = centered.iter().sum::<f64>() / samples.len() as f64;
    for value in &mut *centered {
        *value -= mean;
    }

    let direct_sum = |lag| {
        centered
            .iter()
            .zip(centered.iter().skip(lag))
            .map(|(a, b)| a * b)
            .sum()
    };
    let work = samples.len().saturating_mul(count);
    let fft_work = samples
        .len()
        .checked_add(max_lag)
        .and_then(usize::checked_next_power_of_two)
        .map_or(usize::MAX, |len| {
            len.saturating_mul(len.ilog2() as usize).saturating_mul(4)
        });
    if work > fft_work {
        let mut sums = autocorrelate_prefix(&centered, max_lag).ok_or(SignalError::OutputSize)?;
        if let Some(variance) = sums.first_mut() {
            // Use the direct norm for the critical normalization coefficient.
            *variance = direct_sum(0);
        }
        return Ok((sums, scale));
    }
    let sums = if work >= 262_144 {
        par::map_indexed(count, 8, direct_sum)
    } else {
        (0..count).map(direct_sum).collect()
    };
    Ok((sums, scale))
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "recursion slices follow the current lag and errors are checked before division"
)]
pub(super) fn partial_autocorrelation(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    let correlation = autocorrelation(samples, max_lag)?;
    let mut partial = Vec::with_capacity(correlation.len());
    partial.push(1.0);
    let mut coefficients = f64_scratch(max_lag);
    let mut previous = f64_scratch(max_lag);
    let mut prediction_variance = 1.0;
    for (order, &observed) in correlation.iter().enumerate().skip(1) {
        let predicted = coefficients
            .iter()
            .zip(correlation[1..order].iter().rev())
            .map(|(coefficient, lag)| coefficient * lag)
            .sum::<f64>();
        let reflection = (observed - predicted) / prediction_variance;
        let magnitude = reflection.abs();
        if !magnitude.is_finite() || magnitude >= 1.0 {
            return Err(SignalError::DegenerateVariance);
        }
        previous.clear();
        previous.extend(
            coefficients
                .iter()
                .zip(coefficients.iter().rev())
                .map(|(forward, backward)| forward - reflection * backward),
        );
        previous.push(reflection);
        std::mem::swap(&mut coefficients, &mut previous);
        partial.push(reflection);
        prediction_variance *= (1.0 - magnitude) * (1.0 + magnitude);
    }
    Ok(partial)
}
