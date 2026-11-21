use std::f64::consts::TAU;
use std::sync::Arc;

use num_complex::Complex64;
use rustfft::Fft;
use rustfft::FftDirection;

use super::SignalError;
use crate::kernels::reuse::PlanCache;
use crate::kernels::reuse::fft_plan;
use crate::kernels::reuse::fft_scratch_for;

#[allow(
    clippy::arithmetic_side_effects,
    reason = "scaling finite FFT input keeps its intermediate magnitudes representable"
)]
pub(super) fn analytic_signal(samples: &[f64]) -> Result<Vec<Complex64>, SignalError> {
    let peak = samples
        .iter()
        .try_fold(0.0_f64, |peak, &sample| {
            sample.is_finite().then(|| peak.max(sample.abs()))
        })
        .ok_or(SignalError::NonFiniteInput)?;
    let len = samples.len();
    if len <= 2 || peak <= 0.0 {
        return Ok(samples
            .iter()
            .map(|&value| Complex64::new(value, 0.0))
            .collect());
    }
    // Use a bounded input scale to reduce overflow and underflow in the
    // butterflies, including for constant signals near the limits of f64.
    let amplitude = if (1e-150..=1e150).contains(&peak) {
        1.0
    } else {
        peak
    };
    let mut spectrum: Vec<_> = samples
        .iter()
        .map(|&sample| Complex64::new(sample / amplitude, 0.0))
        .collect();
    let forward = fft_plan(len, FftDirection::Forward);
    let inverse = fft_plan(len, FftDirection::Inverse);
    let scratch_len = forward
        .get_inplace_scratch_len()
        .max(inverse.get_inplace_scratch_len());
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex64::ZERO);
    forward.process_with_scratch(&mut spectrum, &mut scratch);
    // DC and Nyquist contribute only to the real part, which is copied
    // directly from the input. Drop them before the inverse FFT.
    let scale = 2.0 / len as f64;
    for (bin, value) in spectrum.iter_mut().enumerate() {
        if bin > 0 && bin <= (len - 1) / 2 {
            *value *= scale;
        } else {
            *value = Complex64::ZERO;
        }
    }
    inverse.process_with_scratch(&mut spectrum, &mut scratch);
    for (value, &sample) in spectrum.iter_mut().zip(samples) {
        let hilbert = value.im * amplitude;
        if !hilbert.is_finite() {
            return Err(SignalError::NonFiniteResult);
        }
        *value = Complex64::new(sample, hilbert);
    }
    Ok(spectrum)
}

pub(super) fn hilbert(samples: &[f64]) -> Result<Vec<f64>, SignalError> {
    Ok(analytic_signal(samples)?
        .into_iter()
        .map(|value| value.im)
        .collect())
}

pub(super) fn analytic_envelope(samples: &[f64]) -> Result<Vec<f64>, SignalError> {
    analytic_signal(samples)?
        .into_iter()
        .map(|value| {
            let magnitude = value.norm();
            if magnitude.is_finite() {
                Ok(magnitude)
            } else {
                Err(SignalError::NonFiniteResult)
            }
        })
        .collect()
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::float_cmp,
    reason = "FFT work indices are bounded by the checked convolution extent"
)]
pub(super) fn zoom_fft(
    samples: &[Complex64],
    start: f64,
    step: f64,
    bins: usize,
) -> Result<Vec<Complex64>, SignalError> {
    if !start.is_finite() || !step.is_finite() {
        return Err(SignalError::InvalidFrequency);
    }
    if samples
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(SignalError::NonFiniteInput);
    }
    if bins == 0 {
        return Ok(Vec::new());
    }
    if samples.is_empty() {
        return complex_zeros(bins);
    }
    let start = start % 1.0;
    let step = step % 1.0;
    if bins == 1 || step == 0.0 || samples.len() == 1 {
        let mut result = complex_zeros(bins)?;
        let value: Complex64 = if start == 0.0 || samples.len() == 1 {
            samples.iter().sum()
        } else {
            samples
                .iter()
                .enumerate()
                .map(|(index, sample)| {
                    let angle = -TAU * product_cycles(start, index as f64);
                    sample * Complex64::from_polar(1.0, angle)
                })
                .sum()
        };
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(SignalError::NonFiniteResult);
        }
        result.fill(value);
        return Ok(result);
    }
    let size = samples
        .len()
        .checked_add(bins - 1)
        .and_then(usize::checked_next_power_of_two)
        .ok_or(SignalError::OutputSize)?;
    let half_step = step * 0.5;
    let kernel = zoom_kernel(samples.len(), bins, size, half_step)?;
    let mut input = complex_zeros(size)?;
    for (index, (value, &sample)) in input.iter_mut().zip(samples).enumerate() {
        let cycles = product_cycles(start, index as f64) + chirp_cycles(half_step, index);
        *value = sample * Complex64::from_polar(1.0, -TAU * cycles);
    }
    let mut scratch = fft_scratch_for(kernel.scratch_len);
    scratch.resize(kernel.scratch_len, Complex64::ZERO);
    kernel
        .forward
        .process_with_scratch(&mut input, &mut scratch);
    let reciprocal_size = 1.0 / size as f64;
    for (data, response) in input.iter_mut().zip(&kernel.response) {
        *data *= response * reciprocal_size;
    }
    let inverse = Arc::clone(&kernel.inverse);
    drop(kernel);
    inverse.process_with_scratch(&mut input, &mut scratch);
    drop(scratch);
    input
        .iter()
        .take(bins)
        .enumerate()
        .map(|(bin, value)| {
            let output = value * Complex64::from_polar(1.0, -TAU * chirp_cycles(half_step, bin));
            if output.re.is_finite() && output.im.is_finite() {
                Ok(output)
            } else {
                Err(SignalError::NonFiniteResult)
            }
        })
        .collect()
}

#[derive(PartialEq, Eq)]
struct ZoomKey {
    sample_count: usize,
    bins: usize,
    half_step: u64,
}

struct ZoomKernel {
    response: Vec<Complex64>,
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
    scratch_len: usize,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the caller checks the convolution length and handles degenerate grids early"
)]
fn zoom_kernel(
    sample_count: usize,
    bins: usize,
    size: usize,
    half_step: f64,
) -> Result<Arc<ZoomKernel>, SignalError> {
    static KERNELS: PlanCache<ZoomKey, ZoomKernel> = PlanCache::new(8);
    let build = || {
        let mut response = complex_zeros(size)?;
        for (bin, value) in response.iter_mut().take(bins).enumerate() {
            *value = Complex64::from_polar(1.0, TAU * chirp_cycles(half_step, bin));
        }
        for (offset, value) in response.iter_mut().rev().take(sample_count - 1).enumerate() {
            *value = Complex64::from_polar(1.0, TAU * chirp_cycles(half_step, offset + 1));
        }

        let forward = fft_plan(size, FftDirection::Forward);
        let inverse = fft_plan(size, FftDirection::Inverse);
        let scratch_len = forward
            .get_inplace_scratch_len()
            .max(inverse.get_inplace_scratch_len());
        let mut scratch = fft_scratch_for(scratch_len);
        scratch.resize(scratch_len, Complex64::ZERO);
        forward.process_with_scratch(&mut response, &mut scratch);
        Ok(Arc::new(ZoomKernel {
            response,
            forward,
            inverse,
            scratch_len,
        }))
    };
    if size <= 65_536 {
        KERNELS.try_get_or_insert_with(
            ZoomKey {
                sample_count,
                bins,
                half_step: half_step.to_bits(),
            },
            build,
        )
    } else {
        build()
    }
}

fn chirp_cycles(
    half_step: f64,
    index: usize,
) -> f64 {
    product_cycles(half_step, (index as f64).powi(2))
}

fn complex_zeros(len: usize) -> Result<Vec<Complex64>, SignalError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| SignalError::OutputSize)?;
    values.resize(len, Complex64::ZERO);
    Ok(values)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the finite operands give a representable product"
)]
fn product_cycles(
    value: f64,
    count: f64,
) -> f64 {
    // Recover the product's round-off before reducing its phase. A direct
    // multiplication can lose fractional cycles on long chirps.
    let product = value * count;
    (product % 1.0) + value.mul_add(count, -product)
}
