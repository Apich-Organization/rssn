use std::fmt;

use rustfft::FftPlanner;
use rustfft::num_complex::Complex;

/// Error in a checked signal transform, spectral estimate or resampler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignalError {
    /// The window is empty, all zero, or contains a non-finite coefficient.
    InvalidWindow,
    /// The hop size must be positive.
    InvalidHop,
    /// The sample rate must be finite and positive.
    InvalidSampleRate,
    /// The input does not contain a complete frame.
    InsufficientSamples,
    /// A spectrum has a different length than its analysis window.
    FrameShape,
    /// A signal, spectrum or filter coefficient is non-finite.
    NonFiniteInput,
    /// An output or intermediate arithmetic result is not representable.
    NonFiniteResult,
    /// No nonzero window weight covers this output sample.
    UncoveredSample(usize),
    /// An upsampling or downsampling factor is zero.
    InvalidRate,
    /// The output size overflowed or could not be allocated.
    OutputSize,
}

impl fmt::Display for SignalError {
    fn fmt(
        &self,
        formatter: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::InvalidWindow => formatter.write_str("window must be finite and nonzero"),
            | Self::InvalidHop => formatter.write_str("hop size must be positive"),
            | Self::InvalidSampleRate => {
                formatter.write_str("sample rate must be finite and positive")
            },
            | Self::InsufficientSamples => formatter.write_str("no complete frame"),
            | Self::FrameShape => formatter.write_str("spectrum and window lengths differ"),
            | Self::NonFiniteInput => formatter.write_str("signal inputs must be finite"),
            | Self::NonFiniteResult => formatter.write_str("non-finite signal result"),
            | Self::UncoveredSample(index) => {
                write!(formatter, "window does not cover sample {index}")
            },
            | Self::InvalidRate => formatter.write_str("resampling factors must be positive"),
            | Self::OutputSize => formatter.write_str("signal output size is not allocatable"),
        }
    }
}

impl std::error::Error for SignalError {}

/// One-sided power spectral density estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerSpectrum {
    /// Frequencies from DC through the highest nonnegative DFT bin.
    pub frequencies: Vec<f64>,
    /// Signal power per unit frequency, in the same order as `frequencies`.
    pub density: Vec<f64>,
}

#[allow(
    clippy::float_cmp,
    reason = "an all-zero window is specifically invalid"
)]
fn window_scale(window: &[f64]) -> Result<f64, SignalError> {
    let scale = window.iter().map(|value| value.abs()).fold(0.0, f64::max);
    if scale == 0.0 || window.iter().any(|value| !value.is_finite()) {
        return Err(SignalError::InvalidWindow);
    }
    Ok(scale)
}

/// Estimates the one-sided power spectral density of real samples by Welch's method.
///
/// Averages the periodograms of complete windowed frames, whose starts are
/// `hop` samples apart. Pass the window coefficients explicitly, for example
/// with [`super::hann_window`]. The transform uses the window length and
/// performs no detrending. Frequencies have the same unit as `sample_rate`.
/// Only bins with a negative-frequency partner are doubled. Multiplying each
/// bin by the frequency spacing and summing gives the mean window-weighted
/// frame power. Trailing samples that do not complete a frame are ignored.
///
/// # Errors
///
/// Returns an error for an invalid window, hop or sample rate, insufficient
/// samples, non-finite data in a complete frame or non-finite computed power.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "the spectral formulas use nonzero frame lengths and checked finite results"
)]
pub fn welch_psd(
    samples: &[f64],
    window: &[f64],
    hop: usize,
    sample_rate: f64,
) -> Result<PowerSpectrum, SignalError> {
    if hop == 0 {
        return Err(SignalError::InvalidHop);
    }
    if !sample_rate.is_finite() || sample_rate <= 0.0 {
        return Err(SignalError::InvalidSampleRate);
    }
    let scale = window_scale(window)?;
    let frame_len = window.len();
    if samples.len() < frame_len {
        return Err(SignalError::InsufficientSamples);
    }

    let mut weights: Vec<_> = window.iter().map(|value| value / scale).collect();
    let energy_scale = weights
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    for weight in &mut weights {
        *weight /= energy_scale;
    }

    let rate_scale = sample_rate.sqrt();
    let fft = FftPlanner::new().plan_fft_forward(frame_len);
    let mut spectrum = vec![Complex::ZERO; frame_len];
    let mut scratch = vec![Complex::ZERO; fft.get_inplace_scratch_len()];
    let bins = frame_len / 2 + 1;
    let mut density = vec![0.0; bins];
    for (index, frame) in samples.windows(frame_len).step_by(hop).enumerate() {
        for ((value, &sample), &weight) in spectrum.iter_mut().zip(frame).zip(&weights) {
            if !sample.is_finite() {
                return Err(SignalError::NonFiniteInput);
            }
            *value = Complex::new(sample * weight, 0.0);
        }
        fft.process_with_scratch(&mut spectrum, &mut scratch);
        for (bin, (mean, value)) in density.iter_mut().zip(&spectrum).enumerate() {
            let fold = if bin == 0 || (frame_len.is_multiple_of(2) && bin == bins - 1) {
                1.0
            } else {
                std::f64::consts::SQRT_2
            };
            let amplitude = (value.norm() / rate_scale) * fold;
            let power = amplitude * amplitude;
            if !power.is_finite() {
                return Err(SignalError::NonFiniteResult);
            }
            *mean += (power - *mean) / (index + 1) as f64;
        }
    }

    let frequencies = (0..bins)
        .map(|bin| (bin as f64 / frame_len as f64) * sample_rate)
        .collect();
    Ok(PowerSpectrum { frequencies, density })
}

#[derive(Clone, Default)]
struct Overlap {
    sum: f64,
    weight_squared: f64,
    scale: f64,
}

/// Reconstructs a real signal from full [`super::stft`] spectra.
///
/// Applies the same window on synthesis and divides each sample by its sum of
/// squared window weights. Returns the interval covered by the frames, of
/// length `(frames.len() - 1) * hop + window.len()`. The imaginary part after
/// the inverse transform is discarded. Use a window and hop that have
/// nonzero overlap weights across that interval, or pad the analysis when
/// endpoint samples must be preserved.
///
/// # Errors
///
/// Returns an error for invalid parameters, inconsistent or non-finite frames,
/// uncovered samples, output size overflow or non-finite reconstruction.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "overlap-add indices are bounded by the checked reconstruction length"
)]
pub fn inverse_stft(
    frames: &[Vec<Complex<f64>>],
    window: &[f64],
    hop: usize,
) -> Result<Vec<f64>, SignalError> {
    if hop == 0 {
        return Err(SignalError::InvalidHop);
    }
    window_scale(window)?;
    let frame_len = window.len();
    let offsets = frames
        .len()
        .checked_sub(1)
        .ok_or(SignalError::InsufficientSamples)?;
    let len = offsets
        .checked_mul(hop)
        .and_then(|last| last.checked_add(frame_len))
        .ok_or(SignalError::OutputSize)?;
    for frame in frames {
        if frame.len() != frame_len {
            return Err(SignalError::FrameShape);
        }
        if frame
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(SignalError::NonFiniteInput);
        }
    }

    if offsets != 0 && hop > frame_len {
        return Err(SignalError::UncoveredSample(frame_len));
    }
    let mut overlapping = Vec::new();
    overlapping
        .try_reserve_exact(len)
        .map_err(|_| SignalError::OutputSize)?;
    overlapping.resize(len, Overlap::default());
    let fft = FftPlanner::new().plan_fft_inverse(frame_len);
    let mut time = vec![Complex::ZERO; frame_len];
    let mut scratch = vec![Complex::ZERO; fft.get_inplace_scratch_len()];
    let inverse_len = 1.0 / frame_len as f64;
    for (index, frame) in frames.iter().enumerate() {
        let offset = index * hop;
        // Normalize before transforming to avoid overflow in an otherwise
        // representable inverse DFT.
        for (value, &source) in time.iter_mut().zip(frame) {
            *value = source * inverse_len;
        }
        fft.process_with_scratch(&mut time, &mut scratch);
        for ((overlap, sample), &weight) in
            overlapping.iter_mut().skip(offset).zip(&time).zip(window)
        {
            if weight == 0.0 {
                continue;
            }
            let magnitude = weight.abs();
            if magnitude > overlap.scale {
                // Scale by the largest window coefficient at this sample, so
                // even very small nonzero contributions can be normalized.
                let ratio = overlap.scale / magnitude;
                overlap.sum *= ratio;
                overlap.weight_squared *= ratio * ratio;
                overlap.scale = magnitude;
            }

            let weight = weight / overlap.scale;
            overlap.sum += sample.re * weight;
            overlap.weight_squared += weight * weight;
        }
    }

    let mut output = zeroed_output(len)?;
    for (index, (value, overlap)) in output.iter_mut().zip(overlapping).enumerate() {
        if overlap.weight_squared == 0.0 {
            return Err(SignalError::UncoveredSample(index));
        }
        *value = (overlap.sum / overlap.weight_squared) / overlap.scale;
        if !value.is_finite() {
            return Err(SignalError::NonFiniteResult);
        }
    }
    Ok(output)
}

fn zeroed_output(len: usize) -> Result<Vec<f64>, SignalError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(len)
        .map_err(|_| SignalError::OutputSize)?;
    output.resize(len, 0.0);
    Ok(output)
}

/// Upsamples, applies a finite impulse response filter and downsamples.
///
/// This implements zero insertion by `up`, linear convolution with `taps`,
/// then retention of every `down`th sample starting at index zero. It works
/// directly on filter phases without allocating the upsampled signal or
/// multiplying by its inserted zeros. The boundaries are zero extended. The
/// supplied taps determine gain, passband, alias rejection and delay.
///
/// If either input slice is empty the result is empty. Otherwise the output
/// length is `floor(((samples.len() - 1) * up + taps.len() - 1) / down) + 1`.
///
/// # Errors
///
/// Returns an error for a zero factor, non-finite data, non-finite output
/// or an output extent that overflows or cannot be allocated.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "resampling indices are bounded by checked support and there are no zero sampling factors"
)]
pub fn upfirdn(
    samples: &[f64],
    taps: &[f64],
    up: usize,
    down: usize,
) -> Result<Vec<f64>, SignalError> {
    if up == 0 || down == 0 {
        return Err(SignalError::InvalidRate);
    }
    if samples.is_empty() || taps.is_empty() {
        return Ok(Vec::new());
    }
    if samples.iter().chain(taps).any(|value| !value.is_finite()) {
        return Err(SignalError::NonFiniteInput);
    }

    // On 32- and 64-bit targets, the sum of these products fits u128 even
    // when the upsampled signal would be much too large to materialize.
    let end = (samples.len() as u128 - 1) * up as u128 + taps.len() as u128 - 1;
    let len = usize::try_from(end / down as u128 + 1).map_err(|_| SignalError::OutputSize)?;
    let mut output = zeroed_output(len)?;
    for (index, result) in output.iter_mut().enumerate() {
        let output_time = index as u128 * down as u128;
        let head = output_time / up as u128;
        let phase = (output_time % up as u128) as usize;
        let first = head.saturating_sub(samples.len() as u128 - 1);
        let first_tap = phase as u128 + first * up as u128;
        if first_tap >= taps.len() as u128 {
            continue;
        }
        let last = head.min((taps.len().saturating_sub(1).saturating_sub(phase) / up) as u128);
        let count = (last - first + 1) as usize;
        let last_sample = (head - first) as usize;
        let filter = taps.iter().skip(first_tap as usize).step_by(up).take(count);
        let inputs = samples.iter().rev().skip(samples.len() - 1 - last_sample);
        *result = filter.zip(inputs).map(|(tap, sample)| tap * sample).sum();
        if !result.is_finite() {
            return Err(SignalError::NonFiniteResult);
        }
    }
    Ok(output)
}
