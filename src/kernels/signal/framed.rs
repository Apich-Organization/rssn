use std::fmt;

use rustfft::FftDirection;
use rustfft::num_complex::Complex;

use crate::kernels::reuse::fft_plan;
use crate::kernels::reuse::fft_scratch;

/// Error in a checked signal transform, spectral estimate or resampler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignalError {
    /// The window is empty, all zero, or contains a non-finite coefficient.
    InvalidWindow,
    /// The hop size must be positive.
    InvalidHop,
    /// The sample rate must be finite and positive.
    InvalidSampleRate,
    /// The start frequency and frequency step must be finite.
    InvalidFrequency,
    /// The input does not contain a complete frame.
    InsufficientSamples,
    /// A spectrum has a different length than its analysis window.
    FrameShape,
    /// Paired signals must have the same number of samples.
    SampleShape,
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
            | Self::InvalidFrequency => formatter.write_str("frequency and step must be finite"),
            | Self::InsufficientSamples => formatter.write_str("no complete frame"),
            | Self::FrameShape => formatter.write_str("spectrum and window lengths differ"),
            | Self::SampleShape => formatter.write_str("paired signals have different lengths"),
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

/// One-sided cross-spectrum and magnitude-squared coherence of two real signals.
#[derive(Debug, Clone, PartialEq)]
pub struct CrossSpectrum {
    /// Frequencies from DC through the highest nonnegative DFT bin.
    pub frequencies: Vec<f64>,
    /// Mean `conj(X) * Y` per unit frequency after window and one-sided scaling.
    pub density: Vec<Complex<f64>>,
    /// Magnitude-squared coherence at each bin, in `[0, 1]`.
    /// Bins with zero power in either signal have coherence zero.
    pub coherence: Vec<f64>,
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

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the checked window has nonzero finite energy"
)]
fn normalized_window(window: &[f64]) -> Result<Vec<f64>, SignalError> {
    let scale = window_scale(window)?;
    let mut weights: Vec<_> = window.iter().map(|value| value / scale).collect();
    let energy = weights
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    for value in &mut weights {
        *value /= energy;
    }
    Ok(weights)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "a frequency is a fraction of a checked sample rate"
)]
fn frequencies(
    frame_len: usize,
    sample_rate: f64,
) -> Vec<f64> {
    (0..=frame_len / 2)
        .map(|bin| (bin as f64 / frame_len as f64) * sample_rate)
        .collect()
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the spectral formulas use nonzero frame lengths and checked finite results"
)]
pub(super) fn welch_psd(
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
    let weights = normalized_window(window)?;
    let frame_len = window.len();
    if samples.len() < frame_len {
        return Err(SignalError::InsufficientSamples);
    }

    let rate_scale = sample_rate.sqrt();
    let fft = fft_plan(frame_len, FftDirection::Forward);
    let mut spectrum = vec![Complex::ZERO; frame_len];
    let mut scratch = fft_scratch();
    scratch.resize(fft.get_inplace_scratch_len(), Complex::ZERO);
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

    let frequencies = frequencies(frame_len, sample_rate);
    Ok(PowerSpectrum { frequencies, density })
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "spectral statistics check every arithmetic result used for output"
)]
pub(super) fn welch_cross_spectrum(
    x: &[f64],
    y: &[f64],
    window: &[f64],
    hop: usize,
    sample_rate: f64,
) -> Result<CrossSpectrum, SignalError> {
    if x.len() != y.len() {
        return Err(SignalError::SampleShape);
    }
    if hop == 0 {
        return Err(SignalError::InvalidHop);
    }
    if !sample_rate.is_finite() || sample_rate <= 0.0 {
        return Err(SignalError::InvalidSampleRate);
    }
    let weights = normalized_window(window)?;
    let frame_len = window.len();
    if x.len() < frame_len {
        return Err(SignalError::InsufficientSamples);
    }
    let plan = fft_plan(frame_len, FftDirection::Forward);
    let mut first = vec![Complex::ZERO; frame_len];
    let mut second = vec![Complex::ZERO; frame_len];
    let mut scratch = fft_scratch();
    scratch.resize(plan.get_inplace_scratch_len(), Complex::ZERO);
    let mut statistics = vec![(0.0, 0.0, Complex::ZERO); frame_len / 2 + 1];
    let rate_scale = sample_rate.sqrt();
    let frames = x.windows(frame_len).zip(y.windows(frame_len)).step_by(hop);
    for (index, (x, y)) in frames.enumerate() {
        for ((a, b), ((&u, &v), &weight)) in first
            .iter_mut()
            .zip(&mut second)
            .zip(x.iter().zip(y).zip(&weights))
        {
            if !u.is_finite() || !v.is_finite() {
                return Err(SignalError::NonFiniteInput);
            }
            *a = Complex::new(u * weight, 0.0);
            *b = Complex::new(v * weight, 0.0);
        }
        plan.process_with_scratch(&mut first, &mut scratch);
        plan.process_with_scratch(&mut second, &mut scratch);
        let n = (index + 1) as f64;
        for (bin, ((xx, yy, xy), (a, b))) in statistics
            .iter_mut()
            .zip(first.iter().zip(&second))
            .enumerate()
        {
            let fold = if bin == 0 || (frame_len.is_multiple_of(2) && bin == frame_len / 2) {
                1.0
            } else {
                std::f64::consts::SQRT_2
            };
            let a = (a / rate_scale) * fold;
            let b = (b / rate_scale) * fold;
            let power_x = a.norm_sqr();
            let power_y = b.norm_sqr();
            let cross = a.conj() * b;
            if !power_x.is_finite()
                || !power_y.is_finite()
                || !cross.re.is_finite()
                || !cross.im.is_finite()
            {
                return Err(SignalError::NonFiniteResult);
            }
            *xx += (power_x - *xx) / n;
            *yy += (power_y - *yy) / n;
            // The sign can differ between frames. Scale the terms separately
            // so subtracting two large, opposite-signed values cannot overflow.
            *xy += cross / n - *xy / n;
        }
    }
    let (density, coherence) = statistics
        .into_iter()
        .map(|(xx, yy, cross)| {
            let magnitude = if xx > 0.0 && yy > 0.0 {
                (cross.norm() / xx.sqrt() / yy.sqrt()).min(1.0)
            } else {
                0.0
            };
            (cross, magnitude * magnitude)
        })
        .unzip();
    Ok(CrossSpectrum {
        frequencies: frequencies(frame_len, sample_rate),
        density,
        coherence,
    })
}

#[derive(Clone, Default)]
struct Overlap {
    sum: f64,
    weight_squared: f64,
    scale: f64,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "overlap-add indices are bounded by the checked reconstruction length"
)]
pub(super) fn inverse_stft(
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
    let fft = fft_plan(frame_len, FftDirection::Inverse);
    let mut time = vec![Complex::ZERO; frame_len];
    let mut scratch = fft_scratch();
    scratch.resize(fft.get_inplace_scratch_len(), Complex::ZERO);
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

#[allow(
    clippy::arithmetic_side_effects,
    reason = "resampling indices are bounded by checked support and there are no zero sampling factors"
)]
pub(super) fn upfirdn(
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
