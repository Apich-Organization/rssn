use std::fmt;
use std::mem::size_of;
use std::ops::Range;

use rustfft::FftDirection;
use rustfft::num_complex::Complex;

use crate::kernels::numeric_common::scale_product;
use crate::kernels::par;
use crate::kernels::reuse::fft_plan;
use crate::kernels::reuse::fft_scratch_for;

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
    /// The input does not contain enough samples for the requested calculation.
    InsufficientSamples,
    /// The maximum lag must be less than the sample count.
    InvalidLag,
    /// The sample or prediction error variance is zero or numerically singular.
    DegenerateVariance,
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
            | Self::InsufficientSamples => formatter.write_str("insufficient samples"),
            | Self::InvalidLag => {
                formatter.write_str("maximum lag must be less than the sample count")
            },
            | Self::DegenerateVariance => {
                formatter.write_str("zero or numerically singular variance")
            },
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

const WELCH_MIN_PARALLEL_WORK: usize = 262_144;
const WELCH_WORK_PER_BLOCK: usize = 32_768;
const WELCH_MAX_PARTIAL_BYTES: usize = 8 * 1024 * 1024;
const WELCH_MAX_BLOCKS: usize = 16;

struct WelchFrames {
    weights: Vec<f64>,
    hop: usize,
    count: usize,
    inverse_rate_root: f64,
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "frame indices and bin counts are bounded by the input and window lengths"
)]
impl WelchFrames {
    fn new(
        sample_count: usize,
        window: &[f64],
        hop: usize,
        sample_rate: f64,
    ) -> Result<Self, SignalError> {
        if hop == 0 {
            return Err(SignalError::InvalidHop);
        }
        if !sample_rate.is_finite() || sample_rate <= 0.0 {
            return Err(SignalError::InvalidSampleRate);
        }
        let weights = normalized_window(window)?;
        if sample_count < weights.len() {
            return Err(SignalError::InsufficientSamples);
        }
        Ok(Self {
            count: (sample_count - weights.len()) / hop + 1,
            weights,
            hop,
            inverse_rate_root: sample_rate.sqrt().recip(),
        })
    }

    fn frame<'a>(
        &self,
        samples: &'a [f64],
        index: usize,
    ) -> &'a [f64] {
        let start = index * self.hop;
        &samples[start..start + self.weights.len()]
    }

    const fn bins(&self) -> usize {
        self.weights.len() / 2 + 1
    }

    fn covered_regions(&self) -> impl Iterator<Item = Range<usize>> {
        let length = self.weights.len();
        let (count, width) = if self.hop > length {
            (self.count, length)
        } else {
            (1, (self.count - 1) * self.hop + length)
        };
        (0..count).map(move |index| {
            let start = index * self.hop;
            start..start + width
        })
    }

    fn scale_bin(
        &self,
        value: Complex<f64>,
        bin: usize,
        inverse_rate_root: f64,
    ) -> Complex<f64> {
        let fold = if bin == 0 || (self.weights.len().is_multiple_of(2) && bin + 1 == self.bins()) {
            1.0
        } else {
            std::f64::consts::SQRT_2
        };
        value * (inverse_rate_root * fold)
    }

    fn average_bins<T: Send>(
        &self,
        calculate: impl Fn(Range<usize>) -> Result<Vec<T>, SignalError> + Sync,
        merge: impl Fn(&mut T, T, f64),
    ) -> Result<Vec<T>, SignalError> {
        let frame_len = self.weights.len();
        if self.count.saturating_mul(frame_len) < WELCH_MIN_PARALLEL_WORK {
            return calculate(0..self.count);
        }
        let bytes_per_block = self.bins().saturating_mul(size_of::<T>());
        let max_blocks =
            (WELCH_MAX_PARTIAL_BYTES / bytes_per_block.max(1)).clamp(1, WELCH_MAX_BLOCKS);
        let block_len = WELCH_WORK_PER_BLOCK
            .div_ceil(frame_len)
            .max(self.count.div_ceil(max_blocks));
        let block_count = self.count.div_ceil(block_len);
        if block_count == 1 {
            return calculate(0..self.count);
        }

        // Set block boundaries from the problem size so one worker and many
        // workers perform reductions in exactly the same order.
        let mut blocks = par::map_indexed(block_count, 2, |block| {
            let start = block * block_len;
            let end = start + block_len.min(self.count - start);
            calculate(start..end).map(|values| (end - start, values))
        })
        .into_iter();
        let Some(first) = blocks.next() else {
            return Err(SignalError::InsufficientSamples);
        };
        let (mut seen, mut average) = first?;
        for block in blocks {
            let (count, values) = block?;
            seen += count;
            let weight = count as f64 / seen as f64;
            for (mean, value) in average.iter_mut().zip(values) {
                merge(mean, value, weight);
            }
        }
        Ok(average)
    }
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "normalized windows keep finite input products bounded"
)]
fn window_frame(
    output: &mut [Complex<f64>],
    samples: &[f64],
    weights: &[f64],
    normalize: impl Fn(f64) -> Result<f64, SignalError>,
) -> Result<(), SignalError> {
    for ((value, &sample), &weight) in output.iter_mut().zip(samples).zip(weights) {
        *value = Complex::new(normalize(sample)? * weight, 0.0);
    }
    Ok(())
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the spectral formulas use checked finite estimates and known sample counts"
)]
pub(super) fn welch_psd(
    samples: &[f64],
    window: &[f64],
    hop: usize,
    sample_rate: f64,
) -> Result<PowerSpectrum, SignalError> {
    let frames = WelchFrames::new(samples.len(), window, hop, sample_rate)?;
    let plan = fft_plan(window.len(), FftDirection::Forward);
    let scratch_len = plan.get_inplace_scratch_len();
    let density = frames.average_bins(
        |range| {
            let mut spectrum = fft_scratch_for(window.len());
            spectrum.resize(window.len(), Complex::ZERO);
            let mut scratch = fft_scratch_for(scratch_len);
            scratch.resize(scratch_len, Complex::ZERO);
            let mut density = vec![0.0; frames.bins()];
            for (index, frame) in range.enumerate() {
                window_frame(
                    &mut spectrum,
                    frames.frame(samples, frame),
                    &frames.weights,
                    |value| {
                        value
                            .is_finite()
                            .then_some(value)
                            .ok_or(SignalError::NonFiniteInput)
                    },
                )?;
                plan.process_with_scratch(&mut spectrum, &mut scratch);
                let weight = 1.0 / (index + 1) as f64;
                for (bin, (mean, &value)) in density.iter_mut().zip(&*spectrum).enumerate() {
                    let value = frames.scale_bin(value, bin, frames.inverse_rate_root);
                    let mut power = value.norm_sqr();
                    if power < f64::MIN_POSITIVE && value != Complex::ZERO {
                        // Individual squares can round to zero before their
                        // sum reaches the representable subnormal range.
                        power = value.norm().powi(2);
                    }
                    if !power.is_finite() {
                        return Err(SignalError::NonFiniteResult);
                    }
                    *mean += (power - *mean) * weight;
                }
            }
            Ok(density)
        },
        |mean, power, weight| *mean += (power - *mean) * weight,
    )?;

    let frequencies = frequencies(window.len(), sample_rate);
    Ok(PowerSpectrum { frequencies, density })
}

#[derive(Clone, Copy, Default)]
struct CrossBin {
    xx: f64,
    yy: f64,
    xy: Complex<f64>,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "nonnegative powers and weighted complex cross products stay representable"
)]
impl CrossBin {
    fn merge(
        &mut self,
        other: Self,
        weight: f64,
    ) {
        self.xx += (other.xx - self.xx) * weight;
        self.yy += (other.yy - self.yy) * weight;
        // Scale before subtracting because cross products can have either sign.
        self.xy += other.xy * weight - self.xy * weight;
    }
}

/// Finds independent binary scales after an unscaled spectral product loses
/// range. Only samples used by complete frames inform the scales.
#[allow(
    clippy::indexing_slicing,
    reason = "the frame plan confines every region to both inputs"
)]
#[allow(
    clippy::single_call_fn,
    reason = "choose binary scales only for the recovery path"
)]
fn cross_scales(
    x: &[f64],
    y: &[f64],
    frames: &WelchFrames,
) -> Result<[f64; 2], SignalError> {
    let peaks = frames
        .covered_regions()
        .flat_map(|region| x[region.clone()].iter().zip(&y[region]))
        .try_fold([0.0_f64; 2], |[first, second], (&x, &y)| {
            (x.is_finite() && y.is_finite()).then(|| [first.max(x.abs()), second.max(y.abs())])
        })
        .ok_or(SignalError::NonFiniteInput)?;
    Ok(peaks.map(|peak| {
        if peak <= 0.0 {
            1.0
        } else {
            // A binary scale keeps normalization exact whenever its result
            // remains normal, with a finite inverse for subnormal peaks.
            f64::from_bits(peak.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE)
        }
    }))
}

// Leave room in the exponent range for the FFT sum, the square-root sample
// rate factor and the squared spectral products.
const WELCH_FAST_MIN_AMPLITUDE: f64 = 1e-90;
const WELCH_FAST_MAX_AMPLITUDE: f64 = 1e90;

#[allow(
    clippy::needless_bitwise_bool,
    reason = "keeps the contiguous sample check branch-free"
)]
fn ordinary_cross_samples(samples: &[f64]) -> bool {
    const MIN_BITS: u64 = WELCH_FAST_MIN_AMPLITUDE.to_bits();
    const MAX_BITS: u64 = WELCH_FAST_MAX_AMPLITUDE.to_bits();
    const AMPLITUDE_BITS: u64 = u64::MAX >> 1;
    // Positive float encodings are ordered by magnitude. Use a reduction
    // over the contiguous input so overlapping frames do not repeat checks.
    samples.iter().fold(true, |valid, sample| {
        let magnitude = sample.to_bits() & AMPLITUDE_BITS;
        valid
            & ((magnitude.wrapping_sub(MIN_BITS) <= MAX_BITS.saturating_sub(MIN_BITS))
                | (magnitude == 0))
    })
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
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
    let frames = WelchFrames::new(x.len(), window, hop, sample_rate)?;
    let plan = fft_plan(window.len(), FftDirection::Forward);
    let scratch_len = plan.get_inplace_scratch_len();
    let estimate =
        |inverse_scales: Option<[f64; 2]>| {
            let rate_root = if inverse_scales.is_some() {
                1.0
            } else {
                frames.inverse_rate_root
            };
            frames.average_bins(
                |range| {
                    let mut first = fft_scratch_for(window.len());
                    first.resize(window.len(), Complex::ZERO);
                    let mut second = fft_scratch_for(window.len());
                    second.resize(window.len(), Complex::ZERO);
                    let mut scratch = fft_scratch_for(scratch_len);
                    scratch.resize(scratch_len, Complex::ZERO);

                    let mut statistics = vec![CrossBin::default(); frames.bins()];
                    for (index, frame) in range.enumerate() {
                        let (x, y) = (frames.frame(x, frame), frames.frame(y, frame));
                        if let Some([x_inverse, y_inverse]) = inverse_scales {
                            window_frame(&mut first, x, &frames.weights, |value| {
                                Ok(value * x_inverse)
                            })?;
                            window_frame(&mut second, y, &frames.weights, |value| {
                                Ok(value * y_inverse)
                            })?;
                        } else {
                            window_frame(&mut first, x, &frames.weights, Ok)?;
                            window_frame(&mut second, y, &frames.weights, Ok)?;
                        }
                        plan.process_with_scratch(&mut first, &mut scratch);
                        plan.process_with_scratch(&mut second, &mut scratch);
                        let weight = 1.0 / (index + 1) as f64;
                        for (bin, (mean, (&raw_x, &raw_y))) in statistics
                            .iter_mut()
                            .zip(first.iter().zip(&*second))
                            .enumerate()
                        {
                            let a = frames.scale_bin(raw_x, bin, rate_root);
                            let b = frames.scale_bin(raw_y, bin, rate_root);
                            let xx = a.norm_sqr();
                            let yy = b.norm_sqr();
                            let xy = a.conj() * b;
                            if !xx.is_finite()
                                || !yy.is_finite()
                                || !xy.re.is_finite()
                                || !xy.im.is_finite()
                            {
                                return Err(SignalError::NonFiniteResult);
                            }
                            mean.merge(CrossBin { xx, yy, xy }, weight);
                        }
                    }
                    Ok(statistics)
                },
                CrossBin::merge,
            )
        };
    let ordinary = (WELCH_FAST_MIN_AMPLITUDE..=WELCH_FAST_MAX_AMPLITUDE).contains(&sample_rate)
        && frames.covered_regions().all(|region| {
            ordinary_cross_samples(&x[region.clone()]) && ordinary_cross_samples(&y[region])
        });
    let direct = if ordinary {
        estimate(None)
    } else {
        Err(SignalError::NonFiniteResult)
    };
    let (statistics, scales) = match direct {
        | Err(SignalError::NonFiniteResult) => {
            let scales = cross_scales(x, y, &frames)?;
            (estimate(Some(scales.map(f64::recip)))?, Some(scales))
        },
        | result => (result?, None),
    };
    let (mut density, coherence): (Vec<Complex<f64>>, Vec<f64>) = statistics
        .into_iter()
        .map(|CrossBin { xx, yy, xy }| {
            let magnitude = if xx > 0.0 && yy > 0.0 {
                (xy.norm() / xx.sqrt() / yy.sqrt()).min(1.0)
            } else {
                0.0
            };
            (xy, magnitude * magnitude)
        })
        .unzip();
    if let Some(scales) = scales {
        let restore = |value| {
            scale_product(value, &scales, &[sample_rate]).map_err(|_| SignalError::NonFiniteResult)
        };
        for value in &mut density {
            *value = Complex::new(restore(value.re)?, restore(value.im)?);
        }
    }
    Ok(CrossSpectrum {
        frequencies: frequencies(window.len(), sample_rate),
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
    let mut time = fft_scratch_for(frame_len);
    time.resize(frame_len, Complex::ZERO);
    let scratch_len = fft.get_inplace_scratch_len();
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex::ZERO);
    let inverse_len = 1.0 / frame_len as f64;
    let tiny_bound = f64::MIN_POSITIVE * frame_len as f64;
    for (index, frame) in frames.iter().enumerate() {
        let offset = index * hop;
        let tiny_frame = frame
            .iter()
            .all(|value| value.re.abs() <= tiny_bound && value.im.abs() <= tiny_bound);
        if tiny_frame {
            // Tiny coefficients can add to a representable subnormal even
            // when scaling each coefficient would round it to zero.
            time.copy_from_slice(frame);
        } else {
            // Normalizing before the inverse FFT leaves headroom for its
            // intermediate sums.
            for (value, &source) in time.iter_mut().zip(frame) {
                *value = source * inverse_len;
            }
        }
        fft.process_with_scratch(&mut time, &mut scratch);
        if tiny_frame {
            for value in &mut *time {
                *value *= inverse_len;
            }
        }
        for ((overlap, sample), &weight) in
            overlapping.iter_mut().skip(offset).zip(&*time).zip(window)
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
    evaluate_polyphase(samples, taps, up, down, 0, len)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "each output is confined to the support of the two input slices"
)]
pub(super) fn evaluate_polyphase(
    samples: &[f64],
    taps: &[f64],
    up: usize,
    down: usize,
    delay: usize,
    len: usize,
) -> Result<Vec<f64>, SignalError> {
    let mut output = zeroed_output(len)?;
    par::for_each_mut(&mut output, 16_384, |index, result| {
        let output_time = index as u128 * down as u128 + delay as u128;
        let head = output_time / up as u128;
        let phase = (output_time % up as u128) as usize;
        let first = head.saturating_sub(samples.len() as u128 - 1);
        let first_tap = phase as u128 + first * up as u128;
        if first_tap < taps.len() as u128 {
            let last_sample = (head - first) as usize;
            let filter = taps[first_tap as usize..].iter().step_by(up);
            let inputs = samples[..=last_sample].iter().rev();
            *result = filter.zip(inputs).map(|(tap, sample)| tap * sample).sum();
        }
    });
    if output.iter().any(|value| !value.is_finite()) {
        return Err(SignalError::NonFiniteResult);
    }
    Ok(output)
}
