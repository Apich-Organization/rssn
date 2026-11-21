//! # FFT and fast convolution
//!
//! Arbitrary-length complex FFT, real convolution and correlation, and exact
//! integer polynomial products.

use std::sync::Arc;

use num_bigint::BigInt;
use num_complex::Complex64;
use rustfft::Fft;
use rustfft::FftDirection;

use super::exact_eq;
use super::numeric_common::scale_product;
use super::par;
use super::reuse::fft_plan;
use super::reuse::fft_scratch_for;
use super::signal::FilterError;
use super::signal::SignalError;

/// Complex number as `(re, im)`.
pub type Complex = (f64, f64);

/// Discrete Fourier transform of any length (`inverse` includes the `1/n` scaling).
#[must_use]
pub fn fft(
    x: &[Complex],
    inverse: bool,
) -> Vec<Complex> {
    if x.len() <= 1 {
        return x.to_vec();
    }
    let samples: Vec<_> = x.iter().map(|&(re, im)| Complex64::new(re, im)).collect();
    let transformed = if inverse {
        super::fft::ifft_owned(samples)
    } else {
        super::fft::fft_owned(samples)
    };
    transformed.into_iter().map(|z| (z.re, z.im)).collect()
}

fn convolution_fft(
    a: &[f64],
    b: &[f64],
    period: usize,
    used: usize,
) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return vec![0.0; used];
    }
    let mut spectrum = fft_scratch_for(period);
    spectrum.resize(period, Complex64::ZERO);
    for (slot, &sample) in spectrum.iter_mut().zip(a) {
        slot.re = sample;
    }
    for (slot, &sample) in spectrum.iter_mut().zip(b) {
        slot.im = sample;
    }
    let forward = fft_plan(period, FftDirection::Forward);
    let inverse = fft_plan(period, FftDirection::Inverse);
    let scratch_len = forward
        .get_inplace_scratch_len()
        .max(inverse.get_inplace_scratch_len());
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex64::ZERO);

    // Two real sequences share one forward FFT. Their conjugate mirror
    // pairs recover both individual transforms.
    forward.process_with_scratch(&mut spectrum, &mut scratch);
    let normalization = 1.0 / period as f64;
    let multiply_real_bin = |bin: &mut Complex64| {
        *bin = Complex64::new(bin.re * (bin.im * normalization), 0.0);
    };
    multiply_real_bin(&mut spectrum[0]);
    for index in 1..period.div_ceil(2) {
        let mirror = period - index;
        let a = spectrum[index] * 0.5;
        let b = spectrum[mirror].conj() * 0.5;
        let left = a + b;
        let right = Complex64::new(0.0, -normalization) * (a - b);
        let product = left * right;
        spectrum[index] = product;
        spectrum[mirror] = product.conj();
    }
    if period.is_multiple_of(2) {
        multiply_real_bin(&mut spectrum[period / 2]);
    }
    inverse.process_with_scratch(&mut spectrum, &mut scratch);
    spectrum.iter().take(used).map(|v| v.re).collect()
}

/// Linear convolution of two real sequences (`len = a.len() + b.len() - 1`).
#[must_use]
pub fn convolve(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let (signal, taps) = if a.len() >= b.len() {
        (a, b)
    } else {
        (b, a)
    };
    if taps.len() <= 16 {
        return convolve_direct(a, b);
    }

    // Packing both channels is economical at comparable magnitudes but can
    // lose the smaller input to round-off. The ordinary two-spectrum path
    // is retained as a fallback for signals at very different scales.
    let n = a.len() + b.len() - 1;
    if signal.len() / taps.len() >= 16
        && let Ok(output) = convolve_overlap_save(signal, taps)
    {
        return output;
    }
    convolution_fft_scaled(a, b, n.next_power_of_two(), n)
}

fn convolution_fft_scaled(
    a: &[f64],
    b: &[f64],
    period: usize,
    used: usize,
) -> Vec<f64> {
    let peak = |values: &[f64]| values.iter().fold(0.0_f64, |peak, val| peak.max(val.abs()));
    let left_peak = peak(a);
    let right_peak = peak(b);
    if left_peak == 0.0 || right_peak == 0.0 {
        // Preserve non-finite propagation for inputs containing NaN or infinity.
        return convolution_fft(a, b, period, used);
    }
    let ratio = left_peak.min(right_peak) / left_peak.max(right_peak);
    if ratio >= 0.125 {
        return convolution_fft(a, b, period, used);
    }

    let mut left = fft_scratch_for(period);
    let mut right = fft_scratch_for(period);
    for (buffer, input) in [(&mut left, a), (&mut right, b)] {
        buffer.extend(input.iter().map(|&value| Complex64::new(value, 0.0)));
        buffer.resize(period, Complex64::ZERO);
    }
    let forward = fft_plan(period, FftDirection::Forward);
    let inverse = fft_plan(period, FftDirection::Inverse);
    let scratch_len = forward
        .get_inplace_scratch_len()
        .max(inverse.get_inplace_scratch_len());
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex64::ZERO);
    forward.process_with_scratch(&mut left, &mut scratch);
    forward.process_with_scratch(&mut right, &mut scratch);
    let normalization = 1.0 / period as f64;
    for (value, response) in left.iter_mut().zip(right.iter()) {
        *value *= response * normalization;
    }
    inverse.process_with_scratch(&mut left, &mut scratch);
    left.iter().take(used).map(|value| value.re).collect()
}

pub(super) fn convolve_direct(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut output = vec![0.0; a.len() + b.len() - 1];
    for (offset, &x) in a.iter().enumerate() {
        for (coefficient, &y) in output.iter_mut().skip(offset).zip(b) {
            *coefficient += x * y;
        }
    }
    output
}

/// Circular convolution of two real sequences of period
/// `max(a.len(), b.len())`.
///
/// A shorter sequence is zero-padded to the period.
#[must_use]
pub fn circular_convolve(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    let period = a.len().max(b.len());
    convolution_fft_scaled(a, b, period, period)
}

/// Cross-correlation `c[k] = sum_i a[i + k - (len(b) - 1)] * b[i]`
/// (full mode, length `a.len() + b.len() - 1`).
#[must_use]
pub fn correlate(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let rb: Vec<f64> = b.iter().rev().copied().collect();
    convolve(a, &rb)
}

/// Full linear convolution using block FFTs and a reused filter spectrum.
///
/// Uses zero extension and independent overlap-save blocks. The result has
/// `a.len() + taps.len() - 1` samples, or zero samples if either input is empty.
/// This is especially useful for a recording much longer than its filter. Use
/// [`OverlapSavePlan`] when applying the same filter to several recordings.
///
/// # Errors
///
/// Errors on non-finite inputs or output, or an unallocatable extent.
pub fn convolve_overlap_save(
    a: &[f64],
    taps: &[f64],
) -> Result<Vec<f64>, FilterError> {
    if a.is_empty() || taps.is_empty() {
        return Ok(Vec::new());
    }
    OverlapSavePlan::new(taps)?.apply(a)
}

fn needs_fft_amplitude_scaling(magnitude: f64) -> bool {
    // This range leaves headroom for intermediate FFT sums and the spectral
    // product without scaling ordinary recordings and filters.
    magnitude > 0.0 && !(1e-130..=1e130).contains(&magnitude)
}

fn fft_amplitude_scale(peak: f64) -> f64 {
    if needs_fft_amplitude_scaling(peak) {
        // A power of two preserves the significand when scaling. Use the
        // smallest normal scale for subnormal peaks so division is safe.
        f64::from_bits(peak.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE)
    } else {
        1.0
    }
}

/// Reusable FFT filter for full linear convolution by overlap-save.
///
/// The immutable plan can be shared across recordings or embedding threads.
/// The FFT length is the next power of two at or above four times the number
/// of taps. Each block uses zero padding where it extends beyond the recording.
/// Extreme filter and block amplitudes are rescaled before transformation.
pub struct OverlapSavePlan {
    overlap: usize,
    spectrum: Vec<Complex64>,
    tap_scale: f64,
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
    scratch_len: usize,
}

impl OverlapSavePlan {
    /// Creates a plan for a real finite impulse response filter.
    ///
    /// # Errors
    ///
    /// Errors for an empty or non-finite filter, a non-finite transform or
    /// an unallocatable FFT extent.
    pub fn new(taps: &[f64]) -> Result<Self, FilterError> {
        if taps.is_empty() {
            return Err(FilterError::InvalidParameters);
        }
        let mut peak: f64 = 0.0;
        for &tap in taps {
            if !tap.is_finite() {
                return Err(SignalError::NonFiniteInput.into());
            }
            peak = peak.max(tap.abs());
        }
        let length = taps
            .len()
            .checked_mul(4)
            .and_then(usize::checked_next_power_of_two)
            .ok_or(SignalError::OutputSize)?;
        let mut spectrum = Vec::new();
        spectrum
            .try_reserve_exact(length)
            .map_err(|_| SignalError::OutputSize)?;
        let normalization = 1.0 / length as f64;
        let tap_scale = fft_amplitude_scale(peak);
        let gain = tap_scale.recip();
        spectrum.extend(
            taps.iter()
                .map(|&tap| Complex64::new((tap * gain) * normalization, 0.0)),
        );
        spectrum.resize(length, Complex64::ZERO);
        let forward = fft_plan(length, FftDirection::Forward);
        let inverse = fft_plan(length, FftDirection::Inverse);
        let scratch_len = forward
            .get_inplace_scratch_len()
            .max(inverse.get_inplace_scratch_len());
        let mut scratch = fft_scratch_for(scratch_len);
        scratch.resize(scratch_len, Complex64::ZERO);
        forward.process_with_scratch(&mut spectrum, &mut scratch);
        if spectrum
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(SignalError::NonFiniteResult.into());
        }
        Ok(Self {
            overlap: taps.len() - 1,
            spectrum,
            tap_scale,
            forward,
            inverse,
            scratch_len,
        })
    }

    /// Applies the filter with zero extension and returns the full convolution.
    ///
    /// # Errors
    ///
    /// Errors on non-finite inputs or output, or an unallocatable output extent.
    pub fn apply(
        &self,
        samples: &[f64],
    ) -> Result<Vec<f64>, FilterError> {
        if samples.is_empty() {
            return Ok(Vec::new());
        }
        let mut non_finite = false;
        let mut scale_samples = false;
        for &sample in samples {
            non_finite |= !sample.is_finite();
            scale_samples |= needs_fft_amplitude_scaling(sample.abs());
        }
        if non_finite {
            return Err(SignalError::NonFiniteInput.into());
        }
        let len = self.output_len(samples.len())?;
        let mut output = Vec::new();
        output
            .try_reserve_exact(len)
            .map_err(|_| SignalError::OutputSize)?;
        output.resize(len, 0.0);

        let n = self.spectrum.len();
        let valid = n - self.overlap;
        let batch = valid * (16_384 / valid).max(1);
        par::for_chunks_mut(&mut output, batch, 4, |index, result| {
            let mut buffer = fft_scratch_for(n);
            buffer.resize(n, Complex64::ZERO);
            let mut scratch = fft_scratch_for(self.scratch_len);
            scratch.resize(self.scratch_len, Complex64::ZERO);

            for (offset, segment) in result.chunks_mut(valid).enumerate() {
                let start = index * batch + offset * valid;
                let first = start.saturating_sub(self.overlap);
                let left_pad = self.overlap.saturating_sub(start);
                let count = (n - left_pad).min(samples.len() - first);
                let input = &samples[first..first + count];
                let sample_scale = if scale_samples {
                    let peak = input.iter().fold(0.0_f64, |peak, val| peak.max(val.abs()));
                    fft_amplitude_scale(peak)
                } else {
                    1.0
                };
                let (left, block) = buffer.split_at_mut(left_pad);
                let (block, right) = block.split_at_mut(count);
                left.fill(Complex64::ZERO);
                right.fill(Complex64::ZERO);
                if exact_eq(sample_scale, 1.0) {
                    for (value, &sample) in block.iter_mut().zip(input) {
                        *value = Complex64::new(sample, 0.0);
                    }
                } else {
                    let gain = sample_scale.recip();
                    for (value, &sample) in block.iter_mut().zip(input) {
                        *value = Complex64::new(sample * gain, 0.0);
                    }
                }
                self.forward.process_with_scratch(&mut buffer, &mut scratch);
                for (value, response) in buffer.iter_mut().zip(&self.spectrum) {
                    *value *= response;
                }
                self.inverse.process_with_scratch(&mut buffer, &mut scratch);
                let factor = self.tap_scale * sample_scale;
                let values = segment.iter_mut().zip(&buffer[self.overlap..]);
                if exact_eq(factor, 1.0) {
                    for (value, transformed) in values {
                        *value = transformed.re;
                    }
                } else if factor.is_normal() {
                    for (value, transformed) in values {
                        *value = transformed.re * factor;
                    }
                } else {
                    for (value, transformed) in values {
                        *value =
                            scale_product(transformed.re, &[self.tap_scale, sample_scale], &[])
                                .unwrap_or_else(|_| f64::INFINITY.copysign(transformed.re));
                    }
                }
            }
        });
        if output.iter().any(|value| !value.is_finite()) {
            return Err(SignalError::NonFiniteResult.into());
        }
        Ok(output)
    }

    /// Applies this prepared filter to several recordings in input order.
    ///
    /// Each result has the same length and zero extension as [`Self::apply`].
    /// Empty recordings give empty results. Independent recordings can run in
    /// parallel under the active execution policy.
    ///
    /// # Errors
    ///
    /// Returns the error for the first recording that fails [`Self::apply`].
    pub fn apply_batch(
        &self,
        recordings: &[&[f64]],
    ) -> Result<Vec<Vec<f64>>, FilterError> {
        par::map_slice(recordings, 4, |samples| self.apply(samples))
            .into_iter()
            .collect()
    }

    pub(crate) fn output_len(
        &self,
        input_len: usize,
    ) -> Result<usize, SignalError> {
        if input_len == 0 {
            return Ok(0);
        }
        input_len
            .checked_add(self.overlap)
            .ok_or(SignalError::OutputSize)
    }

    pub(crate) const fn retained_budget(&self) -> usize {
        // Charge separately for the owned spectrum and the retained pair of
        // power-of-two FFT plans, which can outlive their cache entries.
        let entries = self
            .spectrum
            .capacity()
            .saturating_add(self.spectrum.len().saturating_mul(4));
        size_of::<Self>().saturating_add(entries.saturating_mul(size_of::<Complex64>()))
    }
}

/// Linear autocorrelation at nonnegative lags, from zero through
/// `samples.len() - 1`.
///
/// Uses zero padding so the sums do not wrap. The zero-lag coefficient
/// is the sum of squares. An empty input gives an empty output.
#[must_use]
pub fn autocorrelate(samples: &[f64]) -> Vec<f64> {
    let Some(max_lag) = samples.len().checked_sub(1) else {
        return Vec::new();
    };
    autocorrelate_prefix(samples, max_lag).unwrap_or_default()
}

pub(super) fn autocorrelate_prefix(
    samples: &[f64],
    max_lag: usize,
) -> Option<Vec<f64>> {
    // The requested lags cannot see values wrapped around by the inverse FFT.
    let count = max_lag.checked_add(1)?;
    let length = samples
        .len()
        .checked_add(max_lag)?
        .checked_next_power_of_two()?;
    std::alloc::Layout::array::<Complex64>(length).ok()?;
    let mut spectrum = fft_scratch_for(length);
    spectrum.extend(samples.iter().map(|&sample| Complex64::new(sample, 0.0)));
    spectrum.resize(length, Complex64::ZERO);

    let forward = fft_plan(length, FftDirection::Forward);
    let inverse = fft_plan(length, FftDirection::Inverse);
    let scratch_len = forward
        .get_inplace_scratch_len()
        .max(inverse.get_inplace_scratch_len());
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex64::ZERO);

    forward.process_with_scratch(&mut spectrum, &mut scratch);
    let normalization = 1.0 / length as f64;
    for value in &mut *spectrum {
        *value = Complex64::new(
            value.re * (value.re * normalization) + value.im * (value.im * normalization),
            0.0,
        );
    }
    inverse.process_with_scratch(&mut spectrum, &mut scratch);
    Some(spectrum.iter().take(count).map(|value| value.re).collect())
}

/// Multiplies two integer polynomials (ascending coefficients) exactly.
///
/// Each output is the corresponding integer coefficient modulo `2^64`,
/// interpreted as an `i64`. Use [`poly_mul_i64_checked`] to detect overflow.
#[must_use]
pub fn poly_mul_i64(
    a: &[i64],
    b: &[i64],
) -> Vec<i64> {
    integer_convolution(a, b, false).unwrap_or_default()
}

/// Multiplies two integer polynomials (ascending coefficients) exactly.
/// Returns `None` if any output coefficient is outside the range of `i64`.
#[must_use]
pub fn poly_mul_i64_checked(
    a: &[i64],
    b: &[i64],
) -> Option<Vec<i64>> {
    integer_convolution(a, b, true)
}

fn integer_convolution(
    a: &[i64],
    b: &[i64],
    checked: bool,
) -> Option<Vec<i64>> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let len = a.len().checked_add(b.len())? - 1;
    let overlap = a.len().min(b.len());
    if !checked && overlap <= 16 {
        let mut out = vec![0_i64; len];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                out[i + j] = out[i + j].wrapping_add(x.wrapping_mul(y));
            }
        }
        return Some(out);
    }
    let max_a = a.iter().map(|x| x.unsigned_abs()).max()?;
    let max_b = b.iter().map(|x| x.unsigned_abs()).max()?;
    if max_a == 0 || max_b == 0 {
        return Some(vec![0; len]);
    }
    let magnitude_bits = (max_a.ilog2() + 1) + (max_b.ilog2() + 1) + (overlap.ilog2() + 1);
    if checked && overlap <= 16 && magnitude_bits < 128 {
        let mut exact = vec![0_i128; len];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                exact[i + j] += i128::from(x) * i128::from(y);
            }
        }
        return exact.into_iter().map(|v| i64::try_from(v).ok()).collect();
    }

    // Kronecker substitution. Choose a power-of-two base wide enough for
    // every signed convolution coefficient, multiply the packed integers,
    // then decode balanced digits without floating-point rounding.
    let width = (magnitude_bits + 1).div_ceil(8) as usize;
    let product = pack_coefficients(a, width)? * pack_coefficients(b, width)?;
    let mut bytes = product.to_signed_bytes_le();
    let extension = if product.sign() == num_bigint::Sign::Minus {
        0xff
    } else {
        0
    };
    let size = width.checked_mul(len)?;
    bytes
        .try_reserve_exact(size.saturating_sub(bytes.len()))
        .ok()?;
    bytes.resize(size, extension);

    let mut carry = 0_i64;
    bytes
        .chunks_exact(width)
        .map(|digit| {
            let sign = digit.last().is_some_and(|b| b & 0x80 != 0);
            // Restore the borrow from the previous negative coefficient.
            // With wrapping results only the low 64 bits are needed.
            let mut low = [if sign { 0xff } else { 0 }; 8];
            let take = width.min(low.len());
            low[..take].copy_from_slice(&digit[..take]);
            let word = i64::from_le_bytes(low);
            let value = if !checked {
                word.wrapping_add(carry)
            } else if width <= low.len() {
                word.checked_add(carry)?
            } else {
                let value = BigInt::from_signed_bytes_le(digit) + carry;
                i64::try_from(value).ok()?
            };
            carry = i64::from(sign);
            Some(value)
        })
        .collect()
}

fn pack_coefficients(
    values: &[i64],
    width: usize,
) -> Option<BigInt> {
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(values.len().checked_mul(width)?)
        .ok()?;
    let mut borrow = 0_i128;
    for &value in values {
        let adjusted = i128::from(value) + borrow;
        borrow = if adjusted < 0 { -1 } else { 0 };
        let digits = adjusted.to_le_bytes();
        bytes.extend_from_slice(&digits[..width.min(digits.len())]);
        bytes.extend(std::iter::repeat_n(
            borrow as u8,
            width.saturating_sub(digits.len()),
        ));
    }
    Some(BigInt::from_signed_bytes_le(&bytes))
}
