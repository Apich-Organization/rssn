//! # FFT and fast convolution
//!
//! Arbitrary-length complex FFT, real convolution and correlation, and exact
//! integer polynomial products.

use num_bigint::BigInt;
use num_complex::Complex64;
use rustfft::FftDirection;

use super::reuse::fft_plan;
use super::reuse::fft_scratch;

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
    let padded = |x: &[f64]| {
        let mut result = Vec::with_capacity(period);
        result.extend(x.iter().map(|&v| Complex64::new(v, 0.0)));
        result.resize(period, Complex64::default());
        result
    };
    let mut left = padded(a);
    let mut right = padded(b);
    let forward = fft_plan(period, FftDirection::Forward);
    let inverse = fft_plan(period, FftDirection::Inverse);
    let scratch_len = forward
        .get_inplace_scratch_len()
        .max(inverse.get_inplace_scratch_len());
    let mut scratch = fft_scratch();
    scratch.resize(scratch_len, Complex64::default());
    forward.process_with_scratch(&mut left, &mut scratch);
    forward.process_with_scratch(&mut right, &mut scratch);
    for (x, y) in left.iter_mut().zip(right) {
        *x *= y;
    }
    inverse.process_with_scratch(&mut left, &mut scratch);
    let normalization = 1.0 / period as f64;
    left.into_iter()
        .take(used)
        .map(|v| v.re * normalization)
        .collect()
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
    let n = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 16 {
        return convolve_direct(a, b);
    }
    convolution_fft(a, b, n.next_power_of_two(), n)
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
    convolution_fft(a, b, period, period)
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
