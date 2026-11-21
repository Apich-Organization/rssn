use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::fmt;

use num_complex::Complex64;

use super::SignalError;
use super::framed::evaluate_polyphase;
use super::kaiser_window;
use crate::kernels::par;

/// Errors from checked filter design, filtering and convolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterError {
    /// Filter coefficients or design parameters are invalid.
    InvalidParameters,
    /// The input, output or work dimensions violate the signal contract.
    Signal(SignalError),
}

impl From<SignalError> for FilterError {
    fn from(error: SignalError) -> Self {
        Self::Signal(error)
    }
}

impl fmt::Display for FilterError {
    fn fmt(
        &self,
        formatter: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::InvalidParameters => formatter.write_str("invalid filter parameters"),
            | Self::Signal(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for FilterError {}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "cutoffs and buffer lengths are checked before designing FIR taps"
)]
pub(super) fn firwin(
    num_taps: usize,
    low: f64,
    high: f64,
    beta: f64,
) -> Result<Vec<f64>, FilterError> {
    if num_taps == 0 || !(low >= 0.0 && low < high && high < 1.0 && beta.is_finite()) {
        return Err(FilterError::InvalidParameters);
    }
    let mut taps = Vec::new();
    taps.try_reserve_exact(num_taps)
        .map_err(|_| SignalError::OutputSize)?;
    let middle = (num_taps - 1) as f64 * 0.5;
    // The bandwidth factor cancels with the gain normalization. Omitting it
    // also avoids subtracting almost identical sinusoids for narrow bands.
    let reference = if low > 0.0 {
        low.midpoint(high)
    } else {
        0.0
    };
    let width = if low > 0.0 {
        (high - low) * 0.5
    } else {
        high
    };
    let mut gain = 0.0;
    for (index, window) in kaiser_window(num_taps, beta).into_iter().enumerate() {
        let offset = index as f64 - middle;
        let angle = PI * width * offset;
        let sinc = if angle.abs() <= f64::EPSILON {
            1.0
        } else {
            angle.sin() / angle
        };
        let carrier = (PI * reference * offset).cos();
        let tap = window * sinc * carrier;
        taps.push(tap);
        gain += tap * carrier;
    }
    if !gain.is_finite() || gain.abs() <= 0.0 {
        return Err(FilterError::InvalidParameters);
    }
    for tap in &mut taps {
        *tap /= gain;
        if !tap.is_finite() {
            return Err(SignalError::NonFiniteResult.into());
        }
    }
    Ok(taps)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "frequencies are reduced before computing unit circle phases"
)]
pub(super) fn frequency_response(
    b: &[f64],
    a: &[f64],
    frequencies: &[f64],
) -> Result<Vec<Complex64>, FilterError> {
    if a.is_empty() || b.is_empty() {
        return Err(FilterError::InvalidParameters);
    }
    if a.iter()
        .chain(b)
        .chain(frequencies)
        .any(|value| !value.is_finite())
    {
        return Err(SignalError::NonFiniteInput.into());
    }
    // The Horner cost scales with the number of coefficients. A short
    // frequency grid can still supply enough work for the active pool.
    let coefficient_count = a.len().saturating_add(b.len());
    let minimum_points = 65_536_usize.div_ceil(coefficient_count).clamp(8, 1_024);
    let result = par::map_slice(frequencies, minimum_points, |&frequency| {
        let delay = Complex64::from_polar(1.0, -TAU * (frequency % 1.0));
        let horner = |coefficients: &[f64]| {
            coefficients
                .iter()
                .rev()
                .fold(Complex64::ZERO, |sum, &coefficient| {
                    sum * delay + coefficient
                })
        };
        let numerator = horner(b);
        let denominator = horner(a);
        // Scale before dividing to avoid squaring an extremely large or
        // small denominator inside complex division.
        let scale = denominator.re.abs().max(denominator.im.abs());
        if scale <= 0.0 || !scale.is_finite() {
            return Complex64::new(f64::NAN, f64::NAN);
        }
        (numerator / scale) / (denominator / scale)
    });
    if result
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(SignalError::NonFiniteResult.into());
    }
    Ok(result)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "filter normalization and all second-order section states are checked"
)]
pub(super) fn sos_filter(
    x: &[f64],
    sections: &[[f64; 6]],
) -> Result<Vec<f64>, FilterError> {
    if x.iter()
        .chain(sections.iter().flatten())
        .any(|value| !value.is_finite())
    {
        return Err(SignalError::NonFiniteInput.into());
    }
    let mut filtered = x.to_vec();
    for &[b0, b1, b2, a0, a1, a2] in sections {
        if a0.abs() <= 0.0 {
            return Err(FilterError::InvalidParameters);
        }
        let [b0, b1, b2, a1, a2] = [b0, b1, b2, a1, a2].map(|value| value / a0);
        if [b0, b1, b2, a1, a2].iter().any(|value| !value.is_finite()) {
            return Err(SignalError::NonFiniteResult.into());
        }
        let (mut first, mut second) = (0.0, 0.0);
        for sample in &mut filtered {
            let input = *sample;
            let output = b0.mul_add(input, first);
            first = b1.mul_add(input, second) - a1 * output;
            second = b2 * input - a2 * output;
            *sample = output;
        }
        if !first.is_finite()
            || !second.is_finite()
            || filtered.iter().any(|value| !value.is_finite())
        {
            return Err(SignalError::NonFiniteResult.into());
        }
    }
    Ok(filtered)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the reduced rates are bounded before allocating and applying the FIR"
)]
pub(super) fn resample_poly(
    x: &[f64],
    up: usize,
    down: usize,
) -> Result<Vec<f64>, FilterError> {
    use num_integer::Integer;

    if up == 0 || down == 0 {
        return Err(SignalError::InvalidRate.into());
    }
    if x.iter().any(|value| !value.is_finite()) {
        return Err(SignalError::NonFiniteInput.into());
    }
    if x.is_empty() || up == down {
        return Ok(x.to_vec());
    }
    let common = up.gcd(&down);
    let up = up / common;
    let down = down / common;
    let len = usize::try_from((x.len() as u128 * up as u128).div_ceil(down as u128))
        .map_err(|_| SignalError::OutputSize)?;
    let phases = up.max(down);
    let delay = phases.checked_mul(10).ok_or(SignalError::OutputSize)?;
    let num_taps = delay
        .checked_mul(2)
        .and_then(|width| width.checked_add(1))
        .ok_or(SignalError::OutputSize)?;
    let mut taps = firwin(num_taps, 0.0, 1.0 / phases as f64, 5.0)?;
    for value in &mut taps {
        *value *= up as f64;
    }
    evaluate_polyphase(x, &taps, up, down, delay, len).map_err(FilterError::from)
}
