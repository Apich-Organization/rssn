use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::fmt;

use num_complex::Complex64;

use super::SignalError;
use super::framed::evaluate_polyphase;
use super::kaiser_window;
use crate::kernels::par;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::f64_scratch;

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

struct DirectForm {
    numerator: Scratch<'static, f64>,
    denominator: Scratch<'static, f64>,
}

impl DirectForm {
    fn new(
        b: &[f64],
        a: &[f64],
    ) -> Self {
        let width = b.len().max(a.len());
        let leading = a[0];
        let normalize = |values: &[f64]| {
            let mut normalized = f64_scratch(width);
            normalized.clear();
            normalized.extend(values.iter().map(|value| value / leading));
            normalized.resize(width, 0.0 / leading);
            normalized
        };
        Self {
            numerator: normalize(b),
            denominator: normalize(a),
        }
    }

    fn apply<'a>(
        &self,
        samples: impl Iterator<Item = &'a mut f64>,
        state: &mut [f64],
    ) {
        let b0 = self.numerator[0];
        let Some(last) = state.len().checked_sub(1) else {
            for sample in samples {
                *sample *= b0;
            }
            return;
        };
        for sample in samples {
            let input = *sample;
            let output = b0.mul_add(input, state[0]);
            *sample = output;
            for k in 0..last {
                state[k] = self.numerator[k + 1].mul_add(
                    input,
                    (-self.denominator[k + 1]).mul_add(output, state[k + 1]),
                );
            }
            state[last] =
                self.numerator[last + 1].mul_add(input, -self.denominator[last + 1] * output);
        }
    }

    fn steady_state(
        &self,
        state: &mut [f64],
    ) -> bool {
        if state.is_empty() {
            return true;
        }
        let width = self.numerator.len();
        let sum_a = par::sum_kahan_fixed(width, 256, usize::MAX, |i| self.denominator[i]);
        let sum_b = par::sum_kahan_fixed(width, 256, usize::MAX, |i| self.numerator[i]);
        let gain = sum_b / sum_a;
        // At unit input and output H(1), the transposed filter state is a
        // reverse cumulative sum of b[k] - H(1)*a[k]. This solves the
        // companion system in linear time without a dense matrix.
        let mut next = 0.0;
        for (k, slot) in state.iter_mut().enumerate().rev() {
            next += (-self.denominator[k + 1]).mul_add(gain, self.numerator[k + 1]);
            *slot = next;
        }
        gain.is_finite() && state.iter().all(|value| value.is_finite())
    }
}

pub(super) fn lfilter(
    b: &[f64],
    a: &[f64],
    x: &[f64],
) -> Vec<f64> {
    if b.is_empty() || a.is_empty() || x.is_empty() {
        return Vec::new();
    }
    if a[0] == 0.0 {
        return vec![f64::NAN; x.len()];
    }
    let filter = DirectForm::new(b, a);
    let mut result = x.to_vec();
    let state_len = filter.numerator.len().saturating_sub(1);
    let mut state = f64_scratch(state_len);
    state.clear();
    state.resize(state_len, 0.0);
    filter.apply(result.iter_mut(), &mut state);
    result
}

pub(super) fn filtfilt(
    b: &[f64],
    a: &[f64],
    x: &[f64],
    checked: bool,
) -> Result<Vec<f64>, FilterError> {
    if b.is_empty() || a.is_empty() || a[0] == 0.0 {
        return Err(FilterError::InvalidParameters);
    }
    if checked && b.iter().chain(a).chain(x).any(|value| !value.is_finite()) {
        return Err(SignalError::NonFiniteInput.into());
    }
    if x.is_empty() {
        return Ok(Vec::new());
    }
    let filter = DirectForm::new(b, a);
    if checked
        && filter
            .numerator
            .iter()
            .chain(filter.denominator.iter())
            .any(|value| !value.is_finite())
    {
        return Err(SignalError::NonFiniteResult.into());
    }
    let state_len = filter.numerator.len().saturating_sub(1);
    let mut steady = f64_scratch(state_len);
    steady.clear();
    steady.resize(state_len, 0.0);
    if !filter.steady_state(&mut steady) {
        if checked {
            return Err(FilterError::InvalidParameters);
        }
        steady.fill(0.0);
    }
    let pad = state_len.saturating_mul(3).min(x.len() - 1);
    let extended_len = pad
        .checked_mul(2)
        .and_then(|extra| extra.checked_add(x.len()))
        .ok_or(SignalError::OutputSize)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(extended_len)
        .map_err(|_| SignalError::OutputSize)?;
    let first = x[0];
    let last = x[x.len() - 1];
    result.extend(
        x[..=pad]
            .iter()
            .skip(1)
            .rev()
            .map(|&sample| first.mul_add(2.0, -sample)),
    );
    result.extend_from_slice(x);
    result.extend(
        x[x.len() - 1 - pad..x.len() - 1]
            .iter()
            .rev()
            .map(|&sample| last.mul_add(2.0, -sample)),
    );

    let mut state = f64_scratch(state_len);
    state.clear();
    let forward_scale = result[0];
    state.extend(steady.iter().map(|value| value * forward_scale));
    filter.apply(result.iter_mut(), &mut state);
    if checked
        && (result.iter().any(|value| !value.is_finite())
            || state.iter().any(|value| !value.is_finite()))
    {
        return Err(SignalError::NonFiniteResult.into());
    }
    let backward_scale = result[extended_len - 1];
    state
        .iter_mut()
        .zip(steady.iter())
        .for_each(|(slot, &value)| *slot = value * backward_scale);
    filter.apply(result.iter_mut().rev(), &mut state);
    if checked
        && (result.iter().any(|value| !value.is_finite())
            || state.iter().any(|value| !value.is_finite()))
    {
        return Err(SignalError::NonFiniteResult.into());
    }
    if pad != 0 {
        result.copy_within(pad..pad + x.len(), 0);
        result.truncate(x.len());
    }
    Ok(result)
}
