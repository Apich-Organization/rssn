//! Signal processing: filters, windows and spectra.
//!
//! # Signal processing kernels
//!
//! Complex Fourier transforms, time-frequency analysis, windows and digital
//! filters. Use [`welch_psd`] to estimate one-sided power spectral density
//! and [`inverse_stft`] for overlap-add reconstruction. [`upfirdn`] applies
//! a caller-supplied FIR filter while changing the sample rate.
//!
//! The Fourier transforms reuse plans and work buffers by transform length.
//! Direct [`convolve`] and [`cross_correlation`] cost `O(n m)`; for FFT-based
//! convolution see [`crate::kernels::fftconv`]. The built-in Hann, Hamming,
//! Blackman and Kaiser windows are symmetric. Filtered frames use the exact
//! window or filter taps supplied by the caller.

use std::f64::consts::PI;

use num_complex::Complex64;
use rustfft::FftDirection;
use rustfft::num_complex::Complex;

use super::fftconv::convolve_direct;
use super::par;
use super::reuse::fft_plan;
use super::reuse::fft_scratch_for;

mod filters;
mod framed;
mod spectral;
mod time_series;

pub use filters::FilterError;
pub use framed::CrossSpectrum;
pub use framed::PowerSpectrum;
pub use framed::SignalError;

/// Estimates the one-sided power spectral density of real samples by Welch's method.
///
/// Averages the periodograms of complete windowed frames, whose starts are
/// `hop` samples apart. Pass the window coefficients explicitly, for example
/// with [`hann_window`]. The transform uses the window length and
/// performs no detrending. Frequencies have the same unit as `sample_rate`.
/// Only bins with a negative-frequency partner are doubled. Multiplying each
/// bin by the frequency spacing and summing gives the mean window-weighted
/// frame power. Trailing samples that do not complete a frame are ignored.
/// Large batches average fixed blocks in order so the result is independent
/// of the worker count.
///
/// # Errors
///
/// Returns an error for an invalid window, hop or sample rate, insufficient
/// samples, non-finite data in a complete frame or non-finite computed power.
pub fn welch_psd(
    samples: &[f64],
    window: &[f64],
    hop: usize,
    sample_rate: f64,
) -> Result<PowerSpectrum, SignalError> {
    framed::welch_psd(samples, window, hop, sample_rate)
}

/// Estimates one-sided cross-spectral density and magnitude-squared coherence.
///
/// Applies the same complete-frame Welch convention as [`welch_psd`] to
/// both inputs. For the Fourier bins `X` and `Y`, it averages `conj(X) * Y`
/// across frames to obtain the cross-spectral density `Pxy`, and computes
/// `abs(Pxy)^2 / (Pxx * Pyy)`. Auto-spectra and cross-spectra share the
/// same window energy and one-sided scaling.
///
/// # Errors
/// Returns an error on unequal lengths or on the invalid input domains
/// described by [`welch_psd`].
pub fn welch_cross_spectrum(
    x: &[f64],
    y: &[f64],
    window: &[f64],
    hop: usize,
    sample_rate: f64,
) -> Result<CrossSpectrum, SignalError> {
    framed::welch_cross_spectrum(x, y, window, hop, sample_rate)
}

/// Computes biased sample autocovariances at lags from zero through `max_lag`.
///
/// For lag `k`, divides the sum of
/// `(samples[t] - mean) * (samples[t + k] - mean)` by the full sample count.
/// Lag zero is the population variance of the input. A constant series gives
/// zero at all lags. Computes short lag ranges directly and uses a zero-padded
/// FFT for sufficiently long ranges.
///
/// # Errors
///
/// Returns an error for empty or non-finite samples, a lag at least as large as
/// the sample count, or an unrepresentable autocovariance or work extent.
pub fn autocovariance(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    time_series::autocovariance(samples, max_lag)
}

/// Computes sample autocorrelations at lags from zero through `max_lag`.
///
/// Divides each biased sample autocovariance by the lag-zero
/// autocovariance. The lag-zero result is one. Input scaling cancels before
/// the division, so the correlation remains computable when the variance
/// itself would overflow or underflow.
///
/// # Errors
///
/// Returns an error for empty or non-finite samples, a lag at least as large as
/// the sample count, a constant series, or an unrepresentable work extent.
pub fn autocorrelation(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    time_series::autocorrelation(samples, max_lag)
}

/// Computes Yule–Walker partial autocorrelations at lags zero through `max_lag`.
///
/// Uses the Durbin–Levinson recursion on the biased sample autocorrelations.
/// For lag `k > 0`, the result is the final coefficient of the order-`k`
/// autoregression. The lag-zero result is one. After calculating the
/// autocorrelations, the recursion takes `O(max_lag²)` time and `O(max_lag)`
/// space.
///
/// # Errors
///
/// Returns the same input errors as [autocorrelation] and also rejects an
/// autoregression with numerically singular prediction-error variance.
pub fn partial_autocorrelation(
    samples: &[f64],
    max_lag: usize,
) -> Result<Vec<f64>, SignalError> {
    time_series::partial_autocorrelation(samples, max_lag)
}

/// Reconstructs a real signal from full [`stft`] spectra.
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
pub fn inverse_stft(
    frames: &[Vec<Complex64>],
    window: &[f64],
    hop: usize,
) -> Result<Vec<f64>, SignalError> {
    framed::inverse_stft(frames, window, hop)
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
pub fn upfirdn(
    samples: &[f64],
    taps: &[f64],
    up: usize,
    down: usize,
) -> Result<Vec<f64>, SignalError> {
    framed::upfirdn(samples, taps, up, down)
}

/// Forms the complex analytic signal of real periodic samples.
///
/// The real part is the original signal and the imaginary part is its
/// discrete Hilbert transform. Positive-frequency bins are doubled and
/// negative bins removed. DC and the Nyquist bin of an even-length signal
/// are preserved. An empty input gives an empty result.
///
/// # Errors
/// Returns an error for non-finite samples or an unrepresentable transform.
pub fn analytic_signal(samples: &[f64]) -> Result<Vec<Complex64>, SignalError> {
    spectral::analytic_signal(samples)
}

/// Returns the discrete Hilbert transform of real periodic samples.
///
/// A positive-frequency cosine becomes a sine, and DC and Nyquist components
/// contribute zero. See [`analytic_signal`] for the boundary convention.
///
/// # Errors
/// Returns an error for non-finite samples or an unrepresentable transform.
pub fn hilbert(samples: &[f64]) -> Result<Vec<f64>, SignalError> {
    spectral::hilbert(samples)
}

/// Returns the amplitude envelope `hypot(x, hilbert(x))` of real samples.
///
/// # Errors
/// Returns an error for non-finite samples or an unrepresentable envelope.
pub fn analytic_envelope(samples: &[f64]) -> Result<Vec<f64>, SignalError> {
    spectral::analytic_envelope(samples)
}

/// Computes arbitrary equally spaced Fourier samples using a chirp convolution.
///
/// Bin `k` is `sum_n x[n] exp(-2 pi i (start + k * step) n)` with `start`
/// and `step` measured in cycles per input sample. Frequencies are periodic
/// modulo one. This is the unit-circle chirp z-transform. It can resolve a
/// sub-band without evaluating the rest of the spectrum. A zero step repeats
/// one frequency. Zero bins give an empty result and an empty input gives
/// `bins` zeros. Time and storage scale with the input and output extents,
/// using an FFT of length `next_power_of_two(x.len() + bins - 1)`.
/// Repeated moderate transforms with the same input length, bin count, and
/// step can share the transformed chirp, even when the start varies.
///
/// # Errors
/// Returns an error for non-finite frequencies or input, unrepresentable
/// intermediate arithmetic, or an output or work extent that is too large.
pub fn zoom_fft(
    samples: &[Complex64],
    start: f64,
    step: f64,
    bins: usize,
) -> Result<Vec<Complex64>, SignalError> {
    spectral::zoom_fft(samples, start, step, bins)
}

/// Computes the one-dimensional discrete Fourier Transform.
///
/// This function uses the `rustfft` library to perform the FFT.
///
/// # Arguments
/// * `input` - A slice of complex samples. Its contents are preserved.
///
/// # Returns
/// A vector of complex numbers representing the FFT of the input.
///
/// # Example
/// ```rust
/// use rssn::kernels::signal::fft;
/// use rustfft::num_complex::Complex;
///
/// let mut input = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
/// ];
///
/// let output = fft(&mut input);
///
/// assert!((output[0].re - 4.0).abs() < 1e-9);
/// ```
pub fn fft(input: &mut [Complex<f64>]) -> Vec<Complex<f64>> {
    super::fft::fft(input)
}

/// Computes the one-dimensional inverse discrete Fourier Transform.
///
/// This is the positive-sign transform followed by the `1/n`
/// normalisation, so `ifft(fft(x)) == x` (up to round-off) for every length.
///
/// # Arguments
/// * `input` - A mutable slice of complex frequency-domain samples.
///
/// # Returns
/// A vector of complex numbers representing the inverse transform.
///
/// # Example
/// ```rust
/// use rssn::kernels::signal::fft;
/// use rssn::kernels::signal::ifft;
/// use rustfft::num_complex::Complex;
///
/// let mut input = vec![Complex::new(1.0, 0.5), Complex::new(-2.0, 0.0)];
///
/// let spectrum = fft(&mut input);
///
/// let back = ifft(&mut spectrum.clone());
///
/// assert!((back[0] - input[0]).norm() < 1e-12);
/// ```
#[must_use]
pub fn ifft(input: &mut [Complex<f64>]) -> Vec<Complex<f64>> {
    super::fft::ifft(input)
}

/// Computes the one-dimensional discrete linear convolution of two sequences.
///
/// Convolution is a mathematical operation that blends two functions to produce a third.
/// In signal processing, it is used to describe the effect of a linear time-invariant system
/// on an input signal.
///
/// # Arguments
/// * `a` - The first input sequence.
/// * `v` - The second input sequence.
///
/// # Returns
/// The discrete linear convolution of `a` and `v`.
///
/// # Example
/// ```rust
/// use rssn::kernels::signal::convolve;
///
/// let a = vec![1.0, 2.0, 3.0];
///
/// let v = vec![0.0, 1.0, 0.5];
///
/// let res = convolve(&a, &v);
///
/// assert_eq!(res, vec![0.0, 1.0, 2.5, 4.0, 1.5]);
/// ```
#[must_use]
pub fn convolve(
    a: &[f64],
    v: &[f64],
) -> Vec<f64> {
    convolve_direct(a, v)
}

/// Computes the discrete cross-correlation of two sequences.
///
/// Cross-correlation is a measure of similarity of two series as a function of the displacement of one relative to the other.
///
/// # Arguments
/// * `a` - The first input sequence.
/// * `v` - The second input sequence.
///
/// # Returns
/// The discrete cross-correlation of `a` and `v`.
///
/// # Example
/// ```rust
/// use rssn::kernels::signal::cross_correlation;
///
/// let a = vec![1.0, 2.0, 3.0];
///
/// let v = vec![0.0, 1.0, 0.5];
///
/// let res = cross_correlation(&a, &v);
/// // correlation(a, v)[k] = sum_i a[i] * v[i-k]
/// ```
#[must_use]
pub fn cross_correlation(
    a: &[f64],
    v: &[f64],
) -> Vec<f64> {
    if a.is_empty() || v.is_empty() {
        return Vec::new();
    }
    let v_rev: Vec<_> = v.iter().rev().copied().collect();
    convolve(a, &v_rev)
}

/// Generates a Hann window of length `n`.
///
/// # Arguments
/// * `n` - The number of points in the output window.
#[must_use]
pub fn hann_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }

    if n == 1 {
        return vec![1.0];
    }

    (0..n)
        .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos()))
        .collect()
}

/// Generates a Hamming window of length `n`.
///
/// # Arguments
/// * `n` - The number of points in the output window.
#[must_use]
pub fn hamming_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }

    if n == 1 {
        return vec![1.0];
    }

    (0..n)
        .map(|i| 0.46f64.mul_add(-(2.0 * PI * i as f64 / (n - 1) as f64).cos(), 0.54))
        .collect()
}

/// Generates a Blackman window of length `n`.
///
/// `w[i] = 0.42 - 0.5 cos(2 pi i / (n - 1)) + 0.08 cos(4 pi i / (n - 1))`,
/// symmetric with both endpoints at zero (in exact arithmetic).
///
/// # Arguments
/// * `n` - The number of points in the output window.
#[must_use]
pub fn blackman_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }

    if n == 1 {
        return vec![1.0];
    }

    (0..n)
        .map(|i| {
            let phase = 2.0 * PI * i as f64 / (n - 1) as f64;
            0.42 - 0.5 * phase.cos() + 0.08 * (2.0 * phase).cos()
        })
        .collect()
}

/// Generates a symmetric Kaiser window with shape parameter `beta`.
///
/// Either sign of a finite `beta` produces the same coefficients. For
/// `n > 1`, a non-finite shape parameter returns NaNs.
#[must_use]
pub fn kaiser_window(
    n: usize,
    beta: f64,
) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    if !beta.is_finite() {
        return vec![f64::NAN; n];
    }
    let beta = beta.abs();
    let denom = crate::kernels::special::bessel_i0e(beta);
    let half = 0.5 * (n.saturating_sub(1) as f64);
    (0..n)
        .map(|i| {
            let r = (i as f64 - half) / half;
            let r_squared = r * r;
            let shape = (1.0 - r_squared).max(0.0).sqrt();
            let scaled = crate::kernels::special::bessel_i0e(beta * shape) / denom;
            // beta * (shape - 1) loses the attenuation near the window center.
            let decay = beta * (r_squared / (1.0 + shape));
            (scaled.ln() - decay).exp()
        })
        .collect()
}

/// Orthonormal Type-II Discrete Cosine Transform (DCT-II) in `O(n log n)`.
#[must_use]
pub fn dct(x: &[f64]) -> Vec<f64> {
    cosine_transform(x, false)
}

/// Inverse of [`dct`] (orthonormal Type-III DCT).
#[must_use]
pub fn idct(x: &[f64]) -> Vec<f64> {
    inverse_cosine_transform(x.iter().copied(), x.len(), false)
}

/// Orthonormal Type-II Discrete Sine Transform (DST-II) in `O(n log n)`.
#[must_use]
pub fn dst(x: &[f64]) -> Vec<f64> {
    let mut output = cosine_transform(x, true);
    output.reverse();
    output
}

/// Inverse of [`dst`] (orthonormal Type-III DST).
#[must_use]
pub fn idst(x: &[f64]) -> Vec<f64> {
    inverse_cosine_transform(x.iter().rev().copied(), x.len(), true)
}

fn cosine_transform(
    input: &[f64],
    alternating: bool,
) -> Vec<f64> {
    let n = input.len();
    if n <= 1 {
        return input.to_vec();
    }
    // Arrange the even-indexed samples first, then the odd-indexed samples
    // in reverse. A phase rotation of this n-point FFT yields the DCT.
    // Negating odd samples and reversing the resulting bins yields the DST.
    let scale = (2.0 / n as f64).sqrt();
    let odd_scale = if alternating { -scale } else { scale };
    let mut reordered = Vec::with_capacity(n);
    reordered.extend(
        input
            .iter()
            .step_by(2)
            .map(|&value| Complex::new(value * scale, 0.0)),
    );
    reordered.extend(
        input
            .iter()
            .skip(1)
            .step_by(2)
            .rev()
            .map(|&value| Complex::new(value * odd_scale, 0.0)),
    );
    let phase = std::f64::consts::FRAC_PI_2 / n as f64;
    super::fft::fft_owned(reordered)
        .into_iter()
        .enumerate()
        .map(|(bin, value)| {
            if bin == 0 {
                value.re * std::f64::consts::FRAC_1_SQRT_2
            } else {
                (value * Complex::from_polar(1.0, -phase * bin as f64)).re
            }
        })
        .collect()
}

fn inverse_cosine_transform(
    input: impl DoubleEndedIterator<Item = f64> + Clone,
    n: usize,
    alternating: bool,
) -> Vec<f64> {
    if n <= 1 {
        return input.collect();
    }
    // Undo the phase rotation. The cosine coefficient of bin n-k contains
    // the negative imaginary component of the rotated bin k.
    let mut spectrum = Vec::with_capacity(n);
    let dc = input.clone().next().unwrap_or_default() / (n as f64).sqrt();
    spectrum.push(Complex::new(dc, 0.0));
    let scale = (0.5 / n as f64).sqrt();
    let phase = std::f64::consts::FRAC_PI_2 / n as f64;
    for (offset, (re, neg_im)) in input.clone().skip(1).zip(input.rev()).enumerate() {
        let angle = phase * (offset + 1) as f64;
        spectrum.push(Complex::new(re * scale, -neg_im * scale) * Complex::from_polar(1.0, angle));
    }
    let plan = fft_plan(n, FftDirection::Inverse);
    let scratch_len = plan.get_inplace_scratch_len();
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex::ZERO);
    plan.process_with_scratch(&mut spectrum, &mut scratch);
    let (even, odd) = spectrum.split_at(n.div_ceil(2));
    let mut result = Vec::with_capacity(n);
    let sign = if alternating { -1.0 } else { 1.0 };
    let mut odd = odd.iter().rev();
    for value in even {
        result.push(value.re);
        if let Some(value) = odd.next() {
            result.push(sign * value.re);
        }
    }
    result
}

/// Short-time Fourier transform of complete frames with starts `hop` samples apart.
///
/// Returns full unnormalised complex spectra. The number of bins per frame is
/// the window length. Incomplete trailing frames are discarded. Empty windows,
/// zero hops and signals shorter than a window yield an empty result.
#[must_use]
pub fn stft(
    x: &[f64],
    window: &[f64],
    hop: usize,
) -> Vec<Vec<Complex<f64>>> {
    let frame_len = window.len();
    if frame_len == 0 || hop == 0 || x.len() < frame_len {
        return Vec::new();
    }
    let plan = fft_plan(frame_len, FftDirection::Forward);
    let count = (x.len() - frame_len) / hop + 1;
    let mut output = vec![Vec::new(); count];
    let transform_frames = |block: usize, frames: &mut [Vec<Complex<f64>>]| {
        let scratch_len = plan.get_inplace_scratch_len();
        let mut scratch = fft_scratch_for(scratch_len);
        scratch.resize(scratch_len, Complex::ZERO);
        for (offset, frame) in frames.iter_mut().enumerate() {
            let start = (block + offset) * hop;
            *frame = x[start..start + frame_len]
                .iter()
                .zip(window)
                .map(|(sample, weight)| Complex::new(sample * weight, 0.0))
                .collect();
            plan.process_with_scratch(frame, &mut scratch);
        }
    };
    let block_size = (16_384 / frame_len).max(1);
    if count.saturating_mul(frame_len) >= 262_144 {
        par::for_chunks_mut(&mut output, block_size, 2, |block, frames| {
            transform_frames(block * block_size, frames);
        });
    } else {
        transform_frames(0, &mut output);
    }
    output
}

/// Power spectrogram `|STFT|^2` of `x` with `window` and `hop`.
#[must_use]
pub fn spectrogram(
    x: &[f64],
    window: &[f64],
    hop: usize,
) -> Vec<Vec<f64>> {
    stft(x, window, hop)
        .into_iter()
        .map(|frame| frame.into_iter().map(|c| c.norm_sqr()).collect())
        .collect()
}

/// Evaluates a digital transfer function at frequencies in cycles per sample.
///
/// Uses the coefficient convention of [`lfilter`]. Frequencies are periodic
/// modulo one, with Nyquist at `0.5`. Returns a complex response per frequency.
///
/// # Errors
///
/// Errors for empty coefficients, non-finite inputs or a singular or
/// unrepresentable response at a requested frequency.
pub fn frequency_response(
    b: &[f64],
    a: &[f64],
    frequencies: &[f64],
) -> Result<Vec<Complex64>, FilterError> {
    filters::frequency_response(b, a, frequencies)
}

/// Designs a symmetric, unit-DC-gain low-pass FIR using a Kaiser window.
///
/// The cutoff is normalized to Nyquist and must lie strictly between zero
/// and one. A larger absolute `beta` improves sidelobe attenuation at the
/// expense of transition width.
///
/// # Errors
///
/// Errors for invalid parameters or unrepresentable coefficients.
pub fn firwin_lowpass(
    num_taps: usize,
    cutoff: f64,
    beta: f64,
) -> Result<Vec<f64>, FilterError> {
    filters::firwin(num_taps, 0.0, cutoff, beta)
}

/// Designs a symmetric band-pass FIR with unit gain at the passband center.
///
/// The cutoffs are ordered and normalized to Nyquist, with
/// `0 < low < high < 1`. The window is Kaiser, parameterized by `beta`.
///
/// # Errors
///
/// Errors for invalid parameters or unrepresentable coefficients.
pub fn firwin_bandpass(
    num_taps: usize,
    low: f64,
    high: f64,
    beta: f64,
) -> Result<Vec<f64>, FilterError> {
    if low <= 0.0 {
        return Err(FilterError::InvalidParameters);
    }
    filters::firwin(num_taps, low, high, beta)
}

/// Applies a zero-state cascade of second-order digital sections.
///
/// Each row of `sections` is `[b0, b1, b2, a0, a1, a2]`. The coefficients
/// are normalized by `a0` and evaluated with transposed direct form II.
/// An empty cascade is the identity.
///
/// # Errors
///
/// Errors for a zero leading denominator, non-finite input, or
/// unrepresentable state or output.
pub fn sos_filter(
    samples: &[f64],
    sections: &[[f64; 6]],
) -> Result<Vec<f64>, FilterError> {
    filters::sos_filter(samples, sections)
}

/// Resamples a real signal by a rational factor using a centered polyphase FIR.
///
/// Reduces the factors by their greatest common divisor, then designs a
/// Kaiser low-pass filter with 10 taps on each side per rate phase
/// (`beta = 5`). The input is zero extended at the boundaries. Output
/// time zero aligns with input time zero, and the output contains
/// `ceil(samples.len() * up / down)` samples.
///
/// # Errors
///
/// Errors for a zero factor, non-finite input or output, or unallocatable
/// filter or output dimensions.
pub fn resample_poly(
    samples: &[f64],
    up: usize,
    down: usize,
) -> Result<Vec<f64>, FilterError> {
    filters::resample_poly(samples, up, down)
}

/// Linear FIR/IIR digital filter (`a[0] * y[n] = sum_k b[k]*x[n-k] - sum_{k>=1} a[k]*y[n-k]`)
/// evaluated via Direct Form II Transposed with pooled scratch state.
#[must_use]
pub fn lfilter(
    b: &[f64],
    a: &[f64],
    x: &[f64],
) -> Vec<f64> {
    lfilter_with_zi(b, a, x, None)
}

fn lfilter_with_zi(
    b: &[f64],
    a: &[f64],
    x: &[f64],
    zi: Option<&[f64]>,
) -> Vec<f64> {
    if x.is_empty() || b.is_empty() || a.is_empty() {
        return Vec::new();
    }
    let a0 = a.first().copied().unwrap_or(1.0);
    if a0 == 0.0 {
        return vec![f64::NAN; x.len()];
    }
    let order = b.len().max(a.len());
    let mut b_norm = super::reuse::f64_scratch(order);
    b_norm.clear();
    for i in 0..order {
        b_norm.push(b.get(i).copied().unwrap_or(0.0) / a0);
    }
    let mut a_norm = super::reuse::f64_scratch(order);
    a_norm.clear();
    for i in 0..order {
        a_norm.push(a.get(i).copied().unwrap_or(0.0) / a0);
    }
    let state_len = order.saturating_sub(1);
    let mut z = super::reuse::f64_scratch(state_len.max(1));
    z.clear();
    if let Some(init) = zi {
        for i in 0..state_len {
            z.push(init.get(i).copied().unwrap_or(0.0));
        }
    } else {
        z.resize(state_len, 0.0);
    }
    let b0 = b_norm[0];

    let mut y = vec![0.0; x.len()];
    if state_len == 0 {
        for (dst, &sample) in y.iter_mut().zip(x) {
            *dst = b0 * sample;
        }
        return y;
    }
    for (dst, &sample) in y.iter_mut().zip(x) {
        let yi = b0.mul_add(sample, z[0]);
        *dst = yi;
        for k in 0..state_len - 1 {
            z[k] = b_norm[k + 1].mul_add(sample, (-a_norm[k + 1]).mul_add(yi, z[k + 1]));
        }
        z[state_len - 1] = b_norm[state_len].mul_add(sample, -a_norm[state_len] * yi);
    }
    y
}

/// Steady-state initial conditions `zi` for Direct Form II Transposed `lfilter`
/// such that a constant unit step input produces zero startup transient.
fn lfilter_zi(
    b: &[f64],
    a: &[f64],
) -> Option<Vec<f64>> {
    let a0 = a.first().copied()?;
    if a0 == 0.0 {
        return None;
    }
    let order = b.len().max(a.len());
    let m = order.saturating_sub(1);
    if m == 0 {
        return Some(Vec::new());
    }
    let b_norm: Vec<f64> = (0..order)
        .map(|i| b.get(i).copied().unwrap_or(0.0) / a0)
        .collect();
    let a_norm: Vec<f64> = (0..order)
        .map(|i| a.get(i).copied().unwrap_or(0.0) / a0)
        .collect();
    // Solve (I - A_companion) * zi = b[1..] - a[1..] * b[0].
    let mut mat = super::dense::Mat::identity(m);
    for i in 0..m {
        mat.set(i, 0, mat.at(i, 0) + a_norm[i + 1]);
        if i + 1 < m {
            mat.set(i, i + 1, mat.at(i, i + 1) - 1.0);
        }
    }
    let rhs: Vec<f64> = (0..m)
        .map(|i| b_norm[i + 1] - a_norm[i + 1] * b_norm[0])
        .collect();
    let lu = super::dense::lu_factor(&mat).ok()?;
    let zi = lu.solve(&rhs);
    zi.iter().all(|v| v.is_finite()).then_some(zi)
}

/// Zero-phase forward-backward digital filter (`filtfilt`) using odd-reflection
/// edge padding (`padlen = 3 * (max(b.len(), a.len()) - 1)`) and steady-state
/// initial conditions.
///
/// Applies [`lfilter`] forward and then backward to eliminate phase distortion.
/// For signals shorter than or equal to the padding length, pads up to `x.len() - 1`.
#[must_use]
pub fn filtfilt(
    b: &[f64],
    a: &[f64],
    x: &[f64],
) -> Vec<f64> {
    if x.is_empty() || b.is_empty() || a.is_empty() {
        return Vec::new();
    }
    if x.len() == 1 {
        let fwd = lfilter(b, a, x);
        return lfilter(b, a, &fwd);
    }
    let ntaps = b.len().max(a.len());
    let padlen = (3 * ntaps.saturating_sub(1)).min(x.len() - 1);
    let mut ext = super::reuse::f64_scratch(x.len() + 2 * padlen);
    ext.clear();
    // Odd reflection at left boundary: 2 * x[0] - x[padlen..=1].
    let x0 = x[0];
    for i in (1..=padlen).rev() {
        ext.push(2.0 * x0 - x[i]);
    }
    ext.extend_from_slice(x);
    // Odd reflection at right boundary: 2 * x[n-1] - x[n-2..n-1-padlen].
    let xn = x[x.len() - 1];
    for i in 1..=padlen {
        ext.push(2.0 * xn - x[x.len() - 1 - i]);
    }
    let zi_unit = lfilter_zi(b, a);
    let zi_fwd: Option<Vec<f64>> = zi_unit
        .as_ref()
        .map(|z| z.iter().map(|&v| v * ext[0]).collect());
    let mut fwd = lfilter_with_zi(b, a, &ext, zi_fwd.as_deref());
    fwd.reverse();
    let zi_bwd: Option<Vec<f64>> = zi_unit
        .as_ref()
        .map(|z| z.iter().map(|&v| v * fwd[0]).collect());
    let mut bwd = lfilter_with_zi(b, a, &fwd, zi_bwd.as_deref());
    bwd.reverse();
    bwd[padlen..padlen + x.len()].to_vec()
}

fn poly_mul(
    p: &[f64],
    q: &[f64],
) -> Vec<f64> {
    if p.is_empty() || q.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0.0; p.len().saturating_add(q.len()).saturating_sub(1)];
    for (i, &pi) in p.iter().enumerate() {
        for (j, &qj) in q.iter().enumerate() {
            if let Some(slot) = out.get_mut(i.saturating_add(j)) {
                *slot += pi * qj;
            }
        }
    }
    out
}

/// Digital low-pass Butterworth filter `(b, a)` of order `order` with
/// normalized cutoff `cutoff` in `(0, 1)` (where `1.0` is Nyquist).
#[must_use]
pub fn butterworth_lowpass(
    order: usize,
    cutoff: f64,
) -> (Vec<f64>, Vec<f64>) {
    if order == 0 || !(cutoff > 0.0 && cutoff < 1.0) {
        return (Vec::new(), Vec::new());
    }
    let omega = (0.5 * PI * cutoff).tan();
    let mut b = vec![1.0];
    let mut a = vec![1.0];
    let pairs = order / 2;
    for k in 0..pairs {
        let angle = PI * (2.0 * (k as f64) + 1.0 + (order as f64)) / (2.0 * (order as f64));
        let sigma = omega * angle.cos();
        let om_k = omega * angle.sin();
        let norm_sq = sigma * sigma + om_k * om_k;
        // Bilinear transform of H(s) = norm_sq / (s^2 - 2*sigma*s + norm_sq)
        let d0 = 1.0 - 2.0 * sigma + norm_sq;
        let d1 = 2.0 * (norm_sq - 1.0);
        let d2 = 1.0 + 2.0 * sigma + norm_sq;
        b = poly_mul(&b, &[norm_sq, 2.0 * norm_sq, norm_sq]);
        a = poly_mul(&a, &[d0, d1, d2]);
    }
    if !order.is_multiple_of(2) {
        // Real pole at s = -omega
        let d0 = 1.0 + omega;
        let d1 = omega - 1.0;
        b = poly_mul(&b, &[omega, omega]);
        a = poly_mul(&a, &[d0, d1]);
    }
    let a0 = a.first().copied().unwrap_or(1.0);
    if a0 != 0.0 {
        for v in &mut b {
            *v /= a0;
        }
        for v in &mut a {
            *v /= a0;
        }
    }
    (b, a)
}

/// Savitzky-Golay smoothing filter coefficients for odd `window_len` and
/// polynomial degree `poly_order < window_len`.
#[must_use]
pub fn savgol_coeffs(
    window_len: usize,
    poly_order: usize,
) -> Vec<f64> {
    if window_len == 0 || window_len.is_multiple_of(2) || poly_order >= window_len {
        return Vec::new();
    }
    let half = (window_len / 2) as isize;
    let ncols = poly_order.saturating_add(1);
    let mut v = crate::kernels::dense::Mat::zeros(window_len, ncols);
    for (r, m) in (-half..=half).enumerate() {
        let mut val = 1.0;
        for c in 0..ncols {
            v.set(r, c, val);
            val *= m as f64 / (half.max(1) as f64);
        }
    }
    // The first row of the Moore-Penrose inverse evaluates the fitted
    // polynomial at the window center. Scaling the abscissas avoids huge
    // Vandermonde powers; SVD avoids squaring its condition number.
    let Ok(inverse) = crate::kernels::dense::pinv(&v) else {
        return Vec::new();
    };
    inverse.data.into_iter().take(window_len).collect()
}

/// Applies a Savitzky-Golay smoothing filter of odd `window_len` and degree
/// `poly_order` to `x` with clamped boundary padding.
#[must_use]
pub fn savgol_filter(
    x: &[f64],
    window_len: usize,
    poly_order: usize,
) -> Vec<f64> {
    let coeffs = savgol_coeffs(window_len, poly_order);
    if coeffs.is_empty() || x.is_empty() {
        return Vec::new();
    }
    let half = window_len / 2;
    let n = x.len();
    let mut out = vec![0.0; n];
    let left_end = half.min(n);
    for i in 0..left_end {
        let mut sum = 0.0;
        for (k, &ck) in coeffs.iter().enumerate() {
            let idx = (i as isize + k as isize - half as isize)
                .clamp(0, n.saturating_sub(1) as isize) as usize;
            sum += ck * x[idx];
        }
        out[i] = sum;
    }
    if n > 2 * half {
        for i in half..n - half {
            out[i] = super::simd::dot_f64(&coeffs, &x[i - half..i + half + 1]);
        }
    }
    let right_start = n.saturating_sub(half).max(left_end);
    for i in right_start..n {
        let mut sum = 0.0;
        for (k, &ck) in coeffs.iter().enumerate() {
            let idx = (i as isize + k as isize - half as isize)
                .clamp(0, n.saturating_sub(1) as isize) as usize;
            sum += ck * x[idx];
        }
        out[i] = sum;
    }
    out
}

#[cfg(test)]
mod filter_tests {
    use super::*;

    #[test]
    fn filtfilt_preserves_constant_and_symmetry() {
        // FIR unit-DC-gain filter preserves constants exactly.
        let b = [0.25, 0.5, 0.25];
        let a = [1.0];
        let constant = vec![3.5; 32];
        let out_const = filtfilt(&b, &a, &constant);
        assert_eq!(out_const.len(), constant.len());
        for &v in &out_const {
            assert!((v - 3.5).abs() < 1e-12);
        }

        // IIR unit-DC-gain filter also preserves constants via steady-state zi.
        let out_iir_const = filtfilt(&[0.2, 0.2, 0.2], &[1.0, -0.4], &constant);
        for &v in &out_iir_const {
            assert!((v - 3.5).abs() < 1e-11);
        }

        // Zero-phase property: symmetric pulse stays symmetric around its center.
        let mut pulse = vec![0.0; 61];
        pulse[30] = 1.0;
        let y = filtfilt(&[0.2, 0.2, 0.2], &[1.0, -0.4], &pulse);
        for i in 0..30 {
            assert!((y[i] - y[60 - i]).abs() < 1e-11);
        }
    }
}
