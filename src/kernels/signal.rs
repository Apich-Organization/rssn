//! # Signal processing kernels
//!
//! A small signal-processing toolbox built on `rustfft`: the complex forward
//! and inverse DFT, linear convolution, cross-correlation, and symmetric
//! Hann/Hamming/Blackman analysis windows.
//!
//! # Limitations and coverage
//!
//! * `fft` / `ifft` are complex transforms: `fft` is unnormalised and `ifft`
//!   divides by `n`. A fresh planner is built on every call, so repeated short
//!   transforms pay rustfft's planning cost; there is no plan reuse, batched,
//!   multi-dimensional or in-place transform here. For a dependency-free FFT
//!   with real-input packing see [`crate::kernels::transforms`].
//! * Convolution and cross-correlation are direct `O(n m)`; there is no
//!   FFT-based fast convolution in this module, no 2-D or circular convolution
//!   and no boundary-mode choice.
//! * The windows are symmetric (not periodic / DFT-even) Hann, Hamming and
//!   Blackman; there is no Kaiser, flat-top, Tukey or custom window, and the
//!   windows are provided for analysis, not for filter design.
//! * There is no filtering, resampling, spectral estimation, STFT or wavelet
//!   helper; the transform accuracy is rustfft's `O(eps log n)` round-off.

use std::f64::consts::PI;

use rustfft::FftPlanner;
use rustfft::num_complex::Complex;

/// Computes the one-dimensional discrete Fourier Transform.
///
/// This function uses the `rustfft` library to perform the FFT.
///
/// # Arguments
/// * `input` - A mutable slice of complex numbers.
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
    let mut planner = FftPlanner::new();

    let fft = planner.plan_fft_forward(input.len());

    let mut buffer = input.to_vec();

    fft.process(&mut buffer);

    buffer
}

/// Computes the one-dimensional inverse discrete Fourier Transform.
///
/// This is the unnormalised forward transform followed by the `1/n`
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
    let n = input.len();

    if n == 0 {
        return vec![];
    }

    let mut planner = FftPlanner::new();

    let inverse = planner.plan_fft_inverse(n);

    let mut buffer = input.to_vec();

    inverse.process(&mut buffer);

    let scale = 1.0 / n as f64;

    for value in &mut buffer {
        *value *= scale;
    }

    buffer
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
    let n = a.len();

    let m = v.len();

    if n == 0 || m == 0 {
        return vec![];
    }

    let mut out = vec![0.0; n + m - 1];

    for i in 0..n {
        for j in 0..m {
            out[i + j] += a[i] * v[j];
        }
    }

    out
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
    let mut v_rev = v.to_vec();

    v_rev.reverse();

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
