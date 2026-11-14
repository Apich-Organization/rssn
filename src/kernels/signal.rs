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

/// Generates a Kaiser window of length `n` with shape parameter `beta`.
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
    let denom = crate::kernels::special::bessel_i0(beta);
    let half = 0.5 * (n.saturating_sub(1) as f64);
    (0..n)
        .map(|i| {
            let r = (i as f64 - half) / half;
            let arg = beta * (1.0 - r * r).max(0.0).sqrt();
            crate::kernels::special::bessel_i0(arg) / denom
        })
        .collect()
}

/// Orthonormal Type-II Discrete Cosine Transform (DCT-II).
#[must_use]
pub fn dct(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;
    let s0 = (1.0 / nf).sqrt();
    let sk = (2.0 / nf).sqrt();
    (0..n)
        .map(|k| {
            let scale = if k == 0 { s0 } else { sk };
            let kf = k as f64;
            let sum: f64 = x
                .iter()
                .enumerate()
                .map(|(i, &xi)| xi * (PI * kf * (i as f64 + 0.5) / nf).cos())
                .sum();
            scale * sum
        })
        .collect()
}

/// Inverse of [`dct`] (orthonormal Type-III DCT).
#[must_use]
pub fn idct(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;
    let s0 = (1.0 / nf).sqrt();
    let sk = (2.0 / nf).sqrt();
    (0..n)
        .map(|i| {
            let if64 = i as f64 + 0.5;
            x.iter()
                .enumerate()
                .map(|(k, &xk)| {
                    let scale = if k == 0 { s0 } else { sk };
                    scale * xk * (PI * (k as f64) * if64 / nf).cos()
                })
                .sum()
        })
        .collect()
}

/// Orthonormal Type-II Discrete Sine Transform (DST-II).
#[must_use]
pub fn dst(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;
    let sk = (2.0 / nf).sqrt();
    let sn = (1.0 / nf).sqrt();
    (0..n)
        .map(|k| {
            let scale = if k.saturating_add(1) == n {
                sn
            } else {
                sk
            };
            let kp1 = k.saturating_add(1) as f64;
            let sum: f64 = x
                .iter()
                .enumerate()
                .map(|(i, &xi)| xi * (PI * kp1 * (i as f64 + 0.5) / nf).sin())
                .sum();
            scale * sum
        })
        .collect()
}

/// Inverse of [`dst`] (orthonormal Type-III DST).
#[must_use]
pub fn idst(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let nf = n as f64;
    let sk = (2.0 / nf).sqrt();
    let sn = (1.0 / nf).sqrt();
    (0..n)
        .map(|i| {
            let if64 = i as f64 + 0.5;
            x.iter()
                .enumerate()
                .map(|(k, &xk)| {
                    let scale = if k.saturating_add(1) == n {
                        sn
                    } else {
                        sk
                    };
                    scale * xk * (PI * (k.saturating_add(1) as f64) * if64 / nf).sin()
                })
                .sum()
        })
        .collect()
}

/// Short-Time Fourier Transform (STFT) of `x` with window `window` and hop size `hop`.
#[must_use]
pub fn stft(
    x: &[f64],
    window: &[f64],
    hop: usize,
) -> Vec<Vec<Complex<f64>>> {
    let win_len = window.len();
    if win_len == 0 || hop == 0 || x.len() < win_len {
        return Vec::new();
    }
    let mut frames = Vec::new();
    let mut start = 0usize;
    while start.saturating_add(win_len) <= x.len() {
        let mut frame: Vec<Complex<f64>> = (0..win_len)
            .map(|j| {
                let xi = x.get(start.saturating_add(j)).copied().unwrap_or(0.0);
                let wi = window.get(j).copied().unwrap_or(0.0);
                Complex::new(xi * wi, 0.0)
            })
            .collect();
        crate::kernels::transforms::fft_slice(&mut frame);
        frames.push(frame);
        start = start.saturating_add(hop);
    }
    frames
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

/// Linear FIR/IIR digital filter (`a[0] * y[n] = sum_k b[k]*x[n-k] - sum_{k>=1} a[k]*y[n-k]`).
#[must_use]
pub fn lfilter(
    b: &[f64],
    a: &[f64],
    x: &[f64],
) -> Vec<f64> {
    if x.is_empty() || b.is_empty() || a.is_empty() {
        return Vec::new();
    }
    let a0 = a.first().copied().unwrap_or(1.0);
    if a0 == 0.0 {
        return vec![f64::NAN; x.len()];
    }
    let mut y = vec![0.0; x.len()];
    for n in 0..x.len() {
        let mut acc = 0.0;
        for (k, &bk) in b.iter().enumerate() {
            if n >= k {
                acc += (bk / a0) * x.get(n.saturating_sub(k)).copied().unwrap_or(0.0);
            }
        }
        for (k, &ak) in a.iter().enumerate().skip(1) {
            if n >= k {
                acc -= (ak / a0) * y.get(n.saturating_sub(k)).copied().unwrap_or(0.0);
            }
        }
        if let Some(slot) = y.get_mut(n) {
            *slot = acc;
        }
    }
    y
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
    (0..n)
        .map(|i| {
            let mut sum = 0.0;
            for (k, &ck) in coeffs.iter().enumerate() {
                let idx = (i as isize + k as isize - half as isize)
                    .clamp(0, n.saturating_sub(1) as isize) as usize;
                sum += ck * x.get(idx).copied().unwrap_or(0.0);
            }
            sum
        })
        .collect()
}
