//! # Numerical Integral Transforms
//!
//! The discrete Fourier transform of `n` complex samples in the
//! `X_k = sum_j x_j exp(-2 pi i j k / n)` convention, in place:
//!
//! * [`fft`] / [`ifft`] on a `Vec`, keeping its length;
//! * [`fft_slice`] / [`ifft_slice`] on a slice, for row/column passes.
//!
//! **Power-of-two lengths** take the iterative radix-2 Cooley-Tukey kernel
//! ([`fft_cooley_tukey_in_place`]) directly on the buffer: `log2(n)` stages of
//! `n/2` butterflies, no allocation.
//!
//! **Every other length** takes Bluestein's chirp-z algorithm. With
//! `w_j = exp(-i pi j^2 / n)` and `jk = (j^2 + k^2 - (k - j)^2) / 2`,
//!
//! ```text
//! X_k = w_k * sum_j (x_j w_j) conj(w_{k-j}),
//! ```
//!
//! which is a linear convolution of length `2n - 1`, evaluated by zero-padding
//! to `m = next_power_of_two(2n)` and using three radix-2 transforms of size
//! `m`. Only the first `n` outputs of that convolution are kept, so the result
//! is the `n`-point DFT itself — the padding is internal to the convolution and
//! never reaches the caller, who sees `n` coefficients for `n` samples.
//!
//! ## Limitations
//!
//! Length is no longer a restriction: every `n >= 1` is a true `n`-point
//! transform, `ifft(fft(x)) = x` at every length, and the returned vector has
//! the same length as the input. What remains true:
//!
//! * **Cost.** A power-of-two `n` costs `log2(n)` stages of `n/2` butterflies.
//!   Any other `n` costs three transforms of the next power of two `m >= 2n`,
//!   i.e. `3 m log2(m)` butterfly stages against `n log2(n)` for the
//!   power-of-two path of the same length — about `6.6x` at `n = 1000`
//!   (`m = 2048`) — plus two `O(m)` scratch buffers.
//! * **No prime special case.** Prime `n` goes through Bluestein like any other
//!   non-power-of-two length; there is no Rader transform.
//! * **Complex input only.** Real sequences pay for a complex transform; there
//!   is no real-to-complex packing and no half-length trick.
//! * **Single-threaded.** Each call is serial; callers parallelise over
//!   independent rows/columns with the slice variants.
//! * **Normalisation.** [`fft`] is unnormalised, [`ifft`] divides by `n`; the
//!   `1/n` is applied to the complex conjugate, matching the radix-2 path
//!   exactly.
//! * **Accuracy.** Accumulation is `f64`. The round trip and Parseval hold to
//!   about `1e-14 * n` relative, and Bluestein adds an `O(eps * n)` chirp-phase
//!   error on top of the radix-2 round-off.

use std::f64::consts::PI;

use num_complex::Complex;

/// The optimized, in-place, iterative Cooley-Tukey Fast Fourier Transform (FFT) algorithm.
///
/// This function operates directly on the slice, avoiding memory allocations and copies.
/// It assumes the input slice length is a power of two; other lengths are
/// handled by the public wrappers through Bluestein's algorithm.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
pub(crate) fn fft_cooley_tukey_in_place(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    let mut j = 0;

    for i in 0..n {
        if i < j {
            data.swap(i, j);
        }

        let mut bit = n >> 1;

        while j & bit != 0 {
            j ^= bit;

            bit >>= 1;
        }

        j ^= bit;
    }

    let mut len = 2;

    while len <= n {
        let half_len = len / 2;

        let w_m = Complex::from_polar(1.0, -2.0 * PI / (len as f64));

        for i in (0..n).step_by(len) {
            let mut w = Complex::new(1.0, 0.0);

            for k in 0..half_len {
                let even_idx = i + k;

                let odd_idx = i + k + half_len;

                let t = w * data[odd_idx];

                let u = data[even_idx];

                data[even_idx] = u + t;

                data[odd_idx] = u - t;

                w *= w_m;
            }
        }

        len *= 2;
    }
}

/// The inverse radix-2 transform: `conj(fft(conj(x))) / n`, in place.
///
/// The input length must be a power of two.
fn inverse_pow2_in_place(data: &mut [Complex<f64>]) {
    for val in data.iter_mut() {
        *val = val.conj();
    }

    fft_cooley_tukey_in_place(data);

    let n_f64 = data.len() as f64;

    for val in data.iter_mut() {
        *val = val.conj() / n_f64;
    }
}

/// Bluestein's chirp-z transform of an arbitrary length, returning `n` samples.
///
/// `X_k = w_k * sum_j (x_j w_j) conj(w_{k-j})` with `w_j = exp(-i pi j^2 / n)`;
/// the sum is a linear convolution of length `2n - 1` done cyclically at
/// `m = next_power_of_two(2n)`, where no wrap-around can reach the first `n`
/// outputs. The chirp angle uses `j^2 mod 2n`, so the phase stays exact for
/// `j` up to the square root of the machine word.
fn bluestein_forward(data: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = data.len();

    let two_n = n.saturating_mul(2);

    let chirp: Vec<Complex<f64>> = (0..n)
        .map(|j| {
            let jj = j.wrapping_mul(j) % two_n;

            Complex::from_polar(1.0, -PI * (jj as f64) / (n as f64))
        })
        .collect();

    let m = two_n.next_power_of_two();

    let mut fa: Vec<Complex<f64>> = data.iter().zip(chirp.iter()).map(|(x, w)| x * w).collect();

    fa.resize(m, Complex::new(0.0, 0.0));

    // b_d = conj(w_d) and b is even in d, so the mirror image fills m - d.
    let mut fb: Vec<Complex<f64>> = (0..m)
        .map(|p| {
            let d = p.min(m.saturating_sub(p));

            chirp.get(d).map_or(Complex::new(0.0, 0.0), Complex::conj)
        })
        .collect();

    fft_cooley_tukey_in_place(&mut fa);

    fft_cooley_tukey_in_place(&mut fb);

    for (a, b) in fa.iter_mut().zip(fb.iter()) {
        *a *= *b;
    }

    inverse_pow2_in_place(&mut fa);

    fa.iter().zip(chirp.iter()).map(|(c, w)| c * w).collect()
}

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
///
/// Works at any length and leaves the length unchanged: a power-of-two `n` uses
/// the in-place radix-2 kernel, any other `n` uses Bluestein's chirp-z
/// algorithm. Either way the result is the exact `n`-point DFT, `O(n log n)`.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(0.0, 0.0),
///     Complex::new(0.0, 0.0),
/// ];
///
/// fft(&mut data);
///
/// assert!((data[0].re - 2.0).abs() < 1e-9);
/// ```
///
/// A length that is not a power of two keeps its length:
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
///
/// let mut data = vec![Complex::new(1.0, 0.0); 3];
///
/// fft(&mut data);
///
/// assert_eq!(data.len(), 3);
/// // [3, 0, 0] for a constant signal.
/// assert!((data[0].re - 3.0).abs() < 1e-12);
/// ```
pub fn fft(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    if n.is_power_of_two() {
        fft_cooley_tukey_in_place(data);
    } else {
        data.copy_from_slice(&bluestein_forward(data));
    }
}

/// Computes the Inverse Fast Fourier Transform (IFFT).
///
/// Works at any length and leaves the length unchanged, so `ifft(fft(x))`
/// returns `x` at every length. The `1/n` normalisation is included.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
/// use rssn::kernels::transforms::ifft;
///
/// let mut data = vec![Complex::new(1.0, 2.0), Complex::new(3.0, 4.0)];
///
/// let original = data.clone();
///
/// fft(&mut data);
///
/// ifft(&mut data);
///
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
///
/// assert!((data[0].im - original[0].im).abs() < 1e-9);
/// ```
pub fn ifft(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    if n.is_power_of_two() {
        inverse_pow2_in_place(data);

        return;
    }

    let conjugated: Vec<Complex<f64>> = data.iter().map(Complex::conj).collect();

    let n_f64 = n as f64;

    let transformed = bluestein_forward(&conjugated);

    for (val, out) in data.iter_mut().zip(transformed.iter()) {
        *val = out.conj() / n_f64;
    }
}

/// Computes the Forward FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// Works at any length and operates in place; a non-power-of-two length is
/// routed through Bluestein's algorithm and copied back into the slice.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft_slice;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(0.0, 0.0),
///     Complex::new(0.0, 0.0),
/// ];
///
/// fft_slice(&mut data);
///
/// assert!((data[0].re - 2.0).abs() < 1e-9);
/// ```
pub fn fft_slice(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    if n.is_power_of_two() {
        fft_cooley_tukey_in_place(data);
    } else {
        let transformed = bluestein_forward(data);

        for (val, out) in data.iter_mut().zip(transformed.iter()) {
            *val = *out;
        }
    }
}

/// Computes the Inverse FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// Works at any length and operates in place; the `1/n` normalisation is
/// included.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft_slice;
/// use rssn::kernels::transforms::ifft_slice;
///
/// let mut data = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
///
/// let original = data.clone();
///
/// fft_slice(&mut data);
///
/// ifft_slice(&mut data);
///
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
/// ```
pub fn ifft_slice(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    if n.is_power_of_two() {
        inverse_pow2_in_place(data);

        return;
    }

    let conjugated: Vec<Complex<f64>> = data.iter().map(Complex::conj).collect();

    let n_f64 = n as f64;

    let transformed = bluestein_forward(&conjugated);

    for (val, out) in data.iter_mut().zip(transformed.iter()) {
        *val = out.conj() / n_f64;
    }
}
