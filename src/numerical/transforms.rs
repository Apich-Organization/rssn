//! # Numerical Integral Transforms
//!
//! This module provides numerical implementations of integral transforms,
//! specifically the Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT).
//! It includes an optimized in-place Cooley-Tukey algorithm and convenience functions
//! for `Vec<Complex<f64>>` inputs.

use std::f64::consts::PI;

use num_complex::Complex;

/// The optimized, in-place, iterative Cooley-Tukey Fast Fourier Transform (FFT) algorithm.
///
/// This function operates directly on the slice, avoiding memory allocations and copies.
/// It assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.

pub(crate) fn fft_cooley_tukey_in_place(
    data: &mut [Complex<f64>]
) {

    let n = data.len();

    if n <= 1 {

        return;
    }

    let mut j = 0;

    for i in 0 .. n {

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

        let w_m = Complex::from_polar(
            1.0,
            -2.0 * PI / (len as f64),
        );

        for i in (0 .. n).step_by(len) {

            let mut w =
                Complex::new(1.0, 0.0);

            for k in 0 .. half_len {

                let even_idx = i + k;

                let odd_idx =
                    i + k + half_len;

                let t =
                    w * data[odd_idx];

                let u = data[even_idx];

                data[even_idx] = u + t;

                data[odd_idx] = u - t;

                w *= w_m;
            }
        }

        len *= 2;
    }
}

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
///
/// Handles inputs of arbitrary length by padding with zeros to the next power of two.
/// This signature maintains compatibility with existing modules.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input sequence.
///
/// # Example
/// ```rust
/// 
/// use num_complex::Complex;
/// use rssn::numerical::transforms::fft;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
/// ];
///
/// fft(&mut data);
///
/// assert!((data[0].re - 4.0).abs() < 1e-9);
/// ```

pub fn fft(
    data: &mut Vec<Complex<f64>>
) {

    let n = data.len();

    if n <= 1 {

        return;
    }

    let next_pow_of_2 =
        n.next_power_of_two();

    if n != next_pow_of_2 {

        data.resize(
            next_pow_of_2,
            Complex::new(0.0, 0.0),
        );
    }

    fft_cooley_tukey_in_place(data);
}

/// Computes the Inverse Fast Fourier Transform (IFFT).
///
/// Handles inputs of arbitrary length by padding with zeros to the next power of two.
/// This signature maintains compatibility with existing modules.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// 
/// use num_complex::Complex;
/// use rssn::numerical::transforms::{fft, ifft};
///
/// let mut data = vec![
///     Complex::new(1.0, 2.0),
///     Complex::new(3.0, 4.0),
/// ];
///
/// let original = data.clone();
///
/// fft(&mut data);
///
/// ifft(&mut data);
///
/// // Note: ifft might return a padded version if the original length was not a power of 2
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
///
/// assert!((data[0].im - original[0].im).abs() < 1e-9);
/// ```

pub fn ifft(
    data: &mut Vec<Complex<f64>>
) {

    let n = data.len();

    if n <= 1 {

        return;
    }

    let next_pow_of_2 =
        n.next_power_of_two();

    if n != next_pow_of_2 {

        data.resize(
            next_pow_of_2,
            Complex::new(0.0, 0.0),
        );
    }

    for val in data.iter_mut() {

        *val = val.conj();
    }

    fft_cooley_tukey_in_place(data);

    let n_f64 = data.len() as f64;

    for val in data.iter_mut() {

        *val = val.conj() / n_f64;
    }
}

/// Computes the Forward FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// This function operates in-place and assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
///
/// # Example
/// ```rust
/// 
/// use num_complex::Complex;
/// use rssn::numerical::transforms::fft_slice;
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

pub fn fft_slice(
    data: &mut [Complex<f64>]
) {

    fft_cooley_tukey_in_place(data);
}

/// Computes the Inverse FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// This function operates in-place and assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// 
/// use num_complex::Complex;
/// use rssn::numerical::transforms::fft_slice;
/// use rssn::numerical::transforms::ifft_slice;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(0.0, 0.0),
/// ];
///
/// let original = data.clone();
///
/// fft_slice(&mut data);
///
/// ifft_slice(&mut data);
///
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
/// ```

pub fn ifft_slice(
    data: &mut [Complex<f64>]
) {

    let n = data.len();

    if n <= 1 {

        return;
    }

    for val in data.iter_mut() {

        *val = val.conj();
    }

    fft_cooley_tukey_in_place(data);

    let n_f64 = n as f64;

    for val in data.iter_mut() {

        *val = val.conj() / n_f64;
    }
}
