//! # Numerical Integral Transforms
//!
//! This module provides numerical implementations of integral transforms,
//! specifically the Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT).
//! It includes an optimized in-place Cooley-Tukey algorithm and convenience functions
//! for `Vec<Complex<f64>>` inputs.

use num_complex::Complex;
use std::f64::consts::PI;

// --- CORE IN-PLACE FFT ALGORITHM ---

/// The optimized, in-place, iterative Cooley-Tukey Fast Fourier Transform (FFT) algorithm.
///
/// This function operates directly on the slice, avoiding memory allocations and copies.
/// It assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
pub(crate) fn fft_cooley_tukey_in_place(data: &mut [Complex<f64>]) {
    let n = data.len();
    if n <= 1 {
        return;
    }

    // 1. Bit-Reversal Permutation
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

    // 2. Butterfly Operations
    let mut len = 2; // Initial butterfly length
    while len <= n {
        let half_len = len / 2;
        // Twiddle Factor (w_m) calculation for the current stage
        let w_m = Complex::from_polar(1.0, -2.0 * PI / (len as f64));

        for i in (0..n).step_by(len) {
            let mut w = Complex::new(1.0, 0.0); // Current twiddle factor

            for k in 0..half_len {
                let even_idx = i + k;
                let odd_idx = i + k + half_len;

                let t = w * data[odd_idx];
                let u = data[even_idx];

                // Butterfly computation: In-place update
                data[even_idx] = u + t;
                data[odd_idx] = u - t;

                w = w * w_m;
            }
        }
        len *= 2;
    }
}

// --- PUBLIC VEC-BASED API (COMPATIBILITY) ---

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
///
/// Handles inputs of arbitrary length by padding with zeros to the next power of two.
/// This signature maintains compatibility with existing modules.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input sequence.
pub fn fft(data: &mut Vec<Complex<f64>>) {
    let n = data.len();
    if n <= 1 {
        return;
    }

    // Vec-specific: Pad with zeros to the next power of two using Vec::resize()
    let next_pow_of_2 = n.next_power_of_two();
    if n != next_pow_of_2 {
        data.resize(next_pow_of_2, Complex::new(0.0, 0.0));
    }

    // Use the optimized in-place function on the Vec's underlying data.
    fft_cooley_tukey_in_place(data);
}

/// Computes the Inverse Fast Fourier Transform (IFFT).
///
/// Handles inputs of arbitrary length by padding with zeros to the next power of two.
/// This signature maintains compatibility with existing modules.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input frequency-domain sequence.
pub fn ifft(data: &mut Vec<Complex<f64>>) {
    let n = data.len();
    if n <= 1 {
        return;
    }

    // Vec-specific: Pad with zeros to the next power of two using Vec::resize()
    let next_pow_of_2 = n.next_power_of_two();
    if n != next_pow_of_2 {
        data.resize(next_pow_of_2, Complex::new(0.0, 0.0));
    }

    // 1. First, take the conjugate of the input data
    for val in data.iter_mut() {
        *val = val.conj();
    }

    // 2. Apply the forward FFT (using the optimized in-place logic)
    fft_cooley_tukey_in_place(data);

    // 3. Then, take the conjugate of the result and scale by 1/N
    let n_f64 = data.len() as f64;
    for val in data.iter_mut() {
        *val = val.conj() / n_f64;
    }
}

// --- NEW SLICE-BASED API (FOR PARALLEL PROCESSING) ---

/// Computes the Forward FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// This function operates in-place and assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
pub fn fft_slice(data: &mut [Complex<f64>]) {
    fft_cooley_tukey_in_place(data);
}

/// Computes the Inverse FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// This function operates in-place and assumes the input slice length is a power of two.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input frequency-domain sequence.
pub fn ifft_slice(data: &mut [Complex<f64>]) {
    let n = data.len();
    if n <= 1 {
        return;
    }

    // 1. First, take the conjugate of the input data
    for val in data.iter_mut() {
        *val = val.conj();
    }

    // 2. Apply the forward FFT (using the optimized in-place logic)
    fft_cooley_tukey_in_place(data);

    // 3. Then, take the conjugate of the result and scale by 1/N
    let n_f64 = n as f64;
    for val in data.iter_mut() {
        *val = val.conj() / n_f64;
    }
}
