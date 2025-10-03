//! # Numerical Signal Processing
//!
//! This module provides numerical signal processing algorithms.
//! It includes implementations for the Fast Fourier Transform (FFT) and convolution,
//! which are fundamental operations in digital signal processing.

use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

/// Computes the one-dimensional discrete Fourier Transform.
///
/// This function uses the `rustfft` library to perform the FFT.
///
/// # Arguments
/// * `input` - A mutable slice of complex numbers.
///
/// # Returns
/// A vector of complex numbers representing the FFT of the input.
pub fn fft(input: &mut [Complex<f64>]) -> Vec<Complex<f64>> {
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(input.len());

    let mut buffer = input.to_vec();
    fft.process(&mut buffer);

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
pub fn convolve(a: &[f64], v: &[f64]) -> Vec<f64> {
    let n = a.len();
    let m = v.len();
    let mut out = vec![0.0; n + m - 1];

    for i in 0..n {
        for j in 0..m {
            out[i + j] += a[i] * v[j];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfft::num_complex::Complex;

    #[test]
    pub(crate) fn test_fft() {
        let mut input = vec![
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
        ];

        let output = fft(&mut input);

        let expected = vec![
            Complex::new(4.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];

        for (o, e) in output.iter().zip(expected.iter()) {
            assert!((o.re - e.re).abs() < 1e-9);
            assert!((o.im - e.im).abs() < 1e-9);
        }
    }

    #[test]
    pub(crate) fn test_convolve() {
        let a = vec![1.0, 2.0, 3.0];
        let v = vec![0.0, 1.0, 0.5];
        let result = convolve(&a, &v);
        let expected = vec![0.0, 1.0, 2.5, 4.0, 1.5];

        assert_eq!(result.len(), expected.len());
        for i in 0..result.len() {
            assert!((result[i] - expected[i]).abs() < 1e-9);
        }
    }
}
