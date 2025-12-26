use std::f64::consts::PI;

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
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::signal::fft;
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
/// assert!(
///     (output[0].re - 4.0).abs()
///         < 1e-9
/// );
/// ```

pub fn fft(
    input: &mut [Complex<f64>]
) -> Vec<Complex<f64>> {

    let mut planner = FftPlanner::new();

    let fft = planner
        .plan_fft_forward(input.len());

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
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::signal::convolve;
///
/// let a = vec![1.0, 2.0, 3.0];
///
/// let v = vec![0.0, 1.0, 0.5];
///
/// let res = convolve(&a, &v);
///
/// assert_eq!(
///     res,
///     vec![0.0, 1.0, 2.5, 4.0, 1.5]
/// );
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

    for i in 0 .. n {

        for j in 0 .. m {

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
/// 
/// use rssn::numerical::signal::cross_correlation;
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

pub fn hann_window(
    n: usize
) -> Vec<f64> {

    if n == 0 {

        return vec![];
    }

    if n == 1 {

        return vec![1.0];
    }

    (0 .. n)
        .map(|i| {

            0.5 * (1.0
                - (2.0 * PI * i as f64
                    / (n - 1) as f64)
                    .cos())
        })
        .collect()
}

/// Generates a Hamming window of length `n`.
///
/// # Arguments
/// * `n` - The number of points in the output window.
#[must_use]

pub fn hamming_window(
    n: usize
) -> Vec<f64> {

    if n == 0 {

        return vec![];
    }

    if n == 1 {

        return vec![1.0];
    }

    (0 .. n)
        .map(|i| {

            0.46f64.mul_add(
                -(2.0 * PI * i as f64
                    / (n - 1) as f64)
                    .cos(),
                0.54,
            )
        })
        .collect()
}
