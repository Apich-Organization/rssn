//! Handle-based FFI API for numerical signal processing.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::matrix::Matrix;
use crate::numerical::signal;
use rustfft::num_complex::Complex;
use std::ptr;

/// Computes the FFT and returns a Matrix<Complex<f64>> as a Matrix<f64> (real, imag interleaved).
#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_fft(
    real: *const f64,
    imag: *const f64,
    len: usize,
) -> *mut Matrix<f64> {

    if real.is_null() || imag.is_null() {

        update_last_error("Null pointer passed to rssn_num_signal_fft".to_string());

        return ptr::null_mut();
    }

    let mut input: Vec<Complex<f64>> = (0..len)
        .map(|i| {

            Complex::new(
                *real.add(i),
                *imag.add(i),
            )
        })
        .collect();

    let output = signal::fft(&mut input);

    let mut flat = Vec::with_capacity(output.len() * 2);

    for c in output {

        flat.push(c.re);

        flat.push(c.im);
    }

    Box::into_raw(Box::new(
        Matrix::new(len, 2, flat),
    ))
}

/// Computes the convolution of two sequences.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_convolve(
    a: *const f64,
    a_len: usize,
    v: *const f64,
    v_len: usize,
) -> *mut Matrix<f64> {

    if a.is_null() || v.is_null() {

        update_last_error("Null pointer passed to rssn_num_signal_convolve".to_string());

        return ptr::null_mut();
    }

    let a_slice = std::slice::from_raw_parts(a, a_len);

    let v_slice = std::slice::from_raw_parts(v, v_len);

    let result = signal::convolve(a_slice, v_slice);

    let n = result.len();

    Box::into_raw(Box::new(
        Matrix::new(1, n, result),
    ))
}

/// Computes the cross-correlation of two sequences.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_cross_correlation(
    a: *const f64,
    a_len: usize,
    v: *const f64,
    v_len: usize,
) -> *mut Matrix<f64> {

    if a.is_null() || v.is_null() {

        update_last_error("Null pointer passed to rssn_num_signal_cross_correlation".to_string());

        return ptr::null_mut();
    }

    let a_slice = std::slice::from_raw_parts(a, a_len);

    let v_slice = std::slice::from_raw_parts(v, v_len);

    let result = signal::cross_correlation(a_slice, v_slice);

    let n = result.len();

    Box::into_raw(Box::new(
        Matrix::new(1, n, result),
    ))
}

/// Generates a Hann window.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_hann_window(n: usize) -> *mut Matrix<f64> {

    let window = signal::hann_window(n);

    Box::into_raw(Box::new(
        Matrix::new(1, n, window),
    ))
}

/// Generates a Hamming window.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_hamming_window(n: usize) -> *mut Matrix<f64> {

    let window = signal::hamming_window(n);

    Box::into_raw(Box::new(
        Matrix::new(1, n, window),
    ))
}
