//! Handle-based FFI API for numerical transforms (FFT/IFFT).

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::transforms;
use num_complex::Complex;
use std::ptr;

/// Computes the Fast Fourier Transform (FFT) in-place.
///
/// # Arguments
/// * `real` - Pointer to the real parts of the input/output sequence.
/// * `imag` - Pointer to the imaginary parts of the input/output sequence.
/// * `len` - Length of the sequence. Must be a power of two for optimal performance.
///
/// # Returns
/// 0 on success, -1 on error.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_fft_inplace(
    real: *mut f64,
    imag: *mut f64,
    len: usize,
) -> i32 {
    if real.is_null() || imag.is_null() {
        update_last_error("Null pointer passed to rssn_num_fft_inplace".to_string());
        return -1;
    }
    if !len.is_power_of_two() {
        update_last_error("FFT length must be a power of two for in-place operation.".to_string());
        return -1;
    }

    let mut data: Vec<Complex<f64>> = (0..len)
        .map(|i| Complex::new(*real.add(i), *imag.add(i)))
        .collect();

    transforms::fft_slice(&mut data);

    for i in 0..len {
        *real.add(i) = data[i].re;
        *imag.add(i) = data[i].im;
    }

    0
}

/// Computes the Inverse Fast Fourier Transform (IFFT) in-place.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_ifft_inplace(
    real: *mut f64,
    imag: *mut f64,
    len: usize,
) -> i32 {
    if real.is_null() || imag.is_null() {
        update_last_error("Null pointer passed to rssn_num_ifft_inplace".to_string());
        return -1;
    }
    if !len.is_power_of_two() {
        update_last_error("IFFT length must be a power of two for in-place operation.".to_string());
        return -1;
    }

    let mut data: Vec<Complex<f64>> = (0..len)
        .map(|i| Complex::new(*real.add(i), *imag.add(i)))
        .collect();

    transforms::ifft_slice(&mut data);

    for i in 0..len {
        *real.add(i) = data[i].re;
        *imag.add(i) = data[i].im;
    }

    0
}
