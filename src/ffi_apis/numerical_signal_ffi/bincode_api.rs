//! Bincode-based FFI API for numerical signal processing.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::signal;
use rustfft::num_complex::Complex;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct FftInput {
    data: Vec<Complex<f64>>,
}

#[derive(Deserialize)]

struct ConvolveInput {
    a: Vec<f64>,
    v: Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_fft_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let mut input: FftInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<Complex<f64>>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let result = signal::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_convolve_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: ConvolveInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let result = signal::convolve(&input.a, &input.v);

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_cross_correlation_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: ConvolveInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let result = signal::cross_correlation(&input.a, &input.v);

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
