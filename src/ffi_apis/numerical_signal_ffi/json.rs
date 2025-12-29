//! JSON-based FFI API for numerical signal processing.

use std::os::raw::c_char;

use rustfft::num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::signal;

#[derive(Deserialize)]

struct FftInput {
    data: Vec<Complex<f64>>,
}

#[derive(Deserialize)]

struct ConvolveInput {
    a: Vec<f64>,
    v: Vec<f64>,
}

/// Computes the Fast Fourier Transform (FFT) of complex data using JSON serialization.
///
/// The FFT converts a signal from time/space domain to frequency domain using the
/// Cooley-Tukey algorithm in O(N log N) time.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `data`: Array of complex numbers to transform
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Complex<f64>>, String>` with
/// the frequency domain representation.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_signal_fft_json(
    input_json: *const c_char
) -> *mut c_char {

    let mut input : FftInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Complex<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        signal::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

/// Computes the discrete convolution of two signals using JSON serialization.
///
/// The convolution is defined as (a * v)[n] = Σ a[k]v[n-k], representing the
/// combined effect of two systems or filtering operation.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `a`: First signal array
///   - `v`: Second signal array (kernel)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the convolution result.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_signal_convolve_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : ConvolveInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = signal::convolve(
        &input.a,
        &input.v,
    );

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

/// Computes the cross-correlation of two signals using JSON serialization.
///
/// Cross-correlation measures similarity between signals as a function of lag:
/// (a ⋆ v)[n] = Σ a[k]v[n+k].
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `a`: First signal array
///   - `v`: Second signal array
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the cross-correlation result.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_signal_cross_correlation_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : ConvolveInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        signal::cross_correlation(
            &input.a,
            &input.v,
        );

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
