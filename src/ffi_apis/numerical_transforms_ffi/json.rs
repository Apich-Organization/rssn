//! JSON-based FFI API for numerical transforms.

use std::os::raw::c_char;

use num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::transforms;

#[derive(Deserialize)]

struct TransformInput {
    data: Vec<Complex<f64>>,
}

/// Computes the Fast Fourier Transform (FFT) in-place via JSON serialization.
///
/// The FFT converts a sequence from the time/space domain to the frequency domain,
/// computing X(k) = Σx(n)e^(-2πikn/N) for k = 0, ..., N-1.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `data`: Array of complex numbers to transform
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Complex<f64>>, String>` with
/// FFT-transformed data in frequency domain.
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

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_fft_json(
    input_json: *const c_char
) -> *mut c_char {

    let mut input : TransformInput = match from_json_string(input_json) {
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

    transforms::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

/// Computes the Inverse Fast Fourier Transform (IFFT) in-place via JSON serialization.
///
/// The IFFT converts a sequence from the frequency domain back to the time/space domain,
/// computing x(n) = (1/N) ΣX(k)e^(2πikn/N) for n = 0, ..., N-1.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `data`: Array of complex numbers in frequency domain
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Complex<f64>>, String>` with
/// IFFT-transformed data in time/space domain.
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

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_ifft_json(
    input_json: *const c_char
) -> *mut c_char {

    let mut input : TransformInput = match from_json_string(input_json) {
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

    transforms::ifft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
