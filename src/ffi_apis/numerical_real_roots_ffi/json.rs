//! JSON-based FFI API for numerical real root finding.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;

#[derive(Deserialize)]

struct FindRootsInput {
    coeffs: Vec<f64>,
    tolerance: f64,
}

/// Finds all real roots of a polynomial using numerical methods and JSON serialization.
///
/// Uses root-finding algorithms to locate all real zeros of the polynomial.
///
/// # Arguments
///
/// * `json_ptr` - A JSON string pointer containing:
///   - `coeffs`: Polynomial coefficients [a₀, a₁, ..., aₙ] for a₀ + a₁x + ... + aₙxⁿ
///   - `tolerance`: Convergence tolerance for root finding
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// an array of real roots found.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

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

pub unsafe extern "C" fn rssn_real_roots_find_roots_json(
    json_ptr: *const c_char
) -> *mut c_char {

    unsafe {

        let json_str = match CStr::from_ptr(
        json_ptr,
    )
    .to_str()
    {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut()
        },
    };

        let input: FindRootsInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

        let poly = Polynomial::new(
            input.coeffs,
        );

        let result =
            real_roots::find_roots(
                &poly,
                input.tolerance,
            );

        let res = match result {
            | Ok(roots) => {
                FfiResult {
                    ok: Some(roots),
                    err: None::<String>,
                }
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

        CString::new(
            serde_json::to_string(&res)
                .unwrap(),
        )
        .unwrap()
        .into_raw()
    }
}
