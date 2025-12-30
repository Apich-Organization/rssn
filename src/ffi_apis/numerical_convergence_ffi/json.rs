//! JSON-based FFI API for numerical convergence operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::convergence;

#[derive(Deserialize)]

struct SeqInput {
    sequence: Vec<f64>,
}

/// JSON FFI for Aitken acceleration.
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

pub unsafe extern "C" fn rssn_convergence_aitken_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: SeqInput =
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

    let res = FfiResult {
        ok : Some(convergence::aitken_acceleration(&input.sequence)),
        err : None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for Richardson extrapolation.
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

pub unsafe extern "C" fn rssn_convergence_richardson_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: SeqInput =
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

    let res = FfiResult {
        ok : Some(convergence::richardson_extrapolation(&input.sequence)),
        err : None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for Wynn's epsilon algorithm.
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

pub unsafe extern "C" fn rssn_convergence_wynn_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: SeqInput =
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

    let res = FfiResult {
        ok: Some(
            convergence::wynn_epsilon(
                &input.sequence,
            ),
        ),
        err: None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}
