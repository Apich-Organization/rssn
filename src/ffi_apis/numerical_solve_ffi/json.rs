//! JSON-based FFI API for numerical equation solvers.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::matrix::Matrix;
use crate::numerical::solve::{
    self,
};

#[derive(Deserialize)]

struct SolveLinearInput {
    matrix: Matrix<f64>,
    vector: Vec<f64>,
}

/// JSON FFI for solving linear systems.
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

pub unsafe extern "C" fn rssn_solve_linear_system_json(
    json_ptr: *const c_char
) -> *mut c_char {

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

    let input: SolveLinearInput =
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

    let result =
        solve::solve_linear_system(
            &input.matrix,
            &input.vector,
        );

    let res = match result {
        | Ok(sol) => {
            FfiResult {
                ok: Some(sol),
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
