//! JSON-based FFI API for numerical integration.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::integrate::QuadratureMethod;
use crate::numerical::integrate::{
    self,
};
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct QuadratureInput {
    expr: Expr,
    var: String,
    a: f64,
    b: f64,
    n_steps: usize,
    method: QuadratureMethod,
}

/// Performs numerical integration via JSON.
///
/// Input JSON format:
/// {
///   "expr": <Expr object>,
///   "var": "x",
///   "a": 0.0,
///   "b": 1.0,
///   "`n_steps"`: 100,
///   "method": "Simpson"
/// }
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

pub unsafe extern "C" fn rssn_numerical_quadrature_json(
    json_ptr: *const c_char
) -> *mut c_char {

    unsafe {

        if json_ptr.is_null() {

            return std::ptr::null_mut(
            );
        }

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

        let input: QuadratureInput =
            match serde_json::from_str(
                json_str,
            ) {
                | Ok(v) => v,
                | Err(e) => {

                    let res : FfiResult<f64, String> = FfiResult {
                ok : None,
                err : Some(format!(
                    "JSON deserialization error: {e}"
                )),
            };

                    return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
                },
            };

        let result =
            integrate::quadrature(
                &input.expr,
                &input.var,
                (input.a, input.b),
                input.n_steps,
                &input.method,
            );

        let res = match result {
            | Ok(val) => {
                FfiResult {
                    ok: Some(val),
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
