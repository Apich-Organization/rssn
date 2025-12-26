//! JSON-based FFI API for numerical integration.

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::integrate::{
    self,
    QuadratureMethod,
};
use crate::symbolic::core::Expr;
use serde::Deserialize;
use std::ffi::{
    CStr,
    CString,
};
use std::os::raw::c_char;

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
///   "n_steps": 100,
///   "method": "Simpson"
/// }
#[no_mangle]

pub unsafe extern "C" fn rssn_numerical_quadrature_json(json_ptr: *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match CStr::from_ptr(json_ptr).to_str() {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let input: QuadratureInput = match serde_json::from_str(json_str) {
        | Ok(v) => v,
        | Err(e) => {

            let res: FfiResult<f64, String> = FfiResult {
                ok: None,
                err: Some(format!(
                    "JSON deserialization error: {}",
                    e
                )),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    let result = integrate::quadrature(
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

    CString::new(serde_json::to_string(&res).unwrap())
        .unwrap()
        .into_raw()
}
