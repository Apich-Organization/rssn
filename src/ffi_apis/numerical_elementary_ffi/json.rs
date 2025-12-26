//! JSON-based FFI API for numerical elementary operations.

use std::collections::HashMap;
use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::elementary;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct EvalRequest {
    expr : Expr,
    vars : HashMap<String, f64>,
}

/// Evaluates an expression from a JSON string.
///
/// Input JSON format: `{"expr": <Expr>, "vars": {"x": 1.0, "y": 2.0}}`
#[no_mangle]

pub unsafe extern "C" fn rssn_num_eval_json(json_ptr : *const c_char) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return std::ptr::null_mut(),
    };

    let req : EvalRequest = match serde_json::from_str(json_str) {
        | Ok(r) => r,
        | Err(e) => {

            let res : FfiResult<f64, String> = FfiResult {
                ok : None,
                err : Some(e.to_string()),
            };

            return CString::new(serde_json::to_string(&res).unwrap())
                .unwrap()
                .into_raw();
        },
    };

    let result = elementary::eval_expr(&req.expr, &req.vars);

    let ffi_res = match result {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    CString::new(serde_json::to_string(&ffi_res).unwrap())
        .unwrap()
        .into_raw()
}
