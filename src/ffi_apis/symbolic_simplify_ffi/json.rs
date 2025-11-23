//! JSON-based FFI API for symbolic simplify functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify;
use std::os::raw::c_char;

/// Simplifies an expression using the heuristic simplifier (JSON input/output).
#[no_mangle]
pub extern "C" fn rssn_json_heuristic_simplify(expr_json: *const c_char) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {
        #[allow(deprecated)]
        let result = simplify::heuristic_simplify(e);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Simplifies an expression using the legacy simplifier (JSON input/output).
#[no_mangle]
pub extern "C" fn rssn_json_simplify(expr_json: *const c_char) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {
        #[allow(deprecated)]
        let result = simplify::simplify(e);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}
