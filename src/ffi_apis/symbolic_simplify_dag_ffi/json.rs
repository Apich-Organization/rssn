//! JSON-based FFI API for symbolic simplify_dag functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag;

/// Simplifies an expression using the DAG-based simplifier (JSON input/output).
#[no_mangle]

pub extern "C" fn rssn_json_simplify_dag(
    expr_json: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    if let Some(e) = expr {

        let result =
            simplify_dag::simplify(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
