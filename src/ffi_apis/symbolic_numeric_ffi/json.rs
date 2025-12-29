use std::ffi::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

/// Numerically evaluates a symbolic expression.

///

/// Takes a JSON string representing an `Expr` as input,

/// and returns a JSON string representing the numerical evaluation of that expression.

#[no_mangle]

pub extern "C" fn rssn_json_evaluate_numerical(
    expr_json: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    if let Some(e) = expr {

        if let Some(result) =
            evaluate_numerical(&e)
        {

            to_json_string(&result)
        } else {

            std::ptr::null_mut()
        }
    } else {

        std::ptr::null_mut()
    }
}
