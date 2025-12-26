use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::integration::*;

/// Integrates an expression using the Risch-Norman algorithm (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_risch_norman_integrate(
    expr_json : *const c_char,
    x_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let x : Option<String> = from_json_string(x_json);

    if let (Some(e), Some(var)) = (expr, x) {

        let result = risch_norman_integrate(&e, &var);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Integrates a rational function (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_integrate_rational_function(
    expr_json : *const c_char,
    x_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let x : Option<String> = from_json_string(x_json);

    if let (Some(e), Some(var)) = (expr, x) {

        match integrate_rational_function_expr(&e, &var) {
            | Ok(result) => to_json_string(&result),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}
