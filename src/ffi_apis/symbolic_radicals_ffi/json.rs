use std::ffi::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::radicals::denest_sqrt;
use crate::symbolic::radicals::simplify_radicals;

/// Simplifies radical expressions (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_simplify_radicals(
    expr_json: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    if let Some(e) = expr {

        let result =
            simplify_radicals(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Denests a nested square root (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_denest_sqrt(
    expr_json: *const c_char
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    if let Some(e) = expr {

        let result = denest_sqrt(&e);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
