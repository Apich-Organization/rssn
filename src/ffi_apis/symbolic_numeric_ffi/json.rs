use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

#[no_mangle]

pub extern "C" fn rssn_json_evaluate_numerical(expr_json : *const c_char) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {

        if let Some(result) = evaluate_numerical(&e) {

            to_json_string(&result)
        } else {

            std::ptr::null_mut()
        }
    } else {

        std::ptr::null_mut()
    }
}
