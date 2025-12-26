//! JSON-based FFI for polynomial operations

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::*;
use std::os::raw::c_char;

/// Checks if an expression is a polynomial in the given variable (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_is_polynomial(
    expr_json: *const c_char,
    var: *const c_char,
) -> bool {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {

        is_polynomial(&e, v)
    } else {

        false
    }
}

/// Computes the degree of a polynomial (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_degree(
    expr_json: *const c_char,
    var: *const c_char,
) -> i64 {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {

        polynomial_degree(&e, v)
    } else {

        -1
    }
}

/// Performs polynomial long division (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_long_division(
    dividend_json: *const c_char,
    divisor_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let dividend: Option<Expr> = from_json_string(dividend_json);

    let divisor: Option<Expr> = from_json_string(divisor_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(d), Some(div), Some(v)) = (dividend, divisor, var_str) {

        let (quotient, remainder) = polynomial_long_division(&d, &div, v);

        let result = serde_json::json!({
            "quotient": quotient,
            "remainder": remainder
        });

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Finds the leading coefficient of a polynomial (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_leading_coefficient(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {

        let result = leading_coefficient(&e, v);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Converts polynomial to coefficient vector (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_to_coeffs_vec(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {

        let coeffs = to_polynomial_coeffs_vec(&e, v);

        to_json_string(&coeffs)
    } else {

        std::ptr::null_mut()
    }
}

/// Checks if an expression contains a variable (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_polynomial_contains_var(
    expr_json: *const c_char,
    var: *const c_char,
) -> bool {

    let expr: Option<Expr> = from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {

        contains_var(&e, v)
    } else {

        false
    }
}
