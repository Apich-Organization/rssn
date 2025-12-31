//! JSON-based FFI for polynomial operations

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::contains_var;
use crate::symbolic::polynomial::is_polynomial;
use crate::symbolic::polynomial::leading_coefficient;
use crate::symbolic::polynomial::polynomial_degree;
use crate::symbolic::polynomial::polynomial_long_division;
use crate::symbolic::polynomial::to_polynomial_coeffs_vec;

/// Checks if an expression is a polynomial in the given variable (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_is_polynomial(
    expr_json: *const c_char,
    var: *const c_char,
) -> bool {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str) {
        | (Some(e), Some(v)) => {
            is_polynomial(&e, v)
        },
        | _ => false,
    }
}

/// Computes the degree of a polynomial (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_degree(
    expr_json: *const c_char,
    var: *const c_char,
) -> i64 {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str) {
        | (Some(e), Some(v)) => {
            polynomial_degree(&e, v)
        },
        | _ => -1,
    }
}

/// Performs polynomial long division (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_long_division(
    dividend_json: *const c_char,
    divisor_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let dividend: Option<Expr> =
        from_json_string(dividend_json);

    let divisor: Option<Expr> =
        from_json_string(divisor_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        dividend,
        divisor,
        var_str,
    ) {
        | (
            Some(d),
            Some(div),
            Some(v),
        ) => {

            let (quotient, remainder) =
                polynomial_long_division(
                    &d, &div, v,
                );

            let result = serde_json::json!({
                "quotient": quotient,
                "remainder": remainder
            });

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Finds the leading coefficient of a polynomial (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_leading_coefficient(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str) {
        | (Some(e), Some(v)) => {

            let result =
                leading_coefficient(
                    &e, v,
                );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Converts polynomial to coefficient vector (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_to_coeffs_vec(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str) {
        | (Some(e), Some(v)) => {

            let coeffs =
            to_polynomial_coeffs_vec(
                &e, v,
            );

            to_json_string(&coeffs)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Checks if an expression contains a variable (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_polynomial_contains_var(
    expr_json: *const c_char,
    var: *const c_char,
) -> bool {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str) {
        | (Some(e), Some(v)) => {
            contains_var(&e, v)
        },
        | _ => false,
    }
}
