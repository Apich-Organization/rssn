//! Handle-based FFI for polynomial operations

use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::*;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Checks if an expression is a polynomial in the given variable (handle-based)
#[no_mangle]

pub extern "C" fn polynomial_is_polynomial_handle(
    expr_handle: *const Expr,
    var: *const c_char,
) -> bool {

    let expr = unsafe {

        &*expr_handle
    };

    let var_str = unsafe {

        CStr::from_ptr(var)
            .to_str()
            .unwrap()
    };

    is_polynomial(expr, var_str)
}

/// Computes the degree of a polynomial (handle-based)
#[no_mangle]

pub extern "C" fn polynomial_degree_handle(expr_handle: *const Expr, var: *const c_char) -> i64 {

    let expr = unsafe {

        &*expr_handle
    };

    let var_str = unsafe {

        CStr::from_ptr(var)
            .to_str()
            .unwrap()
    };

    polynomial_degree(expr, var_str)
}

/// Performs polynomial long division (handle-based)
#[no_mangle]

pub extern "C" fn polynomial_long_division_handle(
    dividend_handle: *const Expr,
    divisor_handle: *const Expr,
    var: *const c_char,
    quotient_out: *mut *mut Expr,
    remainder_out: *mut *mut Expr,
) {

    let dividend = unsafe {

        &*dividend_handle
    };

    let divisor = unsafe {

        &*divisor_handle
    };

    let var_str = unsafe {

        CStr::from_ptr(var)
            .to_str()
            .unwrap()
    };

    let (quotient, remainder) = polynomial_long_division(dividend, divisor, var_str);

    unsafe {

        *quotient_out = Box::into_raw(Box::new(quotient));

        *remainder_out = Box::into_raw(Box::new(remainder));
    }
}

/// Finds the leading coefficient of a polynomial (handle-based)
#[no_mangle]

pub extern "C" fn polynomial_leading_coefficient_handle(
    expr_handle: *const Expr,
    var: *const c_char,
) -> *mut Expr {

    let expr = unsafe {

        &*expr_handle
    };

    let var_str = unsafe {

        CStr::from_ptr(var)
            .to_str()
            .unwrap()
    };

    let result = leading_coefficient(expr, var_str);

    Box::into_raw(Box::new(result))
}

/// Checks if an expression contains a variable (handle-based)
#[no_mangle]

pub extern "C" fn polynomial_contains_var_handle(
    expr_handle: *const Expr,
    var: *const c_char,
) -> bool {

    let expr = unsafe {

        &*expr_handle
    };

    let var_str = unsafe {

        CStr::from_ptr(var)
            .to_str()
            .unwrap()
    };

    contains_var(expr, var_str)
}

/// Frees an Expr handle
#[no_mangle]

pub extern "C" fn polynomial_free_expr_handle(expr_handle: *mut Expr) {

    if !expr_handle.is_null() {

        unsafe {

            let _ = Box::from_raw(expr_handle);
        }
    }
}
