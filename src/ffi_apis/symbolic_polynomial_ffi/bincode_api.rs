//! Bincode-based FFI for polynomial operations

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::*;

/// Checks if an expression is a polynomial in the given variable (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_is_polynomial(
    expr_buf : BincodeBuffer,
    var : *const c_char,
) -> bool {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

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

/// Computes the degree of a polynomial (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_degree(
    expr_buf : BincodeBuffer,
    var : *const c_char,
) -> i64 {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

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

/// Performs polynomial long division (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_long_division(
    dividend_buf : BincodeBuffer,
    divisor_buf : BincodeBuffer,
    var : *const c_char,
) -> BincodeBuffer {

    let dividend : Option<Expr> = from_bincode_buffer(&dividend_buf);

    let divisor : Option<Expr> = from_bincode_buffer(&divisor_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(d), Some(div), Some(v)) = (
        dividend,
        divisor,
        var_str,
    ) {

        let (quotient, remainder) = polynomial_long_division(&d, &div, v);

        let result = (quotient, remainder);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Finds the leading coefficient of a polynomial (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_leading_coefficient(
    expr_buf : BincodeBuffer,
    var : *const c_char,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

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

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Converts polynomial to coefficient vector (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_to_coeffs_vec(
    expr_buf : BincodeBuffer,
    var : *const c_char,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

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

        to_bincode_buffer(&coeffs)
    } else {

        BincodeBuffer::empty()
    }
}

/// Checks if an expression contains a variable (bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_polynomial_contains_var(
    expr_buf : BincodeBuffer,
    var : *const c_char,
) -> bool {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

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
