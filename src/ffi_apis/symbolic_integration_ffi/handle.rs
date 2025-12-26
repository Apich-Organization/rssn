use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::integration::*;

/// Integrates an expression using the Risch-Norman algorithm (Handle)
#[no_mangle]

pub extern "C" fn rssn_risch_norman_integrate_handle(
    expr : *const Expr,
    x : *const c_char,
) -> *mut Expr {

    let expr_ref = unsafe {

        if expr.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*expr
    };

    let x_str = unsafe {

        if x.is_null() {

            return std::ptr::null_mut(
            );
        }

        CStr::from_ptr(x)
            .to_string_lossy()
            .into_owned()
    };

    let result = risch_norman_integrate(
        expr_ref,
        &x_str,
    );

    Box::into_raw(Box::new(result))
}

/// Integrates a rational function (Handle)
#[no_mangle]

pub extern "C" fn rssn_integrate_rational_function_handle(
    expr : *const Expr,
    x : *const c_char,
) -> *mut Expr {

    let expr_ref = unsafe {

        if expr.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*expr
    };

    let x_str = unsafe {

        if x.is_null() {

            return std::ptr::null_mut(
            );
        }

        CStr::from_ptr(x)
            .to_string_lossy()
            .into_owned()
    };

    match integrate_rational_function_expr(expr_ref, &x_str) {
        | Ok(result) => Box::into_raw(Box::new(result)),
        | Err(_) => std::ptr::null_mut(),
    }
}
