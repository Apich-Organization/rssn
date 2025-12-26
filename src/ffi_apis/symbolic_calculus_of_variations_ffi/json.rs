//! JSON-based FFI API for calculus of variations functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::calculus_of_variations;
use crate::symbolic::core::Expr;

/// Computes the Euler-Lagrange equation using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_euler_lagrange(
    lagrangian_json : *const c_char,
    func : *const c_char,
    var : *const c_char,
) -> *mut c_char {

    let lagrangian : Option<Expr> =
        from_json_string(
            lagrangian_json,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

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

    if let (Some(l), Some(f), Some(v)) = (
        lagrangian,
        func_str,
        var_str,
    ) {

        let result = calculus_of_variations::euler_lagrange(
            &l, f, v,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Generates and attempts to solve the Euler-Lagrange equation using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_euler_lagrange(
    lagrangian_json : *const c_char,
    func : *const c_char,
    var : *const c_char,
) -> *mut c_char {

    let lagrangian : Option<Expr> =
        from_json_string(
            lagrangian_json,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

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

    if let (Some(l), Some(f), Some(v)) = (
        lagrangian,
        func_str,
        var_str,
    ) {

        let result =
            calculus_of_variations::solve_euler_lagrange(
                &l, f, v,
            );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Applies Hamilton's Principle using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_hamiltons_principle(
    lagrangian_json : *const c_char,
    func : *const c_char,
    var : *const c_char,
) -> *mut c_char {

    let lagrangian : Option<Expr> =
        from_json_string(
            lagrangian_json,
        );

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

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

    if let (Some(l), Some(f), Some(v)) = (
        lagrangian,
        func_str,
        var_str,
    ) {

        let result =
            calculus_of_variations::hamiltons_principle(
                &l, f, v,
            );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
