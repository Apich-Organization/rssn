//! JSON-based FFI API for symbolic ODE functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::ode;
use std::os::raw::c_char;

/// Solves an ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_ode(
    ode_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let ode_expr: Option<Expr> = from_json_string(ode_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(ode), Some(f), Some(v)) = (
        ode_expr, func_str, var_str,
    ) {

        let result = ode::solve_ode(&ode, f, v, None);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a separable ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_separable_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v)) = (
        equation, func_str, var_str,
    ) {

        match ode::solve_separable_ode(&eq, f, v) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a first-order linear ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_first_order_linear_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v)) = (
        equation, func_str, var_str,
    ) {

        match ode::solve_first_order_linear_ode(&eq, f, v) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a Bernoulli ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_bernoulli_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v)) = (
        equation, func_str, var_str,
    ) {

        match ode::solve_bernoulli_ode(&eq, f, v) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a Riccati ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_riccati_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
    y1_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let y1: Option<Expr> = from_json_string(y1_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v), Some(y)) = (
        equation, func_str, var_str, y1,
    ) {

        match ode::solve_riccati_ode(&eq, f, v, &y) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves a Cauchy-Euler ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_cauchy_euler_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v)) = (
        equation, func_str, var_str,
    ) {

        match ode::solve_cauchy_euler_ode(&eq, f, v) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves an exact ODE using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_exact_ode(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v)) = (
        equation, func_str, var_str,
    ) {

        match ode::solve_exact_ode(&eq, f, v) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Solves by reduction of order using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_solve_by_reduction_of_order(
    equation_json: *const c_char,
    func: *const c_char,
    var: *const c_char,
    y1_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> = from_json_string(equation_json);

    let y1: Option<Expr> = from_json_string(y1_json);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(func)
                .to_str()
                .ok()
        }
    };

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(eq), Some(f), Some(v), Some(y)) = (
        equation, func_str, var_str, y1,
    ) {

        match ode::solve_by_reduction_of_order(&eq, f, v, &y) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}
