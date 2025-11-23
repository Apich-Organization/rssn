//! JSON-based FFI API for symbolic calculus functions.

use crate::ffi_apis::common::*;
use crate::symbolic::calculus;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

/// Differentiates an expression using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_differentiate(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {
        let result = calculus::differentiate(&e, v);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Integrates an expression using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_integrate(expr_json: *const c_char, var: *const c_char) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {
        let result = calculus::integrate(&e, v, None, None);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Computes definite integral using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_definite_integrate(
    expr_json: *const c_char,
    var: *const c_char,
    lower_json: *const c_char,
    upper_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let lower: Option<Expr> = from_json_string(lower_json);
    let upper: Option<Expr> = from_json_string(upper_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(l), Some(u)) = (expr, var_str, lower, upper) {
        let result = calculus::definite_integrate(&e, v, &l, &u);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Checks analytic using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_check_analytic(expr_json: *const c_char, var: *const c_char) -> bool {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {
        calculus::check_analytic(&e, v)
    } else {
        false
    }
}

/// Computes limit using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_limit(
    expr_json: *const c_char,
    var: *const c_char,
    point_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let point: Option<Expr> = from_json_string(point_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(p)) = (expr, var_str, point) {
        let result = calculus::limit(&e, v, &p);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Finds poles using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_find_poles(
    expr_json: *const c_char,
    var: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v)) = (expr, var_str) {
        let result = calculus::find_poles(&e, v);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Calculates residue using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_calculate_residue(
    expr_json: *const c_char,
    var: *const c_char,
    pole_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let pole: Option<Expr> = from_json_string(pole_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(p)) = (expr, var_str, pole) {
        let result = calculus::calculate_residue(&e, v, &p);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Finds pole order using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_find_pole_order(
    expr_json: *const c_char,
    var: *const c_char,
    pole_json: *const c_char,
) -> usize {
    let expr: Option<Expr> = from_json_string(expr_json);
    let pole: Option<Expr> = from_json_string(pole_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(p)) = (expr, var_str, pole) {
        calculus::find_pole_order(&e, v, &p)
    } else {
        0
    }
}

/// Substitutes using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_substitute(
    expr_json: *const c_char,
    var: *const c_char,
    replacement_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let replacement: Option<Expr> = from_json_string(replacement_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(r)) = (expr, var_str, replacement) {
        let result = calculus::substitute(&e, v, &r);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Gets real and imaginary parts using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_get_real_imag_parts(expr_json: *const c_char) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);

    if let Some(e) = expr {
        let (re, im) = calculus::get_real_imag_parts(&e);
        to_json_string(&vec![re, im])
    } else {
        std::ptr::null_mut()
    }
}

/// Computes path integral using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_path_integrate(
    expr_json: *const c_char,
    var: *const c_char,
    contour_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let contour: Option<Expr> = from_json_string(contour_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(c)) = (expr, var_str, contour) {
        let result = calculus::path_integrate(&e, v, &c);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Evaluates at point using JSON.
#[no_mangle]
pub extern "C" fn rssn_json_evaluate_at_point(
    expr_json: *const c_char,
    var: *const c_char,
    value_json: *const c_char,
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let value: Option<Expr> = from_json_string(value_json);
    let var_str = unsafe {
        if var.is_null() {
            None
        } else {
            std::ffi::CStr::from_ptr(var).to_str().ok()
        }
    };

    if let (Some(e), Some(v), Some(val)) = (expr, var_str, value) {
        let result = calculus::evaluate_at_point(&e, v, &val);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}
