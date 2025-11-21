//! Handle-based FFI API for symbolic calculus functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::calculus;
use std::os::raw::c_char;
use std::ffi::CStr;

unsafe fn c_str_to_str<'a>(s: *const c_char) -> Option<&'a str> {
    if s.is_null() {
        None
    } else {
        CStr::from_ptr(s).to_str().ok()
    }
}

/// Differentiates an expression: d/d(var) expr.
/// 
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer and `var` is a valid C string.
#[no_mangle]
pub unsafe extern "C" fn rssn_differentiate(expr: *const Expr, var: *const c_char) -> *mut Expr {
    if expr.is_null() || var.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::differentiate(expr_ref, var_str)))
}

/// Integrates an expression: int(expr) d(var).
/// 
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer and `var` is a valid C string.
#[no_mangle]
pub unsafe extern "C" fn rssn_integrate(expr: *const Expr, var: *const c_char) -> *mut Expr {
    if expr.is_null() || var.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::integrate(expr_ref, var_str, None, None)))
}

/// Checks if an expression is analytic with respect to a variable.
/// 
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer and `var` is a valid C string.
#[no_mangle]
pub unsafe extern "C" fn rssn_check_analytic(expr: *const Expr, var: *const c_char) -> bool {
    if expr.is_null() || var.is_null() {
        return false;
    }
    let expr_ref = &*expr;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return false,
    };
    calculus::check_analytic(expr_ref, var_str)
}

/// Computes the limit of an expression: limit(expr, var -> point).
/// 
/// # Safety
/// The caller must ensure `expr` and `point` are valid Expr pointers and `var` is a valid C string.
#[no_mangle]
pub unsafe extern "C" fn rssn_limit(expr: *const Expr, var: *const c_char, point: *const Expr) -> *mut Expr {
    if expr.is_null() || var.is_null() || point.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let point_ref = &*point;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::limit(expr_ref, var_str, point_ref)))
}

/// Computes the definite integral of an expression.
#[no_mangle]
pub unsafe extern "C" fn rssn_definite_integrate(
    expr: *const Expr,
    var: *const c_char,
    lower: *const Expr,
    upper: *const Expr
) -> *mut Expr {
    if expr.is_null() || var.is_null() || lower.is_null() || upper.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let lower_ref = &*lower;
    let upper_ref = &*upper;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::definite_integrate(expr_ref, var_str, lower_ref, upper_ref)))
}

/// Evaluates an expression at a given point.
#[no_mangle]
pub unsafe extern "C" fn rssn_evaluate_at_point(
    expr: *const Expr,
    var: *const c_char,
    value: *const Expr
) -> *mut Expr {
    if expr.is_null() || var.is_null() || value.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let value_ref = &*value;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::evaluate_at_point(expr_ref, var_str, value_ref)))
}

/// Finds poles of an expression.
#[no_mangle]
pub unsafe extern "C" fn rssn_find_poles(expr: *const Expr, var: *const c_char) -> *mut Vec<Expr> {
    if expr.is_null() || var.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::find_poles(expr_ref, var_str)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_poles_len(poles: *const Vec<Expr>) -> usize {
    if poles.is_null() { 0 } else { (*poles).len() }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_poles_get(poles: *const Vec<Expr>, index: usize) -> *mut Expr {
    if poles.is_null() { return std::ptr::null_mut(); }
    let poles_ref = &*poles;
    if index >= poles_ref.len() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(poles_ref[index].clone()))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_free_poles(poles: *mut Vec<Expr>) {
    if !poles.is_null() {
        let _ = Box::from_raw(poles);
    }
}

/// Calculates the residue of a complex function at a given pole.
#[no_mangle]
pub unsafe extern "C" fn rssn_calculate_residue(expr: *const Expr, var: *const c_char, pole: *const Expr) -> *mut Expr {
    if expr.is_null() || var.is_null() || pole.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let pole_ref = &*pole;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::calculate_residue(expr_ref, var_str, pole_ref)))
}

/// Finds the order of a pole.
#[no_mangle]
pub unsafe extern "C" fn rssn_find_pole_order(expr: *const Expr, var: *const c_char, pole: *const Expr) -> usize {
    if expr.is_null() || var.is_null() || pole.is_null() {
        return 0;
    }
    let expr_ref = &*expr;
    let pole_ref = &*pole;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return 0,
    };
    calculus::find_pole_order(expr_ref, var_str, pole_ref)
}

/// Substitutes a variable with an expression.
#[no_mangle]
pub unsafe extern "C" fn rssn_substitute(expr: *const Expr, var: *const c_char, replacement: *const Expr) -> *mut Expr {
    if expr.is_null() || var.is_null() || replacement.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let replacement_ref = &*replacement;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::substitute(expr_ref, var_str, replacement_ref)))
}

/// Gets real and imaginary parts of an expression.
/// Returns a pointer to a tuple (Expr, Expr) - represented as Vec<Expr> of size 2 for simplicity?
/// Or return two out pointers?
/// I'll return a Vec<Expr> of size 2.
#[no_mangle]
pub unsafe extern "C" fn rssn_get_real_imag_parts(expr: *const Expr) -> *mut Vec<Expr> {
    if expr.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let (re, im) = calculus::get_real_imag_parts(expr_ref);
    Box::into_raw(Box::new(vec![re, im]))
}

/// Computes a path integral.
#[no_mangle]
pub unsafe extern "C" fn rssn_path_integrate(expr: *const Expr, var: *const c_char, contour: *const Expr) -> *mut Expr {
    if expr.is_null() || var.is_null() || contour.is_null() {
        return std::ptr::null_mut();
    }
    let expr_ref = &*expr;
    let contour_ref = &*contour;
    let var_str = match c_str_to_str(var) {
        Some(s) => s,
        None => return std::ptr::null_mut(),
    };
    Box::into_raw(Box::new(calculus::path_integrate(expr_ref, var_str, contour_ref)))
}
