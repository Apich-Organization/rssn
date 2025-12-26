//! Handle-based FFI API for calculus of variations functions.

use crate::symbolic::calculus_of_variations;
use crate::symbolic::core::Expr;
use std::ffi::CStr;
use std::os::raw::c_char;

unsafe fn c_str_to_str<'a>(s: *const c_char) -> Option<&'a str> {

    if s.is_null() {

        None
    } else {

        CStr::from_ptr(s)
            .to_str()
            .ok()
    }
}

/// Computes the Euler-Lagrange equation for a given Lagrangian.
///
/// # Safety
/// The caller must ensure `lagrangian` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[no_mangle]

pub unsafe extern "C" fn rssn_euler_lagrange(
    lagrangian: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    if lagrangian.is_null() || func.is_null() || var.is_null() {

        return std::ptr::null_mut();
    }

    let lagrangian_ref = &*lagrangian;

    let func_str = match c_str_to_str(func) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let var_str = match c_str_to_str(var) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    Box::into_raw(Box::new(
        calculus_of_variations::euler_lagrange(
            lagrangian_ref,
            func_str,
            var_str,
        ),
    ))
}

/// Generates and attempts to solve the Euler-Lagrange equation.
///
/// # Safety
/// The caller must ensure `lagrangian` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[no_mangle]

pub unsafe extern "C" fn rssn_solve_euler_lagrange(
    lagrangian: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    if lagrangian.is_null() || func.is_null() || var.is_null() {

        return std::ptr::null_mut();
    }

    let lagrangian_ref = &*lagrangian;

    let func_str = match c_str_to_str(func) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let var_str = match c_str_to_str(var) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    Box::into_raw(Box::new(
        calculus_of_variations::solve_euler_lagrange(
            lagrangian_ref,
            func_str,
            var_str,
        ),
    ))
}

/// Applies Hamilton's Principle to derive the equations of motion.
///
/// # Safety
/// The caller must ensure `lagrangian` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[no_mangle]

pub unsafe extern "C" fn rssn_hamiltons_principle(
    lagrangian: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    if lagrangian.is_null() || func.is_null() || var.is_null() {

        return std::ptr::null_mut();
    }

    let lagrangian_ref = &*lagrangian;

    let func_str = match c_str_to_str(func) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let var_str = match c_str_to_str(var) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    Box::into_raw(Box::new(
        calculus_of_variations::hamiltons_principle(
            lagrangian_ref,
            func_str,
            var_str,
        ),
    ))
}
