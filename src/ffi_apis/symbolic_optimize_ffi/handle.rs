use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::os::raw::c_int;

use crate::symbolic::core::Expr;
use crate::symbolic::optimize::*;

unsafe fn parse_c_str_array(
    arr: *const *const c_char,
    len: usize,
) -> Option<Vec<String>> {

    if arr.is_null() {

        return None;
    }

    let mut vars =
        Vec::with_capacity(len);

    for i in 0..len {

        let ptr = *arr.add(i);

        if ptr.is_null() {

            return None;
        }

        let c_str = CStr::from_ptr(ptr);

        match c_str.to_str() {
            | Ok(s) => {
                vars.push(s.to_string())
            },
            | Err(_) => return None,
        }
    }

    Some(vars)
}

/// Finds extrema of a function (Handle)
#[no_mangle]

pub extern "C" fn rssn_find_extrema_handle(
    expr_ptr: *const Expr,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut Vec<CriticalPoint> {

    if expr_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let expr = &*expr_ptr;

        let vars_strings = match parse_c_str_array(
            vars_ptr,
            vars_len as usize,
        ) {
            | Some(v) => v,
            | None => return std::ptr::null_mut(),
        };

        let vars_refs: Vec<&str> =
            vars_strings
                .iter()
                .map(|s| s.as_str())
                .collect();

        match find_extrema(
            expr,
            &vars_refs,
        ) {
            | Ok(points) => {
                Box::into_raw(Box::new(
                    points,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Computes Hessian matrix (Handle)
#[no_mangle]

pub extern "C" fn rssn_hessian_matrix_handle(
    expr_ptr: *const Expr,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut Expr {

    if expr_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let expr = &*expr_ptr;

        let vars_strings = match parse_c_str_array(
            vars_ptr,
            vars_len as usize,
        ) {
            | Some(v) => v,
            | None => return std::ptr::null_mut(),
        };

        let vars_refs: Vec<&str> =
            vars_strings
                .iter()
                .map(|s| s.as_str())
                .collect();

        let hessian = hessian_matrix(
            expr,
            &vars_refs,
        );

        Box::into_raw(Box::new(hessian))
    }
}

/// Finds constrained extrema (Handle)
#[no_mangle]

pub extern "C" fn rssn_find_constrained_extrema_handle(
    expr_ptr: *const Expr,
    constraints_ptr: *const Vec<Expr>,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut Vec<HashMap<Expr, Expr>> {

    if expr_ptr.is_null()
        || constraints_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let expr = &*expr_ptr;

        let constraints =
            &*constraints_ptr;

        let vars_strings = match parse_c_str_array(
            vars_ptr,
            vars_len as usize,
        ) {
            | Some(v) => v,
            | None => return std::ptr::null_mut(),
        };

        let vars_refs: Vec<&str> =
            vars_strings
                .iter()
                .map(|s| s.as_str())
                .collect();

        match find_constrained_extrema(
            expr,
            constraints,
            &vars_refs,
        ) {
            | Ok(solutions) => {
                Box::into_raw(Box::new(
                    solutions,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Frees a Vec<CriticalPoint> handle
#[no_mangle]

pub extern "C" fn rssn_free_critical_point_vec_handle(
    ptr: *mut Vec<CriticalPoint>
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Frees a Vec<HashMap<Expr, Expr>> handle
#[no_mangle]

pub extern "C" fn rssn_free_solution_vec_handle(
    ptr: *mut Vec<HashMap<Expr, Expr>>
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}
