//! JSON-based FFI API for symbolic integral transforms.
//!
//! This module provides JSON string-based FFI functions for Fourier, Laplace, and Z-transforms,
//! allowing for language-agnostic integration through JSON serialization/deserialization.

use crate::symbolic::core::Expr;
use crate::symbolic::transforms;
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

/// Computes the symbolic Fourier transform via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_fourier_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::fourier_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_inverse_fourier_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_laplace_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::laplace_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_inverse_laplace_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_z_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::z_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_inverse_z_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_json_string(&transforms::inverse_z_transform(&e, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_convolution_fourier(
    f_json: *const c_char,
    g_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let f: Option<Expr> = from_json_string(f_json);
    let g: Option<Expr> = from_json_string(g_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(f), Some(g), Some(iv), Some(ov)) = (f, g, in_var, out_var) {
        to_json_string(&transforms::convolution_fourier(&f, &g, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_convolution_laplace(
    f_json: *const c_char,
    g_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char
) -> *mut c_char {
    let f: Option<Expr> = from_json_string(f_json);
    let g: Option<Expr> = from_json_string(g_json);
    let in_var: Option<String> = from_json_string(in_var_json);
    let out_var: Option<String> = from_json_string(out_var_json);
    if let (Some(f), Some(g), Some(iv), Some(ov)) = (f, g, in_var, out_var) {
        to_json_string(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } else {
        std::ptr::null_mut()
    }
}
