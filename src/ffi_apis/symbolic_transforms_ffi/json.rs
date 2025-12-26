//! JSON-based FFI API for symbolic integral transforms.
//!
//! This module provides JSON string-based FFI functions for Fourier, Laplace, and Z-transforms,
//! allowing for language-agnostic integration through JSON serialization/deserialization.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::transforms;

// --- Fourier Transform ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_fourier_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::fourier_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_inverse_fourier_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_fourier_time_shift(
    f_omega_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_omega_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::fourier_time_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_fourier_frequency_shift(
    f_omega_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_omega_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::fourier_frequency_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_fourier_scaling(
    f_omega_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_omega_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::fourier_scaling(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_fourier_differentiation(
    f_omega_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_omega_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_json_string(&transforms::fourier_differentiation(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Laplace Transform ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::laplace_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_inverse_laplace_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_time_shift(
    f_s_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_s_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::laplace_time_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_frequency_shift(
    f_s_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_s_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::laplace_frequency_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_scaling(
    f_s_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_s_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::laplace_scaling(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_differentiation(
    f_s_json : *const c_char,
    out_var_json : *const c_char,
    f_zero_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_s_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    let f_zero : Option<Expr> = from_json_string(f_zero_json);

    if let (Some(f), Some(ov), Some(fz)) = (f, out_var, f_zero) {

        to_json_string(&transforms::laplace_differentiation(&f, &ov, &fz))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_laplace_integration(
    f_s_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_s_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_json_string(&transforms::laplace_integration(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Z-Transform ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_z_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::z_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_inverse_z_transform(
    expr_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_z_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_z_time_shift(
    f_z_json : *const c_char,
    k_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_z_json);

    let k : Option<Expr> = from_json_string(k_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(k), Some(ov)) = (f, k, out_var) {

        to_json_string(&transforms::z_time_shift(&f, &k, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_z_scaling(
    f_z_json : *const c_char,
    a_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_z_json);

    let a : Option<Expr> = from_json_string(a_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_json_string(&transforms::z_scaling(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_z_differentiation(
    f_z_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_z_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_json_string(&transforms::z_differentiation(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Utils ---

#[no_mangle]

pub unsafe extern "C" fn rssn_json_convolution_fourier(
    f_json : *const c_char,
    g_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_json);

    let g : Option<Expr> = from_json_string(g_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(g), Some(iv), Some(ov)) = (
        f,
        g,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::convolution_fourier(&f, &g, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_convolution_laplace(
    f_json : *const c_char,
    g_json : *const c_char,
    in_var_json : *const c_char,
    out_var_json : *const c_char,
) -> *mut c_char {

    let f : Option<Expr> = from_json_string(f_json);

    let g : Option<Expr> = from_json_string(g_json);

    let in_var : Option<String> = from_json_string(in_var_json);

    let out_var : Option<String> = from_json_string(out_var_json);

    if let (Some(f), Some(g), Some(iv), Some(ov)) = (
        f,
        g,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_partial_fraction_decomposition(
    expr_json : *const c_char,
    var_json : *const c_char,
) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(expr_json);

    let var : Option<String> = from_json_string(var_json);

    if let (Some(expr), Some(var)) = (expr, var) {

        if let Some(result) = transforms::partial_fraction_decomposition(&expr, &var) {

            return to_json_string(&result);
        }
    }

    std::ptr::null_mut()
}
