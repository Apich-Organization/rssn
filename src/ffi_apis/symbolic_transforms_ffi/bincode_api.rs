//! Bincode-based FFI API for symbolic integral transforms.
//!
//! This module provides binary serialization-based FFI functions for Fourier, Laplace, and Z-transforms,
//! offering efficient binary data interchange for high-performance applications.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::transforms;

// --- Fourier Transform ---

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::fourier_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_inverse_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_time_shift(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_omega_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::fourier_time_shift(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_frequency_shift(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_omega_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::fourier_frequency_shift(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_scaling(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_omega_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::fourier_scaling(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_differentiation(
    f_omega_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_omega_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_bincode_buffer(&transforms::fourier_differentiation(&f, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

// --- Laplace Transform ---

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::laplace_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_inverse_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_time_shift(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::laplace_time_shift(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_frequency_shift(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::laplace_frequency_shift(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_scaling(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::laplace_scaling(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_differentiation(
    f_s_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
    f_zero_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_s_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    let f_zero: Option<Expr> = from_bincode_buffer(&f_zero_buf);

    if let (Some(f), Some(ov), Some(fz)) = (f, out_var, f_zero) {

        to_bincode_buffer(&transforms::laplace_differentiation(&f, &ov, &fz))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laplace_integration(
    f_s_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_s_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_bincode_buffer(&transforms::laplace_integration(&f, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

// --- Z-Transform ---

#[no_mangle]

pub extern "C" fn rssn_bincode_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::z_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_inverse_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(e), Some(iv), Some(ov)) = (
        expr, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::inverse_z_transform(&e, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_z_time_shift(
    f_z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_z_buf);

    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(k), Some(ov)) = (f, k, out_var) {

        to_bincode_buffer(&transforms::z_time_shift(&f, &k, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_z_scaling(
    f_z_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_z_buf);

    let a: Option<Expr> = from_bincode_buffer(&a_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(a), Some(ov)) = (f, a, out_var) {

        to_bincode_buffer(&transforms::z_scaling(&f, &a, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_z_differentiation(
    f_z_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_z_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(ov)) = (f, out_var) {

        to_bincode_buffer(&transforms::z_differentiation(&f, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

// --- Utils ---

#[no_mangle]

pub extern "C" fn rssn_bincode_convolution_fourier(
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_buf);

    let g: Option<Expr> = from_bincode_buffer(&g_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(g), Some(iv), Some(ov)) = (
        f, g, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::convolution_fourier(&f, &g, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_convolution_laplace(
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> = from_bincode_buffer(&f_buf);

    let g: Option<Expr> = from_bincode_buffer(&g_buf);

    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);

    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);

    if let (Some(f), Some(g), Some(iv), Some(ov)) = (
        f, g, in_var, out_var,
    ) {

        to_bincode_buffer(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_partial_fraction_decomposition(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let var: Option<String> = from_bincode_buffer(&var_buf);

    if let (Some(expr), Some(var)) = (expr, var) {

        if let Some(result) = transforms::partial_fraction_decomposition(&expr, &var) {

            return to_bincode_buffer(&result);
        }
    }

    BincodeBuffer::empty()
}
