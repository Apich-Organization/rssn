//! Bincode-based FFI API for symbolic integral transforms.
//!
//! This module provides binary serialization-based FFI functions for Fourier, Laplace, and Z-transforms,
//! offering efficient binary data interchange for high-performance applications.

use crate::symbolic::core::Expr;
use crate::symbolic::transforms;
use crate::ffi_apis::common::*;

/// Computes the symbolic Fourier transform via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::fourier_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_inverse_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::laplace_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_inverse_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::z_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_inverse_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(e), Some(iv), Some(ov)) = (expr, in_var, out_var) {
        to_bincode_buffer(&transforms::inverse_z_transform(&e, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_convolution_fourier(
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let f: Option<Expr> = from_bincode_buffer(&f_buf);
    let g: Option<Expr> = from_bincode_buffer(&g_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(f), Some(g), Some(iv), Some(ov)) = (f, g, in_var, out_var) {
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
    out_var_buf: BincodeBuffer
) -> BincodeBuffer {
    let f: Option<Expr> = from_bincode_buffer(&f_buf);
    let g: Option<Expr> = from_bincode_buffer(&g_buf);
    let in_var: Option<String> = from_bincode_buffer(&in_var_buf);
    let out_var: Option<String> = from_bincode_buffer(&out_var_buf);
    if let (Some(f), Some(g), Some(iv), Some(ov)) = (f, g, in_var, out_var) {
        to_bincode_buffer(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } else {
        BincodeBuffer::empty()
    }
}
