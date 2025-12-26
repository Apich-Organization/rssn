//! Handle-based FFI API for finite field (Galois field) operations.
//!
//! This module provides C-compatible FFI functions for arithmetic in GF(2^8) and
//! general finite fields, including polynomial operations over these fields.

use crate::symbolic::core::Expr;
use crate::symbolic::error_correction_helper::{
    gf256_add,
    gf256_div,
    gf256_exp,
    gf256_inv,
    gf256_log,
    gf256_mul,
    gf256_pow,
    poly_add_gf,
    poly_add_gf256,
    poly_derivative_gf256,
    poly_div_gf,
    poly_div_gf256,
    poly_eval_gf256,
    poly_gcd_gf256,
    poly_mul_gf,
    poly_mul_gf256,
    poly_scale_gf256,
    FiniteField,
};
use std::sync::Arc;

/// Performs addition in GF(2^8) (XOR operation).
#[no_mangle]

pub extern "C" fn rssn_gf256_add(
    a: u8,
    b: u8,
) -> u8 {

    gf256_add(a, b)
}

/// Performs multiplication in GF(2^8).
#[no_mangle]

pub extern "C" fn rssn_gf256_mul(
    a: u8,
    b: u8,
) -> u8 {

    gf256_mul(a, b)
}

/// Computes the exponentiation (anti-logarithm) in GF(2^8).
#[no_mangle]

pub extern "C" fn rssn_gf256_exp(
    log_val: u8
) -> u8 {

    gf256_exp(log_val)
}

/// Computes the discrete logarithm in GF(2^8).
/// Returns 0 if input is 0 (error case, as log(0) is undefined).
#[no_mangle]

pub extern "C" fn rssn_gf256_log(
    a: u8
) -> u8 {

    gf256_log(a).unwrap_or(0)
}

/// Computes a^exp in GF(2^8).
#[no_mangle]

pub extern "C" fn rssn_gf256_pow(
    a: u8,
    exp: u8,
) -> u8 {

    gf256_pow(a, exp)
}

/// Computes the multiplicative inverse in GF(2^8).
/// Returns 0 if input is 0 (error case).
#[no_mangle]

pub extern "C" fn rssn_gf256_inv(
    a: u8
) -> u8 {

    gf256_inv(a).unwrap_or(0)
}

/// Performs division in GF(2^8).
/// Returns 0 if divisor is 0 (error case).
#[no_mangle]

pub extern "C" fn rssn_gf256_div(
    a: u8,
    b: u8,
) -> u8 {

    gf256_div(a, b).unwrap_or(0)
}

/// Evaluates a polynomial over GF(2^8) at point x.
///
/// # Safety
/// Caller must ensure `poly` is a valid pointer to an array of `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_eval_gf256(
    poly: *const u8,
    len: usize,
    x: u8,
) -> u8 {

    if poly.is_null() || len == 0 {

        return 0;
    }

    let slice =
        std::slice::from_raw_parts(
            poly, len,
        );

    poly_eval_gf256(slice, x)
}

/// Adds two polynomials over GF(2^8).
///
/// # Safety
/// Caller must ensure pointers are valid. Result is allocated and must be freed.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_add_gf256(
    p1: *const u8,
    p1_len: usize,
    p2: *const u8,
    p2_len: usize,
    out_len: *mut usize,
) -> *mut u8 {

    if p1.is_null()
        || p2.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    let s1 = std::slice::from_raw_parts(
        p1, p1_len,
    );

    let s2 = std::slice::from_raw_parts(
        p2, p2_len,
    );

    let result = poly_add_gf256(s1, s2);

    *out_len = result.len();

    let boxed =
        result.into_boxed_slice();

    Box::into_raw(boxed) as *mut u8
}

/// Multiplies two polynomials over GF(2^8).
///
/// # Safety
/// Caller must ensure pointers are valid. Result is allocated and must be freed.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_mul_gf256(
    p1: *const u8,
    p1_len: usize,
    p2: *const u8,
    p2_len: usize,
    out_len: *mut usize,
) -> *mut u8 {

    if p1.is_null()
        || p2.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    let s1 = std::slice::from_raw_parts(
        p1, p1_len,
    );

    let s2 = std::slice::from_raw_parts(
        p2, p2_len,
    );

    let result = poly_mul_gf256(s1, s2);

    *out_len = result.len();

    let boxed =
        result.into_boxed_slice();

    Box::into_raw(boxed) as *mut u8
}

/// Scales a polynomial by a constant in GF(2^8).
///
/// # Safety
/// Caller must ensure pointer is valid. Result is allocated and must be freed.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_scale_gf256(
    poly: *const u8,
    len: usize,
    scalar: u8,
    out_len: *mut usize,
) -> *mut u8 {

    if poly.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    let slice =
        std::slice::from_raw_parts(
            poly, len,
        );

    let result =
        poly_scale_gf256(slice, scalar);

    *out_len = result.len();

    let boxed =
        result.into_boxed_slice();

    Box::into_raw(boxed) as *mut u8
}

/// Computes the formal derivative of a polynomial in GF(2^8).
///
/// # Safety
/// Caller must ensure pointer is valid. Result is allocated and must be freed.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_derivative_gf256(
    poly: *const u8,
    len: usize,
    out_len: *mut usize,
) -> *mut u8 {

    if poly.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    let slice =
        std::slice::from_raw_parts(
            poly, len,
        );

    let result =
        poly_derivative_gf256(slice);

    *out_len = result.len();

    let boxed =
        result.into_boxed_slice();

    Box::into_raw(boxed) as *mut u8
}

/// Computes the GCD of two polynomials over GF(2^8).
///
/// # Safety
/// Caller must ensure pointers are valid. Result is allocated and must be freed.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_gcd_gf256(
    p1: *const u8,
    p1_len: usize,
    p2: *const u8,
    p2_len: usize,
    out_len: *mut usize,
) -> *mut u8 {

    if p1.is_null()
        || p2.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    let s1 = std::slice::from_raw_parts(
        p1, p1_len,
    );

    let s2 = std::slice::from_raw_parts(
        p2, p2_len,
    );

    let result = poly_gcd_gf256(s1, s2);

    *out_len = result.len();

    let boxed =
        result.into_boxed_slice();

    Box::into_raw(boxed) as *mut u8
}

/// Creates a new finite field GF(modulus).
///
/// Returns an opaque handle to the field.
#[no_mangle]

pub extern "C" fn rssn_finite_field_new(
    modulus: i64
) -> *mut Arc<FiniteField> {

    Box::into_raw(Box::new(
        FiniteField::new(modulus),
    ))
}

/// Frees a finite field handle.
///
/// # Safety
/// Caller must ensure `field` is a valid pointer returned by `rssn_finite_field_new`.
#[no_mangle]

pub unsafe extern "C" fn rssn_finite_field_free(
    field: *mut Arc<FiniteField>
) {

    if !field.is_null() {

        drop(Box::from_raw(field));
    }
}

/// Adds two polynomials over a general finite field.
///
/// # Safety
/// Caller must ensure all pointers are valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_add_gf(
    p1: *const Expr,
    p2: *const Expr,
    field: *const Arc<FiniteField>,
) -> *mut Expr {

    if p1.is_null()
        || p2.is_null()
        || field.is_null()
    {

        return std::ptr::null_mut();
    }

    match poly_add_gf(
        &*p1, &*p2, &*field,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Multiplies two polynomials over a general finite field.
///
/// # Safety
/// Caller must ensure all pointers are valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_poly_mul_gf(
    p1: *const Expr,
    p2: *const Expr,
    field: *const Arc<FiniteField>,
) -> *mut Expr {

    if p1.is_null()
        || p2.is_null()
        || field.is_null()
    {

        return std::ptr::null_mut();
    }

    match poly_mul_gf(
        &*p1, &*p2, &*field,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}
