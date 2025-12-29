//! JSON-based FFI API for finite field (Galois field) operations.
//!
//! This module provides JSON string-based FFI functions for GF(2^8) and
//! general finite field arithmetic operations.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::error_correction_helper::gf256_add;
use crate::symbolic::error_correction_helper::gf256_inv;
use crate::symbolic::error_correction_helper::gf256_mul;
use crate::symbolic::error_correction_helper::poly_add_gf;
use crate::symbolic::error_correction_helper::poly_add_gf256;
use crate::symbolic::error_correction_helper::poly_eval_gf256;
use crate::symbolic::error_correction_helper::poly_mul_gf;
use crate::symbolic::error_correction_helper::poly_mul_gf256;
use crate::symbolic::error_correction_helper::FiniteField;

/// Performs addition in GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gf256_add(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<u8> =
        from_json_string(a_json);

    let b: Option<u8> =
        from_json_string(b_json);

    if let (Some(va), Some(vb)) = (a, b)
    {

        to_json_string(&gf256_add(
            va, vb,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Performs multiplication in GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gf256_mul(
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let a: Option<u8> =
        from_json_string(a_json);

    let b: Option<u8> =
        from_json_string(b_json);

    if let (Some(va), Some(vb)) = (a, b)
    {

        to_json_string(&gf256_mul(
            va, vb,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes inverse in GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gf256_inv(
    a_json: *const c_char
) -> *mut c_char {

    let a: Option<u8> =
        from_json_string(a_json);

    if let Some(va) = a {

        match gf256_inv(va) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Evaluates a polynomial over GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_poly_eval_gf256(
    poly_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let poly: Option<Vec<u8>> =
        from_json_string(poly_json);

    let x: Option<u8> =
        from_json_string(x_json);

    if let (Some(p), Some(vx)) =
        (poly, x)
    {

        to_json_string(
            &poly_eval_gf256(&p, vx),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Adds two polynomials over GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_poly_add_gf256(
    p1_json: *const c_char,
    p2_json: *const c_char,
) -> *mut c_char {

    let p1: Option<Vec<u8>> =
        from_json_string(p1_json);

    let p2: Option<Vec<u8>> =
        from_json_string(p2_json);

    if let (Some(v1), Some(v2)) =
        (p1, p2)
    {

        to_json_string(&poly_add_gf256(
            &v1, &v2,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Multiplies two polynomials over GF(2^8) via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_poly_mul_gf256(
    p1_json: *const c_char,
    p2_json: *const c_char,
) -> *mut c_char {

    let p1: Option<Vec<u8>> =
        from_json_string(p1_json);

    let p2: Option<Vec<u8>> =
        from_json_string(p2_json);

    if let (Some(v1), Some(v2)) =
        (p1, p2)
    {

        to_json_string(&poly_mul_gf256(
            &v1, &v2,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Adds two polynomials over a general finite field via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_poly_add_gf(
    p1_json: *const c_char,
    p2_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let p1: Option<Expr> =
        from_json_string(p1_json);

    let p2: Option<Expr> =
        from_json_string(p2_json);

    let modulus: Option<i64> =
        from_json_string(modulus_json);

    if let (
        Some(v1),
        Some(v2),
        Some(m),
    ) = (p1, p2, modulus)
    {

        let field = FiniteField::new(m);

        match poly_add_gf(
            &v1,
            &v2,
            &field,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Multiplies two polynomials over a general finite field via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_poly_mul_gf(
    p1_json: *const c_char,
    p2_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let p1: Option<Expr> =
        from_json_string(p1_json);

    let p2: Option<Expr> =
        from_json_string(p2_json);

    let modulus: Option<i64> =
        from_json_string(modulus_json);

    if let (
        Some(v1),
        Some(v2),
        Some(m),
    ) = (p1, p2, modulus)
    {

        let field = FiniteField::new(m);

        match poly_mul_gf(
            &v1,
            &v2,
            &field,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}
