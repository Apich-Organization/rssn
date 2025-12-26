//! Bincode-based FFI API for finite field (Galois field) operations.
//!
//! This module provides binary serialization-based FFI functions for GF(2^8) and
//! general finite field arithmetic operations.

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

/// Performs addition in GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_gf256_add(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<u8> =
        from_bincode_buffer(&a_buf);

    let b: Option<u8> =
        from_bincode_buffer(&b_buf);

    if let (Some(va), Some(vb)) = (a, b)
    {

        to_bincode_buffer(&gf256_add(
            va, vb,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Performs multiplication in GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_gf256_mul(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<u8> =
        from_bincode_buffer(&a_buf);

    let b: Option<u8> =
        from_bincode_buffer(&b_buf);

    if let (Some(va), Some(vb)) = (a, b)
    {

        to_bincode_buffer(&gf256_mul(
            va, vb,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes inverse in GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_gf256_inv(
    a_buf: BincodeBuffer
) -> BincodeBuffer {

    let a: Option<u8> =
        from_bincode_buffer(&a_buf);

    if let Some(va) = a {

        match gf256_inv(va) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Evaluates a polynomial over GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_eval_gf256(
    poly_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let poly: Option<Vec<u8>> =
        from_bincode_buffer(&poly_buf);

    let x: Option<u8> =
        from_bincode_buffer(&x_buf);

    if let (Some(p), Some(vx)) =
        (poly, x)
    {

        to_bincode_buffer(
            &poly_eval_gf256(&p, vx),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Adds two polynomials over GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_add_gf256(
    p1_buf: BincodeBuffer,
    p2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p1: Option<Vec<u8>> =
        from_bincode_buffer(&p1_buf);

    let p2: Option<Vec<u8>> =
        from_bincode_buffer(&p2_buf);

    if let (Some(v1), Some(v2)) =
        (p1, p2)
    {

        to_bincode_buffer(
            &poly_add_gf256(&v1, &v2),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Multiplies two polynomials over GF(2^8) via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_mul_gf256(
    p1_buf: BincodeBuffer,
    p2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p1: Option<Vec<u8>> =
        from_bincode_buffer(&p1_buf);

    let p2: Option<Vec<u8>> =
        from_bincode_buffer(&p2_buf);

    if let (Some(v1), Some(v2)) =
        (p1, p2)
    {

        to_bincode_buffer(
            &poly_mul_gf256(&v1, &v2),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Adds two polynomials over a general finite field via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_add_gf(
    p1_buf: BincodeBuffer,
    p2_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p1: Option<Expr> =
        from_bincode_buffer(&p1_buf);

    let p2: Option<Expr> =
        from_bincode_buffer(&p2_buf);

    let modulus: Option<i64> =
        from_bincode_buffer(
            &modulus_buf,
        );

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
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Multiplies two polynomials over a general finite field via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_poly_mul_gf(
    p1_buf: BincodeBuffer,
    p2_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p1: Option<Expr> =
        from_bincode_buffer(&p1_buf);

    let p2: Option<Expr> =
        from_bincode_buffer(&p2_buf);

    let modulus: Option<i64> =
        from_bincode_buffer(
            &modulus_buf,
        );

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
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}
