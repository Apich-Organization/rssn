//! Bincode-based FFI API for cryptographic operations.
//!
//! This module provides binary serialization-based FFI functions for elliptic curve cryptography,
//! offering efficient binary data interchange for high-performance applications.

use crate::symbolic::cryptography::{
    EllipticCurve, CurvePoint,
};
use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use crate::ffi_apis::common::*;
use num_bigint::BigInt;
use std::sync::Arc;

/// Performs point addition on elliptic curve via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_add(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer,
    p1_x_buf: BincodeBuffer,
    p1_y_buf: BincodeBuffer,
    p2_x_buf: BincodeBuffer,
    p2_y_buf: BincodeBuffer
) -> BincodeBuffer {
    let a: Option<i64> = from_bincode_buffer(&a_buf);
    let b: Option<i64> = from_bincode_buffer(&b_buf);
    let modulus: Option<i64> = from_bincode_buffer(&modulus_buf);
    let p1_x: Option<i64> = from_bincode_buffer(&p1_x_buf);
    let p1_y: Option<i64> = from_bincode_buffer(&p1_y_buf);
    let p2_x: Option<i64> = from_bincode_buffer(&p2_x_buf);
    let p2_y: Option<i64> = from_bincode_buffer(&p2_y_buf);
    
    if let (Some(a), Some(b), Some(m), Some(x1), Some(y1), Some(x2), Some(y2)) = 
        (a, b, modulus, p1_x, p1_y, p2_x, p2_y) {
        let field = PrimeField::new(BigInt::from(m));
        let curve = EllipticCurve {
            a: PrimeFieldElement::new(BigInt::from(a), field.clone()),
            b: PrimeFieldElement::new(BigInt::from(b), field.clone()),
            field: field.clone(),
        };
        let p1 = CurvePoint::Affine {
            x: PrimeFieldElement::new(BigInt::from(x1), field.clone()),
            y: PrimeFieldElement::new(BigInt::from(y1), field.clone()),
        };
        let p2 = CurvePoint::Affine {
            x: PrimeFieldElement::new(BigInt::from(x2), field.clone()),
            y: PrimeFieldElement::new(BigInt::from(y2), field),
        };
        let result = curve.add(&p1, &p2);
        match result {
            CurvePoint::Infinity => to_bincode_buffer(&(true, 0i64, 0i64)), // (is_infinity, x, y)
            CurvePoint::Affine { x, y } => {
                // Convert to i64 for simplicity (assuming small modulus)
                let x_val: i64 = x.value.to_string().parse().unwrap_or(0);
                let y_val: i64 = y.value.to_string().parse().unwrap_or(0);
                to_bincode_buffer(&(false, x_val, y_val))
            }
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Performs scalar multiplication on elliptic curve via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_scalar_mult(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    p_x_buf: BincodeBuffer,
    p_y_buf: BincodeBuffer
) -> BincodeBuffer {
    let a: Option<i64> = from_bincode_buffer(&a_buf);
    let b: Option<i64> = from_bincode_buffer(&b_buf);
    let modulus: Option<i64> = from_bincode_buffer(&modulus_buf);
    let k: Option<i64> = from_bincode_buffer(&k_buf);
    let p_x: Option<i64> = from_bincode_buffer(&p_x_buf);
    let p_y: Option<i64> = from_bincode_buffer(&p_y_buf);
    
    if let (Some(a), Some(b), Some(m), Some(k), Some(px), Some(py)) = 
        (a, b, modulus, k, p_x, p_y) {
        let field = PrimeField::new(BigInt::from(m));
        let curve = EllipticCurve {
            a: PrimeFieldElement::new(BigInt::from(a), field.clone()),
            b: PrimeFieldElement::new(BigInt::from(b), field.clone()),
            field: field.clone(),
        };
        let p = CurvePoint::Affine {
            x: PrimeFieldElement::new(BigInt::from(px), field.clone()),
            y: PrimeFieldElement::new(BigInt::from(py), field),
        };
        let result = curve.scalar_mult(&BigInt::from(k), &p);
        match result {
            CurvePoint::Infinity => to_bincode_buffer(&(true, 0i64, 0i64)),
            CurvePoint::Affine { x, y } => {
                let x_val: i64 = x.value.to_string().parse().unwrap_or(0);
                let y_val: i64 = y.value.to_string().parse().unwrap_or(0);
                to_bincode_buffer(&(false, x_val, y_val))
            }
        }
    } else {
        BincodeBuffer::empty()
    }
}
