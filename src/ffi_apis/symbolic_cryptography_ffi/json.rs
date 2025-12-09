//! JSON-based FFI API for cryptographic operations.
//!
//! This module provides JSON string-based FFI functions for elliptic curve cryptography,
//! enabling language-agnostic integration for ECC operations.

use crate::symbolic::cryptography::{
    EllipticCurve, CurvePoint, generate_keypair, generate_shared_secret,
};
use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use crate::ffi_apis::common::*;
use num_bigint::BigInt;
use std::os::raw::c_char;
use std::sync::Arc;

/// Creates an elliptic curve and performs point addition via JSON interface.
/// Input: {"a": int, "b": int, "modulus": int, "p1": {x, y}, "p2": {x, y}}
#[no_mangle]
pub unsafe extern "C" fn rssn_json_curve_add(
    a_json: *const c_char,
    b_json: *const c_char,
    modulus_json: *const c_char,
    p1_x_json: *const c_char,
    p1_y_json: *const c_char,
    p2_x_json: *const c_char,
    p2_y_json: *const c_char
) -> *mut c_char {
    let a: Option<i64> = from_json_string(a_json);
    let b: Option<i64> = from_json_string(b_json);
    let modulus: Option<i64> = from_json_string(modulus_json);
    let p1_x: Option<i64> = from_json_string(p1_x_json);
    let p1_y: Option<i64> = from_json_string(p1_y_json);
    let p2_x: Option<i64> = from_json_string(p2_x_json);
    let p2_y: Option<i64> = from_json_string(p2_y_json);
    
    if let (Some(a), Some(b), Some(m), Some(x1), Some(y1), Some(x2), Some(y2)) = 
        (a, b, modulus, p1_x, p1_y, p2_x, p2_y) {
        let field = Arc::new(PrimeField::new(BigInt::from(m)));
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
            CurvePoint::Infinity => to_json_string(&"infinity"),
            CurvePoint::Affine { x, y } => {
                let obj = serde_json::json!({
                    "x": x.value.to_string(),
                    "y": y.value.to_string()
                });
                to_json_string(&obj)
            }
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Performs scalar multiplication via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_curve_scalar_mult(
    a_json: *const c_char,
    b_json: *const c_char,
    modulus_json: *const c_char,
    k_json: *const c_char,
    p_x_json: *const c_char,
    p_y_json: *const c_char
) -> *mut c_char {
    let a: Option<i64> = from_json_string(a_json);
    let b: Option<i64> = from_json_string(b_json);
    let modulus: Option<i64> = from_json_string(modulus_json);
    let k: Option<i64> = from_json_string(k_json);
    let p_x: Option<i64> = from_json_string(p_x_json);
    let p_y: Option<i64> = from_json_string(p_y_json);
    
    if let (Some(a), Some(b), Some(m), Some(k), Some(px), Some(py)) = 
        (a, b, modulus, k, p_x, p_y) {
        let field = Arc::new(PrimeField::new(BigInt::from(m)));
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
            CurvePoint::Infinity => to_json_string(&"infinity"),
            CurvePoint::Affine { x, y } => {
                let obj = serde_json::json!({
                    "x": x.value.to_string(),
                    "y": y.value.to_string()
                });
                to_json_string(&obj)
            }
        }
    } else {
        std::ptr::null_mut()
    }
}
