//! Handle-based FFI API for cryptographic operations.
//!
//! This module provides C-compatible FFI functions for elliptic curve cryptography (ECC),
//! including point addition, scalar multiplication, key pair generation, and ECDH shared
//! secret derivation.

use crate::symbolic::cryptography::{
    EllipticCurve, CurvePoint, EcdhKeyPair, generate_keypair, generate_shared_secret,
};
use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use num_bigint::BigInt;
use std::sync::Arc;

/// Creates a new elliptic curve over a prime field.
///
/// # Safety
/// All BigInt pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_elliptic_curve_new(
    a: i64,
    b: i64,
    modulus: i64
) -> *mut EllipticCurve {
    let field = PrimeField::new(BigInt::from(modulus));
    let curve = EllipticCurve {
        a: PrimeFieldElement::new(BigInt::from(a), field.clone()),
        b: PrimeFieldElement::new(BigInt::from(b), field.clone()),
        field,
    };
    Box::into_raw(Box::new(curve))
}

/// Frees an elliptic curve handle.
///
/// # Safety
/// Caller must ensure `curve` was returned by `rssn_elliptic_curve_new`.
#[no_mangle]
pub unsafe extern "C" fn rssn_elliptic_curve_free(curve: *mut EllipticCurve) {
    if !curve.is_null() {
        drop(Box::from_raw(curve));
    }
}

/// Creates an affine curve point.
///
/// # Safety
/// `field` must match the curve's field.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_point_affine(
    x: i64,
    y: i64,
    modulus: i64
) -> *mut CurvePoint {
    let field = PrimeField::new(BigInt::from(modulus));
    let point = CurvePoint::Affine {
        x: PrimeFieldElement::new(BigInt::from(x), field.clone()),
        y: PrimeFieldElement::new(BigInt::from(y), field),
    };
    Box::into_raw(Box::new(point))
}

/// Creates the point at infinity.
#[no_mangle]
pub extern "C" fn rssn_curve_point_infinity() -> *mut CurvePoint {
    Box::into_raw(Box::new(CurvePoint::Infinity))
}

/// Frees a curve point handle.
///
/// # Safety
/// Caller must ensure `point` was returned by a curve point function.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_point_free(point: *mut CurvePoint) {
    if !point.is_null() {
        drop(Box::from_raw(point));
    }
}

/// Adds two curve points.
///
/// # Safety
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_add(
    curve: *const EllipticCurve,
    p1: *const CurvePoint,
    p2: *const CurvePoint
) -> *mut CurvePoint {
    if curve.is_null() || p1.is_null() || p2.is_null() {
        return std::ptr::null_mut();
    }
    let result = (*curve).add(&*p1, &*p2);
    Box::into_raw(Box::new(result))
}

/// Performs scalar multiplication k * P on elliptic curve.
///
/// # Safety
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_scalar_mult(
    curve: *const EllipticCurve,
    k: i64,
    p: *const CurvePoint
) -> *mut CurvePoint {
    if curve.is_null() || p.is_null() {
        return std::ptr::null_mut();
    }
    let result = (*curve).scalar_mult(&BigInt::from(k), &*p);
    Box::into_raw(Box::new(result))
}

/// Generates an ECDH key pair.
///
/// # Safety
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_generate_keypair(
    curve: *const EllipticCurve,
    generator: *const CurvePoint
) -> *mut EcdhKeyPair {
    if curve.is_null() || generator.is_null() {
        return std::ptr::null_mut();
    }
    let keypair = generate_keypair(&*curve, &*generator);
    Box::into_raw(Box::new(keypair))
}

/// Frees an ECDH key pair.
///
/// # Safety
/// Caller must ensure `keypair` was returned by `rssn_generate_keypair`.
#[no_mangle]
pub unsafe extern "C" fn rssn_keypair_free(keypair: *mut EcdhKeyPair) {
    if !keypair.is_null() {
        drop(Box::from_raw(keypair));
    }
}

/// Generates a shared secret using ECDH.
///
/// # Safety
/// All pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_generate_shared_secret(
    curve: *const EllipticCurve,
    private_key: i64,
    other_public_key: *const CurvePoint
) -> *mut CurvePoint {
    if curve.is_null() || other_public_key.is_null() {
        return std::ptr::null_mut();
    }
    let result = generate_shared_secret(&*curve, &BigInt::from(private_key), &*other_public_key);
    Box::into_raw(Box::new(result))
}
