//! Handle-based FFI API for cryptographic operations.
//!
//! This module provides C-compatible FFI functions for elliptic curve cryptography (ECC),
//! including point addition, scalar multiplication, key pair generation, ECDH shared
//! secret derivation, and ECDSA digital signatures.

use crate::symbolic::cryptography::{
    EllipticCurve, CurvePoint, EcdhKeyPair, EcdsaSignature,
    generate_keypair, generate_shared_secret, point_compress, point_decompress,
    ecdsa_sign, ecdsa_verify,
};
use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use num_bigint::BigInt;
use std::sync::Arc;

/// Creates a new elliptic curve over a prime field using the new() constructor.
#[no_mangle]
pub extern "C" fn rssn_elliptic_curve_new(a: i64, b: i64, modulus: i64) -> *mut EllipticCurve {
    let curve = EllipticCurve::new(BigInt::from(a), BigInt::from(b), BigInt::from(modulus));
    Box::into_raw(Box::new(curve))
}

/// Frees an elliptic curve handle.
#[no_mangle]
pub unsafe extern "C" fn rssn_elliptic_curve_free(curve: *mut EllipticCurve) {
    if !curve.is_null() {
        drop(Box::from_raw(curve));
    }
}

/// Creates an affine curve point.
#[no_mangle]
pub extern "C" fn rssn_curve_point_affine(x: i64, y: i64, modulus: i64) -> *mut CurvePoint {
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

/// Checks if a point is the point at infinity.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_point_is_infinity(point: *const CurvePoint) -> bool {
    if point.is_null() { return false; }
    (*point).is_infinity()
}

/// Frees a curve point handle.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_point_free(point: *mut CurvePoint) {
    if !point.is_null() {
        drop(Box::from_raw(point));
    }
}

/// Checks if a point is on the curve.
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_is_on_curve(
    curve: *const EllipticCurve,
    point: *const CurvePoint
) -> bool {
    if curve.is_null() || point.is_null() { return false; }
    (*curve).is_on_curve(&*point)
}

/// Negates a point on the curve (P -> -P).
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_negate(
    curve: *const EllipticCurve,
    point: *const CurvePoint
) -> *mut CurvePoint {
    if curve.is_null() || point.is_null() { return std::ptr::null_mut(); }
    let result = (*curve).negate(&*point);
    Box::into_raw(Box::new(result))
}

/// Doubles a point on the curve (2P).
#[no_mangle]
pub unsafe extern "C" fn rssn_curve_double(
    curve: *const EllipticCurve,
    point: *const CurvePoint
) -> *mut CurvePoint {
    if curve.is_null() || point.is_null() { return std::ptr::null_mut(); }
    let result = (*curve).double(&*point);
    Box::into_raw(Box::new(result))
}

/// Adds two curve points.
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
#[no_mangle]
pub unsafe extern "C" fn rssn_keypair_free(keypair: *mut EcdhKeyPair) {
    if !keypair.is_null() {
        drop(Box::from_raw(keypair));
    }
}

/// Generates a shared secret using ECDH.
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

/// Signs a message using ECDSA.
#[no_mangle]
pub unsafe extern "C" fn rssn_ecdsa_sign(
    message_hash: i64,
    private_key: i64,
    curve: *const EllipticCurve,
    generator: *const CurvePoint,
    order: i64
) -> *mut EcdsaSignature {
    if curve.is_null() || generator.is_null() {
        return std::ptr::null_mut();
    }
    match ecdsa_sign(
        &BigInt::from(message_hash),
        &BigInt::from(private_key),
        &*curve,
        &*generator,
        &BigInt::from(order)
    ) {
        Some(sig) => Box::into_raw(Box::new(sig)),
        None => std::ptr::null_mut(),
    }
}

/// Verifies an ECDSA signature.
#[no_mangle]
pub unsafe extern "C" fn rssn_ecdsa_verify(
    message_hash: i64,
    signature: *const EcdsaSignature,
    public_key: *const CurvePoint,
    curve: *const EllipticCurve,
    generator: *const CurvePoint,
    order: i64
) -> bool {
    if signature.is_null() || public_key.is_null() || curve.is_null() || generator.is_null() {
        return false;
    }
    ecdsa_verify(
        &BigInt::from(message_hash),
        &*signature,
        &*public_key,
        &*curve,
        &*generator,
        &BigInt::from(order)
    )
}

/// Frees an ECDSA signature.
#[no_mangle]
pub unsafe extern "C" fn rssn_ecdsa_signature_free(sig: *mut EcdsaSignature) {
    if !sig.is_null() {
        drop(Box::from_raw(sig));
    }
}

