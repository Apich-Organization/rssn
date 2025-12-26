//! Handle-based FFI API for cryptographic operations.
//!
//! This module provides C-compatible FFI functions for elliptic curve cryptography (ECC).
//! Large integers are passed as strings (decimal).

use crate::ffi_apis::common::{
    c_str_to_str,
    to_c_string,
};
use crate::symbolic::cryptography::{
    ecdsa_sign,
    ecdsa_verify,
    generate_keypair,
    generate_shared_secret,
    point_compress,
    point_decompress,
    CurvePoint,
    EcdhKeyPair,
    EcdsaSignature,
    EllipticCurve,
};
use crate::symbolic::finite_field::{
    PrimeField,
    PrimeFieldElement,
};
use num_bigint::BigInt;
use std::os::raw::c_char;
use std::str::FromStr;
use std::sync::Arc;

/// Use standard decimal string parsing for BigInts.

unsafe fn parse_bigint(s: *const c_char) -> Option<BigInt> {

    if let Some(str_slice) = c_str_to_str(s) {

        BigInt::from_str(str_slice).ok()
    } else {

        None
    }
}

/// Helper to convert BigInt to C string (decimal).

fn bigint_to_string(b: &BigInt) -> *mut c_char {

    to_c_string(b.to_string())
}

// --- EllipticCurve ---

/// Creates a new elliptic curve from decimal strings.
#[no_mangle]

pub unsafe extern "C" fn rssn_elliptic_curve_new(
    a_str: *const c_char,
    b_str: *const c_char,
    modulus_str: *const c_char,
) -> *mut EllipticCurve {

    let a = parse_bigint(a_str);

    let b = parse_bigint(b_str);

    let m = parse_bigint(modulus_str);

    if let (Some(a), Some(b), Some(m)) = (a, b, m) {

        let curve = EllipticCurve::new(a, b, m);

        Box::into_raw(Box::new(curve))
    } else {

        std::ptr::null_mut()
    }
}

/// Frees an elliptic curve handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_elliptic_curve_free(curve: *mut EllipticCurve) {

    if !curve.is_null() {

        drop(Box::from_raw(curve));
    }
}

// --- CurvePoint ---

/// Creates an affine curve point from decimal strings.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_point_affine(
    x_str: *const c_char,
    y_str: *const c_char,
    modulus_str: *const c_char,
) -> *mut CurvePoint {

    let x = parse_bigint(x_str);

    let y = parse_bigint(y_str);

    let m = parse_bigint(modulus_str);

    if let (Some(x), Some(y), Some(m)) = (x, y, m) {

        let field = PrimeField::new(m);

        let point = CurvePoint::Affine {
            x: PrimeFieldElement::new(x, field.clone()),
            y: PrimeFieldElement::new(y, field),
        };

        Box::into_raw(Box::new(point))
    } else {

        std::ptr::null_mut()
    }
}

/// Creates the point at infinity.
#[no_mangle]

pub extern "C" fn rssn_curve_point_infinity() -> *mut CurvePoint {

    Box::into_raw(Box::new(
        CurvePoint::Infinity,
    ))
}

/// Checks if a point is the point at infinity.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_point_is_infinity(point: *const CurvePoint) -> bool {

    if point.is_null() {

        return false;
    }

    (*point).is_infinity()
}

/// Gets the x-coordinate of an affine point as a string. Returns NULL if infinity.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_point_get_x(point: *const CurvePoint) -> *mut c_char {

    if point.is_null() {

        return std::ptr::null_mut();
    }

    if let Some(x) = (*point).x() {

        bigint_to_string(&x.value)
    } else {

        std::ptr::null_mut()
    }
}

/// Gets the y-coordinate of an affine point as a string. Returns NULL if infinity.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_point_get_y(point: *const CurvePoint) -> *mut c_char {

    if point.is_null() {

        return std::ptr::null_mut();
    }

    if let Some(y) = (*point).y() {

        bigint_to_string(&y.value)
    } else {

        std::ptr::null_mut()
    }
}

/// Frees a curve point handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_point_free(point: *mut CurvePoint) {

    if !point.is_null() {

        drop(Box::from_raw(point));
    }
}

// --- Curve Operations ---

#[no_mangle]

pub unsafe extern "C" fn rssn_curve_is_on_curve(
    curve: *const EllipticCurve,
    point: *const CurvePoint,
) -> bool {

    if curve.is_null() || point.is_null() {

        return false;
    }

    (*curve).is_on_curve(&*point)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_curve_negate(
    curve: *const EllipticCurve,
    point: *const CurvePoint,
) -> *mut CurvePoint {

    if curve.is_null() || point.is_null() {

        return std::ptr::null_mut();
    }

    let result = (*curve).negate(&*point);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_curve_double(
    curve: *const EllipticCurve,
    point: *const CurvePoint,
) -> *mut CurvePoint {

    if curve.is_null() || point.is_null() {

        return std::ptr::null_mut();
    }

    let result = (*curve).double(&*point);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_curve_add(
    curve: *const EllipticCurve,
    p1: *const CurvePoint,
    p2: *const CurvePoint,
) -> *mut CurvePoint {

    if curve.is_null() || p1.is_null() || p2.is_null() {

        return std::ptr::null_mut();
    }

    let result = (*curve).add(&*p1, &*p2);

    Box::into_raw(Box::new(result))
}

/// Scalar multiplication. k is a string.
#[no_mangle]

pub unsafe extern "C" fn rssn_curve_scalar_mult(
    curve: *const EllipticCurve,
    k_str: *const c_char,
    p: *const CurvePoint,
) -> *mut CurvePoint {

    let k = parse_bigint(k_str);

    if let (Some(curve), Some(k), Some(p)) = (
        curve.as_ref(),
        k,
        p.as_ref(),
    ) {

        let result = curve.scalar_mult(&k, p);

        Box::into_raw(Box::new(result))
    } else {

        std::ptr::null_mut()
    }
}

// --- Key Management ---

#[no_mangle]

pub unsafe extern "C" fn rssn_generate_keypair(
    curve: *const EllipticCurve,
    generator: *const CurvePoint,
) -> *mut EcdhKeyPair {

    if curve.is_null() || generator.is_null() {

        return std::ptr::null_mut();
    }

    let keypair = generate_keypair(&*curve, &*generator);

    Box::into_raw(Box::new(keypair))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_keypair_get_private_key(kp: *const EcdhKeyPair) -> *mut c_char {

    if let Some(k) = kp.as_ref() {

        bigint_to_string(&k.private_key)
    } else {

        std::ptr::null_mut()
    }
}

/// Returns a NEW handle to the public key point (must be freed).
#[no_mangle]

pub unsafe extern "C" fn rssn_keypair_get_public_key(kp: *const EcdhKeyPair) -> *mut CurvePoint {

    if let Some(k) = kp.as_ref() {

        Box::into_raw(Box::new(
            k.public_key.clone(),
        ))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_keypair_free(keypair: *mut EcdhKeyPair) {

    if !keypair.is_null() {

        drop(Box::from_raw(
            keypair,
        ));
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_generate_shared_secret(
    curve: *const EllipticCurve,
    private_key_str: *const c_char,
    other_public_key: *const CurvePoint,
) -> *mut CurvePoint {

    let pk = parse_bigint(private_key_str);

    if let (Some(c), Some(pk), Some(opub)) = (
        curve.as_ref(),
        pk,
        other_public_key.as_ref(),
    ) {

        let result = generate_shared_secret(c, &pk, opub);

        Box::into_raw(Box::new(result))
    } else {

        std::ptr::null_mut()
    }
}

// --- ECDSA ---

#[no_mangle]

pub unsafe extern "C" fn rssn_ecdsa_sign(
    message_hash_str: *const c_char,
    private_key_str: *const c_char,
    curve: *const EllipticCurve,
    generator: *const CurvePoint,
    order_str: *const c_char,
) -> *mut EcdsaSignature {

    let h = parse_bigint(message_hash_str);

    let pk = parse_bigint(private_key_str);

    let order = parse_bigint(order_str);

    if let (Some(h), Some(pk), Some(c), Some(g), Some(o)) = (
        h,
        pk,
        curve.as_ref(),
        generator.as_ref(),
        order,
    ) {

        match ecdsa_sign(&h, &pk, c, g, &o) {
            | Some(sig) => Box::into_raw(Box::new(sig)),
            | None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_ecdsa_verify(
    message_hash_str: *const c_char,
    signature: *const EcdsaSignature,
    public_key: *const CurvePoint,
    curve: *const EllipticCurve,
    generator: *const CurvePoint,
    order_str: *const c_char,
) -> bool {

    let h = parse_bigint(message_hash_str);

    let order = parse_bigint(order_str);

    if let (Some(h), Some(sig), Some(pk), Some(c), Some(g), Some(o)) = (
        h,
        signature.as_ref(),
        public_key.as_ref(),
        curve.as_ref(),
        generator.as_ref(),
        order,
    ) {

        ecdsa_verify(
            &h, sig, pk, c, g, &o,
        )
    } else {

        false
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_ecdsa_signature_get_r(sig: *const EcdsaSignature) -> *mut c_char {

    if let Some(s) = sig.as_ref() {

        bigint_to_string(&s.r)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_ecdsa_signature_get_s(sig: *const EcdsaSignature) -> *mut c_char {

    if let Some(s) = sig.as_ref() {

        bigint_to_string(&s.s)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_ecdsa_signature_free(sig: *mut EcdsaSignature) {

    if !sig.is_null() {

        drop(Box::from_raw(sig));
    }
}

// --- Compression ---

/// Compresses a point. Returns the x-coordinate string. sets *is_odd to the parity.
#[no_mangle]

pub unsafe extern "C" fn rssn_point_compress(
    point: *const CurvePoint,
    is_odd_out: *mut bool,
) -> *mut c_char {

    if let Some(p) = point.as_ref() {

        if let Some((x, is_odd)) = point_compress(p) {

            if !is_odd_out.is_null() {

                *is_odd_out = is_odd;
            }

            return bigint_to_string(&x);
        }
    }

    std::ptr::null_mut()
}

/// Decompresses a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_point_decompress(
    x_str: *const c_char,
    is_odd: bool,
    curve: *const EllipticCurve,
) -> *mut CurvePoint {

    let x = parse_bigint(x_str);

    if let (Some(x), Some(c)) = (x, curve.as_ref()) {

        if let Some(p) = point_decompress(x, is_odd, c) {

            return Box::into_raw(Box::new(p));
        }
    }

    std::ptr::null_mut()
}

///// Frees an elliptic curve handle.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_elliptic_curve_free(curve: *mut EllipticCurve) {
//     if !curve.is_null() {
//         drop(Box::from_raw(curve));
//     }
// }

// /// Creates an affine curve point.
// #[no_mangle]
// pub extern "C" fn rssn_curve_point_affine(x: i64, y: i64, modulus: i64) -> *mut CurvePoint {
//     let field = PrimeField::new(BigInt::from(modulus));
//     let point = CurvePoint::Affine {
//         x: PrimeFieldElement::new(BigInt::from(x), field.clone()),
//         y: PrimeFieldElement::new(BigInt::from(y), field),
//     };
//     Box::into_raw(Box::new(point))
// }

// /// Creates the point at infinity.
// #[no_mangle]
// pub extern "C" fn rssn_curve_point_infinity() -> *mut CurvePoint {
//     Box::into_raw(Box::new(CurvePoint::Infinity))
// }

// /// Checks if a point is the point at infinity.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_point_is_infinity(point: *const CurvePoint) -> bool {
//     if point.is_null() { return false; }
//     (*point).is_infinity()
// }

// /// Frees a curve point handle.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_point_free(point: *mut CurvePoint) {
//     if !point.is_null() {
//         drop(Box::from_raw(point));
//     }
// }

// /// Checks if a point is on the curve.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_is_on_curve(
//     curve: *const EllipticCurve,
//     point: *const CurvePoint
// ) -> bool {
//     if curve.is_null() || point.is_null() { return false; }
//     (*curve).is_on_curve(&*point)
// }

// /// Negates a point on the curve (P -> -P).
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_negate(
//     curve: *const EllipticCurve,
//     point: *const CurvePoint
// ) -> *mut CurvePoint {
//     if curve.is_null() || point.is_null() { return std::ptr::null_mut(); }
//     let result = (*curve).negate(&*point);
//     Box::into_raw(Box::new(result))
// }

// /// Doubles a point on the curve (2P).
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_double(
//     curve: *const EllipticCurve,
//     point: *const CurvePoint
// ) -> *mut CurvePoint {
//     if curve.is_null() || point.is_null() { return std::ptr::null_mut(); }
//     let result = (*curve).double(&*point);
//     Box::into_raw(Box::new(result))
// }

// /// Adds two curve points.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_add(
//     curve: *const EllipticCurve,
//     p1: *const CurvePoint,
//     p2: *const CurvePoint
// ) -> *mut CurvePoint {
//     if curve.is_null() || p1.is_null() || p2.is_null() {
//         return std::ptr::null_mut();
//     }
//     let result = (*curve).add(&*p1, &*p2);
//     Box::into_raw(Box::new(result))
// }

// /// Performs scalar multiplication k * P on elliptic curve.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_curve_scalar_mult(
//     curve: *const EllipticCurve,
//     k: i64,
//     p: *const CurvePoint
// ) -> *mut CurvePoint {
//     if curve.is_null() || p.is_null() {
//         return std::ptr::null_mut();
//     }
//     let result = (*curve).scalar_mult(&BigInt::from(k), &*p);
//     Box::into_raw(Box::new(result))
// }

// /// Generates an ECDH key pair.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_generate_keypair(
//     curve: *const EllipticCurve,
//     generator: *const CurvePoint
// ) -> *mut EcdhKeyPair {
//     if curve.is_null() || generator.is_null() {
//         return std::ptr::null_mut();
//     }
//     let keypair = generate_keypair(&*curve, &*generator);
//     Box::into_raw(Box::new(keypair))
// }

// /// Frees an ECDH key pair.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_keypair_free(keypair: *mut EcdhKeyPair) {
//     if !keypair.is_null() {
//         drop(Box::from_raw(keypair));
//     }
// }

// /// Generates a shared secret using ECDH.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_generate_shared_secret(
//     curve: *const EllipticCurve,
//     private_key: i64,
//     other_public_key: *const CurvePoint
// ) -> *mut CurvePoint {
//     if curve.is_null() || other_public_key.is_null() {
//         return std::ptr::null_mut();
//     }
//     let result = generate_shared_secret(&*curve, &BigInt::from(private_key), &*other_public_key);
//     Box::into_raw(Box::new(result))
// }

//// Signs a message using ECDSA.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_ecdsa_sign(
//     message_hash: i64,
//     private_key: i64,
//     curve: *const EllipticCurve,
//     generator: *const CurvePoint,
//     order: i64
// ) -> *mut EcdsaSignature {
//     if curve.is_null() || generator.is_null() {
//         return std::ptr::null_mut();
//     }
//     match ecdsa_sign(
//         &BigInt::from(message_hash),
//         &BigInt::from(private_key),
//         &*curve,
//         &*generator,
//         &BigInt::from(order)
//     ) {
//         Some(sig) => Box::into_raw(Box::new(sig)),
//         None => std::ptr::null_mut(),
//     }
// }

// Verifies an ECDSA signature.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_ecdsa_verify(
//     message_hash: i64,
//     signature: *const EcdsaSignature,
//     public_key: *const CurvePoint,
//     curve: *const EllipticCurve,
//     generator: *const CurvePoint,
//     order: i64
// ) -> bool {
//     if signature.is_null() || public_key.is_null() || curve.is_null() || generator.is_null() {
//         return false;
//     }
//     ecdsa_verify(
//         &BigInt::from(message_hash),
//         &*signature,
//         &*public_key,
//         &*curve,
//         &*generator,
//         &BigInt::from(order)
//     )
// }

// /// Frees an ECDSA signature.
// #[no_mangle]
// pub unsafe extern "C" fn rssn_ecdsa_signature_free(sig: *mut EcdsaSignature) {
//     if !sig.is_null() {
//         drop(Box::from_raw(sig));
//     }
// }
