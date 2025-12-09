//! Bincode-based FFI API for cryptographic operations.
//!
//! This module provides binary serialization-based FFI functions for elliptic curve cryptography.

use crate::symbolic::cryptography::{
    EllipticCurve, CurvePoint, EcdhKeyPair, EcdsaSignature,
    generate_keypair, generate_shared_secret, point_compress, point_decompress,
    ecdsa_sign, ecdsa_verify,
};
use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use crate::ffi_apis::common::*;
use num_bigint::BigInt;
use std::sync::Arc;

/// Creates a new elliptic curve via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_elliptic_curve_new(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer
) -> BincodeBuffer {
    let a: Option<BigInt> = from_bincode_buffer(&a_buf);
    let b: Option<BigInt> = from_bincode_buffer(&b_buf);
    let modulus: Option<BigInt> = from_bincode_buffer(&modulus_buf);

    if let (Some(a), Some(b), Some(m)) = (a, b, modulus) {
        to_bincode_buffer(&EllipticCurve::new(a, b, m))
    } else {
        BincodeBuffer::empty()
    }
}

/// Creates an affine curve point via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_point_affine(
    x_buf: BincodeBuffer,
    y_buf: BincodeBuffer,
    modulus_buf: BincodeBuffer
) -> BincodeBuffer {
    let x: Option<BigInt> = from_bincode_buffer(&x_buf);
    let y: Option<BigInt> = from_bincode_buffer(&y_buf);
    let modulus: Option<BigInt> = from_bincode_buffer(&modulus_buf);

    if let (Some(x), Some(y), Some(m)) = (x, y, modulus) {
        let field = PrimeField::new(m);
        let point = CurvePoint::Affine {
            x: PrimeFieldElement::new(x, field.clone()),
            y: PrimeFieldElement::new(y, field)
        };
        to_bincode_buffer(&point)
    } else {
        BincodeBuffer::empty()
    }
}

/// Creates a point at infinity via Bincode interface.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_point_infinity() -> BincodeBuffer {
    to_bincode_buffer(&CurvePoint::Infinity)
}

/// Checks if a point is on the curve.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_is_on_curve(
    curve_buf: BincodeBuffer,
    point_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let point: Option<CurvePoint> = from_bincode_buffer(&point_buf);
    
    if let (Some(c), Some(p)) = (curve, point) {
        to_bincode_buffer(&c.is_on_curve(&p))
    } else {
        BincodeBuffer::empty()
    }
}

/// Negates a point.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_negate(
    curve_buf: BincodeBuffer,
    point_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let point: Option<CurvePoint> = from_bincode_buffer(&point_buf);
    
    if let (Some(c), Some(p)) = (curve, point) {
        to_bincode_buffer(&c.negate(&p))
    } else {
        BincodeBuffer::empty()
    }
}

/// Doubles a point.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_double(
    curve_buf: BincodeBuffer,
    point_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let point: Option<CurvePoint> = from_bincode_buffer(&point_buf);
    
    if let (Some(c), Some(p)) = (curve, point) {
        to_bincode_buffer(&c.double(&p))
    } else {
        BincodeBuffer::empty()
    }
}

/// Adds two points.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_add(
    curve_buf: BincodeBuffer,
    p1_buf: BincodeBuffer,
    p2_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let p1: Option<CurvePoint> = from_bincode_buffer(&p1_buf);
    let p2: Option<CurvePoint> = from_bincode_buffer(&p2_buf);
    
    if let (Some(c), Some(p1), Some(p2)) = (curve, p1, p2) {
        to_bincode_buffer(&c.add(&p1, &p2))
    } else {
        BincodeBuffer::empty()
    }
}

/// Scalar multiplication.
#[no_mangle]
pub extern "C" fn rssn_bincode_curve_scalar_mult(
    curve_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    p_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let k: Option<BigInt> = from_bincode_buffer(&k_buf);
    let p: Option<CurvePoint> = from_bincode_buffer(&p_buf);
    
    if let (Some(c), Some(k), Some(p)) = (curve, k, p) {
        to_bincode_buffer(&c.scalar_mult(&k, &p))
    } else {
        BincodeBuffer::empty()
    }
}

/// Generates a key pair.
#[no_mangle]
pub extern "C" fn rssn_bincode_generate_keypair(
    curve_buf: BincodeBuffer,
    generator_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let gen: Option<CurvePoint> = from_bincode_buffer(&generator_buf);
    
    if let (Some(c), Some(g)) = (curve, gen) {
        to_bincode_buffer(&generate_keypair(&c, &g))
    } else {
        BincodeBuffer::empty()
    }
}

/// Generates a shared secret.
#[no_mangle]
pub extern "C" fn rssn_bincode_generate_shared_secret(
    curve_buf: BincodeBuffer,
    private_key_buf: BincodeBuffer,
    other_public_key_buf: BincodeBuffer
) -> BincodeBuffer {
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let pk: Option<BigInt> = from_bincode_buffer(&private_key_buf);
    let other_pub: Option<CurvePoint> = from_bincode_buffer(&other_public_key_buf);
    
    if let (Some(c), Some(pk), Some(opub)) = (curve, pk, other_pub) {
        to_bincode_buffer(&generate_shared_secret(&c, &pk, &opub))
    } else {
        BincodeBuffer::empty()
    }
}

/// Signs a message.
#[no_mangle]
pub extern "C" fn rssn_bincode_ecdsa_sign(
    message_hash_buf: BincodeBuffer,
    private_key_buf: BincodeBuffer,
    curve_buf: BincodeBuffer,
    generator_buf: BincodeBuffer,
    order_buf: BincodeBuffer
) -> BincodeBuffer {
    let hash: Option<BigInt> = from_bincode_buffer(&message_hash_buf);
    let pk: Option<BigInt> = from_bincode_buffer(&private_key_buf);
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let gen: Option<CurvePoint> = from_bincode_buffer(&generator_buf);
    let order: Option<BigInt> = from_bincode_buffer(&order_buf);
    
    if let (Some(h), Some(pk), Some(c), Some(g), Some(o)) = 
        (hash, pk, curve, gen, order) {
        if let Some(sig) = ecdsa_sign(&h, &pk, &c, &g, &o) {
            to_bincode_buffer(&sig)
        } else {
            BincodeBuffer::empty()
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Verifies a signature.
#[no_mangle]
pub extern "C" fn rssn_bincode_ecdsa_verify(
    message_hash_buf: BincodeBuffer,
    signature_buf: BincodeBuffer,
    public_key_buf: BincodeBuffer,
    curve_buf: BincodeBuffer,
    generator_buf: BincodeBuffer,
    order_buf: BincodeBuffer
) -> BincodeBuffer {
    let hash: Option<BigInt> = from_bincode_buffer(&message_hash_buf);
    let sig: Option<EcdsaSignature> = from_bincode_buffer(&signature_buf);
    let pub_key: Option<CurvePoint> = from_bincode_buffer(&public_key_buf);
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    let gen: Option<CurvePoint> = from_bincode_buffer(&generator_buf);
    let order: Option<BigInt> = from_bincode_buffer(&order_buf);
    
    if let (Some(h), Some(sig), Some(pk), Some(c), Some(g), Some(o)) = 
        (hash, sig, pub_key, curve, gen, order) {
        to_bincode_buffer(&ecdsa_verify(&h, &sig, &pk, &c, &g, &o))
    } else {
        BincodeBuffer::empty()
    }
}

/// Compresses a point.
#[no_mangle]
pub extern "C" fn rssn_bincode_point_compress(
    point_buf: BincodeBuffer
) -> BincodeBuffer {
    let point: Option<CurvePoint> = from_bincode_buffer(&point_buf);
    if let Some(p) = point {
        if let Some(res) = point_compress(&p) {
            to_bincode_buffer(&res)
        } else {
            BincodeBuffer::empty()
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Decompresses a point.
#[no_mangle]
pub extern "C" fn rssn_bincode_point_decompress(
    x_buf: BincodeBuffer,
    is_odd_buf: BincodeBuffer,
    curve_buf: BincodeBuffer
) -> BincodeBuffer {
    let x: Option<BigInt> = from_bincode_buffer(&x_buf);
    let is_odd: Option<bool> = from_bincode_buffer(&is_odd_buf);
    let curve: Option<EllipticCurve> = from_bincode_buffer(&curve_buf);
    
    if let (Some(x), Some(io), Some(c)) = (x, is_odd, curve) {
        if let Some(p) = point_decompress(x, io, &c) {
            to_bincode_buffer(&p)
        } else {
            BincodeBuffer::empty()
        }
    } else {
        BincodeBuffer::empty()
    }
}
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
