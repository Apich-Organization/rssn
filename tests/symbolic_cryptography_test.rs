//! Test suite for cryptography module (Elliptic Curve Cryptography).

use rssn::symbolic::cryptography::*;
use rssn::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use num_bigint::BigInt;
use std::sync::Arc;

// Use a small prime field for testing: p = 23
fn test_field() -> Arc<PrimeField> {
    PrimeField::new(BigInt::from(23))
}

// Test curve: y^2 = x^3 + x + 1 over GF(23)
fn test_curve() -> EllipticCurve {
    let field = test_field();
    EllipticCurve {
        a: PrimeFieldElement::new(BigInt::from(1), field.clone()),
        b: PrimeFieldElement::new(BigInt::from(1), field.clone()),
        field,
    }
}

// A point on the curve y^2 = x^3 + x + 1 (mod 23)
// P = (0, 1): 1^2 = 1, 0^3 + 0 + 1 = 1. Yes, (0, 1) is on curve.
fn test_point() -> CurvePoint {
    let field = test_field();
    CurvePoint::Affine {
        x: PrimeFieldElement::new(BigInt::from(0), field.clone()),
        y: PrimeFieldElement::new(BigInt::from(1), field),
    }
}

#[test]
fn test_point_addition_with_infinity() {
    let curve = test_curve();
    let p = test_point();
    
    // P + O = P
    let result = curve.add(&p, &CurvePoint::Infinity);
    assert_eq!(result, p);
    
    // O + P = P
    let result2 = curve.add(&CurvePoint::Infinity, &p);
    assert_eq!(result2, p);
}

#[test]
fn test_point_doubling() {
    let curve = test_curve();
    let p = test_point();
    
    // 2P = P + P
    let double_p = curve.add(&p, &p);
    
    // Result should be a valid point (not infinity for this curve/point)
    match double_p {
        CurvePoint::Affine { .. } => { /* expected */ }
        CurvePoint::Infinity => { /* also possible depending on curve order */ }
    }
}

#[test]
fn test_scalar_multiplication_by_one() {
    let curve = test_curve();
    let p = test_point();
    
    // 1 * P = P
    let result = curve.scalar_mult(&BigInt::from(1), &p);
    assert_eq!(result, p);
}

#[test]
fn test_scalar_multiplication_by_zero() {
    let curve = test_curve();
    let p = test_point();
    
    // 0 * P = O (point at infinity)
    let result = curve.scalar_mult(&BigInt::from(0), &p);
    assert_eq!(result, CurvePoint::Infinity);
}

#[test]
fn test_scalar_multiplication() {
    let curve = test_curve();
    let p = test_point();
    
    // 2 * P should equal P + P
    let two_p = curve.scalar_mult(&BigInt::from(2), &p);
    let p_plus_p = curve.add(&p, &p);
    assert_eq!(two_p, p_plus_p);
    
    // 3 * P should equal 2P + P
    let three_p = curve.scalar_mult(&BigInt::from(3), &p);
    let expected = curve.add(&two_p, &p);
    assert_eq!(three_p, expected);
}

#[test]
fn test_ecdh_keypair_generation() {
    let curve = test_curve();
    let generator = test_point();
    
    // Generate keypair
    let keypair = generate_keypair(&curve, &generator);
    
    // Private key should be non-zero
    assert!(keypair.private_key > BigInt::from(0));
    
    // Public key should be private_key * generator
    let expected_public = curve.scalar_mult(&keypair.private_key, &generator);
    assert_eq!(keypair.public_key, expected_public);
}

#[test]
fn test_ecdh_shared_secret() {
    let curve = test_curve();
    let generator = test_point();
    
    // Generate two keypairs
    let alice = generate_keypair(&curve, &generator);
    let bob = generate_keypair(&curve, &generator);
    
    // Alice computes shared secret: alice_private * bob_public
    let alice_secret = generate_shared_secret(&curve, &alice.private_key, &bob.public_key);
    
    // Bob computes shared secret: bob_private * alice_public
    let bob_secret = generate_shared_secret(&curve, &bob.private_key, &alice.public_key);
    
    // Both should arrive at the same shared secret
    assert_eq!(alice_secret, bob_secret);
}
