//! # Cryptography Module
//!
//! This module provides implementations for cryptographic primitives and algorithms,
//! particularly focusing on elliptic curve cryptography (ECC). It includes structures
//! for elliptic curves over finite fields, curve points, and functions for key generation
//! and shared secret derivation using ECDH (Elliptic Curve Diffie-Hellman).

use crate::symbolic::finite_field::{PrimeField, PrimeFieldElement};
use num_bigint::{BigInt, RandBigInt};
use num_traits::{One, Zero};
use rand::Rng;

use std::sync::Arc;

/// Represents an elliptic curve over a prime field: y^2 = x^3 + ax + b.
#[derive(Clone)]
pub struct EllipticCurve {
    pub a: PrimeFieldElement,
    pub b: PrimeFieldElement,
    pub field: Arc<PrimeField>,
}

/// Represents a point on an elliptic curve, including the point at infinity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CurvePoint {
    Infinity,
    Affine {
        x: PrimeFieldElement,
        y: PrimeFieldElement,
    },
}

#[derive(Debug, Clone)]
pub struct EcdhKeyPair {
    pub private_key: BigInt,
    pub public_key: CurvePoint,
}

impl EllipticCurve {
    /// Adds two points on the curve.
    ///
    /// This function implements the elliptic curve point addition rules.
    /// It handles cases for adding a point to the point at infinity, adding
    /// two distinct points, and doubling a point.
    ///
    /// # Arguments
    /// * `p1` - The first `CurvePoint`.
    /// * `p2` - The second `CurvePoint`.
    ///
    /// # Returns
    /// A new `CurvePoint` representing the sum of `p1` and `p2`.
    pub fn add(&self, p1: &CurvePoint, p2: &CurvePoint) -> CurvePoint {
        match (p1, p2) {
            (CurvePoint::Infinity, p) => p.clone(),
            (p, CurvePoint::Infinity) => p.clone(),
            (CurvePoint::Affine { x: x1, y: y1 }, CurvePoint::Affine { x: x2, y: y2 }) => {
                if x1 == x2 && *y1 != *y2 {
                    // Points are inverses of each other
                    return CurvePoint::Infinity;
                }

                let m = if x1 == x2 && y1 == y2 {
                    // Point doubling
                    let three = PrimeFieldElement::new(BigInt::from(3), self.field.clone());
                    let two = PrimeFieldElement::new(BigInt::from(2), self.field.clone());
                    let num = three * x1.clone() * x1.clone() + self.a.clone();
                    let den = two * y1.clone();
                    num / den
                } else {
                    // Point addition
                    (y2.clone() - y1.clone()) / (x2.clone() - x1.clone())
                };

                let x3 = m.clone() * m.clone() - x1.clone() - x2.clone();
                let y3 = m * (x1.clone() - x3.clone()) - y1.clone();
                CurvePoint::Affine { x: x3, y: y3 }
            }
        }
    }

    /// Performs scalar multiplication (`k * P`) using the double-and-add algorithm.
    ///
    /// This algorithm efficiently computes `k` times a point `P` on the elliptic curve
    /// by repeatedly doubling `P` and adding `P` based on the binary representation of `k`.
    ///
    /// # Arguments
    /// * `k` - The scalar `BigInt`.
    /// * `p` - The `CurvePoint` to multiply.
    ///
    /// # Returns
    /// A new `CurvePoint` representing `k * P`.
    pub fn scalar_mult(&self, k: &BigInt, p: &CurvePoint) -> CurvePoint {
        let mut res = CurvePoint::Infinity;
        let mut app = p.clone();
        let mut k_clone = k.clone();

        while k_clone > Zero::zero() {
            if &k_clone % 2 != Zero::zero() {
                res = self.add(&res, &app);
            }
            app = self.add(&app, &app);
            k_clone >>= 1;
        }
        res
    }
}

/// Generates a new ECDH (Elliptic Curve Diffie-Hellman) key pair.
///
/// This function randomly selects a private key (a scalar) and computes the
/// corresponding public key (a point on the elliptic curve) by scalar multiplication
/// of the curve's generator point.
///
/// # Arguments
/// * `curve` - The `EllipticCurve` parameters.
/// * `generator` - The base `CurvePoint` (generator point) of the curve.
///
/// # Returns
/// An `EcdhKeyPair` containing the generated private and public keys.
pub fn generate_keypair(curve: &EllipticCurve, generator: &CurvePoint) -> EcdhKeyPair {
    let mut rng = rand::thread_rng();
    // In a real scenario, the private key should be chosen from a specific subgroup order.
    let private_key = rng.gen_bigint_range(&BigInt::one(), &curve.field.modulus);
    let public_key = curve.scalar_mult(&private_key, generator);
    EcdhKeyPair {
        private_key,
        public_key,
    }
}

/// Generates a shared secret using one's own private key and the other party's public key.
///
/// In ECDH, the shared secret is derived by performing scalar multiplication of the
/// other party's public key with one's own private key. This results in a common
/// `CurvePoint` that only both parties can compute.
///
/// # Arguments
/// * `curve` - The `EllipticCurve` parameters.
/// * `own_private_key` - Your own private key (`BigInt`).
/// * `other_public_key` - The other party's public key (`CurvePoint`).
///
/// # Returns
/// A `CurvePoint` representing the shared secret.
pub fn generate_shared_secret(
    curve: &EllipticCurve,
    own_private_key: &BigInt,
    other_public_key: &CurvePoint,
) -> CurvePoint {
    curve.scalar_mult(own_private_key, other_public_key)
}
