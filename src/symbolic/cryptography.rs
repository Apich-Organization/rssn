//! # Cryptography Module
//!
//! This module provides implementations for cryptographic primitives and algorithms,
//! particularly focusing on elliptic curve cryptography (ECC). It includes structures
//! for elliptic curves over finite fields, curve points, and functions for key generation
//! and shared secret derivation using ECDH (Elliptic Curve Diffie-Hellman).
//!
//! ## Features
//! - Elliptic curve operations (point addition, scalar multiplication)
//! - ECDH key exchange
//! - ECDSA digital signatures
//! - Point compression/decompression

use std::sync::Arc;

use num_bigint::BigInt;
use num_bigint::RandBigInt;
use num_traits::One;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::finite_field::PrimeField;
use crate::symbolic::finite_field::PrimeFieldElement;

/// Represents an elliptic curve over a prime field: y^2 = x^3 + ax + b.
#[derive(
    Clone, Serialize, Deserialize,
)]

pub struct EllipticCurve {
    /// Coefficient 'a' in the curve equation
    pub a: PrimeFieldElement,
    /// Coefficient 'b' in the curve equation
    pub b: PrimeFieldElement,
    /// The underlying prime field
    pub field: Arc<PrimeField>,
}

/// Represents a point on an elliptic curve, including the point at infinity.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub enum CurvePoint {
    /// The point at infinity (identity element)
    Infinity,
    /// An affine point with coordinates (x, y)
    Affine {
        /// The x-coordinate of the point.
        x: PrimeFieldElement,
        /// The y-coordinate of the point.
        y: PrimeFieldElement,
    },
}

/// ECDH key pair containing private and public keys.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub struct EcdhKeyPair {
    /// The private key (a scalar)
    pub private_key: BigInt,
    /// The public key (a curve point)
    pub public_key: CurvePoint,
}

/// ECDSA signature containing r and s components.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub struct EcdsaSignature {
    /// The r component of the signature
    pub r: BigInt,
    /// The s component of the signature
    pub s: BigInt,
}

impl CurvePoint {
    /// Returns true if this point is the point at infinity.
    #[must_use]

    pub const fn is_infinity(
        &self
    ) -> bool {

        matches!(self, Self::Infinity)
    }

    /// Returns the x-coordinate if this is an affine point, None if infinity.
    #[must_use]

    pub const fn x(
        &self
    ) -> Option<&PrimeFieldElement>
    {

        match self {
            | Self::Affine {
                x,
                ..
            } => Some(x),
            | Self::Infinity => None,
        }
    }

    /// Returns the y-coordinate if this is an affine point, None if infinity.
    #[must_use]

    pub const fn y(
        &self
    ) -> Option<&PrimeFieldElement>
    {

        match self {
            | Self::Affine {
                y,
                ..
            } => Some(y),
            | Self::Infinity => None,
        }
    }
}

impl EllipticCurve {
    /// Creates a new elliptic curve y^2 = x^3 + ax + b over GF(p).
    ///
    /// # Arguments
    /// * `a` - Coefficient 'a' as `BigInt`
    /// * `b` - Coefficient 'b' as `BigInt`
    /// * `modulus` - The prime modulus defining the field
    ///
    /// # Returns
    /// A new `EllipticCurve` instance
    #[must_use]

    pub fn new(
        a: BigInt,
        b: BigInt,
        modulus: BigInt,
    ) -> Self {

        let field =
            PrimeField::new(modulus);

        Self {
            a: PrimeFieldElement::new(
                a,
                field.clone(),
            ),
            b: PrimeFieldElement::new(
                b,
                field.clone(),
            ),
            field,
        }
    }

    /// Checks if a point is on the curve.
    ///
    /// Verifies that y^2 = x^3 + ax + b for the given point.
    ///
    /// # Arguments
    /// * `point` - The point to verify
    ///
    /// # Returns
    /// `true` if the point is on the curve, `false` otherwise
    #[must_use]

    pub fn is_on_curve(
        &self,
        point: &CurvePoint,
    ) -> bool {

        match point {
            | CurvePoint::Infinity => {
                true
            },
            | CurvePoint::Affine {
                x,
                y,
            } => {

                let lhs = y.clone()
                    * y.clone();

                let rhs = x.clone()
                    * x.clone()
                    * x.clone()
                    + self.a.clone()
                        * x.clone()
                    + self.b.clone();

                lhs == rhs
            },
        }
    }

    /// Negates a point on the curve (P -> -P).
    ///
    /// For an affine point (x, y), the negation is (x, -y).
    ///
    /// # Arguments
    /// * `point` - The point to negate
    ///
    /// # Returns
    /// The negated point
    #[must_use]

    pub fn negate(
        &self,
        point: &CurvePoint,
    ) -> CurvePoint {

        match point {
            | CurvePoint::Infinity => {
                CurvePoint::Infinity
            },
            | CurvePoint::Affine {
                x,
                y,
            } => {
                CurvePoint::Affine {
                    x: x.clone(),
                    y: -y.clone(),
                }
            },
        }
    }

    /// Doubles a point on the curve (2P).
    ///
    /// This is an optimized version of add(P, P).
    ///
    /// # Arguments
    /// * `point` - The point to double
    ///
    /// # Returns
    /// The doubled point
    #[must_use]

    pub fn double(
        &self,
        point: &CurvePoint,
    ) -> CurvePoint {

        match point {
            | CurvePoint::Infinity => {
                CurvePoint::Infinity
            },
            | CurvePoint::Affine {
                x,
                y,
            } => {

                if y.value.is_zero() {

                    return CurvePoint::Infinity;
                }

                let three = PrimeFieldElement::new(
                    BigInt::from(3),
                    self.field.clone(),
                );

                let two = PrimeFieldElement::new(
                    BigInt::from(2),
                    self.field.clone(),
                );

                let m = (three
                    * x.clone()
                    * x.clone()
                    + self.a.clone())
                    / (two * y.clone());

                let x3 = m.clone()
                    * m.clone()
                    - x.clone()
                    - x.clone();

                let y3 = m
                    * (x.clone()
                        - x3.clone())
                    - y.clone();

                CurvePoint::Affine {
                    x: x3,
                    y: y3,
                }
            },
        }
    }

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
    #[must_use]

    pub fn add(
        &self,
        p1: &CurvePoint,
        p2: &CurvePoint,
    ) -> CurvePoint {

        match (p1, p2) {
            | (CurvePoint::Infinity, p) | (p, CurvePoint::Infinity) => p.clone(),
            | (
                CurvePoint::Affine {
                    x: x1,
                    y: y1,
                },
                CurvePoint::Affine {
                    x: x2,
                    y: y2,
                },
            ) => {

                if x1 == x2
                    && *y1 != *y2
                {

                    return CurvePoint::Infinity;
                }

                if x1 == x2 && y1 == y2
                {

                    return self
                        .double(p1);
                }

                let m = (y2.clone()
                    - y1.clone())
                    / (x2.clone()
                        - x1.clone());

                let x3 = m.clone()
                    * m.clone()
                    - x1.clone()
                    - x2.clone();

                let y3 = m
                    * (x1.clone()
                        - x3.clone())
                    - y1.clone();

                CurvePoint::Affine {
                    x: x3,
                    y: y3,
                }
            },
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
    #[must_use]

    pub fn scalar_mult(
        &self,
        k: &BigInt,
        p: &CurvePoint,
    ) -> CurvePoint {

        let mut res =
            CurvePoint::Infinity;

        let mut app = p.clone();

        let mut k_clone = k.clone();

        while k_clone > Zero::zero() {

            if &k_clone % 2
                != Zero::zero()
            {

                res = self
                    .add(&res, &app);
            }

            app = self.double(&app);

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
#[must_use]

pub fn generate_keypair(
    curve: &EllipticCurve,
    generator: &CurvePoint,
) -> EcdhKeyPair {

    let mut rng = rand::thread_rng();

    let private_key = rng
        .gen_bigint_range(
            &BigInt::one(),
            &curve.field.modulus,
        );

    let public_key = curve.scalar_mult(
        &private_key,
        generator,
    );

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
#[must_use]

pub fn generate_shared_secret(
    curve: &EllipticCurve,
    own_private_key: &BigInt,
    other_public_key: &CurvePoint,
) -> CurvePoint {

    curve.scalar_mult(
        own_private_key,
        other_public_key,
    )
}

/// Compresses a curve point to its x-coordinate and a sign bit.
///
/// Point compression reduces storage by only keeping the x-coordinate
/// and the parity of the y-coordinate.
///
/// # Arguments
/// * `point` - The point to compress
///
/// # Returns
/// `Some((x, is_y_odd))` for affine points, `None` for infinity
#[must_use]

pub fn point_compress(
    point: &CurvePoint
) -> Option<(BigInt, bool)> {

    match point {
        | CurvePoint::Infinity => None,
        | CurvePoint::Affine {
            x,
            y,
        } => {

            let is_y_odd = &y.value % 2
                != BigInt::zero();

            Some((
                x.value.clone(),
                is_y_odd,
            ))
        },
    }
}

/// Decompresses a point from its x-coordinate and sign bit.
///
/// Reconstructs the full point by computing y from the curve equation
/// and selecting the correct root based on the sign bit.
///
/// # Arguments
/// * `x` - The x-coordinate
/// * `is_y_odd` - Whether the y-coordinate is odd
/// * `curve` - The elliptic curve
///
/// # Returns
/// `Some(CurvePoint)` if successful, `None` if x is not on the curve
#[must_use]

pub fn point_decompress(
    x: BigInt,
    is_y_odd: bool,
    curve: &EllipticCurve,
) -> Option<CurvePoint> {

    let x_elem = PrimeFieldElement::new(
        x,
        curve.field.clone(),
    );

    // Compute y^2 = x^3 + ax + b
    let y_squared = x_elem.clone()
        * x_elem.clone()
        * x_elem.clone()
        + curve.a.clone()
            * x_elem.clone()
        + curve.b.clone();

    // Compute modular square root using Tonelli-Shanks would be ideal,
    // but for simplicity we use a brute-force approach for small fields
    // In production, use proper modular sqrt
    let modulus = &curve.field.modulus;

    let y_squared_val =
        &y_squared.value;

    // For p ≡ 3 (mod 4), y = y_squared^((p+1)/4)
    let exp = (modulus + 1) / 4;

    let y_val = y_squared_val
        .modpow(&exp, modulus);

    // Verify
    if (&y_val * &y_val) % modulus
        != y_squared_val % modulus
    {

        return None; // Not a quadratic residue
    }

    let y_is_odd =
        &y_val % 2 != BigInt::zero();

    let y_final =
        if y_is_odd == is_y_odd {

            y_val
        } else {

            modulus - &y_val
        };

    Some(CurvePoint::Affine {
        x: x_elem,
        y: PrimeFieldElement::new(
            y_final,
            curve.field.clone(),
        ),
    })
}

/// Signs a message hash using ECDSA.
///
/// # Arguments
/// * `message_hash` - The hash of the message as `BigInt`
/// * `private_key` - The signer's private key
/// * `curve` - The elliptic curve parameters
/// * `generator` - The curve's generator point
/// * `order` - The order of the generator point
///
/// # Returns
/// An `EcdsaSignature` containing (r, s) components
#[must_use]

pub fn ecdsa_sign(
    message_hash: &BigInt,
    private_key: &BigInt,
    curve: &EllipticCurve,
    generator: &CurvePoint,
    order: &BigInt,
) -> Option<EcdsaSignature> {

    let mut rng = rand::thread_rng();

    // Generate random k
    let k = rng.gen_bigint_range(
        &BigInt::one(),
        order,
    );

    // R = k * G
    let r_point = curve
        .scalar_mult(&k, generator);

    let r = match &r_point {
        | CurvePoint::Affine {
            x,
            ..
        } => x.value.clone() % order,
        | CurvePoint::Infinity => {
            return None
        },
    };

    if r.is_zero() {

        return None;
    }

    // s = k^(-1) * (hash + r * private_key) mod order
    let k_inv = mod_inverse(&k, order)?;

    let s = (&k_inv
        * (message_hash
            + &r * private_key))
        % order;

    if s.is_zero() {

        return None;
    }

    Some(EcdsaSignature {
        r,
        s,
    })
}

/// Verifies an ECDSA signature.
///
/// # Arguments
/// * `message_hash` - The hash of the message as `BigInt`
/// * `signature` - The signature to verify
/// * `public_key` - The signer's public key
/// * `curve` - The elliptic curve parameters
/// * `generator` - The curve's generator point
/// * `order` - The order of the generator point
///
/// # Returns
/// `true` if the signature is valid, `false` otherwise
#[must_use]

pub fn ecdsa_verify(
    message_hash: &BigInt,
    signature: &EcdsaSignature,
    public_key: &CurvePoint,
    curve: &EllipticCurve,
    generator: &CurvePoint,
    order: &BigInt,
) -> bool {

    // Check r and s are in valid range
    if signature.r <= BigInt::zero()
        || signature.r >= *order
    {

        return false;
    }

    if signature.s <= BigInt::zero()
        || signature.s >= *order
    {

        return false;
    }

    // w = s^(-1) mod order
    let w = match mod_inverse(
        &signature.s,
        order,
    ) {
        | Some(w) => w,
        | None => return false,
    };

    // u1 = (hash * w) mod order
    let u1 =
        (message_hash * &w) % order;

    // u2 = (r * w) mod order
    let u2 =
        (&signature.r * &w) % order;

    // R' = u1*G + u2*Q
    let point1 = curve
        .scalar_mult(&u1, generator);

    let point2 = curve
        .scalar_mult(&u2, public_key);

    let r_prime =
        curve.add(&point1, &point2);

    match r_prime {
        | CurvePoint::Infinity => false,
        | CurvePoint::Affine {
            x,
            ..
        } => {

            let v = x.value % order;

            v == signature.r
        },
    }
}

/// Computes modular inverse using extended Euclidean algorithm.

fn mod_inverse(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigInt> {

    let (g, x, _) = extended_gcd(a, m);

    if g != BigInt::one() {

        return None;
    }

    Some(((x % m) + m) % m)
}

/// Extended Euclidean algorithm.

fn extended_gcd(
    a: &BigInt,
    b: &BigInt,
) -> (
    BigInt,
    BigInt,
    BigInt,
) {

    if b.is_zero() {

        (
            a.clone(),
            BigInt::one(),
            BigInt::zero(),
        )
    } else {

        let (g, x, y) = extended_gcd(
            b,
            &(a % b),
        );

        (
            g,
            y.clone(),
            x - (a / b) * y,
        )
    }
}
