//! # Geometric Algebra
//!
//! This module provides tools for computations in Clifford and Geometric Algebra.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use num_bigint::BigInt;
use num_traits::One;
use std::collections::BTreeMap;
use std::ops::{Add, Mul, Sub};

/// Represents a multivector in a Clifford algebra.
/// The basis blades are represented by a bitmask. E.g., in 3D:
/// 001 (1) -> e1, 010 (2) -> e2, 100 (4) -> e3
/// 011 (3) -> e12, 101 (5) -> e13, 110 (6) -> e23
/// 111 (7) -> e123 (pseudoscalar)
#[derive(Clone, Debug, PartialEq)]
pub struct Multivector {
    /// A map from the basis blade bitmask to its coefficient.
    pub terms: BTreeMap<u32, Expr>,
    /// The signature of the algebra, e.g., (p, q, r) for (e_i^2 = +1, e_j^2 = -1, e_k^2 = 0)
    pub signature: (u32, u32, u32),
}

impl Multivector {
    /// Creates a new, empty multivector for a given algebra signature.
    ///
    /// # Arguments
    /// * `signature` - A tuple `(p, q, r)` defining the metric of the algebra, where:
    ///   - `p` is the number of basis vectors that square to +1.
    ///   - `q` is the number of basis vectors that square to -1.
    ///   - `r` is the number of basis vectors that square to 0.
    pub fn new(signature: (u32, u32, u32)) -> Self {
        Multivector {
            terms: BTreeMap::new(),
            signature,
        }
    }

    /// Creates a new multivector representing a scalar value.
    ///
    /// A scalar is a grade-0 element of the algebra.
    ///
    /// # Arguments
    /// * `signature` - The signature of the algebra `(p, q, r)`.
    /// * `value` - The scalar value as an `Expr`.
    ///
    /// # Returns
    /// A `Multivector` with a single term for the scalar part (grade 0).
    pub fn scalar(signature: (u32, u32, u32), value: Expr) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(0, value);
        Multivector { terms, signature }
    }

    /// Computes the geometric product of this multivector with another.
    ///
    /// The geometric product is the fundamental product of geometric algebra, combining
    /// the properties of the inner and outer products. It is associative and distributive
    /// but not generally commutative.
    ///
    /// The product of two basis blades `e_A` and `e_B` is computed by considering
    /// commutation rules (swaps) and contractions based on the algebra's metric signature.
    ///
    /// # Arguments
    /// * `other` - The `Multivector` to multiply with.
    ///
    /// # Returns
    /// A new `Multivector` representing the geometric product.
    pub fn geometric_product(&self, other: &Multivector) -> Multivector {
        let mut result = Multivector::new(self.signature);
        for (blade1, coeff1) in &self.terms {
            for (blade2, coeff2) in &other.terms {
                let (sign, metric_scalar, result_blade) = self.blade_product(*blade1, *blade2);

                let new_coeff = simplify(Expr::Mul(
                    Box::new(coeff1.clone()),
                    Box::new(coeff2.clone()),
                ));
                let signed_coeff = simplify(Expr::Mul(
                    Box::new(Expr::Constant(sign)),
                    Box::new(new_coeff),
                ));
                let final_coeff =
                    simplify(Expr::Mul(Box::new(signed_coeff), Box::new(metric_scalar)));

                if let Some(existing_coeff) = result.terms.get_mut(&result_blade) {
                    *existing_coeff = simplify(Expr::Add(
                        Box::new(existing_coeff.clone()),
                        Box::new(final_coeff),
                    ));
                } else {
                    result.terms.insert(result_blade, final_coeff);
                }
            }
        }
        result
    }

    /// Helper to compute the product of two basis blades.
    /// Returns (sign, metric_scalar, resulting_blade)
    pub(crate) fn blade_product(&self, b1: u32, b2: u32) -> (f64, Expr, u32) {
        let b1_mut = b1;
        let mut sign = 1.0;
        // Commutation sign
        for i in 0..32 {
            if (b2 >> i) & 1 == 1 {
                let swaps = (b1_mut >> (i + 1)).count_ones();
                if !swaps.is_multiple_of(2) {
                    sign *= -1.0;
                }
            }
        }

        let common_blades = b1 & b2;
        let mut metric_scalar = Expr::BigInt(BigInt::one());
        for i in 0..32 {
            if (common_blades >> i) & 1 == 1 {
                let (p, q, _r) = self.signature;
                let metric = if i < p {
                    1i64
                } else if i < p + q {
                    -1i64
                } else {
                    0i64
                };
                metric_scalar = simplify(Expr::Mul(
                    Box::new(metric_scalar),
                    Box::new(Expr::BigInt(BigInt::from(metric))),
                ));
            }
        }

        (sign, metric_scalar, b1 ^ b2)
    }

    /// Extracts all terms of a specific grade from the multivector.
    ///
    /// A multivector is a sum of blades of different grades (scalars are grade 0,
    /// vectors are grade 1, bivectors are grade 2, etc.). This function filters
    /// the multivector to keep only the terms corresponding to the desired grade.
    ///
    /// # Arguments
    /// * `grade` - The grade to project onto (e.g., 0 for scalar, 1 for vector).
    ///
    /// # Returns
    /// A new `Multivector` containing only the terms of the specified grade.
    pub fn grade_projection(&self, grade: u32) -> Multivector {
        let mut result = Multivector::new(self.signature);
        for (blade, coeff) in &self.terms {
            if blade.count_ones() == grade {
                result.terms.insert(*blade, coeff.clone());
            }
        }
        result
    }

    /// Computes the outer (or wedge) product of this multivector with another.
    ///
    /// The outer product `A ∧ B` produces a new blade representing the subspace
    /// spanned by the subspaces of A and B. It is grade-increasing: `grade(A ∧ B) = grade(A) + grade(B)`.
    /// It is defined in terms of the geometric product as the grade-sum part:
    /// `A ∧ B = <A B>_{r+s}` where `r=grade(A)` and `s=grade(B)`.
    ///
    /// # Arguments
    /// * `other` - The `Multivector` to compute the outer product with.
    ///
    /// # Returns
    /// A new `Multivector` representing the outer product.
    pub fn outer_product(&self, other: &Multivector) -> Multivector {
        let mut result = Multivector::new(self.signature);
        for r in 0..=self.signature.0 + self.signature.1 {
            for s in 0..=other.signature.0 + other.signature.1 {
                if r + s > self.signature.0 + self.signature.1 {
                    continue;
                }
                let term = self
                    .grade_projection(r)
                    .geometric_product(&other.grade_projection(s));
                result = result + term.grade_projection(r + s);
            }
        }
        result
    }

    /// Computes the inner (or left contraction) product of this multivector with another.
    ///
    /// The inner product `A . B` is a grade-decreasing operation. It is defined in terms
    /// of the geometric product as the grade-difference part:
    /// `A . B = <A B>_{s-r}` where `r=grade(A)` and `s=grade(B)`.
    ///
    /// # Arguments
    /// * `other` - The `Multivector` to compute the inner product with.
    ///
    /// # Returns
    /// A new `Multivector` representing the inner product.
    pub fn inner_product(&self, other: &Multivector) -> Multivector {
        let mut result = Multivector::new(self.signature);
        for r in 0..=self.signature.0 + self.signature.1 {
            for s in 0..=other.signature.0 + other.signature.1 {
                if s < r {
                    continue;
                }
                let term = self
                    .grade_projection(r)
                    .geometric_product(&other.grade_projection(s));
                result = result + term.grade_projection(s - r);
            }
        }
        result
    }

    /// Computes the reverse of the multivector.
    ///
    /// The reverse operation is found by reversing the order of the vectors in each basis blade.
    /// This results in a sign change for any blade `B` depending on its grade `k`:
    /// `reverse(B) = (-1)^(k*(k-1)/2) * B`.
    ///
    /// # Returns
    /// A new `Multivector` representing the reversed multivector.
    pub fn reverse(&self) -> Multivector {
        let mut result = Multivector::new(self.signature);
        for (blade, coeff) in &self.terms {
            let grade = blade.count_ones();
            let sign = if (grade * (grade - 1) / 2) % 2 == 0 {
                1i64
            } else {
                -1i64
            };
            result.terms.insert(
                *blade,
                simplify(Expr::Mul(
                    Box::new(Expr::BigInt(BigInt::from(sign))),
                    Box::new(coeff.clone()),
                )),
            );
        }
        result
    }
}

impl Add for Multivector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.clone();
        for (blade, coeff) in rhs.terms {
            if let Some(existing_coeff) = result.terms.get_mut(&blade) {
                *existing_coeff =
                    simplify(Expr::Add(Box::new(existing_coeff.clone()), Box::new(coeff)));
            } else {
                result.terms.insert(blade, coeff);
            }
        }
        result
    }
}

impl Sub for Multivector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut result = self.clone();
        for (blade, coeff) in rhs.terms {
            if let Some(existing_coeff) = result.terms.get_mut(&blade) {
                *existing_coeff =
                    simplify(Expr::Sub(Box::new(existing_coeff.clone()), Box::new(coeff)));
            } else {
                result
                    .terms
                    .insert(blade, simplify(Expr::Neg(Box::new(coeff))));
            }
        }
        result
    }
}

impl Mul<Expr> for Multivector {
    type Output = Self;
    fn mul(self, scalar: Expr) -> Self {
        let mut result = self.clone();
        for (_, coeff) in result.terms.iter_mut() {
            *coeff = simplify(Expr::Mul(Box::new(coeff.clone()), Box::new(scalar.clone())));
        }
        result
    }
}
