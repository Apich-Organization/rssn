//! # Geometric Algebra
//!
//! This module provides tools for computations in Clifford and Geometric Algebra.

use std::collections::BTreeMap;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

/// Represents a multivector in a Clifford algebra.
///
/// The basis blades are represented by a bitmask. E.g., in 3D:
/// 001 (1) -> e1, 010 (2) -> e2, 100 (4) -> e3
/// 011 (3) -> e12, 101 (5) -> e13, 110 (6) -> e23
/// 111 (7) -> e123 (pseudoscalar)
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub struct Multivector {
    /// A map from the basis blade bitmask to its coefficient.
    pub terms: BTreeMap<u32, Expr>,
    /// The signature of the algebra, e.g., (p, q, r) for (`e_i^2` = +1, `e_j^2` = -1, `e_k^2` = 0)
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
    #[must_use]

    pub const fn new(
        signature: (u32, u32, u32)
    ) -> Self {

        Self {
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
    #[must_use]

    pub fn scalar(
        signature: (u32, u32, u32),
        value: Expr,
    ) -> Self {

        let mut terms = BTreeMap::new();

        terms.insert(0, value);

        Self {
            terms,
            signature,
        }
    }

    /// Creates a new multivector representing a vector (grade-1 element).
    ///
    /// # Arguments
    /// * `signature` - The signature of the algebra `(p, q, r)`.
    /// * `components` - A vector of coefficients for each basis vector.
    ///
    /// # Returns
    /// A `Multivector` representing the vector.
    #[must_use]

    pub fn vector(
        signature: (u32, u32, u32),
        components: Vec<Expr>,
    ) -> Self {

        let mut terms = BTreeMap::new();

        for (i, coeff) in components
            .into_iter()
            .enumerate()
        {

            terms.insert(1 << i, coeff);
        }

        Self {
            terms,
            signature,
        }
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
    #[must_use]

    pub fn geometric_product(
        &self,
        other: &Self,
    ) -> Self {

        let mut result =
            Self::new(self.signature);

        for (blade1, coeff1) in
            &self.terms
        {

            for (blade2, coeff2) in
                &other.terms
            {

                let (
                    sign,
                    metric_scalar,
                    result_blade,
                ) = self.blade_product(
                    *blade1,
                    *blade2,
                );

                let new_coeff =
                    simplify(
                        &Expr::new_mul(
                            coeff1
                                .clone(
                                ),
                            coeff2
                                .clone(
                                ),
                        ),
                    );

                let signed_coeff = simplify(&Expr::new_mul(
                    Expr::Constant(sign),
                    new_coeff,
                ));

                let final_coeff = simplify(&Expr::new_mul(
                    signed_coeff,
                    metric_scalar,
                ));

                if let Some(
                    existing_coeff,
                ) = result
                    .terms
                    .get_mut(
                        &result_blade,
                    )
                {

                    *existing_coeff = simplify(&Expr::new_add(
                        existing_coeff.clone(),
                        final_coeff,
                    ));
                } else {

                    result
                        .terms
                        .insert(
                        result_blade,
                        final_coeff,
                    );
                }
            }
        }

        result.prune_zeros();

        result
    }

    /// Helper to prune terms with zero coefficients.

    fn prune_zeros(&mut self) {

        self.terms
            .retain(|_, coeff| {

                match coeff {
                    | Expr::Constant(c) => c.abs() > f64::EPSILON,
                    | Expr::BigInt(b) => !b.is_zero(),
                    | Expr::Rational(r) => !r.is_zero(),
                    | Expr::Dag(node) => {
                        match node.to_expr() { Ok(expr) => {

                            match expr {
                                | Expr::Constant(c) => c.abs() > f64::EPSILON,
                                | Expr::BigInt(b) => !b.is_zero(),
                                | Expr::Rational(r) => !r.is_zero(),
                                | _ => true,
                            }
                        } _ => {

                            true
                        }}
                    },
                    | _ => true, // Keep symbolic terms
                }
            });
    }

    /// Helper to compute the product of two basis blades.
    /// Returns (sign, `metric_scalar`, `resulting_blade`)

    pub(crate) fn blade_product(
        &self,
        b1: u32,
        b2: u32,
    ) -> (f64, Expr, u32) {

        let b1_mut = b1;

        let mut sign = 1.0;

        for i in 0 .. 32 {

            if (b2 >> i) & 1 == 1 {

                let swaps = (b1_mut
                    >> (i + 1))
                    .count_ones();

                if !swaps
                    .is_multiple_of(2)
                {

                    sign *= -1.0;
                }
            }
        }

        let common_blades = b1 & b2;

        let mut metric_scalar =
            Expr::BigInt(BigInt::one());

        for i in 0 .. 32 {

            if (common_blades >> i) & 1
                == 1
            {

                let (p, q, _r) =
                    self.signature;

                let metric = if i < p {

                    1i64
                } else if i < p + q {

                    -1i64
                } else {

                    0i64
                };

                metric_scalar = simplify(&Expr::new_mul(
                    metric_scalar,
                    Expr::BigInt(BigInt::from(metric)),
                ));
            }
        }

        (
            sign,
            metric_scalar,
            b1 ^ b2,
        )
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
    #[must_use]

    pub fn grade_projection(
        &self,
        grade: u32,
    ) -> Self {

        let mut result =
            Self::new(self.signature);

        for (blade, coeff) in
            &self.terms
        {

            if blade.count_ones()
                == grade
            {

                result.terms.insert(
                    *blade,
                    coeff.clone(),
                );
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
    #[must_use]

    pub fn outer_product(
        &self,
        other: &Self,
    ) -> Self {

        let mut result =
            Self::new(self.signature);

        for r in 0 ..= self.signature.0
            + self.signature.1
        {

            for s in
                0 ..= other.signature.0
                    + other.signature.1
            {

                if r + s
                    > self.signature.0
                        + self
                            .signature
                            .1
                {

                    continue;
                }

                let term = self
                    .grade_projection(r)
                    .geometric_product(&other.grade_projection(s));

                result = result + term.grade_projection(r + s);
            }
        }

        result.prune_zeros();

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
    #[must_use]

    pub fn inner_product(
        &self,
        other: &Self,
    ) -> Self {

        let mut result =
            Self::new(self.signature);

        for r in 0 ..= self.signature.0
            + self.signature.1
        {

            for s in
                0 ..= other.signature.0
                    + other.signature.1
            {

                if s < r {

                    continue;
                }

                let term = self
                    .grade_projection(r)
                    .geometric_product(&other.grade_projection(s));

                result = result + term.grade_projection(s - r);
            }
        }

        result.prune_zeros();

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
    #[must_use]

    pub fn reverse(&self) -> Self {

        let mut result =
            Self::new(self.signature);

        for (blade, coeff) in
            &self.terms
        {

            let grade = i64::from(
                blade.count_ones(),
            );

            let sign = if (grade
                * (grade - 1)
                / 2)
                % 2
                == 0
            {

                1i64
            } else {

                -1i64
            };

            result.terms.insert(
                *blade,
                simplify(&Expr::new_mul(
                    Expr::BigInt(BigInt::from(sign)),
                    coeff.clone(),
                )),
            );
        }

        result
    }

    /// Computes the magnitude (norm) of the multivector.
    ///
    /// The magnitude is defined as `sqrt(M * reverse(M))` where the result
    /// should be a scalar.
    ///
    /// # Returns
    /// An `Expr` representing the magnitude.
    #[must_use]

    pub fn magnitude(&self) -> Expr {

        let product = self
            .geometric_product(
                &self.reverse(),
            );

        let scalar_part =
            product.grade_projection(0);

        scalar_part.terms.get(&0).map_or_else(
            || Expr::Constant(0.0),
            |scalar_coeff| {
                simplify(&Expr::new_sqrt(scalar_coeff.clone()))
            },
        )
    }

    /// Computes the dual of the multivector with respect to the pseudoscalar.
    ///
    /// The dual is defined as `M * I^(-1)` where `I` is the pseudoscalar.
    ///
    /// # Returns
    /// A new `Multivector` representing the dual.
    #[must_use]

    pub fn dual(&self) -> Self {

        let dimension =
            self.signature.0
                + self.signature.1
                + self.signature.2;

        let pseudoscalar_blade =
            (1 << dimension) - 1;

        // Create pseudoscalar multivector
        let mut pseudoscalar =
            Self::new(self.signature);

        pseudoscalar
            .terms
            .insert(
                pseudoscalar_blade,
                Expr::Constant(1.0),
            );

        // Compute dual as M * I^(-1)
        // For simplicity, we use M * I (which works for many cases)
        self.geometric_product(
            &pseudoscalar,
        )
    }

    /// Normalizes the multivector to unit magnitude.
    ///
    /// # Returns
    /// A new `Multivector` with unit magnitude.
    #[must_use]

    pub fn normalize(&self) -> Self {

        let mag = self.magnitude();

        let inv_mag = Expr::new_div(
            Expr::Constant(1.0),
            mag,
        );

        self.clone() * inv_mag
    }
}

impl Add for Multivector {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        let mut result = self;

        for (blade, coeff) in rhs.terms
        {

            if let Some(
                existing_coeff,
            ) = result
                .terms
                .get_mut(&blade)
            {

                *existing_coeff = simplify(&Expr::new_add(
                    existing_coeff.clone(),
                    coeff,
                ));
            } else {

                result.terms.insert(
                    blade, coeff,
                );
            }
        }

        result.prune_zeros();

        result
    }
}

impl Sub for Multivector {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        let mut result = self;

        for (blade, coeff) in rhs.terms
        {

            if let Some(
                existing_coeff,
            ) = result
                .terms
                .get_mut(&blade)
            {

                *existing_coeff = simplify(&Expr::new_sub(
                    existing_coeff.clone(),
                    coeff,
                ));
            } else {

                result.terms.insert(
                    blade,
                    simplify(
                        &Expr::new_neg(
                            coeff,
                        ),
                    ),
                );
            }
        }

        result.prune_zeros();

        result
    }
}

impl Mul<Expr> for Multivector {
    type Output = Self;

    fn mul(
        self,
        scalar: Expr,
    ) -> Self {

        let mut result = self;

        for coeff in result
            .terms
            .values_mut()
        {

            *coeff = simplify(
                &Expr::new_mul(
                    coeff.clone(),
                    scalar.clone(),
                ),
            );
        }

        result.prune_zeros();

        result
    }
}
