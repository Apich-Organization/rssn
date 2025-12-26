//! # Numerical Geometric Algebra (3D)
//!
//! This module provides a `Multivector3D` struct for numerical computations
//! in 3D Geometric Algebra (`G_3`). It implements the geometric product and
//! standard arithmetic operations for multivectors with `f64` components.

use std::ops::Add;
use std::ops::Neg;
use std::ops::Sub;

use serde::Deserialize;
use serde::Serialize;

/// Represents a multivector in 3D Geometric Algebra (`G_3`).
/// Components are: 1 (scalar), e1, e2, e3 (vectors), e12, e23, e31 (bivectors), e123 (pseudoscalar)
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]

pub struct Multivector3D {
    pub s : f64,
    pub v1 : f64,
    pub v2 : f64,
    pub v3 : f64,
    pub b12 : f64,
    pub b23 : f64,
    pub b31 : f64,
    pub pss : f64,
}

impl Add for Multivector3D {
    type Output = Self;

    /// Performs multivector addition.
    ///
    /// Addition is performed component-wise.

    fn add(
        self,
        rhs : Self,
    ) -> Self {

        Self {
            s : self.s + rhs.s,
            v1 : self.v1 + rhs.v1,
            v2 : self.v2 + rhs.v2,
            v3 : self.v3 + rhs.v3,
            b12 : self.b12 + rhs.b12,
            b23 : self.b23 + rhs.b23,
            b31 : self.b31 + rhs.b31,
            pss : self.pss + rhs.pss,
        }
    }
}

impl Sub for Multivector3D {
    type Output = Self;

    /// Performs multivector subtraction.
    ///
    /// Subtraction is performed component-wise.

    fn sub(
        self,
        rhs : Self,
    ) -> Self {

        Self {
            s : self.s - rhs.s,
            v1 : self.v1 - rhs.v1,
            v2 : self.v2 - rhs.v2,
            v3 : self.v3 - rhs.v3,
            b12 : self.b12 - rhs.b12,
            b23 : self.b23 - rhs.b23,
            b31 : self.b31 - rhs.b31,
            pss : self.pss - rhs.pss,
        }
    }
}

impl Neg for Multivector3D {
    type Output = Self;

    /// Performs multivector negation.
    ///
    /// Negation is performed component-wise.

    fn neg(self) -> Self {

        Self {
            s : -self.s,
            v1 : -self.v1,
            v2 : -self.v2,
            v3 : -self.v3,
            b12 : -self.b12,
            b23 : -self.b23,
            b31 : -self.b31,
            pss : -self.pss,
        }
    }
}

/// Implements the geometric product for `Multivector3D`.
///
/// The geometric product is the fundamental product in geometric algebra.
/// It combines the inner (dot) and outer (wedge) products.
/// This implementation uses the full multiplication table for `G_3`,
/// based on `e_i*e_j = -e_j*e_i` for `i != j` and `e_i*e_i = 1`.

impl Multivector3D {
    /// Creates a new `Multivector3D` with all components.
    #[allow(clippy::too_many_arguments)]

    #[must_use] 
    pub const fn new(
        s : f64,
        v1 : f64,
        v2 : f64,
        v3 : f64,
        b12 : f64,
        b23 : f64,
        b31 : f64,
        pss : f64,
    ) -> Self {

        Self {
            s,
            v1,
            v2,
            v3,
            b12,
            b23,
            b31,
            pss,
        }
    }

    /// Returns the reverse of the multivector.
    ///
    /// The reverse operation reverses the order of products of basis vectors.
    /// For G3: s, v1, v2, v3 are unchanged; b12, b23, b31, pss are negated.
    /// Wait, `reverse(e_i` `e_j`) = `e_j` `e_i` = -`e_i` `e_j`. So bivectors are negated.
    /// `reverse(e_1` `e_2` `e_3`) = `e_3` `e_2` `e_1` = -`e_1` `e_2` `e_3`. So pseudoscalar is negated.
    #[must_use]

    pub fn reverse(self) -> Self {

        Self {
            s : self.s,
            v1 : self.v1,
            v2 : self.v2,
            v3 : self.v3,
            b12 : -self.b12,
            b23 : -self.b23,
            b31 : -self.b31,
            pss : -self.pss,
        }
    }

    /// Returns the Clifford conjugate of the multivector.
    ///
    /// Conjugation combines reversal and grade involution.
    #[must_use]

    pub fn conjugate(self) -> Self {

        Self {
            s : self.s,
            v1 : -self.v1,
            v2 : -self.v2,
            v3 : -self.v3,
            b12 : -self.b12,
            b23 : -self.b23,
            b31 : -self.b31,
            pss : self.pss,
        }
    }

    /// Returns the squared norm of the multivector (A * reverse(A))_s.
    #[must_use]

    pub fn norm_sq(self) -> f64 {

        self.s.mul_add(self.s, self.v1 * self.v1)
            + self.v2 * self.v2
            + self.v3 * self.v3
            + self.b12 * self.b12
            + self.b23 * self.b23
            + self.b31 * self.b31
            + self.pss * self.pss
    }

    /// Returns the norm of the multivector.
    #[must_use]

    pub fn norm(self) -> f64 {

        self.norm_sq()
            .sqrt()
    }

    /// Returns the inverse of the multivector, if it exists.
    ///
    /// For simple multivectors (like vectors or blades), this is `A_rev` / |A|^2.
    #[must_use]

    pub fn inv(self) -> Option<Self> {

        let n2 = self.norm_sq();

        if n2.abs() < f64::EPSILON {

            None
        } else {

            let rev = self.reverse();

            Some(Self {
                s : rev.s / n2,
                v1 : rev.v1 / n2,
                v2 : rev.v2 / n2,
                v3 : rev.v3 / n2,
                b12 : rev.b12 / n2,
                b23 : rev.b23 / n2,
                b31 : rev.b31 / n2,
                pss : rev.pss / n2,
            })
        }
    }

    /// Performs the outer (wedge) product.
    #[must_use]

    pub fn wedge(
        self,
        rhs : Self,
    ) -> Self {

        // The wedge product is the grade-increasing part of the geometric product.
        // A ^ B = sum_{r,s} <<a>_r <b>_s>_{r+s}
        Self {
            s : self.s * rhs.s,
            v1 : self.s.mul_add(rhs.v1, self.v1 * rhs.s),
            v2 : self.s.mul_add(rhs.v2, self.v2 * rhs.s),
            v3 : self.s.mul_add(rhs.v3, self.v3 * rhs.s),
            b12 : self.s.mul_add(rhs.b12, self.b12 * rhs.s) + self.v1 * rhs.v2 - self.v2 * rhs.v1,
            b23 : self.s.mul_add(rhs.b23, self.b23 * rhs.s) + self.v2 * rhs.v3 - self.v3 * rhs.v2,
            b31 : self.s.mul_add(rhs.b31, self.b31 * rhs.s) + self.v3 * rhs.v1 - self.v1 * rhs.v3,
            pss : self.v1.mul_add(rhs.b23, self.s.mul_add(
                rhs.pss,
                self.pss * rhs.s,
            ))
                + self.v2 * rhs.b31
                + self.v3 * rhs.b12
                + self.b12 * rhs.v3
                + self.b23 * rhs.v1
                + self.b31 * rhs.v2,
        }
    }

    /// Performs the inner (dot) product.
    #[must_use]

    pub fn dot(
        self,
        rhs : Self,
    ) -> Self {

        // The inner product is the grade-decreasing part of the geometric product.
        // A . B = sum_{r,s} <<a>_r <b>_s>_{|r-s|}
        Self {
            s : self.v2.mul_add(rhs.v2, self.s.mul_add(
                rhs.s,
                self.v1 * rhs.v1,
            ))
                + self.v3 * rhs.v3
                - self.b12 * rhs.b12
                - self.b23 * rhs.b23
                - self.b31 * rhs.b31
                - self.pss * rhs.pss,
            v1 : self.s.mul_add(rhs.v1, self.v1 * rhs.s) - self.v2 * rhs.b12
                + self.v3 * rhs.b31
                + self.b12 * rhs.v2
                - self.b31 * rhs.v3
                - self.b23 * rhs.pss
                - self.pss * rhs.b23,
            v2 : self.s.mul_add(rhs.v2, self.v1 * rhs.b12) + self.v2 * rhs.s
                - self.v3 * rhs.b23
                - self.b12 * rhs.v1
                + self.b23 * rhs.v3
                - self.b31 * rhs.pss
                - self.pss * rhs.b31,
            v3 : self.s.mul_add(rhs.v3, -(self.v1 * rhs.b31)) + self.v2 * rhs.b23 + self.v3 * rhs.s
                - self.b12 * rhs.pss
                - self.b23 * rhs.v2
                + self.b31 * rhs.v1
                - self.pss * rhs.b12,
            b12 : self.s.mul_add(rhs.b12, self.b12 * rhs.s) - self.b23 * rhs.b31 + self.b31 * rhs.b23,
            b23 : self.s.mul_add(rhs.b23, self.b23 * rhs.s) + self.b12 * rhs.b31 - self.b31 * rhs.b12,
            b31 : self.s.mul_add(rhs.b31, self.b31 * rhs.s) - self.b12 * rhs.b23 + self.b23 * rhs.b12,
            pss : self.s.mul_add(rhs.pss, self.pss * rhs.s),
        }
    }
}

impl std::ops::Mul for Multivector3D {
    type Output = Self;

    fn mul(
        self,
        rhs : Self,
    ) -> Self::Output {

        Self {
            s : self.v2.mul_add(rhs.v2, self.s.mul_add(
                rhs.s,
                self.v1 * rhs.v1,
            ))
                + self.v3 * rhs.v3
                - self.b12 * rhs.b12
                - self.b23 * rhs.b23
                - self.b31 * rhs.b31
                - self.pss * rhs.pss,
            v1 : self.v2.mul_add(-rhs.b12, self.s.mul_add(
                rhs.v1,
                self.v1 * rhs.s,
            ))
                + self.v3 * rhs.b31
                + self.b12 * rhs.v2
                - self.b31 * rhs.v3
                - self.b23 * rhs.pss
                - self.pss * rhs.b23,
            v2 : self.v1.mul_add(rhs.b12, self.s.mul_add(
                rhs.v2,
                self.v2 * rhs.s,
            ))
                - self.v3 * rhs.b23
                - self.b12 * rhs.v1
                + self.b23 * rhs.v3
                - self.b31 * rhs.pss
                - self.pss * rhs.b31,
            v3 : self.v1.mul_add(-rhs.b31, self.s.mul_add(
                rhs.v3,
                self.v3 * rhs.s,
            ))
                + self.v2 * rhs.b23
                + self.b31 * rhs.v1
                - self.b23 * rhs.v2
                - self.b12 * rhs.pss
                - self.pss * rhs.b12,
            b12 : self.v1.mul_add(rhs.v2, self.s.mul_add(
                rhs.b12,
                self.b12 * rhs.s,
            ))
                - self.v2 * rhs.v1
                + self.v3 * rhs.pss
                + self.pss * rhs.v3
                - self.b23 * rhs.b31
                + self.b31 * rhs.b23,
            b23 : self.v2.mul_add(rhs.v3, self.s.mul_add(
                rhs.b23,
                self.b23 * rhs.s,
            ))
                - self.v3 * rhs.v2
                + self.v1 * rhs.pss
                + self.pss * rhs.v1
                - self.b31 * rhs.b12
                + self.b12 * rhs.b31,
            b31 : self.v3.mul_add(rhs.v1, self.s.mul_add(
                rhs.b31,
                self.b31 * rhs.s,
            ))
                - self.v1 * rhs.v3
                + self.v2 * rhs.pss
                + self.pss * rhs.v2
                - self.b12 * rhs.b23
                + self.b23 * rhs.b12,
            pss : self.v1.mul_add(rhs.b23, self.s.mul_add(
                rhs.pss,
                self.pss * rhs.s,
            ))
                + self.v2 * rhs.b31
                + self.v3 * rhs.b12
                + self.b12 * rhs.v3
                + self.b23 * rhs.v1
                + self.b31 * rhs.v2,
        }
    }
}
