//! # Numerical Geometric Algebra (3D)
//!
//! This module provides a `Multivector3D` struct for numerical computations
//! in 3D Geometric Algebra (G_3). It implements the geometric product and
//! standard arithmetic operations for multivectors with `f64` components.
use std::ops::{Add, Neg, Sub};
/// Represents a multivector in 3D Geometric Algebra (G_3).
/// Components are: 1 (scalar), e1, e2, e3 (vectors), e12, e23, e31 (bivectors), e123 (pseudoscalar)
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Multivector3D {
    pub s: f64,
    pub v1: f64,
    pub v2: f64,
    pub v3: f64,
    pub b12: f64,
    pub b23: f64,
    pub b31: f64,
    pub pss: f64,
}
impl Add for Multivector3D {
    type Output = Self;
    /// Performs multivector addition.
    ///
    /// Addition is performed component-wise.
    fn add(self, rhs: Self) -> Self {
        Self {
            s: self.s + rhs.s,
            v1: self.v1 + rhs.v1,
            v2: self.v2 + rhs.v2,
            v3: self.v3 + rhs.v3,
            b12: self.b12 + rhs.b12,
            b23: self.b23 + rhs.b23,
            b31: self.b31 + rhs.b31,
            pss: self.pss + rhs.pss,
        }
    }
}
impl Sub for Multivector3D {
    type Output = Self;
    /// Performs multivector subtraction.
    ///
    /// Subtraction is performed component-wise.
    fn sub(self, rhs: Self) -> Self {
        Self {
            s: self.s - rhs.s,
            v1: self.v1 - rhs.v1,
            v2: self.v2 - rhs.v2,
            v3: self.v3 - rhs.v3,
            b12: self.b12 - rhs.b12,
            b23: self.b23 - rhs.b23,
            b31: self.b31 - rhs.b31,
            pss: self.pss - rhs.pss,
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
            s: -self.s,
            v1: -self.v1,
            v2: -self.v2,
            v3: -self.v3,
            b12: -self.b12,
            b23: -self.b23,
            b31: -self.b31,
            pss: -self.pss,
        }
    }
}
/// Implements the geometric product for Multivector3D.
///
/// The geometric product is the fundamental product in geometric algebra.
/// It combines the inner (dot) and outer (wedge) products.
/// This implementation uses the full multiplication table for G_3,
/// based on `e_i*e_j = -e_j*e_i` for `i != j` and `e_i*e_i = 1`.
impl std::ops::Mul for Multivector3D {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            s: self.s * rhs.s + self.v1 * rhs.v1 + self.v2 * rhs.v2 + self.v3 * rhs.v3
                - self.b12 * rhs.b12
                - self.b23 * rhs.b23
                - self.b31 * rhs.b31
                - self.pss * rhs.pss,
            v1: self.s * rhs.v1 + self.v1 * rhs.s - self.v2 * rhs.b12
                + self.v3 * rhs.b31
                + self.b12 * rhs.v2
                - self.b31 * rhs.v3
                + self.b23 * rhs.pss
                - self.pss * rhs.b23,
            v2: self.s * rhs.v2 + self.v2 * rhs.s + self.v1 * rhs.b12
                - self.v3 * rhs.b23
                - self.b12 * rhs.v1
                + self.b23 * rhs.v3
                - self.b31 * rhs.pss
                + self.pss * rhs.b31,
            v3: self.s * rhs.v3 + self.v3 * rhs.s - self.v1 * rhs.b31
                + self.v2 * rhs.b23
                + self.b31 * rhs.v1
                - self.b23 * rhs.v2
                + self.b12 * rhs.pss
                - self.pss * rhs.b12,
            b12: self.s * rhs.b12 + self.b12 * rhs.s + self.v1 * rhs.v2 - self.v2 * rhs.v1
                + self.v3 * rhs.pss
                + self.pss * rhs.v3
                - self.b23 * rhs.b31
                + self.b31 * rhs.b23,
            b23: self.s * rhs.b23 + self.b23 * rhs.s + self.v2 * rhs.v3
                - self.v3 * rhs.v2
                - self.v1 * rhs.pss
                - self.pss * rhs.v1
                - self.b31 * rhs.b12
                + self.b12 * rhs.b31,
            b31: self.s * rhs.b31 + self.b31 * rhs.s + self.v3 * rhs.v1 - self.v1 * rhs.v3
                + self.v2 * rhs.pss
                + self.pss * rhs.v2
                - self.b12 * rhs.b23
                + self.b23 * rhs.b12,
            pss: self.s * rhs.pss
                + self.pss * rhs.s
                + self.v1 * rhs.b23
                + self.v2 * rhs.b31
                + self.v3 * rhs.b12
                + self.b12 * rhs.v3
                + self.b23 * rhs.v1
                + self.b31 * rhs.v2,
        }
    }
}
