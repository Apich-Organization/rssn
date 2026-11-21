//! # Numerical Geometric Algebra (3D)
//!
//! This module provides a `Multivector3D` struct for numerical computations
//! in 3D Geometric Algebra (`G_3`). It implements the geometric product and
//! standard arithmetic operations for multivectors with `f64` components.
//!
//! # Limitations and coverage
//!
//! * Fixed to the 3-D algebra `G_3`: one scalar, three vectors, three
//!   bivectors and one pseudoscalar. There is no general `Cl(p, q)` or
//!   higher-dimensional Clifford algebra and no dimension parameter.
//! * Coefficients are `f64` and the struct is `Copy`. Arithmetic uses fixed
//!   component formulas, with exact intermediates available for difficult inversions.
//! * [`Multivector3D::inv`] returns `None` for a singular element or an
//!   inverse outside the finite `f64` range. There is no exponential,
//!   logarithm, rotor interpolation or trigonometric function of a
//!   multivector, and no outermorphism/matrix representation.
//! * No SIMD or batch evaluation: the type is intended for single-element
//!   numerical work, not for large geometric-algebra arrays.

use std::ops::Add;
use std::ops::Neg;
use std::ops::Sub;

use num_bigint::BigInt;
use num_complex::Complex;
use num_rational::BigRational;
use num_traits::Float;
use num_traits::ToPrimitive;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

/// Represents a multivector in 3D Geometric Algebra (`G_3`).
/// Components are: 1 (scalar), e1, e2, e3 (vectors), e12, e23, e31 (bivectors), e123 (pseudoscalar)
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Multivector3D {
    /// Scalar component.
    pub s: f64,
    /// Vector e1 component.
    pub v1: f64,
    /// Vector e2 component.
    pub v2: f64,
    /// Vector e3 component.
    pub v3: f64,
    /// Bivector e12 component.
    pub b12: f64,
    /// Bivector e23 component.
    pub b23: f64,
    /// Bivector e31 component.
    pub b31: f64,
    /// Pseudoscalar e123 component.
    pub pss: f64,
}

impl Add for Multivector3D {
    type Output = Self;

    /// Performs multivector addition.
    ///
    /// Addition is performed component-wise.
    fn add(
        self,
        rhs: Self,
    ) -> Self {
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
    fn sub(
        self,
        rhs: Self,
    ) -> Self {
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
        s: f64,
        v1: f64,
        v2: f64,
        v3: f64,
        b12: f64,
        b23: f64,
        b31: f64,
        pss: f64,
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
            s: self.s,
            v1: self.v1,
            v2: self.v2,
            v3: self.v3,
            b12: -self.b12,
            b23: -self.b23,
            b31: -self.b31,
            pss: -self.pss,
        }
    }

    /// Returns the Clifford conjugate of the multivector.
    ///
    /// Conjugation combines reversal and grade involution.
    #[must_use]
    pub fn conjugate(self) -> Self {
        Self {
            s: self.s,
            v1: -self.v1,
            v2: -self.v2,
            v3: -self.v3,
            b12: -self.b12,
            b23: -self.b23,
            b31: -self.b31,
            pss: self.pss,
        }
    }

    /// Returns the squared norm of the multivector (A * reverse(A))_s.
    #[must_use]
    pub fn norm_sq(self) -> f64 {
        self.pss.mul_add(
            self.pss,
            self.b31.mul_add(
                self.b31,
                self.b23.mul_add(
                    self.b23,
                    self.b12.mul_add(
                        self.b12,
                        self.v3.mul_add(
                            self.v3,
                            self.v2
                                .mul_add(self.v2, self.s.mul_add(self.s, self.v1 * self.v1)),
                        ),
                    ),
                ),
            ),
        )
    }

    /// Returns the norm of the multivector.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Returns the inverse of a general `Cl(3,0)` multivector.
    ///
    /// The product with its Clifford conjugate is central (scalar plus
    /// pseudoscalar). Components of the input must be finite. Returns `None`
    /// if the determinant is exactly zero or an inverse component overflows.
    #[must_use]
    pub fn inv(self) -> Option<Self> {
        // Within these bounds the squared determinant and small components
        // of its reciprocal can be computed before applying the adjugate.
        const DIRECT_MIN: u64 = 895_u64 << 52; // 2^-128
        const DIRECT_MAX: u64 = 1151_u64 << 52; // 2^128
        const MIN_NORMALIZED: f64 = f64::from_bits(623_u64 << 52); // 2^-400
        const OVERFLOW_NEIGHBORHOOD: f64 = f64::from_bits(2046_u64 << 52); // 2^1023

        let mut largest = 0;
        let mut smallest = u64::MAX;
        for value in self.coefficients() {
            // Finite nonnegative binary64 encodings are ordered by magnitude.
            let bits = value.to_bits() & 0x7fff_ffff_ffff_ffff;
            largest = largest.max(bits);
            if bits != 0 {
                smallest = smallest.min(bits);
            }
        }
        if largest == 0 || largest >= f64::INFINITY.to_bits() {
            return None;
        }
        let max = f64::from_bits(largest);
        if smallest >= DIRECT_MIN && largest <= DIRECT_MAX {
            return self.inverse_fast(max).or_else(|| self.inverse_exact());
        }

        // Normalize by a power of two, even for an entirely subnormal input.
        let exponent = (largest & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE.to_bits());
        let factor = f64::from_bits(exponent).recip();
        let relative_min = f64::from_bits(smallest) * factor;
        if relative_min >= MIN_NORMALIZED
            && let Some(inverse) = self.scale_components(factor).inverse_fast(max * factor)
        {
            let inverse = inverse.scale_components(factor);
            // Approximation in the top binade can round down to a finite
            // value even when the exact inverse rounds to infinity.
            if inverse
                .coefficients()
                .iter()
                .all(|value| value.abs() < OVERFLOW_NEIGHBORHOOD)
            {
                return Some(inverse);
            }
        }
        self.inverse_exact()
    }

    /// Multiplies by the central factor `scalar + pseudoscalar * e123`.
    ///
    /// The `e123` blade commutes with every multivector and squares to `-1`,
    /// so this operation combines four independent pairs of coefficients.
    #[must_use]
    pub fn mul_central(
        self,
        scalar: f64,
        pseudoscalar: f64,
    ) -> Self {
        let multiply = |a: f64, b: f64| {
            (
                (-b).mul_add(pseudoscalar, a * scalar),
                a.mul_add(pseudoscalar, b * scalar),
            )
        };
        let (s, pss) = multiply(self.s, self.pss);
        let (v1, b23) = multiply(self.v1, self.b23);
        let (v2, b31) = multiply(self.v2, self.b31);
        let (v3, b12) = multiply(self.v3, self.b12);
        Self::new(s, v1, v2, v3, b12, b23, b31, pss)
    }

    const fn coefficients(self) -> [f64; 8] {
        [
            self.s, self.v1, self.v2, self.v3, self.b12, self.b23, self.b31, self.pss,
        ]
    }

    fn scale_components(
        self,
        scale: f64,
    ) -> Self {
        let [s, v1, v2, v3, b12, b23, b31, pss] = self.coefficients().map(|x| x * scale);
        Self::new(s, v1, v2, v3, b12, b23, b31, pss)
    }

    fn inverse_fast(
        self,
        max: f64,
    ) -> Option<Self> {
        let positive = self.b31.mul_add(
            self.b31,
            self.b23
                .mul_add(self.b23, self.b12.mul_add(self.b12, self.s * self.s)),
        );
        let negative = self.pss.mul_add(
            self.pss,
            self.v3
                .mul_add(self.v3, self.v2.mul_add(self.v2, self.v1 * self.v1)),
        );
        let real = positive - negative;
        let imaginary = 2.0
            * (-self.v3).mul_add(
                self.b12,
                (-self.v2).mul_add(self.b31, (-self.v1).mul_add(self.b23, self.s * self.pss)),
            );
        // Small determinants may be lost to cancellation. Leave a margin
        // above roundoff before taking the reciprocal in binary64.
        if real.abs().max(imaginary.abs()) <= max * max * (1.0 / 1_048_576.0) {
            let central_and_bivectors = self.pss.abs().to_bits()
                | self.b12.abs().to_bits()
                | self.b23.abs().to_bits()
                | self.b31.abs().to_bits();
            if central_and_bivectors == 0
                && [self.v1, self.v2, self.v3]
                    .into_iter()
                    .filter(|x| x.abs().to_bits() != 0)
                    .count()
                    == 1
            {
                let v = self.v1 + self.v2 + self.v3;
                // For one axis, factoring s²-v² computes the small
                // difference directly. Near the cancellation it is exact.
                let determinant = (self.s - v) * (self.s + v);
                if determinant.abs().to_bits() != 0 {
                    return Some(self.conjugate().scale_components(determinant.recip()));
                }
            }
            return None;
        }
        let norm = real.mul_add(real, imaginary * imaginary);
        Some(self.conjugate().mul_central(real / norm, -imaginary / norm))
    }

    fn inverse_exact(self) -> Option<Self> {
        // Express each input as an integer times a common power of two.
        // This retains exact cancellation and avoids a rational reduction
        // at each multiply. Only the final quotients require rounding.
        let parts = self.coefficients().map(|value| {
            let (mantissa, exponent, sign) = Float::integer_decode(value);
            let zeros = mantissa.trailing_zeros().min(53);
            (mantissa >> zeros, i32::from(exponent) + zeros as i32, sign)
        });
        let power = parts
            .iter()
            .filter(|part| part.0 != 0)
            .map(|part| part.1)
            .min()?;
        let [s, v1, v2, v3, b12, b23, b31, pss] = parts.map(|(mantissa, exponent, sign)| {
            if mantissa == 0 {
                BigInt::zero()
            } else {
                (BigInt::from(mantissa) << ((exponent - power) as usize)) * sign
            }
        });
        let (z, x, y, w) = (
            Complex::new(s, pss),
            Complex::new(v1, b23),
            Complex::new(v2, b31),
            Complex::new(v3, b12),
        );
        let determinant = &z * &z - &x * &x - &y * &y - &w * &w;
        let norm = determinant.norm_sqr();
        if norm.is_zero() {
            return None;
        }
        let denominator = norm << (power.max(0) as usize);
        let numerator_shift = (-power).max(0) as usize;
        let scale = determinant.conj();
        let quotient = |value: BigInt| {
            BigRational::new_raw(value << numerator_shift, denominator.clone())
                .to_f64()
                .filter(|x| x.is_finite())
        };
        let [a, b, c, d] = [z, -x, -y, -w].map(|part| part * &scale);
        Some(Self::new(
            quotient(a.re)?,
            quotient(b.re)?,
            quotient(c.re)?,
            quotient(d.re)?,
            quotient(d.im)?,
            quotient(b.im)?,
            quotient(c.im)?,
            quotient(a.im)?,
        ))
    }

    /// Performs the outer (wedge) product.
    #[must_use]
    pub fn wedge(
        self,
        rhs: Self,
    ) -> Self {
        // The wedge product is the grade-increasing part of the geometric product.
        // A ^ B = sum_{r,s} <<a>_r <b>_s>_{r+s}
        Self {
            s: self.s * rhs.s,
            v1: self.s.mul_add(rhs.v1, self.v1 * rhs.s),
            v2: self.s.mul_add(rhs.v2, self.v2 * rhs.s),
            v3: self.s.mul_add(rhs.v3, self.v3 * rhs.s),
            b12: self.v2.mul_add(
                -rhs.v1,
                self.v1
                    .mul_add(rhs.v2, self.s.mul_add(rhs.b12, self.b12 * rhs.s)),
            ),
            b23: self.v3.mul_add(
                -rhs.v2,
                self.v2
                    .mul_add(rhs.v3, self.s.mul_add(rhs.b23, self.b23 * rhs.s)),
            ),
            b31: self.v1.mul_add(
                -rhs.v3,
                self.v3
                    .mul_add(rhs.v1, self.s.mul_add(rhs.b31, self.b31 * rhs.s)),
            ),
            pss: self.b31.mul_add(
                rhs.v2,
                self.b23.mul_add(
                    rhs.v1,
                    self.b12.mul_add(
                        rhs.v3,
                        self.v3.mul_add(
                            rhs.b12,
                            self.v2.mul_add(
                                rhs.b31,
                                self.v1
                                    .mul_add(rhs.b23, self.s.mul_add(rhs.pss, self.pss * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
        }
    }

    /// Performs the inner (dot) product.
    #[must_use]
    pub fn dot(
        self,
        rhs: Self,
    ) -> Self {
        // The inner product is the grade-decreasing part of the geometric product.
        // A . B = sum_{r,s} <<a>_r <b>_s>_{|r-s|}
        Self {
            s: (-self.pss).mul_add(
                rhs.pss,
                (-self.b31).mul_add(
                    rhs.b31,
                    (-self.b23).mul_add(
                        rhs.b23,
                        self.b12.mul_add(
                            -rhs.b12,
                            self.v3.mul_add(
                                rhs.v3,
                                self.v2
                                    .mul_add(rhs.v2, self.s.mul_add(rhs.s, self.v1 * rhs.v1)),
                            ),
                        ),
                    ),
                ),
            ),
            v1: (-self.pss).mul_add(
                rhs.b23,
                (-self.b23).mul_add(
                    rhs.pss,
                    (-self.b31).mul_add(
                        rhs.v3,
                        self.b12.mul_add(
                            rhs.v2,
                            self.v3.mul_add(
                                rhs.b31,
                                self.v2
                                    .mul_add(-rhs.b12, self.s.mul_add(rhs.v1, self.v1 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            v2: (-self.pss).mul_add(
                rhs.b31,
                (-self.b31).mul_add(
                    rhs.pss,
                    self.b23.mul_add(
                        rhs.v3,
                        (-self.b12).mul_add(
                            rhs.v1,
                            self.v3.mul_add(
                                -rhs.b23,
                                self.v2
                                    .mul_add(rhs.s, self.s.mul_add(rhs.v2, self.v1 * rhs.b12)),
                            ),
                        ),
                    ),
                ),
            ),
            v3: (-self.pss).mul_add(
                rhs.b12,
                (-self.b12).mul_add(
                    rhs.pss,
                    self.b31.mul_add(
                        rhs.v1,
                        (-self.b23).mul_add(
                            rhs.v2,
                            self.v3.mul_add(
                                rhs.s,
                                self.v2
                                    .mul_add(rhs.b23, self.s.mul_add(rhs.v3, -(self.v1 * rhs.b31))),
                            ),
                        ),
                    ),
                ),
            ),
            b12: self.b31.mul_add(
                rhs.b23,
                self.b23
                    .mul_add(-rhs.b31, self.s.mul_add(rhs.b12, self.b12 * rhs.s)),
            ),
            b23: self.b31.mul_add(
                -rhs.b12,
                self.b12
                    .mul_add(rhs.b31, self.s.mul_add(rhs.b23, self.b23 * rhs.s)),
            ),
            b31: self.b23.mul_add(
                rhs.b12,
                self.b12
                    .mul_add(-rhs.b23, self.s.mul_add(rhs.b31, self.b31 * rhs.s)),
            ),
            pss: self.s.mul_add(rhs.pss, self.pss * rhs.s),
        }
    }
}

impl std::ops::Mul for Multivector3D {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self::Output {
        Self {
            s: (-self.pss).mul_add(
                rhs.pss,
                (-self.b31).mul_add(
                    rhs.b31,
                    (-self.b23).mul_add(
                        rhs.b23,
                        self.b12.mul_add(
                            -rhs.b12,
                            self.v3.mul_add(
                                rhs.v3,
                                self.v2
                                    .mul_add(rhs.v2, self.s.mul_add(rhs.s, self.v1 * rhs.v1)),
                            ),
                        ),
                    ),
                ),
            ),
            v1: (-self.pss).mul_add(
                rhs.b23,
                (-self.b23).mul_add(
                    rhs.pss,
                    (-self.b31).mul_add(
                        rhs.v3,
                        self.b12.mul_add(
                            rhs.v2,
                            self.v3.mul_add(
                                rhs.b31,
                                self.v2
                                    .mul_add(-rhs.b12, self.s.mul_add(rhs.v1, self.v1 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            v2: (-self.pss).mul_add(
                rhs.b31,
                (-self.b31).mul_add(
                    rhs.pss,
                    self.b23.mul_add(
                        rhs.v3,
                        (-self.b12).mul_add(
                            rhs.v1,
                            self.v3.mul_add(
                                -rhs.b23,
                                self.v1
                                    .mul_add(rhs.b12, self.s.mul_add(rhs.v2, self.v2 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            v3: (-self.pss).mul_add(
                rhs.b12,
                (-self.b12).mul_add(
                    rhs.pss,
                    self.b31.mul_add(
                        rhs.v1,
                        (-self.b23).mul_add(
                            rhs.v2,
                            self.v2.mul_add(
                                rhs.b23,
                                self.v1
                                    .mul_add(-rhs.b31, self.s.mul_add(rhs.v3, self.v3 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            b12: self.b31.mul_add(
                rhs.b23,
                (-self.b23).mul_add(
                    rhs.b31,
                    self.pss.mul_add(
                        rhs.v3,
                        self.v3.mul_add(
                            rhs.pss,
                            self.v2.mul_add(
                                -rhs.v1,
                                self.v1
                                    .mul_add(rhs.v2, self.s.mul_add(rhs.b12, self.b12 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            b23: self.b12.mul_add(
                rhs.b31,
                (-self.b31).mul_add(
                    rhs.b12,
                    self.pss.mul_add(
                        rhs.v1,
                        self.v1.mul_add(
                            rhs.pss,
                            self.v3.mul_add(
                                -rhs.v2,
                                self.v2
                                    .mul_add(rhs.v3, self.s.mul_add(rhs.b23, self.b23 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            b31: self.b23.mul_add(
                rhs.b12,
                (-self.b12).mul_add(
                    rhs.b23,
                    self.pss.mul_add(
                        rhs.v2,
                        self.v2.mul_add(
                            rhs.pss,
                            self.v1.mul_add(
                                -rhs.v3,
                                self.v3
                                    .mul_add(rhs.v1, self.s.mul_add(rhs.b31, self.b31 * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
            pss: self.b31.mul_add(
                rhs.v2,
                self.b23.mul_add(
                    rhs.v1,
                    self.b12.mul_add(
                        rhs.v3,
                        self.v3.mul_add(
                            rhs.b12,
                            self.v2.mul_add(
                                rhs.b31,
                                self.v1
                                    .mul_add(rhs.b23, self.s.mul_add(rhs.pss, self.pss * rhs.s)),
                            ),
                        ),
                    ),
                ),
            ),
        }
    }
}
