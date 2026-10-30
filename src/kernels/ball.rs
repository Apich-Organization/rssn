//! # Complex balls with outward-rounded real interval arithmetic
//!
//! A [`Cball`] is a rectangle `[re.lo, re.hi] + i [im.lo, im.hi]` of the
//! complex plane. Its arithmetic is the rectangle arithmetic of
//! [`Interval`](crate::backend::tape::Interval), which rounds every
//! elementary operation outward by one ulp, so the computed rectangle
//! encloses the exact image of every point of the operands. That is what
//! makes the Krawczyk test of [`crate::kernels::homotopy_cert`] a proof
//! rather than a numerical check.
//!
//! A system is a list of [`CTerm`]s (a complex coefficient and an exponent
//! vector per term). [`real_terms`] lifts the real `PolyTerm` systems of the
//! homotopy kernel; [`eval_terms`], [`eval_jacobian_terms`] and
//! [`values_at_terms`] evaluate them over balls and points.

use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::kernels::rootfind::PolyTerm;
use num_complex::Complex64 as Cx;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

/// A complex rectangle `[re] + i [im]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cball {
    /// Real part interval.
    pub re: Interval,
    /// Imaginary part interval.
    pub im: Interval,
}

impl Cball {
    /// The point `re + i im`.
    #[must_use]
    pub const fn point(
        re: f64,
        im: f64,
    ) -> Self {
        Self {
            re: Interval::point(re),
            im: Interval::point(im),
        }
    }

    /// A [`Cx`] as a point ball.
    #[must_use]
    pub const fn of(value: Cx) -> Self {
        Self::point(value.re, value.im)
    }

    /// The square `|z - (re + i im)|_inf <= radius` (one radius for both
    /// parts), widened outward.
    #[must_use]
    pub fn square(
        re: f64,
        im: f64,
        radius: f64,
    ) -> Self {
        let (lo, hi) = ((re - radius).next_down(), (re + radius).next_up());
        let (ilo, ihi) = ((im - radius).next_down(), (im + radius).next_up());
        Self {
            re: Interval { lo, hi },
            im: Interval { lo: ilo, hi: ihi },
        }
    }

    /// A ball with separate end points.
    #[must_use]
    pub const fn rectangle(
        re: (f64, f64),
        im: (f64, f64),
    ) -> Self {
        Self {
            re: Interval { lo: re.0, hi: re.1 },
            im: Interval { lo: im.0, hi: im.1 },
        }
    }

    /// `0`.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            re: Interval::point(0.0),
            im: Interval::point(0.0),
        }
    }

    /// `1`.
    #[must_use]
    pub const fn one() -> Self {
        Self {
            re: Interval::point(1.0),
            im: Interval::point(0.0),
        }
    }

    /// The centre.
    #[must_use]
    pub fn centre(&self) -> Cx {
        Cx::new(
            0.5 * self.re.lo + 0.5 * self.re.hi,
            0.5 * self.im.lo + 0.5 * self.im.hi,
        )
    }

    /// The larger of the two half-widths.
    #[must_use]
    pub fn half_width(&self) -> f64 {
        0.5 * (self.re.hi - self.re.lo).max(self.im.hi - self.im.lo)
    }

    /// Whether both parts are finite.
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        self.re.hi.is_finite()
            && self.re.lo.is_finite()
            && self.im.hi.is_finite()
            && self.im.lo.is_finite()
    }

    /// Whether `z` lies in the rectangle.
    #[must_use]
    pub fn contains(
        &self,
        z: Cx,
    ) -> bool {
        self.re.contains(z.re) && self.im.contains(z.im)
    }

    /// A lower bound of `|z|` for `z` in the rectangle.
    #[must_use]
    pub fn modulus_lower(&self) -> f64 {
        let re = if self.re.lo <= 0.0 && self.re.hi >= 0.0 {
            0.0
        } else {
            self.re.lo.abs().min(self.re.hi.abs())
        };
        let im = if self.im.lo <= 0.0 && self.im.hi >= 0.0 {
            0.0
        } else {
            self.im.lo.abs().min(self.im.hi.abs())
        };
        re.hypot(im)
    }

    /// An upper bound of `|z|` over the rectangle.
    #[must_use]
    pub fn modulus_upper(&self) -> f64 {
        self.re
            .lo
            .abs()
            .max(self.re.hi.abs())
            .hypot(self.im.lo.abs().max(self.im.hi.abs()))
    }

    /// The intersection, when the rectangles meet.
    #[must_use]
    pub fn intersect(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let (rlo, rhi) = (self.re.lo.max(other.re.lo), self.re.hi.min(other.re.hi));
        let (ilo, ihi) = (self.im.lo.max(other.im.lo), self.im.hi.min(other.im.hi));
        (rlo <= rhi && ilo <= ihi).then_some(Self {
            re: Interval { lo: rlo, hi: rhi },
            im: Interval { lo: ilo, hi: ihi },
        })
    }

    /// Whether this rectangle lies in the interior of `outer` (strictly in
    /// both coordinates, which is the hypothesis of the Krawczyk theorem).
    #[must_use]
    pub fn inside_strictly(
        &self,
        outer: &Self,
    ) -> bool {
        outer.re.lo < self.re.lo
            && self.re.hi < outer.re.hi
            && outer.im.lo < self.im.lo
            && self.im.hi < outer.im.hi
    }

    /// `z^n` for `n >= 0` by repeated multiplication (one-ball enclosure).
    #[must_use]
    pub fn powi(
        &self,
        n: u32,
    ) -> Self {
        let mut r = Self::one();
        for _ in 0..n {
            r = r * *self;
        }
        r
    }

    /// `self / other`, `None` when the denominator may vanish.
    #[must_use]
    pub fn div(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let denominator = other.re.mul(other.re).add(other.im.mul(other.im));
        if denominator.lo <= 0.0 {
            return None;
        }
        let inverse = Self {
            re: other.re.div(denominator),
            im: other.im.neg().div(denominator),
        };
        Some(*self * inverse)
    }
}

impl Add for Cball {
    type Output = Self;

    fn add(
        self,
        other: Self,
    ) -> Self {
        Self {
            re: self.re.add(other.re),
            im: self.im.add(other.im),
        }
    }
}

impl Sub for Cball {
    type Output = Self;

    fn sub(
        self,
        other: Self,
    ) -> Self {
        Self {
            re: self.re.sub(other.re),
            im: self.im.sub(other.im),
        }
    }
}

impl Mul for Cball {
    type Output = Self;

    fn mul(
        self,
        other: Self,
    ) -> Self {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i, every product an
        // outward-rounded interval product.
        let re = self.re.mul(other.re).sub(self.im.mul(other.im));
        let im = self.re.mul(other.im).add(self.im.mul(other.re));
        Self { re, im }
    }
}

impl Neg for Cball {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }
}

/// One term of a complex polynomial: `coef * z_0^e_0 * z_1^e_1 * ...`.
#[derive(Clone, Debug, PartialEq)]
pub struct CTerm {
    /// Complex coefficient.
    pub coef: Cx,
    /// Exponents, one per variable.
    pub exps: Vec<u32>,
}

impl CTerm {
    /// A real-coefficient term.
    #[must_use]
    pub const fn real(
        coef: f64,
        exps: Vec<u32>,
    ) -> Self {
        Self {
            coef: Cx::new(coef, 0.0),
            exps,
        }
    }
}

/// Lifts a real system to complex terms.
#[must_use]
pub fn real_terms(system: &[Vec<PolyTerm>]) -> Vec<Vec<CTerm>> {
    system
        .iter()
        .map(|equation| {
            equation
                .iter()
                .map(|term| CTerm {
                    coef: Cx::new(term.coef, 0.0),
                    exps: term.exps.clone(),
                })
                .collect()
        })
        .collect()
}

/// The value of one term over the box point `x`.
fn term_value(
    term: &CTerm,
    x: &[Cball],
) -> Cball {
    let mut value = Cball::of(term.coef);
    for (&exponent, &variable) in term.exps.iter().zip(x) {
        value = value * variable.powi(exponent);
    }
    value
}

/// `f_i(X)` for every equation of the system (an enclosure of the range
/// over the box).
#[must_use]
pub fn eval_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
) -> Vec<Cball> {
    system
        .iter()
        .map(|equation| {
            equation
                .iter()
                .fold(Cball::zero(), |sum, term| sum + term_value(term, x))
        })
        .collect()
}

/// The complex Jacobian `∂f_i/∂z_j` over the box.
///
/// Each entry is an enclosure; the derivative of a monomial is
/// `c e_j z_j^{e_j - 1} Π_{k ≠ j} z_k^{e_k}` (zero when `e_j = 0`).
#[must_use]
pub fn eval_jacobian_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
) -> Vec<Vec<Cball>> {
    system
        .iter()
        .map(|equation| {
            (0..x.len())
                .map(|column| {
                    equation
                        .iter()
                        .filter(|term| term.exps.get(column).copied().unwrap_or(0) > 0)
                        .fold(Cball::zero(), |sum, term| {
                            let exponent = term.exps.get(column).copied().unwrap_or(0);
                            let mut value =
                                Cball::of(term.coef) * Cball::point(f64::from(exponent), 0.0);
                            for (k, (&power, &variable)) in term.exps.iter().zip(x).enumerate() {
                                let power = if k == column {
                                    power - 1
                                } else {
                                    power
                                };
                                value = value * variable.powi(power);
                            }
                            sum + value
                        })
                })
                .collect()
        })
        .collect()
}

/// The values `f_i(x)` at a point, in plain `f64` arithmetic.
#[must_use]
pub fn values_at_terms(
    system: &[Vec<CTerm>],
    x: &[Cx],
) -> Vec<Cx> {
    system
        .iter()
        .map(|equation| {
            equation
                .iter()
                .map(|term| {
                    let mut value = term.coef;
                    for (&power, &variable) in term.exps.iter().zip(x) {
                        value *= variable.powi(i32::try_from(power).unwrap_or(0));
                    }
                    value
                })
                .sum()
        })
        .collect()
}

/// The complex Jacobian at a point.
#[must_use]
pub fn jacobian_at_terms(
    system: &[Vec<CTerm>],
    x: &[Cx],
) -> Vec<Vec<Cx>> {
    system
        .iter()
        .map(|equation| {
            (0..x.len())
                .map(|column| {
                    equation
                        .iter()
                        .filter(|term| term.exps.get(column).copied().unwrap_or(0) > 0)
                        .map(|term| {
                            let exponent = f64::from(term.exps.get(column).copied().unwrap_or(0));
                            let mut value = term.coef * exponent;
                            for (k, (&power, &variable)) in term.exps.iter().zip(x).enumerate() {
                                let power = i32::try_from(if k == column {
                                    power - 1
                                } else {
                                    power
                                })
                                .unwrap_or(0);
                                value *= variable.powi(power);
                            }
                            value
                        })
                        .sum()
                })
                .collect()
        })
        .collect()
}

/// `f_i(X)` of a real system over a box.
#[must_use]
pub fn eval_polys(
    system: &[Vec<PolyTerm>],
    x: &[Cball],
) -> Vec<Cball> {
    eval_terms(&real_terms(system), x)
}

/// The Jacobian of a real system over a box.
#[must_use]
pub fn eval_jacobian(
    system: &[Vec<PolyTerm>],
    x: &[Cball],
) -> Vec<Vec<Cball>> {
    eval_jacobian_terms(&real_terms(system), x)
}

/// The values of a real system at a point.
#[must_use]
pub fn values_at(
    system: &[Vec<PolyTerm>],
    x: &[Cx],
) -> Vec<Cx> {
    values_at_terms(&real_terms(system), x)
}

/// The Jacobian of a real system at a point.
#[must_use]
pub fn jacobian_at(
    system: &[Vec<PolyTerm>],
    x: &[Cx],
) -> Vec<Vec<Cx>> {
    jacobian_at_terms(&real_terms(system), x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ball(
        re: (f64, f64),
        im: (f64, f64),
    ) -> Cball {
        Cball::rectangle(re, im)
    }

    #[test]
    fn arithmetic_encloses_the_pointwise_image() {
        let a = ball((1.0, 2.0), (-1.0, 0.5));
        let b = ball((-3.0, 0.5), (0.25, 2.0));
        let sum = a + b;
        let difference = a - b;
        let product = a * b;
        for &xr in &[1.0_f64, 1.5, 2.0] {
            for &xi in &[-1.0_f64, -0.5, 0.5] {
                for &yr in &[-3.0_f64, -1.0, 0.5] {
                    for &yi in &[0.25_f64, 1.0, 2.0] {
                        let (x, y) = (Cx::new(xr, xi), Cx::new(yr, yi));
                        assert!(sum.contains(x + y));
                        assert!(difference.contains(x - y));
                        assert!(product.contains(x * y));
                    }
                }
            }
        }
        assert!(sum.contains(Cx::new(2.0, -0.5)));
        assert!(!sum.contains(Cx::new(3.0, 1.0)));
    }

    #[test]
    fn division_and_powers() {
        let z = ball((1.0, 2.0), (0.5, 1.0));
        let w = ball((3.0, 4.0), (-2.0, -1.0));
        let q = z.div(&w).unwrap();
        for &zr in &[1.0_f64, 2.0] {
            for &zi in &[0.5_f64, 1.0] {
                for &wr in &[3.0_f64, 4.0] {
                    for &wi in &[-2.0_f64, -1.0] {
                        let value = Cx::new(zr, zi) / Cx::new(wr, wi);
                        assert!(q.contains(value), "{value} not in {q:?}");
                    }
                }
            }
        }
        // a denominator that may vanish gives no information
        assert!(z.div(&Cball::square(0.0, 0.0, 1.0)).is_none());
        let cube = ball((1.0, 2.0), (0.0, 1.0)).powi(3);
        for &re in &[1.0_f64, 2.0] {
            for &im in &[0.0_f64, 1.0] {
                assert!(cube.contains(Cx::new(re, im).powi(3)));
            }
        }
        assert!(z.powi(0).contains(Cx::new(1.0, 0.0)));
    }

    #[test]
    fn rectangles_and_norms() {
        let a = ball((1.0, 2.0), (3.0, 4.0));
        let b = ball((1.5, 3.0), (3.5, 9.0));
        let c = a.intersect(&b).unwrap();
        assert!((c.re.lo - 1.5).abs() < 1e-12 && (c.im.hi - 4.0).abs() < 1e-12);
        assert!(a.intersect(&ball((5.0, 6.0), (0.0, 1.0))).is_none());
        assert!(a.inside_strictly(&ball((0.0, 3.0), (2.0, 5.0))));
        assert!(!a.inside_strictly(&a));
        let d = ball((3.0, 4.0), (0.0, 0.0));
        assert!((d.modulus_lower() - 3.0).abs() < 1e-12 && (d.modulus_upper() - 4.0).abs() < 1e-12);
        let through = ball((-1.0, 1.0), (2.0, 3.0));
        assert!((through.modulus_lower() - 2.0).abs() < 1e-12);
        let s = Cball::square(0.5, -0.25, 0.1);
        assert!(s.contains(Cx::new(0.55, -0.2)));
        assert!(!s.contains(Cx::new(0.7, -0.2)));
        assert!((s.centre() - Cx::new(0.5, -0.25)).norm() < 1e-15);
        assert!((s.half_width() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn polynomial_evaluation_encloses_values() {
        // f(x, y) = x^2 y - 3 y + 1, g(x, y) = x + y^2.
        let system = vec![
            vec![
                PolyTerm {
                    coef: 1.0,
                    exps: vec![2, 1],
                },
                PolyTerm {
                    coef: -3.0,
                    exps: vec![0, 1],
                },
                PolyTerm {
                    coef: 1.0,
                    exps: vec![0, 0],
                },
            ],
            vec![
                PolyTerm {
                    coef: 1.0,
                    exps: vec![1, 0],
                },
                PolyTerm {
                    coef: 1.0,
                    exps: vec![0, 2],
                },
            ],
        ];
        let x = ball((0.5, 1.5), (-0.5, 0.5));
        let y = ball((1.0, 2.0), (0.0, 1.0));
        let values = eval_polys(&system, &[x, y]);
        let jacobian = eval_jacobian(&system, &[x, y]);
        for &xr in &[0.5_f64, 1.5] {
            for &xi in &[-0.5_f64, 0.5] {
                for &yr in &[1.0_f64, 2.0] {
                    for &yi in &[0.0_f64, 1.0] {
                        let (p, q) = (Cx::new(xr, xi), Cx::new(yr, yi));
                        assert!(values[0].contains(p * p * q - 3.0 * q + 1.0));
                        assert!(values[1].contains(p + q * q));
                        assert!(jacobian[0][0].contains(2.0 * p * q));
                        assert!(jacobian[0][1].contains(p * p - 3.0));
                        assert!(jacobian[1][0].contains(Cx::new(1.0, 0.0)));
                        assert!(jacobian[1][1].contains(2.0 * q));
                    }
                }
            }
        }
        let point = [Cx::new(0.7, -0.2), Cx::new(1.3, 0.4)];
        let values = values_at(&system, &point);
        assert!(
            (values[0] - (point[0] * point[0] * point[1] - 3.0 * point[1] + 1.0)).norm() < 1e-14
        );
        let jacobian = jacobian_at(&system, &point);
        assert!((jacobian[0][0] - 2.0 * point[0] * point[1]).norm() < 1e-14);
        assert!((jacobian[1][1] - 2.0 * point[1]).norm() < 1e-14);
    }

    #[test]
    fn complex_coefficients_and_the_real_lift_agree() {
        let terms = vec![vec![
            CTerm {
                coef: Cx::new(1.0, 1.0),
                exps: vec![2],
            },
            CTerm {
                coef: Cx::new(-2.0, 0.5),
                exps: vec![0],
            },
        ]];
        let x = [ball((0.4, 0.6), (-0.1, 0.1))];
        let values = eval_terms(&terms, &x);
        let jacobian = eval_jacobian_terms(&terms, &x);
        for &re in &[0.4_f64, 0.6] {
            for &im in &[-0.1_f64, 0.1] {
                let z = Cx::new(re, im);
                assert!(values[0].contains(Cx::new(1.0, 1.0) * z * z + Cx::new(-2.0, 0.5)));
                assert!(jacobian[0][0].contains(Cx::new(2.0, 2.0) * z));
            }
        }
        let real = vec![vec![PolyTerm {
            coef: 3.0,
            exps: vec![1],
        }]];
        let lifted = real_terms(&real);
        assert_eq!(lifted.len(), 1);
        assert_eq!(lifted[0][0].coef, Cx::new(3.0, 0.0));
        let at = values_at_terms(&lifted, &[Cx::new(2.0, 1.0)]);
        assert!((at[0] - Cx::new(6.0, 3.0)).norm() < 1e-15);
    }
}
