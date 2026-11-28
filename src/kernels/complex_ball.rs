//! Complex ball arithmetic (`ComplexBall`: center + radius) with outward rounding.
//!
//! Includes first-order forward AD (`ComplexBallJet`), complex ball linear algebra,
//! Krawczyk root certification, and certified Rouché circle verification.
//!
//! Unlike rectangular [`crate::kernels::ball::Cball`] arithmetic, a
//! [`ComplexBall`] represents a closed Euclidean disk
//! `|z - mid| <= rad`, which avoids the rotation wrapping effect under
//! repeated complex multiplications.

use std::f64::consts::FRAC_PI_2;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

use num_complex::Complex64;

use crate::backend::tape::Interval;
use crate::kernels::ball::Cball;

#[path = "complex_ball_analytic.rs"]
mod analytic;
#[path = "complex_ball_special.rs"]
mod special;

/// Conservative relative ulp slack for elementary complex operations.
const ELEM_ULP_SCALE: f64 = 16.0 * f64::EPSILON;
/// Conservative relative ulp slack for basic complex arithmetic.
const ARITH_ULP_SCALE: f64 = 6.0 * f64::EPSILON;
/// Floor slack to cover subnormal rounding.
const SUBNORMAL_SLACK: f64 = 16.0 * f64::MIN_POSITIVE;

/// Outward-rounded upper bound for `a + b` (non-negative finite or infinity).
#[inline]
fn up_add(
    a: f64,
    b: f64,
) -> f64 {
    let s = a + b;
    if s.is_nan() {
        f64::INFINITY
    } else {
        s.next_up()
    }
}

/// Outward-rounded upper bound for `a * b`.
#[inline]
fn up_mul(
    a: f64,
    b: f64,
) -> f64 {
    let p = a * b;
    if p.is_nan() {
        f64::INFINITY
    } else {
        p.next_up()
    }
}

/// Outward-rounded lower bound for non-negative `a * b`.
#[inline]
fn down_mul(
    a: f64,
    b: f64,
) -> f64 {
    let p = a * b;
    if p.is_nan() || p <= 0.0 {
        0.0
    } else {
        p.next_down().max(0.0)
    }
}

/// Lower bound on the distance of two finite centers, including subtraction rounding.
fn distance_lower(
    a: Complex64,
    b: Complex64,
) -> f64 {
    let axis = |x: f64, y: f64| (x - y).abs().next_down().max(0.0);
    norm_lower(Complex64::new(axis(a.re, b.re), axis(a.im, b.im)))
}

/// Outward-rounded upper bound for `|z|` of a point `z`.
#[inline]
fn norm_upper(z: Complex64) -> f64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return f64::INFINITY;
    }
    let a = z.re.abs();
    let b = z.im.abs();
    let far = a.max(b);
    let near = a.min(b);
    if near == 0.0 {
        return far;
    }
    let r = (near / far).next_up();
    let s = (1.0 + (r * r).next_up()).next_up();
    let root = s.sqrt().next_up();
    (far * root).next_up().max(far)
}

/// Outward-rounded lower bound for `|z|` of a point `z`.
#[inline]
fn norm_lower(z: Complex64) -> f64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return 0.0;
    }
    let a = z.re.abs();
    let b = z.im.abs();
    let far = a.max(b);
    let near = a.min(b);
    if near == 0.0 {
        return far;
    }
    let r = (near / far).next_down().max(0.0);
    let s = (1.0 + (r * r).next_down().max(0.0)).next_down();
    let root = s.sqrt().next_down();
    (far * root).next_down().max(far)
}

/// Conservative rounding error slack for a computed complex center `mid`.
#[inline]
fn center_slack(
    mid: Complex64,
    scale: f64,
) -> f64 {
    if !mid.re.is_finite() || !mid.im.is_finite() {
        return f64::INFINITY;
    }
    let mag = (mid.re.abs() + mid.im.abs()).next_up();
    up_add(up_mul(mag, scale), SUBNORMAL_SLACK)
}

/// A closed complex disk `{ z in C : |z - mid| <= rad }` with `rad >= 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexBall {
    /// Center of the disk.
    pub mid: Complex64,
    /// Non-negative radius of the disk.
    pub rad: f64,
}

impl ComplexBall {
    /// Constructs a ball `{ z : |z - mid| <= rad }`.
    /// Returns `None` if `rad < 0.0` or `rad.is_nan()` or `mid` has NaN components.
    #[must_use]
    pub fn new(
        mid: Complex64,
        rad: f64,
    ) -> Option<Self> {
        if rad.is_nan() || rad < 0.0 || mid.re.is_nan() || mid.im.is_nan() {
            return None;
        }
        Some(Self { mid, rad })
    }

    /// The exact point ball `re + i*im` with radius `0.0`.
    #[must_use]
    pub const fn point(
        re: f64,
        im: f64,
    ) -> Self {
        Self {
            mid: Complex64::new(re, im),
            rad: 0.0,
        }
    }

    /// A [`Complex64`] as a point ball with radius `0.0`.
    #[must_use]
    pub const fn of(z: Complex64) -> Self {
        Self { mid: z, rad: 0.0 }
    }

    /// The zero ball `0 + 0i`.
    #[must_use]
    pub const fn zero() -> Self {
        Self::point(0.0, 0.0)
    }

    /// The unit ball `1 + 0i`.
    #[must_use]
    pub const fn one() -> Self {
        Self::point(1.0, 0.0)
    }

    /// The imaginary unit point ball `0 + 1i`.
    #[must_use]
    pub const fn i() -> Self {
        Self::point(0.0, 1.0)
    }

    /// Converts a rectangular [`Cball`] into an enclosing [`ComplexBall`].
    #[must_use]
    pub fn from_cball(box_ball: Cball) -> Self {
        if !box_ball.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let re_mid = f64::midpoint(box_ball.re.lo, box_ball.re.hi);
        let im_mid = f64::midpoint(box_ball.im.lo, box_ball.im.hi);
        let hw_re = (box_ball.re.hi - re_mid)
            .abs()
            .max((re_mid - box_ball.re.lo).abs())
            .next_up();
        let hw_im = (box_ball.im.hi - im_mid)
            .abs()
            .max((im_mid - box_ball.im.lo).abs())
            .next_up();
        let rad = norm_upper(Complex64::new(hw_re, hw_im));
        Self {
            mid: Complex64::new(re_mid, im_mid),
            rad,
        }
    }

    /// Converts this disk into an enclosing rectangular [`Cball`].
    #[must_use]
    pub fn to_cball(&self) -> Cball {
        if !self.is_finite() {
            return Cball {
                re: Interval::ENTIRE,
                im: Interval::ENTIRE,
            };
        }
        Cball {
            re: Interval {
                lo: (self.mid.re - self.rad).next_down(),
                hi: (self.mid.re + self.rad).next_up(),
            },
            im: Interval {
                lo: (self.mid.im - self.rad).next_down(),
                hi: (self.mid.im + self.rad).next_up(),
            },
        }
    }

    /// Center of the closed disk.
    #[must_use]
    pub const fn center(&self) -> Complex64 {
        self.mid
    }

    /// Nonnegative radius of the closed disk.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.rad
    }

    /// Whether center and radius are finite.
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        self.mid.re.is_finite() && self.mid.im.is_finite() && self.rad.is_finite()
    }

    /// Whether a finite `z` may lie in this closed disk.
    ///
    /// A roundoff-sized uncertainty at the boundary counts as possible membership.
    #[must_use]
    pub fn contains(
        &self,
        z: Complex64,
    ) -> bool {
        z.is_finite() && (!self.is_finite() || distance_lower(z, self.mid) <= self.rad)
    }

    /// Whether `0` may lie in this disk.
    #[must_use]
    pub fn contains_zero(&self) -> bool {
        if !self.is_finite() {
            return true;
        }
        norm_lower(self.mid) <= self.rad
    }

    /// Lower bound on `|z|` for all `z` in this disk.
    #[must_use]
    pub fn modulus_lower(&self) -> f64 {
        if !self.is_finite() {
            return 0.0;
        }
        (norm_lower(self.mid) - self.rad).next_down().max(0.0)
    }

    /// Upper bound on `|z|` for all `z` in this disk.
    #[must_use]
    pub fn modulus_upper(&self) -> f64 {
        if !self.is_finite() {
            return f64::INFINITY;
        }
        up_add(norm_upper(self.mid), self.rad)
    }

    /// Whether two balls may intersect.
    #[must_use]
    pub fn overlaps(
        &self,
        other: &Self,
    ) -> bool {
        if !self.is_finite() || !other.is_finite() {
            return true;
        }
        distance_lower(self.mid, other.mid) <= up_add(self.rad, other.rad)
    }

    /// Whether `self` lies strictly in the interior of `outer`.
    #[must_use]
    pub fn inside_strictly(
        &self,
        outer: &Self,
    ) -> bool {
        if !self.is_finite() || !outer.is_finite() {
            return false;
        }
        let axis = |x: f64, y: f64| {
            let distance = (x - y).abs();
            if distance > 0.0 {
                distance.next_up()
            } else {
                0.0
            }
        };
        let dist_up = norm_upper(Complex64::new(
            axis(self.mid.re, outer.mid.re),
            axis(self.mid.im, outer.mid.im),
        ));
        up_add(dist_up, self.rad) < outer.rad
    }

    /// Complex conjugate of this ball.
    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            mid: self.mid.conj(),
            rad: self.rad,
        }
    }

    /// Negation `-self`.
    #[must_use]
    #[allow(
        clippy::arithmetic_side_effects,
        reason = "negation is exact for finite complex coordinates"
    )]
    pub fn neg_ball(&self) -> Self {
        Self {
            mid: -self.mid,
            rad: self.rad,
        }
    }

    /// Addition `self + other`.
    #[must_use]
    #[allow(
        clippy::arithmetic_side_effects,
        reason = "the disk radius includes a bound for midpoint rounding"
    )]
    pub fn add_ball(
        &self,
        other: &Self,
    ) -> Self {
        if !self.is_finite() || !other.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid + other.mid;
        let rad = up_add(
            up_add(self.rad, other.rad),
            center_slack(mid, ARITH_ULP_SCALE),
        );
        Self { mid, rad }
    }

    /// Subtraction `self - other`.
    #[must_use]
    pub fn sub_ball(
        &self,
        other: &Self,
    ) -> Self {
        self.add_ball(&other.neg_ball())
    }

    /// Multiplication `self * other`.
    #[must_use]
    #[allow(
        clippy::arithmetic_side_effects,
        reason = "the disk radius includes a bound for midpoint rounding"
    )]
    pub fn mul_ball(
        &self,
        other: &Self,
    ) -> Self {
        if !self.is_finite() || !other.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid * other.mid;
        let a_up = norm_upper(self.mid);
        let b_up = norm_upper(other.mid);
        let term1 = up_mul(a_up, other.rad);
        let term2 = up_mul(b_up, self.rad);
        let term3 = up_mul(self.rad, other.rad);
        // Each component forms two rounded products, which can cancel. Bound
        // their combined error using the operand norms before cancellation.
        let product_slack = up_add(up_mul(up_mul(a_up, ARITH_ULP_SCALE), b_up), SUBNORMAL_SLACK);
        let rad = up_add(up_add(up_add(term1, term2), term3), product_slack);
        Self { mid, rad }
    }

    /// Multiplication by a real scalar `s`.
    #[must_use]
    #[allow(
        clippy::arithmetic_side_effects,
        reason = "the disk radius includes a bound for midpoint rounding"
    )]
    pub fn scale_real(
        &self,
        s: f64,
    ) -> Self {
        if !self.is_finite() || !s.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid * s;
        let rad = up_add(
            up_mul(s.abs(), self.rad),
            center_slack(mid, ARITH_ULP_SCALE),
        );
        Self { mid, rad }
    }

    /// Square `self * self`.
    #[must_use]
    pub fn sqr(&self) -> Self {
        self.mul_ball(self)
    }

    /// Reciprocal `1 / self`. Returns `None` if `self` cannot be separated from zero.
    #[must_use]
    pub fn recip(&self) -> Option<Self> {
        let d_min = self.modulus_lower();
        if d_min <= 0.0 || !self.is_finite() {
            return None;
        }
        let m_norm_lo = norm_lower(self.mid);
        if m_norm_lo <= 0.0 {
            return None;
        }
        // The outward rectangular point division scales extreme binary64
        // magnitudes before squaring. Enclose its center-rounding error first.
        let mut inverse = Self::from_cball(Cball::one().div(&Cball::of(self.mid))?);
        if self.rad > 0.0 {
            // |1/z - 1/mid| <= rad / (|mid| * min |z|). Divide in two
            // stages to avoid underflow or overflow of the denominator product.
            let rad_prop = ((self.rad / d_min).next_up() / m_norm_lo).next_up();
            inverse.rad = up_add(inverse.rad, rad_prop);
        }
        Some(inverse)
    }

    /// Division `self / other`. Returns `None` if `other` cannot be separated from zero.
    #[must_use]
    pub fn div_ball(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let inv = other.recip()?;
        Some(self.mul_ball(&inv))
    }

    /// Integer power `self^n`. Returns `None` if `!self.is_finite()` or if `n < 0` and `self` may contain zero.
    #[must_use]
    pub fn powi(
        &self,
        n: i32,
    ) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        if n == 0 {
            return Some(Self::one());
        }
        let (mut base, mut exp) = if n < 0 {
            (self.recip()?, n.unsigned_abs())
        } else {
            (*self, n.unsigned_abs())
        };
        let mut acc = Self::one();
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc.mul_ball(&base);
            }
            exp >>= 1;
            if exp > 0 {
                base = base.sqr();
            }
        }
        Some(acc)
    }

    /// Complex exponential `exp(z)`.
    #[must_use]
    pub fn exp(&self) -> Self {
        if !self.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid.exp();
        let m_up = norm_upper(mid);
        let growth = self.rad.exp_m1().next_up();
        let rad = up_add(up_mul(m_up, growth), center_slack(mid, ELEM_ULP_SCALE));
        Self { mid, rad }
    }

    /// Principal complex logarithm `ln(z)`.
    /// Returns `None` if the ball intersects the non-positive real axis `(-inf, 0]`.
    #[must_use]
    pub fn ln(&self) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        let dist_to_cut = if self.mid.re <= 0.0 {
            self.mid.im.abs()
        } else {
            norm_lower(self.mid)
        };
        if dist_to_cut <= self.rad {
            return None;
        }
        let d_min = self.modulus_lower();
        if d_min <= 0.0 {
            return None;
        }
        let mid = self.mid.ln();
        let rad = up_add(
            (self.rad / d_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Principal square root `sqrt(z)`.
    /// Returns `None` if the ball crosses the negative real axis `(-inf, 0)` away from `0`.
    /// If the ball contains `0`, returns an enclosure centered at `0` covering all branches.
    #[must_use]
    pub fn sqrt(&self) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        if self.contains_zero() {
            let rad = self.modulus_upper().sqrt().next_up();
            return Some(Self {
                mid: Complex64::new(0.0, 0.0),
                rad,
            });
        }
        if self.mid.re <= 0.0 && self.mid.im.abs() <= self.rad {
            return None;
        }
        let d_min = self.modulus_lower();
        if d_min <= 0.0 {
            return None;
        }
        let mid = self.mid.sqrt();
        let denom = down_mul(2.0, d_min.sqrt().next_down()).max(f64::MIN_POSITIVE);
        let rad = up_add(
            (self.rad / denom).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Complex power `self^exponent = exp(exponent * ln(self))`.
    #[must_use]
    pub fn powc(
        &self,
        exponent: &Self,
    ) -> Option<Self> {
        let l = self.ln()?;
        Some(exponent.mul_ball(&l).exp())
    }

    /// Complex sine `sin(z)`.
    #[must_use]
    pub fn sin(&self) -> Self {
        if !self.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid.sin();
        let y_max = up_add(self.mid.im.abs(), self.rad);
        let lip = y_max.cosh().next_up();
        let rad = up_add(up_mul(lip, self.rad), center_slack(mid, ELEM_ULP_SCALE));
        Self { mid, rad }
    }

    /// Complex cosine `cos(z)`.
    #[must_use]
    pub fn cos(&self) -> Self {
        if !self.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid.cos();
        let y_max = up_add(self.mid.im.abs(), self.rad);
        let lip = y_max.cosh().next_up();
        let rad = up_add(up_mul(lip, self.rad), center_slack(mid, ELEM_ULP_SCALE));
        Self { mid, rad }
    }

    /// Complex tangent `tan(z) = sin(z) / cos(z)`.
    #[must_use]
    pub fn tan(&self) -> Option<Self> {
        let s = self.sin();
        let c = self.cos();
        s.div_ball(&c)
    }

    /// Complex hyperbolic sine `sinh(z)`.
    #[must_use]
    pub fn sinh(&self) -> Self {
        if !self.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid.sinh();
        let x_max = up_add(self.mid.re.abs(), self.rad);
        let lip = x_max.cosh().next_up();
        let rad = up_add(up_mul(lip, self.rad), center_slack(mid, ELEM_ULP_SCALE));
        Self { mid, rad }
    }

    /// Complex hyperbolic cosine `cosh(z)`.
    #[must_use]
    pub fn cosh(&self) -> Self {
        if !self.is_finite() {
            return Self {
                mid: Complex64::new(0.0, 0.0),
                rad: f64::INFINITY,
            };
        }
        let mid = self.mid.cosh();
        let x_max = up_add(self.mid.re.abs(), self.rad);
        let lip = x_max.cosh().next_up();
        let rad = up_add(up_mul(lip, self.rad), center_slack(mid, ELEM_ULP_SCALE));
        Self { mid, rad }
    }

    /// Complex hyperbolic tangent `tanh(z) = sinh(z) / cosh(z)`.
    #[must_use]
    pub fn tanh(&self) -> Option<Self> {
        let s = self.sinh();
        let c = self.cosh();
        s.div_ball(&c)
    }

    /// Distance lower bound from `mid` to the real cuts `(-inf, -1] U [1, +inf)`.
    fn avoids_real_outer_cuts(&self) -> bool {
        let dist = if self.mid.re.abs() >= 1.0 {
            self.mid.im.abs()
        } else {
            norm_lower(Complex64::new(1.0 - self.mid.re.abs(), self.mid.im))
        };
        dist > self.rad
    }

    /// Distance lower bound from `mid` to the imaginary cuts `i(-inf, -1] U i[1, +inf)`.
    fn avoids_imag_outer_cuts(&self) -> bool {
        let dist = if self.mid.im.abs() >= 1.0 {
            self.mid.re.abs()
        } else {
            norm_lower(Complex64::new(self.mid.re, 1.0 - self.mid.im.abs()))
        };
        dist > self.rad
    }

    /// Inverse sine `asin(z)`. Returns `None` if the ball meets `(-inf, -1] U [1, inf)`.
    #[must_use]
    pub fn asin(&self) -> Option<Self> {
        if !self.is_finite() || !self.avoids_real_outer_cuts() {
            return None;
        }
        let one_minus_z2 = Self::one().sub_ball(&self.sqr());
        let denom_min = one_minus_z2.modulus_lower().sqrt().next_down();
        if denom_min <= 0.0 {
            return None;
        }
        let mid = self.mid.asin();
        let rad = up_add(
            (self.rad / denom_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Inverse cosine `acos(z) = pi/2 - asin(z)`.
    #[must_use]
    pub fn acos(&self) -> Option<Self> {
        let a = self.asin()?;
        Some(Self::point(FRAC_PI_2, 0.0).sub_ball(&a))
    }

    /// Inverse tangent `atan(z)`. Returns `None` if the ball meets `+-i [1, inf)`.
    #[must_use]
    pub fn atan(&self) -> Option<Self> {
        if !self.is_finite() || !self.avoids_imag_outer_cuts() {
            return None;
        }
        let one_plus_z2 = Self::one().add_ball(&self.sqr());
        let denom_min = one_plus_z2.modulus_lower();
        if denom_min <= 0.0 {
            return None;
        }
        let mid = self.mid.atan();
        let rad = up_add(
            (self.rad / denom_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Inverse hyperbolic sine `asinh(z)`. Returns `None` if the ball meets `+-i [1, inf)`.
    #[must_use]
    pub fn asinh(&self) -> Option<Self> {
        if !self.is_finite() || !self.avoids_imag_outer_cuts() {
            return None;
        }
        let one_plus_z2 = Self::one().add_ball(&self.sqr());
        let denom_min = one_plus_z2.modulus_lower().sqrt().next_down();
        if denom_min <= 0.0 {
            return None;
        }
        let mid = self.mid.asinh();
        let rad = up_add(
            (self.rad / denom_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Inverse hyperbolic cosine `acosh(z)`. Returns `None` if the ball meets `(-inf, 1]`.
    #[must_use]
    pub fn acosh(&self) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        let dist = if self.mid.re <= 1.0 {
            self.mid.im.abs()
        } else {
            norm_lower(Complex64::new(self.mid.re - 1.0, self.mid.im))
        };
        if dist <= self.rad {
            return None;
        }
        let z2_minus_1 = self.sqr().sub_ball(&Self::one());
        let denom_min = z2_minus_1.modulus_lower().sqrt().next_down();
        if denom_min <= 0.0 {
            return None;
        }
        let mid = self.mid.acosh();
        let rad = up_add(
            (self.rad / denom_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Inverse hyperbolic tangent `atanh(z)`. Returns `None` if the ball meets `(-inf, -1] U [1, inf)`.
    #[must_use]
    pub fn atanh(&self) -> Option<Self> {
        if !self.is_finite() || !self.avoids_real_outer_cuts() {
            return None;
        }
        let one_minus_z2 = Self::one().sub_ball(&self.sqr());
        let denom_min = one_minus_z2.modulus_lower();
        if denom_min <= 0.0 {
            return None;
        }
        let mid = self.mid.atanh();
        let rad = up_add(
            (self.rad / denom_min).next_up(),
            center_slack(mid, ELEM_ULP_SCALE),
        );
        Some(Self { mid, rad })
    }

    /// Modulus ball `|z| + 0i` enclosing all `|w|` for `w in self`.
    #[must_use]
    pub fn abs_ball(&self) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        let lo = self.modulus_lower();
        let hi = self.modulus_upper();
        let mid_re = f64::midpoint(lo, hi);
        let rad = (hi - mid_re).abs().max((mid_re - lo).abs()).next_up();
        Some(Self {
            mid: Complex64::new(mid_re, 0.0),
            rad: up_add(
                rad,
                center_slack(Complex64::new(mid_re, 0.0), ARITH_ULP_SCALE),
            ),
        })
    }

    /// Principal argument `arg(z) + 0i` in `(-pi, pi]`.
    /// Returns `None` if the ball meets the non-positive real axis `(-inf, 0]`.
    #[must_use]
    pub fn arg_ball(&self) -> Option<Self> {
        let l = self.ln()?;
        Some(Self {
            mid: Complex64::new(l.mid.im, 0.0),
            rad: l.rad,
        })
    }

    /// Complex secant `sec(z) = 1 / cos(z)`.
    #[must_use]
    pub fn sec(&self) -> Option<Self> {
        analytic::sin_cos(*self).1.recip()
    }

    /// Complex cosecant `csc(z) = 1 / sin(z)`.
    #[must_use]
    pub fn csc(&self) -> Option<Self> {
        analytic::sin_cos(*self).0.recip()
    }

    /// Complex cotangent `cot(z) = cos(z) / sin(z)`.
    #[must_use]
    pub fn cot(&self) -> Option<Self> {
        let (sin, cos) = analytic::sin_cos(*self);
        cos.div_ball(&sin)
    }

    /// Complex hyperbolic secant `sech(z) = 1 / cosh(z)`.
    #[must_use]
    pub fn sech(&self) -> Option<Self> {
        analytic::sinh_cosh(*self).1.recip()
    }

    /// Complex hyperbolic cosecant `csch(z) = 1 / sinh(z)`.
    #[must_use]
    pub fn csch(&self) -> Option<Self> {
        analytic::sinh_cosh(*self).0.recip()
    }

    /// Complex hyperbolic cotangent `coth(z) = cosh(z) / sinh(z)`.
    #[must_use]
    pub fn coth(&self) -> Option<Self> {
        let (sinh, cosh) = analytic::sinh_cosh(*self);
        cosh.div_ball(&sinh)
    }

    /// Encloses the complex error function using a bounded series or an integral bound.
    #[must_use]
    pub fn erf(&self) -> Option<Self> {
        special::erf(*self)
    }

    /// Encloses the principal log-gamma function away from the nonpositive real cut.
    #[must_use]
    pub fn ln_gamma(&self) -> Option<Self> {
        special::ln_gamma(*self)
    }

    /// Encloses the complex gamma function through the principal log-gamma enclosure.
    #[must_use]
    pub fn gamma(&self) -> Option<Self> {
        let result = analytic::exp(special::ln_gamma(*self)?);
        result.is_finite().then_some(result)
    }

    /// Encloses the principal Bessel function J for nonnegative real or negative integral order.
    #[must_use]
    pub fn bessel_j(
        &self,
        nu: f64,
    ) -> Option<Self> {
        special::bessel_j(nu, *self)
    }
}

impl Add for ComplexBall {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        self.add_ball(&rhs)
    }
}

impl Sub for ComplexBall {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        self.sub_ball(&rhs)
    }
}

impl Mul for ComplexBall {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {
        self.mul_ball(&rhs)
    }
}

impl Neg for ComplexBall {
    type Output = Self;

    fn neg(self) -> Self {
        self.neg_ball()
    }
}

/// First-order forward-mode automatic differentiation over [`ComplexBall`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexBallJet {
    /// Function value enclosure.
    pub val: ComplexBall,
    /// Derivative enclosure.
    pub deriv: ComplexBall,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "forward derivatives use outward-rounded disk arithmetic"
)]
impl ComplexBallJet {
    /// Independent variable `z` with derivative `1`.
    #[must_use]
    pub const fn variable(z: ComplexBall) -> Self {
        Self {
            val: z,
            deriv: ComplexBall::one(),
        }
    }

    /// Constant `c` with derivative `0`.
    #[must_use]
    pub const fn constant(c: ComplexBall) -> Self {
        Self {
            val: c,
            deriv: ComplexBall::zero(),
        }
    }

    /// Addition.
    #[must_use]
    pub fn add(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            val: self.val + other.val,
            deriv: self.deriv + other.deriv,
        }
    }

    /// Subtraction.
    #[must_use]
    pub fn sub(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            val: self.val - other.val,
            deriv: self.deriv - other.deriv,
        }
    }

    /// Multiplication.
    #[must_use]
    pub fn mul(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            val: self.val * other.val,
            deriv: self.deriv * other.val + self.val * other.deriv,
        }
    }

    /// Division.
    #[must_use]
    pub fn div(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let val = self.val.div_ball(&other.val)?;
        let num = self.deriv * other.val - self.val * other.deriv;
        let deriv = num.div_ball(&other.val.sqr())?;
        Some(Self { val, deriv })
    }

    /// Integer power.
    #[must_use]
    pub fn powi(
        &self,
        n: i32,
    ) -> Option<Self> {
        if n == 0 {
            return Some(Self::constant(ComplexBall::one()));
        }
        let val = self.val.powi(n)?;
        let prev = self.val.powi(n.checked_sub(1)?)?;
        let deriv = prev.scale_real(f64::from(n)) * self.deriv;
        Some(Self { val, deriv })
    }

    /// Exponential.
    #[must_use]
    pub fn exp(&self) -> Self {
        let val = analytic::exp(self.val);
        Self {
            val,
            deriv: val * self.deriv,
        }
    }

    /// Principal logarithm.
    #[must_use]
    pub fn ln(&self) -> Option<Self> {
        let val = analytic::ln(self.val)?;
        let deriv = self.deriv.div_ball(&self.val)?;
        Some(Self { val, deriv })
    }

    /// Principal square root.
    #[must_use]
    pub fn sqrt(&self) -> Option<Self> {
        let val = analytic::sqrt(self.val)?;
        let two_val = val.scale_real(2.0);
        let deriv = self.deriv.div_ball(&two_val)?;
        Some(Self { val, deriv })
    }

    /// Sine.
    #[must_use]
    pub fn sin(&self) -> Self {
        let (sin, cos) = analytic::sin_cos(self.val);
        Self {
            val: sin,
            deriv: cos * self.deriv,
        }
    }

    /// Cosine.
    #[must_use]
    pub fn cos(&self) -> Self {
        let (sin, cos) = analytic::sin_cos(self.val);
        Self {
            val: cos,
            deriv: -sin * self.deriv,
        }
    }

    /// Hyperbolic sine.
    #[must_use]
    pub fn sinh(&self) -> Self {
        let (sinh, cosh) = analytic::sinh_cosh(self.val);
        Self {
            val: sinh,
            deriv: cosh * self.deriv,
        }
    }

    /// Hyperbolic cosine.
    #[must_use]
    pub fn cosh(&self) -> Self {
        let (sinh, cosh) = analytic::sinh_cosh(self.val);
        Self {
            val: cosh,
            deriv: sinh * self.deriv,
        }
    }

    /// Tangent.
    #[must_use]
    pub fn tan(&self) -> Option<Self> {
        let (numerator, denominator) = analytic::sin_cos(self.val);
        let inverse = denominator.recip()?;
        Some(Self {
            val: numerator * inverse,
            deriv: inverse.sqr() * self.deriv,
        })
    }

    /// Hyperbolic tangent.
    #[must_use]
    pub fn tanh(&self) -> Option<Self> {
        let (numerator, denominator) = analytic::sinh_cosh(self.val);
        let inverse = denominator.recip()?;
        Some(Self {
            val: numerator * inverse,
            deriv: inverse.sqr() * self.deriv,
        })
    }

    /// Complex power `self^exponent` for a constant exponent `c` or general jet.
    #[must_use]
    pub fn powc(
        &self,
        exponent: &Self,
    ) -> Option<Self> {
        let l = self.ln()?;
        let prod = exponent.mul(&l);
        Some(prod.exp())
    }

    /// Inverse sine.
    #[must_use]
    pub fn asin(&self) -> Option<Self> {
        let (val, denom) = analytic::asin(self.val)?;
        let deriv = self.deriv.div_ball(&denom)?;
        Some(Self { val, deriv })
    }

    /// Inverse cosine.
    #[must_use]
    pub fn acos(&self) -> Option<Self> {
        let (val, denom) = analytic::acos(self.val)?;
        let deriv = (-self.deriv).div_ball(&denom)?;
        Some(Self { val, deriv })
    }

    /// Inverse tangent.
    #[must_use]
    pub fn atan(&self) -> Option<Self> {
        let (val, denom) = analytic::atan(self.val)?;
        let deriv = self.deriv.div_ball(&denom)?;
        Some(Self { val, deriv })
    }

    /// Inverse hyperbolic sine.
    #[must_use]
    pub fn asinh(&self) -> Option<Self> {
        let (val, denom) = analytic::asinh(self.val)?;
        let deriv = self.deriv.div_ball(&denom)?;
        Some(Self { val, deriv })
    }

    /// Inverse hyperbolic cosine.
    #[must_use]
    pub fn acosh(&self) -> Option<Self> {
        let (val, denom) = analytic::acosh(self.val)?;
        let deriv = self.deriv.div_ball(&denom)?;
        Some(Self { val, deriv })
    }

    /// Inverse hyperbolic tangent.
    #[must_use]
    pub fn atanh(&self) -> Option<Self> {
        let (val, denom) = analytic::atanh(self.val)?;
        let deriv = self.deriv.div_ball(&denom)?;
        Some(Self { val, deriv })
    }
}

/// Matrix-vector product over [`ComplexBall`].
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "matrix products use outward-rounded disk arithmetic"
)]
pub fn mat_vec_ball(
    a: &[Vec<ComplexBall>],
    x: &[ComplexBall],
) -> Option<Vec<ComplexBall>> {
    let n = x.len();
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        if row.len() != n {
            return None;
        }
        let mut acc = ComplexBall::zero();
        for (entry, xj) in row.iter().zip(x.iter()) {
            acc = acc + (*entry * *xj);
        }
        out.push(acc);
    }
    Some(out)
}

/// Matrix-matrix product over [`ComplexBall`].
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "matrix products use outward-rounded disk arithmetic"
)]
pub fn mat_mul_ball(
    a: &[Vec<ComplexBall>],
    b: &[Vec<ComplexBall>],
) -> Option<Vec<Vec<ComplexBall>>> {
    let m = a.len();
    let k = b.len();
    let n = b.first().map_or(0, Vec::len);
    for row in a {
        if row.len() != k {
            return None;
        }
    }
    for row in b {
        if row.len() != n {
            return None;
        }
    }
    let mut out = vec![vec![ComplexBall::zero(); n]; m];
    for (i, a_row) in a.iter().enumerate() {
        let out_row = out.get_mut(i)?;
        for (j, cell) in out_row.iter_mut().enumerate() {
            let mut acc = ComplexBall::zero();
            for (l, a_il) in a_row.iter().enumerate() {
                let b_lj = b.get(l)?.get(j)?;
                acc = acc + (*a_il * *b_lj);
            }
            *cell = acc;
        }
    }
    Some(out)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "this approximate inverse is subsequently enclosed and checked by disk arithmetic"
)]
fn midpoint_inverse(a: &[Vec<ComplexBall>]) -> Option<Vec<Vec<Complex64>>> {
    let n = a.len();
    if n == 0 || a.iter().any(|row| row.len() != n) {
        return None;
    }
    // Factor the midpoint matrix once.
    let mut lu: Vec<Vec<Complex64>> = a
        .iter()
        .map(|row| row.iter().map(|c| c.mid).collect())
        .collect();
    let mut inv_mid = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    for (i, row) in inv_mid.iter_mut().enumerate() {
        if let Some(cell) = row.get_mut(i) {
            *cell = Complex64::new(1.0, 0.0);
        }
    }
    // Gauss-Jordan with partial pivoting on midpoint matrix to obtain approximate inverse R.
    for col in 0..n {
        let mut pivot_row = col;
        let mut best_norm = 0.0_f64;
        for (r, row) in lu.iter().enumerate().skip(col) {
            let val = row.get(col)?.norm();
            if val > best_norm {
                best_norm = val;
                pivot_row = r;
            }
        }
        if best_norm == 0.0 || !best_norm.is_finite() {
            return None;
        }
        if pivot_row != col {
            lu.swap(col, pivot_row);
            inv_mid.swap(col, pivot_row);
        }
        let inv_piv = Complex64::new(1.0, 0.0) / *lu.get(col)?.get(col)?;
        for j in 0..n {
            if let Some(v) = lu.get_mut(col).and_then(|r| r.get_mut(j)) {
                *v *= inv_piv;
            }
            if let Some(v) = inv_mid.get_mut(col).and_then(|r| r.get_mut(j)) {
                *v *= inv_piv;
            }
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = *lu.get(r)?.get(col)?;
            if factor == Complex64::new(0.0, 0.0) {
                continue;
            }
            for j in 0..n {
                let piv_lu = *lu.get(col)?.get(j)?;
                if let Some(v) = lu.get_mut(r).and_then(|row| row.get_mut(j)) {
                    *v -= factor * piv_lu;
                }
                let piv_inv = *inv_mid.get(col)?.get(j)?;
                if let Some(v) = inv_mid.get_mut(r).and_then(|row| row.get_mut(j)) {
                    *v -= factor * piv_inv;
                }
            }
        }
    }

    Some(inv_mid)
}

/// Solves `A x = b` over [`ComplexBall`] with verified enclosure via midpoint
/// Gaussian elimination and Krawczyk/Rump residual certification.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "matrix contractions use outward-rounded disk arithmetic and finite denominators"
)]
pub fn solve_linear_ball(
    a: &[Vec<ComplexBall>],
    b: &[ComplexBall],
) -> Option<Vec<ComplexBall>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|r| r.len() != n) {
        return None;
    }
    let inv_mid = midpoint_inverse(a)?;

    // Approximate solution x_0 = R * b.mid.
    let mut x0 = Vec::with_capacity(n);
    for row in &inv_mid {
        let mut s = Complex64::new(0.0, 0.0);
        for (rij, bj) in row.iter().zip(b.iter()) {
            s += *rij * bj.mid;
        }
        x0.push(ComplexBall::of(s));
    }

    // Residual enclosure: r = b - A * x0, and C = I - R * A.
    let ax0 = mat_vec_ball(a, &x0)?;
    let r_vec: Vec<ComplexBall> = b
        .iter()
        .zip(ax0.iter())
        .map(|(bi, axi)| *bi - *axi)
        .collect();
    let r_ball: Vec<Vec<ComplexBall>> = inv_mid
        .iter()
        .map(|row| row.iter().map(|c| ComplexBall::of(*c)).collect())
        .collect();
    let z_vec = mat_vec_ball(&r_ball, &r_vec)?;
    let ra = mat_mul_ball(&r_ball, a)?;
    let mut norm_c = 0.0_f64;
    for (i, row) in ra.iter().enumerate() {
        let mut row_sum = 0.0_f64;
        for (j, cell) in row.iter().enumerate() {
            let entry = if i == j {
                ComplexBall::one() - *cell
            } else {
                -*cell
            };
            row_sum = up_add(row_sum, entry.modulus_upper());
        }
        norm_c = norm_c.max(row_sum);
    }
    if norm_c >= 1.0 || !norm_c.is_finite() {
        return None;
    }
    let z_max = z_vec
        .iter()
        .map(ComplexBall::modulus_upper)
        .fold(0.0_f64, f64::max);
    let denom = (1.0 - norm_c).next_down();
    if denom <= 0.0 {
        return None;
    }
    let delta = (z_max / denom).next_up();
    let mut out = Vec::with_capacity(n);
    for xi in &x0 {
        out.push(ComplexBall {
            mid: xi.mid,
            rad: up_add(xi.rad, delta),
        });
    }
    Some(out)
}

/// Certifies existence and uniqueness of a simple complex root of `f` inside `ball`.
///
/// Uses the 1D complex Krawczyk operator:
/// `K(B) = m - f(m)/f'(m) + (1 - f'(B)/f'(m)) * (B - m)`.
/// Returns `Some(K(B))` if `K(B)` lies strictly inside `ball`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "Krawczyk containment is evaluated with outward-rounded disk arithmetic"
)]
pub fn certify_root_1d<F>(
    f_jet: F,
    ball: ComplexBall,
) -> Option<ComplexBall>
where
    F: Fn(ComplexBallJet) -> Option<ComplexBallJet>,
{
    if crate::cancellation::is_cancelled() || !ball.is_finite() || ball.rad <= 0.0 {
        return None;
    }
    let m = ComplexBall::of(ball.mid);
    let jm = f_jet(ComplexBallJet::variable(m))?;
    let jb = f_jet(ComplexBallJet::variable(ball))?;
    let step = jm.val.div_ball(&jm.deriv)?;
    let ratio = jb.deriv.div_ball(&jm.deriv)?;
    let contraction = ComplexBall::one() - ratio;
    let delta = ComplexBall {
        mid: Complex64::new(0.0, 0.0),
        rad: ball.rad,
    };
    let k = (m - step) + contraction * delta;
    if k.inside_strictly(&ball) {
        Some(k)
    } else {
        None
    }
}

/// Multi-dimensional complex Krawczyk root certification for `F(z) = 0`.
///
/// Given a box `balls` (`[B_0, ..., B_{n-1}]`), midpoint evaluation `f_mid = F(m)`,
/// and Jacobian enclosure `jac_box = J_F(B)` (evaluated over `balls`), computes
/// `K(B) = m - R * f_mid + (I - R * J_F(B)) * (B - m)` where `R ~ J_F(m)^{-1}`.
/// Returns `Some(K(B))` if every component `K(B)_i` lies strictly inside `B_i`.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "Krawczyk containment is evaluated with outward-rounded disk arithmetic"
)]
pub fn certify_root_nd<F, J>(
    f_mid: F,
    jac_box: J,
    balls: &[ComplexBall],
) -> Option<Vec<ComplexBall>>
where
    F: Fn(&[ComplexBall]) -> Option<Vec<ComplexBall>>,
    J: Fn(&[ComplexBall]) -> Option<Vec<Vec<ComplexBall>>>,
{
    let n = balls.len();
    if crate::cancellation::is_cancelled()
        || n == 0
        || balls.iter().any(|b| !b.is_finite() || b.rad <= 0.0)
    {
        return None;
    }
    let m: Vec<ComplexBall> = balls.iter().map(|b| ComplexBall::of(b.mid)).collect();
    let fm = f_mid(&m)?;
    if fm.len() != n {
        return None;
    }
    let jm = jac_box(&m)?;
    let jb = jac_box(balls)?;
    if jm.len() != n || jb.len() != n {
        return None;
    }
    // An approximate inverse is enough. Strict Krawczyk inclusion
    // bounds the contraction in the box's weighted maximum norm.
    // Factor once, and use each computed binary64 entry as an exact point.
    let r_mat: Vec<Vec<ComplexBall>> = midpoint_inverse(&jm)?
        .into_iter()
        .map(|row| row.into_iter().map(ComplexBall::of).collect())
        .collect();
    let step = mat_vec_ball(&r_mat, &fm)?;
    let rj = mat_mul_ball(&r_mat, &jb)?;
    let delta: Vec<ComplexBall> = balls
        .iter()
        .map(|b| {
            ComplexBall {
                mid: Complex64::new(0.0, 0.0),
                rad: b.rad,
            }
        })
        .collect();
    let mut c_mat = vec![vec![ComplexBall::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let entry = *rj.get(i)?.get(j)?;
            *c_mat.get_mut(i)?.get_mut(j)? = if i == j {
                ComplexBall::one() - entry
            } else {
                -entry
            };
        }
    }
    let c_delta = mat_vec_ball(&c_mat, &delta)?;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let ki = (*m.get(i)? - *step.get(i)?) + *c_delta.get(i)?;
        if !ki.inside_strictly(balls.get(i)?) {
            return None;
        }
        out.push(ki);
    }
    Some(out)
}

/// Certifies Rouché's theorem on the circle `|z - center| = radius` by covering
/// the circle with `segments` overlapping [`ComplexBall`]s and checking
/// `|g(B)|_upper < |f(B)|_lower` on every segment.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "circle covering bounds use outward-rounded disk arithmetic"
)]
pub fn certify_rouche_circle<F, G>(
    f: F,
    g: G,
    center: Complex64,
    radius: f64,
    segments: usize,
) -> bool
where
    F: Fn(ComplexBall) -> Option<ComplexBall>,
    G: Fn(ComplexBall) -> Option<ComplexBall>,
{
    if !center.re.is_finite()
        || !center.im.is_finite()
        || segments < 4
        || !radius.is_finite()
        || radius <= 0.0
    {
        return false;
    }
    // The midpoint of each arc differs from its endpoints by at most
    // R*pi/segments (arc length). Use one rigorous initial unit direction
    // and advance with its square instead of trusting approximate libm
    // centres independently for every arc.
    let pi = std::f64::consts::PI;
    let pi_disk = ComplexBall {
        mid: Complex64::new(0.0, pi),
        rad: pi.next_up() - pi,
    };
    let count = segments as f64;
    let n = ComplexBall {
        mid: Complex64::new(count, 0.0),
        rad: count.next_up() - count,
    };
    let Some(half_angle) = pi_disk.div_ball(&n) else {
        return false;
    };
    let mut phase = analytic::exp(half_angle);
    let advance = phase.sqr();
    let span = up_mul(radius, half_angle.modulus_upper());
    for _ in 0..segments {
        if crate::cancellation::is_cancelled() {
            return false;
        }
        let mut ball = ComplexBall::of(center) + phase.scale_real(radius);
        ball.rad = up_add(ball.rad, span);
        if !ball.is_finite() {
            return false;
        }
        phase = phase * advance;
        let Some(fb) = f(ball) else {
            return false;
        };
        let Some(gb) = g(ball) else {
            return false;
        };
        if !fb.is_finite() || !gb.is_finite() {
            return false;
        }
        if gb.modulus_upper() >= fb.modulus_lower() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ball_rotations_stay_tight() {
        let angle = std::f64::consts::TAU / 100.0;
        let rot = ComplexBall::of(Complex64::from_polar(1.0, angle));
        let mut z = ComplexBall::new(Complex64::new(1.0, 0.0), 1e-12)
            .unwrap_or_else(|| panic!("valid ball"));
        for _ in 0..100 {
            z = z * rot;
        }
        assert!(z.contains(Complex64::new(1.0, 0.0)));
        assert!(z.rad < 1e-10, "radius grew too much: {}", z.rad);
    }

    #[test]
    fn krawczyk_certifies_cube_root_of_unity() {
        let guess = Complex64::new(-0.5, 0.866_025_4);
        let ball = ComplexBall::new(guess, 1e-3).unwrap_or_else(|| panic!("valid ball"));
        let cert = certify_root_1d(
            |z| {
                let z3 = z.powi(3)?;
                Some(z3.sub(&ComplexBallJet::constant(ComplexBall::one())))
            },
            ball,
        );
        assert!(cert.is_some());
        if let Some(b) = cert {
            let exact = Complex64::new(-0.5, 0.75_f64.sqrt());
            assert!(b.contains(exact));
        }
    }

    #[test]
    fn rouche_certifies_polynomial_on_unit_circle() {
        let ok = certify_rouche_circle(
            |z| Some(z.scale_real(3.0)),
            |z| Some(z.powi(5)? + ComplexBall::one()),
            Complex64::new(0.0, 0.0),
            1.0,
            64,
        );
        assert!(ok);
    }
}
