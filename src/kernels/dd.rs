//! # Double-double arithmetic
//!
//! [`Dd`] is an unevaluated sum `hi + lo` of two `f64` with about 106 bits
//! (31 decimal digits) of precision, built from the error-free transforms
//! `two_sum` and `two_prod` (the latter through fused multiply-add);
//! [`Cdd`] is the corresponding complex type. They are used for the
//! extended-precision endgame of the homotopy solver
//! ([`crate::kernels::homotopy`]) where tightly clustered roots make the
//! `f64` Jacobian numerically singular.

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Double-double real number `hi + lo` with `|lo| <= ulp(hi)/2`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Dd {
    /// Leading part.
    pub hi: f64,
    /// Trailing correction.
    pub lo: f64,
}

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

impl Dd {
    /// Zero.
    pub const ZERO: Self = Self { hi: 0.0, lo: 0.0 };
    /// One.
    pub const ONE: Self = Self { hi: 1.0, lo: 0.0 };

    /// Exact conversion from `f64`.
    #[must_use]
    pub const fn from_f64(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }

    /// Builds a normalised value from an unnormalised pair.
    #[must_use]
    pub fn from_parts(hi: f64, lo: f64) -> Self {
        let (h, l) = two_sum(hi, lo);
        Self { hi: h, lo: l }
    }

    /// Rounds to the nearest `f64`.
    #[must_use]
    pub const fn to_f64(self) -> f64 {
        self.hi
    }

    /// Absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        if self.hi < 0.0 || (self.hi == 0.0 && self.lo < 0.0) { -self } else { self }
    }

    /// Square root (about 100 bits); negative input gives NaN.
    #[must_use]
    pub fn sqrt(self) -> Self {
        if self.hi <= 0.0 {
            return if self.hi == 0.0 { Self::ZERO } else { Self::from_f64(f64::NAN) };
        }
        let x = 1.0 / self.hi.sqrt();
        let ax = self.hi * x;
        let diff = self - Self::from_f64(ax) * Self::from_f64(ax);
        let corr = diff.hi * (x * 0.5);
        Self::from_parts(ax, corr)
    }

    /// Whether both parts are finite.
    #[must_use]
    pub const fn is_finite(self) -> bool {
        self.hi.is_finite() && self.lo.is_finite()
    }
}

impl Add for Dd {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let (s, e) = two_sum(self.hi, o.hi);
        let (t, f) = two_sum(self.lo, o.lo);
        let e = e + t;
        let (s, e) = quick_two_sum(s, e);
        let e = e + f;
        let (h, l) = quick_two_sum(s, e);
        Self { hi: h, lo: l }
    }
}

impl Sub for Dd {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + (-o)
    }
}

impl Neg for Dd {
    type Output = Self;
    fn neg(self) -> Self {
        Self { hi: -self.hi, lo: -self.lo }
    }
}

impl Mul for Dd {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let (p, e) = two_prod(self.hi, o.hi);
        let e = e + (self.hi * o.lo + self.lo * o.hi);
        let (h, l) = quick_two_sum(p, e);
        Self { hi: h, lo: l }
    }
}

impl Div for Dd {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let q1 = self.hi / o.hi;
        let r = self - o * Self::from_f64(q1);
        let q2 = r.hi / o.hi;
        let r = r - o * Self::from_f64(q2);
        let q3 = r.hi / o.hi;
        let (h, l) = quick_two_sum(q1, q2);
        Self { hi: h, lo: l } + Self::from_f64(q3)
    }
}

impl AddAssign for Dd {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}
impl SubAssign for Dd {
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}
impl MulAssign for Dd {
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}
impl DivAssign for Dd {
    fn div_assign(&mut self, o: Self) {
        *self = *self / o;
    }
}

/// Complex number with double-double components.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Cdd {
    /// Real part.
    pub re: Dd,
    /// Imaginary part.
    pub im: Dd,
}

impl Cdd {
    /// Builds a value from two `f64` components.
    #[must_use]
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re: Dd::from_f64(re), im: Dd::from_f64(im) }
    }

    /// Modulus rounded to `f64`.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.re.hi.hypot(self.im.hi)
    }
}

impl Add for Cdd {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self { re: self.re + o.re, im: self.im + o.im }
    }
}
impl Sub for Cdd {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self { re: self.re - o.re, im: self.im - o.im }
    }
}
impl Neg for Cdd {
    type Output = Self;
    fn neg(self) -> Self {
        Self { re: -self.re, im: -self.im }
    }
}
impl Mul for Cdd {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}
impl Div for Cdd {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let den = o.re * o.re + o.im * o.im;
        let num = self * Self { re: o.re, im: -o.im };
        Self { re: num.re / den, im: num.im / den }
    }
}
impl AddAssign for Cdd {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}
impl SubAssign for Cdd {
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}
impl MulAssign for Cdd {
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}
impl DivAssign for Cdd {
    fn div_assign(&mut self, o: Self) {
        *self = *self / o;
    }
}
