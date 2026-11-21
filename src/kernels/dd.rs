//! Double-double real and complex arithmetic.
//!
//! # Double-double arithmetic
//!
//! [`Dd`] is an unevaluated sum `hi + lo` of two `f64` with about 106 bits
//! (31 decimal digits) of precision, built from the error-free transforms
//! `two_sum` and `two_prod` (the latter through fused multiply-add);
//! [`Cdd`] is the corresponding complex type. They are used for the
//! extended-precision endgame of the homotopy solver
//! ([`crate::kernels::homotopy`]) where tightly clustered roots make the
//! `f64` Jacobian numerically singular.

use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

/// Double-double real number `hi + lo` with `|lo| <= ulp(hi)/2`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Dd {
    /// Leading part.
    pub hi: f64,
    /// Trailing correction.
    pub lo: f64,
}

fn two_sum(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn quick_two_sum(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

fn two_prod(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

impl Dd {
    /// One.
    pub const ONE: Self = Self { hi: 1.0, lo: 0.0 };
    /// Zero.
    pub const ZERO: Self = Self { hi: 0.0, lo: 0.0 };

    /// Exact conversion from `f64`.
    #[must_use]
    pub const fn from_f64(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }

    /// Builds a normalised value from an unnormalised pair.
    #[must_use]
    pub fn from_parts(
        hi: f64,
        lo: f64,
    ) -> Self {
        let (h, l) = two_sum(hi, lo);
        Self { hi: h, lo: l }
    }

    /// Rounds to the nearest `f64`.
    #[must_use]
    pub const fn to_f64(self) -> f64 {
        self.hi
    }

    /// Returns the stored high and low components without rounding their sum.
    #[must_use]
    pub const fn parts(self) -> (f64, f64) {
        (self.hi, self.lo)
    }

    /// Absolute value.
    #[must_use]
    pub fn abs(self) -> Self {
        if self.hi < 0.0 || (self.hi == 0.0 && self.lo < 0.0) {
            -self
        } else {
            self
        }
    }

    /// Square root with about 100 bits of precision.
    ///
    /// Negative inputs return NaN. Positive infinity returns infinity.
    #[must_use]
    pub fn sqrt(self) -> Self {
        if self.hi <= 0.0 {
            return if self.hi == 0.0 {
                Self::ZERO
            } else {
                Self::from_f64(f64::NAN)
            };
        }
        if self.hi < f64::from_bits((1023 - 900) << 52) {
            // Restore enough exponent headroom for the exact square residual,
            // including when the original leading term was subnormal.
            const LIFT: f64 = f64::from_bits((1023 + 256) << 52);
            const LIFTED_SUBNORMAL_UNIT: f64 = f64::from_bits((1023 + 256 - 1074) << 52);
            const RESTORE_ROOT: f64 = f64::from_bits((1023 - 128) << 52);
            let lifted = if self.hi.is_subnormal() && self.lo == 0.0 {
                // The positive subnormal's bits count units of 2^-1074.
                Self::from_f64(self.hi.to_bits() as f64 * LIFTED_SUBNORMAL_UNIT)
            } else {
                self.scale(LIFT)
            };
            return lifted.sqrt_residual().scale(RESTORE_ROOT);
        }
        if !self.hi.is_finite() {
            return Self::from_f64(self.hi);
        }
        self.sqrt_residual()
    }

    fn sqrt_residual(self) -> Self {
        let root = self.hi.sqrt();
        // FMA cancels the square before rounding, even near f64::MAX.
        let residual = (-root).mul_add(root, self.hi);
        if self.lo == 0.0 {
            // The hardware root already has the correctly rounded high part.
            // A rounded correction at half an ULP can otherwise move it.
            Self {
                hi: root,
                lo: (0.5 * residual) / root,
            }
        } else {
            // Both contributions are small relative to the leading input.
            let residual = residual + self.lo;
            Self::from_parts(root, (0.5 * residual) / root)
        }
    }

    fn scale(
        self,
        factor: f64,
    ) -> Self {
        Self {
            hi: self.hi * factor,
            lo: self.lo * factor,
        }
    }

    /// Euclidean norm of two double-double components.
    ///
    /// The inputs are scaled before squaring when needed. The result retains
    /// the low part unless the norm rounds to a subnormal or overflows.
    #[must_use]
    pub fn hypot(
        self,
        other: Self,
    ) -> Self {
        const LOW: f64 = f64::from_bits((1023 - 450) << 52);
        const HIGH: f64 = f64::from_bits((1023 + 450) << 52);
        const LIFT: f64 = f64::from_bits((1023 + 650) << 52);
        const LOWER: f64 = f64::from_bits((1023 - 650) << 52);

        let high = self.hi.abs().max(other.hi.abs());
        if !self.hi.is_finite() || !other.hi.is_finite() {
            return Self::from_f64(self.hi.hypot(other.hi));
        }
        if high == 0.0 {
            return Self::ZERO;
        }
        if self.hi == 0.0 {
            return other.abs();
        }
        if other.hi == 0.0 {
            return self.abs();
        }
        if (LOW..=HIGH).contains(&high) {
            return (self * self + other * other).sqrt();
        }
        let (factor, restore) = if high < LOW {
            (LIFT, LOWER)
        } else {
            (LOWER, LIFT)
        };
        let a = self.scale(factor);
        let b = other.scale(factor);
        let root = (a * a + b * b).sqrt();
        let result = root.scale(restore);
        if !result.hi.is_finite() {
            Self::from_f64(result.hi)
        } else if result.hi <= f64::MIN_POSITIVE {
            // Round to the subnormal grid using the scaled low part to break
            // ties. Scaling the two parts independently could round twice.
            let half_spacing = (f64::from_bits(1) * factor) * 0.5;
            let (rounding, tail) = two_sum(root.hi - result.hi * factor, root.lo);
            let odd = result.hi.to_bits() & 1 != 0;
            let rounded = if rounding > half_spacing
                || (rounding.to_bits() == half_spacing.to_bits()
                    && (tail > 0.0 || (tail == 0.0 && odd)))
            {
                result.hi.next_up()
            } else if rounding < -half_spacing
                || (rounding.to_bits() == (-half_spacing).to_bits()
                    && (tail < 0.0 || (tail == 0.0 && odd)))
            {
                result.hi.next_down()
            } else {
                result.hi
            };
            Self::from_f64(rounded)
        } else {
            result
        }
    }

    /// Whether both parts are finite.
    #[must_use]
    pub const fn is_finite(self) -> bool {
        self.hi.is_finite() && self.lo.is_finite()
    }
}

impl Add for Dd {
    type Output = Self;

    fn add(
        self,
        o: Self,
    ) -> Self {
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

    fn sub(
        self,
        o: Self,
    ) -> Self {
        self + (-o)
    }
}

impl Neg for Dd {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl Mul for Dd {
    type Output = Self;

    fn mul(
        self,
        o: Self,
    ) -> Self {
        let (p, e) = two_prod(self.hi, o.hi);
        let e = e + (self.hi * o.lo + self.lo * o.hi);
        let (h, l) = quick_two_sum(p, e);
        Self { hi: h, lo: l }
    }
}

impl Div for Dd {
    type Output = Self;

    fn div(
        self,
        o: Self,
    ) -> Self {
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
    fn add_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self + o;
    }
}
impl SubAssign for Dd {
    fn sub_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self - o;
    }
}
impl MulAssign for Dd {
    fn mul_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self * o;
    }
}
impl DivAssign for Dd {
    fn div_assign(
        &mut self,
        o: Self,
    ) {
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
    pub const fn new(
        re: f64,
        im: f64,
    ) -> Self {
        Self {
            re: Dd::from_f64(re),
            im: Dd::from_f64(im),
        }
    }

    /// Fast `f64` modulus of the leading components.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.re.hi.hypot(self.im.hi)
    }

    /// Modulus computed from both parts of each component.
    #[must_use]
    pub fn norm_dd(self) -> Dd {
        self.re.hypot(self.im)
    }

    /// Multiplicative inverse, using the same scale handling as division.
    ///
    /// The inverse of zero has NaN components.
    #[must_use]
    pub fn recip(self) -> Self {
        Self::new(1.0, 0.0) / self
    }
}

impl Add for Cdd {
    type Output = Self;

    fn add(
        self,
        o: Self,
    ) -> Self {
        Self {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
}
impl Sub for Cdd {
    type Output = Self;

    fn sub(
        self,
        o: Self,
    ) -> Self {
        Self {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }
}
impl Neg for Cdd {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}
impl Mul for Cdd {
    type Output = Self;

    fn mul(
        self,
        o: Self,
    ) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
}
// The leading exponent bits 0111 or 1000 cover exponents -127 through
// 128, where products and quotient corrections have ample headroom.
fn cdd_outside_direct_range(x: f64) -> u64 {
    let magnitude = x.abs().to_bits();
    let changes = (magnitude ^ (magnitude << 1)) & 0x7000_0000_0000_0000;
    let unexpected = changes ^ 0x4000_0000_0000_0000;
    unexpected & 0_u64.wrapping_sub(u64::from(magnitude != 0))
}

#[inline]
fn cdd_divide_direct(
    numerator: Cdd,
    divisor: Cdd,
) -> Cdd {
    let den = divisor.re * divisor.re + divisor.im * divisor.im;
    let num = numerator
        * Cdd {
            re: divisor.re,
            im: -divisor.im,
        };
    Cdd {
        re: num.re / den,
        im: num.im / den,
    }
}

/// A double-double significand with a separate binary exponent.
#[derive(Clone, Copy)]
struct ScaledDd {
    mantissa: Dd,
    exponent: i32,
}

const fn cdd_binary_scale(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

impl ScaledDd {
    fn new(value: Dd) -> Self {
        let bits = value.hi.to_bits();
        let magnitude = value.hi.abs().to_bits();
        if magnitude == 0 {
            return Self {
                mantissa: value,
                exponent: 0,
            };
        }
        if value.hi.is_subnormal() {
            // A normalized double-double with a subnormal high part has no
            // low part. Its bits count units of the smallest subnormal.
            let shift = 63 - magnitude.leading_zeros() as i32;
            let high = (magnitude as f64 * cdd_binary_scale(-shift)).copysign(value.hi);
            return Self {
                mantissa: Dd::from_f64(high),
                exponent: shift - 1074,
            };
        }
        let exponent = (magnitude >> 52) as i32 - 1023;
        let hi = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
        let lo = if exponent == 1023 {
            value.lo * cdd_binary_scale(-1022) * 0.5
        } else {
            value.lo * cdd_binary_scale(-exponent)
        };
        Self {
            mantissa: Dd { hi, lo },
            exponent,
        }
    }

    fn product_sum(
        a: Self,
        b: Self,
        c: Self,
        d: Self,
    ) -> Self {
        let left = a.exponent + b.exponent;
        let right = c.exponent + d.exponent;
        let left_zero = a.mantissa.hi == 0.0 || b.mantissa.hi == 0.0;
        let right_zero = c.mantissa.hi == 0.0 || d.mantissa.hi == 0.0;
        // A gap of 120 bits puts the smaller nonzero product beyond the
        // double-double precision of the sum.
        if left_zero || (!right_zero && right > left + 120) {
            return Self {
                mantissa: c.mantissa * d.mantissa,
                exponent: right,
            };
        }
        if right_zero || left > right + 120 {
            return Self {
                mantissa: a.mantissa * b.mantissa,
                exponent: left,
            };
        }
        let exponent = left.max(right);
        let lhs = a.mantissa.scale(cdd_binary_scale(left - exponent));
        let rhs = c.mantissa.scale(cdd_binary_scale(right - exponent));
        let p = lhs * b.mantissa;
        let q = rhs * d.mantissa;
        let sum = p + q;
        // In a cancellation, the low-low products and the individual
        // rounding errors can become leading terms of the result.
        let mantissa = if sum.hi.abs() < 0.5 * p.hi.abs().max(q.hi.abs()) {
            cdd_canceling_products([(lhs, b.mantissa), (rhs, d.mantissa)])
        } else {
            sum
        };
        Self { mantissa, exponent }
    }

    fn quotient(
        self,
        divisor: Self,
    ) -> Dd {
        if self.mantissa.hi == 0.0 {
            return self.mantissa;
        }
        // The numerator may be tiny after cancellation. Normalize it again
        // before the double-double division forms its residuals.
        let numerator = Self::new(self.mantissa);
        let quotient = Self::new(numerator.mantissa / divisor.mantissa);
        let exponent = self.exponent - divisor.exponent + numerator.exponent + quotient.exponent;
        let value = quotient.mantissa;
        match exponent {
            | 1024.. => Dd::from_f64(f64::INFINITY.copysign(value.hi)),
            | -1022..=1023 => value.scale(cdd_binary_scale(exponent)),
            | -1075..=-1023 => {
                // Round the total to the subnormal lattice, including the
                // tail when the scaled high part is exactly at a midpoint.
                let units = value.abs().scale(cdd_binary_scale(exponent + 1074));
                let lower = units.hi as u64;
                let midpoint = lower as f64 + 0.5;
                let residual = (units.hi - midpoint) + units.lo;
                let rounded =
                    lower + u64::from(residual > 0.0 || (residual == 0.0 && lower & 1 != 0));
                Dd::from_f64(f64::from_bits(rounded).copysign(value.hi))
            },
            | _ => Dd::from_f64(0.0_f64.copysign(value.hi)),
        }
    }
}

/// Sums the sixteen error-free product pieces in a non-overlapping expansion.
#[cold]
fn cdd_canceling_products(products: [(Dd, Dd); 2]) -> Dd {
    let mut expansion = [0.0; 16];
    let mut len = 0;
    for (a, b) in products {
        for x in [a.hi, a.lo] {
            for y in [b.hi, b.lo] {
                let (p, e) = two_prod(x, y);
                for term in [e, p] {
                    let mut high = term;
                    let mut retained = 0;
                    for read in 0..len {
                        let (sum, error) = two_sum(high, expansion[read]);
                        if error != 0.0 {
                            expansion[retained] = error;
                            retained += 1;
                        }
                        high = sum;
                    }
                    if high != 0.0 {
                        expansion[retained] = high;
                        retained += 1;
                    }
                    len = retained;
                }
            }
        }
    }
    expansion[..len]
        .split_last()
        .map_or(Dd::ZERO, |(&hi, lo)| Dd::from_parts(hi, lo.iter().sum()))
}

/// Scales each numerator independently so one component cannot flush the other.
#[cold]
fn cdd_divide_scaled(
    numerator: Cdd,
    divisor: Cdd,
) -> Cdd {
    if ![numerator.re, numerator.im, divisor.re, divisor.im]
        .into_iter()
        .all(Dd::is_finite)
        || (divisor.re.hi == 0.0 && divisor.im.hi == 0.0)
    {
        return cdd_divide_direct(numerator, divisor);
    }
    let [a, b, c, d] = [numerator.re, numerator.im, divisor.re, divisor.im].map(ScaledDd::new);
    let den = ScaledDd::product_sum(c, c, d, d);
    let minus_a = ScaledDd {
        mantissa: -a.mantissa,
        ..a
    };
    Cdd {
        re: ScaledDd::product_sum(a, c, b, d).quotient(den),
        im: ScaledDd::product_sum(b, c, minus_a, d).quotient(den),
    }
}

impl Div for Cdd {
    type Output = Self;

    fn div(
        self,
        o: Self,
    ) -> Self {
        let outside = cdd_outside_direct_range(self.re.hi)
            | cdd_outside_direct_range(self.im.hi)
            | cdd_outside_direct_range(o.re.hi)
            | cdd_outside_direct_range(o.im.hi);
        if outside == 0 {
            cdd_divide_direct(self, o)
        } else {
            cdd_divide_scaled(self, o)
        }
    }
}
impl AddAssign for Cdd {
    fn add_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self + o;
    }
}
impl SubAssign for Cdd {
    fn sub_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self - o;
    }
}
impl MulAssign for Cdd {
    fn mul_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self * o;
    }
}
impl DivAssign for Cdd {
    fn div_assign(
        &mut self,
        o: Self,
    ) {
        *self = *self / o;
    }
}
