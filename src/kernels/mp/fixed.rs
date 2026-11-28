//! Fixed-point working values for the elementary functions.
//!
//! A [`Fx`] is `v * 2^-s` with `v` an integer, so every value it holds is
//! exactly representable as a [`Real`]. Series are evaluated here and the
//! truncation of each division or multiplication is counted in units of
//! `2^-s`; the caller turns that count into the radius of the result
//! [`super::ball::Ball`]. Because every multiplier used by the series in
//! [`super::funcs`] has magnitude below one, an error of `k` units stays
//! below `k` units plus the number of later operations, which is what makes a
//! simple operation count a valid bound.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::One;
use num_traits::Zero;

use super::real::Real;
use super::real::Round;

/// An integer scaled by `2^-s`.
#[derive(Clone, Debug)]
pub struct Fx {
    /// The scaled integer: the value is `v * 2^-s`.
    pub(crate) v: BigInt,
    /// The number of fractional bits.
    pub(crate) s: u32,
}

impl Fx {
    /// Zero at scale `s`.
    pub(crate) fn zero(s: u32) -> Self {
        Self { v: BigInt::zero(), s }
    }

    /// One at scale `s`.
    pub(crate) fn one(s: u32) -> Self {
        Self {
            v: BigInt::one() << usize::try_from(s).unwrap_or(0),
            s,
        }
    }

    /// The integer `n` at scale `s`.
    pub(crate) fn from_i64(
        n: i64,
        s: u32,
    ) -> Self {
        Self {
            v: BigInt::from(n) << usize::try_from(s).unwrap_or(0),
            s,
        }
    }

    /// The exact value of `real` truncated toward zero at scale `s`.
    pub(crate) fn from_real_trunc(
        real: &Real,
        s: u32,
    ) -> Self {
        let shift = real.exponent() + i64::from(s);
        let v = if shift >= 0 {
            real.mantissa() << usize::try_from(shift).unwrap_or(0)
        } else {
            real.mantissa() >> usize::try_from(-shift).unwrap_or(0)
        };
        Self { v, s }
    }

    /// The exact value as a [`Real`] at precision `prec`.
    ///
    /// # Errors
    /// Fails when the exponent leaves the representable range.
    pub(crate) fn to_real(
        &self,
        prec: u32,
    ) -> Result<Real, String> {
        Real::rounded(self.v.clone(), -i64::from(self.s), prec, Round::NearestEven)
    }

    /// The exact value as a [`Real`] with enough precision to be lossless.

    /// The sum, exact.
    pub(crate) fn add(
        &self,
        other: &Self,
    ) -> Self {
        debug_assert_eq!(self.s, other.s);
        Self {
            v: &self.v + &other.v,
            s: self.s,
        }
    }

    /// The difference, exact.
    pub(crate) fn sub(
        &self,
        other: &Self,
    ) -> Self {
        debug_assert_eq!(self.s, other.s);
        Self {
            v: &self.v - &other.v,
            s: self.s,
        }
    }

    /// The product truncated toward zero; adds one unit of error.
    pub(crate) fn mul(
        &self,
        other: &Self,
    ) -> Self {
        debug_assert_eq!(self.s, other.s);
        Self {
            v: (&self.v * &other.v) >> usize::try_from(self.s).unwrap_or(0),
            s: self.s,
        }
    }

    /// The quotient truncated toward zero; adds one unit of error.
    ///
    /// Returns zero for a zero divisor, which no caller relies on: every
    /// division in [`super::funcs`] has a constant, nonzero divisor.
    pub(crate) fn div(
        &self,
        other: &Self,
    ) -> Self {
        debug_assert_eq!(self.s, other.s);
        if other.v.is_zero() {
            return Self::zero(self.s);
        }
        let numerator = &self.v << usize::try_from(self.s).unwrap_or(0);
        Self {
            v: numerator / &other.v,
            s: self.s,
        }
    }

    /// The quotient by a small positive integer, truncated toward zero.
    pub(crate) fn div_i64(
        &self,
        divisor: i64,
    ) -> Self {
        Self {
            v: &self.v / BigInt::from(divisor),
            s: self.s,
        }
    }

    /// The value multiplied by `2^k`, exact.
    pub(crate) fn scale_pow2(
        &self,
        k: i64,
    ) -> Self {
        let v = if k >= 0 {
            &self.v << usize::try_from(k).unwrap_or(0)
        } else {
            &self.v >> usize::try_from(-k).unwrap_or(0)
        };
        Self { v, s: self.s }
    }

    /// The value negated, exact.
    pub(crate) fn neg(&self) -> Self {
        Self {
            v: -self.v.clone(),
            s: self.s,
        }
    }

    /// Whether the value is zero.
    pub(crate) fn is_zero(&self) -> bool {
        self.v.is_zero()
    }

    /// `|self| < |other|`.
    pub(crate) fn abs_lt(
        &self,
        other: &Self,
    ) -> bool {
        self.v.magnitude() < other.v.magnitude()
    }

    /// The square root truncated toward zero; adds one unit of error.
    pub(crate) fn sqrt(&self) -> Self {
        debug_assert_ne!(self.v.sign(), num_bigint::Sign::Minus);
        // sqrt(v * 2^-s) = sqrt(v * 2^s) * 2^-s
        let scaled = self.v.magnitude() << usize::try_from(self.s).unwrap_or(0);
        Self {
            v: BigInt::from_biguint(num_bigint::Sign::Plus, scaled.sqrt()),
            s: self.s,
        }
    }
}

/// `2^-s` as a [`Real`], the unit of the error count.
///
/// # Errors
/// Fails when the exponent leaves the representable range.
pub fn ulp(s: u32) -> Result<Real, String> {
    Real::rounded(BigInt::from(1), -i64::from(s), 1, Round::TowardZero)
}

/// `count * 2^-s` as a [`Real`].
///
/// # Errors
/// Fails when the exponent leaves the representable range.
pub fn units(
    count: u64,
    s: u32,
) -> Result<Real, String> {
    let magnitude = BigUint::from(count);
    let value = Real::rounded(
        BigInt::from_biguint(num_bigint::Sign::Plus, magnitude),
        0,
        64,
        Round::TowardZero,
    )?;
    let unit = ulp(s)?;
    value.mul(&unit, 64, Round::AwayFromZero)
}
