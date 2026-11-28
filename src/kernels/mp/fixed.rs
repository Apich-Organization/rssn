//! Fixed-point working values for the elementary functions.
//!
//! A fixed-point value stores its scaled integer and fractional bit count.
//! The certified series track the error introduced by each truncation.

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use super::real::Real;

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

    /// The value of `real` rounded down at scale `s`.
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

    /// The product rounded down; loses less than one scaled unit.
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

    /// Whether the value is zero.
    pub(crate) fn is_zero(&self) -> bool {
        self.v.is_zero()
    }
}
