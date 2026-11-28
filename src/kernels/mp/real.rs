//! Binary floating point with an arbitrary mantissa precision.
//!
//! A [`Real`] is `mant * 2^exp` with `mant` held in a [`BigInt`] normalised to
//! exactly `prec` bits, so a value carries the precision it was computed at
//! instead of relying on a global setting. Every operation takes the target
//! precision and a [`Round`]ing mode, computes the result **exactly** and
//! rounds once, which is what makes the rounding correct rather than
//! approximately correct: there is no double rounding anywhere in the tower.
//!
//! The exponent is bounded by [`MAX_EXP`] and the precision by [`MAX_BITS`].
//! Both are far beyond what a decimal-digit interface exposes ([`MAX_BITS`] is
//! about 9,800 decimal digits), and they keep every exact intermediate
//! bounded: an addition whose operands are that far apart shifts at most
//! `2 * MAX_EXP` bits, so no operation can be made to allocate unbounded
//! memory by a hostile argument.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_bigint::Sign;
use num_integer::Integer;
use num_traits::Zero;

/// Largest mantissa width, in bits (about 9,800 decimal digits).
pub const MAX_BITS: u32 = 1 << 15;
/// Largest magnitude of the binary exponent.
pub const MAX_EXP: i64 = 1 << 18;

/// How a value that is not representable at the target precision is rounded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Round {
    /// Nearest, ties to an even last bit: the IEEE-754 default.
    NearestEven,
    /// Toward zero (truncation).
    TowardZero,
    /// Away from zero.
    AwayFromZero,
    /// Toward `+inf`.
    TowardPositive,
    /// Toward `-inf`.
    TowardNegative,
}

/// A binary floating-point number `mant * 2^exp` at a fixed precision.
#[derive(Clone, Debug)]
pub struct Real {
    /// Signed mantissa, normalised to exactly `prec` bits when nonzero.
    mant: BigInt,
    /// Binary exponent: the value is `mant * 2^exp`.
    exp: i64,
    /// Mantissa width in bits.
    prec: u32,
}

impl Real {
    /// The value zero at precision `prec`.
    #[must_use]
    pub fn zero(prec: u32) -> Self {
        Self {
            mant: BigInt::from(0),
            exp: 0,
            prec: prec.clamp(1, MAX_BITS),
        }
    }

    /// The value one at precision `prec`.
    #[must_use]
    pub fn one(prec: u32) -> Self {
        Self::from_i64(1, prec).unwrap_or_else(|_| Self::zero(prec))
    }

    /// The integer `value` at precision `prec`, exactly when it fits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn from_i64(
        value: i64,
        prec: u32,
    ) -> Result<Self, String> {
        Self::from_bigint(BigInt::from(value), prec, Round::NearestEven)
    }

    /// The integer `value` at precision `prec`, rounded when it does not fit.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn from_bigint(
        value: BigInt,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        Self::rounded(value, 0, prec, mode)
    }

    /// `mant * 2^exp` rounded to `prec` bits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn rounded(
        mant: BigInt,
        exp: i64,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        Self::rounded_with_inexact(mant, exp, prec, mode).map(|(value, _)| value)
    }

    fn rounded_with_inexact(
        mant: BigInt,
        exp: i64,
        prec: u32,
        mode: Round,
    ) -> Result<(Self, bool), String> {
        let prec = prec.clamp(1, MAX_BITS);
        if mant.is_zero() {
            return Ok((Self { mant, exp: 0, prec }, false));
        }
        let negative = mant.sign() == Sign::Minus;
        let mut magnitude = mant.magnitude().clone();
        let mut exponent = exp;
        let bits = magnitude.bits();
        let width = u64::from(prec);
        let mut inexact = false;
        if bits > width {
            let shift = bits - width;
            let kept = &magnitude >> shift;
            let rest = &magnitude - (&kept << shift);
            inexact = !rest.is_zero();
            let half = BigUint::from(1_u32) << usize::try_from(shift - 1).unwrap_or(0);
            let up = match mode {
                | Round::NearestEven => rest > half || (rest == half && kept.bit(0)),
                | Round::TowardZero => false,
                | Round::AwayFromZero => inexact,
                | Round::TowardPositive => !negative && inexact,
                | Round::TowardNegative => negative && inexact,
            };
            magnitude = if up {
                kept + BigUint::from(1_u32)
            } else {
                kept
            };
            exponent = exponent
                .checked_add(i64::try_from(shift).unwrap_or(i64::MAX))
                .ok_or_else(|| "exponent out of range".to_string())?;
            if magnitude.bits() > width {
                magnitude >>= 1_u32;
                exponent = exponent
                    .checked_add(1)
                    .ok_or_else(|| "exponent out of range".to_string())?;
            }
        } else if bits < width {
            let left = width - bits;
            magnitude <<= left;
            exponent = exponent
                .checked_sub(i64::try_from(left).unwrap_or(i64::MAX))
                .ok_or_else(|| "exponent out of range".to_string())?;
        }
        if exponent.unsigned_abs() > MAX_EXP.unsigned_abs() {
            return Err("exponent out of range".to_string());
        }
        let value = Self {
            mant: BigInt::from_biguint(
                if negative {
                    Sign::Minus
                } else {
                    Sign::Plus
                },
                magnitude,
            ),
            exp: exponent,
            prec,
        };
        Ok((value, inexact))
    }

    /// The exact value of an `f64` at precision `prec`.
    ///
    /// # Errors
    /// Fails for an infinity or a NaN, or when the exponent leaves the range.
    pub fn from_f64(
        value: f64,
        prec: u32,
    ) -> Result<Self, String> {
        if !value.is_finite() {
            return Err("a non-finite number has no exact representation".to_string());
        }
        if value == 0.0 {
            return Ok(Self::zero(prec));
        }
        let bits = value.to_bits();
        let negative = bits >> 63 == 1;
        let exponent_bits = i64::try_from((bits >> 52) & 0x7ff).unwrap_or(0);
        let fraction = bits & ((1_u64 << 52) - 1);
        let (magnitude, exponent) = if exponent_bits == 0 {
            (BigInt::from(fraction), -1074)
        } else {
            (BigInt::from(fraction | (1_u64 << 52)), exponent_bits - 1075)
        };
        let signed = if negative {
            -magnitude
        } else {
            magnitude
        };
        Self::rounded(signed, exponent, prec, Round::NearestEven)
    }

    /// The mantissa, normalised to [`Self::precision`] bits.
    #[must_use]
    pub fn mantissa(&self) -> &BigInt {
        &self.mant
    }

    /// The binary exponent: the value is `mantissa * 2^exponent`.
    #[must_use]
    pub fn exponent(&self) -> i64 {
        self.exp
    }

    /// The mantissa width in bits.
    #[must_use]
    pub fn precision(&self) -> u32 {
        self.prec
    }

    /// Whether the value is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.mant.is_zero()
    }

    /// `-1`, `0` or `1`.
    #[must_use]
    pub fn sign(&self) -> i8 {
        match self.mant.sign() {
            | Sign::Minus => -1,
            | Sign::NoSign => 0,
            | Sign::Plus => 1,
        }
    }

    /// The value negated.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            mant: -self.mant.clone(),
            exp: self.exp,
            prec: self.prec,
        }
    }

    /// The absolute value.
    #[must_use]
    pub fn abs(&self) -> Self {
        Self {
            mant: BigInt::from_biguint(Sign::Plus, self.mant.magnitude().clone()),
            exp: self.exp,
            prec: self.prec,
        }
    }

    /// The same value rounded to `prec` bits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn with_precision(
        &self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        Self::rounded(self.mant.clone(), self.exp, prec, mode)
    }

    /// `floor(log2(|x|))`, or `None` for zero.
    #[must_use]
    pub fn floor_log2(&self) -> Option<i64> {
        if self.mant.is_zero() {
            return None;
        }
        let bits = i64::try_from(self.mant.magnitude().bits()).unwrap_or(i64::MAX);
        Some(self.exp.saturating_add(bits).saturating_sub(1))
    }

    /// The sum, rounded once to `prec` bits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn add(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        self.add_with_inexact(other, prec, mode)
            .map(|(value, _)| value)
    }

    pub(super) fn add_with_inexact(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<(Self, bool), String> {
        if self.mant.is_zero() {
            return Self::rounded_with_inexact(other.mant.clone(), other.exp, prec, mode);
        }
        if other.mant.is_zero() {
            return Self::rounded_with_inexact(self.mant.clone(), self.exp, prec, mode);
        }
        let (high, low) = if self.exp >= other.exp {
            (self, other)
        } else {
            (other, self)
        };
        let shift = high.exp - low.exp;
        if shift > 2 * MAX_EXP {
            return Err("exponent out of range".to_string());
        }
        let aligned = &high.mant << usize::try_from(shift).unwrap_or(0);
        let sum = aligned + &low.mant;
        Self::rounded_with_inexact(sum, low.exp, prec, mode)
    }

    /// The difference, rounded once to `prec` bits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn sub(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        self.add(&other.neg(), prec, mode)
    }

    /// The product, rounded once to `prec` bits.
    ///
    /// # Errors
    /// Fails when the exponent would leave the representable range.
    pub fn mul(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        self.mul_with_inexact(other, prec, mode)
            .map(|(value, _)| value)
    }

    pub(super) fn mul_with_inexact(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<(Self, bool), String> {
        if self.mant.is_zero() || other.mant.is_zero() {
            return Ok((Self::zero(prec), false));
        }
        let product = &self.mant * &other.mant;
        let exponent = self
            .exp
            .checked_add(other.exp)
            .ok_or_else(|| "exponent out of range".to_string())?;
        Self::rounded_with_inexact(product, exponent, prec, mode)
    }

    /// The quotient, rounded once to `prec` bits.
    ///
    /// # Errors
    /// Fails on division by zero or when the exponent leaves the range.
    pub fn div(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        self.div_with_inexact(other, prec, mode)
            .map(|(value, _)| value)
    }

    /// Returns the rounded quotient and whether its exact value differs.
    pub(super) fn div_with_inexact(
        &self,
        other: &Self,
        prec: u32,
        mode: Round,
    ) -> Result<(Self, bool), String> {
        if other.mant.is_zero() {
            return Err("division by zero".to_string());
        }
        if self.mant.is_zero() {
            return Ok((Self::zero(prec), false));
        }
        let negative = self.sign() != other.sign();
        // Scale the numerator so the integer quotient carries `prec + 8`
        // significant bits; the remainder then only has to say "not exact".
        let numerator_bits = i64::try_from(self.mant.magnitude().bits()).unwrap_or(i64::MAX);
        let denominator_bits = i64::try_from(other.mant.magnitude().bits()).unwrap_or(i64::MAX);
        let want = i64::from(prec.clamp(1, MAX_BITS)) + 8;
        let shift = want + denominator_bits - numerator_bits;
        let (numerator, mut exponent, discarded) = if shift >= 0 {
            (
                self.mant.magnitude() << usize::try_from(shift).unwrap_or(0),
                self.exp
                    .checked_sub(shift)
                    .ok_or_else(|| "exponent out of range".to_string())?,
                false,
            )
        } else {
            let right_shift = shift.unsigned_abs();
            let discarded = self
                .mant
                .magnitude()
                .trailing_zeros()
                .is_some_and(|zeros| zeros < right_shift);
            (
                self.mant.magnitude() >> right_shift,
                self.exp
                    .checked_add(-shift)
                    .ok_or_else(|| "exponent out of range".to_string())?,
                discarded,
            )
        };
        exponent = exponent
            .checked_sub(other.exp)
            .ok_or_else(|| "exponent out of range".to_string())?;
        let (quotient, remainder) = numerator.div_rem(other.mant.magnitude());
        // A nonzero remainder or discarded numerator bits put the exact
        // quotient strictly between `q` and `q + 1`. Using `2q + 1` at the next
        // exponent preserves that rounding decision in every mode.
        let (magnitude, exponent) = if remainder.is_zero() && !discarded {
            (quotient, exponent)
        } else {
            (
                (quotient << 1_u32) + BigUint::from(1_u32),
                exponent
                    .checked_sub(1)
                    .ok_or_else(|| "exponent out of range".to_string())?,
            )
        };
        Self::rounded_with_inexact(
            BigInt::from_biguint(
                if negative {
                    Sign::Minus
                } else {
                    Sign::Plus
                },
                magnitude,
            ),
            exponent,
            prec,
            mode,
        )
    }

    /// The square root, rounded once to `prec` bits.
    ///
    /// # Errors
    /// Fails for a negative operand or an exponent outside the range.
    pub fn sqrt(
        &self,
        prec: u32,
        mode: Round,
    ) -> Result<Self, String> {
        if self.sign() < 0 {
            return Err("square root of a negative number".to_string());
        }
        if self.mant.is_zero() {
            return Ok(Self::zero(prec));
        }
        let want = prec.clamp(1, MAX_BITS);
        // An even exponent makes the root's exponent integral; scaling the
        // mantissa to `2 * want + 8` bits gives the integer square root
        // enough guard bits for one correct rounding.
        let mut exponent = self.exp;
        let mut magnitude = self.mant.magnitude().clone();
        if exponent % 2 != 0 {
            magnitude <<= 1_u32;
            exponent -= 1;
        }
        let target = 2 * u64::from(want) + 8;
        // The exponent must stay even for `sqrt(m * 2^e) = sqrt(m) 2^(e/2)`
        // to hold, so the scaling shift is rounded up to an even number.
        let shift = target.saturating_sub(magnitude.bits());
        let shift = shift + (shift % 2);
        magnitude <<= usize::try_from(shift).unwrap_or(0);
        let exponent = exponent
            .checked_sub(i64::try_from(shift).unwrap_or(i64::MAX))
            .ok_or_else(|| "exponent out of range".to_string())?;
        let root = magnitude.sqrt();
        let exact = &root * &root == magnitude;
        let (magnitude, exponent) = if exact {
            (root, exponent / 2)
        } else {
            (
                (root << 1_u32) + BigUint::from(1_u32),
                (exponent / 2)
                    .checked_sub(1)
                    .ok_or_else(|| "exponent out of range".to_string())?,
            )
        };
        Self::rounded(
            BigInt::from_biguint(Sign::Plus, magnitude),
            exponent,
            want,
            mode,
        )
    }

    /// `self * 2^k`, exact.
    ///
    /// # Errors
    /// Fails when the exponent leaves the representable range.
    pub fn scale_pow2(
        &self,
        k: i64,
    ) -> Result<Self, String> {
        let exponent = self
            .exp
            .checked_add(k)
            .ok_or_else(|| "exponent out of range".to_string())?;
        if exponent.unsigned_abs() > MAX_EXP.unsigned_abs() {
            return Err("exponent out of range".to_string());
        }
        Ok(Self {
            mant: self.mant.clone(),
            exp: if self.is_zero() {
                0
            } else {
                exponent
            },
            prec: self.prec,
        })
    }

    /// Total order: `-1`, `0` or `1` as `self` is below, equal to or above
    /// `other`.
    #[must_use]
    pub fn cmp_value(
        &self,
        other: &Self,
    ) -> i8 {
        let (a, b) = (self.sign(), other.sign());
        if a != b {
            return if a < b { -1 } else { 1 };
        }
        if a == 0 {
            return 0;
        }
        let ordering = match (self.floor_log2(), other.floor_log2()) {
            | (Some(x), Some(y)) if x != y => {
                let above = if x > y { 1 } else { -1 };
                return if a < 0 { -above } else { above };
            },
            | (None, None) => return 0,
            | _ => 0,
        };
        let _ = ordering;
        let magnitude_order = match self.exp.cmp(&other.exp) {
            | std::cmp::Ordering::Equal => self.mant.magnitude().cmp(other.mant.magnitude()),
            | std::cmp::Ordering::Greater => {
                let shift = usize::try_from(self.exp - other.exp).unwrap_or(0);
                (self.mant.magnitude() << shift).cmp(other.mant.magnitude())
            },
            | std::cmp::Ordering::Less => {
                let shift = usize::try_from(other.exp - self.exp).unwrap_or(0);
                self.mant
                    .magnitude()
                    .cmp(&(other.mant.magnitude() << shift))
            },
        };
        let above = match magnitude_order {
            | std::cmp::Ordering::Less => -1,
            | std::cmp::Ordering::Equal => 0,
            | std::cmp::Ordering::Greater => 1,
        };
        if a < 0 { -above } else { above }
    }

    /// Whether the value is an exact integer.
    #[must_use]
    pub fn is_integer(&self) -> bool {
        if self.mant.is_zero() || self.exp >= 0 {
            return true;
        }
        let drop = u64::try_from(-self.exp).unwrap_or(u64::MAX);
        if drop >= self.mant.magnitude().bits() {
            return false;
        }
        let mask =
            (BigUint::from(1_u32) << usize::try_from(drop).unwrap_or(0)) - BigUint::from(1_u32);
        (self.mant.magnitude() & mask).is_zero()
    }

    /// The value truncated toward zero to an integer.
    ///
    /// # Errors
    /// Fails when the exponent leaves the representable range.
    pub fn trunc_to_bigint(&self) -> Result<BigInt, String> {
        if self.mant.is_zero() {
            return Ok(BigInt::from(0));
        }
        if self.exp >= 0 {
            return Ok(&self.mant << usize::try_from(self.exp).unwrap_or(0));
        }
        let drop = u64::try_from(-self.exp).unwrap_or(u64::MAX);
        if drop >= self.mant.magnitude().bits() {
            return Ok(BigInt::from(0));
        }
        Ok(BigInt::from_biguint(
            self.mant.sign(),
            self.mant.magnitude() >> drop,
        ))
    }

    /// The value as an `f64`, rounded to nearest-even.
    ///
    /// Subnormals, overflow to infinity and underflow to zero follow IEEE-754,
    /// so this is the correctly rounded conversion of the arbitrary-precision
    /// value rather than a conversion of a `f64` approximation of it.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        if self.mant.is_zero() {
            return 0.0;
        }
        let negative = self.sign() < 0;
        let Some(top) = self.floor_log2() else {
            return 0.0;
        };
        if top > 1023 {
            return if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if top < -1075 {
            return if negative { -0.0 } else { 0.0 };
        }
        let mut grid = if top < -1022 {
            -1074
        } else {
            top - 52
        };
        let Some(mut n) = self.round_to_grid(grid) else {
            return if negative { -0.0 } else { 0.0 };
        };
        if n.is_zero() {
            return if negative { -0.0 } else { 0.0 };
        }
        if n.magnitude().bits() > 53 {
            grid += 1;
            let Some(adjusted) = self.round_to_grid(grid) else {
                return if negative { -0.0 } else { 0.0 };
            };
            n = adjusted;
        }
        let magnitude = n.magnitude().clone();
        let bits = magnitude.bits();
        let (exponent_field, mantissa_field) = if bits >= 53 {
            let exponent = grid + 52 + 1023;
            if exponent >= 2047 {
                return if negative {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                };
            }
            let field = &magnitude - (BigUint::from(1_u32) << 52_u32);
            (u64::try_from(exponent).unwrap_or(2047), field)
        } else {
            // Subnormal, or the smallest normal after a carry.
            let exponent = grid + i64::try_from(bits).unwrap_or(0) - 1 + 1023;
            if exponent > 0 {
                let shift = 53 - bits;
                let scaled = &magnitude << shift;
                let field = &scaled - (BigUint::from(1_u32) << 52_u32);
                (u64::try_from(exponent).unwrap_or(2047), field)
            } else {
                (0, magnitude)
            }
        };
        let low = mantissa_field.to_u64_digits();
        let mantissa = low.first().copied().unwrap_or(0) & ((1_u64 << 52) - 1);
        let sign_bit = u64::from(negative) << 63;
        f64::from_bits(sign_bit | (exponent_field << 52) | mantissa)
    }

    /// `round(self / 2^grid)` to nearest-even, as an integer.
    fn round_to_grid(
        &self,
        grid: i64,
    ) -> Option<BigInt> {
        let shift = self.exp.checked_sub(grid)?;
        if shift >= 0 {
            return Some(&self.mant << usize::try_from(shift).ok()?);
        }
        let drop = u64::try_from(-shift).ok()?;
        if drop > self.mant.magnitude().bits() + 1 {
            return Some(BigInt::from(0));
        }
        let magnitude = self.mant.magnitude().clone();
        let kept = &magnitude >> drop;
        let rest = &magnitude - (&kept << drop);
        let half = BigUint::from(1_u32) << usize::try_from(drop - 1).ok()?;
        let up = rest > half || (rest == half && kept.bit(0));
        let magnitude = if up {
            kept + BigUint::from(1_u32)
        } else {
            kept
        };
        Some(BigInt::from_biguint(
            if self.sign() < 0 {
                Sign::Minus
            } else {
                Sign::Plus
            },
            magnitude,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trunc_to_bigint_discards_fractional_bits_toward_zero() {
        for (mantissa, expected) in [(-15, -3), (-16, -4), (-1, 0), (15, 3)] {
            let Ok(value) = Real::rounded(BigInt::from(mantissa), -2, 8, Round::NearestEven) else {
                panic!("the dyadic test value must be exactly representable");
            };
            assert_eq!(
                value.trunc_to_bigint(),
                Ok(BigInt::from(expected)),
                "{mantissa} / 4"
            );
        }
    }
}
