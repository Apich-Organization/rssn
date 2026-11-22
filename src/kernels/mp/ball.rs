//! Balls: a midpoint with a radius that bounds the error of the midpoint.
//!
//! A [`Ball`] represents every real in `[mid - rad, mid + rad]`. The basic
//! operations propagate the radius by the usual formulas, and the elementary
//! functions either apply a monotone map to the endpoints or a Lipschitz
//! bound to the radius, so the enclosure is preserved without any assumption
//! about the true value. A [`Ball`] is what the `mp` operators return when a
//! host asks for an enclosure rather than a rounded value: the decimal
//! rendering only succeeds when the two endpoints round to the same text,
//! which is exactly the condition for the digits to be certified.

use num_bigint::BigInt;

use super::real::Real;
use super::real::Round;

/// An interval `[mid - rad, mid + rad]` with `rad >= 0`.
#[derive(Clone, Debug)]
pub struct Ball {
    /// The centre of the interval.
    pub mid: Real,
    /// The non-negative radius.
    pub rad: Real,
}

impl Ball {
    /// The ball of radius `rad` around `mid`.
    ///
    /// # Errors
    /// Fails for a negative radius or when an exponent leaves the range.
    pub fn new(
        mid: Real,
        rad: Real,
    ) -> Result<Self, String> {
        if rad.sign() < 0 {
            return Err("a radius cannot be negative".to_string());
        }
        Ok(Self { mid, rad })
    }

    /// The degenerate ball holding exactly `mid`.
    #[must_use]
    pub fn point(mid: Real) -> Self {
        let prec = mid.precision();
        Self {
            mid,
            rad: Real::zero(prec),
        }
    }

    /// The midpoint.
    #[must_use]
    pub fn midpoint(&self) -> &Real {
        &self.mid
    }

    /// The radius.
    #[must_use]
    pub fn radius(&self) -> &Real {
        &self.rad
    }

    /// Whether the interval contains zero.
    #[must_use]
    pub fn contains_zero(&self) -> bool {
        self.mid.abs().cmp_value(&self.rad) <= 0
    }

    /// The width `2 * rad`.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn width(&self) -> Result<Real, String> {
        self.rad
            .add(&self.rad, self.rad.precision(), Round::AwayFromZero)
    }

    /// The sum.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn add(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let mid = self.mid.add(&other.mid, prec, Round::NearestEven)?;
        let rad = self.rad.add(&other.rad, prec, Round::AwayFromZero)?.add(
            &slack(&mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?;
        Ok(Self { mid, rad })
    }

    /// The difference.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn sub(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        self.add(&other.neg(), prec)
    }

    /// The product.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn mul(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let mid = self.mid.mul(&other.mid, prec, Round::NearestEven)?;
        let left = self.mid.abs().mul(&other.rad, prec, Round::AwayFromZero)?;
        let right = other.mid.abs().mul(&self.rad, prec, Round::AwayFromZero)?;
        let corner = self.rad.mul(&other.rad, prec, Round::AwayFromZero)?;
        let rad = left
            .add(&right, prec, Round::AwayFromZero)?
            .add(&corner, prec, Round::AwayFromZero)?
            .add(&slack(&mid, prec)?, prec, Round::AwayFromZero)?;
        Ok(Self { mid, rad })
    }

    /// The quotient, when the divisor's interval excludes zero.
    ///
    /// # Errors
    /// Fails when the divisor contains zero or an exponent leaves the range.
    pub fn div(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        if other.contains_zero() {
            return Err("division by an interval containing zero".to_string());
        }
        let mid = self.mid.div(&other.mid, prec, Round::NearestEven)?;
        let denominator = other.mid.abs().sub(&other.rad, prec, Round::TowardZero)?;
        if denominator.sign() <= 0 {
            return Err("division by an interval containing zero".to_string());
        }
        let numerator = other
            .mid
            .abs()
            .mul(&self.rad, prec, Round::AwayFromZero)?
            .add(
                &self.mid.abs().mul(&other.rad, prec, Round::AwayFromZero)?,
                prec,
                Round::AwayFromZero,
            )?
            .add(
                &self.rad.mul(&other.rad, prec, Round::AwayFromZero)?,
                prec,
                Round::AwayFromZero,
            )?;
        let scale = denominator.mul(&other.mid.abs(), prec, Round::TowardZero)?;
        let rad = numerator.div(&scale, prec, Round::AwayFromZero)?.add(
            &slack(&mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?;
        Ok(Self { mid, rad })
    }

    /// The value negated.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            mid: self.mid.neg(),
            rad: self.rad.clone(),
        }
    }

    /// The absolute value, as the hull of the endpoints' absolute values.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn abs(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let lower = self.mid.abs().sub(&self.rad, prec, Round::TowardZero)?;
        let upper = self.mid.abs().add(&self.rad, prec, Round::AwayFromZero)?;
        if lower.sign() >= 0 {
            return Ok(Self {
                mid: lower
                    .add(&upper, prec, Round::NearestEven)?
                    .scale_pow2(-1)?,
                rad: upper
                    .sub(&lower, prec, Round::AwayFromZero)?
                    .scale_pow2(-1)?,
            });
        }
        Ok(Self {
            mid: Real::zero(prec),
            rad: upper,
        })
    }

    /// Applies a monotone function to the endpoints.
    ///
    /// `f` must be correctly rounded in the given direction, and monotone on
    /// the interval; `increasing` says which way.
    ///
    /// # Errors
    /// Fails when `f` fails, which includes a domain error.
    pub fn map_monotone(
        &self,
        f: &dyn Fn(&Real, u32, Round) -> Result<Real, String>,
        increasing: bool,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(16);
        let lower = self.mid.sub(&self.rad, inner, Round::TowardNegative)?;
        let upper = self.mid.add(&self.rad, inner, Round::AwayFromZero)?;
        let mid = f(&self.mid, inner, Round::NearestEven)?;
        let (low, high) = if increasing {
            (
                f(&lower, inner, Round::TowardNegative)?,
                f(&upper, inner, Round::TowardPositive)?,
            )
        } else {
            (
                f(&upper, inner, Round::TowardNegative)?,
                f(&lower, inner, Round::TowardPositive)?,
            )
        };
        let below = mid.sub(&low, inner, Round::AwayFromZero)?;
        let above = high.sub(&mid, inner, Round::AwayFromZero)?;
        let rad = if below.cmp_value(&above) >= 0 {
            below
        } else {
            above
        };
        Ok(Self {
            mid: mid.with_precision(prec, Round::NearestEven)?,
            rad: rad.add(&slack(&mid, prec)?, prec, Round::AwayFromZero)?,
        })
    }

    /// Widens the radius by a Lipschitz bound of the function applied to the
    /// midpoint.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn lipschitz(
        &self,
        mid: Real,
        bound: &Real,
        prec: u32,
    ) -> Result<Self, String> {
        let rad = self.rad.mul(bound, prec, Round::AwayFromZero)?.add(
            &slack(&mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?;
        Ok(Self { mid, rad })
    }

    /// The two endpoints, `mid - rad` and `mid + rad`.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range.
    pub fn endpoints(
        &self,
        prec: u32,
    ) -> Result<(Real, Real), String> {
        Ok((
            self.mid.sub(&self.rad, prec, Round::TowardNegative)?,
            self.mid.add(&self.rad, prec, Round::AwayFromZero)?,
        ))
    }
}

/// One unit in the last place of `value` at precision `prec`, a bound on the
/// error of rounding `value` to `prec` bits.
fn slack(
    value: &Real,
    prec: u32,
) -> Result<Real, String> {
    // No rounding, no error: a value that is exactly representable at the
    // target precision contributes nothing to the radius. Without this an
    // exact cancellation (`exp(0) - exp(0)`) would carry a spurious radius
    // that stops a later rendering from certifying the digits of zero.
    let width = value.precision();
    if value.is_zero() || width <= prec {
        return Ok(Real::zero(1));
    }
    let dropped = u64::from(width - prec);
    if value
        .mantissa()
        .magnitude()
        .trailing_zeros()
        .is_some_and(|zeros| zeros >= dropped)
    {
        return Ok(Real::zero(1));
    }
    match value.floor_log2() {
        | None => Real::rounded(BigInt::from(1), -i64::from(prec), 1, Round::AwayFromZero),
        | Some(top) => {
            let exponent = top.saturating_sub(i64::from(prec)).saturating_add(1);
            Real::rounded(BigInt::from(1), exponent, 1, Round::AwayFromZero)
        },
    }
}
