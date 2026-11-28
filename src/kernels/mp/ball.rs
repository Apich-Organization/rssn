//! Balls: a midpoint with a radius that bounds the error of the midpoint.
//!
//! A [`Ball`] represents every real in `[mid - rad, mid + rad]`. The basic
//! operations propagate the radius by the usual formulas, and the elementary
//! functions either apply a monotone map to the endpoints or a Lipschitz
//! bound to the radius, so the enclosure is preserved without any assumption
//! about the true value. A [`Ball`] is what the `mp` operators return when a
//! host asks for an enclosure. An enclosure rounds its decimal bounds
//! outward. A rounded value is certified when both endpoints round to the
//! same decimal text.

use num_bigint::BigInt;

use super::real::MAX_EXP;
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
        let (mid, inexact) = self
            .mid
            .add_with_inexact(&other.mid, prec, Round::NearestEven)?;
        let rounding = midpoint_rounding(&mid, inexact)?;
        let rad = self.rad.add(&other.rad, prec, Round::AwayFromZero)?.add(
            &rounding,
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
        let (mid, inexact) = self
            .mid
            .mul_with_inexact(&other.mid, prec, Round::NearestEven)?;
        let rounding = midpoint_rounding(&mid, inexact)?;
        let left = self.mid.abs().mul(&other.rad, prec, Round::AwayFromZero)?;
        let right = other.mid.abs().mul(&self.rad, prec, Round::AwayFromZero)?;
        let corner = self.rad.mul(&other.rad, prec, Round::AwayFromZero)?;
        let rad = left
            .add(&right, prec, Round::AwayFromZero)?
            .add(&corner, prec, Round::AwayFromZero)?
            .add(&rounding, prec, Round::AwayFromZero)?;
        Ok(Self { mid, rad })
    }

    /// Encloses the quotients of the two intervals when the divisor excludes zero.
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
        let (mid, inexact) = self
            .mid
            .div_with_inexact(&other.mid, prec, Round::NearestEven)?;
        let rounding = midpoint_rounding(&mid, inexact)?;
        if self.rad.is_zero() && other.rad.is_zero() {
            return Ok(Self { mid, rad: rounding });
        }

        let minimum_divisor = other
            .mid
            .abs()
            .sub(&other.rad, prec, Round::TowardNegative)?;
        let propagated = if other.rad.is_zero() {
            Real::zero(prec)
        } else {
            let abs_mid = mid.abs();
            let abs_quotient = abs_mid.add(&rounding, prec, Round::TowardPositive);
            match abs_quotient {
                | Ok(bound) => bound.mul(&other.rad, prec, Round::TowardPositive)?,
                | Err(_) => {
                    // Scale the two terms separately if their unscaled sum
                    // would leave the representable range.
                    let central = abs_mid.mul(&other.rad, prec, Round::TowardPositive)?;
                    let error = rounding.mul(&other.rad, prec, Round::TowardPositive)?;
                    central.add(&error, prec, Round::TowardPositive)?
                },
            }
        };
        // For q = a/b, (a + da)/(b + db) - q = (da - q*db)/(b + db).
        // Bound the numerator upward and the positive divisor magnitude downward.
        let variation = propagated
            .add(&self.rad, prec, Round::TowardPositive)?
            .div(&minimum_divisor, prec, Round::TowardPositive)?;
        let rad = rounding.add(&variation, prec, Round::TowardPositive)?;
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
        let midpoint = self.mid.abs();
        if midpoint.cmp_value(&self.rad) >= 0 {
            let mid = midpoint.with_precision(prec, Round::NearestEven)?;
            let rounding = midpoint.sub(&mid, prec, Round::AwayFromZero)?.abs();
            let rad = self.rad.add(&rounding, prec, Round::TowardPositive)?;
            return Ok(Self { mid, rad });
        }
        let upper = midpoint.add(&self.rad, prec, Round::TowardPositive)?;
        if let Ok(mid) = upper.scale_pow2(-1) {
            return Ok(Self {
                rad: mid.clone(),
                mid,
            });
        }
        // At the exponent floor the half-width is unrepresentable.
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
        let (lower, upper) = self.endpoints(inner)?;
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
            self.mid.add(&self.rad, prec, Round::TowardPositive)?,
        ))
    }
}

/// Outward bound for the rounding error of an exact midpoint operation.
fn midpoint_rounding(
    mid: &Real,
    inexact: bool,
) -> Result<Real, String> {
    let prec = mid.precision();
    if !inexact {
        return Ok(Real::zero(prec));
    }
    // The exact value is at most half an ulp from the rounded midpoint.
    // A carry only widens the spacing; at the exponent floor use a full ulp.
    let exponent = mid.exponent().saturating_sub(1).max(-MAX_EXP);
    let unit = Real::rounded(BigInt::from(1_u32), exponent, 1, Round::TowardPositive)?;
    Ok(unit
        .with_precision(prec, Round::TowardPositive)
        .unwrap_or(unit))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_precision_arithmetic_encloses_point_rounding_and_interval_corners() {
        for precision in 1..=4 {
            for left in -8_i64..=8 {
                for right in -8_i64..=8 {
                    for radius in [0, 1] {
                        let a = Ball::new(
                            Real::from_i64(left, 16).unwrap(),
                            Real::from_i64(radius, 16).unwrap(),
                        )
                        .unwrap();
                        let b = Ball::new(
                            Real::from_i64(right, 16).unwrap(),
                            Real::from_i64(radius, 16).unwrap(),
                        )
                        .unwrap();
                        let results = [
                            a.add(&b, precision).unwrap(),
                            a.sub(&b, precision).unwrap(),
                            a.mul(&b, precision).unwrap(),
                        ];
                        for a_end in [left - radius, left + radius] {
                            for b_end in [right - radius, right + radius] {
                                let exact = [a_end + b_end, a_end - b_end, a_end * b_end];
                                for (value, point) in results.iter().zip(exact) {
                                    let point = Real::from_i64(point, 16).unwrap();
                                    let (lo, hi) = value.endpoints(16).unwrap();
                                    assert!(lo.cmp_value(&point) <= 0);
                                    assert!(hi.cmp_value(&point) >= 0);
                                }
                            }
                        }
                    }
                }
            }
        }
        let two = Ball::point(Real::from_i64(2, 16).unwrap());
        assert!(two.add(&two, 1).unwrap().rad.is_zero());
        assert!(two.sub(&two, 1).unwrap().rad.is_zero());
        assert!(two.mul(&two, 1).unwrap().rad.is_zero());
    }

    #[test]
    fn lower_precision_keeps_negative_upper_endpoints() {
        let ball = Ball::new(
            Real::from_i64(-3, 16).unwrap(),
            Real::from_f64(0.25, 16).unwrap(),
        )
        .unwrap();
        let (lower, upper) = ball.endpoints(2).unwrap();
        assert_eq!(lower.to_f64(), -4.0);
        assert_eq!(upper.to_f64(), -2.0);
    }

    #[test]
    fn absolute_value_encloses_point_after_rounding_the_center() {
        let point = Ball::point(Real::from_f64(-2.75, 16).unwrap());
        let (lower, upper) = point.abs(2).unwrap().endpoints(16).unwrap();
        assert!(lower.to_f64() <= 2.75 && upper.to_f64() >= 2.75);

        let crossing = Ball::new(
            Real::from_f64(-0.5, 16).unwrap(),
            Real::from_i64(1, 16).unwrap(),
        )
        .unwrap();
        let (lower, upper) = crossing.abs(2).unwrap().endpoints(16).unwrap();
        assert!(lower.is_zero());
        assert_eq!(upper.to_f64(), 1.5);
    }

    #[test]
    fn absolute_value_encloses_values_at_the_exponent_floor() {
        use super::super::real::MAX_EXP;

        let radius = Real::from_i64(1, 2)
            .unwrap()
            .scale_pow2(MAX_EXP.saturating_neg().saturating_add(1))
            .unwrap();
        let symmetric = Ball::new(Real::zero(2), radius).unwrap();
        let absolute = symmetric.abs(2).unwrap();
        let (lower, upper) = absolute.endpoints(2).unwrap();
        assert!(lower.sign() <= 0);
        assert!(upper.cmp_value(symmetric.radius()) >= 0);
    }
}
