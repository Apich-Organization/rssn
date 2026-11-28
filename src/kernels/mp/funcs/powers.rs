//! Certified powers with exact integer and rational exponent entry points.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;

use super::Ball;
use super::MAX_BITS;
use super::MAX_EXP;
use super::Real;
use super::ScaledBounds;
use super::exp_ball;
use super::ln;
use super::ln_ball;

const GUARD_BITS: u32 = 48;

/// Nonnegative directed integer bounds with a shared binary exponent.
/// Exponents stay wide until the final dyadic can be checked.
struct Positive {
    lower: BigUint,
    upper: BigUint,
    exponent: i128,
}

impl Positive {
    fn enclose_size(
        lower: BigUint,
        upper: BigUint,
        exponent: i128,
        work: u64,
    ) -> Self {
        let drop = upper.bits().saturating_sub(work);
        Self {
            lower: lower >> drop,
            upper: ((upper - 1_u32) >> drop) + 1_u32,
            exponent: exponent + i128::from(drop),
        }
    }

    fn from_base(
        x: &Real,
        reciprocal: bool,
        work: u64,
    ) -> Self {
        let magnitude = x.mantissa().magnitude();
        if reciprocal {
            let shift = work + magnitude.bits();
            let unit = BigUint::from(1_u32) << shift;
            return Self::enclose_size(
                &unit / magnitude,
                unit.div_ceil(magnitude),
                -i128::from(x.exponent()) - i128::from(shift),
                work,
            );
        }
        Self::enclose_size(
            magnitude.clone(),
            magnitude.clone(),
            i128::from(x.exponent()),
            work,
        )
    }

    fn multiply(
        &self,
        other: &Self,
        work: u64,
    ) -> Self {
        Self::enclose_size(
            &self.lower * &other.lower,
            &self.upper * &other.upper,
            self.exponent + other.exponent,
            work,
        )
    }

    fn finish(
        self,
        negative: bool,
    ) -> Result<ScaledBounds, String> {
        let top = self.exponent + i128::from(self.upper.bits()) - 1;
        // A value below the smallest representable dyadic can round up
        // by one binade. The final Real rounding decides boundary cases.
        let min = -i128::from(MAX_EXP) - 1;
        let max = i128::from(MAX_EXP) + i128::from(MAX_BITS);
        if !(min..=max).contains(&top) {
            return Err("power result exponent out of range".into());
        }
        let exponent = i64::try_from(self.exponent).map_err(|error| error.to_string())?;
        let lower = BigInt::from(self.lower);
        let upper = BigInt::from(self.upper);
        let (lower, upper) = if negative {
            (-upper, -lower)
        } else {
            (lower, upper)
        };
        Ok(ScaledBounds {
            lower,
            upper,
            exponent,
        })
    }
}

fn exact_unit(value: i32) -> ScaledBounds {
    ScaledBounds {
        lower: value.into(),
        upper: value.into(),
        exponent: 0,
    }
}

fn integer_bounds(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<ScaledBounds, String> {
    if n == 0 {
        return Ok(exact_unit(1));
    }
    if x.is_zero() {
        return if n < 0 {
            Err("zero to a negative power".into())
        } else {
            Ok(exact_unit(0))
        };
    }
    let negative = x.sign() < 0 && n % 2 != 0;
    if x.floor_log2() == Some(0)
        && x.mantissa().trailing_zeros() == Some(u64::from(x.precision() - 1))
    {
        return Ok(exact_unit(if negative { -1 } else { 1 }));
    }
    let mut power = n.unsigned_abs();
    // Relative error may double at each square. Keep one extra working
    // bit for each exponent bit, in addition to the rounding guard.
    let work = u64::from(prec.clamp(1, MAX_BITS) + GUARD_BITS + power.ilog2() + 1);
    let mut base = Positive::from_base(x, n < 0, work);
    let mut product = Positive {
        lower: 1_u32.into(),
        upper: 1_u32.into(),
        exponent: 0,
    };
    loop {
        if power & 1 == 1 {
            product = product.multiply(&base, work);
        }
        power >>= 1;
        if power == 0 {
            break;
        }
        base = base.multiply(&base, work);
    }
    product.finish(negative)
}

fn ball(
    bounds: ScaledBounds,
    prec: u32,
) -> Result<Ball, String> {
    Ball::from_scaled_bounds(bounds.lower, bounds.upper, bounds.exponent, prec)
}

fn hull(
    left: ScaledBounds,
    right: ScaledBounds,
) -> ScaledBounds {
    let exponent = left.exponent.min(right.exponent);
    let a = (left.exponent - exponent).unsigned_abs();
    let b = (right.exponent - exponent).unsigned_abs();
    ScaledBounds {
        lower: (left.lower << a).min(right.lower << b),
        upper: (left.upper << a).max(right.upper << b),
        exponent,
    }
}

fn working(prec: u32) -> u32 {
    prec.saturating_add(GUARD_BITS).clamp(1, MAX_BITS)
}

pub(super) fn integer(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<Ball, String> {
    ball(integer_bounds(x, n, prec)?, prec)
}

pub(super) fn integer_interval(
    x: &Ball,
    n: i64,
    prec: u32,
) -> Result<Ball, String> {
    if n == 0 || x.rad.is_zero() {
        return integer(&x.mid, n, prec);
    }
    let includes_zero = x.contains_zero();
    if n < 0 && includes_zero {
        return Err("negative integer power of an interval containing zero".into());
    }
    let work = working(prec).max(x.mid.precision()).max(x.rad.precision());
    let (low, high) = x.endpoints(work)?;
    let mut bounds = hull(
        integer_bounds(&low, n, prec)?,
        integer_bounds(&high, n, prec)?,
    );
    if includes_zero && n % 2 == 0 {
        bounds.lower = BigInt::from(0);
    }
    ball(bounds, prec)
}

pub(super) fn point(
    x: &Real,
    y: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if y.is_zero() {
        return ball(exact_unit(1), prec);
    }
    if x.is_zero() {
        return if y.sign() < 0 {
            Err("zero to a negative power".into())
        } else {
            ball(exact_unit(0), prec)
        };
    }
    if y.is_integer() {
        return integer(x, y.integer_exponent()?, prec);
    }
    if x.sign() < 0 {
        return Err("a negative base needs an integer or admissible rational exponent".into());
    }
    let work = working(prec);
    let exponent = ln(x, work)?.mul(&Ball::point(y.clone()), work)?;
    exp_ball(&exponent, prec)
}

pub(super) fn interval(
    x: &Ball,
    y: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    if x.rad.is_zero() && y.rad.is_zero() {
        return point(&x.mid, &y.mid, prec);
    }
    if y.rad.is_zero() && y.mid.is_integer() {
        return integer_interval(x, y.mid.integer_exponent()?, prec);
    }
    let work = working(prec);
    let exponent = ln_ball(x, work)?.mul(y, work)?;
    exp_ball(&exponent, prec)
}

#[derive(Clone, Copy)]
pub(super) struct Rational {
    numerator: i64,
    denominator: u32,
}

impl Rational {
    pub(super) fn reduce(
        numerator: i64,
        denominator: u32,
    ) -> Result<Self, String> {
        if denominator == 0 {
            return Err("rational power denominator must be positive".into());
        }
        let divisor = numerator.unsigned_abs().gcd(&u64::from(denominator));
        let divisor = u32::try_from(divisor).map_err(|error| error.to_string())?;
        Ok(Self {
            numerator: numerator / i64::from(divisor),
            denominator: denominator / divisor,
        })
    }

    pub(super) fn point(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<Ball, String> {
        if self.denominator == 1 {
            return integer(x, self.numerator, prec);
        }
        if x.is_zero() {
            return if self.numerator < 0 {
                Err("zero to a negative power".into())
            } else {
                ball(exact_unit(0), prec)
            };
        }
        if x.sign() < 0 && self.denominator.is_multiple_of(2) {
            return Err("a negative base needs an odd reduced denominator".into());
        }
        let work = working(prec);
        let numerator = Ball::point(Real::from_i64(self.numerator, 64)?);
        let denominator = Ball::point(Real::from_i64(i64::from(self.denominator), 32)?);
        let logarithm = ln(&x.abs(), work)?;
        let exponent = logarithm.mul(&numerator, work)?.div(&denominator, work)?;
        let value = exp_ball(&exponent, prec)?;
        Ok(if x.sign() < 0 && self.numerator % 2 != 0 {
            value.neg()
        } else {
            value
        })
    }

    pub(super) fn interval(
        self,
        x: &Ball,
        prec: u32,
    ) -> Result<Ball, String> {
        if self.denominator == 1 {
            return integer_interval(x, self.numerator, prec);
        }
        if x.rad.is_zero() {
            return self.point(&x.mid, prec);
        }
        if self.numerator < 0 && x.contains_zero() {
            return Err("negative rational power of an interval containing zero".into());
        }
        let work = working(prec).max(x.mid.precision()).max(x.rad.precision());
        let (x_low, x_high) = x.endpoints(work)?;
        let (low_a, high_a) = self.point(&x_low, work)?.endpoints(work)?;
        let (low_b, high_b) = self.point(&x_high, work)?.endpoints(work)?;
        let mut lower = if low_a.cmp_value(&low_b) <= 0 {
            low_a
        } else {
            low_b
        };
        let upper = if high_a.cmp_value(&high_b) >= 0 {
            high_a
        } else {
            high_b
        };
        if x.contains_zero() && lower.sign() > 0 {
            lower = Real::zero(work);
        }
        let exponent = lower.exponent().min(upper.exponent());
        let low = lower.mantissa() << (lower.exponent() - exponent).unsigned_abs();
        let high = upper.mantissa() << (upper.exponent() - exponent).unsigned_abs();
        Ball::from_scaled_bounds(low, high, exponent, prec)
    }
}
