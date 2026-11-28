//! Outward logarithmic bounds for exact dyadic inputs.

use num_bigint::BigInt;
use num_traits::Signed;

use super::ball::Ball;
use super::consts;
use super::fixed::Fx;
use super::real::MAX_BITS;
use super::real::Real;
use super::real::Round;

const LOG_GUARD_BITS: u32 = 40;

struct Bounds {
    low: BigInt,
    high: BigInt,
    scale: u32,
}

impl Bounds {
    fn ball(
        self,
        prec: u32,
    ) -> Result<Ball, String> {
        let inner = prec.saturating_add(LOG_GUARD_BITS);
        let exponent = -i64::from(self.scale);
        let lower = Real::rounded(self.low, exponent, inner, Round::TowardNegative)?;
        let upper = Real::rounded(self.high, exponent, inner, Round::TowardPositive)?;
        Ball::from_bounds(lower, upper, prec)
    }

    fn divide(
        self,
        divisor: &Self,
        prec: u32,
    ) -> Result<Ball, String> {
        let inner = prec.saturating_add(LOG_GUARD_BITS);
        let exponent = i64::from(divisor.scale) - i64::from(self.scale);
        let low_divisor = if self.low.sign() == num_bigint::Sign::Minus {
            &divisor.low
        } else {
            &divisor.high
        };
        let high_divisor = if self.high.sign() == num_bigint::Sign::Minus {
            &divisor.high
        } else {
            &divisor.low
        };
        let lower = Real::rounded_ratio(
            &self.low,
            low_divisor.magnitude(),
            exponent,
            inner,
            Round::TowardNegative,
        )?;
        let upper = Real::rounded_ratio(
            &self.high,
            high_divisor.magnitude(),
            exponent,
            inner,
            Round::TowardPositive,
        )?;
        Ball::from_bounds(lower, upper, prec)
    }
}

/// Bounds 2 atanh(numerator/denominator), where the ratio has magnitude ≤ 1/3.
fn twice_atanh(
    numerator: &BigInt,
    denominator: &BigInt,
    scale: u32,
) -> Result<Bounds, String> {
    let magnitude = Fx {
        v: (numerator.abs() << scale) / denominator,
        s: scale,
    };
    let multiplier = magnitude.mul(&magnitude);
    let mut power = magnitude;
    let mut sum = Fx::zero(scale);
    for k in 0..=i64::from(scale) {
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            // Work in units of 2^-scale. At the truncated argument, the
            // power error is at most two units: a step propagates at most
            // 1/9 of that error and introduces at most 1/3 + 1 units.
            // The term error is thus at most three units. The omitted
            // tail is below four units and truncating the argument adds
            // below two units because atanh'(x) is at most 9/8.
            let lower = &sum.v << 1_u32;
            let upper = (sum.v + (3 * k + 6)) << 1_u32;
            let (low, high) = if numerator.is_negative() {
                (-upper, -lower)
            } else {
                (lower, upper)
            };
            return Ok(Bounds { low, high, scale });
        }
        sum = sum.add(&term);
        power = power.mul(&multiplier);
    }
    Err("logarithm series did not converge".into())
}

fn natural_bounds(
    x: &Real,
    prec: u32,
) -> Result<Bounds, String> {
    if x.sign() <= 0 {
        return Err("logarithm of a non-positive number".into());
    }
    let Some(mut exponent) = x.floor_log2() else {
        return Err("logarithm of a non-positive number".into());
    };
    let magnitude = x.mantissa();
    let mut unit = BigInt::from(1) << (x.precision() - 1);
    // Normalize m to [3/4, 3/2) so near 1 does not cancel a multiple
    // of ln(2). The comparisons and normalization use exact integers.
    if (magnitude << 1_u32) >= &unit * 3 {
        unit <<= 1_u32;
        exponent += 1;
    }
    let numerator = magnitude - &unit;
    let denominator = magnitude + unit;
    let mut scale = prec.clamp(1, MAX_BITS) + LOG_GUARD_BITS;
    if exponent == 0 && numerator.sign() != num_bigint::Sign::NoSign {
        // Near one, k = 0, so use relative working precision for ln(m).
        let extra = denominator.bits().saturating_sub(numerator.bits());
        scale += u32::try_from(extra).map_err(|e| e.to_string())?;
    }
    let mut bounds = if numerator.sign() == num_bigint::Sign::NoSign {
        Bounds {
            low: 0.into(),
            high: 0.into(),
            scale,
        }
    } else {
        twice_atanh(&numerator, &denominator, scale)?
    };
    if exponent != 0 {
        let (low, high) = consts::ln2_at_scale(scale)?;
        let (low, high) = if exponent < 0 {
            (high, low)
        } else {
            (low, high)
        };
        bounds.low += low * exponent;
        bounds.high += high * exponent;
    }
    Ok(bounds)
}

pub(super) fn ln(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    natural_bounds(x, prec)?.ball(prec)
}

pub(super) fn log2(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.sign() > 0 && x.mantissa().trailing_zeros() == Some(u64::from(x.precision() - 1)) {
        let exponent = x.exponent() + i64::from(x.precision() - 1);
        // Every supported dyadic exponent fits exactly in 20 bits.
        let value = Real::from_i64(exponent, 20)?;
        return Ball::from_bounds(value.clone(), value, prec);
    }
    let scale = prec.clamp(1, MAX_BITS) + LOG_GUARD_BITS;
    let (low, high) = consts::ln2_at_scale(scale)?;
    natural_bounds(x, prec)?.divide(&Bounds { low, high, scale }, prec)
}

pub(super) fn log10(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let numerator = natural_bounds(x, prec)?;
    let denominator = natural_bounds(&Real::from_i64(10, 4)?, prec)?;
    numerator.divide(&denominator, prec)
}
