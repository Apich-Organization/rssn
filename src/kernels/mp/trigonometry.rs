//! Certified quarter-turn reduction and alternating trigonometric series.

use num_bigint::BigInt;
use num_bigint::Sign;
use num_traits::Signed;

use super::ball::Ball;
use super::consts;
use super::fixed::Fx;
use super::real::MAX_BITS;
use super::real::Real;
use super::real::Round;

const MAX_ARGUMENT_TOP: i64 = 20;
const GUARD_BITS: u32 = 64;
const TAN_GUARD_BITS: u32 = 32;

struct Reduced {
    magnitude: Fx,
    square: Fx,
    uncertainty: BigInt,
    quadrant: i64,
    negative_input: bool,
    negative_remainder: bool,
}

impl Reduced {
    fn new(
        x: &Real,
        prec: u32,
    ) -> Result<Self, String> {
        let top = x.floor_log2().unwrap_or(0);
        if top > MAX_ARGUMENT_TOP
            || (top == MAX_ARGUMENT_TOP
                && x.mantissa().trailing_zeros() != Some(u64::from(x.precision() - 1)))
        {
            return Err("trigonometric argument magnitude exceeds 2^20".into());
        }
        let extra = u32::try_from(top.unsigned_abs()).map_err(|error| error.to_string())?;
        let scale = prec.clamp(1, MAX_BITS) + GUARD_BITS + extra;
        let (mut lower, mut upper) = scaled_magnitude(x, scale);
        let quadrant = if top >= -1 && !x.is_zero() {
            // Scaling pi by 2^(scale - 1) gives quarter-turn units of 2^-scale.
            let (quarter_low, quarter_high) = consts::pi_at_scale(scale - 1)?;
            let quarter = (&quarter_low + &quarter_high) >> 1_u32;
            let turns = i64::try_from((&lower + (&quarter >> 1_u32)) / &quarter)
                .map_err(|error| error.to_string())?;
            lower -= quarter_high * turns;
            upper -= quarter_low * turns;
            turns % 4
        } else {
            0
        };
        let unit = BigInt::from(1) << scale;
        if lower < -&unit || upper > unit {
            return Err("could not certify the trigonometric reduction".into());
        }
        let centre = (&lower + &upper) >> 1_u32;
        let uncertainty = (&centre - lower).max(upper - &centre);
        let negative_remainder = centre.is_negative();
        let magnitude = Fx {
            v: centre.abs(),
            s: scale,
        };
        let square = magnitude.mul(&magnitude);
        Ok(Self {
            magnitude,
            square,
            uncertainty,
            quadrant,
            negative_input: x.sign() < 0,
            negative_remainder,
        })
    }

    fn component(
        &self,
        sine: bool,
        prec: u32,
    ) -> Result<Ball, String> {
        let even_quadrant = self.quadrant % 2 == 0;
        let use_sine = sine == even_quadrant;
        let negate = if sine {
            (self.quadrant >= 2) ^ self.negative_input
        } else {
            self.quadrant == 1 || self.quadrant == 2
        } ^ (use_sine && self.negative_remainder);
        let (mut value, series_error) = self.series(use_sine)?;
        if negate {
            value = -value;
        }
        // Both functions have absolute derivative at most one.
        let error = &self.uncertainty + series_error;
        Ball::from_scaled_bounds(
            &value - &error,
            value + error,
            -i64::from(self.magnitude.s),
            prec,
        )
    }

    fn series(
        &self,
        sine: bool,
    ) -> Result<(BigInt, i64), String> {
        if self.magnitude.is_zero() {
            let value = if sine {
                BigInt::from(0)
            } else {
                Fx::one(self.magnitude.s).v
            };
            return Ok((value, 0));
        }
        let mut term = if sine {
            self.magnitude.clone()
        } else {
            Fx::one(self.magnitude.s)
        };
        let mut sum = term.clone();
        let offset = i64::from(sine);
        for k in 1..=i64::from(self.magnitude.s) + 1 {
            let denominator = (2 * k - 1 + offset) * (2 * k + offset);
            // Two positive integer floors combine into one. In units of
            // 2^-scale the square error is below one, and the next term's
            // error is bounded by half the preceding error plus 1.5.
            // Thus each term and the first omitted term have error ≤ 3.
            term = term.mul(&self.square).div_i64(denominator);
            if term.is_zero() {
                return Ok((sum.v, 3 * k));
            }
            sum = if k % 2 == 1 {
                sum.sub(&term)
            } else {
                sum.add(&term)
            };
        }
        Err("trigonometric series did not converge".into())
    }
}

fn scaled_magnitude(
    x: &Real,
    scale: u32,
) -> (BigInt, BigInt) {
    let magnitude = x.mantissa().magnitude();
    let shift = x.exponent() + i64::from(scale);
    if shift >= 0 {
        let exact = BigInt::from(magnitude << shift.unsigned_abs());
        return (exact.clone(), exact);
    }
    let dropped = shift.unsigned_abs();
    let lower = BigInt::from_biguint(Sign::Plus, magnitude >> dropped);
    let lost = magnitude
        .trailing_zeros()
        .is_some_and(|zeros| zeros < dropped);
    let upper = &lower + i32::from(lost);
    (lower, upper)
}

pub(super) fn point(
    x: &Real,
    prec: u32,
    sine: bool,
) -> Result<Ball, String> {
    Reduced::new(x, prec)?.component(sine, prec)
}

pub(super) fn pair(
    x: &Real,
    prec: u32,
) -> Result<(Ball, Ball), String> {
    let reduced = Reduced::new(x, prec)?;
    Ok((
        reduced.component(true, prec)?,
        reduced.component(false, prec)?,
    ))
}

fn widen(
    mut value: Ball,
    radius: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if !radius.is_zero() {
        let rad = value
            .rad
            .add(radius, prec, Round::TowardPositive)
            .or_else(|_| value.rad.add(radius, 1, Round::TowardPositive))?;
        // The whole mathematical image always lies in [-1, 1].
        if rad.cmp_value(&Real::from_i64(2, 1)?) >= 0 {
            return Ball::new(Real::zero(prec), Real::one(prec));
        }
        value.rad = rad;
    }
    Ok(value)
}

pub(super) fn interval(
    value: &Ball,
    prec: u32,
    sine: bool,
) -> Result<Ball, String> {
    widen(point(&value.mid, prec, sine)?, &value.rad, prec)
}

pub(super) fn tan(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let working = prec.saturating_add(TAN_GUARD_BITS).min(MAX_BITS);
    let (sine, cosine) = pair(x, working)?;
    sine.div(&cosine, prec)
}

pub(super) fn interval_pair(
    value: &Ball,
    prec: u32,
) -> Result<(Ball, Ball), String> {
    let (sine, cosine) = pair(&value.mid, prec)?;
    Ok((
        widen(sine, &value.rad, prec)?,
        widen(cosine, &value.rad, prec)?,
    ))
}

pub(super) fn tan_interval(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    let working = prec.saturating_add(TAN_GUARD_BITS).min(MAX_BITS);
    let (sine, cosine) = interval_pair(value, working)?;
    sine.div(&cosine, prec)
}
