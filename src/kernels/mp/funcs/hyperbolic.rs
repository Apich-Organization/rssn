//! Exact dyadic bounds for real hyperbolic functions.

use num_bigint::BigInt;
use num_integer::Integer;

use super::Ball;
use super::Fx;
use super::MAX_BITS;
use super::Real;
use super::Round;
use super::ScaledBounds;
use super::exp_bounds;

const GUARD_BITS: u32 = 48;

#[derive(Clone, Copy)]
pub(super) enum Function {
    Sinh,
    Cosh,
    Tanh,
}

impl Function {
    pub(super) fn point(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<Ball, String> {
        let bounds = self.bounds(x, prec)?;
        Ball::from_scaled_bounds(bounds.lower, bounds.upper, bounds.exponent, prec)
    }

    fn bounds(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<ScaledBounds, String> {
        if x.is_zero() {
            return Ok(exact_integer(i32::from(matches!(self, Self::Cosh))));
        }
        let magnitude = x.abs();
        let precision = prec.clamp(1, MAX_BITS);
        let half = Real::rounded(1.into(), -1, 1, Round::NearestEven)?;
        let bounds = if magnitude.cmp_value(&half) <= 0 {
            self.small(&magnitude, precision + GUARD_BITS)?
        } else {
            self.large(&magnitude, precision)?
        };
        if x.sign() < 0 && !matches!(self, Self::Cosh) {
            return Ok(ScaledBounds {
                lower: -bounds.upper,
                upper: -bounds.lower,
                exponent: bounds.exponent,
            });
        }
        Ok(bounds)
    }

    fn small(
        self,
        x: &Real,
        scale: u32,
    ) -> Result<ScaledBounds, String> {
        let shift = 2 * x.exponent() + i64::from(scale);
        let square = x.mantissa() * x.mantissa();
        let square = Fx {
            v: if shift >= 0 {
                square << shift.unsigned_abs()
            } else {
                square >> shift.unsigned_abs()
            },
            s: scale,
        };
        let pair = match self {
            | Self::Sinh => small_series(&square, true)?,
            | Self::Cosh => small_series(&square, false)?,
            | Self::Tanh => {
                let (sin_low, sin_high) = small_series(&square, true)?;
                let (cos_low, cos_high) = small_series(&square, false)?;
                (
                    (sin_low << scale) / cos_high,
                    (sin_high << scale).div_ceil(&cos_low),
                )
            },
        };
        let (low, high) = pair;
        if matches!(self, Self::Cosh) {
            return Ok(ScaledBounds {
                lower: low,
                upper: high,
                exponent: -i64::from(scale),
            });
        }
        // sinh(x)/x and tanh(x)/x stay near one. Keeping x exact avoids
        // cancellation and preserves relative accuracy at tiny arguments.
        Ok(ScaledBounds {
            lower: low * x.mantissa(),
            upper: high * x.mantissa(),
            exponent: x.exponent() - i64::from(scale),
        })
    }

    fn large(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<ScaledBounds, String> {
        if matches!(self, Self::Tanh) {
            return tanh_positive(x, prec);
        }
        let ScaledBounds {
            lower,
            upper,
            exponent,
        } = exp_bounds(x, prec)?;
        let (reciprocal_lower, reciprocal_upper) = if exponent > 0 {
            (BigInt::from(0), BigInt::from(1))
        } else {
            let unit = BigInt::from(1) << (2 * exponent).unsigned_abs();
            (&unit / &upper, unit.div_ceil(&lower))
        };
        let (lower, upper) = if matches!(self, Self::Cosh) {
            (lower + reciprocal_lower, upper + reciprocal_upper)
        } else {
            (lower - reciprocal_upper, upper - reciprocal_lower)
        };
        // Divide by two before converting to a Real. exp(x) can exceed
        // its representable exponent even when sinh(x) and cosh(x) fit.
        Ok(ScaledBounds {
            lower,
            upper,
            exponent: exponent - 1,
        })
    }

    pub(super) fn interval(
        self,
        x: &Ball,
        prec: u32,
    ) -> Result<Ball, String> {
        if x.rad.is_zero() {
            return self.point(&x.mid, prec);
        }
        let working = prec
            .saturating_add(GUARD_BITS)
            .clamp(1, MAX_BITS)
            .max(x.mid.precision())
            .max(x.rad.precision());
        let (input_low, input_high) = x.endpoints(working)?;
        let (lower, upper) = if matches!(self, Self::Cosh) {
            cosh_endpoints(&input_low, &input_high, working)?
        } else {
            (
                self.bounds(&input_low, working)?,
                self.bounds(&input_high, working)?,
            )
        };
        let exponent = lower.exponent.min(upper.exponent);
        Ball::from_scaled_bounds(
            lower.lower << (lower.exponent - exponent).unsigned_abs(),
            upper.upper << (upper.exponent - exponent).unsigned_abs(),
            exponent,
            prec,
        )
    }
}

fn exact_integer(value: i32) -> ScaledBounds {
    ScaledBounds {
        lower: value.into(),
        upper: value.into(),
        exponent: 0,
    }
}

/// Bounds the positive series for cosh(x) or sinh(x)/x at 0 <= x <= 1/2.
fn small_series(
    square: &Fx,
    sinh_quotient: bool,
) -> Result<(BigInt, BigInt), String> {
    let mut term = Fx::one(square.s);
    let mut total = term.clone();
    let offset = i64::from(sinh_quotient);
    for k in 1..=i64::from(square.s) + 1 {
        let divisor = (2 * k - 1 + offset) * (2 * k + offset);
        // In units of 2^-scale, truncating x² loses below one unit.
        // Each next term contracts prior error by at most 1/8 and adds
        // at most 1/2 unit for the square plus one for the combined floors.
        // Therefore each term loses at most two units, and the remaining
        // positive tail is below three units when the working term is zero.
        term = term.mul(square).div_i64(divisor);
        if term.is_zero() {
            let upper = &total.v + (2 * k + 1);
            return Ok((total.v, upper));
        }
        total = total.add(&term);
    }
    Err("hyperbolic series did not converge".into())
}

fn tanh_positive(
    x: &Real,
    prec: u32,
) -> Result<ScaledBounds, String> {
    let cutoff = prec + 32;
    let limit = Real::from_i64(i64::from(cutoff), 20)?;
    if x.cmp_value(&limit) >= 0 {
        // 1 - tanh(x) < 2 exp(-2x) < 2^(1 - 2x).
        let bits = 2 * cutoff - 1;
        let one = BigInt::from(1) << bits;
        return Ok(ScaledBounds {
            lower: &one - 1,
            upper: one,
            exponent: -i64::from(bits),
        });
    }
    let ScaledBounds {
        lower,
        upper,
        exponent,
    } = exp_bounds(x, prec)?;
    let squared_exponent = 2 * exponent;
    let common = squared_exponent.min(0);
    let shift = (squared_exponent - common).unsigned_abs();
    let unity = BigInt::from(1) << common.unsigned_abs();
    let lower = (lower.pow(2) << shift) + &unity;
    let upper = (upper.pow(2) << shift) + &unity;
    let scale = prec + GUARD_BITS;
    let numerator = unity << (scale + 1);
    let correction_lower = &numerator / upper;
    let correction_upper = numerator.div_ceil(&lower);
    let one = BigInt::from(1) << scale;
    // tanh(x) = 1 - 2/(exp(x)^2 + 1). Subtraction reverses the endpoints.
    Ok(ScaledBounds {
        lower: &one - correction_upper,
        upper: one - correction_lower,
        exponent: -i64::from(scale),
    })
}

fn cosh_endpoints(
    low: &Real,
    high: &Real,
    prec: u32,
) -> Result<(ScaledBounds, ScaledBounds), String> {
    if low.sign() >= 0 {
        return Ok((
            Function::Cosh.bounds(low, prec)?,
            Function::Cosh.bounds(high, prec)?,
        ));
    }
    if high.sign() <= 0 {
        return Ok((
            Function::Cosh.bounds(high, prec)?,
            Function::Cosh.bounds(low, prec)?,
        ));
    }
    let left = low.abs();
    let farthest = if left.cmp_value(high) >= 0 {
        &left
    } else {
        high
    };
    Ok((exact_integer(1), Function::Cosh.bounds(farthest, prec)?))
}
