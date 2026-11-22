//! Elementary functions over the arbitrary-precision reals.
//!
//! Every function returns a [`Ball`], not a bare [`Real`]: the midpoint is the
//! computed value and the radius bounds the truncation of the series plus the
//! rounding of each division. The series are the classical ones, evaluated in
//! fixed point at `prec + 32` fractional bits (more when an argument has a
//! large exponent, because the reduction loses that many bits):
//!
//! * `exp` reduces by `k ln 2` and sums `r^j / j!`;
//! * `ln` writes `x = m 2^e` with `m` in `[1, 2)` and sums `2 atanh(z)`,
//!   `z = (m - 1)/(m + 1) <= 1/3`;
//! * `sin` and `cos` reduce modulo `2 pi` and sum their Taylor series;
//! * `atan` halves the angle with `x / (1 + sqrt(1 + x^2))` until the Taylor
//!   series converges quickly.
//!
//! In every series the multiplier has magnitude below one, so an error of `k`
//! units of `2^-s` stays below `k` units: the error count is simply three
//! units per term plus the first omitted term. The cost is quadratic in the
//! precision (no binary splitting), which is why the operator surface caps the
//! digit count: a few thousand digits take about a second, twenty take
//! microseconds.

use num_bigint::BigInt;
use num_traits::Signed;

use super::ball::Ball;
use super::consts;
use super::fixed::Fx;
use super::fixed::units;
use super::real::Real;
use super::real::Round;

/// Largest magnitude whose sine or cosine is reduced exactly enough.
const MAX_TRIG_ARGUMENT: i64 = 1 << 20;

/// Turns a fixed-point result and its error count into a ball.
fn finish(
    value: Fx,
    error: u64,
    prec: u32,
) -> Result<Ball, String> {
    let mid = value.to_real(prec.saturating_add(16))?;
    let rad = units(error.saturating_add(4), value.s)?;
    Ball::new(mid.with_precision(prec, Round::NearestEven)?, rad)
}

/// The guard bits a reduction needs for an argument of magnitude `2^top`.
fn guard(
    prec: u32,
    top: i64,
) -> u32 {
    let extra = u32::try_from(top.max(0)).unwrap_or(u32::MAX).min(1 << 20);
    prec.saturating_add(32).saturating_add(extra)
}

/// `round(x / divisor)` as an integer, corrected so the remainder is small.
fn reduce(
    x: &Fx,
    divisor: &Fx,
) -> BigInt {
    if divisor.v == BigInt::from(0) {
        return BigInt::from(0);
    }
    let numerator = &x.v << usize::try_from(x.s).unwrap_or(0);
    let denominator = &divisor.v << usize::try_from(divisor.s).unwrap_or(0);
    let (quotient, remainder) = num_integer::Integer::div_rem(&numerator, &denominator);
    let twice = &remainder.abs() << 1_u32;
    let magnitude = denominator.abs();
    if twice >= magnitude {
        if x.v.sign() == num_bigint::Sign::Minus {
            quotient - BigInt::from(1)
        } else {
            quotient + BigInt::from(1)
        }
    } else {
        quotient
    }
}

/// `exp(x)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn exp(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::one(prec)));
    }
    let top = x.floor_log2().unwrap_or(0);
    let s = guard(prec, top);
    let ln2 = consts::ln2(s.saturating_add(64))?;
    let x_fx = Fx::from_real_trunc(x, s);
    let ln2_fx = Fx::from_real_trunc(&ln2, s);
    let k = reduce(&x_fx, &ln2_fx);
    let Ok(shift_count) = i64::try_from(k.clone()) else {
        return Err("exponent out of range".to_string());
    };
    let scaled_ln2 = ln2_fx.mul(&Fx::from_i64(shift_count, s));
    let reduced = x_fx.sub(&scaled_ln2);
    let mut term = Fx::one(s);
    let mut sum = term.clone();
    let mut j: i64 = 1;
    loop {
        term = term.mul(&reduced).div_i64(j);
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        j = j.saturating_add(1);
        if j > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    let scaled = sum.scale_pow2(i64::try_from(k).unwrap_or(i64::MAX));
    let count = u64::try_from(j).unwrap_or(u64::MAX);
    finish(Fx { v: scaled.v, s }, 3 * count + 8, prec)
}

/// `ln(x)` as a ball, for `x > 0`.
///
/// # Errors
/// Fails for a non-positive operand or when an exponent leaves the range.
pub fn ln(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.cmp_value(&Real::one(x.precision())) == 0 {
        return Ok(Ball::point(Real::zero(prec)));
    }
    if x.sign() <= 0 {
        return Err("logarithm of a non-positive number".to_string());
    }
    let Some(top) = x.floor_log2() else {
        return Err("logarithm of a non-positive number".to_string());
    };
    let s = guard(prec, top);
    // m in [1, 2): divide by 2^floor(log2 x), not by the mantissa's exponent.
    let mantissa = x.scale_pow2(-top)?;
    let m = Fx::from_real_trunc(&mantissa, s);
    let one = Fx::one(s);
    let z = m.sub(&one).div(&m.add(&one));
    let z2 = z.mul(&z);
    let mut power = z.clone();
    let mut sum = z.clone();
    let mut k: i64 = 1;
    loop {
        power = power.mul(&z2);
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        k = k.saturating_add(1);
        if k > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    let log_m = sum.scale_pow2(1);
    let ln2 = consts::ln2(s.saturating_add(64))?;
    let ln2_fx = Fx::from_real_trunc(&ln2, s);
    let from_exponent = ln2_fx.mul(&Fx::from_i64(top, s));
    let total = log_m.add(&from_exponent);
    let count = u64::try_from(k).unwrap_or(u64::MAX);
    finish(total, 3 * count + 8, prec)
}

/// `sin(x)` as a ball.
///
/// # Errors
/// Fails for an argument beyond `2^20` in magnitude, where the reduction
/// would need more bits than the exponent range allows, or when an exponent
/// leaves the range.
pub fn sin(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let (reduced, _s) = reduce_mod_two_pi(x, prec)?;
    let mut term = reduced.clone();
    let mut sum = reduced.clone();
    let mut j: i64 = 1;
    let square = reduced.mul(&reduced);
    loop {
        term = term.mul(&square).div_i64(2 * j).div_i64(2 * j + 1).neg();
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        j = j.saturating_add(1);
        if j > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    let count = u64::try_from(j).unwrap_or(u64::MAX);
    finish(sum, 3 * count + 8, prec)
}

/// `cos(x)` as a ball.
///
/// # Errors
/// Fails for an argument beyond `2^20` in magnitude or when an exponent
/// leaves the range.
pub fn cos(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::one(prec)));
    }
    let (reduced, s) = reduce_mod_two_pi(x, prec)?;
    let square = reduced.mul(&reduced);
    let mut term = Fx::one(s);
    let mut sum = term.clone();
    let mut j: i64 = 1;
    loop {
        term = term.mul(&square).div_i64(2 * j - 1).div_i64(2 * j).neg();
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        j = j.saturating_add(1);
        if j > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    let count = u64::try_from(j).unwrap_or(u64::MAX);
    finish(sum, 3 * count + 8, prec)
}

/// Reduces `x` modulo `2 pi` into `(-pi, pi]`.
fn reduce_mod_two_pi(
    x: &Real,
    prec: u32,
) -> Result<(Fx, u32), String> {
    let top = x.floor_log2().unwrap_or(0);
    if top > MAX_TRIG_ARGUMENT {
        return Err("trigonometric argument is too large".to_string());
    }
    let s = guard(prec, top);
    let pi = consts::pi(s.saturating_add(64))?;
    let two_pi = pi.scale_pow2(1)?;
    let x_fx = Fx::from_real_trunc(x, s);
    let period = Fx::from_real_trunc(&two_pi, s);
    let k = reduce(&x_fx, &period);
    let reduced = x_fx.sub(&Fx::from_i64(1, s).mul(&period).mul(&Fx { v: k, s }));
    Ok((reduced, s))
}

/// `tan(x)` as a ball.
///
/// # Errors
/// Fails where the cosine interval contains zero, or for an argument beyond
/// `2^20` in magnitude.
pub fn tan(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let sine = sin(x, prec.saturating_add(16))?;
    let cosine = cos(x, prec.saturating_add(16))?;
    sine.div(&cosine, prec)
}

/// `atan(x)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn atan(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let top = x.floor_log2().unwrap_or(0);
    let s = guard(prec, top);
    let (value, count) = atan_fixed(x, s)?;
    finish(value, count, prec)
}

/// `atan` in fixed point at scale `s`, with its error count.
fn atan_fixed(
    x: &Real,
    s: u32,
) -> Result<(Fx, u64), String> {
    let negative = x.sign() < 0;
    let magnitude = x.abs();
    let one = Fx::one(s);
    let mut value = Fx::from_real_trunc(&magnitude, s);
    let mut invert = false;
    if !value.abs_lt(&one) {
        invert = true;
        value = one.div(&value);
    }
    let mut halvings: u32 = 0;
    let threshold = Fx::from_i64(1, s).div_i64(16);
    let mut error: u64 = 4;
    while !value.abs_lt(&threshold) && halvings < 64 {
        let square = value.mul(&value);
        let root = one.add(&square).sqrt();
        value = value.div(&one.add(&root));
        halvings = halvings.saturating_add(1);
        error = error.saturating_add(4);
    }
    let square = value.mul(&value);
    let mut power = value.clone();
    let mut sum = value.clone();
    let mut j: i64 = 1;
    loop {
        power = power.mul(&square);
        let term = power.div_i64(2 * j + 1);
        if term.is_zero() {
            break;
        }
        sum = if j % 2 == 1 {
            sum.sub(&term)
        } else {
            sum.add(&term)
        };
        j = j.saturating_add(1);
        if j > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    error = error.saturating_add(3 * u64::try_from(j).unwrap_or(u64::MAX));
    let mut result = sum.scale_pow2(i64::from(halvings));
    if invert {
        let half_pi = consts::pi(s.saturating_add(64))?.scale_pow2(-1)?;
        let half_pi_fx = Fx::from_real_trunc(&half_pi, s);
        result = half_pi_fx.sub(&result);
        error = error.saturating_add(2);
    }
    if negative {
        result = result.neg();
    }
    Ok((result, error))
}

/// `asin(x)` as a ball, for `|x| <= 1`.
///
/// # Errors
/// Fails outside `[-1, 1]` or when an exponent leaves the range.
pub fn asin(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let one = Real::one(prec.saturating_add(32));
    if x.abs().cmp_value(&one) > 0 {
        return Err("arcsine of a number outside [-1, 1]".to_string());
    }
    if x.abs().cmp_value(&one) == 0 {
        let half_pi = consts::pi(prec.saturating_add(32))?.scale_pow2(-1)?;
        let value = if x.sign() < 0 {
            half_pi.neg()
        } else {
            half_pi
        };
        return Ok(Ball::point(value.with_precision(prec, Round::NearestEven)?));
    }
    let inner = prec.saturating_add(32);
    let square = x.mul(x, inner, Round::NearestEven)?;
    let rest = one.sub(&square, inner, Round::NearestEven)?;
    let root = rest.sqrt(inner, Round::NearestEven)?;
    let ratio = x.div(&root, inner, Round::NearestEven)?;
    atan(&ratio, prec)
}

/// `acos(x)` as a ball, for `|x| <= 1`.
///
/// # Errors
/// Fails outside `[-1, 1]` or when an exponent leaves the range.
pub fn acos(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(32);
    let half_pi = consts::pi(inner)?.scale_pow2(-1)?;
    let arc = asin(x, inner)?;
    let mid = half_pi.sub(&arc.mid, inner, Round::NearestEven)?;
    Ball::new(
        mid.with_precision(prec, Round::NearestEven)?,
        arc.rad.add(&units(4, inner)?, prec, Round::AwayFromZero)?,
    )
}

/// `sinh(x)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn sinh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let inner = prec.saturating_add(16);
    let positive = exp(x, inner)?;
    let negative = exp(&x.neg(), inner)?;
    positive
        .sub(&negative, inner)?
        .div(&Ball::point(Real::from_i64(2, inner)?), prec)
}

/// `cosh(x)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn cosh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::one(prec)));
    }
    let inner = prec.saturating_add(16);
    let positive = exp(x, inner)?;
    let negative = exp(&x.neg(), inner)?;
    positive
        .add(&negative, inner)?
        .div(&Ball::point(Real::from_i64(2, inner)?), prec)
}

/// `tanh(x)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn tanh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let inner = prec.saturating_add(16);
    let doubled = x.scale_pow2(1)?;
    let positive = exp(&doubled, inner)?;
    let one = Ball::point(Real::one(inner));
    let numerator = positive.sub(&one, inner)?;
    let denominator = positive.add(&one, inner)?;
    numerator.div(&denominator, prec)
}

/// `x^y` as a ball.
///
/// # Errors
/// Fails for a negative base with a non-integer exponent, or when an exponent
/// leaves the representable range.
pub fn pow(
    x: &Real,
    y: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if y.is_zero() {
        return Ok(Ball::point(Real::one(prec)));
    }
    if x.is_zero() {
        if y.sign() < 0 {
            return Err("zero to a negative power".to_string());
        }
        return Ok(Ball::point(Real::zero(prec)));
    }
    if y.is_integer() {
        let Some(exponent) = y.trunc_to_bigint().ok().and_then(|v| i64::try_from(v).ok()) else {
            return Err("integer exponent out of range".to_string());
        };
        if exponent.unsigned_abs() <= 1 << 20 {
            return integer_power(x, exponent, prec);
        }
    }
    if x.sign() < 0 {
        return Err("a negative base needs an integer exponent".to_string());
    }
    let inner = prec.saturating_add(32);
    let logarithm = ln(x, inner)?;
    let product = logarithm.mul(&Ball::point(y.clone()), inner)?;
    exp(&product.mid, inner).map(|ball| ball).and_then(|ball| {
        // The radius of the product widens the exponent, which widens the
        // result by exp(radius).
        let widening = exp(&product.rad, inner)?;
        let scaled = ball.mul(&widening, inner)?;
        Ok(Ball {
            mid: scaled.mid.with_precision(prec, Round::NearestEven)?,
            rad: scaled
                .rad
                .add(&units(4, inner)?, prec, Round::AwayFromZero)?,
        })
    })
}

/// `x^n` for an integer `n` by squaring.
fn integer_power(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let mut result = Real::one(inner);
    let mut base = if n < 0 {
        Real::one(inner).div(x, inner, Round::NearestEven)?
    } else {
        x.with_precision(inner, Round::NearestEven)?
    };
    let mut exponent = n.unsigned_abs();
    let mut error: u64 = 2;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = result.mul(&base, inner, Round::NearestEven)?;
            error = error.saturating_add(2);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = base.mul(&base, inner, Round::NearestEven)?;
            error = error.saturating_add(2);
        }
    }
    Ok(Ball {
        mid: result.with_precision(prec, Round::NearestEven)?,
        rad: units(error, inner)?,
    })
}

/// The cube root as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn cbrt(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() {
        return Ok(Ball::point(Real::zero(prec)));
    }
    let negative = x.sign() < 0;
    let target = x.abs();
    let inner = prec.saturating_add(32);
    let guess = target.to_f64().cbrt();
    let mut value = if guess.is_finite() && guess > 0.0 {
        Real::from_f64(guess, inner)?
    } else {
        Real::one(inner)
    };
    let three = Real::from_i64(3, inner)?;
    let two = Real::from_i64(2, inner)?;
    for _ in 0..200 {
        let square = value.mul(&value, inner, Round::NearestEven)?;
        let quotient = target.div(&square, inner, Round::NearestEven)?;
        let doubled = value.mul(&two, inner, Round::NearestEven)?;
        let next = doubled.add(&quotient, inner, Round::NearestEven)?.div(
            &three,
            inner,
            Round::NearestEven,
        )?;
        let delta = next.sub(&value, inner, Round::NearestEven)?.abs();
        value = next;
        if delta.sign() == 0 {
            break;
        }
        if let (Some(a), Some(b)) = (delta.floor_log2(), value.floor_log2())
            && b - a > i64::from(inner)
        {
            break;
        }
    }
    let mid = if negative {
        value.neg()
    } else {
        value
    };
    Ok(Ball {
        mid: mid.with_precision(prec, Round::NearestEven)?,
        rad: units(8, inner)?,
    })
}

/// `sqrt(a^2 + b^2)` as a ball.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn hypot(
    a: &Real,
    b: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(32);
    let left = a.mul(a, inner, Round::NearestEven)?;
    let right = b.mul(b, inner, Round::NearestEven)?;
    let sum = left.add(&right, inner, Round::NearestEven)?;
    let root = sum.sqrt(inner, Round::NearestEven)?;
    Ok(Ball {
        mid: root.with_precision(prec, Round::NearestEven)?,
        rad: units(6, inner)?,
    })
}

/// `log2(x)` as a ball, for `x > 0`.
///
/// # Errors
/// Fails for a non-positive operand or when an exponent leaves the range.
pub fn log2(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(32);
    let numerator = ln(x, inner)?;
    let denominator = Ball::point(consts::ln2(inner)?);
    numerator.div(&denominator, prec)
}

/// `log10(x)` as a ball, for `x > 0`.
///
/// # Errors
/// Fails for a non-positive operand or when an exponent leaves the range.
pub fn log10(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(32);
    let ten = Real::from_i64(10, inner)?;
    let numerator = ln(x, inner)?;
    let denominator = ln(&ten, inner)?;
    numerator.div(&denominator, prec)
}

/// `atan2(y, x)` as a ball, with the quadrant resolved.
///
/// # Errors
/// Fails when both operands are zero or an exponent leaves the range.
pub fn atan2(
    y: &Real,
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    if x.is_zero() && y.is_zero() {
        return Err("atan2 of two zeros".to_string());
    }
    let inner = prec.saturating_add(32);
    if x.sign() == 0 {
        let half_pi = consts::pi(inner)?.scale_pow2(-1)?;
        let value = if y.sign() < 0 {
            half_pi.neg()
        } else {
            half_pi
        };
        return Ok(Ball::point(value.with_precision(prec, Round::NearestEven)?));
    }
    let base = atan(&y.div(x, inner, Round::NearestEven)?, inner)?;
    if x.sign() > 0 {
        return Ok(Ball {
            mid: base.mid.with_precision(prec, Round::NearestEven)?,
            rad: base.rad.add(&units(4, inner)?, prec, Round::AwayFromZero)?,
        });
    }
    let pi = consts::pi(inner)?;
    let shifted = if y.sign() >= 0 {
        base.mid.add(&pi, inner, Round::NearestEven)?
    } else {
        base.mid.sub(&pi, inner, Round::NearestEven)?
    };
    Ok(Ball {
        mid: shifted.with_precision(prec, Round::NearestEven)?,
        rad: base.rad.add(&units(4, inner)?, prec, Round::AwayFromZero)?,
    })
}

/// The ball of a function applied to a ball, for a monotone function.
///
/// # Errors
/// Fails when the function fails on either endpoint.
pub fn map(
    value: &Ball,
    function: &dyn Fn(&Real, u32) -> Result<Ball, String>,
    increasing: bool,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(24);
    value
        .map_monotone(
            &|argument, working, _mode| function(argument, working).map(|ball| ball.mid),
            increasing,
            inner,
        )
        .and_then(|mapped| {
            let widening = function(&value.mid, inner)?;
            Ok(Ball {
                mid: mapped.mid.with_precision(prec, Round::NearestEven)?,
                rad: mapped.rad.add(&widening.rad, inner, Round::AwayFromZero)?,
            })
        })
}

/// `2^x`, a convenience over [`exp`].
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn exp2(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let scaled = x.mul(&consts::ln2(inner)?, inner, Round::NearestEven)?;
    exp(&scaled, prec)
}

/// `x * 10^n` for an integer `n`, exact up to the target precision.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn scale_decimal(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<Real, String> {
    if n == 0 {
        return x.with_precision(prec, Round::NearestEven);
    }
    let factor = if n >= 0 {
        BigInt::from(10_u32).pow(u32::try_from(n).unwrap_or(0))
    } else {
        BigInt::from(1)
    };
    if n >= 0 {
        Real::from_bigint(factor, prec, Round::NearestEven)?.mul(x, prec, Round::NearestEven)
    } else {
        let divisor = Real::from_bigint(
            BigInt::from(10_u32).pow(u32::try_from(-n).unwrap_or(0)),
            prec.saturating_add(8),
            Round::NearestEven,
        )?;
        x.div(&divisor, prec, Round::NearestEven)
    }
}
