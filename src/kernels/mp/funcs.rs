//! Elementary functions over the arbitrary-precision reals.
//!
//! Functions return a computed midpoint and radius as a [`Ball`]. The exact
//! dyadic and interval exponential use directed bounds throughout. The other
//! series estimates use fixed point with at least 32 guard bits:
//!
//! * `exp` scales by an exact power of two, bounds its Taylor series, and
//!   restores the range by outward squaring;
//! * `ln` writes `x = m 2^e` with `m` in `[3/4, 3/2)` and bounds
//!   `2 atanh((m - 1)/(m + 1))`;
//! * `sin` and `cos` share certified quarter-turn reduction and Taylor bounds;
//! * `atan` and `atan2` use exact ratios and outward quadrant corrections.
//!
//! The series use direct summation and the operator surface caps the working
//! precision.

use num_bigint::BigInt;
use num_integer::Integer;

use super::ball::Ball;
use super::consts;
use super::fixed::Fx;
use super::logarithms;
use super::real::MAX_BITS;
use super::real::MAX_EXP;
use super::real::Real;
use super::real::Round;
use super::trigonometry;

mod hyperbolic;
mod inverse_sine;
mod inverse_tangent;
mod powers;

/// Outward integer endpoints with a shared power-of-two exponent.
struct ScaledBounds {
    lower: BigInt,
    upper: BigInt,
    exponent: i64,
}

/// Bounds exp(r) for 0 <= r <= 1/2 using one positive fixed-point series.
fn exp_reduced(
    r: &Real,
    scale: u32,
) -> Result<(BigInt, BigInt), String> {
    let reduced = Fx::from_real_trunc(r, scale);
    let mut term = Fx::one(scale);
    let mut sum = term.clone();
    // Each step at least halves the integer term, so this loop has a fixed bound.
    for j in 1..=i64::from(scale) + 1 {
        // For positive integers the two floors have the same effect as one
        // floor after dividing by 2^scale * j.
        term = term.mul(&reduced).div_i64(j);
        if term.is_zero() {
            // Let q = 2^-scale. For the truncated input each term loses at
            // most q from flooring, and propagates at most half the previous
            // error. Thus every term error is at most 2q. The omitted tail
            // is at most 4q and truncating r contributes at most 2q because
            // exp(r) < 2 throughout [0, 1/2].
            let upper = &sum.v + (2 * j + 4);
            return Ok((sum.v, upper));
        }
        sum = sum.add(&term);
    }
    Err("exponential series did not converge".to_string())
}

fn exp_bounds(
    x: &Real,
    prec: u32,
) -> Result<ScaledBounds, String> {
    let Some(top) = x.floor_log2() else {
        return Ok(ScaledBounds {
            lower: BigInt::from(1),
            upper: BigInt::from(1),
            exponent: 0,
        });
    };
    // e > 2. A magnitude beyond this threshold cannot be exponentiated
    // within the representable range at any supported precision.
    let limit = (MAX_EXP + i64::from(MAX_BITS)).ilog2() + 1;
    if top >= i64::from(limit) {
        return Err("exponential exponent out of range".to_string());
    }
    let squarings = u32::try_from((top + 2).max(0)).map_err(|e| e.to_string())?;
    let reduced = x.abs().scale_pow2(-i64::from(squarings))?;
    // Each squaring can double the relative interval width. Keep those
    // guard bits even when the requested precision is MAX_BITS.
    let scale = prec.clamp(1, MAX_BITS) + squarings + 32;
    let (mut lower, mut upper) = exp_reduced(&reduced, scale)?;
    if x.sign() < 0 {
        let numerator = BigInt::from(1) << (2 * scale);
        (lower, upper) = (&numerator / upper, numerator.div_ceil(&lower));
    }

    let mut exponent = -i64::from(scale);
    for _ in 0..squarings {
        let lower_squared = &lower * &lower;
        let upper_squared = &upper * &upper;
        // Normalize the upper endpoint after each square to keep the
        // integer sizes bounded. Both endpoints use the same power of two.
        let shift = upper_squared.bits().saturating_sub(u64::from(scale) + 1);
        lower = lower_squared >> shift;
        upper = ((upper_squared - 1_i32) >> shift) + 1;
        exponent = 2 * exponent + i64::try_from(shift).map_err(|e| e.to_string())?;
    }
    Ok(ScaledBounds {
        lower,
        upper,
        exponent,
    })
}

fn exp_enclosure(
    bounds: ScaledBounds,
    prec: u32,
) -> Result<Ball, String> {
    let ScaledBounds {
        lower,
        upper,
        exponent,
    } = bounds;
    Ball::from_scaled_bounds(lower, upper, exponent, prec)
}

/// An outward dyadic enclosure of the exponential of an exact real input.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn exp(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    exp_enclosure(exp_bounds(x, prec)?, prec)
}

/// An outward dyadic enclosure of the exponential of a closed interval.
///
/// # Errors
/// Fails when a required endpoint or result exponent is out of range.
pub fn exp_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    if value.rad.is_zero() {
        return exp(&value.mid, prec);
    }
    let working = prec.max(value.mid.precision()).max(value.rad.precision());
    let (left, right) = value.endpoints(working)?;
    let left = exp_bounds(&left, prec)?;
    let right = exp_bounds(&right, prec)?;
    let exponent = left.exponent.min(right.exponent);
    exp_enclosure(
        ScaledBounds {
            lower: left.lower << (left.exponent - exponent).unsigned_abs(),
            upper: right.upper << (right.exponent - exponent).unsigned_abs(),
            exponent,
        },
        prec,
    )
}

/// A certified outward enclosure of the natural logarithm of `x`.
///
/// # Errors
/// Fails for a non-positive operand or when a result exponent is out of range.
pub fn ln(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    logarithms::ln(x, prec)
}

/// A certified natural logarithm enclosure over a positive closed interval.
///
/// # Errors
/// Fails when the interval includes zero or negatives, or on exponent overflow.
pub fn ln_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    map(value, &ln, true, prec)
}

/// A certified outward enclosure of `sin(x)`.
///
/// # Errors
/// Fails for an argument magnitude greater than 2^20 or an unsupported exponent.
pub fn sin(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::point(x, prec, true)
}

/// A certified outward enclosure of `cos(x)`.
///
/// # Errors
/// Fails for an argument magnitude greater than 2^20 or an unsupported exponent.
pub fn cos(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::point(x, prec, false)
}

/// Outward sine and cosine enclosures sharing one argument reduction.
#[derive(Clone, Debug)]
pub struct SinCos {
    /// The sine enclosure.
    pub sin: Ball,
    /// The cosine enclosure.
    pub cos: Ball,
}

/// Computes both sine and cosine with one certified argument reduction.
///
/// # Errors
/// Fails for an argument magnitude greater than 2^20 or an unsupported exponent.
pub fn sin_cos(
    x: &Real,
    prec: u32,
) -> Result<SinCos, String> {
    let (sin, cos) = trigonometry::pair(x, prec)?;
    Ok(SinCos { sin, cos })
}

/// Fused sine and cosine enclosures of a closed interval.
pub(super) fn sin_cos_ball(
    value: &Ball,
    prec: u32,
) -> Result<SinCos, String> {
    let (sin, cos) = trigonometry::interval_pair(value, prec)?;
    Ok(SinCos { sin, cos })
}

/// A certified sine enclosure over a closed interval.
///
/// # Errors
/// Fails if the midpoint magnitude exceeds 2^20 or an exponent is unsupported.
pub fn sin_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::interval(value, prec, true)
}

/// A certified cosine enclosure over a closed interval.
///
/// # Errors
/// Fails if the midpoint magnitude exceeds 2^20 or an exponent is unsupported.
pub fn cos_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::interval(value, prec, false)
}

/// A certified outward enclosure of `tan(x)`.
///
/// # Errors
/// Fails if the cosine enclosure intersects zero, the argument magnitude
/// exceeds 2^20, or a result exponent is unsupported.
pub fn tan(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::tan(x, prec)
}

/// A certified tangent enclosure over a closed interval.
///
/// # Errors
/// Fails if the cosine enclosure intersects zero, the midpoint magnitude
/// exceeds 2^20, or a result exponent is unsupported.
pub fn tan_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    trigonometry::tan_interval(value, prec)
}

/// A certified outward enclosure of the real inverse tangent.
///
/// # Errors
/// Fails when a result exponent is outside the representable range.
pub fn atan(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    inverse_tangent::point(x, &Real::one(1), prec)
}

/// A certified inverse tangent enclosure over a closed dyadic interval.
///
/// # Errors
/// Fails when a required endpoint or result exponent is outside the range.
pub fn atan_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    inverse_tangent::rectangle(x, &Ball::point(Real::one(1)), prec)
}

/// A certified inverse sine enclosure of an exact dyadic.
///
/// # Errors
/// Fails for inputs outside [-1, 1] or when a result exponent is out of range.
pub fn asin(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    inverse_sine::Function::Asin.point(x, prec)
}

/// A certified inverse sine enclosure over a closed real interval.
///
/// # Errors
/// Fails if any part of the interval is outside [-1, 1] or a required exponent
/// is out of range.
pub fn asin_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    inverse_sine::Function::Asin.interval(x, prec)
}

/// A certified inverse cosine enclosure of an exact dyadic.
///
/// # Errors
/// Fails for inputs outside [-1, 1] or when a result exponent is out of range.
pub fn acos(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    inverse_sine::Function::Acos.point(x, prec)
}

/// A certified inverse cosine enclosure over a closed real interval.
///
/// # Errors
/// Fails if any part of the interval is outside [-1, 1] or a required exponent
/// is out of range.
pub fn acos_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    inverse_sine::Function::Acos.interval(x, prec)
}

/// An outward enclosure of the hyperbolic sine of an exact dyadic.
///
/// # Errors
/// Fails when a result exponent is outside the representable range.
pub fn sinh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Sinh.point(x, prec)
}

/// An outward hyperbolic sine enclosure over a closed interval.
///
/// # Errors
/// Fails when an endpoint or result exponent is outside the representable range.
pub fn sinh_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Sinh.interval(x, prec)
}

/// An outward enclosure of the hyperbolic cosine of an exact dyadic.
///
/// # Errors
/// Fails when a result exponent is outside the representable range.
pub fn cosh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Cosh.point(x, prec)
}

/// An outward hyperbolic cosine enclosure over a closed interval.
///
/// # Errors
/// Fails when an endpoint or result exponent is outside the representable range.
pub fn cosh_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Cosh.interval(x, prec)
}

/// An outward enclosure of the hyperbolic tangent of an exact dyadic.
///
/// # Errors
/// Fails when a result exponent is outside the representable range.
pub fn tanh(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Tanh.point(x, prec)
}

/// An outward hyperbolic tangent enclosure over a closed interval.
///
/// # Errors
/// Fails when an endpoint or result exponent is outside the representable range.
pub fn tanh_ball(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    hyperbolic::Function::Tanh.interval(x, prec)
}

/// A certified real power of exact dyadic operands.
///
/// # Errors
/// Fails for negative bases with non-integer exponents, zero with a negative
/// exponent, integer exponents outside signed 64-bit range, or exponent overflow.
pub fn pow(
    x: &Real,
    y: &Real,
    prec: u32,
) -> Result<Ball, String> {
    powers::point(x, y, prec)
}

/// A certified power over closed real base and exponent intervals.
///
/// # Errors
/// Non-integer or varying exponents require a strictly positive base interval.
/// Integral point exponents use the domain of [`powi_ball`].
pub fn pow_ball(
    x: &Ball,
    y: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    powers::interval(x, y, prec)
}

/// An outward power enclosure for an exact signed integer exponent.
///
/// # Errors
/// Fails for zero with a negative exponent or when a result exponent is outside
/// the representable range. Exponent zero returns one.
pub fn powi(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<Ball, String> {
    powers::integer(x, n, prec)
}

/// An outward signed integer power enclosure over a closed real interval.
///
/// # Errors
/// Fails when a negative exponent meets zero or an endpoint exponent is outside
/// the representable range.
pub fn powi_ball(
    x: &Ball,
    n: i64,
    prec: u32,
) -> Result<Ball, String> {
    powers::integer_interval(x, n, prec)
}

/// An outward power enclosure for a reduced exact rational exponent.
///
/// # Errors
/// Fails for a zero denominator, zero base with a negative exponent, negative
/// base with an even reduced denominator, or exponent overflow.
pub fn pow_ratio(
    x: &Real,
    numerator: i64,
    denominator: u32,
    prec: u32,
) -> Result<Ball, String> {
    powers::Rational::reduce(numerator, denominator)?.point(x, prec)
}

/// An outward rational power enclosure over a closed real interval.
///
/// # Errors
/// Fails for a zero denominator, negative input with an even reduced denominator,
/// zero input with a negative exponent, or an exponent outside the range.
pub fn pow_ratio_ball(
    x: &Ball,
    numerator: i64,
    denominator: u32,
    prec: u32,
) -> Result<Ball, String> {
    powers::Rational::reduce(numerator, denominator)?.interval(x, prec)
}

/// An outward dyadic enclosure of the real cube root.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn cbrt(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let (mid, inexact) = x.cbrt_with_inexact(prec, Round::NearestEven)?;
    let rad = super::ball::midpoint_rounding(&mid, inexact)?;
    Ball::new(mid, rad)
}

/// An outward dyadic enclosure of `sqrt(a^2 + b^2)`.
///
/// # Errors
/// Fails when the result exponent leaves the representable range.
pub fn hypot(
    a: &Real,
    b: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let (mid, inexact) = a.hypot_with_inexact(b, prec, Round::NearestEven)?;
    let rad = super::ball::midpoint_rounding(&mid, inexact)?;
    Ball::new(mid, rad)
}

/// A certified outward enclosure of the binary logarithm of `x`.
///
/// # Errors
/// Fails for a non-positive operand or when a result exponent is out of range.
pub fn log2(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    logarithms::log2(x, prec)
}

/// A certified binary logarithm enclosure over a positive closed interval.
///
/// # Errors
/// Fails when the interval includes zero or negatives, or on exponent overflow.
pub fn log2_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    map(value, &log2, true, prec)
}

/// A certified outward enclosure of the decimal logarithm of `x`.
///
/// # Errors
/// Fails for a non-positive operand or when a result exponent is out of range.
pub fn log10(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    logarithms::log10(x, prec)
}

/// A certified decimal logarithm enclosure over a positive closed interval.
///
/// # Errors
/// Fails when the interval includes zero or negatives, or on exponent overflow.
pub fn log10_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    map(value, &log10, true, prec)
}

/// The certified principal angle of `(x, y)` in `(-pi, pi]`.
///
/// # Errors
/// Fails at the origin or when a result exponent is outside the range.
pub fn atan2(
    y: &Real,
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    inverse_tangent::point(y, x, prec)
}

/// The principal-angle enclosure of a closed rectangle with coordinates
/// `x` and `y`. A rectangle lying wholly on the negative real axis has
/// angle +pi under the point convention.
///
/// # Errors
/// Fails if the rectangle contains the origin or intersects the negative
/// real axis without lying wholly on it, or an exponent is out of range.
pub fn atan2_ball(
    y: &Ball,
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    inverse_tangent::rectangle(y, x, prec)
}

/// Encloses the image of a closed interval under a monotone function.
///
/// The caller supplies both the monotonic direction and outward point
/// enclosures for the function.
///
/// # Errors
/// Fails when a point enclosure fails, an interval is invalid, or an
/// exponent leaves the representable range.
pub fn map(
    value: &Ball,
    function: &dyn Fn(&Real, u32) -> Result<Ball, String>,
    increasing: bool,
    prec: u32,
) -> Result<Ball, String> {
    if value.rad.is_zero() {
        let mapped = function(&value.mid, prec)?;
        return Ball::new(mapped.mid, mapped.rad);
    }
    let working = prec
        .saturating_add(40)
        .max(value.mid.precision())
        .max(value.rad.precision())
        .min(MAX_BITS);
    let (lower, upper) = value.outward_endpoints(working)?;
    let (first, last) = if increasing {
        (&lower, &upper)
    } else {
        (&upper, &lower)
    };
    let first = function(first, working)?;
    let last = function(last, working)?;
    if first.rad.sign() < 0 || last.rad.sign() < 0 {
        return Err("a function enclosure radius cannot be negative".into());
    }
    let lower = first
        .mid
        .sub(&first.rad, working, Round::TowardNegative)
        .or_else(|_| first.mid.sub(&first.rad, 1, Round::TowardNegative))?;
    let upper = last
        .mid
        .add(&last.rad, working, Round::TowardPositive)
        .or_else(|_| last.mid.add(&last.rad, 1, Round::TowardPositive))?;
    Ball::from_bounds(lower, upper, prec)
}

pub(super) fn hypot_ball(
    a: &Ball,
    b: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    if a.rad.is_zero() && b.rad.is_zero() {
        return hypot(&a.mid, &b.mid, prec);
    }
    let working = prec.saturating_add(16).min(MAX_BITS);
    let (a_near, a_far) = a.absolute_endpoints(working)?;
    let (b_near, b_far) = b.absolute_endpoints(working)?;
    let lower = a_near
        .hypot_with_inexact(&b_near, working, Round::TowardNegative)?
        .0;
    let upper = a_far
        .hypot_with_inexact(&b_far, working, Round::TowardPositive)?
        .0;
    Ball::from_bounds(lower, upper, prec)
}

/// A certified outward enclosure of `2^x`.
///
/// # Errors
/// Fails when a result exponent is outside the representable range.
pub fn exp2(
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let Some(top) = x.floor_log2() else {
        return Ok(Ball::point(Real::one(prec)));
    };
    if top >= i64::from((MAX_EXP + i64::from(MAX_BITS)).ilog2() + 1) {
        return Err("binary exponential exponent out of range".into());
    }
    if x.is_integer() {
        let exponent = i64::try_from(x.trunc_to_bigint()?).map_err(|e| e.to_string())?;
        return Real::rounded(1.into(), exponent, prec, Round::NearestEven).map(Ball::point);
    }
    if top < -i64::from(prec.clamp(1, MAX_BITS) + 4) {
        // The derivative of 2^t has magnitude at most 2 on [-1, 1].
        let radius = Real::rounded(1.into(), top + 2, 1, Round::TowardPositive)?;
        return Ball::new(Real::one(prec), radius);
    }
    let inner = prec.saturating_add(32);
    let ln2 = consts::named_ball("ln2", inner)?;
    let argument = Ball::point(x.clone()).mul(&ln2, inner)?;
    exp_ball(&argument, prec)
}

/// A certified binary exponential enclosure over a closed interval.
///
/// # Errors
/// Fails when a required endpoint or result exponent is out of range.
pub fn exp2_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    map(value, &exp2, true, prec)
}

/// `x * 10^n` rounded once to the target precision with ties to even.
///
/// # Errors
/// Fails when an exponent is outside the representable range.
pub fn scale_decimal(
    x: &Real,
    n: i64,
    prec: u32,
) -> Result<Real, String> {
    scale_decimal_in(x, n, prec, Round::NearestEven)
}

/// `x * 10^n` rounded once in the requested direction.
///
/// # Errors
/// Fails when an exponent is outside the representable range.
pub fn scale_decimal_in(
    x: &Real,
    n: i64,
    prec: u32,
    mode: Round,
) -> Result<Real, String> {
    if x.is_zero() {
        return Ok(Real::zero(prec));
    }
    if n == 0 {
        return x.with_precision(prec, mode);
    }
    // 10 > 2^3. Beyond this scale no two nonzero representable dyadics
    // can differ by a factor as large as the requested power of ten.
    let scale_limit = (2 * MAX_EXP + i64::from(MAX_BITS)).unsigned_abs() / 3 + 1;
    if n.unsigned_abs() > scale_limit {
        return Err("decimal scaling exponent out of range".into());
    }
    let five = super::decimal::power_of_five(n.unsigned_abs())?;
    let exponent = x.exponent() + n;
    if n >= 0 {
        Real::rounded(x.mantissa() * BigInt::from(five), exponent, prec, mode)
    } else {
        Real::rounded_ratio(x.mantissa(), &five, exponent, prec, mode)
    }
}
