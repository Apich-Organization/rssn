//! Principal inverse tangent from certified quadrant and exact ratio bounds.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::Zero;

use super::Ball;
use super::Fx;
use super::MAX_BITS;
use super::Real;
use super::ScaledBounds;
use super::consts;

const GUARD_BITS: u32 = 56;

fn native_ball(
    bounds: ScaledBounds,
    prec: u32,
) -> Result<Ball, String> {
    Ball::from_scaled_bounds(bounds.lower, bounds.upper, bounds.exponent, prec)
}

pub(super) fn point(
    y: &Real,
    x: &Real,
    prec: u32,
) -> Result<Ball, String> {
    native_ball(atan2_bounds(y, x, prec)?, prec)
}

pub(super) fn atan2_bounds(
    y: &Real,
    x: &Real,
    prec: u32,
) -> Result<ScaledBounds, String> {
    if x.is_zero() && y.is_zero() {
        return Err("atan2 of two zeros".into());
    }
    let (numerator, denominator, quarter_turns, direction) = if x.is_zero() {
        (BigUint::zero(), BigUint::from(1_u32), 2, 1)
    } else if y.is_zero() {
        (
            BigUint::zero(),
            BigUint::from(1_u32),
            if x.sign() < 0 { 4 } else { 0 },
            1,
        )
    } else {
        let abs_x = x.abs();
        let abs_y = y.abs();
        let side = x.sign();
        if abs_x.cmp_value(&abs_y) >= 0 {
            let (n, d) = ratio(&abs_y, &abs_x);
            (n, d, if side < 0 { 4 } else { 0 }, i32::from(side))
        } else {
            let (n, d) = ratio(&abs_x, &abs_y);
            (n, d, 2, -i32::from(side))
        }
    };
    let mut bounds = reduced(numerator, denominator, quarter_turns, direction, prec)?;
    if y.sign() < 0 {
        (bounds.lower, bounds.upper) = (-bounds.upper, -bounds.lower);
    }
    Ok(bounds)
}

fn ratio(
    numerator: &Real,
    denominator: &Real,
) -> (BigUint, BigUint) {
    let n = numerator.mantissa().magnitude();
    let d = denominator.mantissa().magnitude();
    let shift = numerator.exponent() - denominator.exponent();
    if shift >= 0 {
        (n << shift.unsigned_abs(), d.clone())
    } else {
        (n.clone(), d << shift.unsigned_abs())
    }
}

fn reduced(
    numerator: BigUint,
    denominator: BigUint,
    mut quarter_turns: i32,
    mut direction: i32,
    prec: u32,
) -> Result<ScaledBounds, String> {
    let (numerator, denominator) = if (&numerator << 1_u32) >= denominator {
        // atan(t) = pi/4 - atan((1 - t)/(1 + t)), for 1/2 <= t <= 1.
        quarter_turns += direction;
        direction = -direction;
        (&denominator - &numerator, denominator + numerator)
    } else {
        (numerator, denominator)
    };
    let mut scale = prec.clamp(1, MAX_BITS) + GUARD_BITS;
    if quarter_turns == 0 && !numerator.is_zero() {
        // Preserve relative precision near the positive real axis. No
        // multiple of pi is needed on this branch.
        let gap = denominator.bits().saturating_sub(numerator.bits());
        scale += u32::try_from(gap).map_err(|error| error.to_string())?;
    }
    let (mut lower, mut upper) = series(&numerator, &denominator, scale)?;
    if direction < 0 {
        (lower, upper) = (-upper, -lower);
    }
    if quarter_turns != 0 {
        // pi/4 has integral bounds at this scale from the cached pi bounds.
        let (pi_low, pi_high) = consts::pi_at_scale(scale - 2)?;
        lower += pi_low * quarter_turns;
        upper += pi_high * quarter_turns;
    }
    Ok(ScaledBounds {
        lower,
        upper,
        exponent: -i64::from(scale),
    })
}

fn series(
    numerator: &BigUint,
    denominator: &BigUint,
    scale: u32,
) -> Result<(BigInt, BigInt), String> {
    if numerator.is_zero() {
        return Ok((0.into(), 0.into()));
    }
    let value = Fx {
        v: BigInt::from((numerator << scale) / denominator),
        s: scale,
    };
    let square = value.mul(&value);
    let mut power = value;
    let mut total = Fx::zero(scale);
    for k in 0..=i64::from(scale) {
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            // The exact ratio is <= 1/2. In fixed-point units each power
            // contracts the previous error by at most 1/4 and adds at
            // most 1.5. Term error is at most three units. The alternating
            // tail loses at most three more, and |atan'(x)| <= 1 accounts
            // for the initial truncation within one unit.
            let error = 3 * k + 4;
            let lower = (&total.v - error).max(BigInt::from(0));
            return Ok((lower, total.v + error));
        }
        total = if k % 2 == 0 {
            total.add(&term)
        } else {
            total.sub(&term)
        };
        power = power.mul(&square);
    }
    Err("inverse tangent series did not converge".into())
}

fn hull(
    first: ScaledBounds,
    next: ScaledBounds,
) -> ScaledBounds {
    let exponent = first.exponent.min(next.exponent);
    let a = (first.exponent - exponent).unsigned_abs();
    let b = (next.exponent - exponent).unsigned_abs();
    ScaledBounds {
        lower: (first.lower << a).min(next.lower << b),
        upper: (first.upper << a).max(next.upper << b),
        exponent,
    }
}

pub(super) fn rectangle(
    y: &Ball,
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    if y.contains_zero() && x.contains_zero() {
        return Err("atan2 rectangle contains the origin".into());
    }
    if y.rad.is_zero() && x.rad.is_zero() {
        return point(&y.mid, &x.mid, prec);
    }
    let working = prec
        .saturating_add(16)
        .clamp(1, MAX_BITS)
        .max(y.mid.precision())
        .max(y.rad.precision())
        .max(x.mid.precision())
        .max(x.rad.precision());
    let (y_low, y_high) = y.endpoints(working)?;
    let (x_low, x_high) = x.endpoints(working)?;
    if y.contains_zero() && x_low.sign() < 0 {
        if y.mid.is_zero() && y.rad.is_zero() {
            // A rectangle lying wholly on the negative axis has angle +pi.
            return point(&y.mid, &x.mid, prec);
        }
        return Err("atan2 rectangle intersects the negative real axis".into());
    }

    // On a rectangle separated from the origin and branch cut, the angle
    // is monotone or constant along each edge. Its gradient has no
    // interior zero, so the corner hull contains the image.
    let mut bounds = atan2_bounds(&y_low, &x_low, prec)?;
    let distinct_x = x_low.cmp_value(&x_high) != 0;
    let distinct_y = y_low.cmp_value(&y_high) != 0;
    if distinct_x {
        bounds = hull(bounds, atan2_bounds(&y_low, &x_high, prec)?);
    }
    if distinct_y {
        bounds = hull(bounds, atan2_bounds(&y_high, &x_low, prec)?);
        if distinct_x {
            bounds = hull(bounds, atan2_bounds(&y_high, &x_high, prec)?);
        }
    }
    native_ball(bounds, prec)
}
