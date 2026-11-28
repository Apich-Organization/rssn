//! Certified inverse sine and cosine using exact dyadics and the angle kernel.

use super::Ball;
use super::MAX_BITS;
use super::Real;
use super::ScaledBounds;
use super::inverse_tangent;

const GUARD_BITS: u32 = 40;

#[derive(Clone, Copy)]
pub(super) enum Function {
    Asin,
    Acos,
}

impl Function {
    pub(super) fn point(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<Ball, String> {
        enclose(self.bounds(x, prec)?, prec)
    }

    fn bounds(
        self,
        x: &Real,
        prec: u32,
    ) -> Result<ScaledBounds, String> {
        let working = prec
            .clamp(1, MAX_BITS)
            .saturating_add(GUARD_BITS)
            .min(MAX_BITS);
        let (root_low, root_high) = x.unit_complement_sqrt(working)?;
        // Varying the nonnegative root decreases atan2(x, r) when x
        // is positive and decreases atan2(r, x) when x is negative.
        let decreasing = match self {
            | Self::Asin => x.sign() > 0,
            | Self::Acos => x.sign() < 0,
        };
        let (low_argument, high_argument) = if decreasing {
            (root_high, root_low)
        } else {
            (root_low, root_high)
        };
        let angle = |root: &Real| {
            match self {
                | Self::Asin => inverse_tangent::atan2_bounds(x, root, prec),
                | Self::Acos => inverse_tangent::atan2_bounds(root, x, prec),
            }
        };
        let low = angle(&low_argument)?;
        if low_argument.cmp_value(&high_argument) == 0 {
            return Ok(low);
        }
        let high = angle(&high_argument)?;
        Ok(join(low, high))
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
        let (low, high) = x.endpoints(working)?;
        let (a, b) = if matches!(self, Self::Acos) {
            (high, low)
        } else {
            (low, high)
        };
        enclose(join(self.bounds(&a, prec)?, self.bounds(&b, prec)?), prec)
    }
}

fn join(
    lower: ScaledBounds,
    upper: ScaledBounds,
) -> ScaledBounds {
    let exponent = lower.exponent.min(upper.exponent);
    ScaledBounds {
        lower: lower.lower << (lower.exponent - exponent).unsigned_abs(),
        upper: upper.upper << (upper.exponent - exponent).unsigned_abs(),
        exponent,
    }
}

fn enclose(
    bounds: ScaledBounds,
    prec: u32,
) -> Result<Ball, String> {
    Ball::from_scaled_bounds(bounds.lower, bounds.upper, bounds.exponent, prec)
}
