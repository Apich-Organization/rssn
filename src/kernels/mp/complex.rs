//! Arbitrary-precision complex ball arithmetic (`MpComplex`) and
//! Ziv-certified complex expression evaluator.

use super::MAX_BITS;
use super::MAX_DIGITS;
use super::ball::Ball;
use super::consts;
use super::decimal::ball_to_decimal;
use super::decimal::bits_for_digits;
use super::decimal::from_decimal;
use super::funcs;
use super::real::Real;
use super::real::Round;
use num_bigint::BigInt;

/// Maximum expression length accepted by [`eval_complex_text`].
const MAX_TEXT: usize = 4096;
/// Maximum nesting depth accepted by [`eval_complex_text`].
const MAX_DEPTH: usize = 64;
/// Maximum tokens accepted by [`eval_complex_text`].
const MAX_TOKENS: usize = 2048;

/// An arbitrary-precision complex rectangle `re + i * im` where `re` and `im`
/// are rigorous [`Ball`] enclosures.
#[derive(Clone, Debug)]
pub struct MpComplex {
    /// Real part enclosure.
    pub re: Ball,
    /// Imaginary part enclosure.
    pub im: Ball,
}

/// Conservative 1-ulp bound for `value` at precision `prec`.
fn ulp_bound(
    value: &Real,
    prec: u32,
) -> Result<Real, String> {
    let exponent = match value.floor_log2() {
        | None => -i64::from(prec),
        | Some(top) => top.saturating_sub(i64::from(prec)).saturating_add(1),
    };
    Real::rounded(BigInt::from(1), exponent, 1, Round::AwayFromZero)
}

/// Halves an enclosure and rounds outward to the requested precision.
fn half_ball(
    value: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    let half = Real::one(1).scale_pow2(-1)?;
    value.mul(&Ball::point(half), prec)
}

fn half_pi_ball(prec: u32) -> Result<Ball, String> {
    let pi = consts::named_ball("pi", prec.saturating_add(16).min(MAX_BITS))?;
    half_ball(&pi, prec)
}

/// Shared factors in the circular and hyperbolic rectangle identities.
struct RectangleTrig {
    circular: funcs::SinCos,
    sinh: Ball,
    cosh: Ball,
}

impl RectangleTrig {
    fn new(
        angle: &Ball,
        real: &Ball,
        precision: u32,
    ) -> Result<Self, String> {
        let working = precision.saturating_add(16).min(MAX_BITS);
        Ok(Self {
            circular: funcs::sin_cos_ball(angle, working)?,
            sinh: funcs::sinh_ball(real, working)?,
            cosh: funcs::cosh_ball(real, working)?,
        })
    }

    fn sin(
        &self,
        precision: u32,
    ) -> Result<MpComplex, String> {
        Ok(MpComplex {
            re: self.circular.sin.mul(&self.cosh, precision)?,
            im: self.circular.cos.mul(&self.sinh, precision)?,
        })
    }

    fn cos(
        &self,
        precision: u32,
    ) -> Result<MpComplex, String> {
        Ok(MpComplex {
            re: self.circular.cos.mul(&self.cosh, precision)?,
            im: self.circular.sin.mul(&self.sinh, precision)?.neg(),
        })
    }

    fn sinh(
        &self,
        precision: u32,
    ) -> Result<MpComplex, String> {
        let sin = self.sin(precision)?;
        Ok(MpComplex {
            re: sin.im,
            im: sin.re,
        })
    }

    fn cosh(
        &self,
        precision: u32,
    ) -> Result<MpComplex, String> {
        let cos = self.cos(precision)?;
        Ok(MpComplex {
            re: cos.re,
            im: cos.im.neg(),
        })
    }
}

fn is_zero_ball(value: &Ball) -> bool {
    value.mid.is_zero() && value.rad.is_zero()
}

/// A directed bound for sqrt((hypot(x, y) + x) / 2) with x, y >= 0.
fn sqrt_large_component(
    x: &Real,
    y: &Real,
    prec: u32,
    mode: Round,
) -> Result<Real, String> {
    if y.is_zero() {
        return x.sqrt(prec, mode);
    }
    let (lower, upper, scale) = x.scaled_hypot_bounds(y, prec)?;
    let norm = if mode == Round::TowardNegative {
        lower
    } else {
        upper
    };
    norm.sqrt_half_sum_scaled(scale, x, prec, mode)
}

/// A directed bound for sqrt((hypot(x, y) + x) / 2) with y >= 0.
fn sqrt_component(
    x: &Real,
    y: &Real,
    prec: u32,
    mode: Round,
) -> Result<Real, String> {
    if x.sign() >= 0 {
        return sqrt_large_component(x, y, prec, mode);
    }
    if y.is_zero() {
        return Ok(Real::zero(prec));
    }
    // For negative x, use y / (2 sqrt((hypot(x, y) - x) / 2)).
    // This avoids cancellation and reverses the denominator's rounding.
    let denominator_round = if mode == Round::TowardNegative {
        Round::TowardPositive
    } else {
        Round::TowardNegative
    };
    let denominator = sqrt_large_component(&x.neg(), y, prec, denominator_round)?.scale_pow2(1)?;
    y.div(&denominator, prec, mode)
}

/// Squares an interval using its nearest and farthest points from zero.
fn squared_coordinate_bounds(
    value: &Ball,
    prec: u32,
) -> Result<(Real, Real), String> {
    let (near, far) = value.absolute_endpoints(prec)?;
    let low = near.mul(&near, prec, Round::TowardNegative)?;
    let high = far.mul(&far, prec, Round::TowardPositive)?;
    Ok((low, high))
}

/// Evaluates ln(hypot(x, y)) with the norm's binary scale kept separate.
fn log_hypot(
    x: &Real,
    y: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let (lower, upper, scale) = x.scaled_hypot_bounds(y, prec)?;
    let normalized = funcs::ln_ball(&Ball::from_bounds(lower, upper, prec)?, prec)?;
    if scale == 0 {
        return Ok(normalized);
    }
    let ln2 = consts::named_ball("ln2", prec)?;
    let shift = Ball::point(Real::from_i64(scale, 20)?);
    normalized.add(&ln2.mul(&shift, prec)?, prec)
}

/// Divides a ball by positive bounds without rounding the divisor to a ball.
fn divide_by_positive_bounds(
    numerator: &Ball,
    lower: &Real,
    upper: &Real,
    prec: u32,
) -> Result<Ball, String> {
    let (lo, hi) = numerator.endpoints(prec)?;
    let lo_divisor = if lo.sign() < 0 {
        lower
    } else {
        upper
    };
    let hi_divisor = if hi.sign() > 0 {
        lower
    } else {
        upper
    };
    let lo = lo.div(lo_divisor, prec, Round::TowardNegative)?;
    let hi = hi.div(hi_divisor, prec, Round::TowardPositive)?;
    Ball::from_bounds(lo, hi, prec)
}

impl MpComplex {
    /// Constructs an [`MpComplex`] from real and imaginary [`Ball`]s.
    #[must_use]
    pub const fn new(
        re: Ball,
        im: Ball,
    ) -> Self {
        Self { re, im }
    }

    /// Constructs an exact point [`MpComplex`] from two [`Real`] values.
    #[must_use]
    pub fn point(
        re: Real,
        im: Real,
    ) -> Self {
        Self {
            re: Ball::point(re),
            im: Ball::point(im),
        }
    }

    /// Lifts a real [`Ball`] to `re + 0i`.
    #[must_use]
    pub fn from_real(re: Ball) -> Self {
        let prec = re.mid.precision();
        Self {
            re,
            im: Ball::point(Real::zero(prec)),
        }
    }

    /// `0 + 0i` at precision `prec`.
    #[must_use]
    pub fn zero(prec: u32) -> Self {
        Self {
            re: Ball::point(Real::zero(prec)),
            im: Ball::point(Real::zero(prec)),
        }
    }

    /// `1 + 0i` at precision `prec`.
    #[must_use]
    pub fn one(prec: u32) -> Self {
        Self {
            re: Ball::point(Real::one(prec)),
            im: Ball::point(Real::zero(prec)),
        }
    }

    /// `0 + 1i` at precision `prec`.
    #[must_use]
    pub fn i(prec: u32) -> Self {
        Self {
            re: Ball::point(Real::zero(prec)),
            im: Ball::point(Real::one(prec)),
        }
    }

    /// Negation `-self`.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }

    /// Complex conjugate `re - i * im`.
    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: self.im.neg(),
        }
    }

    /// Addition `self + other`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the supported range.
    pub fn add(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        Ok(Self {
            re: self.re.add(&other.re, prec)?,
            im: self.im.add(&other.im, prec)?,
        })
    }

    /// Subtraction `self - other`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the supported range.
    pub fn sub(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        Ok(Self {
            re: self.re.sub(&other.re, prec)?,
            im: self.im.sub(&other.im, prec)?,
        })
    }

    /// Multiplication `self * other`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the supported range.
    pub fn mul(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let ac = self.re.mul(&other.re, prec)?;
        let bd = self.im.mul(&other.im, prec)?;
        let ad = self.re.mul(&other.im, prec)?;
        let bc = self.im.mul(&other.re, prec)?;
        Ok(Self {
            re: ac.sub(&bd, prec)?,
            im: ad.add(&bc, prec)?,
        })
    }

    /// Division `self / other`.
    ///
    /// # Errors
    /// Fails if the divisor contains the origin or an exponent leaves the range.
    pub fn div(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        if other.re.contains_zero() && other.im.contains_zero() {
            return Err("complex division by an enclosure containing zero".to_string());
        }
        let inner = prec.saturating_add(16);
        let num = self.mul(&other.conj(), prec)?;
        if other.re.rad.is_zero() && other.im.rad.is_zero() {
            let denom = other.abs_sq(inner)?;
            return Ok(Self {
                re: num.re.div(&denom, prec)?,
                im: num.im.div(&denom, prec)?,
            });
        }
        // A midpoint and radius can lose a small positive lower bound when the
        // upper bound is large. Keep the directed endpoints through division.
        let (lower, upper) = other.abs_squared_bounds(inner)?;
        Ok(Self {
            re: divide_by_positive_bounds(&num.re, &lower, &upper, prec)?,
            im: divide_by_positive_bounds(&num.im, &lower, &upper, prec)?,
        })
    }

    /// Reciprocal `1 / self`.
    ///
    /// # Errors
    /// Fails if the rectangle contains the origin or an exponent leaves the range.
    pub fn recip(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        Self::one(prec).div(self, prec)
    }

    /// Squared modulus `re^2 + im^2`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn abs_sq(
        &self,
        prec: u32,
    ) -> Result<Ball, String> {
        if self.re.rad.is_zero() && self.im.rad.is_zero() {
            let real = self.re.mul(&self.re, prec)?;
            let imag = self.im.mul(&self.im, prec)?;
            return real.add(&imag, prec);
        }
        let (lower, upper) = self.abs_squared_bounds(prec)?;
        Ball::from_bounds(lower, upper, prec)
    }

    fn abs_squared_bounds(
        &self,
        prec: u32,
    ) -> Result<(Real, Real), String> {
        let (real_low, real_high) = squared_coordinate_bounds(&self.re, prec)?;
        let (imag_low, imag_high) = squared_coordinate_bounds(&self.im, prec)?;
        Ok((
            real_low.add(&imag_low, prec, Round::TowardNegative)?,
            real_high.add(&imag_high, prec, Round::TowardPositive)?,
        ))
    }

    /// Modulus `|self| = sqrt(re^2 + im^2)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn abs(
        &self,
        prec: u32,
    ) -> Result<Ball, String> {
        if self.re.rad.is_zero() && self.im.rad.is_zero() {
            return funcs::hypot(&self.re.mid, &self.im.mid, prec);
        }
        let inner = prec.saturating_add(16);
        let (real_near, real_far) = self.re.absolute_endpoints(inner)?;
        let (imag_near, imag_far) = self.im.absolute_endpoints(inner)?;
        let (lower, _) = real_near.hypot_with_inexact(&imag_near, inner, Round::TowardNegative)?;
        let (upper, _) = real_far.hypot_with_inexact(&imag_far, inner, Round::TowardPositive)?;
        Ball::from_bounds(lower, upper, prec)
    }

    /// Principal argument `arg(self) = atan2(im, re)` in `(-pi, pi]`.
    ///
    /// # Errors
    /// Fails if the enclosure meets the non-positive real axis cut.
    pub fn arg(
        &self,
        prec: u32,
    ) -> Result<Ball, String> {
        if self.im.contains_zero() && self.re.mid.sign() <= 0 {
            return Err("argument of an enclosure meeting the non-positive real axis".into());
        }
        funcs::atan2_ball(&self.im, &self.re, prec)
    }

    /// Complex exponential `exp(x + i*y) = exp(x) * (cos(y) + i*sin(y))`.
    ///
    /// # Errors
    /// Fails when an exponent leaves the range or the imaginary midpoint
    /// lies outside the certified trigonometric range.
    pub fn exp(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(16);
        let ex = funcs::exp_ball(&self.re, inner)?;
        let trig = funcs::sin_cos_ball(&self.im, inner)?;
        Ok(Self {
            re: ex.mul(&trig.cos, prec)?,
            im: ex.mul(&trig.sin, prec)?,
        })
    }

    /// Principal complex logarithm `ln(z) = ln|z| + i*arg(z)`.
    ///
    /// # Errors
    /// Fails if the enclosure meets the branch cut `(-inf, 0]`.
    pub fn ln(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if self.re.contains_zero() && self.im.contains_zero() {
            return Err("logarithm of an enclosure containing zero".into());
        }
        let im = self.arg(prec)?;
        let inner = prec.saturating_add(16);
        let (real_near, real_far) = self.re.absolute_endpoints(inner)?;
        let (imag_near, imag_far) = self.im.absolute_endpoints(inner)?;
        let low = log_hypot(&real_near, &imag_near, inner)?;
        let lower = low.mid.sub(&low.rad, inner, Round::TowardNegative)?;
        let high = if real_near.cmp_value(&real_far) == 0 && imag_near.cmp_value(&imag_far) == 0 {
            low
        } else {
            log_hypot(&real_far, &imag_far, inner)?
        };
        let upper = high.mid.add(&high.rad, inner, Round::TowardPositive)?;
        let re = Ball::from_bounds(lower, upper, prec)?;
        Ok(Self { re, im })
    }

    /// Principal complex square root `sqrt(z)`.
    ///
    /// The negative real axis has the positive imaginary root. A rectangle
    /// straddling the branch cut can return both signs of the imaginary part.
    ///
    /// # Errors
    /// Fails when an intermediate or result exponent is out of range.
    pub fn sqrt(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if self.re.rad.is_zero() && is_zero_ball(&self.im) {
            let root = Ball::point(self.re.mid.abs()).sqrt(prec)?;
            return if self.re.mid.sign() < 0 {
                Ok(Self::new(Ball::point(Real::zero(prec)), root))
            } else {
                Ok(Self::from_real(root))
            };
        }
        let inner = prec.saturating_add(16).min(MAX_BITS);
        let (x_lo, x_hi) = self.re.outward_endpoints(inner)?;
        let (y_lo, y_hi) = self.im.outward_endpoints(inner)?;
        let (near, far) = Ball::absolute_bounds(&y_lo, &y_hi, inner);
        // hypot(x, y) + x increases in x and in |y|. The imaginary
        // magnitude uses the same formula at -x.
        let u_lo = sqrt_component(&x_lo, &near, inner, Round::TowardNegative)?;
        let u_hi = sqrt_component(&x_hi, &far, inner, Round::TowardPositive)?;
        let v_hi = sqrt_component(&x_lo.neg(), &far, inner, Round::TowardPositive)?;
        let (v_lo, v_hi) = if y_lo.sign() >= 0 {
            (
                sqrt_component(&x_hi.neg(), &near, inner, Round::TowardNegative)?,
                v_hi,
            )
        } else if y_hi.sign() < 0 {
            (
                v_hi.neg(),
                sqrt_component(&x_hi.neg(), &near, inner, Round::TowardNegative)?.neg(),
            )
        } else {
            (v_hi.neg(), v_hi)
        };
        Ok(Self::new(
            Ball::from_bounds(u_lo, u_hi, prec)?,
            Ball::from_bounds(v_lo, v_hi, prec)?,
        ))
    }

    /// Integer power `self^n` via binary exponentiation (exact at `n = 0` and for non-positive real bases).
    ///
    /// # Errors
    /// Fails if `n < 0` and `self` cannot be separated from zero, or on exponent overflow.
    pub fn powi(
        &self,
        n: i64,
        prec: u32,
    ) -> Result<Self, String> {
        if n == 0 {
            return Ok(Self::one(prec));
        }
        let inner = prec.saturating_add(16);
        let (mut base, mut exp) = if n < 0 {
            (self.recip(inner)?, n.unsigned_abs())
        } else {
            (self.clone(), n.unsigned_abs())
        };
        let mut acc = Self::one(inner);
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc.mul(&base, inner)?;
            }
            exp >>= 1;
            if exp > 0 {
                base = base.mul(&base, inner)?;
            }
        }
        let re_mid = acc.re.mid.with_precision(prec, Round::NearestEven)?;
        let re_rad = if acc.re.rad.is_zero() && re_mid.cmp_value(&acc.re.mid) == 0 {
            Real::zero(prec)
        } else {
            acc.re
                .rad
                .add(&ulp_bound(&re_mid, prec)?, prec, Round::AwayFromZero)?
        };
        let im_mid = acc.im.mid.with_precision(prec, Round::NearestEven)?;
        let im_rad = if acc.im.rad.is_zero() && im_mid.cmp_value(&acc.im.mid) == 0 {
            Real::zero(prec)
        } else {
            acc.im
                .rad
                .add(&ulp_bound(&im_mid, prec)?, prec, Round::AwayFromZero)?
        };
        Ok(Self {
            re: Ball::new(re_mid, re_rad)?,
            im: Ball::new(im_mid, im_rad)?,
        })
    }

    /// Complex power `self^w`. Uses [`Self::powi`] when `w` is an exact integer point,
    /// otherwise computes `exp(w * ln(self))`.
    ///
    /// # Errors
    /// Fails when an exact integer exponent does not fit a signed 64-bit
    /// integer, a non-integer power requires a logarithm on its branch cut,
    /// or a required exponent leaves the representable range.
    pub fn pow(
        &self,
        w: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&w.im) && w.re.rad.is_zero() && w.re.mid.is_integer() {
            return self.powi(w.re.mid.integer_exponent()?, prec);
        }
        let inner = prec.saturating_add(16);
        let l = self.ln(inner)?;
        let prod = w.mul(&l, inner)?;
        prod.exp(prec)
    }

    /// Complex sine `sin(x + i*y) = sin(x)cosh(y) + i*cos(x)sinh(y)`.
    ///
    /// # Errors
    /// Fails when the real midpoint exceeds the certified trigonometric
    /// range in magnitude or a required exponent is out of range.
    pub fn sin(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        RectangleTrig::new(&self.re, &self.im, prec)?.sin(prec)
    }

    /// Complex cosine `cos(x + i*y) = cos(x)cosh(y) - i*sin(x)sinh(y)`.
    ///
    /// # Errors
    /// Fails when the real midpoint exceeds the certified trigonometric
    /// range in magnitude or a required exponent is out of range.
    pub fn cos(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        RectangleTrig::new(&self.re, &self.im, prec)?.cos(prec)
    }

    /// Complex tangent `tan(z) = sin(z) / cos(z)`.
    ///
    /// # Errors
    /// Fails when the cosine enclosure includes zero, the real midpoint
    /// exceeds the certified trigonometric range in magnitude, or an
    /// intermediate or result exponent is out of range.
    pub fn tan(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.im) {
            return Ok(Self::from_real(funcs::tan_ball(&self.re, prec)?));
        }
        if is_zero_ball(&self.re) {
            return Ok(Self::new(
                Ball::point(Real::zero(prec)),
                funcs::tanh_ball(&self.im, prec)?,
            ));
        }
        let inner = prec.saturating_add(16).min(MAX_BITS);
        let factors = RectangleTrig::new(&self.re, &self.im, inner)?;
        factors.sin(inner)?.div(&factors.cos(inner)?, prec)
    }

    /// Complex hyperbolic sine `sinh(x + i*y) = sinh(x)cos(y) + i*cosh(x)sin(y)`.
    ///
    /// # Errors
    /// Fails when the imaginary midpoint exceeds the certified trigonometric
    /// range in magnitude or a required exponent is out of range.
    pub fn sinh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        RectangleTrig::new(&self.im, &self.re, prec)?.sinh(prec)
    }

    /// Complex hyperbolic cosine `cosh(x + i*y) = cosh(x)cos(y) + i*sinh(x)sin(y)`.
    ///
    /// # Errors
    /// Fails when the imaginary midpoint exceeds the certified trigonometric
    /// range in magnitude or a required exponent is out of range.
    pub fn cosh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        RectangleTrig::new(&self.im, &self.re, prec)?.cosh(prec)
    }

    /// Complex hyperbolic tangent `tanh(z) = sinh(z) / cosh(z)`.
    ///
    /// # Errors
    /// Fails when the hyperbolic cosine enclosure includes zero, the
    /// imaginary midpoint exceeds the certified trigonometric range in
    /// magnitude, or an intermediate or result exponent is out of range.
    pub fn tanh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.im) {
            return Ok(Self::from_real(funcs::tanh_ball(&self.re, prec)?));
        }
        if is_zero_ball(&self.re) {
            return Ok(Self::new(
                Ball::point(Real::zero(prec)),
                funcs::tan_ball(&self.im, prec)?,
            ));
        }
        let inner = prec.saturating_add(16).min(MAX_BITS);
        let factors = RectangleTrig::new(&self.im, &self.re, inner)?;
        factors.sinh(inner)?.div(&factors.cosh(inner)?, prec)
    }

    /// Multiplies a rectangle by i without rounded arithmetic.
    fn times_i(&self) -> Self {
        Self::new(self.im.neg(), self.re.clone())
    }

    fn into_times_neg_i(self) -> Self {
        Self::new(self.im, self.re.neg())
    }

    fn touches_outer_real_cuts(
        &self,
        prec: u32,
    ) -> Result<bool, String> {
        if !self.im.contains_zero() {
            return Ok(false);
        }
        let (lower, upper) = self.re.outward_endpoints(prec)?;
        let one = Real::one(1);
        Ok(lower.cmp_value(&one.neg()) <= 0 || upper.cmp_value(&one) >= 0)
    }

    /// Principal inverse sine of a closed rectangle.
    ///
    /// # Errors
    /// Fails when the rectangle touches the real cuts to the left of -1
    /// or right of 1, including their endpoints for non-real rectangles,
    /// or a required exponent is out of range. Real intervals in [-1, 1]
    /// include both endpoints.
    pub fn asin(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.im) {
            return Ok(Self::from_real(funcs::asin_ball(&self.re, prec)?));
        }
        let inner = prec.saturating_add(24).min(MAX_BITS);
        if self.touches_outer_real_cuts(inner)? {
            return Err("complex inverse sine rectangle meets a real branch cut".into());
        }
        let iz = self.times_i();
        let square = self.mul(self, inner)?;
        let root = Self::one(inner).sub(&square, inner)?.sqrt(inner)?;
        // Both logarithm arguments have positive real part off the cuts
        // and are reciprocals. Choose the sign that adds the magnitude
        // of the imaginary midpoint to the root's real part.
        if self.im.mid.sign() > 0 {
            Ok(root.sub(&iz, inner)?.ln(prec)?.times_i())
        } else {
            Ok(root.add(&iz, inner)?.ln(prec)?.into_times_neg_i())
        }
    }

    /// Principal inverse cosine of a closed rectangle.
    ///
    /// # Errors
    /// Fails on the same cuts and exponent limits as [Self::asin].
    pub fn acos(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.im) {
            return Ok(Self::from_real(funcs::acos_ball(&self.re, prec)?));
        }
        let inner = prec.saturating_add(16).min(MAX_BITS);
        let angle = self.asin(inner)?;
        Self::from_real(half_pi_ball(inner)?).sub(&angle, prec)
    }

    /// Principal inverse tangent of a closed rectangle.
    ///
    /// # Errors
    /// Fails if the rectangle touches the imaginary axis at or beyond
    /// either of ±i, or a required exponent is out of range.
    pub fn atan(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.im) {
            return Ok(Self::from_real(funcs::atan_ball(&self.re, prec)?));
        }
        Ok(self.times_i().atanh(prec)?.into_times_neg_i())
    }

    /// Principal inverse hyperbolic sine of a closed rectangle.
    ///
    /// # Errors
    /// Fails if the rectangle meets an imaginary branch cut or a required
    /// exponent is out of range. Imaginary intervals from -i to i include
    /// the endpoints.
    pub fn asinh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        Ok(self.times_i().asin(prec)?.into_times_neg_i())
    }

    /// Principal inverse hyperbolic cosine of a closed rectangle.
    ///
    /// # Errors
    /// Fails if the rectangle meets the real axis below one or a required
    /// exponent is out of range. One is included.
    pub fn acosh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24).min(MAX_BITS);
        if self.im.contains_zero() {
            let (lower, _) = self.re.outward_endpoints(inner)?;
            if lower.cmp_value(&Real::one(1)) < 0 {
                return Err(
                    "complex inverse hyperbolic cosine rectangle meets its real branch cut".into(),
                );
            }
        }
        let one = Self::one(inner);
        let plus = self.add(&one, inner)?.sqrt(inner)?;
        let minus = self.sub(&one, inner)?.sqrt(inner)?;
        self.add(&plus.mul(&minus, inner)?, inner)?.ln(prec)
    }

    /// Principal inverse hyperbolic tangent of a closed rectangle.
    ///
    /// # Errors
    /// Fails if the rectangle touches either real cut at or beyond ±1,
    /// or a required exponent is out of range.
    pub fn atanh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if is_zero_ball(&self.re) {
            let angle = funcs::atan_ball(&self.im, prec)?;
            return Ok(Self::new(Ball::point(Real::zero(prec)), angle));
        }
        let inner = prec.saturating_add(24).min(MAX_BITS);
        if self.touches_outer_real_cuts(inner)? {
            return Err(
                "complex inverse hyperbolic tangent rectangle meets a real branch cut".into(),
            );
        }
        let one = Self::one(inner);
        // Off the cuts, (1+z)/(1-z) has the sign of Im(z) in its
        // imaginary part. The difference of principal arguments is in
        // (-pi, pi), so this form avoids dividing interval rectangles.
        let plus = one.add(self, inner)?.ln(inner)?;
        let minus = one.sub(self, inner)?.ln(inner)?;
        let difference = plus.sub(&minus, inner)?;
        Ok(Self::new(
            half_ball(&difference.re, prec)?,
            half_ball(&difference.im, prec)?,
        ))
    }
}

/// Lexical token for the complex arbitrary-precision expression parser.
#[derive(Clone, Debug, PartialEq)]
enum Token {
    Number(String),
    Name(String),
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Open,
    Close,
    Comma,
}

fn tokenize(text: &str) -> Result<Vec<Token>, String> {
    if text.len() > MAX_TEXT {
        return Err(format!("expression exceeds {MAX_TEXT} bytes"));
    }
    let mut tokens = Vec::new();
    let chars: Vec<char> = text.chars().collect();
    let mut at = 0;
    while let Some(&ch) = chars.get(at) {
        match ch {
            | ' ' | '\t' | '\n' | '\r' => at += 1,
            | '+' => {
                tokens.push(Token::Plus);
                at += 1;
            },
            | '-' => {
                tokens.push(Token::Minus);
                at += 1;
            },
            | '*' => {
                tokens.push(Token::Star);
                at += 1;
            },
            | '/' => {
                tokens.push(Token::Slash);
                at += 1;
            },
            | '^' => {
                tokens.push(Token::Caret);
                at += 1;
            },
            | '(' => {
                tokens.push(Token::Open);
                at += 1;
            },
            | ')' => {
                tokens.push(Token::Close);
                at += 1;
            },
            | ',' => {
                tokens.push(Token::Comma);
                at += 1;
            },
            | '0'..='9' | '.' => {
                let start = at;
                while chars
                    .get(at)
                    .is_some_and(|c| c.is_ascii_digit() || *c == '.')
                {
                    at += 1;
                }
                if chars.get(at).is_some_and(|c| *c == 'e' || *c == 'E') {
                    let save = at;
                    at += 1;
                    if chars.get(at).is_some_and(|c| *c == '+' || *c == '-') {
                        at += 1;
                    }
                    let dstart = at;
                    while chars.get(at).is_some_and(char::is_ascii_digit) {
                        at += 1;
                    }
                    if at == dstart {
                        at = save;
                    }
                }
                let slice: String = chars.get(start..at).unwrap_or(&[]).iter().collect();
                tokens.push(Token::Number(slice));
            },
            | _ if ch.is_alphabetic() || ch == '_' => {
                let start = at;
                while chars
                    .get(at)
                    .is_some_and(|c| c.is_alphanumeric() || *c == '_')
                {
                    at += 1;
                }
                let slice: String = chars.get(start..at).unwrap_or(&[]).iter().collect();
                tokens.push(Token::Name(slice));
            },
            | _ => return Err(format!("unexpected character `{ch}`")),
        }
        if tokens.len() > MAX_TOKENS {
            return Err(format!("expression exceeds {MAX_TOKENS} tokens"));
        }
    }
    Ok(tokens)
}

struct Parser<'a> {
    tokens: &'a [Token],
    pos: usize,
    prec: u32,
    depth: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a Token> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<&'a Token> {
        let tok = self.tokens.get(self.pos)?;
        self.pos += 1;
        Some(tok)
    }

    fn parse_expr(&mut self) -> Result<MpComplex, String> {
        let mut lhs = self.parse_term()?;
        while let Some(tok) = self.peek() {
            match tok {
                | Token::Plus => {
                    self.pos += 1;
                    let rhs = self.parse_term()?;
                    lhs = lhs.add(&rhs, self.prec)?;
                },
                | Token::Minus => {
                    self.pos += 1;
                    let rhs = self.parse_term()?;
                    lhs = lhs.sub(&rhs, self.prec)?;
                },
                | _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_term(&mut self) -> Result<MpComplex, String> {
        let mut lhs = self.parse_power()?;
        while let Some(tok) = self.peek() {
            match tok {
                | Token::Star => {
                    self.pos += 1;
                    let rhs = self.parse_power()?;
                    lhs = lhs.mul(&rhs, self.prec)?;
                },
                | Token::Slash => {
                    self.pos += 1;
                    let rhs = self.parse_power()?;
                    lhs = lhs.div(&rhs, self.prec)?;
                },
                | _ => break,
            }
        }
        Ok(lhs)
    }

    fn parse_power(&mut self) -> Result<MpComplex, String> {
        let base = self.parse_unary()?;
        if matches!(self.peek(), Some(Token::Caret)) {
            self.pos += 1;
            let exp = self.parse_power()?;
            return base.pow(&exp, self.prec);
        }
        Ok(base)
    }

    fn parse_unary(&mut self) -> Result<MpComplex, String> {
        match self.peek() {
            | Some(Token::Plus) => {
                self.pos += 1;
                self.parse_unary()
            },
            | Some(Token::Minus) => {
                self.pos += 1;
                Ok(self.parse_unary()?.neg())
            },
            | _ => self.parse_primary(),
        }
    }

    fn parse_primary(&mut self) -> Result<MpComplex, String> {
        if self.depth >= MAX_DEPTH {
            return Err("expression nesting limit exceeded".to_string());
        }
        self.depth += 1;
        let result = match self.next() {
            | Some(Token::Number(lit)) => {
                let real_val = from_decimal(lit, self.prec)?;
                let ball = if lit.contains('.') || lit.contains('e') || lit.contains('E') {
                    Ball::new(real_val.clone(), ulp_bound(&real_val, self.prec)?)?
                } else {
                    Ball::point(real_val)
                };
                Ok(MpComplex::from_real(ball))
            },
            | Some(Token::Name(name)) => {
                if matches!(self.peek(), Some(Token::Open)) {
                    self.pos += 1;
                    let mut args = Vec::new();
                    if !matches!(self.peek(), Some(Token::Close)) {
                        loop {
                            args.push(self.parse_expr()?);
                            if matches!(self.peek(), Some(Token::Comma)) {
                                self.pos += 1;
                            } else {
                                break;
                            }
                        }
                    }
                    if !matches!(self.next(), Some(Token::Close)) {
                        return Err(format!("expected `)` after arguments to `{name}`"));
                    }
                    eval_complex_func(name, &args, self.prec)
                } else {
                    eval_complex_const(name, self.prec)
                }
            },
            | Some(Token::Open) => {
                let val = self.parse_expr()?;
                if !matches!(self.next(), Some(Token::Close)) {
                    return Err("expected `)`".to_string());
                }
                Ok(val)
            },
            | other => Err(format!("unexpected token {other:?}")),
        };
        self.depth -= 1;
        result
    }
}

fn eval_complex_const(
    name: &str,
    prec: u32,
) -> Result<MpComplex, String> {
    if name == "I" || name == "i" {
        return Ok(MpComplex::i(prec));
    }
    let val = consts::named(name, prec.saturating_add(16))?;
    let mid = val.with_precision(prec, Round::NearestEven)?;
    let rad = ulp_bound(&mid, prec)?;
    Ok(MpComplex::from_real(Ball::new(mid, rad)?))
}

fn eval_complex_func(
    name: &str,
    args: &[MpComplex],
    prec: u32,
) -> Result<MpComplex, String> {
    let one = |op_name: &str| -> Result<&MpComplex, String> {
        if args.len() != 1 {
            return Err(format!("`{op_name}` expects 1 argument"));
        }
        args.first()
            .ok_or_else(|| format!("missing argument to `{op_name}`"))
    };
    match name {
        | "re" => Ok(MpComplex::from_real(one("re")?.re.clone())),
        | "im" => Ok(MpComplex::from_real(one("im")?.im.clone())),
        | "conj" => Ok(one("conj")?.conj()),
        | "abs" => Ok(MpComplex::from_real(one("abs")?.abs(prec)?)),
        | "arg" => Ok(MpComplex::from_real(one("arg")?.arg(prec)?)),
        | "exp" => one("exp")?.exp(prec),
        | "ln" | "log" => one("ln")?.ln(prec),
        | "sqrt" => one("sqrt")?.sqrt(prec),
        | "sin" => one("sin")?.sin(prec),
        | "cos" => one("cos")?.cos(prec),
        | "tan" => one("tan")?.tan(prec),
        | "sinh" => one("sinh")?.sinh(prec),
        | "cosh" => one("cosh")?.cosh(prec),
        | "tanh" => one("tanh")?.tanh(prec),
        | "asin" => one("asin")?.asin(prec),
        | "acos" => one("acos")?.acos(prec),
        | "atan" => one("atan")?.atan(prec),
        | "asinh" => one("asinh")?.asinh(prec),
        | "acosh" => one("acosh")?.acosh(prec),
        | "atanh" => one("atanh")?.atanh(prec),
        | "erf" => one("erf")?.erf(prec),
        | "erfc" => one("erfc")?.erfc(prec),
        | "lgamma" | "ln_gamma" => one("lgamma")?.ln_gamma(prec),
        | "gamma" => one("gamma")?.gamma(prec),
        | "digamma" => one("digamma")?.digamma(prec),
        | "zeta" => one("zeta")?.riemann_zeta(prec),
        | "besselj" => {
            if args.len() != 2 {
                return Err("`besselj` expects 2 arguments".to_string());
            }
            let nu = args.first().ok_or("missing order")?;
            let z = args.get(1).ok_or("missing argument")?;
            z.bessel_j(nu, prec)
        },
        | "polylog" => {
            if args.len() != 2 {
                return Err("`polylog` expects 2 arguments".to_string());
            }
            let s = args.first().ok_or("missing order")?;
            let z = args.get(1).ok_or("missing argument")?;
            z.polylog(s, prec)
        },
        | "pow" => {
            if args.len() != 2 {
                return Err("`pow` expects 2 arguments".to_string());
            }
            let a = args.first().ok_or("missing base")?;
            let b = args.get(1).ok_or("missing exponent")?;
            a.pow(b, prec)
        },
        | _ => Err(format!("unknown complex function `{name}`")),
    }
}

fn parse_and_eval(
    text: &str,
    prec: u32,
) -> Result<MpComplex, String> {
    let tokens = tokenize(text)?;
    let mut parser = Parser {
        tokens: &tokens,
        pos: 0,
        prec,
        depth: 0,
    };
    let val = parser.parse_expr()?;
    if parser.pos != tokens.len() {
        return Err("unexpected trailing tokens".to_string());
    }
    Ok(val)
}

/// Evaluates a complex expression `text` to `digits` significant decimal digits
/// for both real and imaginary parts, using Ziv's strategy.
///
/// # Errors
/// Fails on syntax errors, branch cuts/poles, or if digits cannot be certified
/// within [`MAX_BITS`].
pub fn eval_complex_text(
    text: &str,
    digits: usize,
) -> Result<(String, String), String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let mut prec = bits_for_digits(digits).saturating_add(96);
    let mut last = String::new();
    for _ in 0..5 {
        let z = parse_and_eval(text, prec)?;
        match (
            ball_to_decimal(&z.re, digits),
            ball_to_decimal(&z.im, digits),
        ) {
            | (Ok(re), Ok(im)) => return Ok((re, im)),
            | (Err(e), _) | (_, Err(e)) => last = e,
        }
        prec = prec.saturating_mul(2).min(MAX_BITS);
    }
    Err(last)
}

/// Evaluates a complex expression `text` to a rigorous enclosure at `digits`
/// significant decimal digits, returning `((re_lo, re_hi), (im_lo, im_hi))`.
///
/// # Errors
/// Fails on syntax errors, branch cuts/poles, or exponent overflow.
pub fn enclosure_complex_text(
    text: &str,
    digits: usize,
) -> Result<((String, String), (String, String)), String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let prec = bits_for_digits(digits).saturating_add(96).min(MAX_BITS);
    let z = parse_and_eval(text, prec)?;
    Ok((
        super::enclosure_bounds_text(&z.re, prec, digits)?,
        super::enclosure_bounds_text(&z.im, prec, digits)?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_identity_to_thirty_digits() {
        let ((re_lo, re_hi), (im_lo, im_hi)) =
            enclosure_complex_text("exp(I * pi)", 30).unwrap_or_else(|e| panic!("{e}"));
        // Directed decimal rounding must enclose the exact value even when
        // neither endpoint can be rendered as that value.
        let prec = 256;
        let decimal = |text: &str| from_decimal(text, prec).unwrap_or_else(|e| panic!("{e}"));
        let width_limit = decimal("2e-29");
        let minus_one = Real::from_i64(-1, prec).unwrap_or_else(|e| panic!("{e}"));
        for (lo, hi, expected) in [(re_lo, re_hi, minus_one), (im_lo, im_hi, Real::zero(prec))] {
            let (lower, upper) = (decimal(&lo), decimal(&hi));
            assert!(
                lower.cmp_value(&expected) <= 0 && upper.cmp_value(&expected) >= 0,
                "{lo} .. {hi}"
            );
            let width = upper
                .sub(&lower, prec, Round::TowardPositive)
                .unwrap_or_else(|e| panic!("{e}"));
            assert!(width.cmp_value(&width_limit) < 0, "{lo} .. {hi}");
        }

        let (re, im) = eval_complex_text("exp(I * pi / 3)", 30).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(re, "0.5");
        assert!(im.starts_with("0.86602540378443864676372317075"));
    }

    #[test]
    fn pythagorean_and_inverse_identities() {
        let (re, im) = eval_complex_text("sin(1 + 2*I)^2 + cos(1 + 2*I)^2 + I", 35)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(re, "1");
        assert_eq!(im, "1");

        let (re2, im2) =
            eval_complex_text("ln(exp(0.5 + 0.7*I))", 30).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(re2, "0.5");
        assert_eq!(im2, "0.7");
    }
}
