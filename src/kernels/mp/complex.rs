//! Arbitrary-precision complex ball arithmetic (`MpComplex`) and
//! Ziv-certified complex expression evaluator.

use super::MAX_BITS;
use super::MAX_DIGITS;
use super::ball::Ball;
use super::consts;
use super::decimal::ball_to_decimal;
use super::decimal::bits_for_digits;
use super::decimal::from_decimal;
use super::decimal::to_decimal_in;
use super::funcs;
use super::real::Real;
use super::real::Round;
use num_bigint::BigInt;
use num_traits::ToPrimitive;

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

/// Real `Ball` square root for a non-negative enclosure.
fn ball_sqrt(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let low = x.mid.sub(&x.rad, inner, Round::TowardNegative)?;
    if low.sign() < 0 {
        return Err("square root of a real ball that may be negative".to_string());
    }
    x.map_monotone(&|v, p, r| v.sqrt(p, r), true, prec)
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

/// Half-pi as a rigorous `Ball` enclosure.
fn half_pi_ball(prec: u32) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let pi = consts::pi(inner)?;
    let half_pi = pi.scale_pow2(-1)?;
    let mid = half_pi.with_precision(prec, Round::NearestEven)?;
    let rad = ulp_bound(&mid, prec)?;
    Ball::new(mid, rad)
}

/// Sound `Ball` enclosure for `sin` or `cos` using `|f'(t)| <= 1`.
fn ball_unit_lipschitz(
    x: &Ball,
    f: &dyn Fn(&Real, u32) -> Result<Ball, String>,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let v = f(&x.mid, inner)?;
    let mid = v.mid.with_precision(prec, Round::NearestEven)?;
    let rad = v.rad.add(&x.rad, prec, Round::AwayFromZero)?.add(
        &ulp_bound(&mid, prec)?,
        prec,
        Round::AwayFromZero,
    )?;
    Ball::new(mid, rad)
}

/// Sound `Ball` enclosure for `cosh` using Mean Value Theorem `|cosh'(t)| <= sinh(|mid| + rad)`.
fn ball_cosh(
    x: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(24);
    let c = funcs::cosh(&x.mid, inner)?;
    let mid = c.mid.with_precision(prec, Round::NearestEven)?;
    let rad = if x.rad.is_zero() {
        c.rad.with_precision(prec, Round::AwayFromZero)?.add(
            &ulp_bound(&mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?
    } else {
        let abs_hi = x.mid.abs().add(&x.rad, inner, Round::AwayFromZero)?;
        let s = funcs::sinh(&abs_hi, inner)?;
        let lip = s.mid.add(&s.rad, inner, Round::AwayFromZero)?;
        let prop = lip.mul(&x.rad, inner, Round::AwayFromZero)?;
        c.rad.add(&prop, prec, Round::AwayFromZero)?.add(
            &ulp_bound(&mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?
    };
    Ball::new(mid, rad)
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
    /// Fails if `|other|^2` contains zero or an exponent leaves the range.
    pub fn div(
        &self,
        other: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let denom = other.abs_sq(prec)?;
        if denom.contains_zero() {
            return Err("complex division by an enclosure containing zero".to_string());
        }
        let num = self.mul(&other.conj(), prec)?;
        Ok(Self {
            re: num.re.div(&denom, prec)?,
            im: num.im.div(&denom, prec)?,
        })
    }

    /// Reciprocal `1 / self`.
    ///
    /// # Errors
    /// Fails if `|self|^2` contains zero or an exponent leaves the range.
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
        let r2 = self.re.mul(&self.re, prec)?;
        let i2 = self.im.mul(&self.im, prec)?;
        let sum = r2.add(&i2, prec)?;
        // Clamp any negative lower endpoint due to outward rounding to >= 0.
        let inner = prec.saturating_add(8);
        let lo = sum.mid.sub(&sum.rad, inner, Round::TowardNegative)?;
        if lo.sign() < 0 {
            let hi = sum.mid.add(&sum.rad, inner, Round::TowardPositive)?;
            let mid_raw = hi.scale_pow2(-1)?;
            let mid = mid_raw.with_precision(prec, Round::NearestEven)?;
            let rad = mid_raw.with_precision(prec, Round::AwayFromZero)?.add(
                &ulp_bound(&mid, prec)?,
                prec,
                Round::AwayFromZero,
            )?;
            return Ball::new(mid, rad);
        }
        Ok(sum)
    }

    /// Modulus `|self| = sqrt(re^2 + im^2)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn abs(
        &self,
        prec: u32,
    ) -> Result<Ball, String> {
        let h = funcs::hypot(&self.re.mid, &self.im.mid, prec)?;
        let extra = self.re.rad.add(&self.im.rad, prec, Round::AwayFromZero)?;
        let rad = h.rad.add(&extra, prec, Round::AwayFromZero)?;
        Ball::new(h.mid, rad)
    }

    /// Principal argument `arg(self) = atan2(im, re)` in `(-pi, pi]`.
    ///
    /// # Errors
    /// Fails if the enclosure meets the non-positive real axis cut.
    pub fn arg(
        &self,
        prec: u32,
    ) -> Result<Ball, String> {
        let inner = prec.saturating_add(16);
        let r = self.abs(inner)?;
        let d_min = r.mid.sub(&r.rad, inner, Round::TowardZero)?;
        if d_min.sign() <= 0 {
            return Err("argument of an enclosure containing zero".to_string());
        }
        let re_hi = self
            .re
            .mid
            .add(&self.re.rad, inner, Round::TowardPositive)?;
        if re_hi.sign() <= 0 && self.im.contains_zero() {
            return Err("argument of an enclosure crossing the branch cut".to_string());
        }
        let a = funcs::atan2(&self.im.mid, &self.re.mid, prec)?;
        let r_sum = self.re.rad.add(&self.im.rad, inner, Round::AwayFromZero)?;
        let lip = r_sum.div(&d_min, prec, Round::AwayFromZero)?;
        let rad = a.rad.add(&lip, prec, Round::AwayFromZero)?;
        Ball::new(a.mid, rad)
    }

    /// Complex exponential `exp(x + i*y) = exp(x) * (cos(y) + i*sin(y))`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn exp(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let ex = funcs::map(&self.re, &funcs::exp, true, prec)?;
        let cy = funcs::map(&self.im, &funcs::cos, false, prec)?;
        let sy = funcs::map(&self.im, &funcs::sin, false, prec)?;
        Ok(Self {
            re: ex.mul(&cy, prec)?,
            im: ex.mul(&sy, prec)?,
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
        let r = self.abs(prec)?;
        if r.contains_zero() {
            return Err("logarithm of an enclosure containing zero".to_string());
        }
        let re = funcs::map(&r, &funcs::ln, true, prec)?;
        let im = self.arg(prec)?;
        Ok(Self { re, im })
    }

    /// Principal complex square root `sqrt(z)`.
    ///
    /// # Errors
    /// Fails if the enclosure crosses the negative real axis `(-inf, 0)`.
    pub fn sqrt(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        if self.re.mid.is_zero()
            && self.re.rad.is_zero()
            && self.im.mid.is_zero()
            && self.im.rad.is_zero()
        {
            return Ok(Self::zero(prec));
        }
        let inner = prec.saturating_add(16);
        let r = self.abs(inner)?;
        if self.re.contains_zero() && self.im.contains_zero() {
            let r_hi = r.mid.add(&r.rad, inner, Round::TowardPositive)?;
            if r_hi.is_zero() {
                return Ok(Self::zero(prec));
            }
            let bound = r_hi.sqrt(prec, Round::AwayFromZero)?;
            let zero = Real::zero(prec);
            return Ok(Self {
                re: Ball::new(zero.clone(), bound.clone())?,
                im: Ball::new(zero, bound)?,
            });
        }
        let re_lo = self
            .re
            .mid
            .sub(&self.re.rad, inner, Round::TowardNegative)?;
        if re_lo.sign() > 0 {
            // Right half-plane: u = sqrt((r + re) / 2), v = im / (2 * u).
            let r_plus_re = r.add(&self.re, inner)?;
            let half_sum = Ball::new(r_plus_re.mid.scale_pow2(-1)?, r_plus_re.rad.scale_pow2(-1)?)?;
            let u = ball_sqrt(&half_sum, prec)?;
            let two_u = Ball::new(u.mid.scale_pow2(1)?, u.rad.scale_pow2(1)?)?;
            let v = self.im.div(&two_u, prec)?;
            return Ok(Self { re: u, im: v });
        }
        let l = self.ln(inner)?;
        let half_l = Self {
            re: Ball::new(l.re.mid.scale_pow2(-1)?, l.re.rad.scale_pow2(-1)?)?,
            im: Ball::new(l.im.mid.scale_pow2(-1)?, l.im.rad.scale_pow2(-1)?)?,
        };
        half_l.exp(prec)
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
    /// Fails if `self` meets the logarithm cut or an exponent overflows.
    pub fn pow(
        &self,
        w: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        if w.im.mid.is_zero() && w.im.rad.is_zero() && w.re.rad.is_zero() {
            if w.re.mid.is_zero() {
                return Ok(Self::one(prec));
            }
            let mant = w.re.mid.mantissa();
            let e = w.re.mid.exponent();
            let int_val = if e >= 0 {
                u32::try_from(e)
                    .ok()
                    .and_then(|shift| mant.to_i64()?.checked_shl(shift))
            } else {
                let shift = e.unsigned_abs();
                if mant.trailing_zeros().is_some_and(|tz| tz >= shift) {
                    (mant >> shift).to_i64()
                } else {
                    None
                }
            };
            if let Some(n) = int_val {
                return self.powi(n, prec);
            }
        }
        let inner = prec.saturating_add(16);
        let l = self.ln(inner)?;
        let prod = w.mul(&l, inner)?;
        prod.exp(prec)
    }

    /// Complex sine `sin(x + i*y) = sin(x)cosh(y) + i*cos(x)sinh(y)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn sin(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let sx = ball_unit_lipschitz(&self.re, &funcs::sin, prec)?;
        let cx = ball_unit_lipschitz(&self.re, &funcs::cos, prec)?;
        let sy = funcs::map(&self.im, &funcs::sinh, true, prec)?;
        let cy = ball_cosh(&self.im, prec)?;
        Ok(Self {
            re: sx.mul(&cy, prec)?,
            im: cx.mul(&sy, prec)?,
        })
    }

    /// Complex cosine `cos(x + i*y) = cos(x)cosh(y) - i*sin(x)sinh(y)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn cos(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let sx = ball_unit_lipschitz(&self.re, &funcs::sin, prec)?;
        let cx = ball_unit_lipschitz(&self.re, &funcs::cos, prec)?;
        let sy = funcs::map(&self.im, &funcs::sinh, true, prec)?;
        let cy = ball_cosh(&self.im, prec)?;
        Ok(Self {
            re: cx.mul(&cy, prec)?,
            im: sx.mul(&sy, prec)?.neg(),
        })
    }

    /// Complex tangent `tan(z) = sin(z) / cos(z)`.
    ///
    /// # Errors
    /// Fails if `cos(z)` contains zero or an exponent leaves the range.
    pub fn tan(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(16);
        self.sin(inner)?.div(&self.cos(inner)?, prec)
    }

    /// Complex hyperbolic sine `sinh(x + i*y) = sinh(x)cos(y) + i*cosh(x)sin(y)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn sinh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let sx = funcs::map(&self.re, &funcs::sinh, true, prec)?;
        let cx = ball_cosh(&self.re, prec)?;
        let sy = ball_unit_lipschitz(&self.im, &funcs::sin, prec)?;
        let cy = ball_unit_lipschitz(&self.im, &funcs::cos, prec)?;
        Ok(Self {
            re: sx.mul(&cy, prec)?,
            im: cx.mul(&sy, prec)?,
        })
    }

    /// Complex hyperbolic cosine `cosh(x + i*y) = cosh(x)cos(y) + i*sinh(x)sin(y)`.
    ///
    /// # Errors
    /// Fails if an exponent leaves the range.
    pub fn cosh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let sx = funcs::map(&self.re, &funcs::sinh, true, prec)?;
        let cx = ball_cosh(&self.re, prec)?;
        let sy = ball_unit_lipschitz(&self.im, &funcs::sin, prec)?;
        let cy = ball_unit_lipschitz(&self.im, &funcs::cos, prec)?;
        Ok(Self {
            re: cx.mul(&cy, prec)?,
            im: sx.mul(&sy, prec)?,
        })
    }

    /// Complex hyperbolic tangent `tanh(z) = sinh(z) / cosh(z)`.
    ///
    /// # Errors
    /// Fails if `cosh(z)` contains zero or an exponent leaves the range.
    pub fn tanh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(16);
        self.sinh(inner)?.div(&self.cosh(inner)?, prec)
    }

    /// Principal inverse sine `asin(z) = -i * ln(i*z + sqrt(1 - z^2))`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn asin(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let i_unit = Self::i(inner);
        let iz = i_unit.mul(self, inner)?;
        let z2 = self.mul(self, inner)?;
        let root = Self::one(inner).sub(&z2, inner)?.sqrt(inner)?;
        let l = iz.add(&root, inner)?.ln(inner)?;
        Self {
            re: l.im,
            im: l.re.neg(),
        }
        .mul(&Self::one(prec), prec)
    }

    /// Principal inverse cosine `acos(z) = pi/2 - asin(z)`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn acos(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(16);
        let a = self.asin(inner)?;
        let hp = Self::from_real(half_pi_ball(inner)?);
        hp.sub(&a, prec)
    }

    /// Principal inverse tangent `atan(z) = (i/2) * ln((1 - i*z) / (1 + i*z))`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn atan(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let i_unit = Self::i(inner);
        let iz = i_unit.mul(self, inner)?;
        let one = Self::one(inner);
        let num = one.sub(&iz, inner)?;
        let den = one.add(&iz, inner)?;
        let l = num.div(&den, inner)?.ln(inner)?;
        // (i/2) * (a + i b) = (-b/2) + i (a/2)
        let re = Ball::new(l.im.neg().mid.scale_pow2(-1)?, l.im.rad.scale_pow2(-1)?)?;
        let im = Ball::new(l.re.mid.scale_pow2(-1)?, l.re.rad.scale_pow2(-1)?)?;
        let re_mid = re.mid.with_precision(prec, Round::NearestEven)?;
        let re_rad = re.rad.with_precision(prec, Round::AwayFromZero)?.add(
            &ulp_bound(&re_mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?;
        let im_mid = im.mid.with_precision(prec, Round::NearestEven)?;
        let im_rad = im.rad.with_precision(prec, Round::AwayFromZero)?.add(
            &ulp_bound(&im_mid, prec)?,
            prec,
            Round::AwayFromZero,
        )?;
        Ok(Self {
            re: Ball::new(re_mid, re_rad)?,
            im: Ball::new(im_mid, im_rad)?,
        })
    }

    /// Principal inverse hyperbolic sine `asinh(z) = ln(z + sqrt(z^2 + 1))`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn asinh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let z2_plus_1 = self.mul(self, inner)?.add(&Self::one(inner), inner)?;
        let root = z2_plus_1.sqrt(inner)?;
        self.add(&root, inner)?.ln(prec)
    }

    /// Principal inverse hyperbolic cosine `acosh(z) = ln(z + sqrt(z - 1) * sqrt(z + 1))`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn acosh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let one = Self::one(inner);
        let zp = self.add(&one, inner)?.sqrt(inner)?;
        let zm = self.sub(&one, inner)?.sqrt(inner)?;
        let prod = zp.mul(&zm, inner)?;
        self.add(&prod, inner)?.ln(prec)
    }

    /// Principal inverse hyperbolic tangent `atanh(z) = 0.5 * ln((1 + z) / (1 - z))`.
    ///
    /// # Errors
    /// Fails if the expression crosses a branch cut.
    pub fn atanh(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let one = Self::one(inner);
        let q = one.add(self, inner)?.div(&one.sub(self, inner)?, inner)?;
        let l = q.ln(inner)?;
        let re_mid =
            l.re.mid
                .scale_pow2(-1)?
                .with_precision(prec, Round::NearestEven)?;
        let re_rad =
            l.re.rad
                .scale_pow2(-1)?
                .with_precision(prec, Round::AwayFromZero)?
                .add(&ulp_bound(&re_mid, prec)?, prec, Round::AwayFromZero)?;
        let im_mid =
            l.im.mid
                .scale_pow2(-1)?
                .with_precision(prec, Round::NearestEven)?;
        let im_rad =
            l.im.rad
                .scale_pow2(-1)?
                .with_precision(prec, Round::AwayFromZero)?
                .add(&ulp_bound(&im_mid, prec)?, prec, Round::AwayFromZero)?;
        Ok(Self {
            re: Ball::new(re_mid, re_rad)?,
            im: Ball::new(im_mid, im_rad)?,
        })
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
    let (re_lo, re_hi) = z.re.endpoints(prec)?;
    let (im_lo, im_hi) = z.im.endpoints(prec)?;
    Ok((
        (
            to_decimal_in(&re_lo, digits, Round::TowardNegative)?,
            to_decimal_in(&re_hi, digits, Round::TowardPositive)?,
        ),
        (
            to_decimal_in(&im_lo, digits, Round::TowardNegative)?,
            to_decimal_in(&im_hi, digits, Round::TowardPositive)?,
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_identity_to_thirty_digits() {
        let ((re_lo, re_hi), (im_lo, im_hi)) =
            enclosure_complex_text("exp(I * pi)", 30).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(re_lo, "-1");
        assert_eq!(re_hi, "-1");
        assert!(im_lo.starts_with('-') || im_lo == "0", "{im_lo} .. {im_hi}");

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
