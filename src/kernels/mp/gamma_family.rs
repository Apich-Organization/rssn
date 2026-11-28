//! Outward Stirling enclosures shared by the gamma family and Bessel J.

use num_bigint::BigInt;
use num_rational::BigRational;
use std::sync::Mutex;
use std::sync::OnceLock;

use super::super::ball::Ball;
use super::super::complex::MpComplex;
use super::super::consts;
use super::super::funcs;
use super::super::real::Real;
use super::super::real::Round;
use super::add_error;
use super::divide_real;
use super::is_zero;

const MAX_BITS: u32 = 1024;
const MAX_SHIFT: i64 = 2304;
const STIRLING_GUARD: u32 = 40;

fn integer(value: i64) -> Result<Ball, String> {
    Ok(Ball::point(Real::from_i64(value, 32)?))
}

pub(super) fn rational(
    value: &BigRational,
    bits: u32,
) -> Result<Ball, String> {
    let denominator = value.denom().magnitude();
    let lower = Real::rounded_ratio(value.numer(), denominator, 0, bits, Round::TowardNegative)?;
    let upper = Real::rounded_ratio(value.numer(), denominator, 0, bits, Round::TowardPositive)?;
    Ball::from_bounds(lower, upper, bits)
}

pub(super) fn bernoulli_even(count: u32) -> Result<Vec<BigRational>, String> {
    // The even Bernoulli prefix is exact, independent of the requested precision.
    static COEFFICIENTS: OnceLock<Mutex<Vec<BigRational>>> = OnceLock::new();
    let mut cache = COEFFICIENTS
        .get_or_init(|| Mutex::new(vec![BigRational::from_integer(1.into())]))
        .lock()
        .map_err(|_| "Stirling coefficient cache poisoned".to_owned())?;
    while cache.len() < count as usize {
        let m = 2 * i64::try_from(cache.len()).map_err(|error| error.to_string())?;
        // Sum binomial(m+1,j) B_j, with B_1=-1/2 and the later odd B_j zero.
        let mut sum = BigRational::new((1 - m).into(), 2.into());
        let mut binomial = BigInt::from(1);
        for (k, coefficient) in (2..m).step_by(2).zip(cache.iter().skip(1)) {
            binomial *= (m + 3 - k) * (m + 2 - k);
            binomial /= (k - 1) * k;
            sum += coefficient * &binomial;
        }
        cache.push(-sum / BigInt::from(m + 1));
    }
    Ok(cache
        .iter()
        .skip(1)
        .take(count.saturating_sub(1) as usize)
        .cloned()
        .collect())
}

struct Stirling {
    shifted: MpComplex,
    shift: i64,
    terms: u32,
    bits: u32,
}

impl Stirling {
    fn new(
        value: &MpComplex,
        bits: u32,
    ) -> Result<Self, String> {
        let bits = bits + STIRLING_GUARD;
        let terms = (bits + 8).div_ceil(8);
        let target = i64::from(8 * terms);
        let (real_low, _) = value.re.outward_endpoints(bits)?;
        let shift = if real_low.cmp_value(&integer(target)?.mid) >= 0 {
            0
        } else {
            if real_low.cmp_value(&integer(target - MAX_SHIFT + 1)?.mid) < 0 {
                return Err("complex gamma family exhausted its recurrence shift budget".into());
            }
            // Compensation for truncation toward zero keeps the whole
            // shifted rectangle in Re(z) >= 8*terms.
            let low =
                i64::try_from(real_low.trunc_to_bigint()?).map_err(|error| error.to_string())?;
            target - low + 1
        };
        let shifted = value.add(&MpComplex::from_real(integer(shift)?), bits)?;
        Ok(Self {
            shifted,
            shift,
            terms,
            bits,
        })
    }

    fn expansion(
        &self,
        derivative: bool,
    ) -> Result<MpComplex, String> {
        let bits = self.bits;
        let z = &self.shifted;
        let half = Ball::point(Real::one(1).scale_pow2(-1)?);
        let inverse = z.recip(bits)?;
        let inverse_squared = inverse.mul(&inverse, bits)?;
        let log_z = z.ln(bits)?;
        let mut result = if derivative {
            let half_inverse =
                MpComplex::new(inverse.re.mul(&half, bits)?, inverse.im.mul(&half, bits)?);
            log_z.sub(&half_inverse, bits)?
        } else {
            let log_two_pi = funcs::ln_ball(
                &consts::named_ball("pi", bits)?.mul(&integer(2)?, bits)?,
                bits,
            )?;
            z.sub(&MpComplex::from_real(half.clone()), bits)?
                .mul(&log_z, bits)?
                .sub(z, bits)?
                .add(&MpComplex::from_real(log_two_pi.mul(&half, bits)?), bits)?
        };
        let mut power = if derivative {
            inverse_squared.clone()
        } else {
            inverse
        };
        for (coefficient, degree) in bernoulli_even(self.terms)?
            .into_iter()
            .zip((2_i64..).step_by(2))
        {
            let divisor = if derivative {
                -degree
            } else {
                degree * (degree - 1)
            };
            let coefficient = rational(&(coefficient / BigInt::from(divisor)), bits)?;
            result = result.add(
                &MpComplex::new(
                    power.re.mul(&coefficient, bits)?,
                    power.im.mul(&coefficient, bits)?,
                ),
                bits,
            )?;
            power = power.mul(&inverse_squared, bits)?;
        }
        // DLMF 5.11(ii), 25.6.2: |B_(2n)| <= 4(2n)!/6^(2n), and in
        // this right half-plane sec(arg(z)/2)^(2n+1) <= 2^(n+1).
        // For Re(z) >= 8n the first omitted term times its sector
        // factor is bounded by 8/(n*288^n) in both expansions.
        let error = Real::one(1).scale_pow2(3 - i64::from(8 * self.terms))?;
        add_error(&mut result, &error, is_zero(&z.im), false, bits)?;
        Ok(result)
    }

    fn factors<'a>(
        &self,
        value: &'a MpComplex,
    ) -> impl Iterator<Item = Result<MpComplex, String>> + 'a {
        let bits = self.bits;
        (0..self.shift).map(move |k| value.add(&MpComplex::from_real(integer(k)?), bits))
    }

    fn product(
        &self,
        value: &MpComplex,
    ) -> Result<MpComplex, String> {
        self.factors(value)
            .try_fold(MpComplex::one(self.bits), |product, factor| {
                let factor = factor?;
                if factor.re.contains_zero() && factor.im.contains_zero() {
                    return Err("complex gamma enclosure meets a non-positive integer pole".into());
                }
                product.mul(&factor, self.bits)
            })
    }
}

fn requested(precision: u32) -> Result<u32, String> {
    if precision > MAX_BITS {
        Err(format!(
            "complex gamma family supports up to {MAX_BITS} requested bits"
        ))
    } else {
        Ok(precision.max(1))
    }
}

fn output(
    value: MpComplex,
    bits: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::new(
        value.re.rounded_midpoint(bits)?,
        value.im.rounded_midpoint(bits)?,
    ))
}

fn positive_integer(value: &MpComplex) -> Option<i64> {
    if value.re.rad.is_zero() && is_zero(&value.im) && value.re.mid.is_integer() {
        return value
            .re
            .mid
            .integer_exponent()
            .ok()
            .filter(|n| (1..=512).contains(n));
    }
    None
}

pub(super) fn log_gamma(
    value: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = requested(precision)?;
    let (re_low, _) = value.re.outward_endpoints(bits + STIRLING_GUARD)?;
    if re_low.sign() <= 0 && value.im.contains_zero() {
        return Err("complex log-gamma enclosure meets the non-positive real cut".into());
    }
    if matches!(positive_integer(value), Some(1 | 2)) {
        return Ok(MpComplex::zero(bits));
    }
    let plan = Stirling::new(value, bits)?;
    let log_shifted = plan.expansion(false)?;
    if re_low.sign() > 0 && is_zero(&value.im) {
        let accumulated = plan.product(value)?.ln(plan.bits)?;
        return log_shifted.sub(&accumulated, bits);
    }
    // Continue the branch real on positive arguments. Taking a single
    // principal log of the recurrence product would lose phase windings.
    let difference = plan.factors(value).try_fold(log_shifted, |total, factor| {
        total.sub(&factor?.ln(plan.bits)?, plan.bits)
    })?;
    output(difference, bits)
}

pub(super) fn gamma(
    value: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = requested(precision)?;
    if let Some(n) = positive_integer(value) {
        let exact: BigInt = (2..n).map(BigInt::from).product();
        return Ok(MpComplex::from_real(rational(
            &BigRational::from_integer(exact),
            bits,
        )?));
    }
    let plan = Stirling::new(value, bits)?;
    let product = plan.product(value)?;
    let logarithm = plan.expansion(false)?;
    if product.re.contains_zero() && product.im.contains_zero() {
        // The interval product can contain the origin through dependency
        // even when every factor excludes its pole. Use factor logarithms
        // in that case and reflect factors strictly in the left half-plane
        // so rectangles crossing the negative axis still admit gamma.
        let mut difference = logarithm;
        let mut negate = false;
        for factor in plan.factors(value) {
            let mut factor = factor?;
            if factor.re.mid.sign() < 0 && !factor.re.contains_zero() {
                factor = factor.neg();
                negate = !negate;
            }
            difference = difference.sub(&factor.ln(plan.bits)?, plan.bits)?;
        }
        let result = difference.exp(bits)?;
        return Ok(if negate {
            result.neg()
        } else {
            result
        });
    }
    let scaled = logarithm.exp(plan.bits)?;
    if is_zero(&value.im) {
        return divide_real(&scaled, &product.re, bits);
    }
    scaled.div(&product, bits)
}

pub(super) fn reciprocal_gamma(
    value: &MpComplex,
    bits: u32,
) -> Result<MpComplex, String> {
    let plan = Stirling::new(value, bits)?;
    let shifted = plan.expansion(false)?.neg().exp(plan.bits)?;
    shifted.mul(&plan.product(value)?, bits)
}

pub(super) fn digamma(
    value: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = requested(precision)?;
    let plan = Stirling::new(value, bits)?;
    let shifted = plan.expansion(true)?;
    let difference = plan.factors(value).try_fold(shifted, |total, factor| {
        total.sub(&factor?.recip(plan.bits)?, plan.bits)
    })?;
    output(difference, bits)
}
