//! Certified Euler–Maclaurin zeta and power-series polylogarithm.

use num_bigint::BigInt;
use num_rational::BigRational;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;

use super::super::ball::Ball;
use super::super::complex::MpComplex;
use super::super::consts;
use super::super::funcs;
use super::super::real::Real;
use super::super::real::Round;
use super::add_error;
use super::gamma_family::bernoulli_even;
use super::gamma_family::rational;
use super::is_zero;

const MAX_BITS: u32 = 512;
const MAX_POLYLOG_TERMS: i64 = 8192;
const LOG_CACHE_CAPACITY: usize = 2048;
const MAX_POLYLOG_NEGATIVE_ORDER: i64 = 256;
const MAX_ZETA_NEGATIVE_ORDER: i64 = 32;
const MAX_ZETA_MAGNITUDE: i64 = 256;

fn integer(n: i64) -> Result<Real, String> {
    Real::from_i64(n, 32)
}

fn scalar(n: i64) -> Result<MpComplex, String> {
    Ok(MpComplex::from_real(Ball::point(integer(n)?)))
}

fn signed_order(value: &MpComplex) -> Option<i64> {
    if value.re.rad.is_zero() && is_zero(&value.im) && value.re.mid.is_integer() {
        return value.re.mid.integer_exponent().ok();
    }
    None
}

fn requested(precision: u32) -> Result<u32, String> {
    if precision > MAX_BITS {
        Err(format!(
            "complex Dirichlet functions support up to {MAX_BITS} requested bits"
        ))
    } else {
        Ok(precision.max(1))
    }
}

fn round(
    value: MpComplex,
    bits: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::new(
        value.re.rounded_midpoint(bits)?,
        value.im.rounded_midpoint(bits)?,
    ))
}

fn magnitude_upper(
    value: &MpComplex,
    bits: u32,
) -> Result<Real, String> {
    let real = value.re.absolute_endpoints(bits)?.1;
    let imag = value.im.absolute_endpoints(bits)?.1;
    let upper = real
        .hypot_with_inexact(&imag, bits, Round::TowardPositive)?
        .0;
    Ok(upper)
}

fn multiply_real(
    value: &MpComplex,
    factor: &Ball,
    bits: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::new(
        value.re.mul(factor, bits)?,
        value.im.mul(factor, bits)?,
    ))
}

fn scale(
    value: &MpComplex,
    exponent: i64,
) -> Result<MpComplex, String> {
    let scaled = |x: &Ball| Ball::new(x.mid.scale_pow2(exponent)?, x.rad.scale_pow2(exponent)?);
    Ok(MpComplex::new(scaled(&value.re)?, scaled(&value.im)?))
}

type CachedLog = RwLock<Option<(u32, Arc<Ball>)>>;

fn integer_log(
    n: i64,
    bits: u32,
) -> Result<Arc<Ball>, String> {
    static LOGS: OnceLock<[CachedLog; LOG_CACHE_CAPACITY]> = OnceLock::new();
    let index = usize::try_from(n).map_err(|error| error.to_string())?;
    if index >= LOG_CACHE_CAPACITY {
        return Ok(Arc::new(funcs::ln(&integer(n)?, bits)?));
    }
    let entry = &LOGS.get_or_init(|| std::array::from_fn(|_| RwLock::new(None)))[index];
    {
        let read = entry
            .read()
            .map_err(|_| "integer logarithm cache poisoned".to_owned())?;
        if let Some((cached_bits, cached)) = read.as_ref() {
            if *cached_bits >= bits {
                return Ok(Arc::clone(cached));
            }
        }
    }
    let computed = Arc::new(funcs::ln(&integer(n)?, bits)?);
    let mut write = entry
        .write()
        .map_err(|_| "integer logarithm cache poisoned".to_owned())?;
    if let Some((cached_bits, cached)) = write.as_ref() {
        if *cached_bits >= bits {
            return Ok(Arc::clone(cached));
        }
    }
    *write = Some((bits, Arc::clone(&computed)));
    drop(write);
    Ok(computed)
}

struct Weights {
    negated_order: MpComplex,
    integer_order: Option<i64>,
    bits: u32,
}

impl Weights {
    fn new(
        order: &MpComplex,
        bits: u32,
    ) -> Self {
        Self {
            negated_order: order.neg(),
            integer_order: signed_order(order).filter(|n| n.unsigned_abs() <= 256),
            bits,
        }
    }

    fn at(
        &self,
        n: i64,
    ) -> Result<MpComplex, String> {
        if n == 1 || self.integer_order == Some(0) {
            return Ok(MpComplex::one(self.bits));
        }
        if let Some(order) = self.integer_order {
            let exact = BigInt::from(n).pow(order.unsigned_abs() as u32);
            let weight = if order < 0 {
                BigRational::from_integer(exact)
            } else {
                BigRational::new(1.into(), exact)
            };
            return Ok(MpComplex::from_real(rational(&weight, self.bits)?));
        }
        let logarithm = integer_log(n, self.bits)?;
        multiply_real(&self.negated_order, &logarithm, self.bits)?.exp(self.bits)
    }
}

fn lower_order(
    value: &MpComplex,
    minimum: i64,
    bits: u32,
) -> Result<(Real, u32), String> {
    let low = value.re.outward_endpoints(bits)?.0;
    if low.cmp_value(&integer(-minimum)?) < 0 {
        return Err(format!(
            "complex Dirichlet real order must be at least {}",
            -minimum
        ));
    }
    let negative = if low.sign() < 0 {
        let truncate = i64::try_from(low.trunc_to_bigint()?).map_err(|error| error.to_string())?;
        let floor = if low.cmp_value(&integer(truncate)?) < 0 {
            truncate - 1
        } else {
            truncate
        };
        u32::try_from(-floor).map_err(|error| error.to_string())?
    } else {
        0
    };
    Ok((low, negative))
}

fn clip_ball(
    ball: &Ball,
    maximum: &Real,
    bits: u32,
) -> Result<Ball, String> {
    let (mut lo, mut hi) = ball.outward_endpoints(bits)?;
    let min = maximum.neg();
    if lo.cmp_value(&min) < 0 {
        lo = min;
    }
    if hi.cmp_value(maximum) > 0 {
        hi = maximum.clone();
    }
    Ball::from_bounds(lo, hi, bits)
}

pub(super) fn polylog(
    argument: &MpComplex,
    order: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = requested(precision)?;
    let inner = bits + 40;
    let modulus = magnitude_upper(argument, inner)?;
    if modulus.is_zero() {
        return Ok(MpComplex::zero(bits));
    }
    let one = Real::one(1);
    if modulus.cmp_value(&one) >= 0 {
        return Err(
            "complex polylog requires the argument rectangle strictly inside the unit disk".into(),
        );
    }
    match signed_order(order) {
        | Some(0) => return argument.div(&MpComplex::one(inner).sub(argument, inner)?, bits),
        | Some(1) => return Ok(MpComplex::one(inner).sub(argument, inner)?.ln(bits)?.neg()),
        | _ => {},
    }

    let (_, negative) = lower_order(order, MAX_POLYLOG_NEGATIVE_ORDER, inner)?;
    let negative = integer(i64::from(negative))?;
    let threshold = modulus.scale_pow2(-i64::from(bits) - 8)?;
    let real_axis = is_zero(&argument.im) && is_zero(&order.im);
    let weights = Weights::new(order, inner);
    let mut power = argument.clone();
    let mut power_upper = modulus.clone();
    let mut sum = MpComplex::zero(inner);
    for k in 1..=MAX_POLYLOG_TERMS {
        let coefficient = weights.at(k)?;
        let k_real = integer(k)?;
        // (1+1/k)^m <= 1/(1-m/k) for an integer m >= -Re(order)
        // and k > m. The upper bound on the ratio decreases with k.
        let denominator = one.sub(
            &negative.div(&k_real, inner, Round::TowardPositive)?,
            inner,
            Round::TowardNegative,
        )?;
        if denominator.cmp_value(&modulus) > 0 {
            let ratio = modulus.div(&denominator, inner, Round::TowardPositive)?;
            let geometric = one.sub(&ratio, inner, Round::TowardNegative)?;
            if geometric.sign() > 0 {
                let tail = power_upper
                    .mul(
                        &magnitude_upper(&coefficient, inner)?,
                        inner,
                        Round::TowardPositive,
                    )?
                    .div(&geometric, inner, Round::TowardPositive)?;
                if tail.cmp_value(&threshold) <= 0 {
                    add_error(&mut sum, &tail, real_axis, false, inner)?;
                    return round(sum, bits);
                }
            }
        }
        sum = sum.add(&power.mul(&coefficient, inner)?, inner)?;
        power_upper = power_upper.mul(&modulus, inner, Round::TowardPositive)?;
        let raw = power.mul(argument, inner)?;
        power = MpComplex::new(
            clip_ball(&raw.re, &power_upper, inner)?,
            clip_ball(&raw.im, &power_upper, inner)?,
        );
    }
    Err("complex polylog series exhausted its certified work limit".into())
}

fn exact_zeta(
    order: i64,
    bits: u32,
) -> Result<Option<MpComplex>, String> {
    if order < 0 && order % 2 == 0 {
        return Ok(Some(MpComplex::zero(bits)));
    }
    if !(-127..=128).contains(&order) || order == 1 || (order > 0 && order % 2 != 0) {
        return Ok(None);
    }
    if order == 0 {
        let half = Ball::point(Real::one(bits).scale_pow2(-1)?.neg());
        return Ok(Some(MpComplex::from_real(half)));
    }
    let degree = if order < 0 {
        1 - order
    } else {
        order
    };
    let bernoulli =
        bernoulli_even(u32::try_from(degree / 2 + 1).map_err(|error| error.to_string())?)?
            .pop()
            .ok_or_else(|| "missing exact Bernoulli coefficient".to_owned())?;
    if order < 0 {
        return Ok(Some(MpComplex::from_real(rational(
            &(-bernoulli / BigInt::from(degree)),
            bits,
        )?)));
    }
    let inner = bits + 40;
    let two_pi = consts::named_ball("pi", inner)?.mul(&Ball::point(integer(2)?), inner)?;
    let factorial: BigInt = (2..=degree).map(BigInt::from).product();
    let coefficient = rational(&(bernoulli / (2 * factorial)), inner)?.abs(inner)?;
    let power = MpComplex::from_real(two_pi).powi(degree, inner)?;
    Ok(Some(MpComplex::from_real(
        power.re.mul(&coefficient, bits)?,
    )))
}

pub(super) fn zeta(
    value: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = requested(precision)?;
    let singularity = value.sub(&MpComplex::one(bits + 40), bits + 40)?;
    if singularity.re.contains_zero() && singularity.im.contains_zero() {
        return Err("complex zeta rectangle meets the pole at one".into());
    }
    if let Some(order) = signed_order(value) {
        if let Some(exact) = exact_zeta(order, bits)? {
            return Ok(exact);
        }
    }
    let (sigma, negative) = lower_order(value, MAX_ZETA_NEGATIVE_ORDER, bits + 40)?;
    let inner = bits + 40 + 12 * negative;
    let size = magnitude_upper(value, inner)?;
    if size.cmp_value(&integer(MAX_ZETA_MAGNITUDE)?) > 0 {
        return Err("complex zeta requires an order rectangle inside the modulus-256 disk".into());
    }
    let size = u32::try_from(size.trunc_to_bigint()?).map_err(|error| error.to_string())?;
    let count = (inner + 16).div_ceil(5);
    let split = (size + 2 * count + 2).next_power_of_two();
    let split_log2 = i64::from(split.ilog2());
    let split = i64::from(split);
    let weights = Weights::new(value, inner);
    let tail_at_split = weights.at(split)?;
    let integral = scale(&tail_at_split, split_log2)?.div(&singularity, inner)?;
    let mut result = integral.add(&scale(&tail_at_split, -1)?, inner)?;
    for k in 1..split {
        result = result.add(&weights.at(k)?, inner)?;
    }
    let coeffs = bernoulli_even(count + 1)?;
    let mut poch = value.clone();
    let mut factorial = BigInt::from(1);
    for (coefficient, degree) in coeffs.into_iter().zip((2_i64..).step_by(2)) {
        factorial *= degree * (degree - 1);
        let coefficient = rational(&(coefficient / &factorial), inner)?;
        let correction = scale(
            &poch.mul(&tail_at_split, inner)?,
            -(degree - 1) * split_log2,
        )?;
        result = result.add(&multiply_real(&correction, &coefficient, inner)?, inner)?;
        poch = poch
            .mul(&value.add(&scalar(degree - 1)?, inner)?, inner)?
            .mul(&value.add(&scalar(degree)?, inner)?, inner)?;
    }

    // DLMF 25.2.9: after the nth even coefficient, bound the odd
    // periodic Bernoulli Fourier series by 4*(2n+1)!/6^(2n+1).
    // Because N >= |s| + 2n, every factor in |(s)_(2n+1)|
    // is <= N. Integrating x^(-Re(s)-2n-1) gives the remaining factor.
    let twice_count = i64::from(2 * count);
    let denominator = sigma.add(&integer(twice_count)?, inner, Round::TowardNegative)?;
    let power_six = BigInt::from(6).pow(2 * count + 1);
    let factor = Real::rounded_ratio(
        &(4 * BigInt::from(split)),
        power_six.magnitude(),
        0,
        inner,
        Round::TowardPositive,
    )?;
    let error = magnitude_upper(&tail_at_split, inner)?
        .mul(&factor, inner, Round::TowardPositive)?
        .div(&denominator, inner, Round::TowardPositive)?;
    add_error(&mut result, &error, is_zero(&value.im), false, inner)?;
    round(result, bits)
}
