//! Certified Bessel J over a bounded pair of dyadic rectangles.

use num_bigint::BigInt;
use num_rational::BigRational;

use super::super::ball::Ball;
use super::super::complex::MpComplex;
use super::super::real::MAX_BITS;
use super::super::real::Real;
use super::super::real::Round;
use super::add_error;
use super::divide_real;
use super::gamma_family::rational;
use super::gamma_family::reciprocal_gamma;
use super::is_zero;

const MAX_INTEGER_ORDER: u64 = 512;
const MAX_GENERAL_BITS: u32 = 512;
const MAX_SERIES_TERMS: i64 = 8192;
const MAX_HALF_ARGUMENT: i64 = 128;

fn integer(n: i64) -> Result<Real, String> {
    Real::from_i64(n, 32)
}

fn modulus_l1_upper(
    value: &MpComplex,
    bits: u32,
) -> Result<Real, String> {
    let (_, real) = value.re.absolute_endpoints(bits)?;
    let (_, imag) = value.im.absolute_endpoints(bits)?;
    real.add(&imag, bits, Round::TowardPositive)
}

fn half_coordinate(value: &Ball) -> Result<Ball, String> {
    Ball::new(value.mid.scale_pow2(-1)?, value.rad.scale_pow2(-1)?)
}

fn normalized_series(
    half_value: &MpComplex,
    order: &MpComplex,
    minimum_order: &Real,
    modulus: &Real,
    bits: u32,
) -> Result<MpComplex, String> {
    let negative_square = half_value.mul(half_value, bits)?.neg();
    let ratio_numerator = modulus
        .mul(modulus, bits, Round::TowardPositive)?
        .scale_pow2(1)?;
    let threshold = Real::one(1).scale_pow2(8 - i64::from(bits))?;
    let real_axis = is_zero(&order.im) && (is_zero(&half_value.re) || is_zero(&half_value.im));
    let mut term = MpComplex::one(bits);
    let mut sum = term.clone();
    for k in 1..=MAX_SERIES_TERMS {
        let k_real = integer(k)?;
        let numerator = term.mul(&negative_square, bits)?;
        let order_plus_k = order.add(&MpComplex::from_real(Ball::point(k_real.clone())), bits)?;
        term = if is_zero(&order.im) {
            divide_real(
                &numerator,
                &order_plus_k.re.mul(&Ball::point(k_real), bits)?,
                bits,
            )?
        } else {
            numerator
                .div(&order_plus_k, bits)
                .and_then(|value| divide_real(&value, &Ball::point(k_real), bits))?
        };

        // DLMF 10.2.2. The ratio after this first omitted term is at
        // most |z/2|^2/((k+1)*(k+1+Re(order)_low)) and decreases in k.
        let next = integer(k + 1)?;
        let denominator = minimum_order.add(&next, bits, Round::TowardNegative)?.mul(
            &next,
            bits,
            Round::TowardNegative,
        )?;
        if ratio_numerator.cmp_value(&denominator) <= 0 {
            let tail = modulus_l1_upper(&term, bits)?.scale_pow2(1)?;
            if tail.cmp_value(&threshold) <= 0 {
                add_error(&mut sum, &tail, real_axis, false, bits)?;
                return Ok(sum);
            }
        }
        sum = sum.add(&term, bits)?;
    }
    Err("complex Bessel J series exhausted its certified work limit".into())
}

pub(super) fn evaluate(
    value: &MpComplex,
    order: &MpComplex,
    precision: u32,
) -> Result<MpComplex, String> {
    let bits = precision.clamp(1, MAX_BITS);
    let integer_order = if is_zero(&order.im) && order.re.rad.is_zero() && order.re.mid.is_integer()
    {
        Some(order.re.mid.integer_exponent()?)
    } else {
        None
    };
    if integer_order.is_some_and(|n| n.unsigned_abs() > MAX_INTEGER_ORDER) {
        return Err(format!(
            "complex Bessel J supports integer order magnitude up to {MAX_INTEGER_ORDER}"
        ));
    }
    if integer_order.is_none() && bits > MAX_GENERAL_BITS {
        return Err(format!(
            "complex Bessel J supports noninteger order up to {MAX_GENERAL_BITS} bits"
        ));
    }
    let order_lower = order.re.outward_endpoints(bits.saturating_add(32))?.0;
    if integer_order.is_none() && order_lower.sign() < 0 {
        return Err("complex Bessel J requires a nonnegative real order enclosure or an exact integer order".into());
    }
    if is_zero(&value.re) && is_zero(&value.im) {
        if integer_order == Some(0) {
            return Ok(MpComplex::one(bits));
        }
        if integer_order.is_some() || order_lower.sign() > 0 {
            return Ok(MpComplex::zero(bits));
        }
        return Err(
            "complex Bessel J at zero requires positive real order or exact integer order".into(),
        );
    }

    let half_value = MpComplex::new(half_coordinate(&value.re)?, half_coordinate(&value.im)?);
    let modulus = modulus_l1_upper(&half_value, bits.saturating_add(32))?;
    if modulus.cmp_value(&integer(MAX_HALF_ARGUMENT)?) > 0 {
        return Err("complex Bessel J requires |Re(value)| + |Im(value)| <= 256".into());
    }
    let size = u32::try_from(modulus.trunc_to_bigint()?).map_err(|error| error.to_string())?;
    let inner = bits.saturating_add(32 + 4 * (size + 1)).min(MAX_BITS);

    if let Some(n) = integer_order {
        let absolute = n.abs();
        let minimum = integer(absolute)?;
        let order = MpComplex::from_real(Ball::point(minimum.clone()));
        let factorial: BigInt = (2..=absolute).map(BigInt::from).product();
        let reciprocal = rational(&BigRational::new(1.into(), factorial), inner)?;
        let factor = half_value.powi(absolute, inner)?;
        let sum = normalized_series(&half_value, &order, &minimum, &modulus, inner)?;
        let product = factor.mul(&sum, inner)?;
        let result = MpComplex::new(
            product.re.mul(&reciprocal, bits)?,
            product.im.mul(&reciprocal, bits)?,
        );
        return Ok(if n < 0 && absolute % 2 != 0 {
            result.neg()
        } else {
            result
        });
    }

    let power = half_value.pow(order, inner)?;
    let reciprocal = reciprocal_gamma(&order.add(&MpComplex::one(inner), inner)?, inner)?;
    let sum = normalized_series(&half_value, order, &order_lower, &modulus, inner)?;
    power.mul(&reciprocal, inner)?.mul(&sum, bits)
}
