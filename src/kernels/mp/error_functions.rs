//! Outward bounds for erf and erfc of closed complex rectangles.

use super::super::ball::Ball;
use super::super::ball::midpoint_rounding;
use super::super::complex::MpComplex;
use super::super::consts;
use super::super::real::MAX_BITS;
use super::super::real::MAX_EXP;
use super::super::real::Real;
use super::super::real::Round;
use super::add_error;
use super::add_radius;
use super::divide_real;
use super::is_zero;

const MAX_SERIES_TERMS: u32 = 8192;
const MAX_SERIES_MODULUS_SQUARED: i64 = 1024;
const GUARD_BITS: u32 = 32;

fn integer(value: i64) -> Result<Real, String> {
    Real::from_i64(value, 16)
}

fn power_of_two(exponent: i64) -> Result<Real, String> {
    Real::rounded(1.into(), exponent.max(-MAX_EXP), 1, Round::TowardPositive)
}

fn sqrt_pi(prec: u32) -> Result<Ball, String> {
    consts::named_ball("pi", prec)?.sqrt(prec)
}

fn square(
    value: &Real,
    prec: u32,
    round: Round,
) -> Result<Real, String> {
    value.mul(value, prec, round)
}

fn complement(
    value: MpComplex,
    constant: i64,
    prec: u32,
) -> Result<MpComplex, String> {
    let re = value.re.neg();
    let (mid, inexact) = re
        .mid
        .add_with_inexact(&integer(constant)?, prec, Round::NearestEven)?;
    let rounding = midpoint_rounding(&mid, inexact)?;
    let mut real = Ball::new(mid, re.rad)?;
    if !rounding.is_zero() {
        add_radius(&mut real, &rounding, prec)?;
    }
    Ok(MpComplex::new(real, value.im.neg().rounded_midpoint(prec)?))
}

fn taylor(
    z: &MpComplex,
    max_squared: &Real,
    prec: u32,
) -> Result<MpComplex, String> {
    if max_squared.cmp_value(&integer(MAX_SERIES_MODULUS_SQUARED)?) > 0 {
        return Err(
            "complex error function requires a tighter rectangle or larger right-sector argument"
                .into(),
        );
    }
    // Reserve for cancellation on the positive real side, including the
    // subtraction in erfc. The analytic remainder is still accounted for
    // if interval dependency requires a wider output.
    let size = u32::try_from(max_squared.trunc_to_bigint()?).map_err(|e| e.to_string())?;
    let inner = prec
        .saturating_add(GUARD_BITS + 4 * (size + 1))
        .min(MAX_BITS);
    let threshold = power_of_two(8 - i64::from(inner))?;
    let negative_square = z.mul(z, inner)?.neg();
    let mut term = z.clone();
    let mut sum = z.clone();
    let max_ratio_test = max_squared.scale_pow2(1)?;
    let real_axis = is_zero(&z.im);
    let imaginary_axis = is_zero(&z.re);
    for n in 1..=MAX_SERIES_TERMS {
        let n = i64::from(n);
        term = divide_real(
            &term.mul(&negative_square, inner)?,
            &Ball::point(integer(n)?),
            inner,
        )?;
        let contribution = divide_real(&term, &Ball::point(integer(2 * n + 1)?), inner)?;

        // NIST DLMF 7.6.1. From term n onward the consecutive ratios
        // of magnitudes are bounded by |z|^2/(n+1). Once this is <= 1/2,
        // twice a bound for the first omitted term encloses the whole tail.
        if max_ratio_test.cmp_value(&integer(n + 1)?) <= 0 {
            let (_, real) = contribution.re.absolute_endpoints(inner)?;
            let (_, imag) = contribution.im.absolute_endpoints(inner)?;
            let error = real
                .add(&imag, inner, Round::TowardPositive)?
                .scale_pow2(1)?;
            if error.cmp_value(&threshold) <= 0 {
                add_error(&mut sum, &error, real_axis, imaginary_axis, inner)?;
                let sqrt_pi = sqrt_pi(inner)?;
                let two = Ball::point(integer(2)?);
                let scale = two.div(&sqrt_pi, inner)?;
                return Ok(MpComplex::new(
                    sum.re.mul(&scale, inner)?,
                    sum.im.mul(&scale, inner)?,
                ));
            }
        }
        sum = sum.add(&contribution, inner)?;
    }
    Err("complex error function series exhausted its certified work limit".into())
}

fn right_sector_erfc(
    z: &MpComplex,
    min_square: &Real,
    prec: u32,
) -> Result<MpComplex, String> {
    let inner = prec.saturating_add(GUARD_BITS).min(MAX_BITS);
    let threshold = power_of_two(-i64::from(prec) - 12)?;
    let inverse = z.recip(inner)?;
    let half_inverse_square = divide_real(
        &inverse.mul(&inverse, inner)?,
        &Ball::point(integer(2)?),
        inner,
    )?
    .neg();
    let largest_ratio =
        Real::one(1).div(&min_square.scale_pow2(1)?, inner, Round::TowardPositive)?;

    let mut term = MpComplex::one(inner);
    let mut sum = term.clone();
    let mut error = Real::one(1);
    for n in 1..=MAX_SERIES_TERMS {
        let odd = integer(2 * i64::from(n) - 1)?;
        let ratio = largest_ratio.mul(&odd, inner, Round::TowardPositive)?;
        if ratio.cmp_value(&Real::one(1)) >= 0 {
            break;
        }

        // NIST DLMF 7.12.1. In |arg z| <= pi/4 the complement's
        // scaled asymptotic remainder is bounded by the first omitted term.
        error = error.mul(&ratio, inner, Round::TowardPositive)?;
        if error.cmp_value(&threshold) <= 0 {
            add_error(&mut sum, &error, is_zero(&z.im), false, inner)?;
            let exponential = z.mul(z, inner)?.neg().exp(inner)?;
            let product = exponential.mul(&inverse, inner)?.mul(&sum, inner)?;
            return divide_real(&product, &sqrt_pi(inner)?, prec);
        }
        let product = term.mul(&half_inverse_square, inner)?;
        let scale = Ball::point(odd);
        term = MpComplex::new(
            product.re.mul(&scale, inner)?,
            product.im.mul(&scale, inner)?,
        );
        sum = sum.add(&term, inner)?;
    }
    Err("complex complementary error function asymptotic budget could not certify the tail".into())
}

fn negligible_erfc(
    real_axis: bool,
    prec: u32,
) -> Result<MpComplex, String> {
    let radius = power_of_two(-MAX_EXP)?;
    Ok(MpComplex::new(
        Ball::new(Real::zero(prec), radius.clone())?,
        if real_axis {
            Ball::point(Real::zero(prec))
        } else {
            Ball::new(Real::zero(prec), radius)?
        },
    ))
}

pub(super) fn evaluate(
    z: &MpComplex,
    erfc: bool,
    precision: u32,
) -> Result<MpComplex, String> {
    let prec = precision.clamp(1, MAX_BITS);
    let inner = prec.saturating_add(GUARD_BITS).min(MAX_BITS);
    let reflect = z.re.mid.sign() < 0;
    let reflected;
    let z = if reflect {
        reflected = z.neg();
        &reflected
    } else {
        z
    };
    let (real_low, _) = z.re.outward_endpoints(inner)?;
    let (_, imag_far) = z.im.absolute_endpoints(inner)?;
    if is_zero(&z.re) && is_zero(&z.im) {
        return Ok(if erfc {
            MpComplex::one(prec)
        } else {
            MpComplex::zero(prec)
        });
    }
    let real_axis = is_zero(&z.im);
    let imaginary_axis = is_zero(&z.re);
    let (_, real_far) = z.re.absolute_endpoints(inner)?;
    let max_coordinate = if real_far.cmp_value(&imag_far) >= 0 {
        &real_far
    } else {
        &imag_far
    };
    if max_coordinate
        .floor_log2()
        .is_some_and(|top| 2 * top < i64::from(inner) - MAX_EXP - 4)
    {
        // Along the straight segment from 0 to z, |exp(-t^2)| < 2
        // for these tiny inputs. The integral then bounds |erf(z)|
        // by 8 times the larger coordinate magnitude without squaring it.
        let radius = max_coordinate.scale_pow2(3)?;
        let mut enclosure = if erfc {
            MpComplex::one(prec)
        } else {
            MpComplex::zero(prec)
        };
        add_error(&mut enclosure, &radius, real_axis, imaginary_axis, prec)?;
        return Ok(enclosure);
    }
    let mut asymptotic = None;
    if real_low.sign() > 0 && real_low.cmp_value(&imag_far) >= 0 {
        // On the right, shifting the defining integration contour gives
        // |erfc(x+iy)| <= exp(y^2-x^2)/(x sqrt(pi)). For x >= 512
        // and |y| <= x/2 this is already below the least positive dyadic.
        if real_low.cmp_value(&integer(512)?) >= 0
            && imag_far.cmp_value(&real_low.scale_pow2(-1)?) <= 0
        {
            asymptotic = Some(negligible_erfc(real_axis, prec)?);
        } else {
            let real_square = square(&real_low, inner, Round::TowardNegative)?;
            let imag_square = square(&imag_far, inner, Round::TowardPositive)?;
            let decay = real_square.sub(&imag_square, inner, Round::TowardNegative)?;
            // ln(2) < 7/10. Above this threshold exp(-decay) is
            // smaller than the least positive dyadic even at one bit.
            let underflow_decay = (7 * MAX_EXP + 9) / 10;
            if real_low.cmp_value(&Real::one(1)) >= 0
                && decay.cmp_value(&Real::from_i64(underflow_decay, 24)?) >= 0
            {
                asymptotic = Some(negligible_erfc(real_axis, prec)?);
            } else if real_square.cmp_value(&integer(i64::from(prec) + 24)?) >= 0 {
                asymptotic = Some(right_sector_erfc(z, &real_square, prec)?);
            }
        }
    }

    let result = if let Some(complementary) = asymptotic {
        if erfc {
            complementary
        } else {
            complement(complementary, 1, prec)?
        }
    } else {
        let max_squared = square(&real_far, inner, Round::TowardPositive)?.add(
            &square(&imag_far, inner, Round::TowardPositive)?,
            inner,
            Round::TowardPositive,
        )?;
        let ordinary = taylor(z, &max_squared, prec)?;
        if erfc {
            complement(ordinary, 1, prec)?
        } else {
            MpComplex::new(
                ordinary.re.rounded_midpoint(prec)?,
                ordinary.im.rounded_midpoint(prec)?,
            )
        }
    };
    if reflect {
        if erfc {
            complement(result, 2, prec)
        } else {
            Ok(result.neg())
        }
    } else {
        Ok(result)
    }
}
