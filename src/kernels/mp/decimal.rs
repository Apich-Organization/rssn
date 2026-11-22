//! Decimal text in and out of the binary tower.
//!
//! A decimal literal is read exactly (the decimal expansion of a finite
//! decimal is a dyadic times a power of five, which a [`Real`] holds without
//! loss up to its precision), and a [`Real`] is written with a requested
//! number of significant digits, correctly rounded. A [`Ball`] is written
//! only when both endpoints round to the same text, which is the certificate
//! the `mp` operators hand back: the digits are then known to be right, not
//! merely plausible.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::Signed;
use num_traits::Zero;

use super::ball::Ball;
use super::real::Real;
use super::real::Round;

/// Largest number of significant decimal digits an operator accepts.
pub const MAX_DIGITS: usize = 10_000;

/// Bits that hold `digits` decimal digits, with guard bits.
#[must_use]
pub fn bits_for_digits(digits: usize) -> u32 {
    let needed = (digits as f64 * 3.321_928_094_f64).ceil() as u64 + 16;
    u32::try_from(needed).unwrap_or(u32::MAX)
}

/// Rounds `numerator / denominator` to nearest, ties to even, for a
/// non-negative numerator and a positive denominator.
fn round_ratio(
    numerator: &BigInt,
    denominator: &BigUint,
    negative: bool,
    mode: Round,
) -> BigInt {
    if denominator.is_zero() {
        return BigInt::zero();
    }
    let magnitude = numerator.magnitude();
    let (quotient, rest) = magnitude.div_rem(denominator);
    let twice = &rest << 1_u32;
    let up = match mode {
        | Round::NearestEven => twice > *denominator || (twice == *denominator && quotient.bit(0)),
        | Round::TowardZero => false,
        | Round::AwayFromZero => !rest.is_zero(),
        | Round::TowardPositive => !negative && !rest.is_zero(),
        | Round::TowardNegative => negative && !rest.is_zero(),
    };
    let quotient = if up {
        quotient + BigUint::from(1_u32)
    } else {
        quotient
    };
    BigInt::from_biguint(
        if negative {
            num_bigint::Sign::Minus
        } else {
            num_bigint::Sign::Plus
        },
        quotient,
    )
}

/// `round(|x| * 10^k)` as an integer, ties to even.
fn round_scaled(
    x: &Real,
    k: i64,
    mode: Round,
) -> Result<BigInt, String> {
    let mantissa = x.mantissa().abs();
    let mut numerator = mantissa.clone();
    let mut denominator = BigUint::from(1_u32);
    if k >= 0 {
        numerator *= BigInt::from(10_u32).pow(u32::try_from(k).unwrap_or(0));
    } else {
        denominator = BigUint::from(10_u32).pow(u32::try_from(-k).unwrap_or(0));
    }
    let exponent = x.exponent();
    if exponent >= 0 {
        numerator <<= usize::try_from(exponent).unwrap_or(0);
    } else {
        denominator <<= usize::try_from(-exponent).unwrap_or(0);
    }
    Ok(round_ratio(&numerator, &denominator, x.sign() < 0, mode))
}

/// `x` as decimal text with `digits` significant digits, correctly rounded.
///
/// # Errors
/// Fails for fewer than one digit or when an exponent leaves the range.
pub fn to_decimal(
    x: &Real,
    digits: usize,
) -> Result<String, String> {
    to_decimal_in(x, digits, Round::NearestEven)
}

/// `x` as decimal text with `digits` significant digits, rounded in `mode`.
///
/// # Errors
/// Fails for fewer than one digit or when an exponent leaves the range.
pub fn to_decimal_in(
    x: &Real,
    digits: usize,
    mode: Round,
) -> Result<String, String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    if x.is_zero() {
        return Ok("0".to_string());
    }
    let negative = x.sign() < 0;
    let magnitude = x.abs();
    let top = magnitude.floor_log2().unwrap_or(0);
    // A first estimate of the number of integer digits, corrected below.
    let mut place = (top as f64 * std::f64::consts::LOG10_2).floor() as i64 + 1;
    let wanted = i64::try_from(digits).unwrap_or(i64::MAX);
    let mut text = String::new();
    for _ in 0..4 {
        let scaled = round_scaled(x, wanted - place, mode)?;
        let rendered = scaled.magnitude().to_string();
        let length = i64::try_from(rendered.len()).unwrap_or(i64::MAX);
        if length > wanted {
            place = place.saturating_add(length - wanted);
            continue;
        }
        if length < wanted {
            place = place.saturating_sub(wanted - length);
            continue;
        }
        text = render(&rendered, place);
        break;
    }
    if text.is_empty() {
        return Err("could not place the decimal point".to_string());
    }
    if negative {
        text.insert(0, '-');
    }
    Ok(trim_zeros(text))
}

/// Drops the trailing zeros of a fractional part, and the point they leave.
///
/// The digits are exact, so `3.0000` and `3` are the same number; keeping the
/// zeros would only pad every exact answer to the requested width.
fn trim_zeros(text: String) -> String {
    if !text.contains('.') {
        return text;
    }
    let trimmed = text.trim_end_matches('0');
    trimmed.strip_suffix('.').unwrap_or(trimmed).to_string()
}

/// Places the decimal point of a digit string whose first digit is at
/// `10^(place - 1)`.
fn render(
    digits: &str,
    place: i64,
) -> String {
    let length = i64::try_from(digits.len()).unwrap_or(0);
    if place <= 0 {
        let zeros = usize::try_from(-place).unwrap_or(0);
        return format!("0.{}{digits}", "0".repeat(zeros));
    }
    if place >= length {
        let zeros = usize::try_from(place - length).unwrap_or(0);
        return format!("{digits}{}", "0".repeat(zeros));
    }
    let split = usize::try_from(place).unwrap_or(0);
    format!("{}.{}", &digits[..split], &digits[split..])
}

/// The exact value of decimal text, rounded to `prec` bits.
///
/// Accepts an optional sign, digits with an optional decimal point, and an
/// optional `e` or `E` exponent. Spaces around the text are ignored.
///
/// # Errors
/// Fails on malformed text, an empty mantissa, or an exponent outside the
/// representable range.
pub fn from_decimal(
    text: &str,
    prec: u32,
) -> Result<Real, String> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err("empty number".to_string());
    }
    let (negative, body) = match trimmed.strip_prefix('-') {
        | Some(rest) => (true, rest),
        | None => (false, trimmed.strip_prefix('+').unwrap_or(trimmed)),
    };
    let (mantissa_text, exponent_text) = match body.find(['e', 'E']) {
        | Some(at) => (&body[..at], Some(&body[at + 1..])),
        | None => (body, None),
    };
    let mut digits = String::new();
    let mut fraction = 0_i64;
    let mut seen_point = false;
    for character in mantissa_text.chars() {
        match character {
            | '.' if !seen_point => seen_point = true,
            | '.' => return Err("more than one decimal point".to_string()),
            | '0'..='9' => {
                digits.push(character);
                if seen_point {
                    fraction = fraction.saturating_add(1);
                }
            },
            | '_' | ' ' => {},
            | _ => return Err(format!("`{character}` is not a digit")),
        }
    }
    if digits.is_empty() {
        return Err("no digits before the exponent".to_string());
    }
    let exponent = match exponent_text {
        | Some(raw) => {
            raw.trim()
                .parse::<i64>()
                .map_err(|_| "malformed exponent".to_string())?
        },
        | None => 0,
    };
    let value = BigInt::parse_bytes(digits.as_bytes(), 10)
        .ok_or_else(|| "malformed mantissa".to_string())?;
    let power = exponent
        .checked_sub(fraction)
        .ok_or_else(|| "exponent out of range".to_string())?;
    if power.unsigned_abs() > u64::try_from(MAX_DIGITS).unwrap_or(u64::MAX) * 4 {
        return Err("exponent out of range".to_string());
    }
    let signed = if negative { -value } else { value };
    // value * 10^power = value * 5^power * 2^power
    if power >= 0 {
        let five = BigInt::from(5_u32).pow(u32::try_from(power).unwrap_or(0));
        Real::rounded(signed * five, power, prec, Round::NearestEven)
    } else {
        let five = BigInt::from(5_u32).pow(u32::try_from(-power).unwrap_or(0));
        let numerator = Real::rounded(signed, power, prec.saturating_add(32), Round::NearestEven)?;
        let denominator = Real::from_bigint(five, prec.saturating_add(32), Round::NearestEven)?;
        numerator.div(&denominator, prec, Round::NearestEven)
    }
}

/// The digits `ball` pins down, or an error when its radius straddles a
/// rounding boundary.
///
/// # Errors
/// Fails when the two endpoints round to different text.
pub fn ball_to_decimal(
    ball: &Ball,
    digits: usize,
) -> Result<String, String> {
    let working = bits_for_digits(digits).saturating_add(64);
    let (low, high) = ball.endpoints(working)?;
    let lower = to_decimal(&low, digits)?;
    let upper = to_decimal(&high, digits)?;
    if lower == upper {
        return Ok(lower);
    }
    Err(format!(
        "the enclosure [{lower}, {upper}] does not determine {digits} digits"
    ))
}

/// The decimal exponent of `x`, `floor(log10(|x|))`, or `None` for zero.
#[must_use]
pub fn decimal_exponent(x: &Real) -> Option<i64> {
    let top = x.floor_log2()?;
    let estimate = (top as f64 * std::f64::consts::LOG10_2).floor() as i64;
    let scaled = round_scaled(&x.abs(), -estimate, Round::NearestEven).ok()?;
    let length = i64::try_from(scaled.magnitude().to_string().len()).unwrap_or(1);
    Some(estimate.saturating_add(length).saturating_sub(1))
}

/// The number of significant decimal digits `x` needs to be exact, capped at
/// [`MAX_DIGITS`].
#[must_use]
pub fn significant_digits(x: &Real) -> usize {
    let Some(top) = x.floor_log2() else {
        return 1;
    };
    let bits = (top as f64 * std::f64::consts::LOG10_2).abs().ceil() as usize + 2;
    bits.min(MAX_DIGITS).max(1)
}
