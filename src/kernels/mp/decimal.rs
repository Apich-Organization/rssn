//! Decimal text in and out of the binary tower.
//!
//! Decimal literals are parsed as integers scaled by powers of ten. A scalar
//! literal is rounded to binary precision and an expression literal can be
//! enclosed using directed rounding. Decimal text can also be rounded directly
//! to significant digits. A [`Real`] is written with a requested
//! number of significant digits, correctly rounded. A [`Ball`] is written
//! only when both endpoints round to the same text, which is the certificate
//! the `mp` operators hand back: the digits are then known to be right, not
//! merely plausible.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::One;
use num_traits::Zero;
use std::cmp::Ordering;

use super::ball::Ball;
use super::real::MAX_EXP;
use super::real::Real;
use super::real::Round;

/// Largest number of significant decimal digits an operator accepts.
pub const MAX_DIGITS: usize = 10_000;

/// The input limits cover fixed-point formatting of the entire dyadic range.
const MAX_DECIMAL_POWER: u64 = 100_000;
pub(crate) const MAX_DECIMAL_TEXT: usize = 100_000;

/// Bits that hold `digits` decimal digits, with guard bits.
#[must_use]
pub fn bits_for_digits(digits: usize) -> u32 {
    let bits = (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32;
    bits.saturating_add(16)
}

/// Whether to increment the magnitude after discarding lower digits.
fn increment_magnitude(
    mode: Round,
    negative: bool,
    inexact: bool,
    half: Ordering,
    odd: bool,
) -> bool {
    match mode {
        | Round::NearestEven => half == Ordering::Greater || (half == Ordering::Equal && odd),
        | Round::TowardZero => false,
        | Round::AwayFromZero => inexact,
        | Round::TowardPositive => !negative && inexact,
        | Round::TowardNegative => negative && inexact,
    }
}

pub(super) fn power_of_five(power: u64) -> Result<BigUint, String> {
    let power = u32::try_from(power).map_err(|_| "exponent out of range".to_string())?;
    Ok(BigUint::from(5_u32).pow(power))
}

/// Exact nonnegative ratio for `|x| * 10^k`.
fn scaled_ratio(
    x: &Real,
    k: i64,
) -> Result<(BigUint, BigUint), String> {
    let mut numerator = x.mantissa().magnitude().clone();
    let mut denominator = BigUint::one();
    let five = power_of_five(k.unsigned_abs())?;
    if k >= 0 {
        numerator *= five;
    } else {
        denominator = five;
    }
    let exponent = x.exponent() + k;
    if exponent >= 0 {
        numerator <<= exponent.unsigned_abs();
    } else {
        denominator <<= exponent.unsigned_abs();
    }
    Ok((numerator, denominator))
}

/// `x` as decimal text with `digits` significant digits, correctly rounded.
///
/// # Errors
/// Fails when an exponent leaves the range.
pub fn to_decimal(
    x: &Real,
    digits: usize,
) -> Result<String, String> {
    to_decimal_in(x, digits, Round::NearestEven)
}

/// `x` as decimal text with `digits` significant digits, rounded in `mode`.
///
/// # Errors
/// Fails when an exponent leaves the range.
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
    let mut place = estimated_exponent(x).unwrap_or(0) + 1;
    let wanted = i64::try_from(digits).unwrap_or(i64::MAX);
    let (mut numerator, mut denominator) = scaled_ratio(x, wanted - place)?;
    loop {
        let (quotient, rest) = numerator.div_rem(&denominator);
        if quotient.is_zero() {
            numerator *= 10_u32;
            place -= 1;
            continue;
        }
        let mut text = quotient.to_str_radix(10);
        let length = i64::try_from(text.len()).unwrap_or(i64::MAX);
        if length != wanted {
            // Correct the initial exponent estimate without recomputing 10^k.
            let difference = length - wanted;
            let factor = BigUint::from(10_u32).pow(
                u32::try_from(difference.unsigned_abs())
                    .map_err(|_| "exponent out of range".to_string())?,
            );
            if difference > 0 {
                denominator *= factor;
            } else {
                numerator *= factor;
            }
            place += difference;
            continue;
        }
        if increment_magnitude(
            mode,
            negative,
            !rest.is_zero(),
            (&rest << 1_u32).cmp(&denominator),
            quotient.bit(0),
        ) {
            text = (quotient + 1_u32).to_str_radix(10);
            if text.len() > digits {
                // Rounding across a decade produces exactly 10^digits.
                text.truncate(digits);
                place += 1;
            }
        }
        return Ok(format_signed(render(&text, place), negative));
    }
}

fn estimated_exponent(x: &Real) -> Option<i64> {
    let top = x.floor_log2()?;
    Some((top as f64 * std::f64::consts::LOG10_2).floor() as i64)
}

fn format_signed(
    mut text: String,
    negative: bool,
) -> String {
    if negative {
        text.insert(0, '-');
    }
    trim_zeros(text)
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
    from_decimal_in(text, prec, Round::NearestEven)
}

/// Rounds an exact decimal to binary precision in the requested mode.
///
/// # Errors
/// Fails on malformed text or when an input or result exceeds the limits.
pub fn from_decimal_in(
    text: &str,
    prec: u32,
    mode: Round,
) -> Result<Real, String> {
    DecimalRatio::parse(text)?.rounded(prec, mode)
}

/// An exact decimal as a signed integer, a positive denominator and a
/// power of two. Keeping the common binary scale avoids rounding operands.
struct DecimalRatio {
    mantissa: BigInt,
    denominator: BigUint,
    exponent: i64,
}

impl DecimalRatio {
    fn parse(text: &str) -> Result<Self, String> {
        let (negative, digits, exponent) = decimal_digits(text)?;
        let digits = digits.trim_start_matches('0');
        if digits.is_empty() {
            return Ok(Self {
                mantissa: BigInt::zero(),
                denominator: BigUint::one(),
                exponent: 0,
            });
        }
        let significant = digits.trim_end_matches('0');
        let trailing = i64::try_from(digits.len() - significant.len())
            .map_err(|_| "exponent out of range".to_string())?;
        let exponent = exponent + trailing;
        let mut mantissa = BigInt::parse_bytes(significant.as_bytes(), 10)
            .ok_or_else(|| "malformed mantissa".to_string())?;
        if negative {
            mantissa = -mantissa;
        }
        let five = power_of_five(exponent.unsigned_abs())?;
        let denominator = if exponent >= 0 {
            mantissa *= BigInt::from(five);
            BigUint::one()
        } else {
            five
        };
        Ok(Self {
            mantissa,
            denominator,
            exponent,
        })
    }

    fn rounded(
        &self,
        prec: u32,
        mode: Round,
    ) -> Result<Real, String> {
        if self.denominator.is_one() {
            return Real::rounded(self.mantissa.clone(), self.exponent, prec, mode);
        }
        Real::rounded_ratio(&self.mantissa, &self.denominator, self.exponent, prec, mode)
    }

    fn distance(
        &self,
        value: &Real,
        prec: u32,
    ) -> Result<Real, String> {
        let exponent = self.exponent.min(value.exponent());
        let left = &self.mantissa << (self.exponent - exponent).unsigned_abs();
        let right = (value.mantissa() << (value.exponent() - exponent).unsigned_abs())
            * BigInt::from(self.denominator.clone());
        let distance = BigInt::from((left - right).magnitude().clone());
        let round = |bits| {
            Real::rounded_ratio(
                &distance,
                &self.denominator,
                exponent,
                bits,
                Round::TowardPositive,
            )
        };
        // The distance from a representable nearest midpoint is at most
        // one ulp there. At the exponent floor, a coarser radius preserves it.
        round(prec)
            .or_else(|_| round(1))
            .or_else(|_| Real::rounded(1.into(), -MAX_EXP, 1, Round::TowardPositive))
    }
}

/// Parses sign, mantissa digits and their decimal power without rounding.
fn decimal_digits(text: &str) -> Result<(bool, String, i64), String> {
    if text.len() > MAX_DECIMAL_TEXT {
        return Err("decimal text is too long".to_string());
    }
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
    let mut digits = String::with_capacity(mantissa_text.len());
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
    let power = exponent
        .checked_sub(fraction)
        .ok_or_else(|| "exponent out of range".to_string())?;
    if power.unsigned_abs() > MAX_DECIMAL_POWER {
        return Err("exponent out of range".to_string());
    }
    Ok((negative, digits, power))
}

/// Rounds the validated decimal digits before any binary conversion.
pub(super) fn round_decimal(
    text: &str,
    digits: usize,
    mode: Round,
) -> Result<String, String> {
    let (negative, source, power) = decimal_digits(text)?;
    // The decimal text-only rounding route has a smaller legacy limit.
    if power.unsigned_abs() > 4 * MAX_DIGITS as u64 {
        return Err("exponent out of range".to_string());
    }
    let significant = source.trim_start_matches('0');
    if significant.is_empty() {
        return Ok("0".to_string());
    }
    let place = i64::try_from(significant.len())
        .ok()
        .and_then(|length| length.checked_add(power))
        .ok_or_else(|| "exponent out of range".to_string())?;
    let digits = digits.clamp(1, MAX_DIGITS);
    let (kept, discarded) = significant.split_at(significant.len().min(digits));
    let Some((&first, tail)) = discarded.as_bytes().split_first() else {
        return Ok(format_signed(render(kept, place), negative));
    };
    let tail_nonzero = tail.iter().any(|digit| *digit != b'0');
    let half = first.cmp(&b'5').then(if tail_nonzero {
        Ordering::Greater
    } else {
        Ordering::Equal
    });
    let odd = kept.as_bytes()[digits - 1] & 1 != 0;
    if !increment_magnitude(mode, negative, first != b'0' || tail_nonzero, half, odd) {
        return Ok(format_signed(render(kept, place), negative));
    }

    // Only the retained prefix is incremented. A carry across every
    // retained digit grows both its length and the decimal point position.
    let mut rounded = kept.as_bytes().to_vec();
    if let Some(last) = rounded.iter().rposition(|digit| *digit != b'9') {
        rounded[last] += 1;
        rounded[last + 1..].fill(b'0');
    } else {
        rounded.fill(b'0');
        rounded.insert(0, b'1');
    }
    let rounded = String::from_utf8(rounded).map_err(|_| "malformed mantissa".to_string())?;
    let place = place
        .checked_add(i64::from(rounded.len() > kept.len()))
        .ok_or_else(|| "exponent out of range".to_string())?;
    Ok(format_signed(render(&rounded, place), negative))
}

/// Encloses an exact decimal in a dyadic ball whose midpoint has `prec` bits.
///
/// # Errors
/// Fails on malformed text or when an input or result exceeds the limits.
pub fn enclose_decimal(
    text: &str,
    prec: u32,
) -> Result<Ball, String> {
    let decimal = DecimalRatio::parse(text)?;
    let mid = decimal.rounded(prec, Round::NearestEven)?;
    let rad = decimal.distance(&mid, prec)?;
    Ball::new(mid, rad)
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
    if ball.rad.is_zero() {
        return to_decimal(&ball.mid, digits);
    }
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
    let mut estimate = estimated_exponent(x)?;
    loop {
        let (numerator, denominator) = scaled_ratio(x, -estimate).ok()?;
        let truncated = numerator / denominator;
        if !truncated.is_zero() {
            let digits = i64::try_from(truncated.to_str_radix(10).len()).ok()?;
            return Some(estimate + digits - 1);
        }
        estimate -= 1;
    }
}

/// The number of significant digits of the terminating decimal for `x`,
/// capped at [`MAX_DIGITS`]. Zero has one significant digit.
#[must_use]
pub fn significant_digits(x: &Real) -> usize {
    let magnitude = x.mantissa().magnitude();
    if magnitude.is_zero() {
        return 1;
    }
    let exponent = x.exponent();
    let decimal_mantissa = if exponent < 0 {
        let denominator_bits = exponent.unsigned_abs();
        let cancelled = magnitude
            .trailing_zeros()
            .unwrap_or(0)
            .min(denominator_bits);
        let power = denominator_bits - cancelled;
        // The numerator is odd when a denominator remains. Multiplication
        // by 5^power contributes no decimal trailing zeros, and 5^2 > 10.
        if power >= 2 * MAX_DIGITS as u64 {
            return MAX_DIGITS;
        }
        let Ok(five) = power_of_five(power) else {
            return MAX_DIGITS;
        };
        (magnitude >> cancelled) * five
    } else {
        magnitude << exponent.unsigned_abs()
    };
    decimal_mantissa
        .to_str_radix(10)
        .trim_end_matches('0')
        .len()
        .min(MAX_DIGITS)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_precision_saturates_when_the_bit_count_cannot_fit() {
        assert_eq!(bits_for_digits(0), 16);
        assert_eq!(bits_for_digits(MAX_DIGITS), 33_236);
        assert_eq!(bits_for_digits(usize::MAX), u32::MAX);
    }
}
