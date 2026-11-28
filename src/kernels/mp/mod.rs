//! Arbitrary-precision floating point and ball arithmetic.
//!
//! # Where this sits
//!
//! The kernel's exact number is a [`crate::graph::number::Number`]: an integer,
//! a rational, or a double. That covers exact algebra and machine arithmetic,
//! and it deliberately has no third case: a transcendental result has no exact
//! form, and adding one would put a precision policy inside every rule.
//!
//! So precision lives here instead, in a tower with three layers:
//!
//! * [`real`] — binary floating point `mant * 2^exp` with the mantissa of any
//!   width, every operation correctly rounded once from the exact result;
//! * [`funcs`] — the elementary functions, each returning a [`Ball`] whose
//!   radius bounds the series truncation and the rounding, so a result is an
//!   enclosure rather than a guess;
//! * [`decimal`] — exact decimal text in, correctly rounded decimal text out,
//!   and a certificate: a ball is rendered only when both endpoints round to
//!   the same digits.
//!
//! [`parse`] ties them together for `mp_eval`, a small total expression
//! language (`sin(1/3) + pi`), because one string per request is a better
//! operator surface than one operator per elementary function.
//!
//! # How it is reached from a term
//!
//! `rules::numeric::mp` exposes the operators `mp_eval`, `mp_const`,
//! `mp_round`, `mp_compare`, `mp_to_float`, `mp_from_float`, `mp_ball` and the
//! four arithmetic pairs. Their operands and results are **decimal strings**
//! carried as ordinary string literals in the graph, which is the whole point
//! of the design: precision travels *with the value*, so it is hash-consed,
//! cached, cloned and shared across worker threads exactly like a number, and
//! no session state, no global rounding mode and no new Aevia value type is
//! needed. `mp_to_float` is then the single, explicitly requested, correctly
//! rounded narrowing step.
//!
//! # Cost
//!
//! The series are evaluated in fixed point without binary splitting, so the
//! cost is quadratic in the precision: twenty digits take microseconds, a few
//! thousand take about a second. [`decimal::MAX_DIGITS`] caps what an operator
//! accepts, and [`real::MAX_BITS`] caps what the tower can hold.
//!
//! [`real`]: crate::kernels::mp::real
//! [`funcs`]: crate::kernels::mp::funcs
//! [`Ball`]: crate::kernels::mp::Ball
//! [`decimal`]: crate::kernels::mp::decimal
//! [`parse`]: crate::kernels::mp::parse
//! [`decimal::MAX_DIGITS`]: crate::kernels::mp::decimal::MAX_DIGITS
//! [`real::MAX_BITS`]: crate::kernels::mp::real::MAX_BITS

pub mod ball;
pub mod complex;
pub mod consts;
pub mod decimal;
pub(crate) mod fixed;
pub mod funcs;
mod logarithms;
pub mod parse;
pub mod real;
pub mod special_complex;
mod trigonometry;

pub use ball::Ball;
pub use complex::MpComplex;
pub use complex::enclosure_complex_text;
pub use complex::eval_complex_text;
pub use decimal::MAX_DIGITS;
pub use decimal::ball_to_decimal;
pub use decimal::bits_for_digits;
pub use decimal::from_decimal;
pub use decimal::to_decimal;
pub use decimal::to_decimal_in;
pub use real::MAX_BITS;
pub use real::MAX_EXP;
pub use real::Real;
pub use real::Round;

/// The rounding modes an operator accepts, by name.
///
/// # Errors
/// Fails for an unknown name.
pub fn round_named(name: &str) -> Result<Round, String> {
    match name {
        | "nearest" | "nearest_even" | "even" => Ok(Round::NearestEven),
        | "zero" | "trunc" | "truncate" => Ok(Round::TowardZero),
        | "away" | "away_from_zero" => Ok(Round::AwayFromZero),
        | "up" | "toward_positive" | "ceil" => Ok(Round::TowardPositive),
        | "down" | "toward_negative" | "floor" => Ok(Round::TowardNegative),
        | _ => Err(format!("unknown rounding mode `{name}`")),
    }
}

/// Evaluates decimal expression `text` to `digits` significant digits.
///
/// Plain decimal literals are rounded directly. Other expressions use
/// progressively higher precision until the enclosure pins the digits.
///
/// # Errors
/// Fails on malformed text, a domain error, or when the digits cannot be
/// pinned within [`MAX_BITS`].
pub fn eval_text(
    text: &str,
    digits: usize,
) -> Result<String, String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    if let Some(value) = parse::round_literal(text, digits)? {
        return Ok(value);
    }
    let mut prec = bits_for_digits(digits).saturating_add(64);
    let mut last = String::new();
    for _ in 0..5 {
        let ball = parse::evaluate(text, prec)?;
        match ball_to_decimal(&ball, digits) {
            | Ok(rendered) => return Ok(rendered),
            | Err(message) => last = message,
        }
        prec = prec.saturating_mul(2).min(MAX_BITS);
    }
    Err(last)
}

/// Evaluates decimal expression `text` to a rigorous enclosure at `digits`
/// significant digits, rounded outward, as `(low, high)`.
///
/// # Errors
/// Fails on malformed text, a domain error, or when an exponent leaves the
/// range.
pub fn enclosure_text(
    text: &str,
    digits: usize,
) -> Result<(String, String), String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let prec = bits_for_digits(digits).saturating_add(96).min(MAX_BITS);
    let ball = parse::evaluate(text, prec)?;
    enclosure_bounds_text(&ball, prec, digits)
}

fn enclosure_bounds_text(
    ball: &Ball,
    precision: u32,
    digits: usize,
) -> Result<(String, String), String> {
    let (low, high) = ball.endpoints(precision)?;
    Ok((
        to_decimal_in(&low, digits, Round::TowardNegative)?,
        to_decimal_in(&high, digits, Round::TowardPositive)?,
    ))
}

/// A constant by name, to `digits` significant digits.
///
/// # Errors
/// Fails for an unknown name or when an exponent leaves the range.
pub fn const_text(
    name: &str,
    digits: usize,
) -> Result<String, String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let prec = bits_for_digits(digits).saturating_add(64);
    let value = consts::named(name, prec)?;
    to_decimal(&value, digits)
}

/// Applies one of the four arithmetic operations to two decimal literals.
///
/// `operation` is `add`, `sub`, `mul` or `div`.
///
/// # Errors
/// Fails on malformed text, division by zero, an unknown operation, or when an
/// exponent leaves the range.
pub fn arithmetic_text(
    operation: &str,
    left: &str,
    right: &str,
    digits: usize,
    mode: Round,
) -> Result<String, String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let prec = bits_for_digits(digits).saturating_add(64);
    let a = from_decimal(left, prec)?;
    let b = from_decimal(right, prec)?;
    let value = match operation {
        | "add" => a.add(&b, prec, mode),
        | "sub" => a.sub(&b, prec, mode),
        | "mul" => a.mul(&b, prec, mode),
        | "div" => a.div(&b, prec, mode),
        | other => return Err(format!("unknown operation `{other}`")),
    }?;
    to_decimal(&value, digits)
}

/// Compares two decimal literals: `-1`, `0` or `1`.
///
/// # Errors
/// Fails on malformed text or when an exponent leaves the range.
pub fn compare_text(
    left: &str,
    right: &str,
) -> Result<i8, String> {
    let prec = bits_for_digits(MAX_DIGITS).min(MAX_BITS);
    let a = from_decimal(left, prec)?;
    let b = from_decimal(right, prec)?;
    Ok(a.cmp_value(&b))
}

/// The correctly rounded `f64` value of a decimal literal.
///
/// # Errors
/// Fails on malformed text or when an exponent leaves the range.
pub fn to_float_text(text: &str) -> Result<f64, String> {
    let prec = bits_for_digits(64).min(MAX_BITS);
    Ok(from_decimal(text, prec)?.to_f64())
}

/// The exact decimal expansion of an `f64`, to `digits` significant digits.
///
/// # Errors
/// Fails for a non-finite value or when an exponent leaves the range.
pub fn from_float_text(
    value: f64,
    digits: usize,
) -> Result<String, String> {
    let digits = digits.clamp(1, MAX_DIGITS);
    let prec = bits_for_digits(digits).saturating_add(64);
    let real = Real::from_f64(value, prec)?;
    to_decimal(&real, digits)
}

/// Rounds decimal text to `digits` significant digits in `mode`.
///
/// # Errors
/// Fails on malformed text or when an exponent leaves the range.
pub fn round_text(
    text: &str,
    digits: usize,
    mode: Round,
) -> Result<String, String> {
    decimal::round_decimal(text, digits, mode)
}

/// The number of significant decimal digits a literal carries.
///
/// # Errors
/// Fails on malformed text.
pub fn digits_of(text: &str) -> Result<usize, String> {
    let trimmed = text.trim();
    let body = trimmed.strip_prefix(['-', '+']).unwrap_or(trimmed);
    let mantissa = body.split(['e', 'E']).next().unwrap_or(body);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    if digits.is_empty() {
        return Err("no digits".to_string());
    }
    let significant = digits.trim_start_matches('0');
    if significant.is_empty() {
        return Ok(1);
    }
    Ok(significant.len())
}

#[cfg(test)]
mod tests {
    use super::Round;

    #[test]
    fn constants_are_the_known_digits() {
        assert_eq!(
            super::const_text("pi", 50).unwrap(),
            "3.1415926535897932384626433832795028841971693993751"
        );
        assert_eq!(super::const_text("e", 20).unwrap(), "2.7182818284590452354");
        assert_eq!(
            super::const_text("ln2", 20).unwrap(),
            "0.69314718055994530942"
        );
        assert_eq!(
            super::const_text("sqrt2", 20).unwrap(),
            "1.4142135623730950488"
        );
        assert_eq!(
            super::const_text("phi", 20).unwrap(),
            "1.6180339887498948482"
        );
    }

    #[test]
    fn elementary_functions_match_their_tables() {
        assert_eq!(
            super::eval_text("sin(1)", 20).unwrap(),
            "0.84147098480789650665"
        );
        assert_eq!(
            super::eval_text("exp(1)", 20).unwrap(),
            "2.7182818284590452354"
        );
        assert_eq!(
            super::eval_text("ln(2)", 20).unwrap(),
            "0.69314718055994530942"
        );
        assert_eq!(
            super::eval_text("atan(1)*4", 20).unwrap(),
            "3.1415926535897932385"
        );
        assert_eq!(
            super::eval_text("sqrt(2)", 20).unwrap(),
            "1.4142135623730950488"
        );
        assert_eq!(super::eval_text("cos(0)", 10).unwrap(), "1");
        assert_eq!(super::eval_text("tanh(0)", 10).unwrap(), "0");
    }

    #[test]
    fn division_rounds_once() {
        assert_eq!(
            super::eval_text("1/3", 20).unwrap(),
            "0.33333333333333333333"
        );
        assert_eq!(
            super::eval_text("2/3", 20).unwrap(),
            "0.66666666666666666667"
        );
        assert_eq!(
            super::arithmetic_text("div", "1", "7", 40, Round::NearestEven).unwrap(),
            "0.1428571428571428571428571428571428571429"
        );
    }

    #[test]
    fn rounding_is_nearest_even_and_directed() {
        assert_eq!(
            super::round_text("2.5", 1, Round::NearestEven).unwrap(),
            "2"
        );
        assert_eq!(
            super::round_text("3.5", 1, Round::NearestEven).unwrap(),
            "4"
        );
        assert_eq!(
            super::round_text("1.25", 2, Round::NearestEven).unwrap(),
            "1.2"
        );
        assert_eq!(
            super::round_text("1.25", 2, Round::TowardPositive).unwrap(),
            "1.3"
        );
        assert_eq!(
            super::round_text("-1.25", 2, Round::TowardNegative).unwrap(),
            "-1.3"
        );
        assert_eq!(
            super::round_text("9.99", 2, Round::NearestEven).unwrap(),
            "10"
        );
    }

    #[test]
    fn decimal_rounding_reads_past_a_binary_halfway_approximation() {
        let above_half = format!("1.25{}1", "0".repeat(120));
        let below_half = format!("1.24{}9", "9".repeat(120));
        assert_eq!(
            super::round_text(&above_half, 2, Round::NearestEven).unwrap(),
            "1.3"
        );
        assert_eq!(
            super::round_text(&below_half, 2, Round::NearestEven).unwrap(),
            "1.2"
        );
        assert_eq!(
            super::round_text(&format!("-{above_half}"), 2, Round::NearestEven).unwrap(),
            "-1.3"
        );

        for (mode, positive, negative) in [
            (Round::NearestEven, "1.2", "-1.2"),
            (Round::TowardZero, "1.2", "-1.2"),
            (Round::AwayFromZero, "1.3", "-1.3"),
            (Round::TowardPositive, "1.3", "-1.2"),
            (Round::TowardNegative, "1.2", "-1.3"),
        ] {
            assert_eq!(super::round_text("1.2500", 2, mode).unwrap(), positive);
            assert_eq!(super::round_text("-1.2500", 2, mode).unwrap(), negative);
        }
        assert_eq!(
            super::round_text("1.35", 2, Round::NearestEven).unwrap(),
            "1.4"
        );
        assert_eq!(
            super::round_text("-1.35", 2, Round::NearestEven).unwrap(),
            "-1.4"
        );
        assert_eq!(
            super::round_text("9.95", 2, Round::NearestEven).unwrap(),
            "10"
        );
        assert_eq!(
            super::round_text("1.2995", 4, Round::NearestEven).unwrap(),
            "1.3"
        );
        assert_eq!(
            super::round_text("000.000125e2", 2, Round::TowardPositive).unwrap(),
            "0.013"
        );
        assert_eq!(
            super::round_text("-0", 1, Round::AwayFromZero).unwrap(),
            "0"
        );
    }

    #[test]
    fn exact_literal_evaluation_obeys_the_expression_grammar() {
        assert_eq!(super::eval_text(" + 0.15 ", 1).unwrap(), "0.2");
        assert_eq!(super::eval_text("-0.15", 1).unwrap(), "-0.2");
        assert_eq!(super::eval_text("0.15e-2", 1).unwrap(), "0.002");
        assert_eq!(super::eval_text("1 + 2", 1).unwrap(), "3");
        for malformed in ["1 2", "1_2", "1.2.3", "0.15e", "1e40001"] {
            assert!(super::eval_text(malformed, 1).is_err(), "{malformed}");
        }
        assert!(super::eval_text(&"1".repeat(4097), 1).is_err());
    }

    #[test]
    fn narrowing_to_a_double_is_correctly_rounded() {
        assert_eq!(super::to_float_text("0.1").unwrap(), 0.1);
        assert_eq!(super::to_float_text("1e-400").unwrap(), 0.0);
        assert!(super::to_float_text("1e400").unwrap().is_infinite());
        // The exact decimal expansion of the double nearest 0.1.
        assert_eq!(
            super::from_float_text(0.1, 30).unwrap(),
            "0.100000000000000005551115123126"
        );
        assert_eq!(super::from_float_text(3.0, 5).unwrap(), "3");
    }

    #[test]
    fn comparison_is_exact_beyond_double_precision() {
        assert_eq!(
            super::compare_text("1.0000000000000000000000000001", "1").unwrap(),
            1
        );
        assert_eq!(
            super::compare_text("0.1", "0.10000000000000001").unwrap(),
            -1
        );
        assert_eq!(super::compare_text("-2", "-2.000").unwrap(), 0);
    }

    #[test]
    fn enclosure_contains_the_value() {
        let (low, high) = super::enclosure_text("pi", 25).unwrap();
        assert!(super::compare_text(&low, "3.141592653589793238462643").unwrap() <= 0);
        assert!(super::compare_text(&high, "3.141592653589793238462643").unwrap() >= 0);
        assert!(super::compare_text(&low, &high).unwrap() <= 0);
    }

    #[test]
    fn decimal_enclosures_round_each_endpoint_outward() {
        for (expression, lower, upper) in [
            ("1/3", "0.3", "0.4"),
            ("-1/3", "-0.4", "-0.3"),
            ("9.99", "9", "10"),
            ("-9.99", "-10", "-9"),
        ] {
            let bounds = super::enclosure_text(expression, 1).unwrap();
            assert_eq!(bounds, (lower.to_owned(), upper.to_owned()), "{expression}");
        }
        assert_eq!(
            super::enclosure_complex_text("(1-I)/3", 1).unwrap(),
            (
                ("0.3".to_owned(), "0.4".to_owned()),
                ("-0.4".to_owned(), "-0.3".to_owned())
            )
        );
    }

    #[test]
    fn digit_count_and_errors() {
        assert_eq!(super::digits_of("0.0012300").unwrap(), 5);
        assert!(super::eval_text("1/0", 10).is_err());
        assert!(super::eval_text("ln(0)", 10).is_err());
        assert!(super::eval_text("sqrt(-1)", 10).is_err());
        assert!(super::eval_text("nonesuch(1)", 10).is_err());
        assert!(super::eval_text("1 +", 10).is_err());
        assert!(super::const_text("nonesuch", 10).is_err());
    }

    #[test]
    fn long_precision_is_still_correct() {
        // 200 digits of pi: the last digits are a fingerprint of the series.
        let text = super::const_text("pi", 200).unwrap();
        assert!(text.starts_with(
            "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
        ));
        assert!(text.len() > 190);
    }
}
