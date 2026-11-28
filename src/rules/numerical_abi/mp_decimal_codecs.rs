// The decimal digit count is distinct from the binary precision of a dyadic.
struct MpDecimalDigits(usize);

impl Arg for MpDecimalDigits {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let digits = usize::try_from(u32::decode(value, index)?)
            .map_err(|_| OpError::Domain("decimal digit count is too large".into()))?;
        if !(1..=mp::MAX_DIGITS).contains(&digits) {
            return Err(OpError::Domain(format!(
                "decimal digits must be in 1..={}",
                mp::MAX_DIGITS
            )));
        }
        Ok(Self(digits))
    }
}

struct MpDecimalText(String);

impl Arg for MpDecimalText {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        if let Value::String(text) = value {
            if text.len() > mp::decimal::MAX_DECIMAL_TEXT {
                return Err(OpError::Domain("decimal text is too long".into()));
            }
            return Ok(Self(text.clone()));
        }
        Err(OpError::Type {
            expected: Box::new(TypeExpr::String),
            got: Box::new(value.type_expr()),
            arg: index,
        })
    }
}

binding!(
    mp_decimal_from_text, "num_mp_decimal_from_text", "numerical_special",
    "`num_mp_decimal_from_text(text, precision)` Parses exact decimal text and rounds to a dyadic with ties to even.",
    "crate::kernels::mp::decimal::from_decimal", DimensionPolicy::Dimensionless,
    (text: super::MpDecimalText, precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |text: super::MpDecimalText, precision: super::MpRealPrecision| { super::mp::decimal::from_decimal(&text.0, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_from_text_in, "num_mp_decimal_from_text_in", "numerical_special",
    "`num_mp_decimal_from_text_in(text, precision, rounding)` Parses exact decimal text and rounds once in the requested direction.",
    "crate::kernels::mp::decimal::from_decimal_in", DimensionPolicy::Dimensionless,
    (text: super::MpDecimalText, precision: super::MpRealPrecision, rounding: super::mp::Round) -> Result<super::mp::Real, OpError>,
    |text: super::MpDecimalText, precision: super::MpRealPrecision, rounding: super::mp::Round| { super::mp::decimal::from_decimal_in(&text.0, precision.0, rounding).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_to_text, "num_mp_decimal_to_text", "numerical_special",
    "`num_mp_decimal_to_text(value, digits)` Formats an exact dyadic to significant digits with ties to even.",
    "crate::kernels::mp::decimal::to_decimal", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, digits: super::MpDecimalDigits) -> Result<String, OpError>,
    |value: super::mp::Real, digits: super::MpDecimalDigits| { super::mp::decimal::to_decimal(&value, digits.0).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_to_text_in, "num_mp_decimal_to_text_in", "numerical_special",
    "`num_mp_decimal_to_text_in(value, digits, rounding)` Formats an exact dyadic to significant digits in the requested direction.",
    "crate::kernels::mp::decimal::to_decimal_in", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, digits: super::MpDecimalDigits, rounding: super::mp::Round) -> Result<String, OpError>,
    |value: super::mp::Real, digits: super::MpDecimalDigits, rounding: super::mp::Round| { super::mp::decimal::to_decimal_in(&value, digits.0, rounding).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_round_text, "num_mp_decimal_round_text", "numerical_special",
    "`num_mp_decimal_round_text(text, digits, rounding)` Rounds a decimal literal directly to significant digits without binary conversion.",
    "crate::kernels::mp::round_text", DimensionPolicy::Dimensionless,
    (text: super::MpDecimalText, digits: super::MpDecimalDigits, rounding: super::mp::Round) -> Result<String, OpError>,
    |text: super::MpDecimalText, digits: super::MpDecimalDigits, rounding: super::mp::Round| { super::mp::round_text(&text.0, digits.0, rounding).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_enclose, "num_mp_decimal_enclose", "numerical_special",
    "`num_mp_decimal_enclose(text, precision)` Encloses an exact decimal in a dyadic ball at the requested midpoint precision.",
    "crate::kernels::mp::decimal::enclose_decimal", DimensionPolicy::Dimensionless,
    (text: super::MpDecimalText, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |text: super::MpDecimalText, precision: super::MpRealPrecision| { super::mp::decimal::enclose_decimal(&text.0, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_format_ball, "num_mp_decimal_format_ball", "numerical_special",
    "`num_mp_decimal_format_ball(interval, digits)` Formats a dyadic ball when both endpoints determine the same rounded digits.",
    "crate::kernels::mp::decimal::ball_to_decimal", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, digits: super::MpDecimalDigits) -> Result<String, OpError>,
    |interval: super::mp::Ball, digits: super::MpDecimalDigits| { super::mp::decimal::ball_to_decimal(&interval, digits.0).map_err(OpError::Domain) }
);

binding!(
    mp_decimal_exponent, "num_mp_decimal_exponent", "numerical_special",
    "`num_mp_decimal_exponent(value)` Returns the exact base-ten order of magnitude, or None for zero.",
    "crate::kernels::mp::decimal::decimal_exponent", DimensionPolicy::Dimensionless,
    (value: super::mp::Real) -> Option<i64>,
    |value: super::mp::Real| { super::mp::decimal::decimal_exponent(&value) }
);

binding!(
    mp_decimal_significant_digits, "num_mp_decimal_significant_digits", "numerical_special",
    "`num_mp_decimal_significant_digits(value)` Counts digits in the exact terminating decimal, capped at 10000.",
    "crate::kernels::mp::decimal::significant_digits", DimensionPolicy::Dimensionless,
    (value: super::mp::Real) -> usize,
    |value: super::mp::Real| { super::mp::decimal::significant_digits(&value) }
);

binding!(
    mp_decimal_working_bits, "num_mp_decimal_working_bits", "numerical_special",
    "`num_mp_decimal_working_bits(digits)` Estimates working binary precision including 16 guard bits for decimal digits.",
    "crate::kernels::mp::decimal::bits_for_digits", DimensionPolicy::Dimensionless,
    (digits: super::MpDecimalDigits) -> u32,
    |digits: super::MpDecimalDigits| { super::mp::decimal::bits_for_digits(digits.0) }
);
