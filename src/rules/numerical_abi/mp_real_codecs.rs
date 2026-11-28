// Portable, exact dyadic records for the finite multiprecision real kernels.
use crate::kernels::mp;

const MP_REAL_FIELDS: [&str; 3] = ["mantissa", "exponent", "precision"];

fn mp_real_mantissa(
    value: &Value,
    index: usize,
) -> Result<&num_bigint::BigInt, OpError> {
    let Value::Int(mantissa) = value else {
        return Err(OpError::Type {
            expected: Box::new(TypeExpr::Int),
            got: Box::new(value.type_expr()),
            arg: index,
        });
    };
    Ok(mantissa)
}

struct MpRealInputMantissa(num_bigint::BigInt);

impl Arg for MpRealInputMantissa {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let mantissa = mp_real_mantissa(value, index)?;
        let limit = mp::MAX_BITS.saturating_mul(2);
        if mantissa.bits() > u64::from(limit) {
            return Err(OpError::Domain(format!(
                "dyadic constructor mantissa must fit {limit} bits"
            )));
        }
        Ok(Self(mantissa.clone()))
    }
}

struct MpRealExponent(i64);

impl Arg for MpRealExponent {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let exponent = i64::decode(value, index)?;
        if exponent.unsigned_abs() > mp::MAX_EXP.unsigned_abs() {
            return Err(OpError::Domain(format!(
                "dyadic exponent must lie within +/-{}",
                mp::MAX_EXP
            )));
        }
        Ok(Self(exponent))
    }
}

struct MpRealPrecision(u32);

impl Arg for MpRealPrecision {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let precision = u32::decode(value, index)?;
        if !(1..=mp::MAX_BITS).contains(&precision) {
            return Err(OpError::Domain(format!(
                "dyadic precision must be in 1..={}",
                mp::MAX_BITS
            )));
        }
        Ok(Self(precision))
    }
}

impl Arg for mp::Round {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::String(mode) = value else {
            return Err(OpError::Type {
                expected: Box::new(TypeExpr::String),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if mode.len() > 32 {
            return Err(OpError::Domain(
                "dyadic rounding mode name is too long".into(),
            ));
        }
        mp::round_named(mode).map_err(OpError::Domain)
    }
}

impl Arg for mp::Real {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &MP_REAL_FIELDS)?;
        let mantissa = mp_real_mantissa(native_field(record, "mantissa")?, index)?;
        let MpRealExponent(exponent) =
            MpRealExponent::decode(native_field(record, "exponent")?, index)?;
        let MpRealPrecision(precision) =
            MpRealPrecision::decode(native_field(record, "precision")?, index)?;
        let bits = mantissa.bits();
        if (bits == 0 && exponent != 0) || (bits != 0 && bits != u64::from(precision)) {
            return Err(OpError::Domain(
                "dyadic records require a mantissa of the stated precision or zero with exponent zero"
                    .into(),
            ));
        }
        // The preceding validation makes this an exact reconstruction in any rounding mode.
        Self::rounded(
            mantissa.clone(),
            exponent,
            precision,
            mp::Round::NearestEven,
        )
        .map_err(OpError::Domain)
    }
}

impl Output for mp::Real {
    fn ty() -> TypeExpr {
        TypeExpr::Record(
            MP_REAL_FIELDS
                .into_iter()
                .map(|name| (name.to_owned(), TypeExpr::Int))
                .collect(),
        )
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Record(RecordValue {
            fields: MP_REAL_FIELDS
                .into_iter()
                .map(str::to_owned)
                .zip([
                    Value::int(self.mantissa().clone()),
                    Value::int(self.exponent()),
                    Value::int(self.precision()),
                ])
                .collect(),
        }))
    }
}

binding!(
    mp_real_round, "num_mp_real_round", "numerical_special",
    "Rounds mantissa times 2^exponent once to the requested bit precision. Returns the integer fields mantissa, exponent and precision. Nonzero mantissas have exactly precision bits and zero has exponent zero. Precision is 1 to 32768 bits, input mantissas may have up to 65536 bits and exponents must lie within +/-262144. Accepts nearest (ties to even), zero, away, up, or down.",
    "crate::kernels::mp::real::Real::rounded", DimensionPolicy::Dimensionless,
    (mantissa: super::MpRealInputMantissa, exponent: super::MpRealExponent,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round)
        -> Result<crate::kernels::mp::Real, OpError>,
    |mantissa: super::MpRealInputMantissa, exponent: super::MpRealExponent,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round| {
        crate::kernels::mp::Real::rounded(mantissa.0, exponent.0, precision.0, mode)
            .map_err(OpError::Domain)
    }
);
binding!(
    mp_real_compare, "num_mp_real_compare", "numerical_special",
    "Compares canonical dyadic records exactly. Returns -1, 0, or 1 and allows different operand precisions.",
    "crate::kernels::mp::real::Real::cmp_value", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Real, right: crate::kernels::mp::Real) -> i32,
    |left: crate::kernels::mp::Real, right: crate::kernels::mp::Real| i32::from(left.cmp_value(&right))
);
binding!(
    mp_real_add, "num_mp_real_add", "numerical_special",
    "Adds two canonical dyadic records and rounds the exact sum once to the requested bit precision and mode.",
    "crate::kernels::mp::real::Real::add", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round)
        -> Result<crate::kernels::mp::Real, OpError>,
    |left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round| {
        left.add(&right, precision.0, mode).map_err(OpError::Domain)
    }
);
binding!(
    mp_real_sub, "num_mp_real_sub", "numerical_special",
    "Subtracts two canonical dyadic records and rounds the exact difference once to the requested bit precision and mode.",
    "crate::kernels::mp::real::Real::sub", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round)
        -> Result<crate::kernels::mp::Real, OpError>,
    |left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round| {
        left.sub(&right, precision.0, mode).map_err(OpError::Domain)
    }
);
binding!(
    mp_real_mul, "num_mp_real_mul", "numerical_special",
    "Multiplies two canonical dyadic records and rounds the exact product once to the requested bit precision and mode.",
    "crate::kernels::mp::real::Real::mul", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round)
        -> Result<crate::kernels::mp::Real, OpError>,
    |left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round| {
        left.mul(&right, precision.0, mode).map_err(OpError::Domain)
    }
);
binding!(
    mp_real_div, "num_mp_real_div", "numerical_special",
    "Divides two canonical dyadic records and rounds the exact quotient once to the requested bit precision and mode. A zero divisor or an out-of-range normalized exponent is a domain error.",
    "crate::kernels::mp::real::Real::div", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round)
        -> Result<crate::kernels::mp::Real, OpError>,
    |left: crate::kernels::mp::Real, right: crate::kernels::mp::Real,
     precision: super::MpRealPrecision, mode: crate::kernels::mp::Round| {
        left.div(&right, precision.0, mode).map_err(OpError::Domain)
    }
);


binding!(
    mp_real_abs, "num_mp_real_abs", "numerical_special",
    "Returns the absolute value of an exact dyadic without rounding.",
    "crate::kernels::mp::real::Real::abs", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> crate::kernels::mp::Real,
    |value: crate::kernels::mp::Real| {
        value.abs()
    }
);


binding!(
    mp_real_neg, "num_mp_real_neg", "numerical_special",
    "Negates an exact dyadic at its existing precision.",
    "crate::kernels::mp::real::Real::neg", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> crate::kernels::mp::Real,
    |value: crate::kernels::mp::Real| {
        value.neg()
    }
);


binding!(
    mp_real_is_zero, "num_mp_real_is_zero", "numerical_special",
    "Tests whether a canonical dyadic is zero.",
    "crate::kernels::mp::real::Real::is_zero", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> bool,
    |value: crate::kernels::mp::Real| {
        value.is_zero()
    }
);


binding!(
    mp_real_is_integer, "num_mp_real_is_integer", "numerical_special",
    "Tests whether a dyadic has no fractional binary digits.",
    "crate::kernels::mp::real::Real::is_integer", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> bool,
    |value: crate::kernels::mp::Real| {
        value.is_integer()
    }
);


binding!(
    mp_real_floor_log2, "num_mp_real_floor_log2", "numerical_special",
    "Returns floor(log2(abs(value))) or None if value is zero.",
    "crate::kernels::mp::real::Real::floor_log2", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> Option<i64>,
    |value: crate::kernels::mp::Real| {
        value.floor_log2()
    }
);


binding!(
    mp_real_trunc_to_bigint, "num_mp_real_trunc_to_bigint", "numerical_special",
    "Truncates a dyadic toward zero and returns an exact integer. The dyadic representation bounds the result to 294912 bits.",
    "crate::kernels::mp::real::Real::trunc_to_bigint", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> Result<num_bigint::BigInt, OpError>,
    |value: crate::kernels::mp::Real| {
        value.trunc_to_bigint().map_err(OpError::Domain)
    }
);


binding!(
    mp_real_with_precision, "num_mp_real_with_precision", "numerical_special",
    "Rounds an existing canonical dyadic to a new bit precision. The rounding modes are nearest, zero, away, up and down.",
    "crate::kernels::mp::real::Real::with_precision", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real, precision: super::MpRealPrecision, rounding: crate::kernels::mp::Round) -> Result<crate::kernels::mp::Real, OpError>,
    |value: crate::kernels::mp::Real, precision: super::MpRealPrecision, rounding: crate::kernels::mp::Round| {
        value.with_precision(precision.0, rounding).map_err(OpError::Domain)
    }
);


binding!(
    mp_real_sqrt, "num_mp_real_sqrt", "numerical_special",
    "Rounds the square root of a nonnegative dyadic at the requested bit precision. A negative input is a domain error.",
    "crate::kernels::mp::real::Real::sqrt", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real, precision: super::MpRealPrecision, rounding: crate::kernels::mp::Round) -> Result<crate::kernels::mp::Real, OpError>,
    |value: crate::kernels::mp::Real, precision: super::MpRealPrecision, rounding: crate::kernels::mp::Round| {
        value.sqrt(precision.0, rounding).map_err(OpError::Domain)
    }
);


binding!(
    mp_real_scale_pow2, "num_mp_real_scale_pow2", "numerical_special",
    "Multiplies a dyadic by an integer power of two exactly. Zero keeps its canonical exponent. A requested exponent outside the representable range is a domain error.",
    "crate::kernels::mp::real::Real::scale_pow2", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real, shift: i64) -> Result<crate::kernels::mp::Real, OpError>,
    |value: crate::kernels::mp::Real, shift: i64| {
        value.scale_pow2(shift).map_err(OpError::Domain)
    }
);
