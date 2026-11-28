// Certified real logarithms, binary exponential and decimal scaling.

binding!(
    mp_real_ln, "num_mp_real_ln", "numerical_special",
    "`num_mp_real_ln(value, precision)` encloses the natural logarithm of an exact dyadic.",
    "crate::kernels::mp::funcs::ln", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::ln(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_ln, "num_mp_ball_ln", "numerical_special",
    "`num_mp_ball_ln(interval, precision)` encloses the natural logarithm over a closed dyadic interval.",
    "crate::kernels::mp::funcs::ln_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::ln_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_log2, "num_mp_real_log2", "numerical_special",
    "`num_mp_real_log2(value, precision)` encloses the base-two logarithm of an exact dyadic.",
    "crate::kernels::mp::funcs::log2", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::log2(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_log2, "num_mp_ball_log2", "numerical_special",
    "`num_mp_ball_log2(interval, precision)` encloses the base-two logarithm over a closed dyadic interval.",
    "crate::kernels::mp::funcs::log2_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::log2_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_log10, "num_mp_real_log10", "numerical_special",
    "`num_mp_real_log10(value, precision)` encloses the base-ten logarithm of an exact dyadic.",
    "crate::kernels::mp::funcs::log10", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::log10(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_log10, "num_mp_ball_log10", "numerical_special",
    "`num_mp_ball_log10(interval, precision)` encloses the base-ten logarithm over a closed dyadic interval.",
    "crate::kernels::mp::funcs::log10_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::log10_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_exp2, "num_mp_real_exp2", "numerical_special",
    "`num_mp_real_exp2(value, precision)` encloses the base-two exponential of an exact dyadic.",
    "crate::kernels::mp::funcs::exp2", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::exp2(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_exp2, "num_mp_ball_exp2", "numerical_special",
    "`num_mp_ball_exp2(interval, precision)` encloses the base-two exponential over a closed dyadic interval.",
    "crate::kernels::mp::funcs::exp2_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::exp2_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_scale_decimal, "num_mp_real_scale_decimal", "numerical_special",
    "`num_mp_real_scale_decimal(value, exponent, precision)` multiplies an exact dyadic by a power of ten and rounds once to nearest-even.",
    "crate::kernels::mp::funcs::scale_decimal", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, exponent: i64, precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |value: super::mp::Real, exponent: i64, precision: super::MpRealPrecision| { super::mp::funcs::scale_decimal(&value, exponent, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_scale_decimal_in, "num_mp_real_scale_decimal_in", "numerical_special",
    "`num_mp_real_scale_decimal_in(value, exponent, precision, rounding)` multiplies an exact dyadic by a power of ten and rounds once in the requested direction.",
    "crate::kernels::mp::funcs::scale_decimal_in", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, exponent: i64, precision: super::MpRealPrecision, rounding: super::mp::Round) -> Result<super::mp::Real, OpError>,
    |value: super::mp::Real, exponent: i64, precision: super::MpRealPrecision, rounding: super::mp::Round| { super::mp::funcs::scale_decimal_in(&value, exponent, precision.0, rounding).map_err(OpError::Domain) }
);
