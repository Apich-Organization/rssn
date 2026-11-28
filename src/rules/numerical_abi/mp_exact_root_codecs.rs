// Checked dyadic accessors, format limits, and outward algebraic roots.
record!(mp::real::RealLimits {
    min_precision: u32,
    max_precision: u32,
    min_exponent: i64,
    max_exponent: i64,
});

binding!(
    mp_real_mantissa, "num_mp_real_mantissa", "numerical_special",
    "`num_mp_real_mantissa(value)` Returns the exact arbitrary-width integer mantissa.",
    "crate::kernels::mp::real::Real::mantissa", DimensionPolicy::Dimensionless,
    (value: super::mp::Real) -> num_bigint::BigInt,
    |value: super::mp::Real| { value.mantissa().clone() }
);

binding!(
    mp_real_exponent, "num_mp_real_exponent", "numerical_special",
    "`num_mp_real_exponent(value)` Returns the binary exponent of the canonical dyadic.",
    "crate::kernels::mp::real::Real::exponent", DimensionPolicy::Dimensionless,
    (value: super::mp::Real) -> i64,
    |value: super::mp::Real| { value.exponent() }
);

binding!(
    mp_real_precision, "num_mp_real_precision", "numerical_special",
    "`num_mp_real_precision(value)` Returns the dyadic significand precision in bits.",
    "crate::kernels::mp::real::Real::precision", DimensionPolicy::Dimensionless,
    (value: super::mp::Real) -> u32,
    |value: super::mp::Real| { value.precision() }
);

binding!(
    mp_ball_midpoint, "num_mp_ball_midpoint", "numerical_special",
    "`num_mp_ball_midpoint(interval)` Returns the midpoint of the closed dyadic interval.",
    "crate::kernels::mp::ball::Ball::midpoint", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball) -> super::mp::Real,
    |interval: super::mp::Ball| { interval.midpoint().clone() }
);

binding!(
    mp_ball_radius, "num_mp_ball_radius", "numerical_special",
    "`num_mp_ball_radius(interval)` Returns the nonnegative radius of the closed dyadic interval.",
    "crate::kernels::mp::ball::Ball::radius", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball) -> super::mp::Real,
    |interval: super::mp::Ball| { interval.radius().clone() }
);

binding!(
    mp_real_limits, "num_mp_real_limits", "numerical_special",
    "`num_mp_real_limits()` Reports the supported dyadic precision and exponent limits.",
    "crate::kernels::mp::real::limits", DimensionPolicy::Dimensionless,
    () -> super::mp::real::RealLimits,
    || { super::mp::real::limits() }
);

binding!(
    mp_real_rounding_mode, "num_mp_real_rounding_mode", "numerical_special",
    "`num_mp_real_rounding_mode(rounding)` Returns the canonical rounding name for a supported spelling.",
    "crate::kernels::mp::round_named", DimensionPolicy::Dimensionless,
    (rounding: super::mp::Round) -> String,
    |rounding: super::mp::Round| { match rounding { super::mp::Round::NearestEven => "nearest", super::mp::Round::TowardZero => "zero", super::mp::Round::AwayFromZero => "away", super::mp::Round::TowardPositive => "up", super::mp::Round::TowardNegative => "down" }.to_owned() }
);

binding!(
    mp_real_cbrt, "num_mp_real_cbrt", "numerical_special",
    "`num_mp_real_cbrt(value, precision)` Encloses the real cube root of a dyadic, including negative inputs.",
    "crate::kernels::mp::funcs::cbrt", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::cbrt(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_cbrt_in, "num_mp_real_cbrt_in", "numerical_special",
    "`num_mp_real_cbrt_in(value, precision, rounding)` Rounds the real cube root of a dyadic in the requested direction.",
    "crate::kernels::mp::real::Real::cbrt", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision, rounding: super::mp::Round) -> Result<super::mp::Real, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision, rounding: super::mp::Round| { value.cbrt(precision.0, rounding).map_err(OpError::Domain) }
);

binding!(
    mp_real_hypot, "num_mp_real_hypot", "numerical_special",
    "`num_mp_real_hypot(left, right, precision)` Encloses the Euclidean norm of two exact dyadics.",
    "crate::kernels::mp::funcs::hypot", DimensionPolicy::Dimensionless,
    (left: super::mp::Real, right: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |left: super::mp::Real, right: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::hypot(&left, &right, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_sqrt, "num_mp_ball_sqrt", "numerical_special",
    "`num_mp_ball_sqrt(interval, precision)` Encloses the square root of a nonnegative closed interval.",
    "crate::kernels::mp::ball::Ball::sqrt", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { interval.sqrt(precision.0).map_err(OpError::Domain) }
);
