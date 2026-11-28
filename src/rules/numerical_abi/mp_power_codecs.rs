// Certified real dyadic, integer and rational powers.

binding!(
    mp_real_pow, "num_mp_real_pow", "numerical_special",
    "`num_mp_real_pow(base, exponent, precision)` encloses the real power of two exact dyadics.",
    "crate::kernels::mp::funcs::pow", DimensionPolicy::Dimensionless,
    (base: super::mp::Real, exponent: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Real, exponent: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::pow(&base, &exponent, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_pow, "num_mp_ball_pow", "numerical_special",
    "`num_mp_ball_pow(base, exponent, precision)` encloses a real power over closed base and exponent intervals.",
    "crate::kernels::mp::funcs::pow_ball", DimensionPolicy::Dimensionless,
    (base: super::mp::Ball, exponent: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Ball, exponent: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::pow_ball(&base, &exponent, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_powi, "num_mp_real_powi", "numerical_special",
    "`num_mp_real_powi(base, exponent, precision)` encloses an exact dyadic raised to a signed integer exponent.",
    "crate::kernels::mp::funcs::powi", DimensionPolicy::Dimensionless,
    (base: super::mp::Real, exponent: i64, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Real, exponent: i64, precision: super::MpRealPrecision| { super::mp::funcs::powi(&base, exponent, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_powi, "num_mp_ball_powi", "numerical_special",
    "`num_mp_ball_powi(base, exponent, precision)` encloses a closed interval raised to a signed integer exponent.",
    "crate::kernels::mp::funcs::powi_ball", DimensionPolicy::Dimensionless,
    (base: super::mp::Ball, exponent: i64, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Ball, exponent: i64, precision: super::MpRealPrecision| { super::mp::funcs::powi_ball(&base, exponent, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_pow_ratio, "num_mp_real_pow_ratio", "numerical_special",
    "`num_mp_real_pow_ratio(base, numerator, denominator, precision)` encloses an exact dyadic raised to the rational exponent numerator / denominator.",
    "crate::kernels::mp::funcs::pow_ratio", DimensionPolicy::Dimensionless,
    (base: super::mp::Real, numerator: i64, denominator: u32, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Real, numerator: i64, denominator: u32, precision: super::MpRealPrecision| { super::mp::funcs::pow_ratio(&base, numerator, denominator, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_pow_ratio, "num_mp_ball_pow_ratio", "numerical_special",
    "`num_mp_ball_pow_ratio(base, numerator, denominator, precision)` encloses a closed interval raised to an exact rational exponent.",
    "crate::kernels::mp::funcs::pow_ratio_ball", DimensionPolicy::Dimensionless,
    (base: super::mp::Ball, numerator: i64, denominator: u32, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |base: super::mp::Ball, numerator: i64, denominator: u32, precision: super::MpRealPrecision| { super::mp::funcs::pow_ratio_ball(&base, numerator, denominator, precision.0).map_err(OpError::Domain) }
);
