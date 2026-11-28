// Directed inverse sine and cosine bounds.

binding!(
    mp_real_asin, "num_mp_real_asin", "numerical_special",
    "`num_mp_real_asin(value, precision)` encloses the inverse sine on the real domain [-1, 1].",
    "crate::kernels::mp::funcs::asin", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::asin(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_asin, "num_mp_ball_asin", "numerical_special",
    "`num_mp_ball_asin(interval, precision)` encloses the inverse sine on the real domain [-1, 1].",
    "crate::kernels::mp::funcs::asin_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::asin_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_acos, "num_mp_real_acos", "numerical_special",
    "`num_mp_real_acos(value, precision)` encloses the inverse cosine on the real domain [-1, 1].",
    "crate::kernels::mp::funcs::acos", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::acos(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_acos, "num_mp_ball_acos", "numerical_special",
    "`num_mp_ball_acos(interval, precision)` encloses the inverse cosine on the real domain [-1, 1].",
    "crate::kernels::mp::funcs::acos_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::acos_ball(&interval, precision.0).map_err(OpError::Domain) }
);
