// Certified point and interval hyperbolic functions.

binding!(
    mp_real_sinh, "num_mp_real_sinh", "numerical_special",
    "`num_mp_real_sinh(value, precision)` encloses the hyperbolic sine of an exact dyadic.",
    "crate::kernels::mp::funcs::sinh", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::sinh(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_sinh, "num_mp_ball_sinh", "numerical_special",
    "`num_mp_ball_sinh(interval, precision)` encloses the hyperbolic sine over a closed dyadic interval.",
    "crate::kernels::mp::funcs::sinh_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::sinh_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_cosh, "num_mp_real_cosh", "numerical_special",
    "`num_mp_real_cosh(value, precision)` encloses the hyperbolic cosine of an exact dyadic.",
    "crate::kernels::mp::funcs::cosh", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::cosh(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_cosh, "num_mp_ball_cosh", "numerical_special",
    "`num_mp_ball_cosh(interval, precision)` encloses the hyperbolic cosine over a closed dyadic interval.",
    "crate::kernels::mp::funcs::cosh_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::cosh_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_tanh, "num_mp_real_tanh", "numerical_special",
    "`num_mp_real_tanh(value, precision)` encloses the hyperbolic tangent of an exact dyadic.",
    "crate::kernels::mp::funcs::tanh", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::tanh(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_tanh, "num_mp_ball_tanh", "numerical_special",
    "`num_mp_ball_tanh(interval, precision)` encloses the hyperbolic tangent over a closed dyadic interval.",
    "crate::kernels::mp::funcs::tanh_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::tanh_ball(&interval, precision.0).map_err(OpError::Domain) }
);
