// Certified inverse tangent and principal argument.

binding!(
    mp_real_atan, "num_mp_real_atan", "numerical_special",
    "`num_mp_real_atan(value, precision)` encloses the real inverse tangent of a dyadic.",
    "crate::kernels::mp::funcs::atan", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::atan(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_atan, "num_mp_ball_atan", "numerical_special",
    "`num_mp_ball_atan(interval, precision)` encloses the real inverse tangent of a closed interval.",
    "crate::kernels::mp::funcs::atan_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::atan_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_atan2, "num_mp_real_atan2", "numerical_special",
    "`num_mp_real_atan2(y, x, precision)` encloses the principal angle of the exact dyadic coordinates (x, y).",
    "crate::kernels::mp::funcs::atan2", DimensionPolicy::Dimensionless,
    (y: super::mp::Real, x: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |y: super::mp::Real, x: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::atan2(&y, &x, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_atan2, "num_mp_ball_atan2", "numerical_special",
    "`num_mp_ball_atan2(y, x, precision)` encloses the principal angle of a closed rectangle away from the origin and wrapping cut.",
    "crate::kernels::mp::funcs::atan2_ball", DimensionPolicy::Dimensionless,
    (y: super::mp::Ball, x: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |y: super::mp::Ball, x: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::atan2_ball(&y, &x, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_arg, "num_mp_complex_arg", "numerical_special",
    "`num_mp_complex_arg(value, precision)` encloses the principal argument of an exact dyadic rectangle off the non-positive real axis.",
    "crate::kernels::mp::complex::MpComplex::arg", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.arg(precision.0).map_err(OpError::Domain) }
);
