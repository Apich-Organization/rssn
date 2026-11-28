// Certified trigonometric enclosures reuse dyadic interval codecs.
record!(mp::funcs::SinCos {
    sin: mp::Ball,
    cos: mp::Ball,
});

binding!(
    mp_real_sin, "num_mp_real_sin", "numerical_special",
    "`num_mp_real_sin(value, precision)` encloses the real sin of an exact dyadic.",
    "crate::kernels::mp::funcs::sin", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::sin(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_sin, "num_mp_ball_sin", "numerical_special",
    "`num_mp_ball_sin(interval, precision)` encloses the real sin over a closed dyadic interval.",
    "crate::kernels::mp::funcs::sin_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::sin_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_cos, "num_mp_real_cos", "numerical_special",
    "`num_mp_real_cos(value, precision)` encloses the real cos of an exact dyadic.",
    "crate::kernels::mp::funcs::cos", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::cos(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_cos, "num_mp_ball_cos", "numerical_special",
    "`num_mp_ball_cos(interval, precision)` encloses the real cos over a closed dyadic interval.",
    "crate::kernels::mp::funcs::cos_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::cos_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_tan, "num_mp_real_tan", "numerical_special",
    "`num_mp_real_tan(value, precision)` encloses the real tan of an exact dyadic.",
    "crate::kernels::mp::funcs::tan", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::tan(&value, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_ball_tan, "num_mp_ball_tan", "numerical_special",
    "`num_mp_ball_tan(interval, precision)` encloses the real tan over a closed dyadic interval.",
    "crate::kernels::mp::funcs::tan_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| { super::mp::funcs::tan_ball(&interval, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_sin_cos, "num_mp_real_sin_cos", "numerical_special",
    "`num_mp_real_sin_cos(value, precision)` encloses sine and cosine with one certified argument reduction.",
    "crate::kernels::mp::funcs::sin_cos", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::funcs::SinCos, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| { super::mp::funcs::sin_cos(&value, precision.0).map_err(OpError::Domain) }
);
