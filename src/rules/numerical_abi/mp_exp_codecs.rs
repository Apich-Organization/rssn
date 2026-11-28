// The real and interval inputs share one certified exponential.
binding!(
    mp_real_exp, "num_mp_real_exp", "numerical_special",
    "`num_mp_real_exp(value, precision)` encloses the exponential of an exact dyadic.",
    "crate::kernels::mp::funcs::exp", DimensionPolicy::Dimensionless,
    (value: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::Real, precision: super::MpRealPrecision| {
        super::mp::funcs::exp(&value, precision.0).map_err(OpError::Domain)
    }
);

binding!(
    mp_ball_exp, "num_mp_ball_exp", "numerical_special",
    "`num_mp_ball_exp(interval, precision)` encloses the exponential over a closed dyadic interval.",
    "crate::kernels::mp::funcs::exp_ball", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, precision: super::MpRealPrecision| {
        super::mp::funcs::exp_ball(&interval, precision.0).map_err(OpError::Domain)
    }
);
