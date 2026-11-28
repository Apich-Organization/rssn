// Certified complex gamma family on dimensionless dyadic rectangles.

binding!(
    mp_complex_log_gamma, "num_mp_complex_log_gamma", "numerical_special",
    "Encloses the analytic log-gamma branch real on positive reals. The closed dyadic rectangle must avoid the nonpositive real cut. Precision is at most 1024 bits and recurrence work is bounded.",
    "crate::kernels::mp::special_complex::MpComplex::ln_gamma", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.ln_gamma(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_gamma, "num_mp_complex_gamma", "numerical_special",
    "Encloses gamma for a closed dyadic rectangle separated from the nonpositive integer poles. Precision is at most 1024 bits and recurrence work is bounded.",
    "crate::kernels::mp::special_complex::MpComplex::gamma", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.gamma(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_digamma, "num_mp_complex_digamma", "numerical_special",
    "Encloses the logarithmic derivative of gamma for a closed dyadic rectangle separated from the poles. Precision is at most 1024 bits and recurrence work is bounded.",
    "crate::kernels::mp::special_complex::MpComplex::digamma", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.digamma(precision.0).map_err(OpError::Domain) }
);
