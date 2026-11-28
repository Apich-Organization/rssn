// Certified dyadic complex rectangle modulus, exponential and logarithm.

binding!(
    mp_complex_abs, "num_mp_complex_abs", "numerical_special",
    "`num_mp_complex_abs(value, precision)` encloses the modulus of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::abs", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.abs(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_exp, "num_mp_complex_exp", "numerical_special",
    "`num_mp_complex_exp(value, precision)` encloses the exponential of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::exp", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.exp(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_ln, "num_mp_complex_ln", "numerical_special",
    "`num_mp_complex_ln(value, precision)` encloses the principal logarithm of a closed dyadic rectangle that avoids the branch cut.",
    "crate::kernels::mp::complex::MpComplex::ln", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.ln(precision.0).map_err(OpError::Domain) }
);
