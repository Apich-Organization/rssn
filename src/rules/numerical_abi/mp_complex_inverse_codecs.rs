// Certified principal inverse functions of complex dyadic rectangles.

binding!(
    mp_complex_asin, "num_mp_complex_asin", "numerical_special",
    "`num_mp_complex_asin(value, precision)` encloses the principal inverse sine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::asin", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.asin(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_acos, "num_mp_complex_acos", "numerical_special",
    "`num_mp_complex_acos(value, precision)` encloses the principal inverse cosine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::acos", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.acos(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_atan, "num_mp_complex_atan", "numerical_special",
    "`num_mp_complex_atan(value, precision)` encloses the principal inverse tangent of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::atan", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.atan(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_asinh, "num_mp_complex_asinh", "numerical_special",
    "`num_mp_complex_asinh(value, precision)` encloses the principal inverse hyperbolic sine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::asinh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.asinh(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_acosh, "num_mp_complex_acosh", "numerical_special",
    "`num_mp_complex_acosh(value, precision)` encloses the principal inverse hyperbolic cosine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::acosh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.acosh(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_atanh, "num_mp_complex_atanh", "numerical_special",
    "`num_mp_complex_atanh(value, precision)` encloses the principal inverse hyperbolic tangent of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::atanh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.atanh(precision.0).map_err(OpError::Domain) }
);
