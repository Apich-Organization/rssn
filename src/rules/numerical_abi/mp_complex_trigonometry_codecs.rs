// Certified circular and hyperbolic functions of complex dyadic rectangles.

binding!(
    mp_complex_sin, "num_mp_complex_sin", "numerical_special",
    "`num_mp_complex_sin(value, precision)` encloses the sine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::sin", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.sin(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_cos, "num_mp_complex_cos", "numerical_special",
    "`num_mp_complex_cos(value, precision)` encloses the cosine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::cos", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.cos(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_tan, "num_mp_complex_tan", "numerical_special",
    "`num_mp_complex_tan(value, precision)` encloses the tangent of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::tan", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.tan(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_sinh, "num_mp_complex_sinh", "numerical_special",
    "`num_mp_complex_sinh(value, precision)` encloses the hyperbolic sine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::sinh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.sinh(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_cosh, "num_mp_complex_cosh", "numerical_special",
    "`num_mp_complex_cosh(value, precision)` encloses the hyperbolic cosine of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::cosh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.cosh(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_tanh, "num_mp_complex_tanh", "numerical_special",
    "`num_mp_complex_tanh(value, precision)` encloses the hyperbolic tangent of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::tanh", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.tanh(precision.0).map_err(OpError::Domain) }
);
