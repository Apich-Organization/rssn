// Certified principal square root and powers for complex dyadic rectangles.

binding!(
    mp_complex_sqrt, "num_mp_complex_sqrt", "numerical_special",
    "`num_mp_complex_sqrt(value, precision)` encloses the principal square root of a closed dyadic rectangle.",
    "crate::kernels::mp::complex::MpComplex::sqrt", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.sqrt(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_pow, "num_mp_complex_pow", "numerical_special",
    "`num_mp_complex_pow(base, exponent, precision)` encloses a signed integer power or the principal complex power of dyadic rectangles.",
    "crate::kernels::mp::complex::MpComplex::pow", DimensionPolicy::Dimensionless,
    (base: super::mp::MpComplex, exponent: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |base: super::mp::MpComplex, exponent: super::mp::MpComplex, precision: super::MpRealPrecision| { base.pow(&exponent, precision.0).map_err(OpError::Domain) }
);
