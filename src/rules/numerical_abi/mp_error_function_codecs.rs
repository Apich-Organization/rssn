// Certified complex error functions on dimensionless dyadic rectangles.

binding!(
    mp_complex_erf, "num_mp_complex_erf", "numerical_special",
    "Encloses the complex error function over a closed dyadic rectangle with a certified series remainder.",
    "crate::kernels::mp::special_complex::MpComplex::erf", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| value.erf(precision.0).map_err(OpError::Domain)
);

binding!(
    mp_complex_erfc, "num_mp_complex_erfc", "numerical_special",
    "Encloses the complex complementary error function over a closed dyadic rectangle. A certified right-sector tail avoids cancellation in the far positive range.",
    "crate::kernels::mp::special_complex::MpComplex::erfc", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| value.erfc(precision.0).map_err(OpError::Domain)
);
