// Certified principal Bessel J on argument and order rectangles.

binding!(
    mp_complex_bessel_j, "num_mp_complex_bessel_j", "numerical_special",
    "Encloses principal Bessel J for dimensionless dyadic argument and order rectangles. Exact integer orders through magnitude 512 admit either side of the cut. General orders require Re(order) >= 0 and up to 512 bits. The argument rectangle must stay within |Re(z)| + |Im(z)| <= 256.",
    "crate::kernels::mp::special_complex::MpComplex::bessel_j", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, order: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, order: super::mp::MpComplex, precision: super::MpRealPrecision| { value.bessel_j(&order, precision.0).map_err(OpError::Domain) }
);
