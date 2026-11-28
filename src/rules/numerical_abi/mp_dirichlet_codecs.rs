// Certified zeta and polylogarithm on dimensionless dyadic rectangles.

binding!(
    mp_complex_riemann_zeta, "num_mp_complex_riemann_zeta", "numerical_special",
    "Encloses the analytic continuation of Riemann zeta. General order rectangles must avoid one and satisfy Re(order) >= -32 and |order| <= 256. Exact negative even integer orders are trivial zeros and the supported exact odd orders extend down to -127. Precision is at most 512 bits.",
    "crate::kernels::mp::special_complex::MpComplex::riemann_zeta", DimensionPolicy::Dimensionless,
    (order: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |order: super::mp::MpComplex, precision: super::MpRealPrecision| { order.riemann_zeta(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_polylog, "num_mp_complex_polylog", "numerical_special",
    "Encloses polylogarithm when the argument rectangle lies strictly within the unit disk and the order rectangle has Re(order) >= -256. Returns an error when a certified tail cannot be reached within the bounded series. Precision is at most 512 bits.",
    "crate::kernels::mp::special_complex::MpComplex::polylog", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, order: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, order: super::mp::MpComplex, precision: super::MpRealPrecision| { value.polylog(&order, precision.0).map_err(OpError::Domain) }
);
