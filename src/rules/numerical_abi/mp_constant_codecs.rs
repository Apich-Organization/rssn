// Certified constants share the same dyadic codecs as the arithmetic operators.

binding!(
    mp_real_const_pi, "num_mp_real_const_pi", "numerical_special",
    "`num_mp_real_const_pi(precision)` returns the circle constant pi, correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::pi", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |precision: super::MpRealPrecision| { super::mp::consts::pi(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_e, "num_mp_real_const_e", "numerical_special",
    "`num_mp_real_const_e(precision)` returns Euler’s number e, correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::e", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |precision: super::MpRealPrecision| { super::mp::consts::e(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_ln2, "num_mp_real_const_ln2", "numerical_special",
    "`num_mp_real_const_ln2(precision)` returns the natural logarithm of two, correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::ln2", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |precision: super::MpRealPrecision| { super::mp::consts::ln2(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_sqrt2, "num_mp_real_const_sqrt2", "numerical_special",
    "`num_mp_real_const_sqrt2(precision)` returns the square root of two, correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::sqrt2", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |precision: super::MpRealPrecision| { super::mp::consts::sqrt2(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_phi, "num_mp_real_const_phi", "numerical_special",
    "`num_mp_real_const_phi(precision)` returns the golden ratio, correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::phi", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |precision: super::MpRealPrecision| { super::mp::consts::phi(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_named, "num_mp_real_const_named", "numerical_special",
    "`num_mp_real_const_named(name, precision)` returns a named constant correctly rounded to the requested bit precision.",
    "crate::kernels::mp::consts::named", DimensionPolicy::Dimensionless,
    (name: String, precision: super::MpRealPrecision) -> Result<super::mp::Real, OpError>,
    |name: String, precision: super::MpRealPrecision| { super::mp::consts::named(&name, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_enclose, "num_mp_real_const_enclose", "numerical_special",
    "`num_mp_real_const_enclose(name, precision)` returns an outward dyadic enclosure of a named constant.",
    "crate::kernels::mp::consts::named_ball", DimensionPolicy::Dimensionless,
    (name: String, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |name: String, precision: super::MpRealPrecision| { super::mp::consts::named_ball(&name, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_names, "num_mp_real_const_names", "numerical_special",
    "`num_mp_real_const_names()` lists the exact names accepted by the dyadic constant functions.",
    "crate::kernels::mp::consts::names", DimensionPolicy::Dimensionless,
    () -> Vec<String>,
    || { super::mp::consts::names().iter().map(|name| (*name).to_owned()).collect() }
);

binding!(
    mp_real_const_pi_scaled, "num_mp_real_const_pi_scaled", "numerical_special",
    "`num_mp_real_const_pi_scaled(bits)` returns the exact integer floor of pi times 2^bits. Accepts 0 to 32768 bits.",
    "crate::kernels::mp::consts::pi_scaled", DimensionPolicy::Dimensionless,
    (bits: u32) -> Result<num_bigint::BigInt, OpError>,
    |bits: u32| { super::mp::consts::pi_scaled(bits).map_err(OpError::Domain) }
);

binding!(
    mp_real_const_functions, "num_mp_real_const_functions", "numerical_special",
    "`num_mp_real_const_functions()` lists the function names recognized by the multiprecision expression parser.",
    "crate::kernels::mp::parse::functions", DimensionPolicy::Dimensionless,
    () -> Vec<String>,
    || { super::mp::parse::functions().iter().map(|name| (*name).to_owned()).collect() }
);
