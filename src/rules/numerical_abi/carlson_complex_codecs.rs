// Typed contracts for the complex Carlson symmetric elliptic integrals.
use crate::kernels::special_complex_ext::elliptic_modular as carlson;

fn checked_carlson(value: Complex64) -> Result<Complex64, OpError> {
    if !value.is_finite() {
        return Err(OpError::Domain(
            "Carlson integral is singular, outside binary64 range or did not converge".into(),
        ));
    }
    Ok(value)
}

fn carlson_rj_branch(
    x: Complex64,
    y: Complex64,
    z: Complex64,
    p: Complex64,
) -> bool {
    // Carlson (1995), section 2.3. The R_D degeneration also accepts
    // other quadrants, provided the arguments select the same branch lip.
    [x, y, z]
        .into_iter()
        .any(|value| carlson::carlson_same_branch(p, value))
        || (x.re >= 0.0 && y.re >= 0.0 && z.re >= 0.0 && p.re > 0.0)
}

binding!(
    carlson_rf_complex,
    "num_special_complex_carlson_rf",
    "numerical_special",
    "Carlson R_F accepts up to one zero and uses principal square roots. Inputs and output are dimensionless. Poles and values outside binary64 range produce a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::carlson_rf_complex",
    DimensionPolicy::Dimensionless,
    (x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64| {
        super::checked_carlson(super::carlson::carlson_rf_complex(x, y, z))
    }
);

binding!(
    carlson_rc_complex,
    "num_special_complex_carlson_rc",
    "numerical_special",
    "Carlson R_C(x, y) is R_F(x, y, y); y must be nonzero. Inputs and output are dimensionless. Poles and values outside binary64 range produce a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::carlson_rc_complex",
    DimensionPolicy::Dimensionless,
    (x: num_complex::Complex64, y: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |x: num_complex::Complex64, y: num_complex::Complex64| {
        super::checked_carlson(super::carlson::carlson_rc_complex(x, y))
    }
);

binding!(
    carlson_rd_complex,
    "num_special_complex_carlson_rd",
    "numerical_special",
    "Carlson R_D uses principal square roots. The last argument must be nonzero, and the first two cannot both vanish. Inputs and output are dimensionless. Poles and values outside binary64 range produce a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::carlson_rd_complex",
    DimensionPolicy::Dimensionless,
    (x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64| {
        super::checked_carlson(super::carlson::carlson_rd_complex(x, y, z))
    }
);

binding!(
    carlson_rg_complex,
    "num_special_complex_carlson_rg",
    "numerical_special",
    "Carlson R_G uses principal square roots and accepts zero arguments. Inputs and output are dimensionless. Poles and values outside binary64 range produce a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::carlson_rg_complex",
    DimensionPolicy::Dimensionless,
    (x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64| {
        super::checked_carlson(super::carlson::carlson_rg_complex(x, y, z))
    }
);

binding!(
    carlson_rj_complex,
    "num_special_complex_carlson_rj",
    "numerical_special",
    "Carlson R_J accepts nonnegative real parts for x, y, z and a positive real part for p. If p equals one of x, y, z, the corresponding R_D branch applies. Inputs and output are dimensionless. Poles and values outside binary64 range produce a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::carlson_rj_complex",
    DimensionPolicy::Dimensionless,
    (x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64, p: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |x: num_complex::Complex64, y: num_complex::Complex64, z: num_complex::Complex64, p: num_complex::Complex64| {
        if !super::carlson_rj_branch(x, y, z, p) {
            return Err(OpError::Domain("R_J duplication requires nonnegative real parts for x, y, z and a positive real part for p, or an R_D degeneration".into()));
        }
        super::checked_carlson(super::carlson::carlson_rj_complex(x, y, z, p))
    }
);
