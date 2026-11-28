// Numeric cut operations share their kernels with the symbolic complex rules.
use crate::kernels::complex_coords as branch_cuts;

binding!(
    complex_coord_arg_branch,
    "num_complex_coords_arg_branch",
    "numerical_special",
    "`num_complex_coords_arg_branch(z, theta0)` returns the argument of z in (theta0, theta0 + 2π], with theta0 in radians. The upper endpoint is used on the cut. At zero the result follows the signed-zero convention of atan2 before wrapping. Both arguments are dimensionless.",
    "crate::kernels::complex_coords::arg_branch",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, theta0: f64) -> f64,
    super::branch_cuts::arg_branch
);
binding!(
    complex_coord_log_cut,
    "num_complex_coords_log_cut",
    "numerical_special",
    "`num_complex_coords_log_cut(z, theta0)` returns ln(|z|) + i·arg_branch(z, theta0) for dimensionless nonzero z. The cut angle theta0 is in radians. The upper endpoint is used on the cut.",
    "crate::kernels::complex_coords::log_cut",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, theta0: f64) -> Result<num_complex::Complex64, OpError>,
    |z: num_complex::Complex64, theta0: f64| {
        if z == num_complex::Complex64::new(0.0, 0.0) {
            return Err(OpError::Domain("the branch logarithm requires nonzero z".into()));
        }
        Ok(super::branch_cuts::log_cut(z, theta0))
    }
);
binding!(
    complex_coord_sqrt_cut,
    "num_complex_coords_sqrt_cut",
    "numerical_special",
    "`num_complex_coords_sqrt_cut(z, theta0)` returns sqrt(|z|)·exp(i·arg_branch(z, theta0)/2), with theta0 in radians. The result has the square root of the input dimension. Zero maps to zero.",
    "crate::kernels::complex_coords::sqrt_cut",
    DimensionPolicy::Root,
    (z: num_complex::Complex64, theta0: f64) -> num_complex::Complex64,
    super::branch_cuts::sqrt_cut
);
binding!(
    complex_coord_pow_cut,
    "num_complex_coords_pow_cut",
    "numerical_special",
    "`num_complex_coords_pow_cut(z, w, theta0)` returns exp(w·log_cut(z, theta0)) for dimensionless nonzero base z and dimensionless exponent w. The cut angle theta0 is in radians. The upper endpoint is used on the cut.",
    "crate::kernels::complex_coords::pow_cut",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, w: num_complex::Complex64, theta0: f64) -> Result<num_complex::Complex64, OpError>,
    |z: num_complex::Complex64, w: num_complex::Complex64, theta0: f64| {
        if z == num_complex::Complex64::new(0.0, 0.0) {
            return Err(OpError::Domain("the branch power requires a nonzero base".into()));
        }
        Ok(super::branch_cuts::pow_cut(z, w, theta0))
    }
);
