// Numerical bounds in normalized polynomial coordinates and on the unit disk.
struct RadialBounds {
    lower: f64,
    upper: f64,
}
record!(RadialBounds {
    lower: f64,
    upper: f64
});

fn polynomial_radius(radius: Option<f64>) -> Result<f64, OpError> {
    radius
        .ok_or_else(|| OpError::Domain("expected a finite nonconstant polynomial".into()))
        .and_then(|radius| {
            if radius.is_finite() {
                Ok(radius)
            } else {
                Err(OpError::Domain(
                    "polynomial root radius exceeds binary64 range".into(),
                ))
            }
        })
}

const fn radial_bounds(bounds: (f64, f64)) -> RadialBounds {
    let (lower, upper) = bounds;
    RadialBounds { lower, upper }
}

fn invalid_disk_radius() -> OpError {
    OpError::Domain("unit disk radius must satisfy 0 <= r < 1".into())
}

binding!(
    cauchy_radius,
    "num_cauchy_root_bound",
    "numerical_vector",
    "Classic Cauchy upper bound from a complex coefficient vector in ascending order. Coefficients may share ordinate units. The polynomial coordinate and the result are dimensionless.",
    "crate::kernels::complex_roots::cauchy_root_bound",
    DimensionPolicy::UnitlessResult,
    (coefficients: Vec<num_complex::Complex64>) -> Result<f64, OpError>,
    |coefficients: Vec<num_complex::Complex64>| super::polynomial_radius(
        crate::kernels::complex_roots::cauchy_root_bound(&coefficients)
    )
);
binding!(
    sharp_cauchy_radius,
    "num_cauchy_sharp_root_bound",
    "numerical_vector",
    "Compute the sharp Cauchy radius from ascending complex coefficients with a common ordinate dimension. Returns zero for a nonconstant monomial. The polynomial coordinate and the result are dimensionless.",
    "crate::kernels::complex_roots::cauchy_sharp_root_bound",
    DimensionPolicy::UnitlessResult,
    (coefficients: Vec<num_complex::Complex64>) -> Result<f64, OpError>,
    |coefficients: Vec<num_complex::Complex64>| super::polynomial_radius(
        crate::kernels::complex_roots::cauchy_sharp_root_bound(&coefficients)
    )
);
binding!(
    fujiwara_radius,
    "num_fujiwara_root_bound",
    "numerical_vector",
    "Fujiwara upper bound from ascending complex coefficients with a common ordinate dimension. The polynomial coordinate and the result are dimensionless.",
    "crate::kernels::complex_roots::fujiwara_root_bound",
    DimensionPolicy::UnitlessResult,
    (coefficients: Vec<num_complex::Complex64>) -> Result<f64, OpError>,
    |coefficients: Vec<num_complex::Complex64>| super::polynomial_radius(
        crate::kernels::complex_roots::fujiwara_root_bound(&coefficients)
    )
);
binding!(
    enestrom_kakeya,
    "num_enestrom_kakeya_bounds",
    "numerical_vector",
    "Root annulus for at least two positive, finite real coefficients in ascending order. Coefficients may share ordinate units. The lower and upper root radii are dimensionless.",
    "crate::kernels::complex_roots::enestrom_kakeya_bounds",
    DimensionPolicy::UnitlessResult,
    (coefficients: Vec<f64>) -> Result<super::RadialBounds, OpError>,
    |coefficients: Vec<f64>| crate::kernels::complex_roots::enestrom_kakeya_bounds(&coefficients)
        .map(super::radial_bounds)
        .ok_or_else(|| OpError::Domain("root annulus requires at least two finite positive coefficients".into()))
);
binding!(
    koebe_growth,
    "num_koebe_growth_bounds",
    "numerical_special",
    "Koebe lower and upper bounds on |f(z)| for normalized univalent functions at |z|=r in [0,1).",
    "crate::kernels::univalent::koebe_growth_bounds",
    DimensionPolicy::Dimensionless,
    (r: f64) -> Result<super::RadialBounds, OpError>,
    |r: f64| crate::kernels::univalent::koebe_growth_bounds(r).map(super::radial_bounds).ok_or_else(super::invalid_disk_radius)
);
binding!(
    koebe_distortion,
    "num_koebe_distortion_bounds",
    "numerical_special",
    "Koebe lower and upper bounds on |f'(z)| for normalized univalent functions at |z|=r in [0,1).",
    "crate::kernels::univalent::koebe_distortion_bounds",
    DimensionPolicy::Dimensionless,
    (r: f64) -> Result<super::RadialBounds, OpError>,
    |r: f64| crate::kernels::univalent::koebe_distortion_bounds(r).map(super::radial_bounds).ok_or_else(super::invalid_disk_radius)
);
binding!(
    koebe_rotation,
    "num_koebe_rotation_bound",
    "numerical_special",
    "Goluzin's sharp upper bound in radians on the absolute continuous argument of f'(z) for |z|=r in [0,1).",
    "crate::kernels::univalent::koebe_rotation_bound",
    DimensionPolicy::Dimensionless,
    (r: f64) -> Result<f64, OpError>,
    |r: f64| crate::kernels::univalent::koebe_rotation_bound(r).ok_or_else(super::invalid_disk_radius)
);
