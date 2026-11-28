// Named portable outputs for the algebraic coordinate transforms.
struct RiemannSpherePoint {
    x: f64,
    y: f64,
    z: f64,
}
record!(RiemannSpherePoint {
    x: f64,
    y: f64,
    z: f64
});

struct JoukowskyBranches {
    plus: num_complex::Complex64,
    minus: num_complex::Complex64,
}
record!(JoukowskyBranches {
    plus: num_complex::Complex64,
    minus: num_complex::Complex64
});

binding!(
    riemann_sphere_project,
    "num_complex_coords_to_riemann_sphere",
    "numerical_special",
    "`num_complex_coords_to_riemann_sphere(point)` projects a finite dimensionless complex point to a unit sphere record with x, y and z fields. Zero maps to the south pole.",
    "crate::kernels::complex_coords::to_riemann_sphere",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64) -> super::RiemannSpherePoint,
    |point: num_complex::Complex64| {
        let (x, y, z) = crate::kernels::complex_coords::to_riemann_sphere(point);
        super::RiemannSpherePoint { x, y, z }
    }
);
binding!(
    riemann_sphere_unproject,
    "num_complex_coords_from_riemann_sphere",
    "numerical_special",
    "`num_complex_coords_from_riemann_sphere(x, y, z)` inverts stereographic projection for a dimensionless unit sphere point. The north pole returns None. Invalid coordinates or an unrepresentable finite image produce a domain error.",
    "crate::kernels::complex_coords::from_riemann_sphere",
    DimensionPolicy::Dimensionless,
    (x: f64, y: f64, z: f64) -> Result<Option<num_complex::Complex64>, OpError>,
    |x: f64, y: f64, z: f64| {
        let length = x.hypot(y).hypot(z);
        if !length.is_finite() || (length - 1.0).abs() > 16.0 * f64::EPSILON {
            return Err(OpError::Domain("stereographic projection requires a point on the unit sphere".into()));
        }
        Ok(crate::kernels::complex_coords::from_riemann_sphere(x, y, z))
    }
);
binding!(
    joukowsky_inverse,
    "num_conformal_joukowsky_inverse",
    "numerical_special",
    "`num_conformal_joukowsky_inverse(point)` returns the two dimensionless inverse branches as named plus and minus fields, ordered by the principal square root of point² - 1. An unrepresentable branch produces a domain error.",
    "crate::kernels::conformal::joukowsky_inverse",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64) -> super::JoukowskyBranches,
    |point: num_complex::Complex64| {
        let (plus, minus) = crate::kernels::conformal::joukowsky_inverse(point);
        super::JoukowskyBranches { plus, minus }
    }
);
