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

struct PolarCoordinates {
    angle: f64,
    radius: f64,
}
record!(PolarCoordinates {
    angle: f64,
    radius: f64
});

binding!(
    complex_coordinate_to_polar,
    "num_complex_coords_to_polar",
    "numerical_special",
    "Returns a measured radius and the principal angle in radians. Signed zero selects the angle on the negative real axis.",
    "crate::kernels::complex_coords::to_polar",
    DimensionPolicy::PolarCoordinates,
    (point: num_complex::Complex64) -> super::PolarCoordinates,
    |point: num_complex::Complex64| {
        let (radius, angle) = crate::kernels::complex_coords::to_polar(point);
        super::PolarCoordinates { angle, radius }
    }
);
binding!(
    complex_coordinate_from_polar,
    "num_complex_coords_from_polar",
    "numerical_special",
    "Reconstructs a complex point in radius units from a nonnegative radius and an angle in radians.",
    "crate::kernels::complex_coords::from_polar",
    DimensionPolicy::PreserveWithUnitlessSecond,
    (radius: f64, angle: f64) -> Result<num_complex::Complex64, OpError>,
    |radius: f64, angle: f64| {
        if radius < 0.0 {
            return Err(OpError::Domain("polar radius must be nonnegative".into()));
        }
        Ok(crate::kernels::complex_coords::from_polar(radius, angle))
    }
);
binding!(
    complex_coordinate_abs_squared,
    "num_complex_coords_abs_sq",
    "numerical_special",
    "Returns the squared modulus in squared coordinate units. An unrepresentable result is a domain error.",
    "crate::kernels::complex_coords::abs_sq",
    DimensionPolicy::Square,
    (point: num_complex::Complex64) -> f64,
    crate::kernels::complex_coords::abs_sq
);
binding!(
    complex_coordinate_sign,
    "num_complex_coords_sign",
    "numerical_special",
    "Returns the dimensionless unit direction of a finite nonzero point; zero maps to zero.",
    "crate::kernels::complex_coords::sign_complex",
    DimensionPolicy::UnitlessResult,
    (point: num_complex::Complex64) -> num_complex::Complex64,
    crate::kernels::complex_coords::sign_complex
);
binding!(
    complex_coordinate_unwrap_phase,
    "num_complex_coords_unwrap_phase",
    "numerical_special",
    "Unwraps phases in radians of nonzero samples with a shared coordinate dimension. Consecutive increments are in (−π, π].",
    "crate::kernels::complex_coords::unwrap_phase",
    DimensionPolicy::UnitlessResult,
    (samples: Vec<num_complex::Complex64>) -> Result<Vec<f64>, OpError>,
    |samples: Vec<num_complex::Complex64>| {
        if samples.contains(&num_complex::Complex64::new(0.0, 0.0)) {
            return Err(OpError::Domain("phase is undefined at a zero sample".into()));
        }
        Ok(crate::kernels::complex_coords::unwrap_phase(&samples))
    }
);
binding!(
    complex_coordinate_continue_root,
    "num_complex_coords_continue_root",
    "numerical_special",
    "Tracks the continuous degree-th root of a dimensionless sampled path away from zero. A positive degree and increments of less than π between samples select a branch.",
    "crate::kernels::complex_coords::analytic_continuation_root",
    DimensionPolicy::Dimensionless,
    (points: Vec<num_complex::Complex64>, degree: u32, initial_branch: i32)
        -> Result<Vec<num_complex::Complex64>, OpError>,
    |points: Vec<num_complex::Complex64>, degree: u32, initial_branch: i32| {
        if degree == 0 || points.contains(&num_complex::Complex64::new(0.0, 0.0)) {
            return Err(OpError::Domain("root continuation needs a positive degree and a nonzero path".into()));
        }
        Ok(crate::kernels::complex_coords::analytic_continuation_root(&points, degree, initial_branch))
    }
);
