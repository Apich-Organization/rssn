// Portable conformal map and polynomial basis exports.
use crate::kernels::conformal::SchwarzChristoffelMap;
use crate::kernels::contour_deform::ContourIntegral;

impl Arg for SchwarzChristoffelMap {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        const FIELDS: &[&str] = &["prevertices", "interior_angles_pi", "scale", "shift"];
        let field = |name| model_field(value, name, FIELDS);
        Ok(Self {
            prevertices: Vec::<f64>::decode(field("prevertices")?, index)?,
            interior_angles_pi: Vec::<f64>::decode(field("interior_angles_pi")?, index)?,
            scale: Complex64::decode(field("scale")?, index)?,
            shift: Complex64::decode(field("shift")?, index)?,
        })
    }
}

record!(SchwarzChristoffelMap {
    prevertices: Vec<f64>,
    interior_angles_pi: Vec<f64>,
    scale: Complex64,
    shift: Complex64
});
record!(ContourIntegral {
    value: Complex64,
    error: f64,
    evaluations: usize,
    converged: bool
});

binding!(
    conformal_schwarz_christoffel_evaluate,
    "num_conformal_schwarz_christoffel_evaluate",
    "numerical_special",
    "Evaluates a normalized Schwarz–Christoffel map on the open upper half-plane. Provide the prepared map, point, positive absolute tolerance and panel budget. The result contains the complex value, absolute error estimate, evaluations and convergence.",
    "crate::kernels::conformal::SchwarzChristoffelMap::eval",
    DimensionPolicy::Dimensionless,
    (map: super::SchwarzChristoffelMap, point: num_complex::Complex64, tolerance: f64, panels: usize) -> Result<super::ContourIntegral, OpError>,
    |map: super::SchwarzChristoffelMap, point: num_complex::Complex64, tolerance: f64, panels: usize| {
        let result = map.eval_checked(point, tolerance, panels);
        crate::cancellation::checkpoint()?;
        result.map_err(OpError::Domain)
    }
);
binding!(
    conformal_bergman_basis,
    "num_conformal_bergman_polynomial_basis",
    "numerical_special",
    "Evaluates a dimensionless Bergman expansion with rows of ascending complex basis coefficients. Zero-pad lower-degree polynomials. The diagonal evaluates the basis once per row.",
    "crate::kernels::conformal::bergman_kernel_polynomial_basis",
    DimensionPolicy::Dimensionless,
    (basis: Vec<Vec<num_complex::Complex64>>, z: num_complex::Complex64, w: num_complex::Complex64) -> num_complex::Complex64,
    |basis: Vec<Vec<num_complex::Complex64>>, z: num_complex::Complex64, w: num_complex::Complex64| {
        crate::kernels::conformal::bergman_kernel_polynomial_basis(&basis, z, w)
    }
);
binding!(
    conformal_mobius_three_points,
    "num_conformal_mobius_three_points",
    "numerical_special",
    "Builds the projective Möbius map sending three distinct finite source points to three distinct finite target points. Coefficients are homogeneous.",
    "crate::kernels::conformal::mobius_from_three_points",
    DimensionPolicy::Dimensionless,
    (z1: num_complex::Complex64, z2: num_complex::Complex64, z3: num_complex::Complex64, w1: num_complex::Complex64, w2: num_complex::Complex64, w3: num_complex::Complex64) -> Result<super::Mobius, OpError>,
    |z1: num_complex::Complex64, z2: num_complex::Complex64, z3: num_complex::Complex64, w1: num_complex::Complex64, w2: num_complex::Complex64, w3: num_complex::Complex64| {
        crate::kernels::conformal::mobius_from_three_points(z1, z2, z3, w1, w2, w3).map_err(OpError::Domain)
    }
);
binding!(
    conformal_schwarz_christoffel_solve,
    "num_conformal_schwarz_christoffel_solve",
    "numerical_special",
    "Fits a reusable Schwarz–Christoffel map to a simple counterclockwise polygon in normalized coordinates. Pass a positive tolerance and Newton iteration budget.",
    "crate::kernels::conformal::solve_sc_uhp",
    DimensionPolicy::Dimensionless,
    (vertices: Vec<num_complex::Complex64>, max_iterations: usize, tolerance: f64) -> Result<super::SchwarzChristoffelMap, OpError>,
    |vertices: Vec<num_complex::Complex64>, max_iterations: usize, tolerance: f64| {
        let solution = crate::kernels::conformal::solve_sc_uhp(&vertices, max_iterations, tolerance);
        crate::cancellation::checkpoint()?;
        solution.map_err(OpError::Domain)
    }
);
