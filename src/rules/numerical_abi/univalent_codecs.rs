// Checked diagnostics on the normalized, dimensionless unit disk.
use crate::kernels::univalent;

fn checked_univalent_result<T>(result: Result<T, String>) -> Result<T, OpError> {
    crate::cancellation::checkpoint()?;
    result.map_err(OpError::Domain)
}

binding!(
    univalent_bieberbach,
    "num_univalent_bieberbach",
    "numerical_special",
    "Checks the necessary Bieberbach–de Branges coefficient inequalities for dimensionless Taylor coefficients [a0,a1,...]. False rules out the normalized class within the nonnegative tolerance.",
    "crate::kernels::univalent::check_bieberbach_de_branges",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, tolerance: f64) -> Result<bool, OpError>,
    |coefficients: Vec<num_complex::Complex64>, tolerance: f64| {
        if tolerance < 0.0 || coefficients.get(1).is_none_or(|v| *v == num_complex::Complex64::default()) {
            return Err(OpError::Domain("Bieberbach check needs a nonzero a1 and nonnegative tolerance".into()));
        }
        Ok(super::univalent::check_bieberbach_de_branges(&coefficients, tolerance))
    }
);
binding!(
    univalent_fekete_szego,
    "num_univalent_fekete_szego",
    "numerical_special",
    "Returns |a3 - mu*a2^2| for normalized dimensionless Taylor coefficients.",
    "crate::kernels::univalent::fekete_szego_functional",
    DimensionPolicy::Dimensionless,
    (a2: num_complex::Complex64, a3: num_complex::Complex64, mu: num_complex::Complex64) -> f64,
    super::univalent::fekete_szego_functional
);
binding!(
    univalent_hankel,
    "num_univalent_hankel",
    "numerical_special",
    "Returns the dimensionless Hankel determinant with entry coefficients[start+i+j]. Order zero returns one. Orders through 256 are supported.",
    "crate::kernels::univalent::hankel_determinant",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, order: usize, start: usize) -> Result<num_complex::Complex64, OpError>,
    |coefficients: Vec<num_complex::Complex64>, order: usize, start: usize| {
        super::checked_univalent_result(super::univalent::hankel_determinant(&coefficients, order, start))
    }
);
binding!(
    univalent_grunsky_matrix,
    "num_univalent_grunsky_matrix",
    "numerical_special",
    "Builds the symmetric normalized complex Grunsky matrix. Requires 2*order+2 dimensionless Taylor coefficients for a positive order through 92.",
    "crate::kernels::univalent::grunsky_matrix",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, order: usize) -> Result<Vec<Vec<num_complex::Complex64>>, OpError>,
    |coefficients: Vec<num_complex::Complex64>, order: usize| {
        super::checked_univalent_result(super::univalent::grunsky_matrix(&coefficients, order))
    }
);
binding!(
    univalent_grunsky_check,
    "num_univalent_grunsky_check",
    "numerical_special",
    "Tests the complex bilinear Grunsky inequality for a supplied vector and checks a power-iteration estimate of the spectral norm. This is a finite-dimensional numerical diagnostic.",
    "crate::kernels::univalent::check_grunsky_inequality",
    DimensionPolicy::Dimensionless,
    (matrix: Vec<Vec<num_complex::Complex64>>, test_vector: Vec<num_complex::Complex64>, tolerance: f64) -> Result<bool, OpError>,
    |matrix: Vec<Vec<num_complex::Complex64>>, test_vector: Vec<num_complex::Complex64>, tolerance: f64| {
        super::checked_univalent_result(super::univalent::check_grunsky_inequality(&matrix, &test_vector, tolerance))
    }
);
binding!(
    univalent_becker,
    "num_univalent_becker",
    "numerical_special",
    "Tests Becker’s pointwise dimensionless inequality at a point in the unit disk with a nonzero first derivative.",
    "crate::kernels::univalent::check_becker_univalence",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64, first: num_complex::Complex64, second: num_complex::Complex64) -> bool,
    super::univalent::check_becker_univalence
);
binding!(
    univalent_loewner_step,
    "num_univalent_loewner_step",
    "numerical_special",
    "Takes one radial Loewner RK4 step in the normalized disk. The driving angle is in radians, the step is nonnegative and each sampled state must stay away from the driving pole.",
    "crate::kernels::univalent::loewner_slit_step",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64, driving_angle: f64, step: f64) -> Result<num_complex::Complex64, OpError>,
    |point: num_complex::Complex64, driving_angle: f64, step: f64| {
        if point.norm() >= 1.0 || step < 0.0 {
            return Err(OpError::Domain("Loewner step needs a point in the open unit disk and nonnegative step".into()));
        }
        super::univalent::loewner_step_checked(point, driving_angle, step)
            .ok_or_else(|| OpError::Domain("Loewner step is undefined at a driving pole or exceeds Float range".into()))
    }
);
