// Portable algebraic covering data with packed ramification cycles.
use crate::kernels::complex_coords::RiemannSurfaceInfo;

record!(RiemannSurfaceInfo {
    degree: usize,
    branch_points: Vec<Complex64>,
    has_branch_at_infinity: bool,
    ramification_indices: Vec<Vec<usize>>,
    monodromy_generators: Vec<Vec<usize>>,
    genus: usize
});

binding!(
    riemann_hyperelliptic,
    "num_complex_surface_hyperelliptic",
    "numerical_special",
    "Returns normalized hyperelliptic covering topology for w^2 with given finite roots. Generators act on zero-based sheets, and ramification cycles are packed with offsets.",
    "crate::kernels::complex_coords::algebraic_riemann_surface_hyperelliptic",
    DimensionPolicy::Dimensionless,
    (roots: Vec<num_complex::Complex64>) -> super::RiemannSurfaceInfo,
    |roots: Vec<num_complex::Complex64>| {
        crate::kernels::complex_coords::algebraic_riemann_surface_hyperelliptic(&roots)
    }
);
binding!(
    riemann_superelliptic,
    "num_complex_surface_superelliptic",
    "numerical_special",
    "Returns normalized superelliptic covering topology for w^degree with given finite roots. The infinity permutation uses a positive loop in the 1/z coordinate.",
    "crate::kernels::complex_coords::algebraic_riemann_surface_superelliptic",
    DimensionPolicy::Dimensionless,
    (roots: Vec<num_complex::Complex64>, degree: usize) -> super::RiemannSurfaceInfo,
    |roots: Vec<num_complex::Complex64>, degree: usize| {
        crate::kernels::complex_coords::algebraic_riemann_surface_superelliptic(&roots,degree)
    }
);
binding!(
    riemann_monodromy,
    "num_complex_surface_monodromy_checked",
    "numerical_special",
    "Tracks polynomial fibers on a normalized closed piecewise-linear loop. Coefficient rows index w degree and columns index z degree in ascending order. The last vertex is connected back to the first. Returns zero-based sheet images or an error for ambiguous fibers.",
    "crate::kernels::complex_coords::monodromy_around_loop",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<Vec<num_complex::Complex64>>, path: Vec<num_complex::Complex64>, subdivisions: usize) -> Result<Vec<usize>, OpError>,
    |coefficients: Vec<Vec<num_complex::Complex64>>, path: Vec<num_complex::Complex64>, subdivisions: usize| {
        use crate::kernels::complex_coords::{MonodromyError,monodromy_checked};
        monodromy_checked(&coefficients,&path,subdivisions).map_err(|error| match error {
            MonodromyError::Budget=>OpError::BudgetExceeded,
            MonodromyError::Cancelled=>crate::cancellation::Cancelled.into(),
            MonodromyError::Inputs=>OpError::Domain("invalid algebraic covering path or polynomial".into()),
            MonodromyError::Resolution=>OpError::Domain("algebraic fiber roots or path are unresolved in Float".into()),
        })
    }
);
