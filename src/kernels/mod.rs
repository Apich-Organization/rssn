//! # Numeric kernels
//!
//! Plain numerical routines with no dependency on the expression graph:
//! they take slices, closures and matrices and return numbers. The rule
//! sets in [`crate::rules`] wrap them as reduction kernels for the numeric
//! phase; the simulations in [`crate::sim`] call them directly.

pub mod calculus;
pub mod combinatorics;
/// Coordinate transformations between Cartesian, polar, cylindrical and spherical systems.
pub mod coordinates;
/// Sequence and series convergence acceleration (Aitken, Richardson, Wynn).
pub mod convergence;
pub mod complex;
pub mod computer_graphics;
pub mod error_correction;
pub mod finite_field;
pub mod fractal_geometry_and_chaos;
/// Curvature of a metric given as a closure: Christoffel symbols, Riemann and Ricci tensors.
pub mod differential_geometry;
pub mod functional_analysis;
pub mod geometric_algebra;
pub mod graph;
/// Indefinite sums and products at non-integer arguments.
pub mod indefinite_sum;
pub mod integrate;
pub mod interpolate;
pub mod matrix;
pub mod number_theory;
pub mod ode;
pub mod optimize;
pub mod pde;
pub mod polynomial;
/// Real root isolation and refinement for polynomials.
pub mod real_roots;
/// Partial sums and accelerated infinite sums.
pub mod series;
/// Signal processing: filters, windows and spectra.
pub mod signal;
pub mod solve;
pub mod sparse;
pub mod special;
pub mod stats;
pub mod tensor;
/// Computational topology: simplicial complexes, Betti numbers, persistence.
pub mod topology;
/// Exact homology with torsion, persistent homology, cubical complexes.
pub mod homology;
/// Exact linear algebra over the rationals.
pub mod qlinalg;
/// Semisimple Lie algebras: root systems, weights, representations; structure theory.
pub mod lie_structure;
pub mod transforms;
pub mod vector;
