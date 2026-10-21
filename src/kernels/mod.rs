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
/// Interpolation and approximation: splines, PCHIP, B-splines, Chebyshev, AAA, Pade.
pub mod approx;
/// Dense linear algebra: LU, Cholesky, QR, SVD, symmetric and general eigenvalues.
pub mod dense;
/// FFT for arbitrary lengths and fast convolution.
pub mod fftconv;
/// Sparse CSR matrices and Krylov solvers with preconditioners.
pub mod krylov;
/// Radau IIA(5/9/13) and variable-order BDF stiff integrators with dense output, events and DAEs.
pub mod ode_stiff;
/// Adaptive, stiff, symplectic, BVP, DAE and DDE integrators.
pub mod ode_adaptive;
/// Unconstrained, global and constrained optimisation and linear programming.
pub mod optim;
/// Adaptive Gauss-Kronrod, double-exponential, Clenshaw-Curtis, Filon and multi-dimensional quadrature.
pub mod quadrature;
/// Seeded deterministic random number generation.
pub mod random;
/// Robust statistics, quantiles and kernel density estimation.
pub mod robust;
/// Scalar and system root finding, polynomial roots, homotopy continuation.
pub mod rootfind;
/// Elliptic integrals, hypergeometric functions (real and complex argument) and integer-order Bessel functions.
pub mod special_ext;
/// Fill-reducing orderings: approximate minimum degree and nested dissection.
pub mod ordering;
/// Structured (banded, sparse) Jacobians, colouring finite differences and banded LU.
pub mod jacobian;
/// Krylov-Schur, block Lanczos and Jacobi-Davidson eigensolvers, generalised problems.
pub mod eigs_ext;
/// Complex-parameter pFq, 2F1, 1F1 and Tricomi U by Taylor continuation of the defining ODE.
pub mod hypergeometric;
/// Algebraic multigrid for non-symmetric systems (Ruge-Stueben, Petrov-Galerkin aggregation).
pub mod amg_general;
/// Multilevel nested dissection (heavy-edge coarsening, FM refinement, vertex-cover separators).
pub mod ordering_ml;
/// MINRES, BiCGSTAB(l), IDR(s) and flexible GMRES.
pub mod krylov_ext;
/// Lanczos and Arnoldi eigensolvers with restarts and shift-invert.
pub mod eigs;
/// Matrix functions: expm, sqrtm, logm, phi functions and exp(tA)v.
pub mod matfun;
/// Geometric multigrid (V/W/FMG) and smoothed-aggregation algebraic multigrid.
pub mod multigrid;
/// Exponential (ETDRK4, exp-Rosenbrock), IMEX, symplectic and variational integrators.
pub mod ode_ext;
/// Stochastic differential equations: Euler-Maruyama, Milstein, SRK, Heun and adaptive steps.
pub mod sde;
/// Homotopy continuation in projective coordinates with Cauchy endgame and path-jumping protection.
pub mod homotopy;

/// Exact IEEE-754 equality (`-0.0 == 0.0`, NaN never equal).
///
/// For the few places where bit-exact equality is the intended predicate
/// (sentinels, "no change after rounding" tests such as `|h| + s == s`).
#[must_use]
pub(crate) fn exact_eq(
    a: f64,
    b: f64,
) -> bool {
    a.partial_cmp(&b) == Some(std::cmp::Ordering::Equal)
}
