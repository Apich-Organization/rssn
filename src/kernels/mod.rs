//! # Numeric kernels
//!
//! Plain numerical routines with no dependency on the expression graph:
//! they take slices, closures and matrices and return numbers. The rule
//! sets in [`crate::rules`] wrap them as reduction kernels for the numeric
//! phase; the simulations in [`crate::sim`] call them directly.

/// Algebraic multigrid for non-symmetric systems (Ruge-Stueben, Petrov-Galerkin aggregation).
pub mod amg_general;
/// Interpolation and approximation: splines, PCHIP, B-splines, Chebyshev, AAA, Pade.
pub mod approx;
/// Complex balls with outward-rounded interval arithmetic (Krawczyk proofs).
pub mod ball;
pub mod calculus;
pub mod combinatorics;
pub mod complex;
pub mod computer_graphics;
/// Sequence and series convergence acceleration (Aitken, Richardson, Wynn).
pub mod convergence;
/// Coordinate transformations between Cartesian, polar, cylindrical and spherical systems.
pub mod coordinates;
/// Double-double real and complex arithmetic.
pub mod dd;
/// Dense linear algebra: LU, Cholesky, QR, SVD, symmetric and general eigenvalues.
pub mod dense;
/// Curvature of a metric given as a closure: Christoffel symbols, Riemann and Ricci tensors.
pub mod differential_geometry;
/// Lanczos and Arnoldi eigensolvers with restarts and shift-invert.
pub mod eigs;
/// Krylov-Schur, block Lanczos and Jacobi-Davidson eigensolvers, generalised problems.
pub mod eigs_ext;
pub mod error_correction;
/// FFT for arbitrary lengths and fast convolution.
pub mod fftconv;
pub mod finite_field;
pub mod fractal_geometry_and_chaos;
pub mod functional_analysis;
pub mod geometric_algebra;
/// Weighted shortest paths, and partition refinement, canonical labelling and
/// automorphism groups (`refine`, `canonical_labelling`, `isomorphic`, all
/// over the `ColourGraph` input type).
pub mod graph;
/// Exact homology with torsion, persistent homology, cubical complexes.
pub mod homology;
/// Homotopy continuation in projective coordinates with Cauchy endgame and path-jumping protection.
pub mod homotopy;
/// Smale alpha-theory certificates for homotopy endpoints.
pub mod homotopy_cert;
/// Complex-parameter pFq, 2F1, 1F1 and Tricomi U by Taylor continuation of the defining ODE.
pub mod hypergeometric;
/// Indefinite sums and products at non-integer arguments.
pub mod indefinite_sum;
pub mod integrate;
pub mod interpolate;
/// Structured (banded, sparse) Jacobians, colouring finite differences and banded LU.
pub mod jacobian;
/// Sparse CSR matrices and Krylov solvers with preconditioners.
pub mod krylov;
/// MINRES, BiCGSTAB(l), IDR(s) and flexible GMRES.
pub mod krylov_ext;
/// Semisimple Lie algebras: root systems, weights, representations; structure theory.
pub mod lie_structure;
/// Matrix functions: expm, sqrtm, logm, phi functions and exp(tA)v.
pub mod matfun;
pub mod matrix;
/// Meijer G-function by Slater's theorem and numerical Mellin-Barnes integration.
pub mod meijer;
/// Geometric multigrid (V/W/FMG) and smoothed-aggregation algebraic multigrid.
pub mod multigrid;
pub mod number_theory;
pub mod ode;
/// Adaptive, stiff, symplectic, BVP, DAE and DDE integrators.
pub mod ode_adaptive;
/// Exponential (ETDRK4, exp-Rosenbrock), IMEX, symplectic and variational integrators.
pub mod ode_ext;
/// Radau IIA(5/9/13) and variable-order BDF stiff integrators with dense output, events and DAEs.
pub mod ode_stiff;
/// Unconstrained, global and constrained optimisation and linear programming.
pub mod optim;
pub mod optimize;
/// Fill-reducing orderings: approximate minimum degree and nested dissection.
pub mod ordering;
/// Multilevel nested dissection (heavy-edge coarsening, FM refinement, vertex-cover separators).
pub mod ordering_ml;
/// Deterministic data-parallel helpers (ordered maps, fixed-tree sums, chunked loops).
pub mod par;
pub mod pde;
pub mod polynomial;
/// Exact linear algebra over the rationals.
pub mod qlinalg;
/// The quadratic sieve and its self-initialising variant for integer factorisation.
pub mod quadratic_sieve;
/// Adaptive Gauss-Kronrod, double-exponential, Clenshaw-Curtis, Filon and multi-dimensional quadrature.
pub mod quadrature;
/// Seeded deterministic random number generation.
pub mod random;
/// Real root isolation and refinement for polynomials.
pub mod real_roots;
/// Robust statistics, quantiles and kernel density estimation.
pub mod robust;
/// Scalar and system root finding, polynomial roots, homotopy continuation.
pub mod rootfind;
/// Stochastic differential equations: Euler-Maruyama, Milstein, SRK, Heun and adaptive steps.
pub mod sde;
/// Partial sums and accelerated infinite sums.
pub mod series;
/// Signal processing: filters, windows and spectra.
pub mod signal;
pub mod solve;
pub mod sparse;
pub mod special;
/// Elliptic integrals, hypergeometric functions (real and complex argument) and integer-order Bessel functions.
pub mod special_ext;
pub mod stats;
pub mod tensor;
/// Computational topology: simplicial complexes, Betti numbers, persistence.
pub mod topology;
pub mod transforms;
pub mod vector;
pub mod wavelets;

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
