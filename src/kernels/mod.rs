//! # Numeric kernels
//!
//! Plain numerical routines with no dependency on the expression graph:
//! they take slices, closures and matrices and return numbers. The rule
//! sets in [`crate::rules`] wrap them as reduction kernels for the numeric
//! phase; the simulations in [`crate::sim`] call them directly.

pub mod amg_general;
pub mod analytic_continuation;
/// Interpolation and approximation: splines, PCHIP, B-splines, Chebyshev, AAA, Pade.
pub mod approx;
pub mod autodiff;
pub mod ball;
pub mod bsr;
pub mod calculus;
pub mod combinatorics;
pub mod complex;
pub mod complex_ball;
pub mod complex_coords;
pub mod complex_dynamics;
pub mod complex_roots;
pub mod complex_transforms;
pub mod computer_graphics;
pub mod conformal;
pub mod contour_deform;
pub mod convergence;
pub mod coordinates;
pub mod dd;
pub mod dense;
/// Curvature of a metric given as a closure: Christoffel symbols, Riemann and Ricci tensors.
pub mod differential_geometry;
/// Univariate and multivariate probability distributions.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::use_self,
    clippy::needless_pass_by_value,
    clippy::too_many_lines
)]
pub mod distributions;
pub mod distributions_ext;
pub mod eigs;
pub mod eigs_ext;
pub mod error_correction;
/// FFT for arbitrary lengths and fast convolution.
pub mod fftconv;
pub mod finite_field;
pub mod fractal_geometry_and_chaos;
pub mod fresnel;
pub mod functional_analysis;
pub mod geometric_algebra;
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    clippy::missing_const_for_fn,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::suboptimal_flops,
    clippy::needless_pass_by_value
)]
pub mod geometry;
pub mod graph;
pub mod holonomic;
pub mod homology;
pub mod homotopy;
pub mod homotopy_cert;
pub mod hypergeometric;
/// Indefinite sums and products at non-integer arguments.
pub mod indefinite_sum;
pub mod integrate;
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    clippy::missing_const_for_fn,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::suboptimal_flops,
    clippy::needless_pass_by_value,
    clippy::doc_markdown,
    clippy::must_use_candidate
)]
pub mod interp_ext;
pub mod interpolate;
pub mod jacobian;
pub mod krylov;
pub mod krylov_ext;
/// Semisimple Lie algebras: root systems, weights, representations; structure theory.
pub mod lie_structure;
pub mod lowrank;
pub mod matfun;
pub mod matrix;
/// Meijer G-function by Slater's theorem and numerical Mellin-Barnes integration.
pub mod meijer;
/// Division-free modular arithmetic (Montgomery and Barrett) shared by every
/// kernel that works modulo a machine-word modulus.
pub mod modular;
pub mod multigrid;
/// Unconstrained, global and constrained optimisation and linear programming.
/// Matrix-free trust-region Newton-CG with Hessian-vector actions.
pub mod newton_cg;
/// L-stable SDIRK3(2) integrator for stiff initial-value problems.
/// Extended number theory and exact combinatorics kernels.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::single_call_fn
)]
pub mod number_ext;
pub mod number_theory;
pub mod ode;
pub use ode::OdeMethod;
pub use ode::OdeStepMode;
pub use ode::solve_ode;
pub use ode::solve_ode_events;
pub mod ode_adaptive;
mod ode_contract;
pub mod ode_ext;
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::type_complexity,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::excessive_precision
)]
pub mod ode_geom;
pub mod ode_sdirk;
pub mod ode_stiff;
pub mod online_stats;
pub mod optim;
pub mod optimize;
pub mod ordering;
pub mod ordering_ml;
pub mod oscillatory;
/// Deterministic data-parallel helpers (ordered maps, fixed-tree sums, chunked loops).
pub mod par;
pub mod pde;
pub mod permgroup;
pub mod polynomial;
pub mod qlinalg;
pub mod quadratic_sieve;
pub mod quadrature;
pub mod random;
/// Real root isolation and refinement for polynomials.
pub mod real_roots;
/// Bounded storage for reusable numerical plans and scratch buffers.
pub mod reuse;
pub mod riemann_hilbert;
/// Robust statistics, quantiles and kernel density estimation.
pub mod robust;
pub mod rootfind;
pub mod schur;
pub mod sde;
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::use_self
)]
pub mod sde_ext;
/// Partial sums and accelerated infinite sums.
pub mod series;
pub mod signal;
/// Runtime-dispatched SIMD BLAS-1 and complex-vector microkernels.
pub mod simd;
pub mod simplex_quadrature;
pub mod solve;
pub mod sparse;
pub mod sparse_exact;
pub mod special;
pub mod special_complex;
pub mod special_complex_ext;
pub mod special_ext;
pub mod stats;
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps
)]
pub mod stats_ext;
pub mod steepest_descent;
pub mod tensor;
pub mod topology;
pub mod transforms;
pub mod univalent;
pub mod vector;
pub mod wavelets;

/// The error type of the stochastic integrators of [`sde`].
///
/// Re-exported at the module root because every caller that uses an SDE scheme
/// matches on it; the schemes themselves stay under [`sde`].
pub use sde::SdeError;

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

/// Tensor-train compression and contractions.
pub mod tensor_train;

pub mod fft;

/// Scalar complex arithmetic and transcendental functions.
pub mod complex_scalar;
/// Descriptive statistics on floating and exact samples.
pub mod descriptive;
/// Exact integer arithmetic and factorization.
pub mod number_exact;
/// Shared errors for low-level numerical routines.
pub mod numeric_common;
/// Native CSR arithmetic, solvers and preconditioners.
pub mod sparse_native;
/// Exact special-function polynomial recurrences.
pub mod special_exact;

/// Generic tensor contractions with caller-provided scalar operations.
pub mod contractions;

/// Geometry import and repair for triangle meshes.
///
/// Includes STL and STEP readers, mesh repair, and Boolean operations.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod cad;
/// Analysis and design of control systems.
///
/// Includes state space, transfer functions, Riccati equations, LQR/LQG/MPC, PID tuning,
/// frequency response, and system identification.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod control;
/// Solvers for convex optimisation.
///
/// Supports dense and sparse quadratic programs, second-order cone programs,
/// and semidefinite programs.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod convex;
/// Generalised linear models fitted by iteratively reweighted least squares.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod glm;
/// Markov chain Monte Carlo samplers and convergence diagnostics.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod mcmc;
/// Mesh generation and quality refinement.
///
/// Includes constrained Delaunay triangulation, tetrahedralisation, and advancing front.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod meshgen;
/// Mixed-integer linear programming by branch and bound, and sparse SQP.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod milp;
/// Arbitrary-precision floating point and ball arithmetic.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod mp;
/// Survival analysis: Kaplan-Meier, Nelson-Aalen, log-rank tests and Cox regression.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod survival;
/// Time-series models: autocorrelation, ARIMA, exponential smoothing and state-space filtering.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod timeseries;
