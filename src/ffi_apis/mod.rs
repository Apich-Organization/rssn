//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(clippy::no_mangle_with_rust_abi)]

#[macro_use]
pub mod macros;
pub mod common;
pub mod compute_cache_ffi;
pub mod compute_state_ffi;
pub mod constant_ffi;
pub mod ffi_api;
pub mod symbolic_calculus_ffi;
pub mod symbolic_cas_foundations_ffi;
pub mod symbolic_complex_analysis_ffi;
pub mod symbolic_convergence_ffi;
pub mod symbolic_coordinates_ffi;
pub mod symbolic_differential_geometry_ffi;
pub mod symbolic_discrete_groups_ffi;
pub mod symbolic_elementary_ffi;
pub mod symbolic_finite_field_ffi;
pub mod symbolic_fractal_geometry_and_chaos_ffi;
pub mod symbolic_functional_analysis_ffi;
pub mod symbolic_geometric_algebra_ffi;
pub mod symbolic_grobner_ffi;
pub mod symbolic_group_theory_ffi;
pub mod symbolic_handles_ffi;
pub mod symbolic_integral_equations_ffi;
pub mod symbolic_integration_ffi;
pub mod symbolic_lie_groups_ffi;
pub mod symbolic_logic_ffi;
pub mod symbolic_matrix_ffi;
pub mod symbolic_multi_valued_ffi;
pub mod symbolic_number_theory_ffi;
pub mod symbolic_numeric_ffi;
pub mod symbolic_ode_ffi;
pub mod symbolic_optimize_ffi;
pub mod symbolic_pde_ffi;
pub mod symbolic_poly_factorization_ffi;
pub mod symbolic_polynomial_ffi;
pub mod symbolic_radicals_ffi;
pub mod symbolic_real_roots_ffi;
pub mod symbolic_rewriting_ffi;
pub mod symbolic_series_ffi;
pub mod symbolic_simplify_dag_ffi;
pub mod symbolic_simplify_ffi;
pub mod symbolic_solve_ffi;
pub mod symbolic_tensor_ffi;
pub mod symbolic_topology_ffi;
pub mod symbolic_unit_unification_ffi;
pub mod symbolic_vector_calculus_ffi;
pub mod symbolic_vector_ffi;
