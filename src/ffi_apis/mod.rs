//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]

#[macro_use]

pub mod macros;

pub mod common;
pub mod compute_cache_ffi;
pub mod compute_state_ffi;
pub mod constant_ffi;
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]
pub mod ffi_api;
pub mod jit_ffi;
pub mod numerical_calculus_ffi;
pub mod numerical_calculus_of_variations_ffi;
pub mod numerical_combinatorics_ffi;
pub mod numerical_complex_analysis_ffi;
pub mod numerical_computer_graphics_ffi;
pub mod numerical_convergence_ffi;
pub mod numerical_coordinates_ffi;
pub mod numerical_differential_geometry_ffi;
pub mod numerical_elementary_ffi;
pub mod numerical_error_correction_ffi;
pub mod numerical_finite_field_ffi;
pub mod numerical_fractal_geometry_and_chaos_ffi;
pub mod numerical_functional_analysis_ffi;
pub mod numerical_geometric_algebra_ffi;
pub mod numerical_graph_ffi;
pub mod numerical_integrate_ffi;
pub mod numerical_interpolate_ffi;
pub mod numerical_matrix_ffi;
pub mod numerical_multi_valued_ffi;
pub mod numerical_number_theory_ffi;
pub mod numerical_ode_ffi;
pub mod numerical_optimize_ffi;
pub mod numerical_physics_cfd_ffi;
pub mod numerical_physics_fea_ffi;
pub mod numerical_physics_ffi;
pub mod numerical_physics_md_ffi;
pub mod numerical_polynomial_ffi;
pub mod numerical_real_roots_ffi;
pub mod numerical_series_ffi;
pub mod numerical_signal_ffi;
pub mod numerical_solve_ffi;
pub mod numerical_sparse_ffi;
pub mod numerical_special_ffi;
pub mod numerical_stats_ffi;
pub mod numerical_tensor_ffi;
pub mod numerical_topology_ffi;
pub mod numerical_transforms_ffi;
pub mod numerical_vector_calculus_ffi;
pub mod numerical_vector_ffi;
pub mod physics_bem_ffi;
pub mod physics_cnm_ffi;
pub mod physics_em_ffi;
pub mod physics_fdm_ffi;
pub mod physics_fem_ffi;
pub mod physics_fvm_ffi;
pub mod physics_mm_ffi;
pub mod physics_mtm_ffi;
pub mod physics_rkm_ffi;
pub mod physics_sim_fdtd_ffi;
pub mod physics_sim_geodesic_ffi;
pub mod physics_sim_gpe_ffi;
pub mod physics_sim_ising_ffi;
pub mod physics_sim_linear_elasticity_ffi;
pub mod physics_sim_navier_stokes_ffi;
pub mod physics_sim_schrodinger_ffi;
pub mod physics_sm_ffi;
pub mod plugins_ffi;
pub mod symbolic_cad_ffi;
pub mod symbolic_calculus_ffi;
pub mod symbolic_calculus_of_variations_ffi;
pub mod symbolic_cas_foundations_ffi;
pub mod symbolic_classical_mechanics_ffi;
pub mod symbolic_combinatorics_ffi;
pub mod symbolic_complex_analysis_ffi;
pub mod symbolic_computer_graphics_ffi;
pub mod symbolic_convergence_ffi;
pub mod symbolic_coordinates_ffi;
pub mod symbolic_cryptography_ffi;
pub mod symbolic_differential_geometry_ffi;
pub mod symbolic_discrete_groups_ffi;
pub mod symbolic_electromagnetism_ffi;
pub mod symbolic_elementary_ffi;
pub mod symbolic_error_correction_ffi;
pub mod symbolic_error_correction_helper_ffi;
pub mod symbolic_finite_field_ffi;
pub mod symbolic_fractal_geometry_and_chaos_ffi;
pub mod symbolic_functional_analysis_ffi;
pub mod symbolic_geometric_algebra_ffi;
pub mod symbolic_graph_algorithms_ffi;
pub mod symbolic_graph_ffi;
pub mod symbolic_graph_isomorphism_and_coloring_ffi;
pub mod symbolic_graph_operations_ffi;
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
pub mod symbolic_proof_ffi;
pub mod symbolic_quantum_field_theory_ffi;
pub mod symbolic_quantum_mechanics_ffi;
pub mod symbolic_radicals_ffi;
pub mod symbolic_real_roots_ffi;
pub mod symbolic_relativity_ffi;
pub mod symbolic_rewriting_ffi;
pub mod symbolic_series_ffi;
pub mod symbolic_simplify_dag_ffi;
pub mod symbolic_simplify_ffi;
pub mod symbolic_solid_state_physics_ffi;
pub mod symbolic_solve_ffi;
pub mod symbolic_special_ffi;
pub mod symbolic_special_functions_ffi;
pub mod symbolic_stats_ffi;
pub mod symbolic_stats_inference_ffi;
pub mod symbolic_stats_information_theory_ffi;
pub mod symbolic_stats_probability_ffi;
pub mod symbolic_stats_regression_ffi;
pub mod symbolic_tensor_ffi;
pub mod symbolic_thermodynamics_ffi;
pub mod symbolic_topology_ffi;
pub mod symbolic_transforms_ffi;
pub mod symbolic_unit_unification_ffi;
pub mod symbolic_vector_calculus_ffi;
pub mod symbolic_vector_ffi;
pub mod nightly_ffi;
