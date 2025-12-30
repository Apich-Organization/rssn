//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]
// In ffi_apis, we use raw pointers to pass data to and from Rust and C.
// clippy::not_unsafe_ptr_arg_deref is triggered by this.
#![allow(
    clippy::not_unsafe_ptr_arg_deref
)]
// This is enforced by clippy::nursery. It has too high false positive rate.
#![allow(
    clippy::option_if_let_else
)]

#[macro_use]
/// FFI macros.

pub mod macros;

/// Common FFI utilities.
pub mod common;
/// FFI for compute cache.
pub mod compute_cache_ffi;
/// FFI for compute state.
pub mod compute_state_ffi;
/// FFI for constants.
pub mod constant_ffi;
/// General FFI API (deprecated).
pub mod ffi_api;
/// FFI for JIT compilation.
pub mod jit_ffi;
/// FFI for nightly features.
pub mod nightly_ffi;
/// FFI for numerical calculus.
pub mod numerical_calculus_ffi;
/// FFI for numerical calculus of variations.
pub mod numerical_calculus_of_variations_ffi;
/// FFI for numerical combinatorics.
pub mod numerical_combinatorics_ffi;
/// FFI for numerical complex analysis.
pub mod numerical_complex_analysis_ffi;
/// FFI for numerical computer graphics.
pub mod numerical_computer_graphics_ffi;
/// FFI for numerical convergence.
pub mod numerical_convergence_ffi;
/// FFI for numerical coordinates.
pub mod numerical_coordinates_ffi;
/// FFI for numerical differential geometry.
pub mod numerical_differential_geometry_ffi;
/// FFI for numerical elementary functions.
pub mod numerical_elementary_ffi;
/// FFI for numerical error correction.
pub mod numerical_error_correction_ffi;
/// FFI for numerical finite field.
pub mod numerical_finite_field_ffi;
/// FFI for numerical fractal geometry and chaos.
pub mod numerical_fractal_geometry_and_chaos_ffi;
/// FFI for numerical functional analysis.
pub mod numerical_functional_analysis_ffi;
/// FFI for numerical geometric algebra.
pub mod numerical_geometric_algebra_ffi;
/// FFI for numerical graph theory.
pub mod numerical_graph_ffi;
/// FFI for numerical integration.
pub mod numerical_integrate_ffi;
/// FFI for numerical interpolation.
pub mod numerical_interpolate_ffi;
/// FFI for numerical matrix operations.
pub mod numerical_matrix_ffi;
/// FFI for numerical multi-valued functions.
pub mod numerical_multi_valued_ffi;
/// FFI for numerical number theory.
pub mod numerical_number_theory_ffi;
/// FFI for numerical ODEs.
pub mod numerical_ode_ffi;
/// FFI for numerical optimization.
pub mod numerical_optimize_ffi;
/// FFI for numerical CFD physics simulations.
pub mod numerical_physics_cfd_ffi;
/// FFI for numerical FEA physics simulations.
pub mod numerical_physics_fea_ffi;
/// FFI for numerical physics.
pub mod numerical_physics_ffi;
/// FFI for numerical molecular dynamics physics simulations.
pub mod numerical_physics_md_ffi;
/// FFI for numerical polynomials.
pub mod numerical_polynomial_ffi;
/// FFI for numerical real roots.
pub mod numerical_real_roots_ffi;
/// FFI for numerical series.
pub mod numerical_series_ffi;
/// FFI for numerical signal processing.
pub mod numerical_signal_ffi;
/// FFI for numerical solvers.
pub mod numerical_solve_ffi;
/// FFI for numerical sparse matrices.
pub mod numerical_sparse_ffi;
/// FFI for numerical special functions.
pub mod numerical_special_ffi;
/// FFI for numerical statistics.
pub mod numerical_stats_ffi;
/// FFI for numerical tensor operations.
pub mod numerical_tensor_ffi;
/// FFI for numerical topology.
pub mod numerical_topology_ffi;
/// FFI for numerical transforms.
pub mod numerical_transforms_ffi;
/// FFI for numerical vector calculus.
pub mod numerical_vector_calculus_ffi;
/// FFI for numerical vector operations.
pub mod numerical_vector_ffi;
/// FFI for BEM physics simulations.
pub mod physics_bem_ffi;
/// FFI for CNM physics simulations.
pub mod physics_cnm_ffi;
/// FFI for EM physics simulations.
pub mod physics_em_ffi;
/// FFI for FDM physics simulations.
pub mod physics_fdm_ffi;
/// FFI for FEM physics simulations.
pub mod physics_fem_ffi;
/// FFI for FVM physics simulations.
pub mod physics_fvm_ffi;
/// FFI for MM physics simulations.
pub mod physics_mm_ffi;
/// FFI for MTM physics simulations.
pub mod physics_mtm_ffi;
/// FFI for RKM physics simulations.
pub mod physics_rkm_ffi;
/// FFI for FDTD physics simulations.
pub mod physics_sim_fdtd_ffi;
/// FFI for geodesic relativity simulations.
pub mod physics_sim_geodesic_ffi;
/// FFI for GPE superfluidity simulations.
pub mod physics_sim_gpe_ffi;
/// FFI for Ising statistical simulations.
pub mod physics_sim_ising_ffi;
/// FFI for linear elasticity simulations.
pub mod physics_sim_linear_elasticity_ffi;
/// FFI for Navier-Stokes fluid simulations.
pub mod physics_sim_navier_stokes_ffi;
/// FFI for Schrodinger quantum simulations.
pub mod physics_sim_schrodinger_ffi;
/// FFI for SM physics simulations.
pub mod physics_sm_ffi;
/// FFI for plugins.
pub mod plugins_ffi;
/// FFI for symbolic CAD operations.
pub mod symbolic_cad_ffi;
/// FFI for symbolic calculus.
pub mod symbolic_calculus_ffi;
/// FFI for symbolic calculus of variations.
pub mod symbolic_calculus_of_variations_ffi;
/// FFI for symbolic CAS foundations.
pub mod symbolic_cas_foundations_ffi;
/// FFI for symbolic classical mechanics.
pub mod symbolic_classical_mechanics_ffi;
/// FFI for symbolic combinatorics.
pub mod symbolic_combinatorics_ffi;
/// FFI for symbolic complex analysis.
pub mod symbolic_complex_analysis_ffi;
/// FFI for symbolic computer graphics.
pub mod symbolic_computer_graphics_ffi;
/// FFI for symbolic convergence.
pub mod symbolic_convergence_ffi;
/// FFI for symbolic coordinates.
pub mod symbolic_coordinates_ffi;
/// FFI for symbolic cryptography.
pub mod symbolic_cryptography_ffi;
/// FFI for symbolic differential geometry.
pub mod symbolic_differential_geometry_ffi;
/// FFI for symbolic discrete groups.
pub mod symbolic_discrete_groups_ffi;
/// FFI for symbolic electromagnetism.
pub mod symbolic_electromagnetism_ffi;
/// FFI for symbolic elementary functions.
pub mod symbolic_elementary_ffi;
/// FFI for symbolic error correction.
pub mod symbolic_error_correction_ffi;
/// FFI for symbolic error correction helpers.
pub mod symbolic_error_correction_helper_ffi;
/// FFI for symbolic finite field.
pub mod symbolic_finite_field_ffi;
/// FFI for symbolic fractal geometry and chaos.
pub mod symbolic_fractal_geometry_and_chaos_ffi;
/// FFI for symbolic functional analysis.
pub mod symbolic_functional_analysis_ffi;
/// FFI for symbolic geometric algebra.
pub mod symbolic_geometric_algebra_ffi;
/// FFI for symbolic graph algorithms.
pub mod symbolic_graph_algorithms_ffi;
/// FFI for symbolic graph theory.
pub mod symbolic_graph_ffi;
/// FFI for symbolic graph isomorphism and coloring.
pub mod symbolic_graph_isomorphism_and_coloring_ffi;
/// FFI for symbolic graph operations.
pub mod symbolic_graph_operations_ffi;
/// FFI for symbolic Grobner bases.
pub mod symbolic_grobner_ffi;
/// FFI for symbolic group theory.
pub mod symbolic_group_theory_ffi;
/// FFI for symbolic handles.
pub mod symbolic_handles_ffi;
/// FFI for symbolic integral equations.
pub mod symbolic_integral_equations_ffi;
/// FFI for symbolic integration.
pub mod symbolic_integration_ffi;
/// FFI for symbolic Lie groups.
pub mod symbolic_lie_groups_ffi;
/// FFI for symbolic logic.
pub mod symbolic_logic_ffi;
/// FFI for symbolic matrix operations.
pub mod symbolic_matrix_ffi;
/// FFI for symbolic multi-valued functions.
pub mod symbolic_multi_valued_ffi;
/// FFI for symbolic number theory.
pub mod symbolic_number_theory_ffi;
/// FFI for symbolic numeric operations.
pub mod symbolic_numeric_ffi;
/// FFI for symbolic ODEs.
pub mod symbolic_ode_ffi;
/// FFI for symbolic optimization.
pub mod symbolic_optimize_ffi;
/// FFI for symbolic PDEs.
pub mod symbolic_pde_ffi;
/// FFI for symbolic polynomial factorization.
pub mod symbolic_poly_factorization_ffi;
/// FFI for symbolic polynomials.
pub mod symbolic_polynomial_ffi;
/// FFI for symbolic proofs.
pub mod symbolic_proof_ffi;
/// FFI for symbolic quantum field theory.
pub mod symbolic_quantum_field_theory_ffi;
/// FFI for symbolic quantum mechanics.
pub mod symbolic_quantum_mechanics_ffi;
/// FFI for symbolic radicals.
pub mod symbolic_radicals_ffi;
/// FFI for symbolic real roots.
pub mod symbolic_real_roots_ffi;
/// FFI for symbolic relativity.
pub mod symbolic_relativity_ffi;
/// FFI for symbolic rewriting.
pub mod symbolic_rewriting_ffi;
/// FFI for symbolic series.
pub mod symbolic_series_ffi;
/// FFI for symbolic DAG simplification.
pub mod symbolic_simplify_dag_ffi;
/// FFI for symbolic simplification.
pub mod symbolic_simplify_ffi;
/// FFI for symbolic solid state physics.
pub mod symbolic_solid_state_physics_ffi;
/// FFI for symbolic solvers.
pub mod symbolic_solve_ffi;
/// FFI for symbolic special functions.
pub mod symbolic_special_ffi;
/// FFI for symbolic special functions.
pub mod symbolic_special_functions_ffi;
/// FFI for symbolic statistics.
pub mod symbolic_stats_ffi;
/// FFI for symbolic statistical inference.
pub mod symbolic_stats_inference_ffi;
/// FFI for symbolic statistical information theory.
pub mod symbolic_stats_information_theory_ffi;
/// FFI for symbolic statistical probability.
pub mod symbolic_stats_probability_ffi;
/// FFI for symbolic statistical regression.
pub mod symbolic_stats_regression_ffi;
/// FFI for symbolic tensor operations.
pub mod symbolic_tensor_ffi;
/// FFI for symbolic thermodynamics.
pub mod symbolic_thermodynamics_ffi;
/// FFI for symbolic topology.
pub mod symbolic_topology_ffi;
/// FFI for symbolic transforms.
pub mod symbolic_transforms_ffi;
/// FFI for symbolic unit unification.
pub mod symbolic_unit_unification_ffi;
/// FFI for symbolic vector calculus.
pub mod symbolic_vector_calculus_ffi;
/// FFI for symbolic vector operations.
pub mod symbolic_vector_ffi;
