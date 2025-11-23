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
pub mod ffi_api;
pub mod constant_ffi;
pub mod compute_cache_ffi;
pub mod compute_state_ffi;
pub mod symbolic_elementary_ffi;
pub mod symbolic_calculus_ffi;
pub mod symbolic_simplify_dag_ffi;
pub mod symbolic_simplify_ffi;
pub mod symbolic_polynomial_ffi;
pub mod symbolic_matrix_ffi;
pub mod symbolic_numeric_ffi;
pub mod symbolic_vector_ffi;
pub mod symbolic_tensor_ffi;
pub mod symbolic_coordinates_ffi;
pub mod symbolic_unit_unification_ffi;
pub mod symbolic_solve_ffi;
