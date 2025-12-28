//! FFI bindings for symbolic convergence operations.
//!
//! This module provides foreign function interface (FFI) bindings for symbolic
//! convergence-related operations, allowing them to be called from other
//! languages. It includes APIs for bincode serialization, handle-based
//! operations, and JSON serialization.
pub mod bincode_api;
pub mod handle;
pub mod json;
