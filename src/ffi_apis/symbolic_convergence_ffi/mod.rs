//! FFI bindings for symbolic convergence operations.
//!
//! This module provides foreign function interface (FFI) bindings for symbolic
//! convergence-related operations, allowing them to be called from other
//! languages. It includes APIs for bincode serialization, handle-based
//! operations, and JSON serialization.

/// Bincode-based FFI bindings for symbolic convergence operations.
pub mod bincode_api;
/// Handle-based FFI bindings exposing opaque convergence-related `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic convergence operations.
pub mod json;
