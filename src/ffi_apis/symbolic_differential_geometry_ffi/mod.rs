//! FFI bindings for symbolic differential geometry operations.
//!
//! This module exposes foreign function interface (FFI) bindings for differential
//! geometric constructions such as differential forms, exterior derivatives, and
//! curvature tensors, allowing them to be used from other languages.

/// Bincode-based FFI bindings for symbolic differential_geometry
pub mod bincode_api;
/// Handle-based FFI bindings for differential geometry using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic differential geometry operations.
pub mod json;
