//! FFI bindings for symbolic integration operations.
//!
//! This module exposes foreign function interface (FFI) bindings for indefinite and
//! definite integration, integration by parts, and related calculus operations.
pub mod bincode_api;
/// Handle-based FFI bindings for integration using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for integration using serialized expressions.
pub mod json;
