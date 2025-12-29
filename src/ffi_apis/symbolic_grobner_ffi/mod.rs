//! FFI bindings for symbolic Gröbner basis computations.
//!
//! This module exposes foreign function interface (FFI) bindings for multivariate
//! polynomial ideal operations and Gröbner basis algorithms.
/// Bincode-based FFI bindings for symbolic Gröbner basis operations.
pub mod bincode_api;
/// Handle-based FFI bindings for Gröbner bases using opaque polynomial ideal handles.
pub mod handle;
/// JSON-based FFI bindings for Gröbner basis computations using serialized polynomial data.
pub mod json;
