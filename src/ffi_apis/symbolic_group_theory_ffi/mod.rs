//! FFI bindings for symbolic group theory operations.
//!
//! This module exposes foreign function interface (FFI) bindings for abstract groups,
//! homomorphisms, and group-theoretic constructions.
/// Bincode-based FFI bindings for symbolic group theory operations.
pub mod bincode_api;
/// Handle-based FFI bindings for group theory using opaque `Group` and related handles.
pub mod handle;
/// JSON-based FFI bindings for group theory using serialized group-theoretic data.
pub mod json;
