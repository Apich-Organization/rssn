//! FFI bindings for symbolic Lie group and Lie algebra operations.
//!
//! This module exposes foreign function interface (FFI) bindings for matrix Lie groups,
//! Lie algebras, exponentials, logarithms, and related structures.
pub mod bincode_api;
/// Handle-based FFI bindings for Lie groups using opaque group and algebra handles.
pub mod handle;
/// JSON-based FFI bindings for Lie groups and algebras using serialized symbolic data.
pub mod json;
