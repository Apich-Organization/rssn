//! FFI bindings for symbolic discrete group operations.
//!
//! This module exposes foreign function interface (FFI) bindings for constructing
//! finite discrete groups such as cyclic, dihedral, symmetric, and Klein four groups,
//! enabling them to be used from other languages.
pub mod bincode_api;
/// Handle-based FFI bindings for discrete groups using opaque `Group` handles.
pub mod handle;
/// JSON-based FFI bindings for discrete groups using serialized group values.
pub mod json;
