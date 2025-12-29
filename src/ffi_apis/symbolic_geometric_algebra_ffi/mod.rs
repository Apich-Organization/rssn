//! FFI bindings for symbolic geometric algebra operations.
//!
//! This module exposes foreign function interface (FFI) bindings for multivectors,
//! geometric products, and other Clifford / geometric algebra constructions.
/// Bincode-based FFI bindings for symbolic geometric algebra operations.
pub mod bincode_api;
/// Handle-based FFI bindings for geometric algebra using opaque multivector handles.
pub mod handle;
/// JSON-based FFI bindings for geometric algebra using serialized symbolic values.
pub mod json;
