//! FFI bindings for symbolic finite field operations.
//!
//! This module exposes foreign function interface (FFI) bindings for arithmetic
//! and algebraic structures over finite fields, enabling them to be used from
//! other languages.
pub mod bincode_api;
/// Handle-based FFI bindings for finite fields using opaque handles to field elements and structures.
pub mod handle;
/// JSON-based FFI bindings for finite field operations using serialized values.
pub mod json;
