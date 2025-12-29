//! FFI API for symbolic CAS foundations (expansion, factorization, normalization).
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (opaque pointers)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for symbolic CAS foundations (expansion, factorization, normalization).
pub mod bincode_api;
/// Handle-based FFI bindings for symbolic CAS foundations using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic CAS foundations using serialized data.
pub mod json;

pub use bincode_api::*;
pub use handle::*;
pub use json::*;
