//! FFI API for symbolic CAS foundations (expansion, factorization, normalization).
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (opaque pointers)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

pub mod handle;
pub mod json;
pub mod bincode_api;

pub use handle::*;
pub use json::*;
pub use bincode_api::*;
