//! FFI API for symbolic integral equation operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// Bincode-based FFI bindings for symbolic integral equation operations.
pub mod bincode_api;
/// Handle-based FFI bindings for integral equations using opaque handles to equations and solutions.
pub mod handle;
/// JSON-based FFI bindings for integral equations using serialized symbolic data.
pub mod json;

pub use bincode_api::*;
pub use handle::*;
pub use json::*;
