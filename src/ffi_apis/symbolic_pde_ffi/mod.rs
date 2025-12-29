//! FFI APIs for the symbolic PDE module.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI bindings for symbolic PDE operations.
pub mod bincode_api;
/// Handle-based FFI bindings for PDE using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for PDE using serialized data.
pub mod json;

pub use bincode_api::*;
pub use handle::*;
pub use json::*;
