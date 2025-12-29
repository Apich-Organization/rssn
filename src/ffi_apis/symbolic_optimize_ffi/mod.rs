//! FFI bindings for symbolic optimization operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI bindings for symbolic optimization operations.
pub mod bincode_api;
/// Handle-based FFI bindings for optimization using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for optimization using serialized data.
pub mod json;
