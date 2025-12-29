//! FFI bindings for symbolic multi-valued operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI bindings for symbolic multi-valued operations.
pub mod bincode_api;
/// Handle-based FFI bindings for multi-valued operations using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for multi-valued operations using serialized data.
pub mod json;
