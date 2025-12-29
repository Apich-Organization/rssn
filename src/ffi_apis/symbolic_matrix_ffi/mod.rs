//! FFI bindings for symbolic matrix operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// Bincode-based FFI bindings for symbolic matrix operations.
pub mod bincode_api;
/// Handle-based FFI bindings for matrix operations using opaque `Matrix` handles.
pub mod handle;
/// JSON-based FFI bindings for matrix operations using serialized matrix data.
pub mod json;
