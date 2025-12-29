//! FFI bindings for symbolic stats operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for symbolic stats operations using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for symbolic stats operations using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic stats operations using serialized `Expr` values.
pub mod json;
