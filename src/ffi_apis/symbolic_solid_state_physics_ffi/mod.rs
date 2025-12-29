//! FFI API for the Solid-State Physics module.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for Solid-State Physics using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for Solid-State Physics using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for Solid-State Physics using serialized `Expr` values.
pub mod json;
