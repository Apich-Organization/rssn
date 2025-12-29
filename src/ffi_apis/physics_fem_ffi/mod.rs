//! FFI API for physics FEM (Finite Element Method) functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for physics FEM functions.
pub mod bincode_api;
/// Handle-based FFI bindings for physics FEM using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for physics FEM using serialized data.
pub mod json;
