//! FFI API for physics MTM (Multigrid Transfer Method) functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for physics MTM functions.
pub mod bincode_api;
/// Handle-based FFI bindings for physics MTM using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for physics MTM using serialized data.
pub mod json;
