//! FFI API for numerical error correction functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for numerical error correction functions.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical error correction using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical error correction using serialized data.
pub mod json;
