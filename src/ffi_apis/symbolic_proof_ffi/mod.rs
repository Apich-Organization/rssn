//! FFI bindings for symbolic proof operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI bindings for symbolic proof operations.
pub mod bincode_api;
/// Handle-based FFI bindings for proofs using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for proofs using serialized data.
pub mod json;
