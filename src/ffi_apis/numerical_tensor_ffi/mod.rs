//! FFI API for numerical tensor operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for numerical tensor operations.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical tensor using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical tensor using serialized data.
pub mod json;
