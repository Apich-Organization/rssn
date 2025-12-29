//! FFI API for numerical stats functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for numerical stats functions.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical stats using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical stats using serialized data.
pub mod json;
