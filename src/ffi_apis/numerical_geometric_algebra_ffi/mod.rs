//! FFI API for numerical geometric algebra operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for numerical geometric algebra operations.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical geometric algebra using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical geometric algebra using serialized data.
pub mod json;
