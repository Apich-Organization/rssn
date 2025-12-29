//! FFI API for numerical sparse matrix operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for numerical sparse matrix operations.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical sparse matrix using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical sparse matrix using serialized data.
pub mod json;
