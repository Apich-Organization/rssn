//! FFI API for numerical coordinate transformations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for numerical coordinate transformations.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical coordinate transformations using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical coordinate transformations using serialized data.
pub mod json;
