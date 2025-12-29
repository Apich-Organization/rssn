//! FFI bindings for symbolic radical operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI API.
pub mod bincode_api;
/// Handle-based FFI API.
pub mod handle;
/// JSON-based FFI API.
pub mod json;
