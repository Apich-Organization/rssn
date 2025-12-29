//! FFI bindings for symbolic quantum mechanics operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// Bincode-based FFI bindings for symbolic quantum mechanics operations.
pub mod bincode_api;
/// Handle-based FFI bindings for quantum mechanics using opaque handles.
pub mod handle;
/// JSON-based FFI API.
pub mod json;
