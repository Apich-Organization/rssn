//! FFI API for the Calculus of Variations module.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// bincode-based FFI bindings for Calculus of Variations.
pub mod bincode_api;
/// Handle-based FFI bindings for Calculus of Variations using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for Calculus of Variations using serialized data.
pub mod json;
