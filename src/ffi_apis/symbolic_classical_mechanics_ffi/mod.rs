//! FFI API for the Classical Mechanics module.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for Classical Mechanics.
pub mod bincode_api;
/// Handle-based FFI bindings for Classical Mechanics using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for Classical Mechanics using serialized data.
pub mod json;
