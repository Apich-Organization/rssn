//! FFI API for symbolic handle manager.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

/// The bincode based FFI API for symbolic handles.
pub mod bincode_api;
/// The handle based FFI API for symbolic handles.
pub mod handle;
/// The JSON based FFI API for symbolic handles.
pub mod json;

pub use bincode_api::*;
pub use handle::*;
pub use json::*;
