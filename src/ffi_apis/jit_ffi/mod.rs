//! JIT FFI module.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Handle-based FFI bindings for JIT using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for JIT using serialized data.
pub mod json;
