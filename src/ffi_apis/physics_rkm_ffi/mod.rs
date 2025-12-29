//! FFI API for physics RKM (Runge-Kutta Methods) functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for physics RKM functions.
pub mod bincode_api;
/// Handle-based FFI bindings for physics RKM using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for physics RKM using serialized data.
pub mod json;
