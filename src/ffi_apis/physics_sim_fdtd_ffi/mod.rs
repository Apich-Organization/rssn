//! FFI API for physics sim FDTD electrodynamics functions.

//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for physics sim FDTD electrodynamics functions.
pub mod bincode_api;
/// Handle-based FFI bindings for physics sim FDTD electrodynamics using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for physics sim FDTD electrodynamics using serialized data.
pub mod json;
