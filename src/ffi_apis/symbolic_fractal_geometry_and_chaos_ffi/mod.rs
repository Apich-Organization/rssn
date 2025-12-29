//! FFI bindings for symbolic fractal geometry and chaos operations.
//!
//! This module exposes foreign function interface (FFI) bindings for iterated
//! maps, fractal sets, and chaotic dynamical systems, allowing them to be
//! explored from other languages.

/// Bincode-based FFI bindings for symbolic fractal geometry and chaos operations.
pub mod bincode_api;
/// Handle-based FFI bindings for fractal and chaos computations using opaque symbolic handles.
pub mod handle;
/// JSON-based FFI bindings for fractal geometry and chaos using serialized data.
pub mod json;
