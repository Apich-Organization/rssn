//! FFI bindings for symbolic functional analysis operations.
//!
//! This module exposes foreign function interface (FFI) bindings for operators on
//! function spaces, norms, inner products, and spectral properties, enabling
//! functional-analytic constructs to be used from other languages.

/// Bincode-based FFI bindings for symbolic functional analysis operations.
pub mod bincode_api;
/// Handle-based FFI bindings for functional analysis using opaque handles to operators and spaces.
pub mod handle;
/// JSON-based FFI bindings for functional analysis using serialized symbolic data.
pub mod json;
