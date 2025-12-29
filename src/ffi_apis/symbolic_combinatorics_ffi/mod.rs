//! FFI bindings for symbolic combinatorics operations.
//!
//! This module exposes Handle-based, JSON, and Bincode APIs for permutations,
//! combinations, Catalan numbers, Stirling numbers, and Bell numbers.
/// Bincode-based FFI bindings for combinatorial functions using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for combinatorial functions using opaque `Expr` pointers.
pub mod handle;
/// JSON-based FFI bindings for combinatorial functions using serialized `Expr` values.
pub mod json;
