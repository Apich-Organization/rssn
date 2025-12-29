//! # CAD FFI Module
//!
//! This module provides FFI bindings for the Cylindrical Algebraic Decomposition (CAD)
//! functionality. It exposes Handle-based, JSON, and Bincode APIs for computing
//! and interacting with CAD cells.

/// Bincode-based FFI bindings for Cylindrical Algebraic Decomposition (CAD) operations.
pub mod bincode_api;
/// Handle-based FFI bindings exposing opaque CAD cell handles.
pub mod handle;
/// JSON-based FFI bindings for interoperable CAD computations.
pub mod json;
