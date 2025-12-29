//! FFI bindings for symbolic integration operations.
//!
//! This module exposes foreign function interface (FFI) bindings for indefinite and
//! definite integration, integration by parts, and related calculus operations.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// Bincode-based FFI bindings for symbolic integration operations.
pub mod bincode_api;
/// Handle-based FFI bindings for integration using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for integration using serialized expressions.
pub mod json;
