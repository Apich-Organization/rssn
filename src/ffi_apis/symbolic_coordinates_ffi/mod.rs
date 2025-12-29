//! FFI bindings for symbolic coordinate operations.
//!
//! This module provides foreign function interface (FFI) bindings for symbolic
//! coordinate-related operations, allowing them to be called from other
//! languages. It includes APIs for bincode serialization, handle-based
//! operations, and JSON serialization.

/// Bincode-based FFI bindings for symbolic coordinate operations.
pub mod bincode_api;
/// Handle-based FFI bindings operating on coordinate `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic coordinate operations.
pub mod json;
