//! FFI bindings for symbolic graph data structures and basic operations.
//!
//! This module exposes foreign function interface (FFI) bindings for constructing
//! and manipulating symbolic graphs, which are then consumed by higher-level
//! graph algorithms.

/// Bincode-based FFI bindings for symbolic graph operations.
pub mod bincode_api;
/// Handle-based FFI bindings for graph structures using opaque graph handles.
pub mod handle;
/// JSON-based FFI bindings for graph structures using serialized graph representations.
pub mod json;
