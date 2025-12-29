//! FFI bindings for symbolic graph isomorphism and coloring operations.
//!
//! This module exposes foreign function interface (FFI) bindings for testing graph
//! isomorphism, computing canonical labelings, and performing graph coloring
//! algorithms on symbolic graphs.

/// Bincode-based FFI bindings for symbolic graph isomorphism and coloring operations.
pub mod bincode_api;
/// Handle-based FFI bindings for graph isomorphism and coloring using opaque graph handles.
pub mod handle;
/// JSON-based FFI bindings for graph isomorphism and coloring using serialized graph data.
pub mod json;
