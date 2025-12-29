//! FFI bindings for symbolic graph operations.
//!
//! This module exposes foreign function interface (FFI) bindings for operations on
//! symbolic graphs, such as union, intersection, complements, and subgraph
//! extraction.
/// Bincode-based FFI bindings for symbolic graph operations.
pub mod bincode_api;
/// Handle-based FFI bindings for graph operations using opaque graph handles.
pub mod handle;
/// JSON-based FFI bindings for graph operations using serialized graph representations.
pub mod json;
