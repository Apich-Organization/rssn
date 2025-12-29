//! FFI bindings for symbolic graph algorithm operations.
//!
//! This module exposes foreign function interface (FFI) bindings for graph algorithms
//! such as shortest paths, connectivity, flows, and centrality measures on symbolic graphs.
pub mod bincode_api;
/// Handle-based FFI bindings for graph algorithms using opaque graph and result handles.
pub mod handle;
/// JSON-based FFI bindings for graph algorithms using serialized graph data.
pub mod json;
