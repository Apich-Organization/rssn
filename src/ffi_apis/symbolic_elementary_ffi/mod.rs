//! FFI APIs for the symbolic elementary module.
//!
//! This module provides three different FFI API versions:
//! - **Handle-based**: Traditional C-style functions with opaque pointers
//! - **JSON-based**: String serialization for easy language interop
//! - **Bincode-based**: Binary serialization for high performance

/// bincode-based FFI bindings for symbolic elementary operations using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for symbolic elementary operations using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic elementary operations using serialized `Expr` values.
pub mod json;

// Re-export all functions for convenience
pub use bincode_api::*;
pub use handle::*;
pub use json::*;
