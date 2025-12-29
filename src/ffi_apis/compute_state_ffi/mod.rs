//! FFI APIs for the compute state module.
//!
//! This module provides three different FFI API versions:
//! - **Handle-based**: Traditional C-style functions with opaque pointers
//! - **JSON-based**: String serialization for easy language interop
//! - **Bincode-based**: Binary serialization for high performance

/// bincode-based FFI bindings for compute state operations.
pub mod bincode_api;
/// Handle-based FFI bindings for compute state using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for compute state using serialized data.
pub mod json;

// Re-export all functions for convenience
pub use bincode_api::*;
pub use handle::*;
pub use json::*;
