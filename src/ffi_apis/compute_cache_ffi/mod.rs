//! FFI APIs for the compute cache module.
//!
//! This module provides three different FFI API versions:
//! - **Handle-based**: Traditional C-style functions
//! - **JSON-based**: String serialization for easy language interop
//! - **Bincode-based**: Binary serialization for high performance

pub mod handle;
pub mod json;
pub mod bincode_api;

// Re-export all functions for convenience
pub use handle::*;
pub use json::*;
pub use bincode_api::*;
