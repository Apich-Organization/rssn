//! FFI API for symbolic handle manager.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (C-style functions)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)

pub mod handle;
pub mod json;
pub mod bincode_api;

pub use handle::*;
pub use json::*;
pub use bincode_api::*;
