//! FFI bindings for symbolic coordinate operations.
//!
//! This module provides foreign function interface (FFI) bindings for symbolic
//! coordinate-related operations, allowing them to be called from other
//! languages. It includes APIs for bincode serialization, handle-based
//! operations, and JSON serialization.

pub mod bincode_api;
pub mod handle;
pub mod json;
