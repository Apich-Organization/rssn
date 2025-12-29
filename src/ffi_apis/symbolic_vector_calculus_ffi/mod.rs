//! FFI API for symbolic vector calculus functions.
//!
//! This module provides three types of FFI interfaces:
//! - Handle-based API (opaque pointers)
//! - JSON-based API (string serialization)
//! - Bincode-based API (binary serialization)
/// bincode-based FFI bindings for symbolic vector calculus using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for symbolic vector calculus using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for symbolic vector calculus using serialized `Expr` values.
pub mod json;

pub use bincode_api::*;
pub use handle::*;
pub use json::*;
