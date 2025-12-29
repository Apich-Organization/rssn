//! # FFI API for Numerical Elementary Operations
//!
//! Provides bindings for evaluating symbolic expressions and direct numerical functions.
/// bincode-based FFI bindings for numerical elementary operations.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical elementary using opaque handles.
pub mod handle;
/// JSON-based FFI bindings for numerical elementary using serialized data.
pub mod json;
