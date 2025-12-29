//! FFI API for finite field (Galois field) operations and error correction helpers.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported operations include GF(2^8) arithmetic (add, mul, inv, div),
//! polynomial operations over GF(2^8), and general finite field arithmetic.

/// bincode-based FFI bindings for error correction helpers using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for error correction helpers using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for error correction helpers using serialized `Expr` values.
pub mod json;
