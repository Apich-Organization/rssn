//! FFI API for error-correcting codes.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported codes include Hamming(7,4) single-error correcting code and
//! Reed-Solomon codes with configurable error correction capability.

pub mod handle;
pub mod json;
pub mod bincode_api;
