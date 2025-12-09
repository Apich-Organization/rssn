//! FFI API for cryptographic operations.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported operations include elliptic curve point addition, scalar multiplication,
//! ECDH key pair generation, and shared secret derivation.

pub mod handle;
pub mod json;
pub mod bincode_api;
