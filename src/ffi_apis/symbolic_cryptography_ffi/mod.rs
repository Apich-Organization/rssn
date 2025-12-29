//! FFI API for cryptographic operations.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported operations include:
//! - Elliptic curve primitives: Point validation, arithmetic (add, double, negate, multiply).
//! - Key Management: ECDH key pair generation.
//! - Protocol Operations: ECDH shared secret derivation.
//! - Digital Signatures: ECDSA signing and verification.
//! - Serialization: Point compression and decompression.
/// bincode-based FFI bindings for cryptographic operations using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for cryptographic operations using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for cryptographic operations using serialized `Expr` values.
pub mod json;
