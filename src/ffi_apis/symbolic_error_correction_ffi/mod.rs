//! FFI API for error-correcting codes.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! ## Supported Codes
//! 
//! ### Hamming(7,4)
//! - Single-error correcting code
//! - Functions: `encode`, `decode`, `distance`, `weight`, `check`
//!
//! ### Reed-Solomon
//! - Multi-symbol error correction with configurable capability
//! - Functions: `encode`, `decode`, `check`, `error_count`
//!
//! ### CRC-32 (IEEE 802.3)
//! - Error detection via cyclic redundancy check
//! - Functions: `compute`, `verify`, `update`, `finalize`

pub mod handle;
pub mod json;
pub mod bincode_api;

