//! FFI API for symbolic integral transforms.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported transforms include Fourier, Laplace, and Z-transforms with their
//! inverse operations and various properties (time shift, frequency shift, 
//! scaling, differentiation, and convolution theorems).

pub mod handle;
pub mod json;
pub mod bincode_api;
