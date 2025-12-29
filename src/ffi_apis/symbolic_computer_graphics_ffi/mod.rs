//! FFI API for computer graphics operations.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Opaque pointer-based API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! Supported operations include 2D/3D transformations (translation, rotation, scaling),
//! perspective and orthographic projections, Bezier curves, B-splines, and polygon meshes.
/// bincode-based FFI bindings for computer graphics operations using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for computer graphics operations using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for computer graphics operations using serialized `Expr` values.
pub mod json;
