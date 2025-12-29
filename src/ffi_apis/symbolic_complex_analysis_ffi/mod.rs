//! FFI for complex analysis.
//!
//! This module exposes FFI bindings for symbolic complex analysis operations,
//! including analytic continuation, singularity classification, and contour
//! integration in the complex plane.

pub mod bincode_api;
/// Handle-based FFI bindings for complex analysis using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for complex analysis using serialized `Expr` values.
pub mod json;
