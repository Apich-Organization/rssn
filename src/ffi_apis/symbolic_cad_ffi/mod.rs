//! # CAD FFI Module
//!
//! This module provides FFI bindings for the Cylindrical Algebraic Decomposition (CAD)
//! functionality. It exposes Handle-based, JSON, and Bincode APIs for computing
//! and interacting with CAD cells.

pub mod bincode_api;
pub mod handle;
pub mod json;
