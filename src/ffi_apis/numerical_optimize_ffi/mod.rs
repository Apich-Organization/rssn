//! FFI API for numerical optimization functions.
//!
//! This module provides three different FFI interfaces for numerical optimization:
//! - Bincode-based serialization API
//! - Handle-based API for managing optimizer state
//! - JSON-based serialization API

/// Bincode serialization API for numerical optimization.
pub mod bincode_api;
/// Handle-based API for managing numerical optimizer instances.
pub mod handle;
/// JSON serialization API for numerical optimization.
pub mod json;
