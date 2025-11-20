//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(clippy::no_mangle_with_rust_abi)]

pub mod common;
pub mod ffi_api;
pub mod constant_ffi;
pub mod compute_cache_ffi;
pub mod compute_state_ffi;
