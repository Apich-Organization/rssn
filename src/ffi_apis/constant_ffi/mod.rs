//! FFI APIs for the constant module.
//!
//! This module provides three different FFI API versions:
//! - **Handle-based**: Traditional C-style functions returning strings
//! - **JSON-based**: String serialization for easy language interop
//! - **Bincode-based**: Binary serialization for high performance
//!
//! # Examples
//!
//! ## Handle-based API (C/C++)
//! ```c
//! char* date = rssn_get_build_date();
//! printf("Build date: %s\n", date);
//! rssn_free_string(date);
//! ```
//!
//! ## JSON-based API (Python, JavaScript, etc.)
//! ```c
//! char* json = rssn_get_build_info_json();
//! // Parse JSON: {"build_date": "...", "commit_sha": "...", ...}
//! rssn_free_string(json);
//! ```
//!
//! ## Bincode-based API (High performance)
//! ```c
//! BincodeBuffer buf = rssn_get_build_info_bincode();
//! // Deserialize binary data
//! rssn_free_bincode_buffer(buf);
//! ```

pub mod handle;
pub mod json;
pub mod bincode_api;

// Re-export all functions for convenience
pub use handle::*;
pub use json::*;
pub use bincode_api::*;
