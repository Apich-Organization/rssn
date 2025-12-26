//! JSON-based FFI API for constant module.
//!
//! This provides string-based serialization for easy language interop.

use crate::ffi_apis::common::to_c_string;
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

/// Build information structure for JSON serialization.
#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct BuildInfo {
    pub build_date: String,
    pub commit_sha: String,
    pub rustc_version: String,
    pub cargo_target_triple: String,
    pub system_info: String,
}

/// Returns all build information as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_info_json() -> *mut c_char {

    let info = BuildInfo {
        build_date: crate::constant::get_build_date().to_string(),
        commit_sha: crate::constant::get_commit_sha().to_string(),
        rustc_version: crate::constant::get_rustc_version().to_string(),
        cargo_target_triple: crate::constant::get_cargo_target_triple().to_string(),
        system_info: crate::constant::get_system_info().to_string(),
    };

    match serde_json::to_string(&info) {
        | Ok(json) => to_c_string(json),
        | Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the build date as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_date_json() -> *mut c_char {

    let date = crate::constant::get_build_date();

    match serde_json::to_string(&date) {
        | Ok(json) => to_c_string(json),
        | Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the commit SHA as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha_json() -> *mut c_char {

    let sha = crate::constant::get_commit_sha();

    match serde_json::to_string(&sha) {
        | Ok(json) => to_c_string(json),
        | Err(_) => std::ptr::null_mut(),
    }
}
