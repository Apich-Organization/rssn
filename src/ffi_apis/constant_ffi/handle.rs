//! Handle-based FFI API for constant module.
//!
//! This provides traditional C-style opaque pointer functions.

use std::ffi::CString;
use std::os::raw::c_char;

/// Returns the build date as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_date() -> *mut c_char {

    let date = crate::constant::get_build_date();

    match CString::new(date) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the commit SHA as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha() -> *mut c_char {

    let sha = crate::constant::get_commit_sha();

    match CString::new(sha) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the rustc version as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_rustc_version() -> *mut c_char {

    let version = crate::constant::get_rustc_version();

    match CString::new(version) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the cargo target triple as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_cargo_target_triple() -> *mut c_char {

    let triple = crate::constant::get_cargo_target_triple();

    match CString::new(triple) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Returns the system info as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_system_info() -> *mut c_char {

    let info = crate::constant::get_system_info();

    match CString::new(info) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}
