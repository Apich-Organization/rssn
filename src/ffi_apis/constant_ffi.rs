use crate::constant;
use std::ffi::CString;
use std::os::raw::c_char;

/// Returns the build date as a C string.
/// The caller is responsible for freeing the memory.
#[no_mangle]
pub extern "C" fn rssn_get_build_date() -> *mut c_char {
    let s = constant::get_build_date();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Returns the commit SHA as a C string.
/// The caller is responsible for freeing the memory.
#[no_mangle]
pub extern "C" fn rssn_get_commit_sha() -> *mut c_char {
    let s = constant::get_commit_sha();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Returns the rustc version as a C string.
/// The caller is responsible for freeing the memory.
#[no_mangle]
pub extern "C" fn rssn_get_rustc_version() -> *mut c_char {
    let s = constant::get_rustc_version();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Returns the cargo target triple as a C string.
/// The caller is responsible for freeing the memory.
#[no_mangle]
pub extern "C" fn rssn_get_cargo_target_triple() -> *mut c_char {
    let s = constant::get_cargo_target_triple();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Returns the system info as a C string.
/// The caller is responsible for freeing the memory.
#[no_mangle]
pub extern "C" fn rssn_get_system_info() -> *mut c_char {
    let s = constant::get_system_info();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Frees a string returned by any of the rssn_get_* functions.
#[no_mangle]
pub extern "C" fn rssn_free_string(s: *mut c_char) {
    if s.is_null() {
        return;
    }
    unsafe {
        let _ = CString::from_raw(s);
    }
}
