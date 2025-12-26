//! Bincode-based FFI API for constant module.
//!
//! This provides binary serialization for high-performance interop.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::constant_ffi::json::BuildInfo;

/// Returns all build information as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_build_info_bincode() -> BincodeBuffer {

    let info = BuildInfo {
        build_date: crate::constant::get_build_date().to_string(),
        commit_sha: crate::constant::get_commit_sha().to_string(),
        rustc_version: crate::constant::get_rustc_version().to_string(),
        cargo_target_triple: crate::constant::get_cargo_target_triple().to_string(),
        system_info: crate::constant::get_system_info().to_string(),
    };

    match bincode_next::serde::encode_to_vec(
        &info,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Returns the build date as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_build_date_bincode() -> BincodeBuffer {

    let date = crate::constant::get_build_date();

    match bincode_next::serde::encode_to_vec(
        &date,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Returns the commit SHA as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha_bincode() -> BincodeBuffer {

    let sha = crate::constant::get_commit_sha();

    match bincode_next::serde::encode_to_vec(
        &sha,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}
