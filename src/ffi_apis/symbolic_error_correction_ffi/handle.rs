//! Handle-based FFI API for error-correcting codes.
//!
//! This module provides C-compatible FFI functions for Hamming codes and Reed-Solomon codes,
//! including encoding and decoding operations with error correction capabilities.

use crate::symbolic::error_correction::{
    hamming_encode, hamming_decode, rs_encode, rs_decode,
};

/// Encodes 4 data bits into a 7-bit Hamming(7,4) codeword.
///
/// # Safety
/// Caller must ensure `data` points to 4 bytes and `out` points to 7 bytes of allocated memory.
#[no_mangle]
pub unsafe extern "C" fn rssn_hamming_encode(data: *const u8, out: *mut u8) -> i32 {
    if data.is_null() || out.is_null() {
        return -1;
    }
    let slice = std::slice::from_raw_parts(data, 4);
    match hamming_encode(slice) {
        Some(codeword) => {
            for (i, &b) in codeword.iter().enumerate() {
                *out.add(i) = b;
            }
            0
        }
        None => -1,
    }
}

/// Decodes a 7-bit Hamming(7,4) codeword, correcting single-bit errors.
///
/// # Safety
/// Caller must ensure `codeword` points to 7 bytes and `data_out` points to 4 bytes.
/// `error_pos` will receive the 1-based error position or 0 if no error.
#[no_mangle]
pub unsafe extern "C" fn rssn_hamming_decode(
    codeword: *const u8,
    data_out: *mut u8,
    error_pos: *mut u8
) -> i32 {
    if codeword.is_null() || data_out.is_null() || error_pos.is_null() {
        return -1;
    }
    let slice = std::slice::from_raw_parts(codeword, 7);
    match hamming_decode(slice) {
        Ok((data, pos)) => {
            for (i, &b) in data.iter().enumerate() {
                *data_out.add(i) = b;
            }
            *error_pos = pos.unwrap_or(0) as u8;
            0
        }
        Err(_) => -1,
    }
}

/// Encodes data using Reed-Solomon code with n_sym error correction symbols.
///
/// # Safety
/// Caller must ensure `data` is valid. Returns allocated memory that must be freed.
#[no_mangle]
pub unsafe extern "C" fn rssn_rs_encode(
    data: *const u8,
    data_len: usize,
    n_sym: usize,
    out_len: *mut usize
) -> *mut u8 {
    if data.is_null() || out_len.is_null() {
        return std::ptr::null_mut();
    }
    let slice = std::slice::from_raw_parts(data, data_len);
    match rs_encode(slice, n_sym) {
        Ok(codeword) => {
            *out_len = codeword.len();
            let boxed = codeword.into_boxed_slice();
            Box::into_raw(boxed) as *mut u8
        }
        Err(_) => std::ptr::null_mut(),
    }
}

/// Decodes a Reed-Solomon codeword, correcting errors if possible.
///
/// # Safety
/// Caller must ensure `codeword` is valid. Returns allocated memory that must be freed.
#[no_mangle]
pub unsafe extern "C" fn rssn_rs_decode(
    codeword: *const u8,
    codeword_len: usize,
    n_sym: usize,
    out_len: *mut usize
) -> *mut u8 {
    if codeword.is_null() || out_len.is_null() {
        return std::ptr::null_mut();
    }
    let slice = std::slice::from_raw_parts(codeword, codeword_len);
    match rs_decode(slice, n_sym) {
        Ok(data) => {
            *out_len = data.len();
            let boxed = data.into_boxed_slice();
            Box::into_raw(boxed) as *mut u8
        }
        Err(_) => std::ptr::null_mut(),
    }
}

/// Frees memory allocated by rs_encode or rs_decode.
///
/// # Safety
/// Caller must ensure `ptr` was returned by rssn_rs_encode or rssn_rs_decode.
#[no_mangle]
pub unsafe extern "C" fn rssn_rs_free(ptr: *mut u8, len: usize) {
    if !ptr.is_null() && len > 0 {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, len));
    }
}
