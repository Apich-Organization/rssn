//! JSON-based FFI API for error-correcting codes.
//!
//! This module provides JSON string-based FFI functions for Hamming codes and Reed-Solomon codes,
//! enabling language-agnostic integration for error correction algorithms.

use crate::symbolic::error_correction::{
    hamming_encode, hamming_decode, rs_encode, rs_decode,
};
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

/// Encodes 4 data bits into a 7-bit Hamming(7,4) codeword via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_hamming_encode(data_json: *const c_char) -> *mut c_char {
    let data: Option<Vec<u8>> = from_json_string(data_json);
    if let Some(d) = data {
        match hamming_encode(&d) {
            Some(codeword) => to_json_string(&codeword),
            None => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Decodes a 7-bit Hamming(7,4) codeword via JSON interface.
/// Returns JSON object with "data" and "error_pos" fields.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_hamming_decode(codeword_json: *const c_char) -> *mut c_char {
    let codeword: Option<Vec<u8>> = from_json_string(codeword_json);
    if let Some(c) = codeword {
        match hamming_decode(&c) {
            Ok((data, error_pos)) => {
                let result = serde_json::json!({
                    "data": data,
                    "error_pos": error_pos
                });
                to_json_string(&result)
            }
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Encodes data using Reed-Solomon code via JSON interface.
/// Input: {"data": [bytes], "n_sym": number}
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rs_encode(
    data_json: *const c_char,
    n_sym_json: *const c_char
) -> *mut c_char {
    let data: Option<Vec<u8>> = from_json_string(data_json);
    let n_sym: Option<usize> = from_json_string(n_sym_json);
    if let (Some(d), Some(n)) = (data, n_sym) {
        match rs_encode(&d, n) {
            Ok(codeword) => to_json_string(&codeword),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

/// Decodes a Reed-Solomon codeword via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rs_decode(
    codeword_json: *const c_char,
    n_sym_json: *const c_char
) -> *mut c_char {
    let codeword: Option<Vec<u8>> = from_json_string(codeword_json);
    let n_sym: Option<usize> = from_json_string(n_sym_json);
    if let (Some(c), Some(n)) = (codeword, n_sym) {
        match rs_decode(&c, n) {
            Ok(data) => to_json_string(&data),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}
