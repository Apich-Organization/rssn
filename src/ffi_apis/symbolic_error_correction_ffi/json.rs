//! JSON-based FFI API for error-correcting codes.
//!
//! This module provides JSON string-based FFI functions for Hamming codes, Reed-Solomon codes,
//! and CRC-32, enabling language-agnostic integration for error correction algorithms.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::error_correction::crc32_compute;
use crate::symbolic::error_correction::crc32_finalize;
use crate::symbolic::error_correction::crc32_update;
use crate::symbolic::error_correction::crc32_verify;
use crate::symbolic::error_correction::hamming_check;
use crate::symbolic::error_correction::hamming_decode;
use crate::symbolic::error_correction::hamming_distance;
use crate::symbolic::error_correction::hamming_encode;
use crate::symbolic::error_correction::hamming_weight;
use crate::symbolic::error_correction::rs_check;
use crate::symbolic::error_correction::rs_decode;
use crate::symbolic::error_correction::rs_encode;
use crate::symbolic::error_correction::rs_error_count;

/// Encodes 4 data bits into a 7-bit Hamming(7,4) codeword via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hamming_encode(
    data_json : *const c_char
) -> *mut c_char {

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    if let Some(d) = data {

        match hamming_encode(&d) {
            | Some(codeword) => {
                to_json_string(
                    &codeword,
                )
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Decodes a 7-bit Hamming(7,4) codeword via JSON interface.
/// Returns JSON object with "data" and "error_pos" fields.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hamming_decode(
    codeword_json : *const c_char
) -> *mut c_char {

    let codeword : Option<Vec<u8>> =
        from_json_string(codeword_json);

    if let Some(c) = codeword {

        match hamming_decode(&c) {
            | Ok((data, error_pos)) => {

                let result = serde_json::json!({
                    "data": data,
                    "error_pos": error_pos
                });

                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Encodes data using Reed-Solomon code via JSON interface.
/// Input: {"data": [bytes], "n_sym": number}
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rs_encode(
    data_json : *const c_char,
    n_sym_json : *const c_char,
) -> *mut c_char {

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    let n_sym : Option<usize> =
        from_json_string(n_sym_json);

    if let (Some(d), Some(n)) =
        (data, n_sym)
    {

        match rs_encode(&d, n) {
            | Ok(codeword) => {
                to_json_string(
                    &codeword,
                )
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Decodes a Reed-Solomon codeword via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rs_decode(
    codeword_json : *const c_char,
    n_sym_json : *const c_char,
) -> *mut c_char {

    let codeword : Option<Vec<u8>> =
        from_json_string(codeword_json);

    let n_sym : Option<usize> =
        from_json_string(n_sym_json);

    if let (Some(c), Some(n)) =
        (codeword, n_sym)
    {

        match rs_decode(&c, n) {
            | Ok(data) => {
                to_json_string(&data)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Hamming Code Extensions
// ============================================================================

/// Computes Hamming distance between two byte slices via JSON interface.
/// Input: {"a": [bytes], "b": [bytes]}
/// Returns: distance as integer, or null on error
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hamming_distance(
    a_json : *const c_char,
    b_json : *const c_char,
) -> *mut c_char {

    let a : Option<Vec<u8>> =
        from_json_string(a_json);

    let b : Option<Vec<u8>> =
        from_json_string(b_json);

    if let (Some(av), Some(bv)) = (a, b)
    {

        match hamming_distance(&av, &bv)
        {
            | Some(dist) => {
                to_json_string(&dist)
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Hamming weight of a byte slice via JSON interface.
/// Input: [bytes]
/// Returns: weight as integer
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hamming_weight(
    data_json : *const c_char
) -> *mut c_char {

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    if let Some(d) = data {

        let weight = hamming_weight(&d);

        to_json_string(&weight)
    } else {

        std::ptr::null_mut()
    }
}

/// Checks if a Hamming(7,4) codeword is valid via JSON interface.
/// Input: [7 bytes]
/// Returns: boolean
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hamming_check(
    codeword_json : *const c_char
) -> *mut c_char {

    let codeword : Option<Vec<u8>> =
        from_json_string(codeword_json);

    if let Some(c) = codeword {

        let valid = hamming_check(&c);

        to_json_string(&valid)
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// Reed-Solomon Enhancements
// ============================================================================

/// Checks if a Reed-Solomon codeword is valid via JSON interface.
/// Returns: boolean
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rs_check(
    codeword_json : *const c_char,
    n_sym_json : *const c_char,
) -> *mut c_char {

    let codeword : Option<Vec<u8>> =
        from_json_string(codeword_json);

    let n_sym : Option<usize> =
        from_json_string(n_sym_json);

    if let (Some(c), Some(n)) =
        (codeword, n_sym)
    {

        let valid = rs_check(&c, n);

        to_json_string(&valid)
    } else {

        std::ptr::null_mut()
    }
}

/// Estimates error count in a Reed-Solomon codeword via JSON interface.
/// Returns: error count as integer
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rs_error_count(
    codeword_json : *const c_char,
    n_sym_json : *const c_char,
) -> *mut c_char {

    let codeword : Option<Vec<u8>> =
        from_json_string(codeword_json);

    let n_sym : Option<usize> =
        from_json_string(n_sym_json);

    if let (Some(c), Some(n)) =
        (codeword, n_sym)
    {

        let count =
            rs_error_count(&c, n);

        to_json_string(&count)
    } else {

        std::ptr::null_mut()
    }
}

// ============================================================================
// CRC-32
// ============================================================================

/// Computes CRC-32 checksum via JSON interface.
/// Input: [bytes]
/// Returns: u32 checksum
#[no_mangle]

pub unsafe extern "C" fn rssn_json_crc32_compute(
    data_json : *const c_char
) -> *mut c_char {

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    if let Some(d) = data {

        let crc = crc32_compute(&d);

        to_json_string(&crc)
    } else {

        std::ptr::null_mut()
    }
}

/// Verifies CRC-32 checksum via JSON interface.
/// Input: data as [bytes], expected_crc as u32
/// Returns: boolean
#[no_mangle]

pub unsafe extern "C" fn rssn_json_crc32_verify(
    data_json : *const c_char,
    expected_crc_json : *const c_char,
) -> *mut c_char {

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    let expected_crc : Option<u32> =
        from_json_string(
            expected_crc_json,
        );

    if let (Some(d), Some(crc)) =
        (data, expected_crc)
    {

        let valid =
            crc32_verify(&d, crc);

        to_json_string(&valid)
    } else {

        std::ptr::null_mut()
    }
}

/// Updates CRC-32 incrementally via JSON interface.
/// Input: current crc as u32, data as [bytes]
/// Returns: updated crc as u32
#[no_mangle]

pub unsafe extern "C" fn rssn_json_crc32_update(
    crc_json : *const c_char,
    data_json : *const c_char,
) -> *mut c_char {

    let crc : Option<u32> =
        from_json_string(crc_json);

    let data : Option<Vec<u8>> =
        from_json_string(data_json);

    if let (Some(c), Some(d)) =
        (crc, data)
    {

        let updated =
            crc32_update(c, &d);

        to_json_string(&updated)
    } else {

        std::ptr::null_mut()
    }
}

/// Finalizes CRC-32 computation via JSON interface.
/// Input: running crc as u32
/// Returns: final crc as u32
#[no_mangle]

pub unsafe extern "C" fn rssn_json_crc32_finalize(
    crc_json : *const c_char
) -> *mut c_char {

    let crc : Option<u32> =
        from_json_string(crc_json);

    if let Some(c) = crc {

        let final_crc =
            crc32_finalize(c);

        to_json_string(&final_crc)
    } else {

        std::ptr::null_mut()
    }
}
