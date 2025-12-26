//! Handle-based FFI API for error-correcting codes.
//!
//! This module provides C-compatible FFI functions for Hamming codes, Reed-Solomon codes,
//! and CRC-32, including encoding, decoding, and verification operations with error
//! detection and correction capabilities.

use crate::symbolic::error_correction::{
    crc32_compute, crc32_finalize, crc32_update, crc32_verify, hamming_check, hamming_decode,
    hamming_distance, hamming_encode, hamming_weight, rs_check, rs_decode, rs_encode,
    rs_error_count,
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
    error_pos: *mut u8,
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
    out_len: *mut usize,
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
    out_len: *mut usize,
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

// ============================================================================
// Hamming Code Extensions
// ============================================================================

/// Computes Hamming distance between two byte slices.
///
/// # Safety
/// Caller must ensure `a` and `b` point to `len` bytes each.
/// Returns -1 on error (null pointers or different lengths).
#[no_mangle]

pub unsafe extern "C" fn rssn_hamming_distance(
    a: *const u8,
    a_len: usize,
    b: *const u8,
    b_len: usize,
) -> i32 {

    if a.is_null() || b.is_null() {

        return -1;
    }

    if a_len != b_len {

        return -1;
    }

    let slice_a = std::slice::from_raw_parts(a, a_len);

    let slice_b = std::slice::from_raw_parts(b, b_len);

    match hamming_distance(slice_a, slice_b) {
        Some(dist) => dist as i32,
        None => -1,
    }
}

/// Computes Hamming weight (number of 1s) of a byte slice.
///
/// # Safety
/// Caller must ensure `data` points to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_hamming_weight(data: *const u8, len: usize) -> i32 {

    if data.is_null() {

        return -1;
    }

    let slice = std::slice::from_raw_parts(data, len);

    hamming_weight(slice) as i32
}

/// Checks if a Hamming(7,4) codeword is valid without correcting.
///
/// # Safety
/// Caller must ensure `codeword` points to 7 bytes.
/// Returns 1 if valid, 0 if invalid, -1 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_hamming_check(codeword: *const u8) -> i32 {

    if codeword.is_null() {

        return -1;
    }

    let slice = std::slice::from_raw_parts(codeword, 7);

    if hamming_check(slice) {

        1
    } else {

        0
    }
}

// ============================================================================
// Reed-Solomon Enhancements
// ============================================================================

/// Checks if a Reed-Solomon codeword is valid without attempting correction.
///
/// # Safety
/// Caller must ensure `codeword` points to `codeword_len` bytes.
/// Returns 1 if valid, 0 if invalid, -1 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_rs_check(
    codeword: *const u8,
    codeword_len: usize,
    n_sym: usize,
) -> i32 {

    if codeword.is_null() {

        return -1;
    }

    let slice = std::slice::from_raw_parts(codeword, codeword_len);

    if rs_check(slice, n_sym) {

        1
    } else {

        0
    }
}

/// Estimates the number of errors in a Reed-Solomon codeword.
///
/// # Safety
/// Caller must ensure `codeword` points to `codeword_len` bytes.
/// Returns error count or -1 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_rs_error_count(
    codeword: *const u8,
    codeword_len: usize,
    n_sym: usize,
) -> i32 {

    if codeword.is_null() {

        return -1;
    }

    let slice = std::slice::from_raw_parts(codeword, codeword_len);

    rs_error_count(slice, n_sym) as i32
}

// ============================================================================
// CRC-32
// ============================================================================

/// Computes CRC-32 checksum of data.
///
/// # Safety
/// Caller must ensure `data` points to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_crc32_compute(data: *const u8, len: usize) -> u32 {

    if data.is_null() {

        return 0;
    }

    let slice = std::slice::from_raw_parts(data, len);

    crc32_compute(slice)
}

/// Verifies CRC-32 checksum of data.
///
/// # Safety
/// Caller must ensure `data` points to `len` bytes.
/// Returns 1 if valid, 0 if invalid.
#[no_mangle]

pub unsafe extern "C" fn rssn_crc32_verify(data: *const u8, len: usize, expected_crc: u32) -> i32 {

    if data.is_null() {

        return 0;
    }

    let slice = std::slice::from_raw_parts(data, len);

    if crc32_verify(slice, expected_crc) {

        1
    } else {

        0
    }
}

/// Updates an existing CRC-32 with additional data (for incremental computation).
///
/// # Safety
/// Caller must ensure `data` points to `len` bytes.
/// Use 0xFFFFFFFF as initial crc for first call.
#[no_mangle]

pub unsafe extern "C" fn rssn_crc32_update(crc: u32, data: *const u8, len: usize) -> u32 {

    if data.is_null() {

        return crc;
    }

    let slice = std::slice::from_raw_parts(data, len);

    crc32_update(crc, slice)
}

/// Finalizes a CRC-32 computation started with crc32_update.
#[no_mangle]

pub extern "C" fn rssn_crc32_finalize(crc: u32) -> u32 {

    crc32_finalize(crc)
}
