//! Bincode-based FFI API for error-correcting codes.
//!
//! This module provides binary serialization-based FFI functions for Hamming codes,
//! Reed-Solomon codes, and CRC-32, offering efficient binary data interchange for
//! high-performance applications.

use crate::ffi_apis::common::*;
use crate::symbolic::error_correction::{
    crc32_compute,
    crc32_finalize,
    crc32_update,
    crc32_verify,
    hamming_check,
    hamming_decode,
    hamming_distance,
    hamming_encode,
    hamming_weight,
    rs_check,
    rs_decode,
    rs_encode,
    rs_error_count,
};

/// Encodes 4 data bits into a 7-bit Hamming(7,4) codeword via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_hamming_encode(data_buf: BincodeBuffer) -> BincodeBuffer {

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        match hamming_encode(&d) {
            Some(codeword) => to_bincode_buffer(&codeword),
            None => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Decodes a 7-bit Hamming(7,4) codeword via Bincode interface.
/// Returns tuple of (data, error_pos).
#[no_mangle]

pub extern "C" fn rssn_bincode_hamming_decode(codeword_buf: BincodeBuffer) -> BincodeBuffer {

    let codeword: Option<Vec<u8>> = from_bincode_buffer(&codeword_buf);

    if let Some(c) = codeword {

        match hamming_decode(&c) {
            Ok((data, error_pos)) => to_bincode_buffer(&(data, error_pos)),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Encodes data using Reed-Solomon code via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rs_encode(
    data_buf: BincodeBuffer,
    n_sym_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    let n_sym: Option<usize> = from_bincode_buffer(&n_sym_buf);

    if let (Some(d), Some(n)) = (data, n_sym) {

        match rs_encode(&d, n) {
            Ok(codeword) => to_bincode_buffer(&codeword),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Decodes a Reed-Solomon codeword via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rs_decode(
    codeword_buf: BincodeBuffer,
    n_sym_buf: BincodeBuffer,
) -> BincodeBuffer {

    let codeword: Option<Vec<u8>> = from_bincode_buffer(&codeword_buf);

    let n_sym: Option<usize> = from_bincode_buffer(&n_sym_buf);

    if let (Some(c), Some(n)) = (codeword, n_sym) {

        match rs_decode(&c, n) {
            Ok(data) => to_bincode_buffer(&data),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Hamming Code Extensions
// ============================================================================

/// Computes Hamming distance between two byte slices via Bincode interface.
/// Input: (a: Vec<u8>, b: Vec<u8>)
/// Returns: Option<usize>
#[no_mangle]

pub extern "C" fn rssn_bincode_hamming_distance(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<Vec<u8>> = from_bincode_buffer(&a_buf);

    let b: Option<Vec<u8>> = from_bincode_buffer(&b_buf);

    if let (Some(av), Some(bv)) = (a, b) {

        let dist = hamming_distance(&av, &bv);

        to_bincode_buffer(&dist)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Hamming weight of a byte slice via Bincode interface.
/// Input: Vec<u8>
/// Returns: usize
#[no_mangle]

pub extern "C" fn rssn_bincode_hamming_weight(data_buf: BincodeBuffer) -> BincodeBuffer {

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        let weight = hamming_weight(&d);

        to_bincode_buffer(&weight)
    } else {

        BincodeBuffer::empty()
    }
}

/// Checks if a Hamming(7,4) codeword is valid via Bincode interface.
/// Input: Vec<u8> (7 bytes)
/// Returns: bool
#[no_mangle]

pub extern "C" fn rssn_bincode_hamming_check(codeword_buf: BincodeBuffer) -> BincodeBuffer {

    let codeword: Option<Vec<u8>> = from_bincode_buffer(&codeword_buf);

    if let Some(c) = codeword {

        let valid = hamming_check(&c);

        to_bincode_buffer(&valid)
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Reed-Solomon Enhancements
// ============================================================================

/// Checks if a Reed-Solomon codeword is valid via Bincode interface.
/// Input: (codeword: Vec<u8>, n_sym: usize)
/// Returns: bool
#[no_mangle]

pub extern "C" fn rssn_bincode_rs_check(
    codeword_buf: BincodeBuffer,
    n_sym_buf: BincodeBuffer,
) -> BincodeBuffer {

    let codeword: Option<Vec<u8>> = from_bincode_buffer(&codeword_buf);

    let n_sym: Option<usize> = from_bincode_buffer(&n_sym_buf);

    if let (Some(c), Some(n)) = (codeword, n_sym) {

        let valid = rs_check(&c, n);

        to_bincode_buffer(&valid)
    } else {

        BincodeBuffer::empty()
    }
}

/// Estimates error count in a Reed-Solomon codeword via Bincode interface.
/// Input: (codeword: Vec<u8>, n_sym: usize)
/// Returns: usize
#[no_mangle]

pub extern "C" fn rssn_bincode_rs_error_count(
    codeword_buf: BincodeBuffer,
    n_sym_buf: BincodeBuffer,
) -> BincodeBuffer {

    let codeword: Option<Vec<u8>> = from_bincode_buffer(&codeword_buf);

    let n_sym: Option<usize> = from_bincode_buffer(&n_sym_buf);

    if let (Some(c), Some(n)) = (codeword, n_sym) {

        let count = rs_error_count(&c, n);

        to_bincode_buffer(&count)
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// CRC-32
// ============================================================================

/// Computes CRC-32 checksum via Bincode interface.
/// Input: Vec<u8>
/// Returns: u32
#[no_mangle]

pub extern "C" fn rssn_bincode_crc32_compute(data_buf: BincodeBuffer) -> BincodeBuffer {

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        let crc = crc32_compute(&d);

        to_bincode_buffer(&crc)
    } else {

        BincodeBuffer::empty()
    }
}

/// Verifies CRC-32 checksum via Bincode interface.
/// Input: (data: Vec<u8>, expected_crc: u32)
/// Returns: bool
#[no_mangle]

pub extern "C" fn rssn_bincode_crc32_verify(
    data_buf: BincodeBuffer,
    expected_crc_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    let expected_crc: Option<u32> = from_bincode_buffer(&expected_crc_buf);

    if let (Some(d), Some(crc)) = (data, expected_crc) {

        let valid = crc32_verify(&d, crc);

        to_bincode_buffer(&valid)
    } else {

        BincodeBuffer::empty()
    }
}

/// Updates CRC-32 incrementally via Bincode interface.
/// Input: (crc: u32, data: Vec<u8>)
/// Returns: u32
#[no_mangle]

pub extern "C" fn rssn_bincode_crc32_update(
    crc_buf: BincodeBuffer,
    data_buf: BincodeBuffer,
) -> BincodeBuffer {

    let crc: Option<u32> = from_bincode_buffer(&crc_buf);

    let data: Option<Vec<u8>> = from_bincode_buffer(&data_buf);

    if let (Some(c), Some(d)) = (crc, data) {

        let updated = crc32_update(c, &d);

        to_bincode_buffer(&updated)
    } else {

        BincodeBuffer::empty()
    }
}

/// Finalizes CRC-32 computation via Bincode interface.
/// Input: u32 (running crc)
/// Returns: u32 (final crc)
#[no_mangle]

pub extern "C" fn rssn_bincode_crc32_finalize(crc_buf: BincodeBuffer) -> BincodeBuffer {

    let crc: Option<u32> = from_bincode_buffer(&crc_buf);

    if let Some(c) = crc {

        let final_crc = crc32_finalize(c);

        to_bincode_buffer(&final_crc)
    } else {

        BincodeBuffer::empty()
    }
}
