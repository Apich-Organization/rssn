//! Bincode-based FFI API for error-correcting codes.
//!
//! This module provides binary serialization-based FFI functions for Hamming codes and
//! Reed-Solomon codes, offering efficient binary data interchange for high-performance applications.

use crate::symbolic::error_correction::{
    hamming_encode, hamming_decode, rs_encode, rs_decode,
};
use crate::ffi_apis::common::*;

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
    n_sym_buf: BincodeBuffer
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
    n_sym_buf: BincodeBuffer
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
