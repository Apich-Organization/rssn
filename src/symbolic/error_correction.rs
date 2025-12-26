//! # Error-Correcting Codes
//!
//! This module provides implementations for various error-correcting codes (ECC).
//! It includes functions for encoding data into codewords and decoding potentially
//! corrupted codewords back into data, with error detection and correction capabilities.
//!
//! ## Supported Codes
//!
//! ### Hamming(7,4) - Single-bit Error Correction
//! - `hamming_encode` - Encode 4-bit data into 7-bit codeword
//! - `hamming_decode` - Decode and correct single-bit errors
//! - `hamming_distance` - Compute distance between two codewords
//! - `hamming_weight` - Count number of 1s in codeword
//! - `hamming_check` - Verify codeword validity without correction
//!
//! ### Reed-Solomon - Multi-symbol Error Correction
//! - `rs_encode` - Encode data with configurable error correction symbols
//! - `rs_decode` - Decode and correct errors in received codeword
//! - `rs_check` - Verify codeword validity without correction attempt
//! - `rs_error_count` - Estimate number of errors in codeword
//!
//! ### CRC-32 (IEEE 802.3) - Error Detection
//! - `crc32_compute` - Compute 32-bit checksum
//! - `crc32_verify` - Verify data against expected checksum
//! - `crc32_update` - Incremental CRC computation for streaming data
//! - `crc32_finalize` - Finalize incremental CRC computation
//!
//! ## Examples
//!
//! ### Hamming Code
//! ```
//! use rssn::symbolic::error_correction::{hamming_encode, hamming_decode, hamming_check};
//!
//! // Encode 4 bits of data
//! let data = vec![1, 0, 1, 1];
//! let codeword = hamming_encode(&data).unwrap();
//!
//! // Verify the codeword is valid
//! assert!(hamming_check(&codeword));
//!
//! // Decode (with optional error correction)
//! let (decoded, error_pos) = hamming_decode(&codeword).unwrap();
//! assert_eq!(decoded, data);
//! ```
//!
//! ### Reed-Solomon Code
//! ```
//! use rssn::symbolic::error_correction::{rs_encode, rs_decode, rs_check};
//!
//! let data = b"Hello".to_vec();
//! let codeword = rs_encode(&data, 4).unwrap();  // 4 error correction symbols
//!
//! assert!(rs_check(&codeword, 4));
//! let decoded = rs_decode(&codeword, 4).unwrap();
//! assert_eq!(decoded, data);
//! ```
//!
//! ### CRC-32
//! ```
//! use rssn::symbolic::error_correction::{crc32_compute, crc32_verify};
//!
//! let data = b"Important data";
//! let checksum = crc32_compute(data);
//!
//! assert!(crc32_verify(data, checksum));
//! ```

use crate::symbolic::error_correction_helper::{
    gf256_add, gf256_div, gf256_exp, gf256_inv, gf256_mul, poly_add_gf256, poly_div_gf256,
    poly_eval_gf256, poly_mul_gf256,
};

// ============================================================================
// Hamming Codes
// ============================================================================

/// Computes the Hamming distance between two byte slices.
///
/// The Hamming distance is the number of positions at which the corresponding
/// bits are different.
///
/// # Arguments
/// * `a` - First byte slice
/// * `b` - Second byte slice
///
/// # Returns
/// The Hamming distance, or `None` if slices have different lengths
#[must_use]
pub fn hamming_distance(a: &[u8], b: &[u8]) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }
    Some(a.iter().zip(b.iter()).filter(|(&x, &y)| x != y).count())
}

/// Computes the Hamming weight (number of 1s) of a byte slice.
///
/// For binary data where each byte is 0 or 1, this counts the number of 1s.
///
/// # Arguments
/// * `data` - Byte slice where each byte represents a bit (0 or 1)
///
/// # Returns
/// The count of non-zero bytes
#[must_use]
pub fn hamming_weight(data: &[u8]) -> usize {
    data.iter().filter(|&&x| x != 0).count()
}

/// Encodes a 4-bit data block into a 7-bit Hamming(7,4) codeword.
///
/// Hamming(7,4) is a single-error correcting code. It takes 4 data bits
/// and adds 3 parity bits to create a 7-bit codeword.
///
/// # Arguments
/// * `data` - A slice of 4 bytes, each representing a bit (0 or 1).
///
/// # Returns
/// A `Option<Vec<u8>>` of 7 bits representing the codeword, or `None` if the input length is not 4.
#[must_use]
pub fn hamming_encode(data: &[u8]) -> Option<Vec<u8>> {
    if data.len() != 4 {
        return None;
    }
    let d3 = data[0];
    let d5 = data[1];
    let d6 = data[2];
    let d7 = data[3];
    let p1 = d3 ^ d5 ^ d7;
    let p2 = d3 ^ d6 ^ d7;
    let p4 = d5 ^ d6 ^ d7;
    Some(vec![
        p1, p2, d3, p4, d5, d6, d7,
    ])
}

/// Checks if a Hamming(7,4) codeword is valid without correcting.
///
/// # Arguments
/// * `codeword` - A 7-bit codeword
///
/// # Returns
/// `true` if the codeword is valid (no errors), `false` otherwise
#[must_use]
pub fn hamming_check(codeword: &[u8]) -> bool {
    if codeword.len() != 7 {
        return false;
    }
    let p1_in = codeword[0];
    let p2_in = codeword[1];
    let d3_in = codeword[2];
    let p4_in = codeword[3];
    let d5_in = codeword[4];
    let d6_in = codeword[5];
    let d7_in = codeword[6];

    let p1_calc = d3_in ^ d5_in ^ d7_in;
    let p2_calc = d3_in ^ d6_in ^ d7_in;
    let p4_calc = d5_in ^ d6_in ^ d7_in;

    p1_in == p1_calc && p2_in == p2_calc && p4_in == p4_calc
}

/// Decodes a 7-bit Hamming(7,4) codeword, correcting a single-bit error if found.
///
/// This function calculates syndrome bits to detect and locate a single-bit error.
/// If an error is found, it flips the erroneous bit to correct the codeword.
///
/// # Arguments
/// * `codeword` - A slice of 7 bytes, each representing a bit (0 or 1).
///
/// # Returns
/// A `Result` containing:
/// - `Ok((data, error_pos))` where `data` is the 4-bit corrected data and `error_pos` is
///   `Some(index)` if an error was corrected at the given 1-based index, or `None` if no error was found.
/// - `Err(String)` if the input length is not 7.
pub fn hamming_decode(codeword: &[u8]) -> Result<(Vec<u8>, Option<usize>), String> {
    if codeword.len() != 7 {
        return Err("Codeword length must be 7".to_string());
    }
    let p1_in = codeword[0];
    let p2_in = codeword[1];
    let d3_in = codeword[2];
    let p4_in = codeword[3];
    let d5_in = codeword[4];
    let d6_in = codeword[5];
    let d7_in = codeword[6];
    let p1_calc = d3_in ^ d5_in ^ d7_in;
    let p2_calc = d3_in ^ d6_in ^ d7_in;
    let p4_calc = d5_in ^ d6_in ^ d7_in;
    let c1 = p1_in ^ p1_calc;
    let c2 = p2_in ^ p2_calc;
    let c4 = p4_in ^ p4_calc;
    let error_pos = (c4 << 2) | (c2 << 1) | c1;
    let mut corrected_codeword = codeword.to_vec();
    let error_index = if error_pos != 0 {
        let index = error_pos as usize - 1;
        if index < corrected_codeword.len() {
            corrected_codeword[index] ^= 1;
        }
        Some(error_pos as usize)
    } else {
        None
    };
    let corrected_data = vec![
        corrected_codeword[2],
        corrected_codeword[4],
        corrected_codeword[5],
        corrected_codeword[6],
    ];
    Ok((corrected_data, error_index))
}

// ============================================================================
// Reed-Solomon Codes
// ============================================================================

/// Computes the generator polynomial for a Reed-Solomon code with `n_sym` error correction symbols.
pub(crate) fn rs_generator_poly(n_sym: usize) -> Result<Vec<u8>, String> {
    if n_sym == 0 {
        return Err("Number of symbols must be positive".to_string());
    }
    let mut g = vec![1];
    for i in 0..n_sym {
        let p = vec![
            1,
            gf256_exp(i as u8),
        ];
        g = poly_mul_gf256(&g, &p);
    }
    Ok(g)
}

/// Encodes a data message using a Reed-Solomon code, adding `n_sym` error correction symbols.
///
/// Reed-Solomon codes are non-binary cyclic error-correcting codes. This function
/// appends `n_sym` zero bytes to the data message, divides the resulting polynomial
/// by the generator polynomial, and appends the remainder as parity symbols.
///
/// # Arguments
/// * `data` - The data message as a slice of `u8` bytes.
/// * `n_sym` - The number of error correction symbols to add.
///
/// # Returns
/// A `Result` containing the encoded codeword as a `Vec<u8>`, or an error string
/// if the message length exceeds the maximum allowed for the chosen code.
pub fn rs_encode(data: &[u8], n_sym: usize) -> Result<Vec<u8>, String> {
    if data.len() + n_sym > 255 {
        return Err("Message length + number of symbols cannot exceed 255".to_string());
    }
    let gen_poly = rs_generator_poly(n_sym)?;
    let mut message_poly = data.to_vec();
    message_poly.extend(vec![0; n_sym]);
    let remainder = poly_div_gf256(message_poly, &gen_poly)?;
    let mut codeword = data.to_vec();
    codeword.extend(remainder);
    Ok(codeword)
}

/// Calculates the syndromes of a received codeword.
pub(crate) fn rs_calc_syndromes(codeword_poly: &[u8], n_sym: usize) -> Vec<u8> {
    let mut syndromes = vec![0; n_sym];
    for (i, syndrome) in syndromes.iter_mut().enumerate().take(n_sym) {
        *syndrome = poly_eval_gf256(codeword_poly, gf256_exp(i as u8));
    }
    syndromes
}

/// Checks if a Reed-Solomon codeword is valid without attempting correction.
///
/// This is faster than full decoding when you only need to verify integrity.
///
/// # Arguments
/// * `codeword` - The received codeword
/// * `n_sym` - Number of error correction symbols
///
/// # Returns
/// `true` if the codeword is valid (all syndromes are zero)
#[must_use]
pub fn rs_check(codeword: &[u8], n_sym: usize) -> bool {
    let syndromes = rs_calc_syndromes(codeword, n_sym);
    syndromes.iter().all(|&s| s == 0)
}

/// Estimates the number of errors in a Reed-Solomon codeword.
///
/// This uses the syndrome values to estimate the error count.
///
/// # Arguments
/// * `codeword` - The received codeword
/// * `n_sym` - Number of error correction symbols
///
/// # Returns
/// Estimated number of errors (0 if codeword is valid)
#[must_use]
pub fn rs_error_count(codeword: &[u8], n_sym: usize) -> usize {
    let syndromes = rs_calc_syndromes(codeword, n_sym);
    if syndromes.iter().all(|&s| s == 0) {
        return 0;
    }
    let sigma = rs_find_error_locator_poly(&syndromes);
    sigma.len() - 1 // degree of error locator polynomial
}

/// Finds the error locator polynomial `sigma` using the Berlekamp-Massey algorithm.
#[allow(clippy::cast_possible_wrap)]
pub(crate) fn rs_find_error_locator_poly(syndromes: &[u8]) -> Vec<u8> {
    let mut sigma = vec![1];
    let mut prev_sigma = vec![1];
    let mut l = 0;
    let mut m = -1;
    let mut b = 1;
    for n in 0..syndromes.len() {
        let mut d = syndromes[n];
        for i in 1..=l {
            d = gf256_add(d, gf256_mul(sigma[sigma.len() - 1 - i], syndromes[n - i]));
        }
        if d != 0 {
            let t = sigma.clone();
            let mut correction = vec![b];
            correction.extend(vec![0; (n as i32 - m) as usize]);
            correction = poly_mul_gf256(&correction, &prev_sigma);
            sigma = poly_add_gf256(&sigma, &correction);
            if 2 * l <= n {
                l = n + 1 - l;
                m = n as i32;
                prev_sigma = t;
                b = d;
            }
        }
    }
    sigma
}

/// Finds the locations of errors by finding the roots of the error locator polynomial.
pub(crate) fn rs_find_error_locations(
    sigma: &[u8],
    codeword_len: usize,
) -> Result<Vec<usize>, String> {
    let mut error_locs = Vec::new();
    let err_poly_degree = sigma.len() - 1;
    for i in 0..codeword_len {
        let x = gf256_exp((255 - i) as u8);
        if poly_eval_gf256(sigma, x) == 0 {
            error_locs.push(i);
        }
    }
    if error_locs.len() != err_poly_degree {
        return Err("Failed to find the correct number of error locations.".to_string());
    }
    Ok(error_locs)
}

/// Decodes a Reed-Solomon codeword, correcting errors if found.
///
/// This function calculates syndromes, uses the Berlekamp-Massey algorithm to find
/// the error locator polynomial, determines error locations, and then calculates
/// error magnitudes to correct the corrupted symbols in the codeword.
///
/// # Arguments
/// * `codeword` - The received codeword as a slice of `u8` bytes.
/// * `n_sym` - The number of error correction symbols used during encoding.
///
/// # Returns
/// A `Result` containing the corrected data message as a `Vec<u8>`, or an error string
/// if error correction fails (e.g., too many errors).
pub fn rs_decode(codeword: &[u8], n_sym: usize) -> Result<Vec<u8>, String> {
    let mut codeword_poly = codeword.to_vec();
    let syndromes = rs_calc_syndromes(&codeword_poly, n_sym);
    if syndromes.iter().all(|&s| s == 0) {
        return Ok(codeword[..codeword.len() - n_sym].to_vec());
    }
    let sigma = rs_find_error_locator_poly(&syndromes);
    let error_locs = rs_find_error_locations(&sigma, codeword.len())?;
    let mut omega = poly_mul_gf256(&syndromes, &sigma);
    omega.truncate(n_sym);
    for &err_loc in &error_locs {
        let x_inv = gf256_inv(gf256_exp((codeword.len() - 1 - err_loc) as u8));
        let mut sigma_prime_eval = 0;
        for i in (1..sigma.len()).step_by(2) {
            sigma_prime_eval = gf256_add(sigma_prime_eval, sigma[i]);
        }
        let y = gf256_div(poly_eval_gf256(&omega, x_inv?), sigma_prime_eval);
        codeword_poly[err_loc] = gf256_add(codeword_poly[err_loc], y?);
    }
    Ok(codeword_poly[..codeword.len() - n_sym].to_vec())
}

// ============================================================================
// CRC-32
// ============================================================================

/// CRC-32 polynomial (IEEE 802.3)
const CRC32_POLYNOMIAL: u32 = 0xEDB88320;

/// Computes CRC-32 checksum of data.
///
/// Uses the standard IEEE 802.3 polynomial.
///
/// # Arguments
/// * `data` - The data to compute checksum for
///
/// # Returns
/// The 32-bit CRC checksum
#[must_use]
pub fn crc32_compute(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFFFFFF;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ CRC32_POLYNOMIAL;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}

/// Verifies CRC-32 checksum of data.
///
/// # Arguments
/// * `data` - The data to verify
/// * `expected_crc` - The expected CRC-32 checksum
///
/// # Returns
/// `true` if the computed CRC matches the expected value
#[must_use]
pub fn crc32_verify(data: &[u8], expected_crc: u32) -> bool {
    crc32_compute(data) == expected_crc
}

/// Updates an existing CRC-32 with additional data.
///
/// This allows computing CRC incrementally for streaming data.
///
/// # Arguments
/// * `crc` - Current CRC value (use 0xFFFFFFFF for initial call)
/// * `data` - Additional data to process
///
/// # Returns
/// Updated CRC value (call `!result` to get final CRC)
#[must_use]
pub fn crc32_update(crc: u32, data: &[u8]) -> u32 {
    let mut crc = crc;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ CRC32_POLYNOMIAL;
            } else {
                crc >>= 1;
            }
        }
    }
    crc
}

/// Finalizes a CRC-32 computation started with `crc32_update`.
///
/// # Arguments
/// * `crc` - The running CRC value from `crc32_update` calls
///
/// # Returns
/// The final CRC-32 checksum
#[must_use]
pub const fn crc32_finalize(crc: u32) -> u32 {
    !crc
}
