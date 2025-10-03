//! # Error-Correcting Codes
//!
//! This module provides implementations for various error-correcting codes (ECC).
//! It includes functions for encoding data into codewords and decoding potentially
//! corrupted codewords back into data, with error detection and correction capabilities.
//! Specific implementations include Hamming codes and Reed-Solomon codes.

use crate::symbolic::error_correction_helper::{
    gf256_add, gf256_div, gf256_exp, gf256_inv, gf256_mul, poly_add_gf256, poly_div_gf256,
    poly_eval_gf256, poly_mul_gf256,
};

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

    Some(vec![p1, p2, d3, p4, d5, d6, d7])
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
            corrected_codeword[index] ^= 1; // Flip the bit
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

// =====================================================================================
// region: Reed-Solomon Codes
// =====================================================================================

/// Computes the generator polynomial for a Reed-Solomon code with `n_sym` error correction symbols.
pub(crate) fn rs_generator_poly(n_sym: usize) -> Result<Vec<u8>, String> {
    if n_sym == 0 {
        return Err("Number of symbols must be positive".to_string());
    }
    let mut g = vec![1];
    for i in 0..n_sym {
        let p = vec![1, gf256_exp(i as u8)];
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

    let remainder = poly_div_gf256(message_poly, &gen_poly);

    let mut codeword = data.to_vec();
    codeword.extend(remainder);

    Ok(codeword)
}

/// Calculates the syndromes of a received codeword.
pub(crate) fn rs_calc_syndromes(codeword_poly: &[u8], n_sym: usize) -> Vec<u8> {
    let mut syndromes = vec![0; n_sym];
    for i in 0..n_sym {
        syndromes[i] = poly_eval_gf256(codeword_poly, gf256_exp(i as u8));
    }
    syndromes
}

/// Finds the error locator polynomial `sigma` using the Berlekamp-Massey algorithm.
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

        let y = gf256_div(poly_eval_gf256(&omega, x_inv), sigma_prime_eval);
        codeword_poly[err_loc] = gf256_add(codeword_poly[err_loc], y);
    }

    Ok(codeword_poly[..codeword.len() - n_sym].to_vec())
}
