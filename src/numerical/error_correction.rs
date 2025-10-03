//! # Numerical Error Correction Codes
//!
//! This module provides numerical implementations of error correction codes,
//! specifically focusing on Reed-Solomon codes over GF(2^8). It includes
//! functions for encoding messages and decoding codewords to correct errors,
//! utilizing polynomial arithmetic over finite fields.

use crate::numerical::finite_field::{gf256_add, gf256_div, gf256_inv, gf256_mul};

/// Represents a polynomial over GF(2^8)
#[derive(Clone)]
pub struct PolyGF256(Vec<u8>);

impl PolyGF256 {
    pub(crate) fn eval(&self, x: u8) -> u8 {
        self.0
            .iter()
            .rfold(0, |acc, &coeff| gf256_add(gf256_mul(acc, x), coeff))
    }
}

impl PolyGF256 {
    pub(crate) fn degree(&self) -> usize {
        if self.0.is_empty() {
            0
        } else {
            self.0.len() - 1
        }
    }

    pub(crate) fn poly_add(&self, other: Self) -> Self {
        let mut result = vec![0; self.0.len().max(other.0.len())];
        //for i in 0..self.0.len() { result[i + result.len() - self.0.len()] = self.0[i]; }
        //for i in 0..other.0.len() { result[i + result.len() - other.0.len()] ^= other.0[i]; }
        // Assuming this block is inside a function where result is defined and mutable.
        //
        // Original code:
        // 20 | for i in 0..self.0.len() { result[i + result.len() - self.0.len()] = self.0[i]; }
        // 21 | for i in 0..other.0.len() { result[i + result.len() - other.0.len()] ^= other.0[i]; }

        // Fixed code:

        // 1. Calculate the length *before* the loops start.
        let result_len = result.len();

        // 2. Use the pre-calculated length in the first loop:
        for i in 0..self.0.len() {
            // result_len is a local primitive, not a borrow of 'result'.
            result[i + result_len - self.0.len()] = self.0[i];
        }

        // 3. Use the pre-calculated length in the second loop:
        for i in 0..other.0.len() {
            // result_len is a local primitive, not a borrow of 'result'.
            result[i + result_len - other.0.len()] ^= other.0[i];
        }
        PolyGF256(result)
    }

    pub(crate) fn poly_sub(&self, other: Self) -> Self {
        self.poly_add(other) // In GF(2^8), add is the same as sub
    }

    pub(crate) fn poly_mul(&self, other: Self) -> Self {
        let mut result = vec![0; self.degree() + other.degree() + 1];
        for i in 0..=self.degree() {
            for j in 0..=other.degree() {
                result[i + j] ^= gf256_mul(self.0[i], other.0[j]);
            }
        }
        PolyGF256(result)
    }

    pub(crate) fn poly_div(&self, divisor: Self) -> (Self, Self) {
        let mut rem = self.0.clone();
        let mut quot = vec![0; self.degree() + 1];
        let divisor_lead_inv = gf256_inv(*divisor.0.first().unwrap());

        while rem.len() >= divisor.0.len() {
            let lead_coeff = *rem.first().unwrap();
            let q_coeff = gf256_mul(lead_coeff, divisor_lead_inv);
            let deg_diff = rem.len() - divisor.0.len();
            quot[deg_diff] = q_coeff;

            for i in 0..divisor.0.len() {
                rem[i] ^= gf256_mul(divisor.0[i], q_coeff);
            }
            rem.remove(0);
        }
        (PolyGF256(quot), PolyGF256(rem))
    }

    pub(crate) fn derivative(&self) -> Self {
        let mut deriv = vec![0; self.degree()];
        for i in 1..=self.degree() {
            if i % 2 != 0 {
                deriv[i - 1] = self.0[i];
            }
        }
        PolyGF256(deriv)
    }
}

/// Encodes a message using Reed-Solomon codes over GF(2^8).
///
/// This function implements a systematic encoding scheme for Reed-Solomon codes.
/// It appends `n_parity` parity symbols to the message, which are computed by
/// evaluating the message polynomial at specific points in the finite field.
///
/// # Arguments
/// * `message` - A slice of bytes representing the message.
/// * `n_parity` - The number of parity symbols to add.
///
/// # Returns
/// A `Result` containing the full codeword (message + parity symbols), or an error
/// if the total length exceeds the field size.
pub fn reed_solomon_encode(message: &[u8], n_parity: usize) -> Result<Vec<u8>, String> {
    if message.len() + n_parity > 255 {
        return Err("Message + parity length cannot exceed 255".to_string());
    }
    // Interpret the message as coefficients of a polynomial
    let msg_poly = PolyGF256(message.to_vec());

    // The generator polynomial for n_parity symbols is (x-α^0)(x-α^1)...(x-α^{n_parity-1})
    // For simplicity, we will evaluate the message polynomial at n_parity points
    // to get the parity symbols directly. This is a systematic encoding.
    let mut codeword = message.to_vec();
    for i in 0..n_parity {
        // Evaluate P(α^i). α is 2, a generator of the field.
        let alpha_i = gf256_pow(2, i as u8);
        let parity_symbol = msg_poly.eval(alpha_i);
        codeword.push(parity_symbol);
    }
    Ok(codeword)
}

/// Decodes a Reed-Solomon codeword, correcting errors.
///
/// This implementation uses the Sugiyama algorithm (based on the Extended Euclidean Algorithm)
/// to find the error locator polynomial, Chien search to find error locations, and Forney's
/// algorithm to find error magnitudes. It corrects errors in-place within the `codeword`.
///
/// # Arguments
/// * `codeword` - A mutable slice of bytes representing the received codeword.
/// * `n_parity` - The number of parity symbols in the original encoding.
///
/// # Returns
/// A `Result` indicating success or an error if decoding fails (e.g., too many errors).
pub fn reed_solomon_decode(codeword: &mut [u8], n_parity: usize) -> Result<(), String> {
    // 1. Calculate Syndromes
    let syndromes = calculate_syndromes(codeword, n_parity);
    if syndromes.iter().all(|&s| s == 0) {
        return Ok(()); // No errors detected
    }

    // 2. Find Error Locator Polynomial using Sugiyama's Algorithm (based on EEA)
    let (sigma, omega) = find_error_locator_poly(&syndromes, n_parity)?;

    // 3. Find error locations using Chien Search
    let error_locations = chien_search(&sigma);
    if error_locations.is_empty() {
        return Err("Failed to find error locations.".to_string());
    }

    // 4. Find error magnitudes using Forney's Algorithm
    let error_magnitudes = forney_algorithm(&omega, &sigma, &error_locations);

    // 5. Correct the errors
    for (i, loc) in error_locations.iter().enumerate() {
        let codeword_pos = codeword.len() - 1 - (*loc as usize);
        codeword[codeword_pos] = gf256_add(codeword[codeword_pos], error_magnitudes[i]);
    }

    Ok(())
}

/// Calculates the syndromes of a received codeword.
pub(crate) fn calculate_syndromes(codeword: &[u8], n_parity: usize) -> Vec<u8> {
    let received_poly = PolyGF256(codeword.to_vec());
    let mut syndromes = Vec::with_capacity(n_parity);
    for i in 0..n_parity {
        let alpha_i = gf256_pow(2, i as u8);
        syndromes.push(received_poly.eval(alpha_i));
    }
    syndromes
}

/// Finds the error locator and evaluator polynomials using the Extended Euclidean Algorithm.
pub(crate) fn find_error_locator_poly(
    syndromes: &[u8],
    n_parity: usize,
) -> Result<(PolyGF256, PolyGF256), String> {
    let s = PolyGF256(syndromes.iter().rev().cloned().collect());
    let mut z_k = vec![0u8; n_parity + 1];
    z_k[n_parity] = 1;
    let z = PolyGF256(z_k);

    let (mut r_prev, mut r_curr) = (z, s);
    let (mut t_prev, mut t_curr) = (PolyGF256(vec![0]), PolyGF256(vec![1]));

    while r_curr.degree() >= n_parity / 2 {
        let (q, r_next) = r_prev.poly_div(r_curr.clone());
        let t_next = t_prev.poly_sub(q.poly_mul(t_curr.clone()));

        r_prev = r_curr;
        r_curr = r_next;

        t_prev = t_curr;
        t_curr = t_next;
    }
    // sigma = t_curr, omega = r_curr
    Ok((t_curr, r_curr))
}

/// Finds the roots of the error locator polynomial to determine error locations.
pub(crate) fn chien_search(sigma: &PolyGF256) -> Vec<u8> {
    let mut error_locs = Vec::new();
    for i in 0..255 {
        let alpha_inv = gf256_inv(gf256_pow(2, i));
        if sigma.eval(alpha_inv) == 0 {
            error_locs.push(i);
        }
    }
    error_locs
}

/// Computes error magnitudes using Forney's algorithm.
pub(crate) fn forney_algorithm(omega: &PolyGF256, sigma: &PolyGF256, error_locs: &[u8]) -> Vec<u8> {
    let sigma_prime = sigma.derivative();
    let mut magnitudes = Vec::new();
    for &loc in error_locs {
        let x_inv = gf256_inv(gf256_pow(2, loc));
        let omega_val = omega.eval(x_inv);
        let sigma_prime_val = sigma_prime.eval(x_inv);
        let magnitude = gf256_div(gf256_mul(omega_val, x_inv), sigma_prime_val);
        magnitudes.push(magnitude);
    }
    magnitudes
}

pub(crate) fn gf256_pow(base: u8, exp: u8) -> u8 {
    let mut res = 1;
    let mut b = base;
    let mut e = exp;
    while e > 0 {
        if e % 2 == 1 {
            res = gf256_mul(res, b);
        }
        b = gf256_mul(b, b);
        e /= 2;
    }
    res
}
