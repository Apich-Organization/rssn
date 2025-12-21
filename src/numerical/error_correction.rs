//! # Numerical Error Correction Codes
//!
//! This module provides numerical implementations of error correction codes,
//! specifically focusing on Reed-Solomon codes over GF(2^8). It includes
//! functions for encoding messages and decoding codewords to correct errors,
//! utilizing polynomial arithmetic over finite fields.
use crate::numerical::finite_field::{gf256_add, gf256_div, gf256_inv, gf256_mul, gf256_pow};
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
    pub(crate) const fn degree(&self) -> usize {
        if self.0.is_empty() {
            0
        } else {
            self.0.len() - 1
        }
    }
    pub(crate) fn poly_add(&self, other: &Self) -> Self {
        let mut result = vec![0; self.0.len().max(other.0.len())];
        let result_len = result.len();
        for i in 0..self.0.len() {
            result[i + result_len - self.0.len()] = self.0[i];
        }
        for i in 0..other.0.len() {
            result[i + result_len - other.0.len()] ^= other.0[i];
        }
        PolyGF256(result)
    }
    pub(crate) fn poly_sub(&self, other: &Self) -> Self {
        self.poly_add(other)
    }
    pub(crate) fn poly_mul(&self, other: &Self) -> Self {
        let mut result = vec![0; self.degree() + other.degree() + 1];
        for i in 0..=self.degree() {
            for j in 0..=other.degree() {
                result[i + j] ^= gf256_mul(self.0[i], other.0[j]);
            }
        }
        PolyGF256(result)
    }
    pub(crate) fn poly_div(&self, divisor: &Self) -> Result<(Self, Self), String> {
        if divisor.0.is_empty() {
            return Err("Division by zero polynomial".to_string());
        }
        let mut rem = self.0.clone();
        let mut quot = vec![0; self.degree() + 1];
        let divisor_lead_inv = gf256_inv(divisor.0[0])?;
        while rem.len() >= divisor.0.len() {
            let lead_coeff = rem[0];
            let q_coeff = gf256_mul(lead_coeff, divisor_lead_inv);
            let deg_diff = rem.len() - divisor.0.len();
            quot[deg_diff] = q_coeff;
            for (i, var) in rem.iter_mut().enumerate().take(divisor.0.len()) {
                *var ^= gf256_mul(divisor.0[i], q_coeff);
            }
            rem.remove(0);
        }
        Ok((PolyGF256(quot), PolyGF256(rem)))
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
    let msg_poly = PolyGF256(message.to_vec());
    let mut codeword = message.to_vec();
    for i in 0..n_parity {
        let alpha_i = gf256_pow(2, i as u64);
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
    let syndromes = calculate_syndromes(codeword, n_parity);
    if syndromes.iter().all(|&s| s == 0) {
        return Ok(());
    }
    let (sigma, omega) = find_error_locator_poly(&syndromes, n_parity)?;
    let error_locations = chien_search(&sigma)?;
    if error_locations.is_empty() {
        return Err("Failed to find error locations.".to_string());
    }
    let error_magnitudes = forney_algorithm(&omega, &sigma, &error_locations)?;
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
        let alpha_i = gf256_pow(2, i as u64);
        syndromes.push(received_poly.eval(alpha_i));
    }
    syndromes
}
/// Finds the error locator and evaluator polynomials using the Extended Euclidean Algorithm.
pub(crate) fn find_error_locator_poly(
    syndromes: &[u8],
    n_parity: usize,
) -> Result<(PolyGF256, PolyGF256), String> {
    let s = PolyGF256(syndromes.iter().rev().copied().collect());
    let mut z_k = vec![0u8; n_parity + 1];
    z_k[n_parity] = 1;
    let z = PolyGF256(z_k);
    let (mut r_prev, mut r_curr) = (z, s);
    let (mut t_prev, mut t_curr) = (PolyGF256(vec![0]), PolyGF256(vec![1]));
    while r_curr.degree() >= n_parity / 2 {
        let (q, r_next) = r_prev.poly_div(&r_curr.clone())?;
        let t_next = t_prev.poly_sub(&q.poly_mul(&t_curr.clone()));
        r_prev = r_curr;
        r_curr = r_next;
        t_prev = t_curr;
        t_curr = t_next;
    }
    Ok((t_curr, r_curr))
}
/// Finds the roots of the error locator polynomial to determine error locations.
pub(crate) fn chien_search(sigma: &PolyGF256) -> Result<Vec<u8>, String> {
    let mut error_locs = Vec::new();
    for i in 0..255u8 {
        let alpha_inv = gf256_inv(gf256_pow(2, i as u64))?;
        if sigma.eval(alpha_inv) == 0 {
            error_locs.push(i);
        }
    }
    Ok(error_locs)
}
/// Computes error magnitudes using Forney's algorithm.
pub(crate) fn forney_algorithm(
    omega: &PolyGF256,
    sigma: &PolyGF256,
    error_locs: &[u8],
) -> Result<Vec<u8>, String> {
    let sigma_prime = sigma.derivative();
    let mut magnitudes = Vec::new();
    for &loc in error_locs {
        let x_inv = gf256_inv(gf256_pow(2, loc as u64))?;
        let omega_val = omega.eval(x_inv);
        let sigma_prime_val = sigma_prime.eval(x_inv);
        let magnitude = gf256_div(gf256_mul(omega_val, x_inv), sigma_prime_val)?;
        magnitudes.push(magnitude);
    }
    Ok(magnitudes)
}
