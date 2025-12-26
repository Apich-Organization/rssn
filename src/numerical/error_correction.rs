//! # Numerical Error Correction Codes
//!
//! This module provides numerical implementations of error correction codes,
//! specifically focusing on Reed-Solomon codes over GF(2^8), Hamming codes,
//! BCH codes, and CRC checksums. It includes functions for encoding messages
//! and decoding codewords to correct errors, utilizing polynomial arithmetic
//! over finite fields.
//!
//! ## Supported Codes
//!
//! ### Reed-Solomon Codes
//! - `reed_solomon_encode` - Encode data with configurable error correction symbols
//! - `reed_solomon_decode` - Decode and correct errors in received codeword
//! - `reed_solomon_check` - Verify codeword validity
//!
//! ### Hamming Codes
//! - `hamming_encode_numerical` - Encode data into Hamming(7,4) codeword
//! - `hamming_decode_numerical` - Decode with single-bit error correction
//! - `hamming_distance_numerical` - Compute distance between codewords
//! - `hamming_weight_numerical` - Count number of 1s in codeword
//!
//! ### BCH Codes
//! - `bch_encode` - BCH encoding with multiple error correction
//! - `bch_decode` - BCH decoding with error correction
//!
//! ### CRC Checksums
//! - `crc32_compute_numerical` - Compute CRC-32 checksum
//! - `crc32_verify_numerical` - Verify data against checksum
//! - `crc16_compute` - Compute CRC-16 checksum
//! - `crc8_compute` - Compute CRC-8 checksum
//!
//! ## Examples
//!
//! ### Reed-Solomon Encoding/Decoding
//! ```
//! 
//! use rssn::numerical::error_correction::{
//!     reed_solomon_decode,
//!     reed_solomon_encode,
//! };
//!
//! let message = vec![
//!     0x01, 0x02, 0x03, 0x04,
//! ];
//!
//! let codeword =
//!     reed_solomon_encode(&message, 4).unwrap();
//!
//! // Introduce an error
//! let mut corrupted = codeword.clone();
//!
//! corrupted[0] ^= 0xFF;
//!
//! // Decode and correct
//! reed_solomon_decode(&mut corrupted, 4).unwrap();
//!
//! assert_eq!(
//!     &corrupted[..4],
//!     &message
//! );
//! ```
//!
//! ### CRC-32 Checksum
//! ```
//! 
//! use rssn::numerical::error_correction::{
//!     crc32_compute_numerical,
//!     crc32_verify_numerical,
//! };
//!
//! let data = b"Hello, World!";
//!
//! let checksum = crc32_compute_numerical(data);
//!
//! assert!(crc32_verify_numerical(data, checksum));
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::finite_field::gf256_add;
use crate::numerical::finite_field::gf256_div;
use crate::numerical::finite_field::gf256_inv;
use crate::numerical::finite_field::gf256_mul;
use crate::numerical::finite_field::gf256_pow;

// ============================================================================
// Polynomial over GF(2^8)
// ============================================================================

/// Represents a polynomial over GF(2^8).
///
/// The polynomial is stored in descending order of powers, i.e., the first
/// element is the coefficient of the highest power term.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]

pub struct PolyGF256(pub Vec<u8>);

impl PolyGF256 {
    /// Creates a new polynomial from coefficients.
    ///
    /// # Arguments
    /// * `coeffs` - Coefficients in descending order of powers.
    #[must_use]

    pub const fn new(coeffs : Vec<u8>) -> Self {

        Self(coeffs)
    }

    /// Returns the degree of the polynomial.
    #[must_use]

    pub const fn degree(&self) -> usize {

        if self.0.is_empty() {

            0
        } else {

            self.0.len() - 1
        }
    }

    /// Evaluates the polynomial at a given point using Horner's method.
    ///
    /// # Arguments
    /// * `x` - The evaluation point in GF(2^8).
    ///
    /// # Returns
    /// The value p(x) in GF(2^8).
    #[must_use]

    pub fn eval(
        &self,
        x : u8,
    ) -> u8 {

        self.0
            .iter()
            .rfold(0, |acc, &coeff| {

                gf256_add(
                    gf256_mul(acc, x),
                    coeff,
                )
            })
    }

    /// Adds two polynomials over GF(2^8).
    ///
    /// In GF(2^8), addition is XOR.
    #[must_use]

    pub fn poly_add(
        &self,
        other : &Self,
    ) -> Self {

        let mut result = vec![
            0;
            self.0
                .len()
                .max(other.0.len())
        ];

        let result_len = result.len();

        for i in 0 .. self.0.len() {

            result[i + result_len - self.0.len()] = self.0[i];
        }

        for i in 0 .. other.0.len() {

            result[i + result_len - other.0.len()] ^= other.0[i];
        }

        Self(result)
    }

    /// Subtracts two polynomials over GF(2^8).
    ///
    /// In GF(2^8), subtraction is the same as addition (XOR).
    #[must_use]

    pub fn poly_sub(
        &self,
        other : &Self,
    ) -> Self {

        self.poly_add(other)
    }

    /// Multiplies two polynomials over GF(2^8).
    #[must_use]

    pub fn poly_mul(
        &self,
        other : &Self,
    ) -> Self {

        let mut result = vec![0; self.degree() + other.degree() + 1];

        for i in 0 ..= self.degree() {

            for j in 0 ..= other.degree() {

                result[i + j] ^= gf256_mul(
                    self.0[i],
                    other.0[j],
                );
            }
        }

        Self(result)
    }

    /// Divides two polynomials over GF(2^8).
    ///
    /// # Arguments
    /// * `divisor` - The divisor polynomial.
    ///
    /// # Returns
    /// A tuple (quotient, remainder), or an error if divisor is zero.

    pub fn poly_div(
        &self,
        divisor : &Self,
    ) -> Result<(Self, Self), String> {

        if divisor.0.is_empty() {

            return Err("Division by zero polynomial".to_string());
        }

        let mut rem = self.0.clone();

        let mut quot = vec![0; self.degree() + 1];

        let divisor_lead_inv = gf256_inv(divisor.0[0])?;

        while rem.len() >= divisor.0.len() {

            let lead_coeff = rem[0];

            let q_coeff = gf256_mul(
                lead_coeff,
                divisor_lead_inv,
            );

            let deg_diff = rem.len() - divisor.0.len();

            quot[deg_diff] = q_coeff;

            for (i, var) in rem
                .iter_mut()
                .enumerate()
                .take(divisor.0.len())
            {

                *var ^= gf256_mul(
                    divisor.0[i],
                    q_coeff,
                );
            }

            rem.remove(0);
        }

        Ok((
            Self(quot),
            Self(rem),
        ))
    }

    /// Returns the derivative of the polynomial over GF(2^8).
    ///
    /// In GF(2^8), the derivative has a special property:
    /// d/dx(x^n) = n * x^(n-1), where n is reduced mod 2.
    /// Thus, only odd-power terms contribute to the derivative.
    #[must_use]

    pub fn derivative(&self) -> Self {

        let mut deriv = vec![0; self.degree()];

        for i in 1 ..= self.degree() {

            if i % 2 != 0 {

                deriv[i - 1] = self.0[i];
            }
        }

        Self(deriv)
    }

    /// Scales the polynomial by a constant in GF(2^8).
    #[must_use]

    pub fn scale(
        &self,
        c : u8,
    ) -> Self {

        Self(
            self.0
                .iter()
                .map(|&coeff| gf256_mul(coeff, c))
                .collect(),
        )
    }

    /// Normalizes the polynomial by removing leading zeros.
    #[must_use]

    pub fn normalize(&self) -> Self {

        let start = self
            .0
            .iter()
            .position(|&x| x != 0)
            .unwrap_or(self.0.len());

        Self(self.0[start ..].to_vec())
    }
}

// ============================================================================
// Reed-Solomon Codes
// ============================================================================

/// Computes the generator polynomial for a Reed-Solomon code with `n_parity` error correction symbols.
///
/// The generator polynomial is the product of (x - α^i) for i = 0 to n_parity-1,
/// where α is the primitive element (2) of GF(2^8).

fn rs_generator_poly(n_parity : usize) -> Result<Vec<u8>, String> {

    if n_parity == 0 {

        return Err("Number of parity symbols must be positive".to_string());
    }

    let mut g = vec![1u8];

    for i in 0 .. n_parity {

        // Multiply by (x - α^i), which in GF(2^8) is (x + α^i) since subtraction = addition
        let root = gf256_pow(2, i as u64);

        let factor = vec![1u8, root];

        g = poly_mul_gf256(&g, &factor);
    }

    Ok(g)
}

/// Multiplies two polynomials over GF(2^8).

fn poly_mul_gf256(
    p1 : &[u8],
    p2 : &[u8],
) -> Vec<u8> {

    if p1.is_empty() || p2.is_empty() {

        return vec![];
    }

    let mut result = vec![0; p1.len() + p2.len() - 1];

    for i in 0 .. p1.len() {

        for j in 0 .. p2.len() {

            result[i + j] ^= gf256_mul(p1[i], p2[j]);
        }
    }

    result
}

/// Performs polynomial long division over GF(2^8).
/// Returns the remainder after dividing dividend by divisor.

fn poly_div_gf256(
    mut dividend : Vec<u8>,
    divisor : &[u8],
) -> Result<Vec<u8>, String> {

    if divisor.is_empty() {

        return Err("Divisor cannot be empty".to_string());
    }

    let divisor_len = divisor.len();

    let lead_divisor = divisor[0];

    let lead_divisor_inv = gf256_inv(lead_divisor)?;

    while dividend.len() >= divisor_len {

        let lead_dividend = dividend[0];

        if lead_dividend == 0 {

            dividend.remove(0);

            continue;
        }

        let coeff = gf256_mul(
            lead_dividend,
            lead_divisor_inv,
        );

        for i in 0 .. divisor_len {

            let term = gf256_mul(coeff, divisor[i]);

            dividend[i] ^= term;
        }

        dividend.remove(0);
    }

    // Pad remainder to have n_parity bytes
    while dividend.len() < divisor_len - 1 {

        dividend.insert(0, 0);
    }

    Ok(dividend)
}

/// Evaluates a polynomial over GF(2^8) at a given point.

fn poly_eval_gf256(
    poly : &[u8],
    x : u8,
) -> u8 {

    let mut y = 0u8;

    for &coeff in poly {

        y = gf256_mul(y, x) ^ coeff;
    }

    y
}

/// Encodes a message using Reed-Solomon codes over GF(2^8).
///
/// This function implements a systematic encoding scheme for Reed-Solomon codes.
/// It appends `n_parity` parity symbols to the message. The parity symbols are
/// computed by dividing the message polynomial (shifted by `n_parity` positions)
/// by the generator polynomial and taking the remainder.
///
/// # Arguments
/// * `message` - A slice of bytes representing the message.
/// * `n_parity` - The number of parity symbols to add.
///
/// # Returns
/// A `Result` containing the full codeword (message + parity symbols), or an error
/// if the total length exceeds the field size.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::reed_solomon_encode;
///
/// let message = vec![0x01, 0x02, 0x03];
///
/// let codeword = reed_solomon_encode(&message, 4).unwrap();
///
/// assert_eq!(codeword.len(), 7); // 3 data + 4 parity
/// ```

pub fn reed_solomon_encode(
    message : &[u8],
    n_parity : usize,
) -> Result<Vec<u8>, String> {

    if message.len() + n_parity > 255 {

        return Err("Message + parity length cannot exceed 255".to_string());
    }

    if n_parity == 0 {

        return Ok(message.to_vec());
    }

    // Generate the generator polynomial
    let gen_poly = rs_generator_poly(n_parity)?;

    // Shift message polynomial by n_parity positions (multiply by x^n_parity)
    let mut message_poly = message.to_vec();

    message_poly.extend(vec![0; n_parity]);

    // Compute remainder = message_poly mod gen_poly
    let remainder = poly_div_gf256(
        message_poly,
        &gen_poly,
    )?;

    // Codeword = message + remainder (systematic encoding)
    let mut codeword = message.to_vec();

    codeword.extend(&remainder);

    Ok(codeword)
}

/// Decodes a Reed-Solomon codeword, correcting errors.
///
/// This implementation uses the Berlekamp-Massey algorithm to find the error locator
/// polynomial, Chien search to find error locations, and Forney's algorithm to find
/// error magnitudes. It corrects errors in-place within the `codeword`.
///
/// # Arguments
/// * `codeword` - A mutable slice of bytes representing the received codeword.
/// * `n_parity` - The number of parity symbols in the original encoding.
///
/// # Returns
/// A `Result` indicating success or an error if decoding fails (e.g., too many errors).
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::{
///     reed_solomon_decode,
///     reed_solomon_encode,
/// };
///
/// let message = vec![0x01, 0x02, 0x03];
///
/// let mut codeword =
///     reed_solomon_encode(&message, 4).unwrap();
///
/// codeword[0] ^= 0xFF; // Introduce error
/// reed_solomon_decode(&mut codeword, 4).unwrap();
///
/// assert_eq!(
///     &codeword[..3],
///     &message
/// );
/// ```

pub fn reed_solomon_decode(
    codeword : &mut [u8],
    n_parity : usize,
) -> Result<(), String> {

    let syndromes = calculate_syndromes(codeword, n_parity);

    if syndromes
        .iter()
        .all(|&s| s == 0)
    {

        return Ok(());
    }

    // Use Berlekamp-Massey to find error locator polynomial
    let sigma = berlekamp_massey(&syndromes);

    // Use Chien search to find error locations
    let error_locations = chien_search_extended(
        &sigma,
        codeword.len(),
    )?;

    if error_locations.is_empty() {

        return Err("Failed to find error locations.".to_string());
    }

    // Compute error evaluator polynomial omega = S(x) * sigma(x) mod x^n_parity
    let mut omega = poly_mul_gf256(&syndromes, &sigma);

    if omega.len() > n_parity {

        omega = omega[omega.len() - n_parity ..].to_vec();
    }

    // Use Forney's algorithm to find error magnitudes
    let error_magnitudes = forney_algorithm_extended(
        &omega,
        &sigma,
        &error_locations,
        codeword.len(),
    )?;

    // Correct errors
    for (i, &loc) in error_locations
        .iter()
        .enumerate()
    {

        codeword[loc] ^= error_magnitudes[i];
    }

    Ok(())
}

/// Berlekamp-Massey algorithm to find the error locator polynomial.

fn berlekamp_massey(syndromes : &[u8]) -> Vec<u8> {

    let n = syndromes.len();

    let mut sigma = vec![1u8]; // Error locator polynomial
    let mut b = vec![1u8]; // Previous error locator polynomial
    let mut l = 0usize; // Number of errors found

    for i in 0 .. n {

        // Compute discrepancy
        let mut delta = syndromes[i];

        for j in 1 ..= l.min(sigma.len() - 1) {

            delta ^= gf256_mul(
                sigma[j],
                syndromes[i - j],
            );
        }

        // Shift b: b = x * b
        b.push(0);

        if delta != 0 {

            let t = sigma.clone();

            // sigma = sigma - delta * x^(i-m) * b
            // where x^(i-m) * b is already stored in b
            let scaled_b : Vec<u8> = b
                .iter()
                .map(|&c| gf256_mul(c, delta))
                .collect();

            sigma = poly_add_gf256(&sigma, &scaled_b);

            if 2 * l <= i {

                l = i + 1 - l;

                // b = t / delta
                let delta_inv = gf256_inv(delta).unwrap_or(0);

                b = t
                    .iter()
                    .map(|&c| gf256_mul(c, delta_inv))
                    .collect();
            }
        }
    }

    sigma
}

/// Add two polynomials over GF(2^8).

fn poly_add_gf256(
    p1 : &[u8],
    p2 : &[u8],
) -> Vec<u8> {

    let max_len = p1
        .len()
        .max(p2.len());

    let mut result = vec![0u8; max_len];

    for (i, &c) in p1
        .iter()
        .enumerate()
    {

        result[max_len - p1.len() + i] ^= c;
    }

    for (i, &c) in p2
        .iter()
        .enumerate()
    {

        result[max_len - p2.len() + i] ^= c;
    }

    result
}

/// Extended Chien search to find error locations.

fn chien_search_extended(
    sigma : &[u8],
    codeword_len : usize,
) -> Result<Vec<usize>, String> {

    let mut error_locs = Vec::new();

    let num_errors = sigma.len() - 1;

    for i in 0 .. codeword_len {

        // Evaluate sigma at alpha^(-i) = alpha^(255-i)
        let x = gf256_pow(
            2,
            ((255 - i) % 255) as u64,
        );

        let eval = poly_eval_gf256(sigma, x);

        if eval == 0 {

            error_locs.push(i);
        }
    }

    if error_locs.len() != num_errors && num_errors > 0 {

        return Err(format!(
            "Found {} error locations, expected {}",
            error_locs.len(),
            num_errors
        ));
    }

    Ok(error_locs)
}

/// Extended Forney's algorithm to compute error magnitudes.

fn forney_algorithm_extended(
    omega : &[u8],
    sigma : &[u8],
    error_locs : &[usize],
    codeword_len : usize,
) -> Result<Vec<u8>, String> {

    // Compute formal derivative of sigma
    let sigma_prime = poly_derivative_gf256(sigma);

    let mut magnitudes = Vec::new();

    for &loc in error_locs {

        // X_i = alpha^loc, X_i^(-1) = alpha^(-loc) = alpha^(255-loc)
        let x_inv = gf256_pow(
            2,
            ((255 - loc) % 255) as u64,
        );

        let omega_val = poly_eval_gf256(omega, x_inv);

        let sigma_prime_val = poly_eval_gf256(&sigma_prime, x_inv);

        if sigma_prime_val == 0 {

            return Err("Division by zero in Forney algorithm".to_string());
        }

        let magnitude = gf256_div(
            omega_val,
            sigma_prime_val,
        )?;

        magnitudes.push(magnitude);
    }

    Ok(magnitudes)
}

/// Compute formal derivative of polynomial in GF(2^8).

fn poly_derivative_gf256(poly : &[u8]) -> Vec<u8> {

    if poly.len() <= 1 {

        return vec![0];
    }

    // In GF(2), derivative of x^n is n*x^(n-1), and n mod 2 matters
    // Only odd-indexed coefficients survive
    let n = poly.len() - 1;

    let mut result = Vec::new();

    // poly = a_n*x^n + a_(n-1)*x^(n-1) + ... + a_1*x + a_0
    // derivative = n*a_n*x^(n-1) + ... + a_1
    for (i, &coeff) in poly
        .iter()
        .enumerate()
    {

        let power = n - i;

        if power > 0 && power % 2 == 1 {

            result.push(coeff);
        } else if power > 0 {

            result.push(0);
        }
    }

    // Remove leading zeros
    while result.len() > 1 && result[0] == 0 {

        result.remove(0);
    }

    if result.is_empty() {

        result.push(0);
    }

    result
}

/// Checks if a Reed-Solomon codeword is valid.
///
/// # Arguments
/// * `codeword` - The codeword to check.
/// * `n_parity` - The number of parity symbols.
///
/// # Returns
/// `true` if the codeword is valid (all syndromes are zero).
#[must_use]

pub fn reed_solomon_check(
    codeword : &[u8],
    n_parity : usize,
) -> bool {

    let syndromes = calculate_syndromes(codeword, n_parity);

    syndromes
        .iter()
        .all(|&s| s == 0)
}

/// Calculates the syndromes of a received codeword.
///
/// The syndrome `S_i` is computed as the codeword polynomial evaluated at α^i.
#[must_use]

pub fn calculate_syndromes(
    codeword : &[u8],
    n_parity : usize,
) -> Vec<u8> {

    let mut syndromes = Vec::with_capacity(n_parity);

    for i in 0 .. n_parity {

        let alpha_i = gf256_pow(2, i as u64);

        syndromes.push(poly_eval_gf256(
            codeword,
            alpha_i,
        ));
    }

    syndromes
}

/// Finds the roots of the error locator polynomial to determine error locations.
///
/// Uses Chien search to efficiently evaluate the polynomial at all field elements.

pub fn chien_search(sigma : &PolyGF256) -> Result<Vec<u8>, String> {

    let mut error_locs = Vec::new();

    for i in 0 .. 255u8 {

        let alpha_inv = gf256_inv(gf256_pow(
            2,
            u64::from(i),
        ))?;

        if sigma.eval(alpha_inv) == 0 {

            error_locs.push(i);
        }
    }

    Ok(error_locs)
}

/// Computes error magnitudes using Forney's algorithm.
///
/// # Arguments
/// * `omega` - The error evaluator polynomial.
/// * `sigma` - The error locator polynomial.
/// * `error_locs` - The positions of errors in the codeword.
///
/// # Returns
/// A vector of error magnitudes for each error location.

pub fn forney_algorithm(
    omega : &PolyGF256,
    sigma : &PolyGF256,
    error_locs : &[u8],
) -> Result<Vec<u8>, String> {

    let sigma_prime = sigma.derivative();

    let mut magnitudes = Vec::new();

    for &loc in error_locs {

        let x_inv = gf256_inv(gf256_pow(
            2,
            u64::from(loc),
        ))?;

        let omega_val = omega.eval(x_inv);

        let sigma_prime_val = sigma_prime.eval(x_inv);

        let magnitude = gf256_div(
            gf256_mul(omega_val, x_inv),
            sigma_prime_val,
        )?;

        magnitudes.push(magnitude);
    }

    Ok(magnitudes)
}

// ============================================================================
// Hamming Codes
// ============================================================================

/// Computes the Hamming distance between two byte slices.
///
/// The Hamming distance is the number of positions at which the corresponding
/// bits differ.
///
/// # Arguments
/// * `a` - First byte slice
/// * `b` - Second byte slice
///
/// # Returns
/// The Hamming distance, or `None` if slices have different lengths.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::hamming_distance_numerical;
///
/// let a = vec![0, 1, 1, 0];
///
/// let b = vec![1, 1, 0, 0];
///
/// assert_eq!(
///     hamming_distance_numerical(&a, &b),
///     Some(2)
/// );
/// ```
#[must_use]

pub fn hamming_distance_numerical(
    a : &[u8],
    b : &[u8],
) -> Option<usize> {

    if a.len() != b.len() {

        return None;
    }

    Some(
        a.iter()
            .zip(b.iter())
            .filter(|(&x, &y)| x != y)
            .count(),
    )
}

/// Computes the Hamming weight (number of 1s) of a byte slice.
///
/// For binary data where each byte is 0 or 1, this counts the number of 1s.
///
/// # Arguments
/// * `data` - Byte slice where each byte represents a bit (0 or 1)
///
/// # Returns
/// The count of non-zero bytes.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::hamming_weight_numerical;
///
/// let data = vec![1, 0, 1, 1, 0, 1];
///
/// assert_eq!(
///     hamming_weight_numerical(&data),
///     4
/// );
/// ```
#[must_use]

pub fn hamming_weight_numerical(data : &[u8]) -> usize {

    data.iter()
        .filter(|&&x| x != 0)
        .count()
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
/// A `Some(Vec<u8>)` of 7 bits representing the codeword, or `None` if the input length is not 4.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::hamming_encode_numerical;
///
/// let data = vec![1, 0, 1, 1];
///
/// let codeword = hamming_encode_numerical(&data).unwrap();
///
/// assert_eq!(codeword.len(), 7);
/// ```
#[must_use]

pub fn hamming_encode_numerical(data : &[u8]) -> Option<Vec<u8>> {

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

/// Decodes a 7-bit Hamming(7,4) codeword, correcting a single-bit error if found.
///
/// # Arguments
/// * `codeword` - A slice of 7 bytes, each representing a bit (0 or 1).
///
/// # Returns
/// A `Result` containing:
/// - `Ok((data, error_pos))` where `data` is the 4-bit corrected data and `error_pos` is
///   `Some(index)` if an error was corrected at the given 1-based index, or `None` if no error was found.
/// - `Err(String)` if the input length is not 7.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::{hamming_decode_numerical, hamming_encode_numerical};
///
/// let data = vec![1, 0, 1, 1];
///
/// let mut codeword = hamming_encode_numerical(&data).unwrap();
///
/// codeword[2] ^= 1; // Introduce error at position 3 (1-indexed)
/// let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();
///
/// assert_eq!(decoded, data);
///
/// assert_eq!(error_pos, Some(3));
/// ```

pub fn hamming_decode_numerical(
    codeword : &[u8]
) -> Result<
    (
        Vec<u8>,
        Option<usize>,
    ),
    String,
> {

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

    Ok((
        corrected_data,
        error_index,
    ))
}

/// Checks if a Hamming(7,4) codeword is valid.
///
/// # Arguments
/// * `codeword` - A 7-bit codeword
///
/// # Returns
/// `true` if the codeword is valid (no errors), `false` otherwise.
#[must_use]

pub fn hamming_check_numerical(codeword : &[u8]) -> bool {

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

// ============================================================================
// BCH Codes
// ============================================================================

/// Encodes data using a BCH code.
///
/// BCH (Bose-Chaudhuri-Hocquenghem) codes are a class of cyclic error-correcting codes.
/// This is a simplified implementation for demonstration.
///
/// # Arguments
/// * `data` - The data bits to encode.
/// * `t` - The error correction capability (can correct up to t errors).
///
/// # Returns
/// The encoded codeword.
#[must_use]

pub fn bch_encode(
    data : &[u8],
    t : usize,
) -> Vec<u8> {

    // Simplified BCH encoding: append parity bits
    let n_parity = 2 * t;

    let mut codeword = data.to_vec();

    // Calculate parity bits using XOR of subsets
    for i in 0 .. n_parity {

        let mut parity = 0u8;

        for (j, &bit) in data
            .iter()
            .enumerate()
        {

            if ((j + 1) >> i) & 1 == 1 {

                parity ^= bit;
            }
        }

        codeword.push(parity);
    }

    codeword
}

/// Decodes a BCH codeword and attempts to correct errors.
///
/// # Arguments
/// * `codeword` - The received codeword.
/// * `t` - The error correction capability.
///
/// # Returns
/// The decoded data, or an error if too many errors are present.

pub fn bch_decode(
    codeword : &[u8],
    t : usize,
) -> Result<Vec<u8>, String> {

    let n_parity = 2 * t;

    if codeword.len() < n_parity {

        return Err("Codeword too short".to_string());
    }

    let data_len = codeword.len() - n_parity;

    let data = &codeword[.. data_len];

    let parity_bits = &codeword[data_len ..];

    // Check parity and find error syndrome
    let mut syndrome = 0usize;

    for i in 0 .. n_parity {

        let mut expected_parity = 0u8;

        for (j, &bit) in data
            .iter()
            .enumerate()
        {

            if ((j + 1) >> i) & 1 == 1 {

                expected_parity ^= bit;
            }
        }

        if expected_parity != parity_bits[i] {

            syndrome |= 1 << i;
        }
    }

    if syndrome == 0 {

        // No errors
        return Ok(data.to_vec());
    }

    // Attempt single error correction
    if syndrome > 0 && syndrome <= data_len {

        let mut corrected = data.to_vec();

        corrected[syndrome - 1] ^= 1;

        return Ok(corrected);
    }

    Err("Unable to correct errors".to_string())
}

// ============================================================================
// CRC Checksums
// ============================================================================

/// CRC-32 polynomial (IEEE 802.3 / ISO 3309 / PKZIP)

const CRC32_POLYNOMIAL : u32 = 0xEDB88320;

/// CRC-16 polynomial (IBM / ANSI)

const CRC16_POLYNOMIAL : u16 = 0xA001;

/// CRC-8 polynomial (ITU)

const CRC8_POLYNOMIAL : u8 = 0x07;

/// Computes the CRC-32 checksum of the given data.
///
/// Uses the standard IEEE 802.3 polynomial (0xEDB88320 in reflected form).
///
/// # Arguments
/// * `data` - The data to compute the checksum for.
///
/// # Returns
/// The 32-bit CRC checksum.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::crc32_compute_numerical;
///
/// let data = b"Hello, World!";
///
/// let checksum = crc32_compute_numerical(data);
///
/// assert_eq!(checksum, 0xEC4AC3D0);
/// ```
#[must_use]

pub fn crc32_compute_numerical(data : &[u8]) -> u32 {

    let mut crc : u32 = 0xFFFFFFFF;

    for byte in data {

        crc ^= u32::from(*byte);

        for _ in 0 .. 8 {

            if crc & 1 != 0 {

                crc = (crc >> 1) ^ CRC32_POLYNOMIAL;
            } else {

                crc >>= 1;
            }
        }
    }

    !crc
}

/// Verifies the CRC-32 checksum of data.
///
/// # Arguments
/// * `data` - The data to verify.
/// * `expected_crc` - The expected CRC-32 checksum.
///
/// # Returns
/// `true` if the computed CRC matches the expected value.
#[must_use]

pub fn crc32_verify_numerical(
    data : &[u8],
    expected_crc : u32,
) -> bool {

    crc32_compute_numerical(data) == expected_crc
}

/// Updates an existing CRC-32 with additional data (for streaming).
///
/// # Arguments
/// * `crc` - Current CRC value (use 0xFFFFFFFF for initial call).
/// * `data` - Additional data to process.
///
/// # Returns
/// Updated CRC value (call `crc32_finalize_numerical` to get final CRC).
#[must_use]

pub fn crc32_update_numerical(
    crc : u32,
    data : &[u8],
) -> u32 {

    let mut crc = crc;

    for byte in data {

        crc ^= u32::from(*byte);

        for _ in 0 .. 8 {

            if crc & 1 != 0 {

                crc = (crc >> 1) ^ CRC32_POLYNOMIAL;
            } else {

                crc >>= 1;
            }
        }
    }

    crc
}

/// Finalizes a CRC-32 computation started with `crc32_update_numerical`.
///
/// # Arguments
/// * `crc` - The running CRC value.
///
/// # Returns
/// The final CRC-32 checksum.
#[must_use]

pub const fn crc32_finalize_numerical(crc : u32) -> u32 {

    !crc
}

/// Computes the CRC-16 checksum of the given data.
///
/// Uses the IBM/ANSI polynomial (0xA001 in reflected form).
///
/// # Arguments
/// * `data` - The data to compute the checksum for.
///
/// # Returns
/// The 16-bit CRC checksum.
///
/// # Example
/// ```
/// 
/// use rssn::numerical::error_correction::crc16_compute;
///
/// let data = b"123456789";
///
/// let checksum = crc16_compute(data);
///
/// assert_eq!(checksum, 0xBB3D);
/// ```
#[must_use]

pub fn crc16_compute(data : &[u8]) -> u16 {

    let mut crc : u16 = 0xFFFF;

    for byte in data {

        crc ^= u16::from(*byte);

        for _ in 0 .. 8 {

            if crc & 1 != 0 {

                crc = (crc >> 1) ^ CRC16_POLYNOMIAL;
            } else {

                crc >>= 1;
            }
        }
    }

    crc
}

/// Computes the CRC-8 checksum of the given data.
///
/// Uses the ITU polynomial (0x07).
///
/// # Arguments
/// * `data` - The data to compute the checksum for.
///
/// # Returns
/// The 8-bit CRC checksum.
#[must_use]

pub fn crc8_compute(data : &[u8]) -> u8 {

    let mut crc : u8 = 0;

    for byte in data {

        crc ^= *byte;

        for _ in 0 .. 8 {

            if crc & 0x80 != 0 {

                crc = (crc << 1) ^ CRC8_POLYNOMIAL;
            } else {

                crc <<= 1;
            }
        }
    }

    crc
}

// ============================================================================
// Interleaving for Burst Error Correction
// ============================================================================

/// Interleaves data to spread burst errors across multiple codewords.
///
/// This technique helps in correcting burst errors by distributing consecutive
/// symbols across different codewords.
///
/// # Arguments
/// * `data` - The data to interleave.
/// * `depth` - The interleaving depth (number of codewords to spread across).
///
/// # Returns
/// The interleaved data.
#[must_use]

pub fn interleave(
    data : &[u8],
    depth : usize,
) -> Vec<u8> {

    if depth == 0 || data.is_empty() {

        return data.to_vec();
    }

    let n = data.len();

    let rows = n.div_ceil(depth);

    let mut result = Vec::with_capacity(n);

    for col in 0 .. depth {

        for row in 0 .. rows {

            let idx = row * depth + col;

            if idx < n {

                result.push(data[idx]);
            }
        }
    }

    result
}

/// De-interleaves data (reverses the interleave operation).
///
/// # Arguments
/// * `data` - The interleaved data.
/// * `depth` - The interleaving depth used during interleaving.
///
/// # Returns
/// The de-interleaved data.
#[must_use]

pub fn deinterleave(
    data : &[u8],
    depth : usize,
) -> Vec<u8> {

    if depth == 0 || data.is_empty() {

        return data.to_vec();
    }

    let n = data.len();

    let rows = n.div_ceil(depth);

    let full_cols = n % depth;

    let full_cols = if full_cols == 0 {

        depth
    } else {

        full_cols
    };

    let mut result = vec![0u8; n];

    let mut idx = 0;

    for col in 0 .. depth {

        let col_len = if col < full_cols {

            rows
        } else {

            rows - 1
        };

        for row in 0 .. col_len {

            let orig_idx = row * depth + col;

            if orig_idx < n {

                result[orig_idx] = data[idx];

                idx += 1;
            }
        }
    }

    result
}

// ============================================================================
// Convolutional Codes (Simplified)
// ============================================================================

/// Encodes data using a simple rate-1/2 convolutional code.
///
/// Uses generators G1 = 1+D^2 (0b101) and G2 = 1+D+D^2 (0b111).
///
/// # Arguments
/// * `data` - The binary data to encode (each byte is 0 or 1).
///
/// # Returns
/// The encoded data (twice the length of input).
#[must_use]

pub fn convolutional_encode(data : &[u8]) -> Vec<u8> {

    let mut state : u8 = 0;

    let mut output = Vec::with_capacity(data.len() * 2);

    for &bit in data {

        state = (state >> 1) | ((bit & 1) << 2);

        // G1 = 0b101: bits 0 and 2
        let g1 = (state & 1) ^ ((state >> 2) & 1);

        // G2 = 0b111: bits 0, 1, and 2
        let g2 = (state & 1) ^ ((state >> 1) & 1) ^ ((state >> 2) & 1);

        output.push(g1);

        output.push(g2);
    }

    // Flush encoder state
    for _ in 0 .. 2 {

        state >>= 1;

        let g1 = (state & 1) ^ ((state >> 2) & 1);

        let g2 = (state & 1) ^ ((state >> 1) & 1) ^ ((state >> 2) & 1);

        output.push(g1);

        output.push(g2);
    }

    output
}

/// Computes the minimum distance of a code given a set of codewords.
///
/// # Arguments
/// * `codewords` - A slice of codewords.
///
/// # Returns
/// The minimum Hamming distance between any two distinct codewords.
#[must_use]

pub fn minimum_distance(codewords : &[Vec<u8>]) -> Option<usize> {

    if codewords.len() < 2 {

        return None;
    }

    let mut min_dist : Option<usize> = None;

    for i in 0 .. codewords.len() {

        for j in (i + 1) .. codewords.len() {

            if let Some(dist) = hamming_distance_numerical(
                &codewords[i],
                &codewords[j],
            ) {

                min_dist = Some(
                    min_dist.map_or(dist, |m| {
                        m.min(dist)
                    }),
                );
            }
        }
    }

    min_dist
}

/// Computes the code rate (number of information bits per coded bit).
///
/// # Arguments
/// * `k` - Number of information bits.
/// * `n` - Total codeword length.
///
/// # Returns
/// The code rate as a floating-point number.
#[must_use]

pub fn code_rate(
    k : usize,
    n : usize,
) -> f64 {

    if n == 0 {

        0.0
    } else {

        k as f64 / n as f64
    }
}

/// Estimates the error correction capability from the minimum distance.
///
/// A code with minimum distance d can correct up to floor((d-1)/2) errors.
///
/// # Arguments
/// * `min_distance` - The minimum distance of the code.
///
/// # Returns
/// The number of errors that can be corrected.
#[must_use]

pub const fn error_correction_capability(min_distance : usize) -> usize {

    if min_distance == 0 {

        0
    } else {

        (min_distance - 1) / 2
    }
}

/// Estimates the error detection capability from the minimum distance.
///
/// A code with minimum distance d can detect up to d-1 errors.
///
/// # Arguments
/// * `min_distance` - The minimum distance of the code.
///
/// # Returns
/// The number of errors that can be detected.
#[must_use]

pub const fn error_detection_capability(min_distance : usize) -> usize {

    if min_distance == 0 {

        0
    } else {

        min_distance - 1
    }
}
