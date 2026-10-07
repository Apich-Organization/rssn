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
//! - `Bch` - narrow-sense primitive binary BCH code over GF(2^m): generator
//!   polynomial from minimal polynomials, systematic encoding, decoding with
//!   syndromes, Berlekamp-Massey and Chien search
//! - `bch_encode` - shortened BCH encoding of a bit vector (wrapper over `Bch`)
//! - `bch_decode` - matching decoder with error correction
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
//! use rssn::kernels::error_correction::reed_solomon_decode;
//! use rssn::kernels::error_correction::reed_solomon_encode;
//!
//! let message = vec![0x01, 0x02, 0x03, 0x04];
//!
//! let codeword = reed_solomon_encode(&message, 4).unwrap();
//!
//! // Introduce an error
//! let mut corrupted = codeword.clone();
//!
//! corrupted[0] ^= 0xFF;
//!
//! // Decode and correct
//! reed_solomon_decode(&mut corrupted, 4).unwrap();
//!
//! assert_eq!(&corrupted[..4], &message);
//! ```
//!
//! ### CRC-32 Checksum
//! ```
//! use rssn::kernels::error_correction::crc32_compute_numerical;
//! use rssn::kernels::error_correction::crc32_verify_numerical;
//!
//! let data = b"Hello, World!";
//!
//! let checksum = crc32_compute_numerical(data);
//!
//! assert!(crc32_verify_numerical(data, checksum));
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::kernels::finite_field::gf256_add;
use crate::kernels::finite_field::gf256_div;
use crate::kernels::finite_field::gf256_inv;
use crate::kernels::finite_field::gf256_mul;
use crate::kernels::finite_field::gf256_pow;

// ============================================================================
// Polynomial over GF(2^8)
// ============================================================================

/// Represents a polynomial over GF(2^8).
///
/// The polynomial is stored in ascending order of powers: `coeffs[i]` is the
/// coefficient of `x^i`, so the first element is the constant term. Every
/// method of this type (`eval`, `poly_add`, `poly_mul`, `poly_div`,
/// `derivative`, `normalize`) uses this one order.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolyGF256(pub Vec<u8>);

impl PolyGF256 {
    /// Creates a new polynomial from coefficients.
    ///
    /// # Arguments
    /// * `coeffs` - Coefficients in ascending order of powers (`coeffs[i]` is the
    ///   coefficient of `x^i`).
    #[must_use]
    pub const fn new(coeffs: Vec<u8>) -> Self {
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
        x: u8,
    ) -> u8 {
        self.0
            .iter()
            .rfold(0, |acc, &coeff| gf256_add(gf256_mul(acc, x), coeff))
    }

    /// Adds two polynomials over GF(2^8).
    ///
    /// In GF(2^8), addition is XOR.
    #[must_use]
    pub fn poly_add(
        &self,
        other: &Self,
    ) -> Self {
        let mut result = vec![0; self.0.len().max(other.0.len())];

        let result_len = result.len();

        // Ascending order: the constant terms line up at index 0.
        result[..self.0.len()].copy_from_slice(&self.0[..]);

        for (r, o) in result.iter_mut().zip(&other.0) {
            *r ^= *o;
        }

        debug_assert_eq!(result.len(), result_len);

        Self(result)
    }

    /// Subtracts two polynomials over GF(2^8).
    ///
    /// In GF(2^8), subtraction is the same as addition (XOR).
    #[must_use]
    pub fn poly_sub(
        &self,
        other: &Self,
    ) -> Self {
        self.poly_add(other)
    }

    /// Multiplies two polynomials over GF(2^8).
    #[must_use]
    pub fn poly_mul(
        &self,
        other: &Self,
    ) -> Self {
        let mut result = vec![0; self.degree() + other.degree() + 1];

        for i in 0..=self.degree() {
            for j in 0..=other.degree() {
                result[i + j] ^= gf256_mul(self.0[i], other.0[j]);
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
    ///
    /// # Errors
    /// Returns an error if the divisor polynomial is empty (zero polynomial).
    pub fn poly_div(
        &self,
        divisor: &Self,
    ) -> Result<(Self, Self), String> {
        if divisor.0.is_empty() {
            return Err("Division by zero \
                 polynomial"
                .to_string());
        }

        // Ignore zero high-order coefficients so the true leading term divides.
        let d = divisor.normalize();

        if d.0.is_empty() {
            return Err("Division by zero \
                 polynomial"
                .to_string());
        }

        let d_len = d.0.len();

        let lead_inv = gf256_inv(d.0[d_len - 1])?;

        let mut rem = self.0.clone();

        if rem.len() < d_len {
            return Ok((Self(vec![0]), Self(rem)));
        }

        let mut quot = vec![0; rem.len() - d_len + 1];

        // Eliminate the highest power first; quotient coefficient k belongs to x^k.
        for k in (0..quot.len()).rev() {
            let q_coeff = gf256_mul(rem[k + d_len - 1], lead_inv);

            if q_coeff == 0 {
                continue;
            }

            quot[k] = q_coeff;

            for (i, &dc) in d.0.iter().enumerate() {
                rem[k + i] ^= gf256_mul(dc, q_coeff);
            }
        }

        rem.truncate(d_len - 1);

        Ok((Self(quot), Self(rem)))
    }

    /// Returns the derivative of the polynomial over GF(2^8).
    ///
    /// In GF(2^8), the derivative has a special property:
    /// d/dx(x^n) = n * x^(n-1), where n is reduced mod 2.
    /// Thus, only odd-power terms contribute to the derivative.
    #[must_use]
    pub fn derivative(&self) -> Self {
        let mut deriv = vec![0; self.degree()];

        for i in 1..=self.degree() {
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
        c: u8,
    ) -> Self {
        Self(self.0.iter().map(|&coeff| gf256_mul(coeff, c)).collect())
    }

    /// Normalizes the polynomial by removing zero coefficients of the highest
    /// powers (the trailing elements in ascending order).
    #[must_use]
    pub fn normalize(&self) -> Self {
        let end = self.0.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);

        Self(self.0[..end].to_vec())
    }
}

// ============================================================================
// Reed-Solomon Codes
// ============================================================================

/// Computes the generator polynomial for a Reed-Solomon code with `n_parity` error correction symbols.
///
/// The generator polynomial is the product of (x - α^i) for i = 0 to n_parity-1,
/// where α is the primitive element (2) of GF(2^8).
pub(crate) fn rs_generator_poly(n_parity: usize) -> Result<Vec<u8>, String> {
    if n_parity == 0 {
        return Err("Number of parity symbols \
             must be positive"
            .to_string());
    }

    let mut g = vec![1u8];

    for i in 0..n_parity {
        // Multiply by (x - α^i), which in GF(2^8) is (x + α^i) since subtraction = addition
        let root = gf256_pow(2, i as u64);

        let factor = vec![1u8, root];

        g = poly_mul_gf256(&g, &factor);
    }

    Ok(g)
}

/// Multiplies two polynomials over GF(2^8).
fn poly_mul_gf256(
    p1: &[u8],
    p2: &[u8],
) -> Vec<u8> {
    if p1.is_empty() || p2.is_empty() {
        return vec![];
    }

    let mut result = vec![0; p1.len() + p2.len() - 1];

    for i in 0..p1.len() {
        for j in 0..p2.len() {
            result[i + j] ^= gf256_mul(p1[i], p2[j]);
        }
    }

    result
}

/// Performs polynomial long division over GF(2^8).
/// Returns the remainder after dividing dividend by divisor.
fn poly_div_gf256(
    mut dividend: Vec<u8>,
    divisor: &[u8],
) -> Result<Vec<u8>, String> {
    if divisor.is_empty() {
        return Err("Divisor cannot \
                    be empty"
            .to_string());
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

        let coeff = gf256_mul(lead_dividend, lead_divisor_inv);

        for i in 0..divisor_len {
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
    poly: &[u8],
    x: u8,
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
/// # Errors
/// Returns an error if:
/// - The total length (message + parity) exceeds 255.
/// - The number of parity symbols is invalid.
///
/// # Example
/// ```
/// use rssn::kernels::error_correction::reed_solomon_encode;
///
/// let message = vec![0x01, 0x02, 0x03];
///
/// let codeword = reed_solomon_encode(&message, 4).unwrap();
///
/// assert_eq!(codeword.len(), 7); // 3 data + 4 parity
/// ```
pub fn reed_solomon_encode(
    message: &[u8],
    n_parity: usize,
) -> Result<Vec<u8>, String> {
    if message.len() + n_parity > 255 {
        return Err("Message + parity length \
             cannot exceed 255"
            .to_string());
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
    let remainder = poly_div_gf256(message_poly, &gen_poly)?;

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
/// # Errors
/// Returns an error if:
/// - The number of errors exceeds the correction capability of the code.
/// - No valid error locations can be found.
/// - The Berlekamp-Massey or Forney algorithms fail to converge.
///
/// # Example
/// ```
/// use rssn::kernels::error_correction::reed_solomon_decode;
/// use rssn::kernels::error_correction::reed_solomon_encode;
///
/// let message = vec![0x01, 0x02, 0x03];
///
/// let mut codeword = reed_solomon_encode(&message, 4).unwrap();
///
/// codeword[0] ^= 0xFF; // Introduce error
/// reed_solomon_decode(&mut codeword, 4).unwrap();
///
/// assert_eq!(&codeword[..3], &message);
/// ```
pub fn reed_solomon_decode(
    codeword: &mut [u8],
    n_parity: usize,
) -> Result<(), String> {
    let syndromes = calculate_syndromes(codeword, n_parity);

    if syndromes.iter().all(|&s| s == 0) {
        return Ok(());
    }

    // Use Berlekamp-Massey to find error locator polynomial
    let sigma = berlekamp_massey(&syndromes);

    // More than t = floor(n_parity / 2) errors cannot be located uniquely;
    // report it instead of "correcting" to a different codeword.
    let sigma_degree = sigma.iter().rposition(|&x| x != 0).unwrap_or(0);

    if sigma_degree > n_parity / 2 {
        return Err("Uncorrectable errors: more than \
                    n_parity / 2 errors."
            .to_string());
    }

    // Use Chien search to find error locations
    let error_locations = chien_search_extended(&sigma, codeword.len())?;

    if error_locations.is_empty() {
        return Err("Failed to find \
                    error locations."
            .to_string());
    }

    // Error evaluator omega(x) = S(x) sigma(x) mod x^n_parity; both are stored
    // lowest power first, so the low-order coefficients are kept.
    let mut omega = poly_mul_gf256(&syndromes, &sigma);

    omega.truncate(n_parity);

    // Use Forney's algorithm to find error magnitudes
    let error_magnitudes =
        forney_algorithm_extended(&omega, &sigma, &error_locations, codeword.len())?;

    // Correct a copy first so that a failed verification leaves the input intact.
    let mut corrected = codeword.to_vec();

    for (i, &loc) in error_locations.iter().enumerate() {
        corrected[loc] ^= error_magnitudes[i];
    }

    if !calculate_syndromes(&corrected, n_parity)
        .iter()
        .all(|&s| s == 0)
    {
        return Err("Uncorrectable errors: \
                    corrected word is not a codeword."
            .to_string());
    }

    codeword.copy_from_slice(&corrected);

    Ok(())
}

/// Berlekamp-Massey algorithm to find the error locator polynomial.
pub(crate) fn berlekamp_massey(syndromes: &[u8]) -> Vec<u8> {
    let mut sigma = vec![1u8];

    let mut b = vec![1u8];

    let mut l = 0;

    for i in 0..syndromes.len() {
        let mut delta = syndromes[i];

        for j in 1..=l {
            if j < sigma.len() {
                delta ^= gf256_mul(sigma[j], syndromes[i - j]);
            }
        }

        b.insert(0, 0);

        // Correctly represents b(x) * x in Little-Endian

        if delta != 0 {
            let t = sigma.clone();

            let scaled_b: Vec<u8> = b.iter().map(|&c| gf256_mul(c, delta)).collect();

            sigma = poly_add_gf256(&sigma, &scaled_b);

            if 2 * l <= i {
                l = i + 1 - l;

                let inv_delta = gf256_inv(delta).unwrap();

                b = t.iter().map(|&c| gf256_mul(c, inv_delta)).collect();
            }
        }
    }

    sigma
}

/// Little-Endian Evaluation (Horner's Method)
fn poly_eval_gf256_le(
    poly: &[u8],
    x: u8,
) -> u8 {
    let mut result = 0;

    // Iterate from the highest power down to the constant term
    // result = (...((an*x + an-1)*x + an-2)...)*x + a0
    for &coeff in poly.iter().rev() {
        result = gf256_mul(result, x) ^ coeff;
    }

    result
}

/// Add two polynomials over GF(2^8).
fn poly_add_gf256(
    p1: &[u8],
    p2: &[u8],
) -> Vec<u8> {
    let mut result = vec![0u8; p1.len().max(p2.len())];

    for (i, &c) in p1.iter().enumerate() {
        result[i] ^= c;
    }

    for (i, &c) in p2.iter().enumerate() {
        result[i] ^= c;
    }

    result
}

/// Extended Chien search to find error locations.
fn chien_search_extended(
    sigma: &[u8],
    codeword_len: usize,
) -> Result<Vec<usize>, String> {
    let mut error_locs = Vec::new();

    // Find the actual degree of the polynomial (last non-zero coefficient)
    let degree = sigma.iter().rposition(|&x| x != 0).unwrap_or(0);

    // If no errors (sigma is just [1]), we are done
    if degree == 0 {
        return Ok(Vec::new());
    }

    for i in 0..codeword_len {
        // Map array index i to the power of x: (n - 1 - i)
        let power = codeword_len - 1 - i;

        // Root we are looking for is alpha^(-power)
        let x_inv = gf256_pow(2, power as u64);

        let x = gf256_inv(x_inv).unwrap_or(0);

        if poly_eval_gf256_le(sigma, x) == 0 {
            error_locs.push(i);
        }
    }

    // CRITICAL: Validation
    // A polynomial of degree L must have exactly L roots in the field
    // to be a valid error locator.
    if error_locs.len() != degree {
        return Err(format!(
            "Uncorrectable errors: \
             found {} roots, but \
             expected {}",
            error_locs.len(),
            degree
        ));
    }

    Ok(error_locs)
}

/// Extended Forney's algorithm to compute error magnitudes.
fn forney_algorithm_extended(
    omega: &[u8],
    sigma: &[u8],
    error_locs: &[usize],
    codeword_len: usize,
) -> Result<Vec<u8>, String> {
    let sigma_prime = poly_derivative_gf256(sigma);

    let mut magnitudes = Vec::new();

    for &loc in error_locs {
        // 1. Calculate the same 'x' (root) used in Chien Search
        let power = codeword_len - 1 - loc;

        let x_inv = gf256_pow(2, power as u64); // alpha^j
        let x = gf256_inv(x_inv).unwrap_or(0); // alpha^-j (the root)

        // 2. Evaluate using Little-Endian logic
        let omega_val = poly_eval_gf256_le(omega, x);

        let sigma_prime_val = poly_eval_gf256_le(&sigma_prime, x);

        if sigma_prime_val == 0 {
            return Err("Division by zero in \
                 Forney"
                .to_string());
        }

        // 3. With roots alpha^0..alpha^(n-1) (first consecutive root b = 0)
        // Forney's formula is Y_i = X_i^(1-b) Omega(X_i^-1) / Sigma'(X_i^-1),
        // i.e. an extra factor X_i = alpha^power.
        let magnitude = gf256_mul(x_inv, gf256_div(omega_val, sigma_prime_val)?);

        magnitudes.push(magnitude);
    }

    Ok(magnitudes)
}

/// Compute formal derivative of polynomial in GF(2^8).
fn poly_derivative_gf256(poly: &[u8]) -> Vec<u8> {
    let mut result = Vec::new();

    // derivative of sum(a_i * x^i) is sum(i * a_i * x^{i-1})
    // In GF(2^8), i * a_i is 0 if i is even, a_i if i is odd.
    for (i, &coeff) in poly.iter().enumerate().skip(1) {
        if i % 2 == 1 {
            result.push(coeff);
        } else {
            result.push(0);
        }
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
    codeword: &[u8],
    n_parity: usize,
) -> bool {
    let syndromes = calculate_syndromes(codeword, n_parity);

    syndromes.iter().all(|&s| s == 0)
}

/// Calculates the syndromes of a received codeword.
///
/// The syndrome `S_i` is computed as the codeword polynomial evaluated at α^i.
#[must_use]
pub fn calculate_syndromes(
    codeword: &[u8],
    n_parity: usize,
) -> Vec<u8> {
    let mut syndromes = Vec::with_capacity(n_parity);

    for i in 0..n_parity {
        let alpha_i = gf256_pow(2, i as u64);

        syndromes.push(poly_eval_gf256(codeword, alpha_i));
    }

    syndromes
}

/// Finds the roots of the error locator polynomial to determine error locations.
///
/// Uses Chien search to efficiently evaluate the polynomial at all field elements.
///
/// # Errors
/// Returns an error if the field element inversion fails.
pub fn chien_search(sigma: &PolyGF256) -> Result<Vec<u8>, String> {
    let mut error_locs = Vec::new();

    for i in 0..255u8 {
        let alpha_inv = gf256_inv(gf256_pow(2, u64::from(i)))?;

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
///
/// # Errors
/// Returns an error if division by zero occurs during magnitude calculation.
pub fn forney_algorithm(
    omega: &PolyGF256,
    sigma: &PolyGF256,
    error_locs: &[u8],
) -> Result<Vec<u8>, String> {
    let sigma_prime = sigma.derivative();

    let mut magnitudes = Vec::new();

    for &loc in error_locs {
        let x_inv = gf256_inv(gf256_pow(2, u64::from(loc)))?;

        let omega_val = omega.eval(x_inv);

        let sigma_prime_val = sigma_prime.eval(x_inv);

        let magnitude = gf256_div(gf256_mul(omega_val, x_inv), sigma_prime_val)?;

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
/// use rssn::kernels::error_correction::hamming_distance_numerical;
///
/// let a = vec![0, 1, 1, 0];
///
/// let b = vec![1, 1, 0, 0];
///
/// assert_eq!(hamming_distance_numerical(&a, &b), Some(2));
/// ```
#[must_use]
pub fn hamming_distance_numerical(
    a: &[u8],
    b: &[u8],
) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    Some(a.iter().zip(b.iter()).filter(|&(&x, &y)| x != y).count())
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
/// use rssn::kernels::error_correction::hamming_weight_numerical;
///
/// let data = vec![1, 0, 1, 1, 0, 1];
///
/// assert_eq!(hamming_weight_numerical(&data), 4);
/// ```
#[must_use]
pub fn hamming_weight_numerical(data: &[u8]) -> usize {
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
/// A `Some(Vec<u8>)` of 7 bits representing the codeword, or `None` if the input length is not 4.
///
/// # Example
/// ```
/// use rssn::kernels::error_correction::hamming_encode_numerical;
///
/// let data = vec![1, 0, 1, 1];
///
/// let codeword = hamming_encode_numerical(&data).unwrap();
///
/// assert_eq!(codeword.len(), 7);
/// ```
#[must_use]
pub fn hamming_encode_numerical(data: &[u8]) -> Option<Vec<u8>> {
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
/// # Arguments
/// * `codeword` - A slice of 7 bytes, each representing a bit (0 or 1).
///
/// # Returns
/// A `Result` containing:
/// - `Ok((data, error_pos))` where `data` is the 4-bit corrected data and `error_pos` is
///   `Some(index)` if an error was corrected at the given 1-based index, or `None` if no error was found.
/// - `Err(String)` if the input length is not 7.
///
/// # Errors
/// Returns an error if the input codeword length is not exactly 7.
///
/// # Example
/// ```
/// use rssn::kernels::error_correction::hamming_decode_numerical;
/// use rssn::kernels::error_correction::hamming_encode_numerical;
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
pub fn hamming_decode_numerical(codeword: &[u8]) -> Result<(Vec<u8>, Option<usize>), String> {
    if codeword.len() != 7 {
        return Err("Codeword length \
                    must be 7"
            .to_string());
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

/// Checks if a Hamming(7,4) codeword is valid.
///
/// # Arguments
/// * `codeword` - A 7-bit codeword
///
/// # Returns
/// `true` if the codeword is valid (no errors), `false` otherwise.
#[must_use]
pub const fn hamming_check_numerical(codeword: &[u8]) -> bool {
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

/// Primitive polynomials of GF(2^m) (bit `i` is the coefficient of `x^i`),
/// indexed by `m`.
const BCH_PRIMITIVE_POLYS: [u32; 17] = [
    0, 0, 0x7, 0xB, 0x13, 0x25, 0x43, 0x89, 0x11D, 0x211, 0x409, 0x805, 0x1053, 0x201B, 0x4443,
    0x8003, 0x1100B,
];

/// Largest field degree `m` supported by [`Bch`] (block length `2^m - 1`).
pub const BCH_MAX_M: usize = 16;

/// A narrow-sense primitive binary BCH code of length `n = 2^m - 1`
/// correcting `t` errors.
///
/// The generator polynomial is the least common multiple of the minimal
/// polynomials of `alpha, alpha^2, ..., alpha^(2t)` over GF(2), where
/// `alpha` is a root of a primitive polynomial of degree `m`; the code
/// dimension is `k = n - deg g`. Encoding is systematic (the message bits
/// come first, followed by the `n - k` parity bits) and decoding computes
/// the syndromes `S_j = r(alpha^j)`, solves the key equation with the
/// Berlekamp-Massey algorithm and locates the errors with a Chien search.
///
/// Bits are `u8` values `0` or `1`; the first bit of a vector is the
/// coefficient of the highest power of `x`. Shorter messages are
/// supported by *shortening*: a message of `l <= k` bits is encoded as if
/// padded with `k - l` leading zeros, which are then not transmitted.
///
/// # Example
/// ```
/// use rssn::kernels::error_correction::Bch;
///
/// let code = Bch::new(4, 2).unwrap(); // the (15, 7) double-error-correcting code
/// assert_eq!((code.n(), code.k()), (15, 7));
/// let msg = [1, 0, 1, 1, 0, 0, 1];
/// let mut cw = code.encode(&msg).unwrap();
/// cw[2] ^= 1;
/// cw[11] ^= 1;
/// assert_eq!(code.decode(&cw).unwrap(), msg);
/// ```
#[derive(Debug, Clone)]
pub struct Bch {
    m: usize,
    n: usize,
    t: usize,
    exp: Vec<u32>,
    log: Vec<u32>,
    /// Generator polynomial over GF(2), `gen[i]` the coefficient of `x^i`.
    genpoly: Vec<u8>,
}

impl Bch {
    /// Builds the BCH code of length `2^m - 1` correcting `t` errors.
    ///
    /// # Errors
    /// Fails if `m` is outside `2..=16`, `t == 0`, or `t` is so large that
    /// the code has no information bits (`k < 1`).
    pub fn new(
        m: usize,
        t: usize,
    ) -> Result<Self, String> {
        if !(2..=BCH_MAX_M).contains(&m) {
            return Err(format!("BCH field degree m must be in 2..={BCH_MAX_M}, got {m}"));
        }

        if t == 0 {
            return Err("BCH error capability t must be at least 1".to_string());
        }

        let n = (1usize << m) - 1;

        if 2 * t >= n {
            return Err(format!("t = {t} too large for length {n}"));
        }

        let poly = BCH_PRIMITIVE_POLYS[m];

        let mut exp = vec![0u32; 2 * n];

        let mut log = vec![0u32; n + 1];

        let mut x = 1u32;

        for (i, slot) in exp.iter_mut().enumerate().take(n) {
            *slot = x;

            log[x as usize] = i as u32;

            x <<= 1;

            if x >> m != 0 {
                x ^= poly;
            }
        }

        for i in n..2 * n {
            exp[i] = exp[i - n];
        }

        let mut code = Self { m, n, t, exp, log, genpoly: Vec::new() };

        code.genpoly = code.build_generator()?;

        if code.genpoly.len() > n {
            return Err(format!("t = {t} too large for length {n}: no information bits"));
        }

        Ok(code)
    }

    /// Builds the smallest-length code of degree `m <= 16` correcting `t`
    /// errors that carries at least `data_bits` information bits.
    ///
    /// # Errors
    /// Fails if no supported code is large enough.
    pub fn for_data_len(
        data_bits: usize,
        t: usize,
    ) -> Result<Self, String> {
        for m in 2..=BCH_MAX_M {
            if let Ok(code) = Self::new(m, t)
                && code.k() >= data_bits
            {
                return Ok(code);
            }
        }

        Err(format!("no BCH code with t = {t} carries {data_bits} data bits"))
    }

    /// Field degree `m` (the block length is `2^m - 1`).
    #[must_use]
    pub const fn m(&self) -> usize {
        self.m
    }

    /// Block length `n = 2^m - 1`.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Dimension `k = n - deg g` (message bits per block).
    #[must_use]
    pub const fn k(&self) -> usize {
        self.n + 1 - self.genpoly.len()
    }

    /// Designed error-correction capability `t`.
    #[must_use]
    pub const fn t(&self) -> usize {
        self.t
    }

    /// Number of parity bits `n - k`.
    #[must_use]
    pub const fn parity_len(&self) -> usize {
        self.genpoly.len() - 1
    }

    /// Generator polynomial coefficients over GF(2), lowest power first.
    #[must_use]
    pub fn generator(&self) -> &[u8] {
        &self.genpoly
    }

    fn mul(
        &self,
        a: u32,
        b: u32,
    ) -> u32 {
        if a == 0 || b == 0 {
            0
        } else {
            self.exp[(self.log[a as usize] + self.log[b as usize]) as usize]
        }
    }

    fn inv(
        &self,
        a: u32,
    ) -> u32 {
        self.exp[self.n - self.log[a as usize] as usize]
    }

    fn build_generator(&self) -> Result<Vec<u8>, String> {
        let n = self.n;

        let mut covered = vec![false; n];

        // g(x) over GF(2^m), lowest power first.
        let mut g: Vec<u32> = vec![1];

        for i in 1..=2 * self.t {
            if covered[i % n] {
                continue;
            }

            // minimal polynomial of alpha^i: prod over the cyclotomic coset.
            let mut minpoly: Vec<u32> = vec![1];

            let mut j = i % n;

            loop {
                covered[j] = true;

                let root = self.exp[j];

                let mut next = vec![0u32; minpoly.len() + 1];

                for (idx, &c) in minpoly.iter().enumerate() {
                    next[idx + 1] ^= c;

                    next[idx] ^= self.mul(c, root);
                }

                minpoly = next;

                j = (j * 2) % n;

                if j == i % n {
                    break;
                }
            }

            let mut prod = vec![0u32; g.len() + minpoly.len() - 1];

            for (a, &ga) in g.iter().enumerate() {
                for (b, &mb) in minpoly.iter().enumerate() {
                    prod[a + b] ^= self.mul(ga, mb);
                }
            }

            g = prod;
        }

        g.iter()
            .map(|&c| match c {
                0 => Ok(0u8),
                1 => Ok(1u8),
                _ => Err("generator polynomial is not binary".to_string()),
            })
            .collect()
    }

    /// Encodes `data` (`1..=k` bits) into a systematic codeword
    /// `data ++ parity` of `data.len() + (n - k)` bits.
    ///
    /// # Errors
    /// Fails for an empty message, a message longer than `k`, or entries
    /// other than `0` and `1`.
    pub fn encode(
        &self,
        data: &[u8],
    ) -> Result<Vec<u8>, String> {
        if data.is_empty() {
            return Err("empty message".to_string());
        }

        if data.len() > self.k() {
            return Err(format!("message of {} bits exceeds k = {}", data.len(), self.k()));
        }

        if data.iter().any(|&b| b > 1) {
            return Err("message entries must be 0 or 1".to_string());
        }

        let p = self.parity_len();

        // Remainder of x^p m(x) modulo g(x) by bitwise long division.
        let mut rem = vec![0u8; p];

        for &bit in data {
            let fb = bit ^ rem[p - 1];

            for i in (1..p).rev() {
                rem[i] = rem[i - 1] ^ (fb & self.genpoly[i]);
            }

            rem[0] = fb & self.genpoly[0];
        }

        let mut cw = data.to_vec();

        cw.extend(rem.iter().rev());

        Ok(cw)
    }

    /// The syndromes `S_1..S_2t` of a received word (shortened words are
    /// treated as zero-padded to length `n`).
    fn syndromes(
        &self,
        word: &[u8],
    ) -> Vec<u32> {
        let len = word.len();

        (1..=2 * self.t)
            .map(|j| {
                let mut s = 0u32;

                for (pos, &bit) in word.iter().enumerate() {
                    if bit == 1 {
                        let e = len - 1 - pos;

                        s ^= self.exp[(j * e) % self.n];
                    }
                }

                s
            })
            .collect()
    }

    /// Whether `word` is a codeword (all syndromes vanish).
    #[must_use]
    pub fn is_codeword(
        &self,
        word: &[u8],
    ) -> bool {
        word.len() <= self.n && word.iter().all(|&b| b <= 1) && self.syndromes(word).iter().all(|&s| s == 0)
    }

    /// Berlekamp-Massey: the error locator polynomial `Lambda(x)` (lowest
    /// power first, `Lambda(0) = 1`) from the syndromes.
    fn berlekamp_massey(
        &self,
        s: &[u32],
    ) -> Vec<u32> {
        let mut c: Vec<u32> = vec![1];

        let mut b: Vec<u32> = vec![1];

        let mut l = 0usize;

        let mut shift = 1usize;

        let mut bd = 1u32;

        for nn in 0..s.len() {
            let mut d = s[nn];

            for i in 1..=l.min(c.len() - 1) {
                d ^= self.mul(c[i], s[nn - i]);
            }

            if d == 0 {
                shift += 1;

                continue;
            }

            let coef = self.mul(d, self.inv(bd));

            let old = c.clone();

            if c.len() < b.len() + shift {
                c.resize(b.len() + shift, 0);
            }

            for (i, &bi) in b.iter().enumerate() {
                c[i + shift] ^= self.mul(coef, bi);
            }

            if 2 * l <= nn {
                l = nn + 1 - l;

                b = old;

                bd = d;

                shift = 1;
            } else {
                shift += 1;
            }
        }

        c.truncate(l + 1);

        c
    }

    /// Corrects up to `t` bit errors in a (possibly shortened) codeword.
    ///
    /// Returns the corrected codeword and the number of corrected bits.
    ///
    /// # Errors
    /// Fails if the word is empty, longer than `n`, not binary, or if the
    /// error pattern is uncorrectable (more than `t` errors that are
    /// detected; patterns of more than `t` errors may also be mis-decoded
    /// to a different codeword, as for every bounded-distance decoder).
    pub fn correct(
        &self,
        word: &[u8],
    ) -> Result<(Vec<u8>, usize), String> {
        if word.len() <= self.parity_len() || word.len() > self.n {
            return Err(format!(
                "codeword length {} outside {}..={}",
                word.len(),
                self.parity_len() + 1,
                self.n
            ));
        }

        if word.iter().any(|&b| b > 1) {
            return Err("codeword entries must be 0 or 1".to_string());
        }

        let s = self.syndromes(word);

        if s.iter().all(|&x| x == 0) {
            return Ok((word.to_vec(), 0));
        }

        let lambda = self.berlekamp_massey(&s);

        let deg = lambda.len() - 1;

        if deg == 0 || deg > self.t {
            return Err("Unable to correct errors".to_string());
        }

        // Chien search over the exponents e = 0..len: an error at exponent
        // e has locator alpha^e, a root of Lambda at alpha^-e.
        let len = word.len();

        let mut corrected = word.to_vec();

        let mut found = 0usize;

        for e in 0..self.n {
            let x = self.exp[(self.n - e) % self.n];

            let mut acc = 0u32;

            let mut xp = 1u32;

            for &c in &lambda {
                acc ^= self.mul(c, xp);

                xp = self.mul(xp, x);
            }

            if acc == 0 {
                if e >= len {
                    return Err("Unable to correct errors".to_string());
                }

                corrected[len - 1 - e] ^= 1;

                found += 1;
            }
        }

        if found != deg || !self.is_codeword(&corrected) {
            return Err("Unable to correct errors".to_string());
        }

        Ok((corrected, found))
    }

    /// Decodes a (possibly shortened) codeword to its message bits,
    /// correcting up to `t` errors.
    ///
    /// # Errors
    /// See [`Bch::correct`].
    pub fn decode(
        &self,
        word: &[u8],
    ) -> Result<Vec<u8>, String> {
        let (cw, _) = self.correct(word)?;

        Ok(cw[..cw.len() - self.parity_len()].to_vec())
    }
}

/// Encodes data bits with a binary BCH code correcting `t` errors.
///
/// The code is the shortened narrow-sense primitive BCH code of the smallest
/// block length `2^m - 1` that can carry `data.len()` bits (see [`Bch`]);
/// the result is the systematic codeword `data ++ parity`. Use
/// [`Bch::encode`] for explicit control of `(n, k)` and for error details.
///
/// # Arguments
/// * `data` - The data bits (`0` or `1`).
/// * `t` - The error correction capability (corrects up to `t` bit errors).
///
/// # Returns
/// The encoded codeword, or an empty vector if the arguments are invalid
/// (empty or non-binary data, `t == 0`, or no supported code is large
/// enough).
#[must_use]
pub fn bch_encode(
    data: &[u8],
    t: usize,
) -> Vec<u8> {
    Bch::for_data_len(data.len(), t)
        .and_then(|code| code.encode(data))
        .unwrap_or_default()
}

/// Decodes a codeword produced by [`bch_encode`], correcting up to `t`
/// bit errors.
///
/// # Arguments
/// * `codeword` - The received codeword.
/// * `t` - The error correction capability.
///
/// # Returns
/// The decoded data bits, or an error if too many errors are present.
///
/// # Errors
/// Returns an error if:
/// - `t` is zero or the codeword length matches no shortened BCH code;
/// - Too many errors are detected to perform reliable correction.
pub fn bch_decode(
    codeword: &[u8],
    t: usize,
) -> Result<Vec<u8>, String> {
    // `bch_encode` picks the smallest m whose code holds the data, so the
    // data length is the unique self-consistent one.
    for m in 2..=BCH_MAX_M {
        let Ok(code) = Bch::new(m, t) else { continue };

        let p = code.parity_len();

        if codeword.len() <= p {
            continue;
        }

        let data_len = codeword.len() - p;

        if data_len > code.k() {
            continue;
        }

        let smaller_fits = (2..m).any(|m2| Bch::new(m2, t).is_ok_and(|c| c.k() >= data_len));

        if smaller_fits {
            continue;
        }

        return code.decode(codeword);
    }

    Err("Codeword too short or invalid for any BCH code with this t".to_string())
}

// ============================================================================
// CRC Checksums
// ============================================================================

/// CRC-32 polynomial (IEEE 802.3 / ISO 3309 / PKZIP)
const CRC32_POLYNOMIAL: u32 = 0xEDB8_8320;

/// CRC-16 polynomial (IBM / ANSI)
const CRC16_POLYNOMIAL: u16 = 0xA001;

/// CRC-8 polynomial (ITU)
const CRC8_POLYNOMIAL: u8 = 0x07;

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
/// use rssn::kernels::error_correction::crc32_compute_numerical;
///
/// let data = b"Hello, World!";
///
/// let checksum = crc32_compute_numerical(data);
///
/// assert_eq!(checksum, 0xEC4AC3D0);
/// ```
#[must_use]
pub fn crc32_compute_numerical(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;

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
    data: &[u8],
    expected_crc: u32,
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
    crc: u32,
    data: &[u8],
) -> u32 {
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

/// Finalizes a CRC-32 computation started with `crc32_update_numerical`.
///
/// # Arguments
/// * `crc` - The running CRC value.
///
/// # Returns
/// The final CRC-32 checksum.
#[must_use]
pub const fn crc32_finalize_numerical(crc: u32) -> u32 {
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
/// use rssn::kernels::error_correction::crc16_compute;
///
/// let data = b"123456789";
///
/// let checksum = crc16_compute(data);
///
/// assert_eq!(checksum, 0xBB3D);
/// ```
#[must_use]
pub fn crc16_compute(data: &[u8]) -> u16 {
    let mut crc: u16 = 0x0000;

    for byte in data {
        crc ^= u16::from(*byte);

        for _ in 0..8 {
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
pub fn crc8_compute(data: &[u8]) -> u8 {
    let mut crc: u8 = 0;

    for byte in data {
        crc ^= *byte;

        for _ in 0..8 {
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
    data: &[u8],
    depth: usize,
) -> Vec<u8> {
    if depth == 0 || data.is_empty() {
        return data.to_vec();
    }

    let n = data.len();

    let rows = n.div_ceil(depth);

    let mut result = Vec::with_capacity(n);

    for col in 0..depth {
        for row in 0..rows {
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
    data: &[u8],
    depth: usize,
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

    for col in 0..depth {
        let col_len = if col < full_cols {
            rows
        } else {
            rows - 1
        };

        for row in 0..col_len {
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
pub fn convolutional_encode(data: &[u8]) -> Vec<u8> {
    let mut state: u8 = 0;

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
    for _ in 0..2 {
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
pub fn minimum_distance(codewords: &[Vec<u8>]) -> Option<usize> {
    if codewords.len() < 2 {
        return None;
    }

    let mut min_dist: Option<usize> = None;

    for i in 0..codewords.len() {
        for j in (i + 1)..codewords.len() {
            if let Some(dist) = hamming_distance_numerical(&codewords[i], &codewords[j]) {
                min_dist = Some(min_dist.map_or(dist, |m| m.min(dist)));
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
    k: usize,
    n: usize,
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
pub const fn error_correction_capability(min_distance: usize) -> usize {
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
pub const fn error_detection_capability(min_distance: usize) -> usize {
    if min_distance == 0 {
        0
    } else {
        min_distance - 1
    }
}
