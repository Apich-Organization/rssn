//! # Finite Fields (Galois Fields)
//!
//! This module provides structures and functions for arithmetic in finite fields.
//! It is a foundational component for advanced algebra, cryptography, and error-correcting codes.

use crate::symbolic::core::Expr;
use crate::symbolic::number_theory::extended_gcd;
// use num_bigint::{BigInt, ToBigInt as OtherToBigInt};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::Arc;

#[derive(Debug, PartialEq, Eq)]
pub struct FiniteField {
    pub modulus: BigInt,
}

impl FiniteField {
    /// Creates a new finite field `GF(modulus)`.
    ///
    /// # Arguments
    /// * `modulus` - The characteristic of the finite field (a prime number).
    ///
    /// # Returns
    /// An `Arc<Self>` pointing to the newly created field structure.
    pub fn new(modulus: i64) -> Arc<Self> {
        Arc::new(FiniteField {
            modulus: BigInt::from(modulus),
        })
    }
}

#[derive(Debug, Clone)]
pub struct FieldElement {
    pub value: BigInt,
    pub field: Arc<FiniteField>,
}

impl FieldElement {
    /// Creates a new element in a finite field.
    ///
    /// The value is reduced modulo the field's characteristic.
    ///
    /// # Arguments
    /// * `value` - The initial `BigInt` value of the element.
    /// * `field` - An `Arc` pointing to the `FiniteField` this element belongs to.
    ///
    /// # Returns
    /// A new `FieldElement`.
    pub fn new(value: BigInt, field: Arc<FiniteField>) -> Self {
        FieldElement {
            value: value % &field.modulus,
            field,
        }
    }

    /// Computes the multiplicative inverse of the element in the finite field.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod modulus)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(FieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    pub fn inverse(&self) -> Option<Self> {
        let (g, x, _) = extended_gcd(
            &Expr::BigInt(self.value.clone()),
            &Expr::BigInt(self.field.modulus.clone()),
        );
        if let Expr::BigInt(g_val) = g {
            if g_val.is_one() {
                let inv = x.to_bigint().unwrap_or_default();
                let modulus = &self.field.modulus;
                return Some(FieldElement::new(
                    (inv % modulus + modulus) % modulus,
                    self.field.clone(),
                ));
            }
        }
        None
    }
}

// --- Trait Implementations for FieldElement ---

impl Add for FieldElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot add elements from different fields.");
        }
        let val = (self.value + rhs.value) % &self.field.modulus;
        FieldElement::new(val, self.field.clone())
    }
}

impl Sub for FieldElement {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot subtract elements from different fields.");
        }
        let val = (self.value - rhs.value + &self.field.modulus) % &self.field.modulus;
        FieldElement::new(val, self.field.clone())
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot multiply elements from different fields.");
        }
        let val = (self.value * rhs.value) % &self.field.modulus;
        FieldElement::new(val, self.field.clone())
    }
}

impl Div for FieldElement {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot divide elements from different fields.");
        }
        let inv_rhs = rhs
            .inverse()
            .expect("Division by zero or non-invertible element.");
        self * inv_rhs
    }
}

impl Neg for FieldElement {
    type Output = Self;
    fn neg(self) -> Self {
        let val = (-self.value + &self.field.modulus) % &self.field.modulus;
        FieldElement::new(val, self.field.clone())
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.field == other.field && self.value == other.value
    }
}
impl Eq for FieldElement {}

// =====================================================================================
// region: Extension Field GF(2^8) for Reed-Solomon
// =====================================================================================

const GF256_GENERATOR_POLY: u16 = 0x11d; // x^8 + x^4 + x^3 + x^2 + 1
const GF256_MODULUS: usize = 256;

static mut GF256_LOG: [u8; GF256_MODULUS] = [0; GF256_MODULUS];
static mut GF256_EXP: [u8; GF256_MODULUS] = [0; GF256_MODULUS];
static mut GF256_TABLES_INITIALIZED: bool = false;

pub(crate) fn init_gf256_tables() {
    unsafe {
        if GF256_TABLES_INITIALIZED {
            return;
        }
        let mut x: u16 = 1;
        for i in 0..255 {
            GF256_EXP[i] = x as u8;
            GF256_LOG[x as usize] = i as u8;
            x <<= 1;
            if x >= 256 {
                x ^= GF256_GENERATOR_POLY;
            }
        }
        GF256_EXP[255] = GF256_EXP[0];
        GF256_TABLES_INITIALIZED = true;
    }
}

/// Computes the exponentiation (anti-logarithm) in GF(2^8).
///
/// This function uses a precomputed lookup table to find the field element
/// corresponding to a given logarithm.
///
/// # Arguments
/// * `log_val` - The logarithm of the element.
///
/// # Returns
/// The field element `alpha ^ log_val`.
pub fn gf256_exp(log_val: u8) -> u8 {
    init_gf256_tables();
    unsafe { GF256_EXP[log_val as usize] }
}

/// Performs addition in the finite field GF(2^8).
///
/// In fields of characteristic 2, addition is equivalent to a bitwise XOR operation.
#[inline]
pub fn gf256_add(a: u8, b: u8) -> u8 {
    a ^ b
}

/// Performs multiplication in GF(2^8) using precomputed lookup tables.
///
/// Multiplication is performed by adding the logarithms of the operands and then
/// finding the anti-logarithm of the result.
#[inline]
pub fn gf256_mul(a: u8, b: u8) -> u8 {
    init_gf256_tables();
    if a == 0 || b == 0 {
        0
    } else {
        unsafe {
            let log_a = GF256_LOG[a as usize] as u16;
            let log_b = GF256_LOG[b as usize] as u16;
            GF256_EXP[((log_a + log_b) % 255) as usize]
        }
    }
}

/// Computes the multiplicative inverse of an element in GF(2^8).
///
/// The inverse is calculated using the logarithm and exponentiation tables.
///
/// # Panics
/// Panics if the input `a` is 0, as 0 has no multiplicative inverse.
#[inline]
pub fn gf256_inv(a: u8) -> u8 {
    init_gf256_tables();
    if a == 0 {
        panic!("Cannot invert 0");
    }
    unsafe { GF256_EXP[(255 - GF256_LOG[a as usize] as u16) as usize] }
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the multiplicative inverse of the divisor.
///
/// # Panics
/// Panics if the divisor `b` is 0.
#[inline]
pub fn gf256_div(a: u8, b: u8) -> u8 {
    if b == 0 {
        panic!("Division by zero");
    }
    if a == 0 {
        return 0;
    }
    init_gf256_tables();
    unsafe {
        let log_a = GF256_LOG[a as usize] as u16;
        let log_b = GF256_LOG[b as usize] as u16;
        GF256_EXP[((log_a + 255 - log_b) % 255) as usize]
    }
}

// =====================================================================================
// region: Polynomial Operations over Finite Fields
// =====================================================================================

/// Evaluates a polynomial over GF(2^8) at a given point `x`.
///
/// This function uses Horner's method for efficient polynomial evaluation.
///
/// # Arguments
/// * `poly` - A slice of `u8` representing the polynomial coefficients.
/// * `x` - The point at which to evaluate the polynomial.
///
/// # Returns
/// The result of the polynomial evaluation as a `u8`.
pub fn poly_eval_gf256(poly: &[u8], x: u8) -> u8 {
    let mut y = 0;
    for coeff in poly.iter() {
        y = gf256_mul(y, x) ^ coeff;
    }
    y
}

/// Adds two polynomials over GF(2^8).
///
/// Polynomial addition in GF(2^8) is performed by XORing corresponding coefficients.
///
/// # Arguments
/// * `p1` - The first polynomial as a slice of `u8` coefficients.
/// * `p2` - The second polynomial as a slice of `u8` coefficients.
///
/// # Returns
/// A `Vec<u8>` representing the sum polynomial.
pub fn poly_add_gf256(p1: &[u8], p2: &[u8]) -> Vec<u8> {
    let mut result = vec![0; std::cmp::max(p1.len(), p2.len())];
    let res_len = result.len();
    for i in 0..p1.len() {
        result[i + res_len - p1.len()] = p1[i];
    }
    for i in 0..p2.len() {
        result[i + res_len - p2.len()] ^= p2[i];
    }
    result
}

/// Multiplies two polynomials over GF(2^8).
///
/// Polynomial multiplication is performed by convolving the coefficients,
/// with each coefficient multiplication and addition done in GF(2^8).
///
/// # Arguments
/// * `p1` - The first polynomial as a slice of `u8` coefficients.
/// * `p2` - The second polynomial as a slice of `u8` coefficients.
///
/// # Returns
/// A `Vec<u8>` representing the product polynomial.
pub fn poly_mul_gf256(p1: &[u8], p2: &[u8]) -> Vec<u8> {
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

/// Divides two polynomials over GF(2^8).
///
/// This function performs polynomial long division. It returns the remainder.
///
/// # Arguments
/// * `dividend` - The dividend polynomial as a `Vec<u8>` (will be consumed).
/// * `divisor` - The divisor polynomial as a slice of `u8` coefficients.
///
/// # Returns
/// A `Vec<u8>` representing the remainder polynomial.
///
/// # Panics
/// Panics if the divisor is empty.
pub fn poly_div_gf256(mut dividend: Vec<u8>, divisor: &[u8]) -> Vec<u8> {
    if divisor.is_empty() {
        panic!("Divisor cannot be empty");
    }

    let divisor_len = divisor.len();
    let lead_divisor = divisor[0];
    let lead_divisor_inv = gf256_inv(lead_divisor);

    while dividend.len() >= divisor_len {
        let lead_dividend = dividend[0];
        let coeff = gf256_mul(lead_dividend, lead_divisor_inv);

        for i in 0..divisor_len {
            let term = gf256_mul(coeff, divisor[i]);
            dividend[i] ^= term;
        }
        dividend.remove(0);
    }
    dividend
}

pub(crate) fn expr_to_field_elements(
    p_expr: &Expr,
    field: &Arc<FiniteField>,
) -> Result<Vec<FieldElement>, String> {
    if let Expr::Polynomial(coeffs) = p_expr {
        coeffs
            .iter()
            .map(|c| {
                c.to_bigint()
                    .map(|val| FieldElement::new(val, field.clone()))
                    .ok_or_else(|| format!("Invalid coefficient in polynomial: {}", c))
            })
            .collect()
    } else {
        Err(format!("Expression is not a polynomial: {}", p_expr))
    }
}

pub(crate) fn field_elements_to_expr(coeffs: &[FieldElement]) -> Expr {
    let expr_coeffs = coeffs
        .iter()
        .map(|c| Expr::BigInt(c.value.clone()))
        .collect();
    Expr::Polynomial(expr_coeffs)
}

/// Adds two polynomials whose coefficients are `FieldElement`s from a given finite field.
///
/// # Arguments
/// * `p1_expr` - The first polynomial as an `Expr::Polynomial`.
/// * `p2_expr` - The second polynomial as an `Expr::Polynomial`.
/// * `field` - The finite field over which the polynomials are defined.
///
/// # Returns
/// * `Ok(Expr::Polynomial)` representing the sum.
/// * `Err(String)` if the input expressions are not valid polynomials or contain invalid coefficients.
pub fn poly_add_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: Arc<FiniteField>,
) -> Result<Expr, String> {
    let c1 = expr_to_field_elements(p1_expr, &field)?;
    let c2 = expr_to_field_elements(p2_expr, &field)?;
    let mut result_coeffs = vec![];

    let len1 = c1.len();
    let len2 = c2.len();
    let max_len = std::cmp::max(len1, len2);

    for i in 0..max_len {
        let val1 = if i < len1 {
            c1[len1 - 1 - i].clone()
        } else {
            FieldElement::new(Zero::zero(), field.clone())
        };
        let val2 = if i < len2 {
            c2[len2 - 1 - i].clone()
        } else {
            FieldElement::new(Zero::zero(), field.clone())
        };
        result_coeffs.push(val1 + val2);
    }
    result_coeffs.reverse();

    Ok(field_elements_to_expr(&result_coeffs))
}

/// Multiplies two polynomials whose coefficients are `FieldElement`s from a given finite field.
///
/// # Arguments
/// * `p1_expr` - The first polynomial as an `Expr::Polynomial`.
/// * `p2_expr` - The second polynomial as an `Expr::Polynomial`.
/// * `field` - The finite field over which the polynomials are defined.
///
/// # Returns
/// * `Ok(Expr::Polynomial)` representing the product.
/// * `Err(String)` if the input expressions are not valid polynomials or contain invalid coefficients.
pub fn poly_mul_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: Arc<FiniteField>,
) -> Result<Expr, String> {
    let c1 = expr_to_field_elements(p1_expr, &field)?;
    let c2 = expr_to_field_elements(p2_expr, &field)?;

    if c1.is_empty() || c2.is_empty() {
        return Ok(Expr::Polynomial(vec![]));
    }

    let deg1 = c1.len() - 1;
    let deg2 = c2.len() - 1;
    let mut result_coeffs = vec![FieldElement::new(Zero::zero(), field.clone()); deg1 + deg2 + 1];

    for i in 0..=deg1 {
        for j in 0..=deg2 {
            let term_mul = c1[i].clone() * c2[j].clone();
            result_coeffs[i + j] = result_coeffs[i + j].clone() + term_mul;
        }
    }

    Ok(field_elements_to_expr(&result_coeffs))
}

/// Divides two polynomials whose coefficients are `FieldElement`s from a given finite field.
///
/// # Arguments
/// * `p1_expr` - The dividend polynomial as an `Expr::Polynomial`.
/// * `p2_expr` - The divisor polynomial as an `Expr::Polynomial`.
/// * `field` - The finite field over which the polynomials are defined.
///
/// # Returns
/// * `Ok((Expr::Polynomial, Expr::Polynomial))` representing the quotient and remainder.
/// * `Err(String)` if the input expressions are not valid polynomials, contain invalid coefficients,
///   or if division by the zero polynomial is attempted.
pub fn poly_div_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: Arc<FiniteField>,
) -> Result<(Expr, Expr), String> {
    let mut num = expr_to_field_elements(p1_expr, &field)?;
    let den = expr_to_field_elements(p2_expr, &field)?;

    if den.iter().all(|c| c.value.is_zero()) {
        return Err("Division by zero polynomial".to_string());
    }

    let mut quotient = vec![FieldElement::new(Zero::zero(), field.clone()); num.len()];

    let lead_den_inv = den
        .first()
        .unwrap()
        .inverse()
        .ok_or("Leading coefficient is not invertible".to_string())?;

    while num.len() >= den.len() {
        let lead_num = num.first().unwrap().clone();
        let coeff = lead_num * lead_den_inv.clone();
        let degree_diff = num.len() - den.len();
        quotient[degree_diff] = coeff.clone();

        for (i, den_coeff) in den.iter().enumerate() {
            let term = coeff.clone() * den_coeff.clone();
            num[i] = num[i].clone() - term;
        }
        num.remove(0);
    }

    let first_non_zero = num
        .iter()
        .position(|c| !c.value.is_zero())
        .unwrap_or(num.len());
    let remainder = &num[first_non_zero..];

    Ok((
        field_elements_to_expr(&quotient),
        field_elements_to_expr(remainder),
    ))
}

// Helper trait to convert Expr to BigInt
trait ToBigInt {
    fn to_bigint(&self) -> Option<BigInt>;
}

impl ToBigInt for Expr {
    fn to_bigint(&self) -> Option<BigInt> {
        match self {
            Expr::BigInt(i) => Some(i.clone()),
            Expr::Constant(_) => None, // Cannot convert f64 to BigInt without loss
            _ => None,
        }
    }
}
