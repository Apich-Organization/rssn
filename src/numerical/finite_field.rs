//! # Finite Field Arithmetic
//!
//! This module provides numerical implementations for arithmetic in finite fields.
//! It includes support for prime fields GF(p) and optimized arithmetic for GF(2^8).

use std::ops::{Add, Div, Mul, Sub};

// =====================================================================================
// region: Prime Fields GF(p) over u64
// =====================================================================================

/// Represents an element in a prime field GF(p), where p is the modulus.
///
/// The value is stored as a `u64`, and all arithmetic operations are performed
/// modulo the specified `modulus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimeFieldElement {
    /// The value of the field element.
    pub value: u64,
    /// The modulus (prime number) of the field.
    pub modulus: u64,
}

impl PrimeFieldElement {
    /// Creates a new `PrimeFieldElement`.
    ///
    /// The initial value is reduced modulo the `modulus`.
    ///
    /// # Arguments
    /// * `value` - The initial value of the element.
    /// * `modulus` - The prime modulus of the field.
    pub fn new(value: u64, modulus: u64) -> Self {
        PrimeFieldElement {
            value: value % modulus,
            modulus,
        }
    }

    /// Computes the multiplicative inverse of the element.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod modulus)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(PrimeFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    pub fn inverse(&self) -> Option<Self> {
        let (g, x, _) = extended_gcd_u64(self.value as i64, self.modulus as i64);
        if g == 1 {
            let inv = (x % self.modulus as i64 + self.modulus as i64) % self.modulus as i64;
            Some(PrimeFieldElement::new(inv as u64, self.modulus))
        } else {
            None
        }
    }
}

/// Implements the Extended Euclidean Algorithm for `i64` integers.
///
/// This function computes the greatest common divisor (GCD) of `a` and `b`
/// and also finds coefficients `x` and `y` such that `a*x + b*y = gcd(a, b)`.
///
/// # Arguments
/// * `a` - The first integer.
/// * `b` - The second integer.
///
/// # Returns
/// A tuple `(g, x, y)` where `g` is the GCD, and `x`, `y` are the Bézout coefficients.
pub(crate) fn extended_gcd_u64(a: i64, b: i64) -> (i64, i64, i64) {
    /// Implements the Extended Euclidean Algorithm for `i64` integers.
    ///
    /// This function computes the greatest common divisor (GCD) of `a` and `b`
    /// and also finds coefficients `x` and `y` such that `a*x + b*y = gcd(a, b)`.
    ///
    /// # Arguments
    /// * `a` - The first integer.
    /// * `b` - The second integer.
    ///
    /// # Returns
    /// A tuple `(g, x, y)` where `g` is the GCD, and `x`, `y` are the Bézout coefficients.
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd_u64(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

impl Add for PrimeFieldElement {
    type Output = Self;
    /// Performs addition in the prime field.
    fn add(self, rhs: Self) -> Self {
        let val = (self.value + rhs.value) % self.modulus;
        Self::new(val, self.modulus)
    }
}

impl Sub for PrimeFieldElement {
    type Output = Self;
    /// Performs subtraction in the prime field.
    fn sub(self, rhs: Self) -> Self {
        let val = (self.value + self.modulus - rhs.value) % self.modulus;
        Self::new(val, self.modulus)
    }
}

impl Mul for PrimeFieldElement {
    type Output = Self;
    /// Performs multiplication in the prime field.
    ///
    /// Uses `u128` for the intermediate multiplication to prevent overflow
    /// before the modulo operation.
    fn mul(self, rhs: Self) -> Self {
        let val = ((self.value as u128 * rhs.value as u128) % self.modulus as u128) as u64;
        Self::new(val, self.modulus)
    }
}

impl Div for PrimeFieldElement {
    type Output = Self;
    /// Performs division in the prime field.
    ///
    /// This is achieved by multiplying by the multiplicative inverse of the divisor.
    ///
    /// # Panics
    /// Panics if the divisor is zero or not invertible.
    fn div(self, rhs: Self) -> Self {
        let inv_rhs = rhs
            .inverse()
            .expect("Division by zero or non-invertible element.");
        self * inv_rhs
    }
}

// =====================================================================================
// region: Optimized GF(2^8) Arithmetic
// =====================================================================================

const GF256_GENERATOR_POLY: u16 = 0x11d; // x^8 + x^4 + x^3 + x^2 + 1
const GF256_MODULUS: usize = 256;

static mut GF256_LOG: [u8; GF256_MODULUS] = [0; GF256_MODULUS];
static mut GF256_EXP: [u8; GF256_MODULUS] = [0; GF256_MODULUS];
static mut GF256_TABLES_INITIALIZED: bool = false;

/// Initializes the logarithm and exponentiation tables for GF(2^8) arithmetic.
///
/// This function precomputes lookup tables to speed up multiplication and division
/// in the field. It is called automatically on the first use of a GF(2^8) arithmetic function.
/// This function is not thread-safe if called concurrently.
#[allow(unsafe_code)]
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

/// Performs addition in GF(2^8).
///
/// In GF(2^8), addition is equivalent to a bitwise XOR operation.
#[inline]
pub fn gf256_add(a: u8, b: u8) -> u8 {
    a ^ b
}

/// Performs multiplication in GF(2^8) using lookup tables.
///
/// This function uses precomputed logarithm and exponentiation tables for efficiency.
///
/// # Panics
/// This function is not thread-safe during the first call that initializes the tables.
#[allow(unsafe_code)]
#[inline]
pub fn gf256_mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }
    init_gf256_tables();
    unsafe {
        let log_a = GF256_LOG[a as usize] as u16;
        let log_b = GF256_LOG[b as usize] as u16;
        GF256_EXP[((log_a + log_b) % 255) as usize]
    }
}

/// Computes the multiplicative inverse in GF(2^8).
///
/// # Panics
/// * Panics if `a` is 0, as 0 has no multiplicative inverse.
/// * This function is not thread-safe during the first call that initializes the tables.
#[allow(unsafe_code)]
#[inline]
pub fn gf256_inv(a: u8) -> u8 {
    if a == 0 {
        panic!("Cannot invert 0");
    }
    init_gf256_tables();
    unsafe { GF256_EXP[(255 - GF256_LOG[a as usize] as u16) as usize] }
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the inverse of the divisor.
///
/// # Panics
/// * Panics if the divisor `b` is 0.
/// * This function is not thread-safe during the first call that initializes the tables.
#[allow(unsafe_code)]
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
