//! # Finite Field Arithmetic
//!
//! This module provides numerical implementations for arithmetic in finite fields.
//! It includes support for prime fields GF(p) and optimized arithmetic for GF(2^8).

use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

/// Represents an element in a prime field GF(p), where p is the modulus.
///
/// The value is stored as a `u64`, and all arithmetic operations are performed
/// modulo the specified `modulus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]

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
    #[must_use]

    pub const fn new(
        value: u64,
        modulus: u64,
    ) -> Self {

        Self {
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
    #[must_use]

    pub fn inverse(&self) -> Option<Self> {

        let (g, x, _) = extended_gcd_u64(
            self.value,
            self.modulus,
        );

        if g == 1 {

            let m = i128::from(self.modulus);

            let inv = (x % m + m) % m;

            Some(Self::new(
                inv as u64,
                self.modulus,
            ))
        } else {

            None
        }
    }

    /// Computes (self^exp) modulo modulus.
    #[must_use]

    pub fn pow(
        &self,
        mut exp: u64,
    ) -> Self {

        let mut res = 1u128;

        let mut base = u128::from(self.value);

        let m = u128::from(self.modulus);

        while exp > 0 {

            if exp % 2 == 1 {

                res = (res * base) % m;
            }

            base = (base * base) % m;

            exp /= 2;
        }

        Self::new(
            res as u64,
            self.modulus,
        )
    }
}

use num_traits::{One, Zero};
use std::ops::Neg;

impl Zero for PrimeFieldElement {
    fn zero() -> Self {

        Self {
            value: 0,
            modulus: 2,
        } // Dummy modulus, should be careful
    }

    fn is_zero(&self) -> bool {

        self.value == 0
    }
}

impl One for PrimeFieldElement {
    fn one() -> Self {

        Self {
            value: 1,
            modulus: 2,
        } // Dummy modulus
    }
}

impl Neg for PrimeFieldElement {
    type Output = Self;

    fn neg(self) -> Self {

        if self.value == 0 {

            self
        } else {

            Self::new(
                self.modulus - self.value,
                self.modulus,
            )
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

pub(crate) fn extended_gcd_u64(
    a: u64,
    b: u64,
) -> (u64, i128, i128) {

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

        (
            g,
            y - (i128::from(b) / i128::from(a)) * x,
            x,
        )
    }
}

impl Add for PrimeFieldElement {
    type Output = Self;

    /// Performs addition in the prime field.

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        let val = (self.value + rhs.value) % self.modulus;

        Self::new(val, self.modulus)
    }
}

impl Sub for PrimeFieldElement {
    type Output = Self;

    /// Performs subtraction in the prime field.

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

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

    fn mul(
        self,
        rhs: Self,
    ) -> Self {

        let val =
            ((u128::from(self.value) * u128::from(rhs.value)) % u128::from(self.modulus)) as u64;

        Self::new(val, self.modulus)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]

impl Div for PrimeFieldElement {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {

        let inv_rhs = match rhs.inverse() {
            Some(inv) => inv,
            None => {

                return Self::new(0, self.modulus);
            }
        };

        self * inv_rhs
    }
}

const GF256_GENERATOR_POLY: u16 = 0x11d;

const GF256_MODULUS: usize = 256;

struct Gf256Tables {
    log: [u8; GF256_MODULUS],
    exp: [u8; GF256_MODULUS],
}

static GF256_TABLES: std::sync::LazyLock<Gf256Tables> = std::sync::LazyLock::new(|| {

    let mut log_table = [0u8; GF256_MODULUS];

    let mut exp_table = [0u8; GF256_MODULUS];

    let mut x: u16 = 1;

    //for i in 0..255 {
    for (i, value) in exp_table
        .iter_mut()
        .enumerate()
        .take(255)
    {

        //exp_table[i] = x as u8;
        *value = x as u8;

        log_table[x as usize] = i as u8;

        x <<= 1;

        if x >= 256 {

            x ^= GF256_GENERATOR_POLY;
        }
    }

    exp_table[255] = exp_table[0];

    Gf256Tables {
        log: log_table,
        exp: exp_table,
    }
});

/// Performs addition in GF(2^8).
///
/// In GF(2^8), addition is equivalent to a bitwise XOR operation.
#[inline]
#[must_use]

pub const fn gf256_add(
    a: u8,
    b: u8,
) -> u8 {

    a ^ b
}

/// Performs multiplication in GF(2^8) using lookup tables.
///
/// This function uses precomputed logarithm and exponentiation tables for efficiency.
#[inline]
#[must_use]

pub fn gf256_mul(
    a: u8,
    b: u8,
) -> u8 {

    if a == 0 || b == 0 {

        return 0;
    }

    let log_a = u16::from(GF256_TABLES.log[a as usize]);

    let log_b = u16::from(GF256_TABLES.log[b as usize]);

    GF256_TABLES.exp[((log_a + log_b) % 255) as usize]
}

/// Computes the multiplicative inverse in GF(2^8).
///
/// # Panics
/// * Panics if `a` is 0, as 0 has no multiplicative inverse.
#[inline]

pub fn gf256_inv(a: u8) -> Result<u8, String> {

    if a == 0 {

        return Err("Cannot invert 0".to_string());
    }

    Ok(GF256_TABLES.exp[(255 - u16::from(GF256_TABLES.log[a as usize])) as usize])
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the inverse of the divisor.
///
/// # Panics
/// * Panics if the divisor `b` is 0.
#[inline]

pub fn gf256_div(
    a: u8,
    b: u8,
) -> Result<u8, String> {

    if b == 0 {

        return Err("Division by zero".to_string());
    }

    if a == 0 {

        return Ok(0);
    }

    let log_a = u16::from(GF256_TABLES.log[a as usize]);

    let log_b = u16::from(GF256_TABLES.log[b as usize]);

    Ok(GF256_TABLES.exp[((log_a + 255 - log_b) % 255) as usize])
}

/// Performs exponentiation in GF(2^8).
#[must_use]

pub fn gf256_pow(
    a: u8,
    exp: u64,
) -> u8 {

    if exp == 0 {

        return 1;
    }

    if a == 0 {

        return 0;
    }

    let log_a = u16::from(GF256_TABLES.log[a as usize]);

    let new_log = (u128::from(log_a) * u128::from(exp)) % 255;

    GF256_TABLES.exp[new_log as usize]
}

impl AddAssign for PrimeFieldElement {
    fn add_assign(
        &mut self,
        rhs: Self,
    ) {

        *self = *self + rhs;
    }
}

impl SubAssign for PrimeFieldElement {
    fn sub_assign(
        &mut self,
        rhs: Self,
    ) {

        *self = *self - rhs;
    }
}

impl MulAssign for PrimeFieldElement {
    fn mul_assign(
        &mut self,
        rhs: Self,
    ) {

        *self = *self * rhs;
    }
}

impl DivAssign for PrimeFieldElement {
    fn div_assign(
        &mut self,
        rhs: Self,
    ) {

        *self = *self / rhs;
    }
}
