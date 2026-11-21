//! # Finite Field Arithmetic
//!
//! This module provides numerical implementations for arithmetic in finite fields.
//! It includes support for prime fields GF(p) and arithmetic for GF(2^8) in the
//! polynomial basis modulo `0x11d`. Use [`Gf256Multiplier`] to reuse a fixed
//! coefficient across byte slices.

use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;

use serde::Deserialize;
use serde::Serialize;

/// Represents an element in a prime field GF(p), where p is the modulus.
///
/// The value is stored as a `u64`, and all arithmetic operations are performed
/// modulo the specified `modulus`.
///
/// # Modulus-agnostic elements
/// A `modulus` of `0` marks a *modulus-agnostic* element: the integer
/// `value` (a non-negative multiple of the field's `1`) that has not yet
/// been bound to a field. [`Zero::zero`] and [`One::one`] return such
/// elements, so they can be created without knowing the field. In a binary
/// operation the result adopts the modulus of the concrete operand, and the
/// agnostic value is reduced modulo it; `x + zero() == x` and
/// `x * one() == x` hold for every concrete `x`. Between two agnostic
/// elements the operations act on non-negative integers (addition and
/// multiplication saturate at `u64::MAX`, subtraction saturates at `0`);
/// negating a non-zero agnostic element is not representable and leaves it
/// unchanged. When both operands are concrete with different moduli the
/// left operand's modulus is used.
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
            value: if modulus == 0 {
                value
            } else {
                value % modulus
            },
            modulus,
        }
    }

    /// Returns `true` for a modulus-agnostic element (see the type docs).
    #[must_use]
    pub const fn is_agnostic(&self) -> bool {
        self.modulus == 0
    }

    /// Modulus the result of a binary operation adopts: the left operand's
    /// unless it is agnostic.
    const fn result_modulus(
        &self,
        rhs: &Self,
    ) -> u64 {
        if self.modulus == 0 {
            rhs.modulus
        } else {
            self.modulus
        }
    }

    /// The value reduced modulo `m` (unchanged for `m == 0`).
    const fn reduced(
        &self,
        m: u64,
    ) -> u64 {
        if m == 0 {
            self.value
        } else {
            self.value % m
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
    ///
    /// # Panics
    /// Panics if the internal modular inverse calculation results in a value
    /// that cannot be converted to `u64`. This is not expected to happen if
    /// the modulus is a valid `u64`.
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        if self.modulus == 0 {
            // Only 1 is invertible among agnostic integers.
            return (self.value == 1).then_some(*self);
        }

        let (g, x, _) = extended_gcd_u64(self.value, self.modulus);

        if g == 1 {
            let m = i128::from(self.modulus);

            let inv = (x % m + m) % m;

            Some(Self::new(u64::try_from(inv).unwrap(), self.modulus))
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
        if self.modulus == 0 {
            // Agnostic integers: saturating power.
            return Self::new(
                self.value
                    .saturating_pow(u32::try_from(exp).unwrap_or(u32::MAX)),
                0,
            );
        }

        let mut res = 1u128 % u128::from(self.modulus);

        let mut base = u128::from(self.value);

        let m = u128::from(self.modulus);

        while exp > 0 {
            if exp % 2 == 1 {
                res = (res * base) % m;
            }

            base = (base * base) % m;

            exp /= 2;
        }

        Self::new(res as u64, self.modulus)
    }
}

use std::ops::Neg;

use num_traits::One;
use num_traits::Zero;

impl Zero for PrimeFieldElement {
    fn zero() -> Self {
        // Modulus-agnostic zero: adopts the modulus of the other operand.
        Self { value: 0, modulus: 0 }
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }
}

impl One for PrimeFieldElement {
    fn one() -> Self {
        // Modulus-agnostic one: adopts the modulus of the other operand.
        Self { value: 1, modulus: 0 }
    }
}

impl Neg for PrimeFieldElement {
    type Output = Self;

    fn neg(self) -> Self {
        if self.value == 0 || self.modulus == 0 {
            self
        } else {
            Self::new(self.modulus - self.value, self.modulus)
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

        (g, y - (i128::from(b) / i128::from(a)) * x, x)
    }
}

impl Add for PrimeFieldElement {
    type Output = Self;

    /// Performs addition in the prime field.
    fn add(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        let sum = u128::from(self.reduced(m)) + u128::from(rhs.reduced(m));

        if m == 0 {
            return Self::new(u64::try_from(sum).unwrap_or(u64::MAX), 0);
        }

        Self::new((sum % u128::from(m)) as u64, m)
    }
}

impl Sub for PrimeFieldElement {
    type Output = Self;

    /// Performs subtraction in the prime field.
    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        if m == 0 {
            return Self::new(self.value.saturating_sub(rhs.value), 0);
        }

        let val = (u128::from(self.reduced(m)) + u128::from(m) - u128::from(rhs.reduced(m)))
            % u128::from(m);

        Self::new(val as u64, m)
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
        let m = self.result_modulus(&rhs);

        let prod = u128::from(self.reduced(m)) * u128::from(rhs.reduced(m));

        if m == 0 {
            return Self::new(u64::try_from(prod).unwrap_or(u64::MAX), 0);
        }

        Self::new((prod % u128::from(m)) as u64, m)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for PrimeFieldElement {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {
        let m = self.result_modulus(&rhs);

        let rhs = Self::new(rhs.value, m);

        let Some(inv_rhs) = rhs.inverse() else {
            return Self::new(0, m);
        };

        self * inv_rhs
    }
}

const GF256_GENERATOR_POLY: u16 = 0x11d;
const GF256_SIZE: usize = 256;
const GF256_ORDER: usize = GF256_SIZE - 1;

struct Gf256Tables {
    log: [u8; GF256_SIZE],
    // Repeat the cycle to index by either a sum or an offset difference of logs.
    // Two spare entries give an array large enough for any two byte indices.
    exp: [u8; 2 * GF256_SIZE],
}

#[expect(
    clippy::indexing_slicing,
    reason = "the loop bounds the exponent index and polynomial reduction leaves a byte"
)]
static GF256_TABLES: Gf256Tables = {
    let mut tables = Gf256Tables {
        log: [0; GF256_SIZE],
        exp: [0; 2 * GF256_SIZE],
    };
    let mut i = 0;
    let mut x: u16 = 1;
    while i < tables.exp.len() {
        tables.exp[i] = x as u8;
        if i < GF256_ORDER {
            tables.log[x as usize] = i as u8;
        }
        x <<= 1;
        if x >= GF256_SIZE as u16 {
            x ^= GF256_GENERATOR_POLY;
        }
        i += 1;
    }
    tables
};

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
/// The compile-time logarithm and exponentiation tables use the `0x11d` polynomial basis.
#[inline]
#[must_use]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "byte log indices fit their table and the sum of two bytes is at most 510"
)]
pub const fn gf256_mul(
    a: u8,
    b: u8,
) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }

    let log_a = GF256_TABLES.log[a as usize] as usize;
    let log_b = GF256_TABLES.log[b as usize] as usize;
    GF256_TABLES.exp[log_a + log_b]
}

/// Computes the multiplicative inverse in GF(2^8).
///
/// # Errors
/// * Returns an error if `a` is 0, as 0 has no multiplicative inverse.
#[inline]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "subtracting a byte log from 255 gives an exponent in 0..=255"
)]
pub fn gf256_inv(a: u8) -> Result<u8, String> {
    if a == 0 {
        return Err("Cannot invert 0".to_string());
    }

    Ok(GF256_TABLES.exp[GF256_ORDER - usize::from(GF256_TABLES.log[usize::from(a)])])
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the inverse of the divisor.
///
/// # Errors
/// * Returns an error if the divisor `b` is 0.
#[inline]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "adding 255 to the log bounds the nonnegative difference by 510"
)]
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

    let log_a = usize::from(GF256_TABLES.log[usize::from(a)]);
    let log_b = usize::from(GF256_TABLES.log[usize::from(b)]);
    Ok(GF256_TABLES.exp[log_a + GF256_ORDER - log_b])
}

/// Performs exponentiation in GF(2^8).
///
/// Exponent zero returns `1`, including for a zero base.
#[must_use]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "the reduced exponent product fits u16 and the residue is at most 254"
)]
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

    let log_a = u16::from(GF256_TABLES.log[usize::from(a)]);
    // Nonzero elements have order dividing 255, so reducing first keeps the product in u16.
    let reduced_exp = (exp % GF256_ORDER as u64) as u16;
    let new_log = (log_a * reduced_exp) % GF256_ORDER as u16;
    GF256_TABLES.exp[usize::from(new_log)]
}

/// A reusable multiplier for a fixed coefficient in GF(2^8).
///
/// It owns a 256-byte lookup table and needs no heap allocation. Reuse it for operations
/// on several slices with the same coefficient, such as erasure-coding stripes.
///
/// ```
/// use rssn::kernels::finite_field::Gf256Multiplier;
///
/// let multiplier = Gf256Multiplier::new(2);
/// let mut parity = [0x10, 0x20];
/// multiplier.mul_add_assign(&[0x80, 1], &mut parity)?;
/// assert_eq!(parity, [0x0d, 0x22]);
/// # Ok::<(), String>(())
/// ```
#[derive(Clone, Debug)]
pub struct Gf256Multiplier {
    coefficient: u8,
    products: [u8; GF256_SIZE],
}

#[expect(
    clippy::indexing_slicing,
    reason = "the constructor and byte lookups index a 256-entry table"
)]
impl Gf256Multiplier {
    /// Builds a multiplier using the same polynomial basis as [`gf256_mul`].
    #[must_use]
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the table loop increments only while the index is below 256"
    )]
    pub const fn new(coefficient: u8) -> Self {
        let mut products = [0; GF256_SIZE];
        if coefficient > 1 {
            let mut value = 1;
            while value < GF256_SIZE {
                products[value] = gf256_mul(coefficient, value as u8);
                value += 1;
            }
        }
        Self {
            coefficient,
            products,
        }
    }

    /// Multiplies each byte of the slice by the coefficient in place.
    pub fn scale_assign(
        &self,
        values: &mut [u8],
    ) {
        match self.coefficient {
            | 0 => values.fill(0),
            | 1 => {},
            | _ => {
                for value in values {
                    *value = self.products[usize::from(*value)];
                }
            },
        }
    }

    /// Adds the scaled source into the destination using field addition (XOR).
    ///
    /// # Errors
    /// Returns an error and leaves the destination unchanged when the lengths differ.
    pub fn mul_add_assign(
        &self,
        source: &[u8],
        destination: &mut [u8],
    ) -> Result<(), String> {
        if source.len() != destination.len() {
            return Err("GF(256) source and destination lengths must match".into());
        }
        match self.coefficient {
            | 0 => {},
            | 1 => {
                for (&src, dst) in source.iter().zip(destination) {
                    *dst ^= src;
                }
            },
            | _ => {
                for (&src, dst) in source.iter().zip(destination) {
                    *dst ^= self.products[usize::from(src)];
                }
            },
        }
        Ok(())
    }
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
