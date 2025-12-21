//! # Symbolic Finite Field Arithmetic
//!
//! This module provides symbolic structures for arithmetic in finite fields (Galois fields).
//! It defines prime fields GF(p) and extension fields GF(p^n), along with the necessary
//! arithmetic operations for their elements and for polynomials over these fields.
#![allow(clippy::should_implement_trait)]
use crate::symbolic::number_theory::extended_gcd_inner;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

// Helper module for serializing Arc<T>
mod arc_serde {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::sync::Arc;

    pub fn serialize<S, T>(arc: &Arc<T>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
        T: Serialize,
    {
        arc.as_ref().serialize(serializer)
    }

    pub fn deserialize<'de, D, T>(deserializer: D) -> Result<Arc<T>, D::Error>
    where
        D: Deserializer<'de>,
        T: Deserialize<'de>,
    {
        T::deserialize(deserializer).map(Arc::new)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrimeField {
    pub modulus: BigInt,
}
impl PrimeField {
    /// Creates a new prime field `GF(p)` with the given modulus.
    ///
    /// # Arguments
    /// * `modulus` - A `BigInt` representing the prime modulus `p` of the field.
    ///
    /// # Returns
    /// An `Arc<PrimeField>` pointing to the newly created field structure.
    #[must_use]
    pub fn new(modulus: BigInt) -> Arc<Self> {
        Arc::new(Self { modulus })
    }
}
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrimeFieldElement {
    pub value: BigInt,
    #[serde(with = "arc_serde")]
    pub field: Arc<PrimeField>,
}
impl PrimeFieldElement {
    /// Creates a new element in a prime field.
    ///
    /// The value is reduced modulo the field's characteristic (the prime modulus).
    ///
    /// # Arguments
    /// * `value` - The initial `BigInt` value of the element.
    /// * `field` - An `Arc` pointing to the `PrimeField` this element belongs to.
    ///
    /// # Returns
    /// A new `PrimeFieldElement`.
    #[must_use]
    pub fn new(value: BigInt, field: Arc<PrimeField>) -> Self {
        let modulus = &field.modulus;
        let mut val = value % modulus;
        if val < Zero::zero() {
            val += modulus;
        }
        Self { value: val, field }
    }
    /// Computes the multiplicative inverse of the element in the prime field.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod p)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(PrimeFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        let (g, x, _) = extended_gcd_inner(self.value.clone(), self.field.modulus.clone());
        if g.is_one() {
            let modulus = &self.field.modulus;
            let mut inv = x % modulus;
            if inv < Zero::zero() {
                inv += modulus;
            }
            Some(Self::new(inv, self.field.clone()))
        } else {
            None
        }
    }
}
impl Add for PrimeFieldElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            return Self::new(Zero::zero(), self.field);
        }
        let val = (self.value + rhs.value) % &self.field.modulus;
        Self::new(val, self.field)
    }
}
impl Sub for PrimeFieldElement {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            return Self::new(Zero::zero(), self.field);
        }
        let val = (self.value - rhs.value + &self.field.modulus) % &self.field.modulus;
        Self::new(val, self.field)
    }
}
impl Mul for PrimeFieldElement {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            return Self::new(Zero::zero(), self.field);
        }
        let val = (self.value * rhs.value) % &self.field.modulus;
        Self::new(val, self.field)
    }
}
impl Div for PrimeFieldElement {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            return Self::new(Zero::zero(), self.field);
        }
        let inv_rhs = match rhs.inverse() {
            Some(inv) => inv,
            None => return Self::new(Zero::zero(), self.field),
        };
        self * inv_rhs
    }
}
impl Neg for PrimeFieldElement {
    type Output = Self;
    fn neg(self) -> Self {
        let val = (-self.value + &self.field.modulus) % &self.field.modulus;
        Self::new(val, self.field)
    }
}
impl PartialEq for PrimeFieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.field == other.field && self.value == other.value
    }
}
impl Eq for PrimeFieldElement {}
impl Zero for PrimeFieldElement {
    fn zero() -> Self {
        let dummy_field = PrimeField::new(BigInt::from(2));
        Self::new(Zero::zero(), dummy_field)
    }
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}
impl One for PrimeFieldElement {
    fn one() -> Self {
        let dummy_field = PrimeField::new(BigInt::from(2));
        Self::new(One::one(), dummy_field)
    }
}
impl AddAssign for PrimeFieldElement {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}
impl SubAssign for PrimeFieldElement {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}
impl MulAssign for PrimeFieldElement {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}
impl DivAssign for PrimeFieldElement {
    fn div_assign(&mut self, rhs: Self) {
        *self = self.clone() / rhs;
    }
}
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct FiniteFieldPolynomial {
    pub coeffs: Vec<PrimeFieldElement>,
    #[serde(with = "arc_serde")]
    pub field: Arc<PrimeField>,
}
impl FiniteFieldPolynomial {
    /// Creates a new polynomial over a prime field.
    ///
    /// This function also removes leading zero coefficients to keep the representation canonical.
    ///
    /// # Arguments
    /// * `coeffs` - A vector of `PrimeFieldElement`s representing the coefficients in descending order of degree.
    /// * `field` - An `Arc` pointing to the `PrimeField` the coefficients belong to.
    ///
    /// # Returns
    /// A new `FiniteFieldPolynomial`.
    #[must_use]
    pub fn new(coeffs: Vec<PrimeFieldElement>, field: Arc<PrimeField>) -> Self {
        let first_non_zero = coeffs
            .iter()
            .position(|c| !c.value.is_zero())
            .unwrap_or(coeffs.len());
        Self {
            coeffs: coeffs[first_non_zero..].to_vec(),
            field,
        }
    }
    /// Returns the degree of the polynomial.
    ///
    /// The degree is the highest power of the variable with a non-zero coefficient.
    /// The degree of the zero polynomial is defined as -1.
    ///
    /// # Returns
    /// An `isize` representing the degree.
    #[must_use]
    pub const fn degree(&self) -> isize {
        if self.coeffs.is_empty() {
            -1
        } else {
            (self.coeffs.len() - 1) as isize
        }
    }
    /// Performs polynomial long division over the prime field.
    ///
    /// # Arguments
    /// * `divisor` - The polynomial to divide by.
    ///
    /// # Returns
    /// A tuple `(quotient, remainder)`.
    ///
    /// # Panics
    /// Panics if the divisor is the zero polynomial.
    pub fn long_division(self, divisor: &Self) -> Result<(Self, Self), String> {
        if divisor.coeffs.is_empty() || divisor.coeffs.iter().all(|c| c.value.is_zero()) {
            return Err("Division by zero polynomial".to_string());
        }
        let mut quotient =
            vec![PrimeFieldElement::new(Zero::zero(), self.field.clone()); self.coeffs.len()];
        let mut remainder = self.coeffs.clone();
        let divisor_deg = divisor.coeffs.len() - 1;
        let lead_divisor_inv = divisor.coeffs[0]
            .inverse()
            .ok_or_else(|| "Leading coefficient of divisor is not invertible.".to_string())?;
        while remainder.len() > divisor_deg && !remainder.is_empty() {
            let lead_rem = remainder[0].clone();
            let coeff = lead_rem * lead_divisor_inv.clone();
            let degree_diff = remainder.len() - divisor.coeffs.len();
            if degree_diff < quotient.len() {
                quotient[degree_diff] = coeff.clone();
            }
            for (i, vars) in remainder.iter_mut().enumerate().take(divisor_deg + 1) {
                let term = coeff.clone() * divisor.coeffs[i].clone();
                *vars = (*vars).clone() - term;
            }
            remainder.remove(0);
        }
        Ok((
            Self::new(quotient, self.field.clone()),
            Self::new(remainder, self.field),
        ))
    }
}
impl Add for FiniteFieldPolynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut result_coeffs =
            vec![PrimeFieldElement::new(Zero::zero(), self.field.clone()); max_len];
        let self_start = max_len - self.coeffs.len();
        result_coeffs[self_start..(self.coeffs.len() + self_start)]
            .clone_from_slice(&self.coeffs[..]);
        let rhs_start = max_len - rhs.coeffs.len();
        for i in 0..rhs.coeffs.len() {
            result_coeffs[rhs_start + i] =
                result_coeffs[rhs_start + i].clone() + rhs.coeffs[i].clone();
        }
        Self::new(result_coeffs, self.field)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for FiniteFieldPolynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let neg_rhs_coeffs = rhs.coeffs.into_iter().map(|c| -c).collect();
        let neg_rhs = Self::new(neg_rhs_coeffs, rhs.field);
        self + neg_rhs
    }
}
impl Mul for FiniteFieldPolynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Self::new(vec![], self.field);
        }
        let deg1 = self.coeffs.len() - 1;
        let deg2 = rhs.coeffs.len() - 1;
        let mut result_coeffs =
            vec![PrimeFieldElement::new(Zero::zero(), self.field.clone()); deg1 + deg2 + 1];
        for i in 0..=deg1 {
            for j in 0..=deg2 {
                let term_mul = self.coeffs[i].clone() * rhs.coeffs[j].clone();
                let existing = result_coeffs[i + j].clone();
                result_coeffs[i + j] = existing + term_mul;
            }
        }
        Self::new(result_coeffs, self.field)
    }
}
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ExtensionField {
    #[serde(with = "arc_serde")]
    pub prime_field: Arc<PrimeField>,
    pub irreducible_poly: FiniteFieldPolynomial,
}
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ExtensionFieldElement {
    pub poly: FiniteFieldPolynomial,
    #[serde(with = "arc_serde")]
    pub field: Arc<ExtensionField>,
}
impl ExtensionFieldElement {
    /// Creates a new element in an extension field.
    ///
    /// The element is represented by a polynomial, which is reduced modulo the
    /// field's irreducible polynomial to keep it in canonical form.
    ///
    /// # Arguments
    /// * `poly` - The `FiniteFieldPolynomial` representing the initial value of the element.
    /// * `field` - An `Arc` pointing to the `ExtensionField` this element belongs to.
    ///
    /// # Returns
    /// A new `ExtensionFieldElement`.
    #[must_use]
    pub fn new(poly: FiniteFieldPolynomial, field: Arc<ExtensionField>) -> Self {
        match poly.long_division(&field.irreducible_poly.clone()) {
            Ok((_, remainder)) => Self {
                poly: remainder,
                field,
            },
            Err(_) => Self {
                poly: FiniteFieldPolynomial::new(vec![], field.prime_field.clone()),
                field,
            },
        }
    }
    /// Computes the multiplicative inverse of the element in the extension field.
    ///
    /// This is done using the Extended Euclidean Algorithm for polynomials.
    ///
    /// # Returns
    /// * `Some(ExtensionFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible.
    #[must_use]
    pub fn inverse(&self) -> Option<Self> {
        let (gcd, _, inv) =
            poly_extended_gcd(self.poly.clone(), self.field.irreducible_poly.clone()).ok()?;
        if gcd.degree() > 0 || gcd.coeffs.is_empty() {
            return None;
        }
        let inv_factor = gcd.coeffs[0].inverse()?;
        Some(Self::new(
            inv * FiniteFieldPolynomial::new(vec![inv_factor], self.poly.field.clone()),
            self.field.clone(),
        ))
    }
}
pub(crate) fn poly_extended_gcd(
    a: FiniteFieldPolynomial,
    b: FiniteFieldPolynomial,
) -> Result<
    (
        FiniteFieldPolynomial,
        FiniteFieldPolynomial,
        FiniteFieldPolynomial,
    ),
    String,
> {
    let zero_poly = FiniteFieldPolynomial::new(vec![], a.field.clone());
    if b.coeffs.is_empty() || b.coeffs.iter().all(|c| c.value.is_zero()) {
        let one_poly = FiniteFieldPolynomial::new(
            vec![PrimeFieldElement::new(One::one(), a.field.clone())],
            a.field.clone(),
        );
        return Ok((a, one_poly, zero_poly));
    }
    let (q, r) = a.long_division(&b)?;
    let (g, x, y) = poly_extended_gcd(b, r)?;
    let t = x - (q * y.clone());
    Ok((g, y, t))
}
impl ExtensionFieldElement {
    pub fn add(self, rhs: Self) -> Result<Self, String> {
        Ok(Self::new(self.poly + rhs.poly, self.field))
    }
    pub fn sub(self, rhs: Self) -> Result<Self, String> {
        Ok(Self::new(self.poly - rhs.poly, self.field))
    }
    pub fn mul(self, rhs: Self) -> Result<Self, String> {
        Ok(Self::new(self.poly * rhs.poly, self.field))
    }
    pub fn div(self, rhs: &Self) -> Result<Self, String> {
        let inv_rhs = rhs
            .inverse()
            .ok_or_else(|| "Division by zero or non-invertible element.".to_string())?;
        self.mul(inv_rhs)
    }
}
impl Neg for ExtensionFieldElement {
    type Output = Self;
    fn neg(self) -> Self {
        let zero_poly = FiniteFieldPolynomial::new(vec![], self.poly.field.clone());
        Self::new(zero_poly - self.poly, self.field)
    }
}
