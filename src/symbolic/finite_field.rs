//! # Symbolic Finite Field Arithmetic
//!
//! This module provides symbolic structures for arithmetic in finite fields (Galois fields).
//! It defines prime fields GF(p) and extension fields GF(p^n), along with the necessary
//! arithmetic operations for their elements and for polynomials over these fields.

use crate::symbolic::number_theory::extended_gcd_inner;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;
// =====================================================================================
// region: Prime Field GF(p)
// =====================================================================================

#[derive(Debug, PartialEq, Eq, Clone)]
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
    pub fn new(modulus: BigInt) -> Arc<Self> {
        Arc::new(PrimeField { modulus })
    }
}

#[derive(Debug, Clone)]
pub struct PrimeFieldElement {
    pub value: BigInt,
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
    pub fn new(value: BigInt, field: Arc<PrimeField>) -> Self {
        let modulus = &field.modulus;
        let mut val = value % modulus;
        if val < Zero::zero() {
            val += modulus;
        }
        PrimeFieldElement { value: val, field }
    }

    /// Computes the multiplicative inverse of the element in the prime field.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod p)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(PrimeFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    pub fn inverse(&self) -> Option<Self> {
        let (g, x, _) = extended_gcd_inner(self.value.clone(), self.field.modulus.clone());
        if g.is_one() {
            let modulus = &self.field.modulus;
            let mut inv = x % modulus;
            if inv < Zero::zero() {
                inv += modulus;
            }
            Some(PrimeFieldElement::new(inv, self.field.clone()))
        } else {
            None
        }
    }
}

impl Add for PrimeFieldElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot add elements from different fields.");
        }
        let val = (self.value + rhs.value) % &self.field.modulus;
        PrimeFieldElement::new(val, self.field.clone())
    }
}

impl Sub for PrimeFieldElement {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot subtract elements from different fields.");
        }
        let val = (self.value - rhs.value + &self.field.modulus) % &self.field.modulus;
        PrimeFieldElement::new(val, self.field.clone())
    }
}

impl Mul for PrimeFieldElement {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.field != rhs.field {
            panic!("Cannot multiply elements from different fields.");
        }
        let val = (self.value * rhs.value) % &self.field.modulus;
        PrimeFieldElement::new(val, self.field.clone())
    }
}

impl Div for PrimeFieldElement {
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

impl Neg for PrimeFieldElement {
    type Output = Self;
    fn neg(self) -> Self {
        let val = (-self.value + &self.field.modulus) % &self.field.modulus;
        PrimeFieldElement::new(val, self.field.clone())
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
        // This is a bit of a hack, as we don't have a global field context.
        // A better solution would be for the Field trait to have a `zero(field)` method.
        // We assume a dummy field here, which is not ideal but works for now.
        let dummy_field = PrimeField::new(BigInt::from(2)); // Modulo 2 is arbitrary
        PrimeFieldElement::new(Zero::zero(), dummy_field)
    }
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}

impl One for PrimeFieldElement {
    fn one() -> Self {
        let dummy_field = PrimeField::new(BigInt::from(2));
        PrimeFieldElement::new(One::one(), dummy_field)
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

// =====================================================================================
// region: Polynomials over Prime Fields
// =====================================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteFieldPolynomial {
    pub coeffs: Vec<PrimeFieldElement>,
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
    pub fn new(coeffs: Vec<PrimeFieldElement>, field: Arc<PrimeField>) -> Self {
        let first_non_zero = coeffs
            .iter()
            .position(|c| !c.value.is_zero())
            .unwrap_or(coeffs.len());
        FiniteFieldPolynomial {
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
    pub fn degree(&self) -> isize {
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
    pub fn long_division(self, divisor: Self) -> (Self, Self) {
        if divisor.coeffs.is_empty() || divisor.coeffs.iter().all(|c| c.value.is_zero()) {
            panic!("Division by zero polynomial");
        }
        let mut quotient =
            vec![PrimeFieldElement::new(Zero::zero(), self.field.clone()); self.coeffs.len()];
        let mut remainder = self.coeffs.clone();
        let divisor_deg = divisor.coeffs.len() - 1;
        let lead_divisor_inv = divisor.coeffs[0].inverse().unwrap();

        while remainder.len() > divisor_deg && !remainder.is_empty() {
            let lead_rem = remainder[0].clone();
            let coeff = lead_rem * lead_divisor_inv.clone();
            let degree_diff = remainder.len() - divisor.coeffs.len();
            if degree_diff < quotient.len() {
                quotient[degree_diff] = coeff.clone();
            }

            for i in 0..=divisor_deg {
                //for (i, _item) in remainder.iter_mut().enumerate().take(divisor_deg + 1) {
                let term = coeff.clone() * divisor.coeffs[i].clone();
                remainder[i] = Sub::sub(remainder[i].clone(), term);
            }
            remainder.remove(0);
        }
        (
            FiniteFieldPolynomial::new(quotient, self.field.clone()),
            FiniteFieldPolynomial::new(remainder, self.field),
        )
    }
}

impl Add for FiniteFieldPolynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result_coeffs = vec![
            PrimeFieldElement::new(Zero::zero(), self.field.clone());
            std::cmp::max(self.coeffs.len(), rhs.coeffs.len())
        ];
        for (i, c) in self.coeffs.iter().rev().enumerate() {
            let new_idx = result_coeffs.len() - 1 - i;
            result_coeffs[new_idx] = Add::add(result_coeffs[new_idx].clone(), c.clone());
        }
        for (i, c) in rhs.coeffs.iter().rev().enumerate() {
            let new_idx = result_coeffs.len() - 1 - i;
            result_coeffs[new_idx] = Add::add(result_coeffs[new_idx].clone(), c.clone());
        }
        result_coeffs.reverse();
        FiniteFieldPolynomial::new(result_coeffs, self.field.clone())
    }
}

impl Sub for FiniteFieldPolynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let neg_rhs_coeffs = rhs.coeffs.into_iter().map(|c| -c).collect();
        let neg_rhs = FiniteFieldPolynomial::new(neg_rhs_coeffs, rhs.field);
        Add::add(self, neg_rhs)
    }
}

impl Mul for FiniteFieldPolynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return FiniteFieldPolynomial::new(vec![], self.field.clone());
        }
        let deg1 = self.coeffs.len() - 1;
        let deg2 = rhs.coeffs.len() - 1;
        let mut result_coeffs =
            vec![PrimeFieldElement::new(Zero::zero(), self.field.clone()); deg1 + deg2 + 1];

        for i in 0..=deg1 {
            for j in 0..=deg2 {
                let term_mul = Mul::mul(self.coeffs[i].clone(), rhs.coeffs[j].clone());
                let existing = result_coeffs[i + j].clone();
                result_coeffs[i + j] = Add::add(existing, term_mul);
            }
        }
        FiniteFieldPolynomial::new(result_coeffs, self.field.clone())
    }
}

// =====================================================================================
// region: Extension Fields GF(p^n)
// =====================================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtensionField {
    pub prime_field: Arc<PrimeField>,
    pub irreducible_poly: FiniteFieldPolynomial,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtensionFieldElement {
    pub poly: FiniteFieldPolynomial,
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
    pub fn new(poly: FiniteFieldPolynomial, field: Arc<ExtensionField>) -> Self {
        let (_, remainder) = poly.long_division(field.irreducible_poly.clone());
        ExtensionFieldElement {
            poly: remainder,
            field,
        }
    }

    /// Computes the multiplicative inverse of the element in the extension field.
    ///
    /// This is done using the Extended Euclidean Algorithm for polynomials.
    ///
    /// # Returns
    /// * `Some(ExtensionFieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible.
    pub fn inverse(&self) -> Option<Self> {
        let (gcd, _, inv) =
            poly_extended_gcd(self.poly.clone(), self.field.irreducible_poly.clone());
        if gcd.degree() > 0 || gcd.coeffs.is_empty() {
            return None;
        }
        let inv_factor = gcd.coeffs[0].inverse()?;
        Some(ExtensionFieldElement::new(
            Mul::mul(
                inv,
                FiniteFieldPolynomial::new(vec![inv_factor], self.poly.field.clone()),
            ),
            self.field.clone(),
        ))
    }
}

pub(crate) fn poly_extended_gcd(
    a: FiniteFieldPolynomial,
    b: FiniteFieldPolynomial,
) -> (
    FiniteFieldPolynomial,
    FiniteFieldPolynomial,
    FiniteFieldPolynomial,
) {
    let zero_poly = FiniteFieldPolynomial::new(vec![], a.field.clone());
    if b.coeffs.is_empty() || b.coeffs.iter().all(|c| c.value.is_zero()) {
        let one_poly = FiniteFieldPolynomial::new(
            vec![PrimeFieldElement::new(One::one(), a.field.clone())],
            a.field.clone(),
        );
        return (a, one_poly, zero_poly);
    }

    let (q, r) = a.clone().long_division(b.clone());
    let (g, x, y) = poly_extended_gcd(b, r);

    let t = Sub::sub(x, Mul::mul(q, y.clone()));
    (g, y, t)
}

impl Add for ExtensionFieldElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        ExtensionFieldElement::new(Add::add(self.poly, rhs.poly), self.field.clone())
    }
}

impl Sub for ExtensionFieldElement {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        ExtensionFieldElement::new(Sub::sub(self.poly, rhs.poly), self.field.clone())
    }
}

impl Mul for ExtensionFieldElement {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        ExtensionFieldElement::new(Mul::mul(self.poly, rhs.poly), self.field.clone())
    }
}

impl Div for ExtensionFieldElement {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let inv_rhs = rhs
            .inverse()
            .expect("Division by zero or non-invertible element.");
        //self * inv_rhs
        self.mul(inv_rhs)
    }
}

impl Neg for ExtensionFieldElement {
    type Output = Self;
    fn neg(self) -> Self {
        let zero_poly = FiniteFieldPolynomial::new(vec![], self.poly.field.clone());
        ExtensionFieldElement::new(Sub::sub(zero_poly, self.poly), self.field.clone())
    }
}
