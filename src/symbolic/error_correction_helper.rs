//! # Finite Fields (Galois Fields)
//!
//! This module provides structures and functions for arithmetic in finite fields.
//! It is a foundational component for advanced algebra, cryptography, and error-correcting codes.
//!
//! ## Features
//! - General finite field arithmetic (`FieldElement`)
//! - GF(2^8) with lookup table optimizations
//! - Polynomial operations over finite fields
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use crate::symbolic::core::Expr;
use crate::symbolic::number_theory::extended_gcd;

/// Represents a finite field GF(p) where p is the modulus.
#[derive(Debug, PartialEq, Eq)]

pub struct FiniteField {
    /// The modulus (characteristic) of the field
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
    #[must_use]

    pub fn new(
        modulus: i64
    ) -> Arc<Self> {

        Arc::new(Self {
            modulus: BigInt::from(
                modulus,
            ),
        })
    }

    /// Creates a finite field from a `BigInt` modulus.
    #[must_use]

    pub fn from_bigint(
        modulus: BigInt
    ) -> Arc<Self> {

        Arc::new(Self {
            modulus,
        })
    }
}

/// Represents an element in a finite field.
#[derive(Debug, Clone)]

pub struct FieldElement {
    /// The value of the element (reduced modulo the field characteristic)
    pub value: BigInt,
    /// The field this element belongs to
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
    #[must_use]

    pub fn new(
        value: BigInt,
        field: Arc<FiniteField>,
    ) -> Self {

        let reduced = ((value
            % &field.modulus)
            + &field.modulus)
            % &field.modulus;

        Self {
            value: reduced,
            field,
        }
    }

    /// Returns true if this element is zero.
    #[must_use]

    pub fn is_zero(&self) -> bool {

        self.value.is_zero()
    }

    /// Returns true if this element is one.
    #[must_use]

    pub fn is_one(&self) -> bool {

        self.value.is_one()
    }

    /// Computes the multiplicative inverse of the element in the finite field.
    ///
    /// The inverse `x` of an element `a` is such that `a * x = 1 (mod modulus)`.
    /// This implementation uses the Extended Euclidean Algorithm.
    ///
    /// # Returns
    /// * `Some(FieldElement)` containing the inverse if it exists.
    /// * `None` if the element is not invertible (i.e., its value is not coprime to the modulus).
    #[must_use]

    pub fn inverse(
        &self
    ) -> Option<Self> {

        let (g, x, _) = extended_gcd(
            &Expr::BigInt(
                self.value.clone(),
            ),
            &Expr::BigInt(
                self.field
                    .modulus
                    .clone(),
            ),
        );

        if let Expr::BigInt(g_val) = g {

            if g_val.is_one() {

                let inv = x
                    .to_bigint()
                    .unwrap_or_default(
                    );

                let modulus =
                    &self.field.modulus;

                return Some(
                    Self::new(
                        (inv % modulus
                            + modulus)
                            % modulus,
                        self.field
                            .clone(),
                    ),
                );
            }
        }

        None
    }

    /// Computes a^exp mod p using binary exponentiation.
    ///
    /// # Arguments
    /// * `exp` - The exponent (non-negative)
    ///
    /// # Returns
    /// A new `FieldElement` representing self^exp
    #[must_use]

    pub fn pow(
        &self,
        exp: u64,
    ) -> Self {

        if exp == 0 {

            return Self::new(
                BigInt::one(),
                self.field.clone(),
            );
        }

        let mut result = Self::new(
            BigInt::one(),
            self.field.clone(),
        );

        let mut base = self.clone();

        let mut e = exp;

        while e > 0 {

            if e & 1 == 1 {

                result = (result
                    * base.clone())
                .unwrap_or_else(|_| {

                    Self::new(
                        BigInt::zero(),
                        self.field
                            .clone(),
                    )
                });
            }

            base = (base.clone()
                * base.clone())
            .unwrap_or_else(|_| {

                Self::new(
                    BigInt::zero(),
                    self.field.clone(),
                )
            });

            e >>= 1;
        }

        result
    }
}

impl Add for FieldElement {
    type Output = Result<Self, String>;
/// # Errors
///
/// This function will return an error if `self` and `rhs` belong to different
/// finite fields.


    fn add(
        self,
        rhs: Self,
    ) -> Self::Output {

        if self.field != rhs.field {

            return Err(
                "Cannot add elements \
                 from different \
                 fields."
                    .to_string(),
            );
        }

        let val = (self.value
            + rhs.value)
            % &self.field.modulus;

        Ok(Self::new(
            val,
            self.field,
        ))
    }
}

impl Sub for FieldElement {
    type Output = Result<Self, String>;
/// # Errors
///
/// This function will return an error if `self` and `rhs` belong to different
/// finite fields.


    fn sub(
        self,
        rhs: Self,
    ) -> Self::Output {

        if self.field != rhs.field {

            return Err(
                "Cannot subtract \
                 elements from \
                 different fields."
                    .to_string(),
            );
        }

        let val = (self.value
            - rhs.value
            + &self.field.modulus)
            % &self.field.modulus;

        Ok(Self::new(
            val,
            self.field,
        ))
    }
}

impl Mul for FieldElement {
    type Output = Result<Self, String>;
/// # Errors
///
/// This function will return an error if `self` and `rhs` belong to different
/// finite fields.


    fn mul(
        self,
        rhs: Self,
    ) -> Self::Output {

        if self.field != rhs.field {

            return Err(
                "Cannot multiply \
                 elements from \
                 different fields."
                    .to_string(),
            );
        }

        let val = (self.value
            * rhs.value)
            % &self.field.modulus;

        Ok(Self::new(
            val,
            self.field,
        ))
    }
}

impl Div for FieldElement {
    type Output = Result<Self, String>;

    fn div(
        self,
        rhs: Self,
    ) -> Self::Output {

        if self.field != rhs.field {

            return Err(
                "Cannot divide \
                 elements from \
                 different fields."
                    .to_string(),
            );
        }

        let inv_rhs = rhs
            .inverse()
            .ok_or_else(|| {

                "Division by zero or \
                 non-invertible \
                 element."
                    .to_string()
            })?;

        self.mul(inv_rhs)
    }
}

impl Neg for FieldElement {
    type Output = Self;

    fn neg(self) -> Self {

        let val = (-self.value
            + &self.field.modulus)
            % &self.field.modulus;

        Self::new(val, self.field)
    }
}

impl PartialEq for FieldElement {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {

        self.field == other.field
            && self.value == other.value
    }
}

impl Eq for FieldElement {
}

const GF256_GENERATOR_POLY: u16 = 0x11d;

const GF256_MODULUS: usize = 256;

struct Gf256Tables {
    log: [u8; GF256_MODULUS],
    exp: [u8; GF256_MODULUS],
}

static GF256_TABLES:
    std::sync::LazyLock<Gf256Tables> =
    std::sync::LazyLock::new(|| {

        let mut log_table =
            [0u8; GF256_MODULUS];

        let mut exp_table =
            [0u8; GF256_MODULUS];

        let mut x: u16 = 1;

        for i in 0 .. 255 {

            exp_table[i] = x as u8;

            log_table[x as usize] =
                i as u8;

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

/// Computes the exponentiation (anti-logarithm) in GF(2^8).
///
/// Returns `alpha^log_val` where alpha is the primitive element.
///
/// # Arguments
/// * `log_val` - The logarithm of the element.
///
/// # Returns
/// The field element `alpha ^ log_val`.
#[must_use]

pub fn gf256_exp(log_val: u8) -> u8 {

    GF256_TABLES.exp[log_val as usize]
}

/// Computes the discrete logarithm in GF(2^8).
///
/// Returns `log_alpha(a)` where alpha is the primitive element.
///
/// # Arguments
/// * `a` - The field element (must be non-zero)
///
/// # Returns
/// `Ok(log)` if a is non-zero, `Err` if a is zero

pub fn gf256_log(
    a: u8
) -> Result<u8, String> {

    if a == 0 {

        return Err(
            "Logarithm of zero is \
             undefined"
                .to_string(),
        );
    }

    Ok(GF256_TABLES.log[a as usize])
}

/// Performs addition in the finite field GF(2^8).
///
/// In fields of characteristic 2, addition is equivalent to a bitwise XOR operation.
#[inline]
#[must_use]

pub const fn gf256_add(
    a: u8,
    b: u8,
) -> u8 {

    a ^ b
}

/// Performs multiplication in GF(2^8) using precomputed lookup tables.
///
/// Multiplication is performed by adding the logarithms of the operands and then
/// finding the anti-logarithm of the result.
#[inline]
#[must_use]

pub fn gf256_mul(
    a: u8,
    b: u8,
) -> u8 {

    if a == 0 || b == 0 {

        0
    } else {

        let log_a = u16::from(
            GF256_TABLES.log
                [a as usize],
        );

        let log_b = u16::from(
            GF256_TABLES.log
                [b as usize],
        );

        GF256_TABLES.exp[((log_a
            + log_b)
            % 255)
            as usize]
    }
}

/// Computes the multiplicative inverse of an element in GF(2^8).
///
/// The inverse is calculated using the logarithm and exponentiation tables.
///
/// # Errors
///
/// This function will return an error if `a` is 0, as 0 has no multiplicative inverse.
#[inline]

pub fn gf256_inv(
    a: u8
) -> Result<u8, String> {

    if a == 0 {

        return Err("Cannot invert 0"
            .to_string());
    }

    Ok(
        GF256_TABLES.exp[(255
            - u16::from(
                GF256_TABLES.log
                    [a as usize],
            ))
            as usize],
    )
}

/// Performs division in GF(2^8).
///
/// Division is implemented as multiplication by the multiplicative inverse of the divisor.
///
/// # Errors
///
/// This function will return an error if `b` (the divisor) is 0, as division by zero is undefined.
#[inline]

pub fn gf256_div(
    a: u8,
    b: u8,
) -> Result<u8, String> {

    if b == 0 {

        return Err("Division by zero"
            .to_string());
    }

    if a == 0 {

        return Ok(0);
    }

    let log_a = u16::from(
        GF256_TABLES.log[a as usize],
    );

    let log_b = u16::from(
        GF256_TABLES.log[b as usize],
    );

    Ok(
        GF256_TABLES.exp[((log_a + 255
            - log_b)
            % 255)
            as usize],
    )
}

/// Computes a^exp in GF(2^8).
///
/// Uses logarithm tables for efficiency.
///
/// # Arguments
/// * `a` - The base element
/// * `exp` - The exponent
///
/// # Returns
/// a^exp in GF(2^8)
#[must_use]

pub fn gf256_pow(
    a: u8,
    exp: u8,
) -> u8 {

    if a == 0 {

        return u8::from(exp == 0);
    }

    if exp == 0 {

        return 1;
    }

    let log_a = u16::from(
        GF256_TABLES.log[a as usize],
    );

    let log_result =
        (log_a * u16::from(exp)) % 255;

    GF256_TABLES.exp
        [log_result as usize]
}

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
#[must_use]

pub fn poly_eval_gf256(
    poly: &[u8],
    x: u8,
) -> u8 {

    let mut y = 0;

    for coeff in poly {

        y = gf256_mul(y, x) ^ coeff;
    }

    y
}

/// Adds two polynomials over GF(2^8).
///
/// Polynomial addition in GF(2^8) is performed by `XORing` corresponding coefficients.
///
/// # Arguments
/// * `p1` - The first polynomial as a slice of `u8` coefficients.
/// * `p2` - The second polynomial as a slice of `u8` coefficients.
///
/// # Returns
/// A `Vec<u8>` representing the sum polynomial.
#[must_use]

pub fn poly_add_gf256(
    p1: &[u8],
    p2: &[u8],
) -> Vec<u8> {

    let mut result = vec![
        0;
        std::cmp::max(
            p1.len(),
            p2.len()
        )
    ];

    let res_len = result.len();

    for i in 0 .. p1.len() {

        result
            [i + res_len - p1.len()] =
            p1[i];
    }

    for i in 0 .. p2.len() {

        result
            [i + res_len - p2.len()] ^=
            p2[i];
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
#[must_use]

pub fn poly_mul_gf256(
    p1: &[u8],
    p2: &[u8],
) -> Vec<u8> {

    if p1.is_empty() || p2.is_empty() {

        return vec![];
    }

    let mut result =
        vec![
            0;
            p1.len() + p2.len() - 1
        ];

    for i in 0 .. p1.len() {

        for j in 0 .. p2.len() {

            result[i + j] ^=
                gf256_mul(p1[i], p2[j]);
        }
    }

    result
}

/// Scales a polynomial by a constant in GF(2^8).
///
/// # Arguments
/// * `poly` - The polynomial coefficients
/// * `scalar` - The scalar to multiply by
///
/// # Returns
/// A new polynomial with scaled coefficients
#[must_use]

pub fn poly_scale_gf256(
    poly: &[u8],
    scalar: u8,
) -> Vec<u8> {

    poly.iter()
        .map(|&c| gf256_mul(c, scalar))
        .collect()
}

/// Computes the formal derivative of a polynomial in GF(2^8).
///
/// For a polynomial `a_n`*x^n + ... + `a_1`*x + `a_0`, the derivative is
/// n*`a_n`*x^(n-1) + ... + `a_1`. In characteristic 2, even-indexed terms vanish.
///
/// # Arguments
/// * `poly` - The polynomial coefficients (highest degree first)
///
/// # Returns
/// The derivative polynomial
#[must_use]

pub fn poly_derivative_gf256(
    poly: &[u8]
) -> Vec<u8> {

    if poly.len() <= 1 {

        return vec![0];
    }

    let n = poly.len() - 1; // degree
    let mut result =
        Vec::with_capacity(n);

    for (i, &coeff) in poly
        .iter()
        .enumerate()
    {

        let power = n - i;

        if power > 0 {

            // In GF(2^8), multiply by power. Only odd powers survive in char 2.
            if power % 2 == 1 {

                result.push(coeff);
            } else {

                result.push(0);
            }
        }
    }

    // Remove leading zeros
    while result.len() > 1
        && result[0] == 0
    {

        result.remove(0);
    }

    result
}

/// Computes the GCD of two polynomials over GF(2^8) using Euclidean algorithm.
///
/// # Arguments
/// * `p1` - First polynomial
/// * `p2` - Second polynomial
///
/// # Returns
/// The GCD polynomial (monic)
#[must_use]

pub fn poly_gcd_gf256(
    p1: &[u8],
    p2: &[u8],
) -> Vec<u8> {

    // Remove leading zeros
    let strip_leading =
        |p: &[u8]| -> Vec<u8> {

            let first_non_zero = p
                .iter()
                .position(|&x| x != 0)
                .unwrap_or(p.len());

            if first_non_zero >= p.len()
            {

                vec![0]
            } else {

                p[first_non_zero ..]
                    .to_vec()
            }
        };

    let mut a = strip_leading(p1);

    let mut b = strip_leading(p2);

    while b.len() > 1
        || (b.len() == 1 && b[0] != 0)
    {

        if b.is_empty()
            || (b.len() == 1
                && b[0] == 0)
        {

            break;
        }

        // Compute remainder of a / b
        let remainder =
            match poly_div_gf256(
                a.clone(),
                &b,
            ) {
                | Ok(r) => {
                    strip_leading(&r)
                },
                | Err(_) => break,
            };

        a = b;

        b = remainder;
    }

    // Make monic (leading coefficient = 1)
    if !a.is_empty() && a[0] != 0 {

        if let Ok(inv) = gf256_inv(a[0])
        {

            a = poly_scale_gf256(
                &a, inv,
            );
        }
    }

    a
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

pub fn poly_div_gf256(
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

    let lead_divisor_inv =
        gf256_inv(lead_divisor)?;

    while dividend.len() >= divisor_len
    {

        let lead_dividend = dividend[0];

        let coeff = gf256_mul(
            lead_dividend,
            lead_divisor_inv,
        );

        for i in 0 .. divisor_len {

            let term = gf256_mul(
                coeff,
                divisor[i],
            );

            dividend[i] ^= term;
        }

        dividend.remove(0);
    }

    Ok(dividend)
}

pub(crate) fn expr_to_field_elements(
    p_expr: &Expr,
    field: &Arc<FiniteField>,
) -> Result<Vec<FieldElement>, String> {

    if let Expr::Polynomial(coeffs) =
        p_expr
    {

        coeffs
            .iter()
            .map(|c| {

                c.to_bigint()
                    .map(|val| FieldElement::new(val, field.clone()))
                    .ok_or_else(|| format!("Invalid coefficient in polynomial: {c}"))
            })
            .collect()
    } else {

        Err(format!(
            "Expression is not a \
             polynomial: {p_expr}"
        ))
    }
}

pub(crate) fn field_elements_to_expr(
    coeffs: &[FieldElement]
) -> Expr {

    let expr_coeffs = coeffs
        .iter()
        .map(|c| {

            Expr::BigInt(
                c.value.clone(),
            )
        })
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
/// * `Err(String)` if the input expressions are not valid polynomials, contain invalid coefficients,
///   or if the underlying `FieldElement::add` operation fails due to field mismatch.

pub fn poly_add_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: &Arc<FiniteField>,
) -> Result<Expr, String> {

    let c1 = expr_to_field_elements(
        p1_expr,
        field,
    )?;

    let c2 = expr_to_field_elements(
        p2_expr,
        field,
    )?;

    let mut result_coeffs = vec![];

    let len1 = c1.len();

    let len2 = c2.len();

    let max_len =
        std::cmp::max(len1, len2);

    for i in 0 .. max_len {

        let val1 = if i < len1 {

            c1[len1 - 1 - i].clone()
        } else {

            FieldElement::new(
                Zero::zero(),
                field.clone(),
            )
        };

        let val2 = if i < len2 {

            c2[len2 - 1 - i].clone()
        } else {

            FieldElement::new(
                Zero::zero(),
                field.clone(),
            )
        };

        result_coeffs
            .push((val1 + val2)?);
    }

    result_coeffs.reverse();

    Ok(
        field_elements_to_expr(
            &result_coeffs,
        ),
    )
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
/// * `Err(String)` if the input expressions are not valid polynomials, contain invalid coefficients,
///   or if the underlying `FieldElement::mul` operation fails due to field mismatch.

pub fn poly_mul_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: &Arc<FiniteField>,
) -> Result<Expr, String> {

    let c1 = expr_to_field_elements(
        p1_expr,
        field,
    )?;

    let c2 = expr_to_field_elements(
        p2_expr,
        field,
    )?;

    if c1.is_empty() || c2.is_empty() {

        return Ok(Expr::Polynomial(
            vec![],
        ));
    }

    let deg1 = c1.len() - 1;

    let deg2 = c2.len() - 1;

    let mut result_coeffs = vec![
        FieldElement::new(
            Zero::zero(),
            field.clone()
        );
        deg1 + deg2 + 1
    ];

    for i in 0 ..= deg1 {

        for j in 0 ..= deg2 {

            let term_mul = (c1[i]
                .clone()
                * c2[j].clone())?;

            result_coeffs[i + j] =
                (result_coeffs[i + j]
                    .clone()
                    + term_mul)?;
        }
    }

    Ok(
        field_elements_to_expr(
            &result_coeffs,
        ),
    )
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
///   if division by the zero polynomial is attempted, or if a leading coefficient is not invertible.

pub fn poly_div_gf(
    p1_expr: &Expr,
    p2_expr: &Expr,
    field: &Arc<FiniteField>,
) -> Result<(Expr, Expr), String> {

    let mut num =
        expr_to_field_elements(
            p1_expr,
            field,
        )?;

    let den = expr_to_field_elements(
        p2_expr,
        field,
    )?;

    if den
        .iter()
        .all(|c| c.value.is_zero())
    {

        return Err("Division by \
                    zero polynomial"
            .to_string());
    }

    let mut quotient = vec![
        FieldElement::new(
            Zero::zero(),
            field.clone()
        );
        num.len()
    ];

    let lead_den_inv = den
        .first()
        .ok_or(
            "Divisor polynomial is \
             empty."
                .to_string(),
        )?
        .inverse()
        .ok_or(
            "Leading coefficient is \
             not invertible"
                .to_string(),
        )?;

    while num.len() >= den.len() {

        let lead_num = match num.first()
        {
            | Some(n) => n.clone(),
            | None => return Err(
                "Dividend became \
                 empty unexpectedly."
                    .to_string(),
            ),
        };

        let coeff = (lead_num
            * lead_den_inv.clone())?;

        let degree_diff =
            num.len() - den.len();

        quotient[degree_diff] =
            coeff.clone();

        for (i, den_coeff) in den
            .iter()
            .enumerate()
        {

            let term = (coeff.clone()
                * den_coeff.clone())?;

            num[i] = (num[i].clone()
                - term)?;
        }

        num.remove(0);
    }

    let first_non_zero = num
        .iter()
        .position(|c| {

            !c.value.is_zero()
        })
        .unwrap_or(num.len());

    let remainder =
        &num[first_non_zero ..];

    Ok((
        field_elements_to_expr(
            &quotient,
        ),
        field_elements_to_expr(
            remainder,
        ),
    ))
}

trait ToBigInt {
    fn to_bigint(
        &self
    ) -> Option<BigInt>;
}

impl ToBigInt for Expr {
    fn to_bigint(
        &self
    ) -> Option<BigInt> {

        match self {
            | Self::BigInt(i) => {
                Some(i.clone())
            },
            | Self::Constant(_) => None,
            | _ => None,
        }
    }
}
