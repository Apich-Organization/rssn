//! # Polynomial Factorization over Finite Fields
//!
//! This module provides implementations of algorithms for factoring polynomials
//! over finite fields. It includes Berlekamp's algorithm for small fields,
//! Cantor-Zassenhaus algorithm for larger fields, and square-free factorization.
//! It also contains a simplified approach to Berlekamp-Zassenhaus for integer polynomials.

use crate::numerical::matrix::Matrix;
use crate::symbolic::finite_field::{FiniteFieldPolynomial, PrimeField, PrimeFieldElement};
// Note: This dependency would need to be added to Cargo.toml
use itertools::Itertools;
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use rand;
use std::ops::Mul;
use std::ops::Sub;
use std::sync::Arc;

// =====================================================================================
// region: Main Factorization Dispatcher
// =====================================================================================

/// Factors a polynomial over a finite field.
///
/// This function acts as a dispatcher, choosing between Berlekamp's algorithm
/// (for small fields) and Cantor-Zassenhaus algorithm (for larger fields)
/// based on the size of the field's modulus.
///
/// # Arguments
/// * `poly` - The polynomial to factor, represented as a `FiniteFieldPolynomial`.
///
/// # Returns
/// A `Vec<FiniteFieldPolynomial>` containing the irreducible factors of the polynomial.
pub fn factor_gf(poly: &FiniteFieldPolynomial) -> Vec<FiniteFieldPolynomial> {
    // A simple heuristic: use Berlekamp for small fields, Cantor-Zassenhaus for large fields.
    if poly.field.modulus.to_u64().unwrap_or(u64::MAX) < 50 {
        berlekamp_factorization(poly)
    } else {
        cantor_zassenhaus(poly)
    }
}

/// Computes the derivative of a polynomial over a prime field.
///
/// This function applies the standard power rule for differentiation to each term
/// of the polynomial, with arithmetic performed modulo the field's characteristic.
///
/// # Arguments
/// * `p` - The polynomial to differentiate.
///
/// # Returns
/// A new `FiniteFieldPolynomial` representing the derivative.
pub fn poly_derivative_gf(p: &FiniteFieldPolynomial) -> FiniteFieldPolynomial {
    if p.coeffs.is_empty() {
        return FiniteFieldPolynomial::new(vec![], p.field.clone());
    }
    let mut deriv_coeffs = Vec::with_capacity(p.coeffs.len().saturating_sub(1));
    let degree = p.degree() as usize;
    for i in 0..degree {
        let original_coeff = &p.coeffs[i];
        let power = degree - i;
        let new_coeff_val = &original_coeff.value * BigInt::from(power);
        deriv_coeffs.push(PrimeFieldElement::new(new_coeff_val, p.field.clone()));
    }
    FiniteFieldPolynomial::new(deriv_coeffs, p.field.clone())
}

// =====================================================================================
// region: Berlekamp & Berlekamp-Zassenhaus
// =====================================================================================

/// Performs square-free factorization of a polynomial over a prime field.
///
/// This algorithm decomposes a polynomial `f(x)` into a product of powers of distinct
/// irreducible polynomials: `f(x) = f_1(x) * f_2(x)^2 * ... * f_k(x)^k`.
/// It is a preprocessing step for many factorization algorithms.
///
/// # Arguments
/// * `f` - The polynomial to factorize.
///
/// # Returns
/// A `Vec<(FiniteFieldPolynomial, usize)>` where each tuple contains a square-free
/// polynomial and its multiplicity.
pub fn square_free_factorization_gf(
    f: FiniteFieldPolynomial,
) -> Vec<(FiniteFieldPolynomial, usize)> {
    // Implementation based on Yun's algorithm for square-free factorization.
    let mut factors = Vec::new();
    let mut i = 1;
    let mut f_i = f;
    while f_i.degree() > 0 {
        let f_prime = poly_derivative_gf(&f_i);
        let g = poly_gcd_gf(f_i.clone(), f_prime);
        let h = f_i.clone().long_division(g.clone()).0;
        if h.degree() > 0 {
            factors.push((h, i));
        }
        f_i = g;
        i += 1;
    }
    factors
}

/*
/// Computes base^exp mod modulus for polynomials over a prime field.
pub(crate) fn poly_pow_mod(
    base: FiniteFieldPolynomial,
    exp: &BigInt,
    modulus: &FiniteFieldPolynomial,
) -> FiniteFieldPolynomial {
    let mut res = FiniteFieldPolynomial::new(
        vec![PrimeFieldElement::new(One::one(), base.field.clone())],
        base.field.clone(),
    );
    let mut b = base;
    let mut e = exp.clone();
    while e > Zero::zero() {
        if &e % 2 == One::one() {
            res = res.mul(b.clone()).long_division(modulus.clone()).1;
        }
        b = b.clone().mul(b.clone()).long_division(modulus.clone()).1;
        e >>= 1;
    }
    res
}
*/

/// Factors a square-free polynomial over a small prime field using Berlekamp's algorithm.

///

/// Berlekamp's algorithm is a classical method for factoring polynomials over finite fields.

/// It constructs a matrix (Berlekamp matrix) whose null space provides information about

/// the factors. The algorithm then uses GCD computations to split the polynomial.

///

/// # Arguments

/// * `f` - The square-free polynomial to factor.

///

/// # Returns

/// A `Vec<FiniteFieldPolynomial>` containing the irreducible factors.

pub fn berlekamp_factorization(f: &FiniteFieldPolynomial) -> Vec<FiniteFieldPolynomial> {
    let p_val = f.field.modulus.to_u64().unwrap();
    let n = f.degree() as usize;
    if n <= 1 {
        return vec![f.clone()];
    }

    let mut q_data = Vec::new();
    let x_poly = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(One::one(), f.field.clone()),
            PrimeFieldElement::new(Zero::zero(), f.field.clone()),
        ],
        f.field.clone(),
    );

    for i in 0..n {
        let exp = BigInt::from(p_val).pow(i as u32);
        let x_pow_mod_f = poly_pow_mod(x_poly.clone(), &exp, f);
        let mut row = vec![PrimeFieldElement::new(Zero::zero(), f.field.clone()); n];
        let offset = n.saturating_sub(x_pow_mod_f.coeffs.len());
        for (j, coeff) in x_pow_mod_f.coeffs.iter().enumerate() {
            row[offset + j] = coeff.clone();
        }
        q_data.extend(row);
    }
    let mut q_matrix = Matrix::new(n, n, q_data);

    for i in 0..n {
        let val = q_matrix.get(i, i).clone();
        *q_matrix.get_mut(i, i) = val - PrimeFieldElement::new(One::one(), f.field.clone());
    }

    let null_space_matrix = q_matrix.null_space();
    let basis_vectors = null_space_matrix.get_cols();
    let r = basis_vectors.len();

    if r == 1 {
        return vec![f.clone()];
    }

    let mut factors = vec![f.clone()];
    for v_coeffs in basis_vectors.iter().skip(1) {
        let v = FiniteFieldPolynomial::new(v_coeffs.clone(), f.field.clone());
        let mut new_factors = Vec::new();
        for s in 0..p_val {
            let s_elem = PrimeFieldElement::new(BigInt::from(s), f.field.clone());
            let v_minus_s = v.clone() - FiniteFieldPolynomial::new(vec![s_elem], f.field.clone());

            for factor in factors.iter() {
                let h = poly_gcd_gf(factor.clone(), v_minus_s.clone());
                if h.degree() > 0 && h.degree() < factor.degree() {
                    new_factors.push(h.clone());
                    new_factors.push(factor.clone().long_division(h).0);
                } else {
                    new_factors.push(factor.clone());
                }
            }
            factors = new_factors;
            new_factors = Vec::new();
            if factors.len() == r {
                break;
            }
        }
        if factors.len() == r {
            break;
        }
    }
    factors
}

// ... (Berlekamp-Zassenhaus and Hensel lifting would go here) ...

/// Factors a polynomial with integer coefficients using the Berlekamp-Zassenhaus algorithm.
///
/// This algorithm combines modular factorization (using Berlekamp's algorithm over `GF(p)`),
/// Hensel lifting (to lift factors from `GF(p)` to `GF(p^k)`), and recombination techniques
/// to find factors over the integers.
///
/// **Note**: This is a simplified implementation. A full implementation would handle
/// content, leading coefficients, square-free factorization, and more robust Hensel lifting.
///
/// # Arguments
/// * `poly` - The polynomial to factor, assumed to have integer coefficients.
///
/// # Returns
/// A `Vec<FiniteFieldPolynomial>` containing the factors over the integers.
pub fn berlekamp_zassenhaus(poly: &FiniteFieldPolynomial) -> Vec<FiniteFieldPolynomial> {
    // For simplicity, assume poly is monic and square-free.
    // A full implementation would handle content, leading coefficients, and square-free factorization.
    // Stage 1: Modular Factorization
    // Choose a prime p that does not divide the leading coefficient (already assumed to be 1).
    // A robust implementation would have a list of primes and check conditions.
    let p = BigInt::from(5);
    let field = PrimeField::new(p.clone());
    let f_mod_p = poly_with_field(poly, field);
    let factors_mod_p = berlekamp_factorization(&f_mod_p);
    if factors_mod_p.len() <= 1 {
        return vec![poly.clone()]; // Already irreducible mod p, likely irreducible over Z.
    }
    // Stage 2: Hensel Lifting
    // Lift the factorization f = g*h mod p to f = g_k*h_k mod p^k
    // For now, we lift the first factor against the product of the rest.
    let g_mod_p = factors_mod_p[0].clone();
    let h_mod_p = factors_mod_p.iter().skip(1).fold(
        FiniteFieldPolynomial::new(
            vec![PrimeFieldElement::new(One::one(), f_mod_p.field.clone())],
            f_mod_p.field.clone(),
        ),
        |acc, factor| acc * factor.clone(),
    );
    // We need to lift to a bound p^k > 2 * B, where B is Mignotte's bound on coefficients.
    // For simplicity, we lift to a fixed power, e.g., p^4.
    let k = 4;
    let (g_lifted, _h_lifted) = match hensel_lift(poly, &g_mod_p, &h_mod_p, &p, k) {
        Some((g, h)) => (g, h),
        None => return vec![poly.clone()], // Lifting failed
    };
    // Stage 3: Recombination of Factors
    // We have a set of factors mod p^k. We need to find which subsets multiply to true factors over Z.
    // This is the most complex part. We will try combinations.
    let mut true_factors = Vec::new();
    let mut remaining_poly = poly.clone();
    let lifted_factors = [g_lifted]; // In a full impl, this would be all lifted factors.
    for i in 1..=lifted_factors.len() {
        for subset in lifted_factors.iter().combinations(i) {
            let mut potential_factor = FiniteFieldPolynomial::new(
                vec![PrimeFieldElement::new(One::one(), poly.field.clone())],
                poly.field.clone(),
            );
            for factor in subset {
                potential_factor = potential_factor * factor.clone();
            }
            // Center the coefficients of the potential factor around 0.
            let p_k = p.pow(k);
            let p_k_half = &p_k / 2;
            let centered_coeffs = potential_factor
                .coeffs
                .into_iter()
                .map(|c| {
                    let mut val = c.value;
                    if val > p_k_half {
                        val -= &p_k;
                    }
                    PrimeFieldElement::new(val, poly.field.clone())
                })
                .collect();
            let centered_factor = FiniteFieldPolynomial::new(centered_coeffs, poly.field.clone());
            // Trial division
            let (quotient, remainder) = remaining_poly
                .clone()
                .long_division(centered_factor.clone());
            if remainder.coeffs.is_empty() || remainder.coeffs.iter().all(|c| c.value.is_zero()) {
                true_factors.push(centered_factor);
                remaining_poly = quotient;
                // A full implementation would remove used factors and restart combinations.
            }
        }
    }
    if !remaining_poly.coeffs.is_empty() {
        true_factors.push(remaining_poly);
    }
    true_factors
}
/// Lifts a factorization f ≡ g*h (mod p) to a factorization f ≡ g_k*h_k (mod p^k).
pub(crate) fn hensel_lift(
    f: &FiniteFieldPolynomial,
    g: &FiniteFieldPolynomial,
    h: &FiniteFieldPolynomial,
    p: &BigInt,
    k: u32,
) -> Option<(FiniteFieldPolynomial, FiniteFieldPolynomial)> {
    let mut g_i = g.clone();
    let mut h_i = h.clone();
    let mut current_p = p.clone();
    for _ in 0..k.ilog2() + 1 {
        let field = PrimeField::new(current_p.clone());
        let f_mod_pi = poly_with_field(f, field.clone());
        let g_i_mod_pi = poly_with_field(&g_i, field.clone());
        let h_i_mod_pi = poly_with_field(&h_i, field.clone());
        let e = f_mod_pi - (g_i_mod_pi.clone() * h_i_mod_pi.clone());
        if e.coeffs.is_empty() {
            current_p = &current_p * &current_p;
            continue;
        }
        let e_prime_coeffs = e
            .coeffs
            .into_iter()
            .map(|c| PrimeFieldElement::new(c.value / &current_p, field.clone()))
            .collect();
        let e_prime = FiniteFieldPolynomial::new(e_prime_coeffs, field.clone());
        let (gcd, s, t) = poly_extended_gcd(g_i_mod_pi.clone(), h_i_mod_pi.clone());
        if gcd.degree() > 0 {
            return None;
        }
        let d_h = (s * e_prime.clone()).long_division(h_i_mod_pi.clone()).1;
        let d_g = (t * e_prime).long_division(g_i_mod_pi).1;
        g_i = g_i + poly_mul_scalar(&d_g, &current_p);
        h_i = h_i + poly_mul_scalar(&d_h, &current_p);
        current_p = &current_p * &current_p;
    }
    Some((g_i, h_i))
}

// =====================================================================================
// region: Cantor-Zassenhaus
// =====================================================================================

/// Factors a square-free polynomial over a large prime field using Cantor-Zassenhaus algorithm.
///
/// The Cantor-Zassenhaus algorithm is a probabilistic algorithm for factoring polynomials
/// over finite fields. It relies on distinct-degree factorization and equal-degree splitting.
///
/// # Arguments
/// * `f` - The square-free polynomial to factor.
///
/// # Returns
/// A `Vec<FiniteFieldPolynomial>` containing the irreducible factors.
pub fn cantor_zassenhaus(f: &FiniteFieldPolynomial) -> Vec<FiniteFieldPolynomial> {
    let ddf_factors = distinct_degree_factorization(f);
    let mut final_factors = Vec::new();

    for (poly_product, degree) in ddf_factors {
        if poly_product.degree() as usize == degree {
            // Already irreducible
            final_factors.push(poly_product);
        } else {
            // Perform Equal-Degree Splitting
            let mut split_factors = equal_degree_splitting(&poly_product, degree);
            final_factors.append(&mut split_factors);
        }
    }
    final_factors
}

/// Performs Distinct-Degree Factorization (DDF) of a polynomial over a finite field.
///
/// DDF groups irreducible factors by their degree. It uses the property that
/// `x^(p^d) - x` is the product of all monic irreducible polynomials over `GF(p)`
/// whose degree divides `d`.
///
/// # Arguments
/// * `f` - The polynomial to factor.
///
/// # Returns
/// A `Vec<(FiniteFieldPolynomial, usize)>` where each tuple contains a polynomial
/// (which is a product of irreducible factors of a certain degree) and that degree.
pub fn distinct_degree_factorization(
    f: &FiniteFieldPolynomial,
) -> Vec<(FiniteFieldPolynomial, usize)> {
    let mut factors = Vec::new();
    let p = &f.field.modulus;
    let x = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(One::one(), f.field.clone()),
            PrimeFieldElement::new(Zero::zero(), f.field.clone()),
        ],
        f.field.clone(),
    );

    let mut h = x.clone();
    let mut f_star = f.clone();
    let mut d = 1;

    while f_star.degree() >= 2 * (d as isize) {
        h = poly_pow_mod(h.clone(), p, &f_star);
        let g_d = poly_gcd_gf(f_star.clone(), h.clone() - x.clone());

        if g_d.degree() > 0 {
            factors.push((g_d.clone(), d));
            f_star = f_star.long_division(g_d).0;
        }
        d += 1;
    }

    if f_star.degree() > 0 {
        factors.push((f_star.clone(), f_star.degree() as usize));
    }
    factors
}

/// Performs Equal-Degree Splitting.
pub(crate) fn equal_degree_splitting(
    f: &FiniteFieldPolynomial,
    d: usize,
) -> Vec<FiniteFieldPolynomial> {
    if f.degree() as usize == d {
        return vec![f.clone()];
    }

    let mut factors = vec![f.clone()];
    let mut result = Vec::new();

    while let Some(current_f) = factors.pop() {
        if (current_f.degree() as usize) == d {
            result.push(current_f);
            continue;
        }

        let p = &current_f.field.modulus;
        let exp = (p.pow(d as u32) - BigInt::one()) / 2;

        // Loop until a non-trivial factor is found
        loop {
            let a = random_poly(current_f.degree() as usize - 1, current_f.field.clone());
            let b = poly_pow_mod(a, &exp, &current_f)
                - FiniteFieldPolynomial::new(
                    vec![PrimeFieldElement::new(One::one(), current_f.field.clone())],
                    current_f.field.clone(),
                );
            let g = poly_gcd_gf(current_f.clone(), b);

            if g.degree() > 0 && g.degree() < current_f.degree() {
                factors.push(g.clone());
                factors.push(current_f.long_division(g).0);
                break; // Found a split, break the random-try loop
            }
            // Otherwise, try another random polynomial
        }
    }
    result
}

/// Generates a random monic polynomial of a given degree.
pub(crate) fn random_poly(degree: usize, field: Arc<PrimeField>) -> FiniteFieldPolynomial {
    let mut coeffs = Vec::with_capacity(degree + 1);
    coeffs.push(PrimeFieldElement::new(One::one(), field.clone())); // Monic
    for _ in 0..degree {
        let random_val = BigInt::from(rand::random::<u64>()) % &field.modulus;
        coeffs.push(PrimeFieldElement::new(random_val, field.clone()));
    }
    FiniteFieldPolynomial::new(coeffs, field)
}

// =====================================================================================
// region: Helpers
// =====================================================================================

/// Computes the greatest common divisor (GCD) of two polynomials over a prime field.
pub fn poly_gcd_gf(a: FiniteFieldPolynomial, b: FiniteFieldPolynomial) -> FiniteFieldPolynomial {
    if b.coeffs.is_empty() || b.coeffs.iter().all(|c| c.value.is_zero()) {
        a
    } else {
        let (_, remainder) = a.clone().long_division(b.clone());
        poly_gcd_gf(b, remainder)
    }
}

/// Computes base^exp mod modulus for polynomials over a prime field.
pub(crate) fn poly_pow_mod(
    base: FiniteFieldPolynomial,
    exp: &BigInt,
    modulus: &FiniteFieldPolynomial,
) -> FiniteFieldPolynomial {
    let mut res = FiniteFieldPolynomial::new(
        vec![PrimeFieldElement::new(One::one(), base.field.clone())],
        base.field.clone(),
    );
    let mut b = base;
    let mut e = exp.clone();
    while e > Zero::zero() {
        if &e % 2 == One::one() {
            res = res.mul(b.clone()).long_division(modulus.clone()).1;
        }
        b = b.clone().mul(b.clone()).long_division(modulus.clone()).1;
        e >>= 1;
    }
    res
}

/// Helper to multiply a polynomial by a scalar BigInt.
pub fn poly_mul_scalar(poly: &FiniteFieldPolynomial, scalar: &BigInt) -> FiniteFieldPolynomial {
    let new_coeffs = poly
        .coeffs
        .iter()
        .map(|c| PrimeFieldElement::new(&c.value * scalar, c.field.clone()))
        .collect();
    FiniteFieldPolynomial::new(new_coeffs, poly.field.clone())
}
/// Helper to change the field of a polynomial's coefficients.
pub(crate) fn poly_with_field(
    poly: &FiniteFieldPolynomial,
    field: Arc<PrimeField>,
) -> FiniteFieldPolynomial {
    let new_coeffs = poly
        .coeffs
        .iter()
        .map(|c| PrimeFieldElement::new(c.value.clone(), field.clone()))
        .collect();
    FiniteFieldPolynomial::new(new_coeffs, field)
}
/// Polynomial Extended Euclidean Algorithm for `a(x)s(x) + b(x)t(x) = gcd(a(x), b(x))`.
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
