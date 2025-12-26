//! # Gröbner Bases
//!
//! This module provides functions for computing Gröbner bases of polynomial ideals.
//! It includes implementations for multivariate polynomial division, S-polynomial computation,
//! and Buchberger's algorithm for generating a Gröbner basis. Monomial orderings are also supported.
//!
//! ## Overview
//!
//! A Gröbner basis is a particular kind of generating set of an ideal in a polynomial ring.
//! Gröbner bases have many applications in computational algebra, including:
//! - Solving systems of polynomial equations
//! - Ideal membership testing
//! - Computing polynomial remainders
//! - Elimination theory
//!
//! ## Examples
//!
//! ### Computing a Gröbner Basis
//! ```
//! use rssn::symbolic::core::{Expr, Monomial, SparsePolynomial};
//! use rssn::symbolic::grobner::{buchberger, MonomialOrder};
//! use std::collections::BTreeMap;
//!
//! // Create polynomials: x^2 - y and xy - 1
//! let mut poly1_terms = BTreeMap::new();
//! let mut mono1 = BTreeMap::new();
//! mono1.insert("x".to_string(), 2);
//! poly1_terms.insert(Monomial(mono1), Expr::new_constant(1.0));
//! let mut mono2 = BTreeMap::new();
//! mono2.insert("y".to_string(), 1);
//! poly1_terms.insert(Monomial(mono2), Expr::new_constant(-1.0));
//!
//! let poly1 = SparsePolynomial { terms: poly1_terms };
//!
//! // Compute Gröbner basis
//! let basis = vec![poly1];
//! let grobner = buchberger(&basis, MonomialOrder::Lexicographical).unwrap();
//! ```

use crate::symbolic::core::{Expr, Monomial, SparsePolynomial};
use crate::symbolic::polynomial::{add_poly, mul_poly};
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// Defines the monomial ordering to be used in polynomial division.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
#[repr(C)]

pub enum MonomialOrder {
    Lexicographical,
    GradedLexicographical,
    GradedReverseLexicographical,
}

/// Compares two monomials based on a given ordering.

pub(crate) fn compare_monomials(
    m1: &Monomial,
    m2: &Monomial,
    order: MonomialOrder,
) -> Ordering {

    match order {
        MonomialOrder::Lexicographical => {
            m1.0.iter()
                .cmp(m2.0.iter())
        }
        _ => {
            m1.0.iter()
                .cmp(m2.0.iter())
        }
    }
}

/// Performs division of a multivariate polynomial `f` by a set of divisors `F`.
///
/// This function implements the multivariate division algorithm, which is a generalization
/// of polynomial long division. It repeatedly subtracts multiples of the divisors from `f`
/// until a remainder is obtained that cannot be further reduced.
///
/// # Arguments
/// * `f` - The dividend polynomial.
/// * `divisors` - A slice of `SparsePolynomial`s representing the divisors.
/// * `order` - The `MonomialOrder` to use for division.
///
/// # Returns
/// A tuple `(quotients, remainder)` where `quotients` is a `Vec<SparsePolynomial>`
/// (one for each divisor) and `remainder` is a `SparsePolynomial`.

pub fn poly_division_multivariate(
    f: &SparsePolynomial,
    divisors: &[SparsePolynomial],
    order: MonomialOrder,
) -> Result<
    (
        Vec<SparsePolynomial>,
        SparsePolynomial,
    ),
    String,
> {

    let mut quotients = vec![
        SparsePolynomial {
            terms: BTreeMap::new()
        };
        divisors.len()
    ];

    let mut remainder = SparsePolynomial {
        terms: BTreeMap::new(),
    };

    let mut p = f.clone();

    while !p.terms.is_empty() {

        let mut division_occurred = false;

        let lead_term_p = match p
            .terms
            .keys()
            .max_by(|a, b| compare_monomials(a, b, order))
        {
            Some(lt) => lt.clone(),
            None => continue,
        };

        for (i, divisor) in divisors
            .iter()
            .enumerate()
        {

            if divisor
                .terms
                .is_empty()
            {

                continue;
            }

            let lead_term_g = match divisor
                .terms
                .keys()
                .max_by(|a, b| compare_monomials(a, b, order))
            {
                Some(lt) => lt.clone(),
                None => unreachable!(),
            };

            if is_divisible(
                &lead_term_p,
                &lead_term_g,
            ) {

                let coeff_p = match p
                    .terms
                    .get(&lead_term_p)
                {
                    Some(c) => c,
                    None => {

                        return Err("Logic error: lead term not in polynomial terms".to_string());
                    }
                };

                let coeff_g = match divisor
                    .terms
                    .get(&lead_term_g)
                {
                    Some(c) => c,
                    None => {

                        return Err(
                            "Logic error: lead term not found in divisor terms".to_string(),
                        );
                    }
                };

                let coeff_ratio = simplify(&Expr::new_div(
                    coeff_p.clone(),
                    coeff_g.clone(),
                ));

                let mono_ratio = subtract_monomials(
                    &lead_term_p,
                    &lead_term_g,
                );

                let mut t_terms = BTreeMap::new();

                t_terms.insert(
                    mono_ratio,
                    coeff_ratio,
                );

                let t = SparsePolynomial { terms: t_terms };

                quotients[i] = add_poly(&quotients[i], &t);

                let t_g = mul_poly(&t, divisor);

                p = subtract_poly(&p, &t_g);

                division_occurred = true;

                break;
            }
        }

        if !division_occurred {

            let coeff = match p
                .terms
                .remove(&lead_term_p)
            {
                Some(c) => c,
                None => {

                    return Err("Logic error: lead term not found for removal".to_string());
                }
            };

            remainder
                .terms
                .insert(lead_term_p, coeff);
        }
    }

    Ok((quotients, remainder))
}

pub(crate) fn is_divisible(
    m1: &Monomial,
    m2: &Monomial,
) -> bool {

    m2.0.iter()
        .all(|(var, exp2)| {

            m1.0.get(var)
                .is_some_and(|exp1| exp1 >= exp2)
        })
}

pub(crate) fn subtract_monomials(
    m1: &Monomial,
    m2: &Monomial,
) -> Monomial {

    let mut result = m1.0.clone();

    for (var, exp2) in &m2.0 {

        let exp1 = result
            .entry(var.clone())
            .or_insert(0);

        *exp1 -= exp2;
    }

    Monomial(
        result
            .into_iter()
            .filter(|(_, exp)| *exp > 0)
            .collect(),
    )
}

#[must_use]

pub fn subtract_poly(
    p1: &SparsePolynomial,
    p2: &SparsePolynomial,
) -> SparsePolynomial {

    let mut result_terms = p1.terms.clone();

    for (mono, coeff) in &p2.terms {

        let entry = result_terms
            .entry(mono.clone())
            .or_insert_with(|| Expr::Constant(0.0));

        *entry = simplify(&Expr::new_sub(
            entry.clone(),
            coeff.clone(),
        ));
    }

    result_terms.retain(|_, v| !is_zero(v));

    SparsePolynomial {
        terms: result_terms,
    }
}

/// Computes the leading term (monomial, coefficient) of a polynomial.

pub(crate) fn leading_term(
    p: &SparsePolynomial,
    order: MonomialOrder,
) -> Option<(Monomial, Expr)> {

    p.terms
        .iter()
        .max_by(|(m1, _), (m2, _)| compare_monomials(m1, m2, order))
        .map(|(m, c)| (m.clone(), c.clone()))
}

/// Computes the least common multiple (LCM) of two monomials.

pub(crate) fn lcm_monomial(
    m1: &Monomial,
    m2: &Monomial,
) -> Monomial {

    let mut lcm_map = m1.0.clone();

    for (var, &exp2) in &m2.0 {

        let exp1 = lcm_map
            .entry(var.clone())
            .or_insert(0);

        *exp1 = std::cmp::max(*exp1, exp2);
    }

    Monomial(lcm_map)
}

/// Computes the S-polynomial of two polynomials.
/// S(f, g) = (lcm(LM(f), LM(g)) / LT(f)) * f - (lcm(LM(f), LM(g)) / LT(g)) * g

pub(crate) fn s_polynomial(
    p1: &SparsePolynomial,
    p2: &SparsePolynomial,
    order: MonomialOrder,
) -> Option<SparsePolynomial> {

    let (lm1, lc1) = leading_term(p1, order)?;

    let (lm2, lc2) = leading_term(p2, order)?;

    let lcm = lcm_monomial(&lm1, &lm2);

    let t1_mono = subtract_monomials(&lcm, &lm1);

    let t1_coeff = simplify(&Expr::new_div(
        Expr::Constant(1.0),
        lc1,
    ));

    let mut t1_terms = BTreeMap::new();

    t1_terms.insert(t1_mono, t1_coeff);

    let t1 = SparsePolynomial { terms: t1_terms };

    let t2_mono = subtract_monomials(&lcm, &lm2);

    let t2_coeff = simplify(&Expr::new_div(
        Expr::Constant(1.0),
        lc2,
    ));

    let mut t2_terms = BTreeMap::new();

    t2_terms.insert(t2_mono, t2_coeff);

    let t2 = SparsePolynomial { terms: t2_terms };

    let term1 = mul_poly(&t1, p1);

    let term2 = mul_poly(&t2, p2);

    Some(subtract_poly(
        &term1, &term2,
    ))
}

/// Computes a Gröbner basis for a polynomial ideal using Buchberger's algorithm.
///
/// Buchberger's algorithm takes a set of generators for a polynomial ideal and
/// produces a Gröbner basis for that ideal. A Gröbner basis has many desirable
/// properties, such as simplifying polynomial division and solving systems of
/// polynomial equations.
///
/// # Arguments
/// * `basis` - A slice of `SparsePolynomial`s representing the initial generators of the ideal.
/// * `order` - The `MonomialOrder` to use for computations.
///
/// # Returns
/// A `Vec<SparsePolynomial>` representing the Gröbner basis.

pub fn buchberger(
    basis: &[SparsePolynomial],
    order: MonomialOrder,
) -> Result<Vec<SparsePolynomial>, String> {

    if basis.is_empty() {

        return Ok(vec![]);
    }

    let mut g = basis.to_vec();

    let mut pairs: Vec<(usize, usize)> = (0..g.len())
        .flat_map(|i| (i + 1..g.len()).map(move |j| (i, j)))
        .collect();

    while let Some((i, j)) = pairs.pop() {

        if let Some(s_poly) = s_polynomial(&g[i], &g[j], order) {

            let (_, remainder) = poly_division_multivariate(&s_poly, &g, order)?;

            if !remainder
                .terms
                .is_empty()
            {

                let new_poly_idx = g.len();

                for k in 0..new_poly_idx {

                    pairs.push((k, new_poly_idx));
                }

                g.push(remainder);
            }
        }
    }

    Ok(g)
}
