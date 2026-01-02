#![allow(clippy::match_same_arms)]

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
//! 
//! use std::collections::BTreeMap;
//!
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::core::Monomial;
//! use rssn::symbolic::core::SparsePolynomial;
//! use rssn::symbolic::grobner::buchberger;
//! use rssn::symbolic::grobner::MonomialOrder;
//!
//! // Create polynomials: x^2 - y and xy - 1
//! let mut poly1_terms = BTreeMap::new();
//!
//! let mut mono1 = BTreeMap::new();
//!
//! mono1.insert("x".to_string(), 2);
//!
//! poly1_terms.insert(
//!     Monomial(mono1),
//!     Expr::new_constant(1.0),
//! );
//!
//! let mut mono2 = BTreeMap::new();
//!
//! mono2.insert("y".to_string(), 1);
//!
//! poly1_terms.insert(
//!     Monomial(mono2),
//!     Expr::new_constant(-1.0),
//! );
//!
//! let poly1 = SparsePolynomial {
//!     terms : poly1_terms,
//! };
//!
//! // Compute Gröbner basis
//! let basis = vec![poly1];
//!
//! let grobner = buchberger(
//!     &basis,
//!     MonomialOrder::Lexicographical,
//! )
//! .unwrap();
//! ```

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;

use crate::symbolic::core::Expr;
use crate::symbolic::core::Monomial;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::polynomial::add_poly;
use crate::symbolic::polynomial::mul_poly;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;

/// Defines the monomial ordering to be used in polynomial division.
#[derive(
    Debug,
    Clone,
    Copy,
    serde::Serialize,
    serde::Deserialize,
)]
#[repr(C)]

pub enum MonomialOrder {
    /// Dictionary order: compares exponents of variables in a fixed sequence.
    Lexicographical,
    /// Compares total degree first, then uses lexicographical order to break ties.
    GradedLexicographical,
    /// Compares total degree first, then uses reverse lexicographical order to break ties.
    GradedReverseLexicographical,
}

/// Compares two monomials based on a given ordering.

pub(crate) fn compare_monomials(
    m1: &Monomial,
    m2: &Monomial,
    order: MonomialOrder,
) -> Ordering {

    match order {
        | MonomialOrder::Lexicographical => {
            compare_lex(m1, m2)
        },
        | MonomialOrder::GradedLexicographical => {
            let deg1: u32 = m1.0.values().sum();
            let deg2: u32 = m2.0.values().sum();
            match deg1.cmp(&deg2) {
                Ordering::Equal => compare_lex(m1, m2),
                other => other,
            }
        }
        | MonomialOrder::GradedReverseLexicographical => {
            let deg1: u32 = m1.0.values().sum();
            let deg2: u32 = m2.0.values().sum();
            match deg1.cmp(&deg2) {
                Ordering::Equal => compare_revlex(m1, m2),
                other => other,
            }
        }
    }
}

pub(crate) fn compare_lex(
    m1: &Monomial,
    m2: &Monomial,
) -> Ordering {

    let all_vars: BTreeSet<&String> =
        m1.0.keys()
            .chain(m2.0.keys())
            .collect();

    for var in all_vars {

        let e1 =
            m1.0.get(var)
                .unwrap_or(&0);

        let e2 =
            m2.0.get(var)
                .unwrap_or(&0);

        match e1.cmp(e2) {
            | Ordering::Equal => {
                continue;
            },
            | other => return other,
        }
    }

    Ordering::Equal
}

pub(crate) fn compare_revlex(
    m1: &Monomial,
    m2: &Monomial,
) -> Ordering {

    let mut all_vars: Vec<&String> =
        m1.0.keys()
            .chain(m2.0.keys())
            .collect();

    all_vars.sort();

    all_vars.dedup();

    for var in all_vars
        .iter()
        .rev()
    {

        let e1 =
            m1.0.get(*var)
                .unwrap_or(&0);

        let e2 =
            m2.0.get(*var)
                .unwrap_or(&0);

        match e1.cmp(e2) {
            Ordering::Equal => continue,
            Ordering::Greater => return Ordering::Less,
            Ordering::Less => return Ordering::Greater,
        }
    }

    Ordering::Equal
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
///
/// # Errors
///
/// This function will return an error if there is a logic error during term processing,
/// such as a leading term not being found in the polynomial or divisor terms.

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
            terms : BTreeMap::new()
        };
        divisors.len()
    ];

    let mut remainder =
        SparsePolynomial {
            terms: BTreeMap::new(),
        };

    let mut p = f.clone();

    while !p.terms.is_empty() {

        let mut division_occurred =
            false;

        let lead_term_p = match p
            .terms
            .keys()
            .max_by(|a, b| {

                compare_monomials(
                    a, b, order,
                )
            }) {
            | Some(lt) => lt.clone(),
            | None => continue,
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
                | Some(lt) => lt.clone(),
                | None => unreachable!(),
            };

            if is_divisible(
                &lead_term_p,
                &lead_term_g,
            ) {

                let coeff_p = match p
                    .terms
                    .get(&lead_term_p)
                {
                    | Some(c) => c,
                    | None => {

                        return Err("Logic error: lead term not in polynomial terms".to_string());
                    },
                };

                let coeff_g =
                    match divisor
                        .terms
                        .get(
                        &lead_term_g,
                    ) {
                        | Some(c) => c,
                        | None => {

                            return Err(
                            "Logic error: lead term not found in divisor terms".to_string(),
                        );
                        },
                    };

                let coeff_ratio =
                    simplify(
                        &Expr::new_div(
                            coeff_p
                                .clone(
                                ),
                            coeff_g
                                .clone(
                                ),
                        ),
                    );

                let mono_ratio =
                    subtract_monomials(
                        &lead_term_p,
                        &lead_term_g,
                    );

                let mut t_terms =
                    BTreeMap::new();

                t_terms.insert(
                    mono_ratio,
                    coeff_ratio,
                );

                let t =
                    SparsePolynomial {
                        terms: t_terms,
                    };

                quotients[i] = add_poly(
                    &quotients[i],
                    &t,
                );

                let t_g = mul_poly(
                    &t,
                    divisor,
                );

                p = subtract_poly(
                    &p, &t_g,
                );

                division_occurred =
                    true;

                break;
            }
        }

        if !division_occurred {

            let coeff = match p
                .terms
                .remove(&lead_term_p)
            {
                | Some(c) => c,
                | None => {

                    return Err("Logic error: lead term not found for removal".to_string());
                },
            };

            remainder
                .terms
                .insert(
                    lead_term_p,
                    coeff,
                );
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
                .is_some_and(|exp1| {

                    exp1 >= exp2
                })
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

/// Subtracts one polynomial from another.
///
/// # Arguments
/// * `p1` - The first polynomial.
/// * `p2` - The second polynomial.
///
/// # Returns
/// A `SparsePolynomial` representing `p1 - p2`.

pub fn subtract_poly(
    p1: &SparsePolynomial,
    p2: &SparsePolynomial,
) -> SparsePolynomial {

    let mut result_terms =
        p1.terms.clone();

    for (mono, coeff) in &p2.terms {

        let entry = result_terms
            .entry(mono.clone())
            .or_insert_with(|| {

                Expr::Constant(0.0)
            });

        *entry =
            simplify(&Expr::new_sub(
                entry.clone(),
                coeff.clone(),
            ));
    }

    result_terms
        .retain(|_, v| !is_zero(v));

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
        .max_by(|(m1, _), (m2, _)| {

            compare_monomials(
                m1, m2, order,
            )
        })
        .map(|(m, c)| {

            (m.clone(), c.clone())
        })
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

        *exp1 =
            std::cmp::max(*exp1, exp2);
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

    let (lm1, lc1) =
        leading_term(p1, order)?;

    let (lm2, lc2) =
        leading_term(p2, order)?;

    let lcm = lcm_monomial(&lm1, &lm2);

    let t1_mono =
        subtract_monomials(&lcm, &lm1);

    let t1_coeff =
        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            lc1,
        ));

    let mut t1_terms = BTreeMap::new();

    t1_terms.insert(t1_mono, t1_coeff);

    let t1 = SparsePolynomial {
        terms: t1_terms,
    };

    let t2_mono =
        subtract_monomials(&lcm, &lm2);

    let t2_coeff =
        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            lc2,
        ));

    let mut t2_terms = BTreeMap::new();

    t2_terms.insert(t2_mono, t2_coeff);

    let t2 = SparsePolynomial {
        terms: t2_terms,
    };

    let term1 = mul_poly(&t1, p1);

    let term2 = mul_poly(&t2, p2);

    Some(subtract_poly(
        &term1,
        &term2,
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
///
/// # Errors
///
/// This function will return an error if `poly_division_multivariate` encounters
/// an error during the reduction of S-polynomials.

pub fn buchberger(
    basis: &[SparsePolynomial],
    order: MonomialOrder,
) -> Result<Vec<SparsePolynomial>, String>
{

    if basis.is_empty() {

        return Ok(vec![]);
    }

    let mut g = basis.to_vec();

    let mut pairs: Vec<(usize, usize)> =
        (0 .. g.len())
            .flat_map(|i| {

                (i + 1 .. g.len()).map(
                    move |j| (i, j),
                )
            })
            .collect();

    while let Some((i, j)) = pairs.pop()
    {

        if let Some(s_poly) =
            s_polynomial(
                &g[i], &g[j], order,
            )
        {

            let (_, remainder) = poly_division_multivariate(&s_poly, &g, order)?;

            if !remainder
                .terms
                .is_empty()
            {

                let new_poly_idx =
                    g.len();

                for k in
                    0 .. new_poly_idx
                {

                    pairs.push((
                        k,
                        new_poly_idx,
                    ));
                }

                g.push(remainder);
            }
        }
    }

    Ok(reduced_basis(
        g, order,
    ))
}

/// Reduces a Gröbner basis to its reduced form.
///
/// A reduced Gröbner basis is a unique (for a given order) basis where:
/// 1. The leading coefficient of each polynomial is 1.
/// 2. For each polynomial p in the basis, no monomial in p is divisible by any LT(g) for g in G \ {p}.

pub fn reduced_basis(
    basis: Vec<SparsePolynomial>,
    order: MonomialOrder,
) -> Vec<SparsePolynomial> {

    if basis.is_empty() {

        return vec![];
    }

    // 1. Get a minimal basis
    let mut minimal = Vec::new();

    let mut sorted_basis = basis;

    sorted_basis.retain(|p| {
        !p.terms.is_empty()
    });

    for i in 0 .. sorted_basis.len() {

        let lt_i = match leading_term(
            &sorted_basis[i],
            order,
        ) {
            | Some((m, _)) => m,
            | None => continue,
        };

        let mut redundant = false;

        for j in 0 .. sorted_basis.len()
        {

            if i == j {

                continue;
            }

            let lt_j =
                match leading_term(
                    &sorted_basis[j],
                    order,
                ) {
                    | Some((m, _)) => m,
                    | None => continue,
                };

            if is_divisible(
                &lt_i, &lt_j,
            ) {

                if lt_i == lt_j && i > j
                {

                    redundant = true;

                    break;
                } else if lt_i != lt_j {

                    redundant = true;

                    break;
                }
            }
        }

        if !redundant {

            minimal.push(
                sorted_basis[i].clone(),
            );
        }
    }

    // 2. Reduce each polynomial by the others
    let mut reduced = Vec::new();

    for i in 0 .. minimal.len() {

        let mut others =
            minimal.clone();

        others.remove(i);

        let (_, rem) =
            poly_division_multivariate(
                &minimal[i],
                &others,
                order,
            )
            .unwrap_or((
                vec![],
                minimal[i].clone(),
            ));

        // 3. Make monic
        if let Some((_, lc)) =
            leading_term(&rem, order)
        {

            let mut monic_terms =
                BTreeMap::new();

            for (m, c) in rem.terms {

                let monic_c = simplify(
                    &Expr::new_div(
                        c,
                        lc.clone(),
                    ),
                );

                if !is_zero(&monic_c) {

                    monic_terms.insert(
                        m,
                        monic_c,
                    );
                }
            }

            reduced.push(
                SparsePolynomial {
                    terms: monic_terms,
                },
            );
        }
    }

    reduced.sort_by(|p1, p2| {

        let lt1 =
            leading_term(p1, order)
                .map(|(m, _)| m);

        let lt2 =
            leading_term(p2, order)
                .map(|(m, _)| m);

        match (lt1, lt2) {
            | (Some(m1), Some(m2)) => {
                compare_monomials(
                    &m1, &m2, order,
                )
            },
            | (Some(_), None) => {
                Ordering::Greater
            },
            | (None, Some(_)) => {
                Ordering::Less
            },
            | (None, None) => {
                Ordering::Equal
            },
        }
    });

    reduced
}
