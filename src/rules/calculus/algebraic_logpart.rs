//! The logarithmic part of an algebraic integral.
//!
//! After Hermite reduction the denominator `C` is square-free.  For an
//! irreducible factor `f` of `C` and a root `alpha` of `f` the residue of
//! `(A + B y)/C` at the place `(alpha, y_alpha)` is
//!
//! ```text
//! r = (A(alpha) + B(alpha) y_alpha) / C'(alpha),
//! ```
//!
//! an element of `Q[x, y] / (f, y^2 - D)`.  A logarithm can be written
//! over `Q(x)` only when the residues of the conjugate places agree, which
//! forces `r` to lie in the base residue field.  This module handles the
//! cases that do:
//!
//! * `r` rational: `r log(f)`;
//! * `r = u + v x` in an irreducible quadratic residue field: the conjugate
//!   pair becomes `(u - v p) log(f) - 2 v m atan((2x + p)/m)` with
//!   `f = x^2 + p x + q` and `m = sqrt(-(p^2 - 4 q))`;
//! * `r` rational for a factor of degree three or more: `r log(f)`.
//!
//! A residue with a `y` component (which, since `C` is co-prime to `B`,
//! happens whenever `B` is not zero) would need a logarithm whose divisor
//! singles out one of the two places above `alpha`.  For a radicand of
//! degree three or four such a function need not exist at all — the curve
//! has positive genus — and constructing it when it does is outside the
//! scope of the Hermite step.  The caller then returns the reduced
//! integral, which is still an identity, instead of a wrong closed form.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

use super::algebraic_field as field;
use super::algebraic_term as term;

/// `∫ q dx` for a polynomial `q`.
fn integrate_poly(q: &[BigRational]) -> QPoly {
    let mut out: QPoly = std::iter::once(BigRational::zero())
        .chain(
            q.iter()
                .enumerate()
                .map(|(k, c)| c / BigRational::from_integer(BigInt::from(k.saturating_add(1)))),
        )
        .collect();
    while out.last().is_some_and(Zero::is_zero) {
        out.pop();
    }
    out
}

/// Converts a factor over `Z` to a monic polynomial over `Q`.
fn monic_factor(f: &[BigInt]) -> QPoly {
    let q: QPoly = f
        .iter()
        .map(|c| BigRational::from_integer(c.clone()))
        .collect();
    field::monic(&q)
}

/// The complete logarithmic part of `∫ (a + b y)/c dx`, including the
/// polynomial part, or `None` if some residue cannot be expressed over
/// `Q(x)`.
pub(super) fn complete(
    graph: &mut Graph,
    x: NodeId,
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> Option<NodeId> {
    let mut terms: Vec<NodeId> = Vec::new();
    // Polynomial part of the integrand: (q_a + q_b y) + (r_a + r_b y)/c.
    let (qa, ra) = field::divrem(a, c)?;
    let (qb, rb) = field::divrem(b, c)?;
    if !qb.is_empty() {
        // ∫ q_b y dx is a hyperelliptic (in general non-elementary)
        // integral; do not pretend to have found a closed form.
        return None;
    }
    if !qa.is_empty() {
        let integrated = integrate_poly(&qa);
        terms.push(term::poly_node(graph, x, &integrated));
    }
    if field::degree(c).is_none_or(|k| k == 0) {
        return Some(term::add_all(graph, &terms));
    }
    let (_, factors) = univariate::factor(c);
    if factors.iter().any(|(_, e)| *e != 1) {
        return None;
    }
    let cp = field::derivative(c);
    for (f, _) in &factors {
        let fq = monic_factor(f);
        let shared = field::gcd(&fq, d);
        if field::degree(&shared).is_none_or(|k| k > 0) {
            return None;
        }
        let cp_mod = field::rem(&cp, &fq)?;
        let inv = field::inverse_mod(&cp_mod, &fq)?;
        let residue_a = field::rem(&field::mul(&ra, &inv), &fq)?;
        let residue_b = field::rem(&field::mul(&rb, &inv), &fq)?;
        let degree = field::degree(&fq);
        let f_node = term::poly_node(graph, x, &fq);
        if !residue_b.is_empty() {
            // A `y` component in the residue: the conjugate places have
            // different residues, so no logarithm of a rational function
            // produces them.
            return None;
        }
        let residue = residue_a;
        let coefficient = residue.first().cloned().unwrap_or_else(Zero::zero);
        match degree {
            | Some(1) => {
                if coefficient.is_zero() {
                    continue;
                }
                let log = term::log_node(graph, f_node)?;
                let scale = term::rat_node(graph, &coefficient);
                terms.push(graph.node(crate::graph::op::core::MUL, &[scale, log]));
            },
            | Some(2) => {
                let v = residue.get(1).cloned().unwrap_or_else(Zero::zero);
                let p = fq.get(1).cloned().unwrap_or_else(Zero::zero);
                let q0 = fq.first().cloned().unwrap_or_else(Zero::zero);
                let four = BigRational::from_integer(BigInt::from(4));
                let p_squared = &p * &p;
                let four_q = &four * &q0;
                let discriminant = p_squared - four_q;
                if v.is_zero() {
                    if coefficient.is_zero() {
                        continue;
                    }
                    let log = term::log_node(graph, f_node)?;
                    let scale = term::rat_node(graph, &coefficient);
                    terms.push(graph.node(crate::graph::op::core::MUL, &[scale, log]));
                    continue;
                }
                if discriminant.is_zero() || discriminant.is_positive() {
                    // A repeated or split quadratic would have been split
                    // by `factor`; be conservative.
                    return None;
                }
                let minus_discriminant = -discriminant;
                let root = term::rat_node(graph, &minus_discriminant);
                let m = term::sqrt_node(graph, root);
                let log = term::log_node(graph, f_node)?;
                // (u - v p/2) log(f) - v m atan((2x + p)/m)
                let half = BigRational::from_integer(BigInt::from(1))
                    / BigRational::from_integer(BigInt::from(2));
                let log_coefficient = &coefficient - &v * &p * &half;
                if !log_coefficient.is_zero() {
                    let scale = term::rat_node(graph, &log_coefficient);
                    terms.push(graph.node(crate::graph::op::core::MUL, &[scale, log]));
                }
                let slope = BigRational::from_integer(BigInt::from(2));
                let linear = term::poly_node(graph, x, &[p, slope]);
                let arg = term::div_node(graph, linear, m);
                let arc = term::atan_node(graph, arg)?;
                let minus_v = -v;
                let scale = term::rat_node(graph, &minus_v);
                let first = graph.node(crate::graph::op::core::MUL, &[scale, m]);
                terms.push(graph.node(crate::graph::op::core::MUL, &[first, arc]));
            },
            | _ => {
                if residue.len() > 1 {
                    return None;
                }
                if coefficient.is_zero() {
                    continue;
                }
                let log = term::log_node(graph, f_node)?;
                let scale = term::rat_node(graph, &coefficient);
                terms.push(graph.node(crate::graph::op::core::MUL, &[scale, log]));
            },
        }
    }
    Some(term::add_all(graph, &terms))
}
