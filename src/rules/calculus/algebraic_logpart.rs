//! The logarithmic part of an algebraic integral.
//!
//! After Hermite reduction the denominator `C` is square-free.  For an
//! irreducible factor `w` of `C` with a root `alpha` the residue of
//! `(A + B y)/C` at the place `(alpha, y_alpha)` is
//!
//! ```text
//! r = (A(alpha) + B(alpha) y_alpha) / C'(alpha),
//! ```
//!
//! an element of the residue field `Q(alpha, y_alpha)` of the curve.  In the
//! algebra `Q[alpha, y]/(w(alpha), y^2 - D(alpha))` this is `a + b y` with
//! `a`, `b` polynomials in `alpha`; `b` is non-zero exactly when `B` is
//! non-zero, because `C` is co-prime to `B`.
//!
//! A logarithm of an *algebraic* function `g` contributes the integer
//! valuation `ord_P(g)` at every place, so the logarithmic part of an
//! elementary integral is `log(g)` for some `g` in `Q(x)(y)` when the
//! residues are integers, and a sum of such logarithms (equivalently a
//! logarithm of a rational function) when they are rational.  This module
//! handles three cases:
//!
//! * a residue `a` in `Q(alpha)` (no `y` component): the conjugates give
//!   `sum_sigma a(alpha_sigma) log(x - alpha_sigma)`, which is expressed
//!   over `Q(x)`; for a linear or quadratic residue field that is the usual
//!   logarithm / arctangent pair, and a constant residue of a factor of
//!   degree three or more is a plain logarithm;
//! * a residue with a `y` component: the full integral is searched for as
//!   `log(g)` with `g = P + Q y`.  The condition `g' = f g` is *linear* in
//!   `P` and `Q`, so it becomes an exact linear system over `Q`; a
//!   non-trivial solution is a genuine (and by differentiation explicitly
//!   checked) antiderivative.  This is what closes positive-genus cases
//!   such as `integral(3 D D' / (2 y (1 + y^3)), x) = log(1 + y^3)`,
//!   `y^2 = D = x^3 + 1`, whose log part is a logarithm of an algebraic
//!   function; no logarithm of a rational function can represent it.
//! * either way the `y`-trace of the answer is checked: with
//!   `N(g) = P^2 - Q^2 D` the norm, `d/dx log N(g) = Tr(g'/g) = 2 A / C`.
//!   The trace `2 a` of the residue (the sum of the two conjugate residues
//!   at the same `x = alpha`) is the base-field projection of the
//!   logarithmic part; it is a *necessary* condition, not the antiderivative
//!   itself.  For a residue with a `y` component the anti-trace part is
//!   represented by the non-rational factor `g`, never by `2 a` alone.
//!
//! A residue with a non-constant `y` component that is *not* the
//! logarithmic derivative of an algebraic function (for an elliptic or
//! hyperelliptic radicand that is the generic situation: the remaining
//! integral is of the third kind and has non-zero periods) is genuinely
//! non-elementary.  The caller then keeps the Hermite-reduced identity
//! rather than inventing a closed form.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

use super::algebraic_field as field;
use super::algebraic_field::Element;
use super::algebraic_field::Frac;
use super::algebraic_hermite as hermite;
use super::algebraic_term as term;

/// Largest denominator degree for which the logarithmic ansatz is solved.
const MAX_ANSATZ_DEGREE: usize = 24;
/// Largest radicand degree for which the ansatz is solved.
const MAX_RADICAND_DEGREE: usize = 12;
/// Extra degree allowed above the denominator degree in the ansatz (a
/// residue of `k` at a simple pole needs a zero of order `k`).
const ANSATZ_MARGIN: usize = 2;

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

/// Trims trailing zero coefficients.
fn trim(mut p: QPoly) -> QPoly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

/// The coefficient of `x^k`, or zero outside the polynomial.
fn coefficient(
    p: &[BigRational],
    k: isize,
) -> BigRational {
    if k < 0 {
        return BigRational::zero();
    }
    usize::try_from(k)
        .ok()
        .and_then(|k| p.get(k))
        .cloned()
        .unwrap_or_else(Zero::zero)
}

/// Reduced row echelon form of a rational matrix, returning the non-zero
/// rows and their pivot columns.
fn rref(
    mut rows: Vec<Vec<BigRational>>,
    columns: usize,
) -> (Vec<Vec<BigRational>>, Vec<usize>) {
    let mut pivots = Vec::new();
    let mut row = 0_usize;
    for column in 0..columns {
        let Some(swap) = (row..rows.len()).find(|&r| !rows[r][column].is_zero()) else {
            continue;
        };
        rows.swap(row, swap);
        let lead = rows[row][column].clone();
        for slot in &mut rows[row] {
            *slot /= &lead;
        }
        let pivot_row = rows[row].clone();
        for (other, row_slice) in rows.iter_mut().enumerate() {
            if other == row {
                continue;
            }
            let factor = row_slice[column].clone();
            if factor.is_zero() {
                continue;
            }
            for (slot, value) in row_slice.iter_mut().zip(&pivot_row) {
                *slot -= value * &factor;
            }
        }
        pivots.push(column);
        row = row.saturating_add(1);
        if row == rows.len() {
            break;
        }
    }
    rows.truncate(row);
    (rows, pivots)
}

/// A basis of the nullspace of the rational matrix `rows`.
fn nullspace(
    rows: Vec<Vec<BigRational>>,
    columns: usize,
) -> Vec<Vec<BigRational>> {
    let (reduced, pivots) = rref(rows, columns);
    let mut basis = Vec::new();
    for free in 0..columns {
        if pivots.contains(&free) {
            continue;
        }
        let mut vector = vec![BigRational::zero(); columns];
        vector[free] = BigRational::from_integer(BigInt::from(1));
        for (r, &pivot) in pivots.iter().enumerate() {
            vector[pivot] = -reduced[r][free].clone();
        }
        basis.push(vector);
    }
    basis
}

/// Solves `g'/g = (a + b y)/c` for `g = P + Q y` with rational polynomials
/// `P`, `Q`.
///
/// `g' = f g` is linear in `(P, Q)`; multiplying by `2 y c` clears the
/// denominators and the `y` in `y' = d'/(2 y)`, and comparing the
/// coefficients of `1` and `y` gives two polynomial identities:
///
/// ```text
/// c P' = a P + b d Q
/// c (2 d Q' + d' Q) = 2 d (a Q + b P)
/// ```
///
/// Every solution of these linear equations satisfies `g' = f g` exactly
/// (`2 y c` is not a zero divisor in `Q(x)(y)`), so the candidate is checked
/// with exact field arithmetic before it is returned.
fn log_derivative(
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let n = field::degree(c)?;
    if n == 0 || n > MAX_ANSATZ_DEGREE {
        return None;
    }
    if field::degree(d).is_none_or(|k| k > MAX_RADICAND_DEGREE) {
        return None;
    }
    let m = n.saturating_add(ANSATZ_MARGIN);
    let columns = 2_usize.saturating_mul(m.saturating_add(1));
    let bd = field::mul(b, d);
    let two = BigRational::from_integer(BigInt::from(2));
    let two_d = field::scale(d, &two);
    let c_two_d = field::mul(c, &two_d);
    let c_d_prime = field::mul(c, &field::derivative(d));
    let two_d_a = field::scale(&field::mul(d, a), &two);
    let two_d_b = field::scale(&field::mul(d, b), &two);
    let max_degree = n
        .saturating_add(m)
        .saturating_add(2)
        .saturating_add(field::degree(d).unwrap_or(0));
    let mut rows: Vec<Vec<BigRational>> = Vec::new();
    for k in 0..=max_degree {
        let k = isize::try_from(k).ok()?;
        let mut row_a = vec![BigRational::zero(); columns];
        let mut row_b = vec![BigRational::zero(); columns];
        for i in 0..=m {
            let is = isize::try_from(i).ok()?;
            let iq = BigRational::from_integer(BigInt::from(i));
            // c P' - a P - b d Q.
            let c_p = coefficient(c, k - is + 1) * &iq;
            let a_p = coefficient(a, k - is);
            row_a[i] += c_p - a_p;
            row_a[m.saturating_add(1).saturating_add(i)] -= coefficient(&bd, k - is);
            // c (2 d Q' + d' Q) - 2 d a Q - 2 d b P.
            let c_two_d_q = coefficient(&c_two_d, k - is + 1) * &iq;
            let c_d_q = coefficient(&c_d_prime, k - is);
            let two_d_a_q = coefficient(&two_d_a, k - is);
            row_b[m.saturating_add(1).saturating_add(i)] += c_two_d_q + c_d_q - two_d_a_q;
            row_b[i] -= coefficient(&two_d_b, k - is);
        }
        if row_a.iter().any(|v| !v.is_zero()) {
            rows.push(row_a);
        }
        if row_b.iter().any(|v| !v.is_zero()) {
            rows.push(row_b);
        }
    }
    if rows.is_empty() {
        return None;
    }
    for vector in nullspace(rows, columns) {
        let first = m.saturating_add(1);
        let p = trim(vector[..first].to_vec());
        let q = trim(vector[first..].to_vec());
        if p.is_empty() && q.is_empty() {
            continue;
        }
        if verifies_log_derivative(d, a, b, c, &p, &q) {
            return Some((p, q));
        }
    }
    None
}

/// Whether `g'/g = (a + b y)/c` for `g = P + Q y`, exactly.
fn verifies_log_derivative(
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
    p: &[BigRational],
    q: &[BigRational],
) -> bool {
    let Some(f) = Element::quotient(a, b, c) else {
        return false;
    };
    let g = Element {
        a: Frac::polynomial(p),
        b: Frac::polynomial(q),
    };
    let Some(derivative) = hermite::element_derivative(&g, d) else {
        return false;
    };
    let Some(product) = g.mul(&f, d) else {
        return false;
    };
    derivative.a == product.a && derivative.b == product.b
}

/// The `y`-trace of the found logarithm: the norm `N = P^2 - Q^2 d` must
/// satisfy `N'/N = 2 a / c`, the trace of the integrand.  This is the
/// base-field projection of the logarithmic part; checking it rejects an
/// ansatz solution that only matches the anti-trace component.
fn trace_matches(
    d: &[BigRational],
    a: &[BigRational],
    c: &[BigRational],
    p: &[BigRational],
    q: &[BigRational],
) -> bool {
    let q_squared_d = field::mul(&field::mul(q, q), d);
    let norm = field::sub(&field::mul(p, p), &q_squared_d);
    let Some(log_derivative) = Frac::new(&field::derivative(&norm), &norm) else {
        return true;
    };
    let two = BigRational::from_integer(BigInt::from(2));
    let Some(trace) = Frac::new(&field::scale(a, &two), c) else {
        return false;
    };
    log_derivative == trace
}

/// Tries to write `∫ (a + b y)/c dx` as `∫ polynomial + log(g)`.
fn logarithmic_antiderivative(
    graph: &mut Graph,
    x: NodeId,
    d: &[BigRational],
    a: &[BigRational],
    b: &[BigRational],
    c: &[BigRational],
) -> Option<NodeId> {
    let (qa, ra) = field::divrem(a, c)?;
    let (qb, rb) = field::divrem(b, c)?;
    if !qb.is_empty() {
        // ∫ q_b y dx is a hyperelliptic integral of the second kind; it has
        // no logarithmic primitive.
        return None;
    }
    let (p, q) = log_derivative(d, &ra, &rb, c)?;
    if !trace_matches(d, &ra, c, &p, &q) {
        return None;
    }
    let mut terms: Vec<NodeId> = Vec::new();
    if !qa.is_empty() {
        terms.push(term::poly_node(graph, x, &integrate_poly(&qa)));
    }
    let p_node = term::poly_node(graph, x, &p);
    let q_node = term::poly_node(graph, x, &q);
    let radicand_node = term::poly_node(graph, x, d);
    let sqrt_node = term::sqrt_node(graph, radicand_node);
    let qy = graph.node(crate::graph::op::core::MUL, &[q_node, sqrt_node]);
    let body = graph.node(crate::graph::op::core::ADD, &[p_node, qy]);
    terms.push(term::log_node(graph, body)?);
    Some(term::add_all(graph, &terms))
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
    if let Some(answer) = logarithmic_antiderivative(graph, x, d, a, b, c) {
        return Some(answer);
    }
    residue_logpart(graph, x, d, a, b, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one() -> BigRational {
        BigRational::from_integer(BigInt::from(1))
    }

    /// `g = 1 + D y`, so `g'/g = 3 D'(y - D^2)/(2 (1 - D^3))`.
    fn logarithmic_example(d: &[BigRational]) -> (QPoly, QPoly, QPoly) {
        let three = BigRational::from_integer(BigInt::from(3));
        let two = BigRational::from_integer(BigInt::from(2));
        let d_prime = field::derivative(d);
        let d_squared = field::mul(d, d);
        let minus_three = -three.clone();
        let a = field::scale(&field::mul(&d_squared, &d_prime), &minus_three);
        let b = field::scale(&d_prime, &three);
        let d_cubed = field::mul(&d_squared, d);
        let one_poly = vec![one()];
        let c = field::scale(&field::sub(&one_poly, &d_cubed), &two);
        (a, b, c)
    }

    #[test]
    fn logarithmic_antiderivative_closes_positive_genus() {
        for d in [
            vec![one(), BigRational::zero(), BigRational::zero(), one()],
            vec![BigRational::zero(), -one(), BigRational::zero(), one()],
            vec![
                one(),
                BigRational::zero(),
                BigRational::zero(),
                BigRational::zero(),
                one(),
            ],
        ] {
            let (a, b, c) = logarithmic_example(&d);
            let (p, q) = log_derivative(&d, &a, &b, &c).expect("no logarithmic antiderivative");
            assert!(verifies_log_derivative(&d, &a, &b, &c, &p, &q));
            assert!(trace_matches(&d, &a, &c, &p, &q));
            // The solution is the expected `g = 1 + D y`, up to a scalar.
            let scale = p.first().cloned().unwrap_or_else(Zero::zero);
            assert!(!scale.is_zero());
            assert_eq!(field::scale(&d, &scale), q, "unexpected y coefficient");
            assert_eq!(field::scale(&vec![one()], &scale), p);
        }
    }

    #[test]
    fn non_logarithmic_residue_is_rejected() {
        // `y/(x - 2)`, `y^2 = x^3 + 1`: the residues are `+-3`, but the
        // remaining differential of the second kind is not exact, so the
        // integral is not a logarithm of an algebraic function.
        let d = vec![one(), BigRational::zero(), BigRational::zero(), one()];
        let c = vec![BigRational::from_integer(BigInt::from(-2)), one()];
        assert!(log_derivative(&d, &[], &[one()], &c).is_none());
        // The trace of a pure `y` integrand is zero, so the norm check is
        // not vacuous.
        let (a, b, c) = logarithmic_example(&d);
        let (p, q) = log_derivative(&d, &a, &b, &c).expect("no logarithmic antiderivative");
        assert!(!trace_matches(&d, &[], &c, &p, &q));
    }
}

/// The residue-by-residue logarithmic part for residues without a `y`
/// component (and the polynomial part of the integrand).
fn residue_logpart(
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
            // A simple pole of `C` at a root of `D` is a ramified place:
            // `y` vanishes there and the residue needs the confluent
            // (square-root) case, which is left to the caller's identity.
            return None;
        }
        let cp_mod = field::rem(&cp, &fq)?;
        let inv = field::inverse_mod(&cp_mod, &fq)?;
        let residue_a = field::rem(&field::mul(&ra, &inv), &fq)?;
        let residue_b = field::rem(&field::mul(&rb, &inv), &fq)?;
        let degree = field::degree(&fq);
        let f_node = term::poly_node(graph, x, &fq);
        if !residue_b.is_empty() {
            // A genuine `y` component in the residue.  Its trace is the
            // base-field projection `2 a`; the anti-trace component is not
            // a logarithm of a rational function, so it is left to the
            // logarithmic antiderivative above (and otherwise reported as
            // non-elementary by the caller).
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
