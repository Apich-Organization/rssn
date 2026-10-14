//! Radical solutions of irreducible polynomials of degree six and more that
//! have a structure.
//!
//! * **Functional decomposition** `f = g ∘ h` (Kozen–Landau): for each
//!   divisor `d` of the degree, `h` of degree `d` is read off the expansion
//!   of `f^(1/m)` (`m = deg f / d`) and `f` is expanded in powers of `h` by
//!   repeated division; the decomposition is accepted when every remainder
//!   is a constant. The real roots of `g` are found recursively (they may be
//!   radicals, or `root_of` objects), and for each real root `u` the real
//!   roots of `h(x) = u` by the symbolic solver; a real root `x` of `f`
//!   forces `u = h(x)` to be real, so nothing is lost.
//! * **Self-reciprocal polynomials** `x^n f(1/x) = f` of degree `2k`: with
//!   `u = x + 1/x` they become a polynomial of degree `k` in `u`
//!   (`x^j + x^-j` is the Dickson polynomial in `u`), and `x` solves
//!   `x^2 - u x + 1 = 0`, real for `|u| >= 2`.
//!
//! When some step has no exact form the whole factor is left to the
//! `root_of` objects.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;

type Q = BigRational;

fn trim(p: &mut Vec<Q>) {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
}

fn mul(
    a: &[Q],
    b: &[Q],
) -> Vec<Q> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![Q::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i + j) {
                *slot += x * y;
            }
        }
    }
    out
}

fn power(
    p: &[Q],
    n: usize,
) -> Vec<Q> {
    let mut out = vec![Q::one()];
    for _ in 0..n {
        out = mul(&out, p);
    }
    out
}

/// Quotient and remainder of `a` by `b`.
fn divrem(
    a: &[Q],
    b: &[Q],
) -> Option<(Vec<Q>, Vec<Q>)> {
    let lead = b.last()?.clone();
    if lead.is_zero() {
        return None;
    }
    let mut rest = a.to_vec();
    let db = b.len() - 1;
    let mut quotient = vec![Q::zero(); a.len().saturating_sub(db).max(1)];
    while rest.len() > db {
        let top = rest.last()?.clone();
        let shift = rest.len() - 1 - db;
        let factor = top / &lead;
        for (i, c) in b.iter().enumerate() {
            let slot = rest.get_mut(shift + i)?;
            *slot -= &factor * c;
        }
        *quotient.get_mut(shift)? = factor;
        rest.pop();
        trim(&mut rest);
        if rest.is_empty() {
            break;
        }
    }
    trim(&mut quotient);
    Some((quotient, rest))
}

/// `f = g ∘ h` with `deg h = d`, as `(g, h)`.
fn decompose_with_degree(
    f: &[Q],
    d: usize,
) -> Option<(Vec<Q>, Vec<Q>)> {
    let n = f.len() - 1;
    if d < 2 || n % d != 0 || n / d < 2 {
        return None;
    }
    let m = n / d;
    let lead = f.last()?.clone();
    let monic: Vec<Q> = f.iter().map(|c| c / &lead).collect();
    let mut h = vec![Q::zero(); d + 1];
    *h.get_mut(d)? = Q::one();
    let m_q = Q::from_integer(BigInt::from(m));
    for k in 1..d {
        let hm = power(&h, m);
        let c = monic.get(n - k)?.clone() - hm.get(n - k).cloned().unwrap_or_else(Q::zero);
        *h.get_mut(d - k)? += c / &m_q;
    }
    // Expand f in powers of h.
    let mut g = Vec::new();
    let mut current = f.to_vec();
    while current.len() > d {
        let (quotient, remainder) = divrem(&current, &h)?;
        if remainder.len() > 1 {
            return None;
        }
        g.push(remainder.first().cloned().unwrap_or_else(Q::zero));
        current = quotient;
    }
    // The last quotient has degree < d, hence must be a constant (if the
    // degrees match, the loop leaves g_m).
    if current.len() > 1 {
        return None;
    }
    g.push(current.first().cloned().unwrap_or_else(Q::zero));
    trim(&mut g);
    (g.len() == m + 1).then_some((g, h))
}

/// `g(u)` with `f(x) / x^k = g(x + 1/x)` for a self-reciprocal `f` of
/// degree `2k`.
fn reciprocal_reduction(f: &[Q]) -> Option<Vec<Q>> {
    let n = f.len() - 1;
    if n % 2 != 0 || n < 4 || (0..=n).any(|i| f.get(i) != f.get(n - i)) {
        return None;
    }
    let k = n / 2;
    // Dickson polynomials D_j(u) = x^j + x^-j, D_0 = 1 (the middle term).
    let mut dickson: Vec<Vec<Q>> = vec![vec![Q::one()], vec![Q::zero(), Q::one()]];
    for j in 2..=k {
        let a = mul(&[Q::zero(), Q::one()], dickson.get(j - 1)?);
        let b = dickson.get(j - 2)?;
        let mut next = a;
        for (i, c) in b.iter().enumerate() {
            let value = next.get(i).cloned().unwrap_or_else(Q::zero) - c;
            if i < next.len() {
                *next.get_mut(i)? = value;
            } else {
                next.push(value);
            }
        }
        // D_2 = u^2 - 2: the recurrence needs D_0 = 2.
        if j == 2 {
            let first = next.first().cloned().unwrap_or_else(Q::zero) - Q::one();
            *next.get_mut(0)? = first;
        }
        dickson.push(next);
    }
    let mut g = vec![Q::zero(); k + 1];
    for j in 0..=k {
        let c = f.get(k + j)?;
        let piece = dickson.get(j)?;
        for (i, p) in piece.iter().enumerate() {
            *g.get_mut(i)? += c * p;
        }
    }
    trim(&mut g);
    Some(g)
}

fn value_of(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.eval(node, &Env::numeric(0.0)).filter(|v| v.is_finite())
}

/// The real roots of the integer polynomial `f` (ascending) in radicals
/// when it is decomposable or self-reciprocal.
pub(super) fn roots(
    graph: &mut Graph,
    f: &[BigInt],
) -> Option<Vec<(f64, NodeId)>> {
    let f: Vec<Q> = f.iter().cloned().map(Q::from_integer).collect();
    let n = f.len() - 1;
    let mut found: Vec<NodeId> = Vec::new();
    if let Some(g) = reciprocal_reduction(&f) {
        for u in super::rational_polynomial_roots(graph, &g)? {
            let value = value_of(graph, u)?;
            if value.abs() <= 2.0 {
                continue;
            }
            let two = graph.int(2);
            let minus_four = graph.int(-4);
            let square = graph.node(core::POW, &[u, two]);
            let inner = graph.node(core::ADD, &[square, minus_four]);
            let half = graph.num(Number::fraction(1, 2)?);
            let root = graph.node(core::POW, &[inner, half]);
            let minus_one = graph.int(-1);
            let minus_root = graph.node(core::MUL, &[minus_one, root]);
            for signed in [root, minus_root] {
                let sum = graph.node(core::ADD, &[u, signed]);
                found.push(graph.node(core::MUL, &[half, sum]));
            }
        }
    } else {
        let (g, h) = (2..=n / 2).find_map(|d| decompose_with_degree(&f, d))?;
        let x_symbol = graph.interner_mut().fresh_symbol("x");
        let x = graph.symbol_node(x_symbol);
        for u in super::rational_polynomial_roots(graph, &g)? {
            let mut terms = Vec::new();
            for (i, c) in h.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                let coefficient = graph.num(Number::rat(c.clone()));
                let exponent = graph.int(i64::try_from(i).ok()?);
                let power = if i == 0 { graph.int(1) } else { graph.node(core::POW, &[x, exponent]) };
                terms.push(graph.node(core::MUL, &[coefficient, power]));
            }
            let minus_one = graph.int(-1);
            let minus_u = graph.node(core::MUL, &[minus_one, u]);
            terms.push(minus_u);
            let equation = graph.node(core::ADD, &terms);
            found.extend(super::solve_for(graph, equation, x, 0)?);
        }
    }
    let mut out = Vec::new();
    for node in found {
        out.push((value_of(graph, node)?, node));
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-12 * (1.0 + a.0.abs()));
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: &[i64]) -> Vec<Q> {
        v.iter().map(|&c| Q::from_integer(BigInt::from(c))).collect()
    }

    #[test]
    fn decomposition_and_reduction() {
        // (x^2 + x + 1)^3 + 2 = g(h) with g = u^3 + 2, h = x^2 + x + 1... h has constant 0 by
        // convention, so g absorbs it: g(u) = (u + 1)^3 + 2, h = x^2 + x.
        let h = q(&[1, 1, 1]);
        let f: Vec<Q> = {
            let mut p = power(&h, 3);
            p[0] += Q::from_integer(BigInt::from(2));
            p
        };
        let (g, hh) = decompose_with_degree(&f, 2).unwrap_or_else(|| panic!("decomposes"));
        assert_eq!(hh, q(&[0, 1, 1]));
        assert_eq!(g, q(&[3, 3, 3, 1]));
        assert!(decompose_with_degree(&q(&[1, 1, 0, 0, 0, 0, 1]), 2).is_none());
        // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1 is self-reciprocal: u^3 + u^2 - 2u - 1.
        let g = reciprocal_reduction(&q(&[1, 1, 1, 1, 1, 1, 1])).unwrap_or_else(|| panic!("reciprocal"));
        assert_eq!(g, q(&[-1, -2, 1, 1]));
        // x^4 + 3x^2 + 1: u^2 + 1.
        assert_eq!(reciprocal_reduction(&q(&[1, 0, 3, 0, 1])).unwrap_or_else(|| panic!("reciprocal")), q(&[1, 0, 1]));
    }
}
