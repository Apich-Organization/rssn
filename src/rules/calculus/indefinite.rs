//! Indefinite sums and products.
//!
//! | operator | value |
//! |---|---|
//! | `indefinite_sum(f, k)` | `T(k)` with `T(k + 1) - T(k) = f(k)` |
//! | `indefinite_product(f, k)` | `P(k)` with `P(k + 1) / P(k) = f(k)` |
//!
//! Indefinite sums are linear and are found term by term: Gosper's
//! algorithm for hypergeometric terms (polynomials, `c^k`, factorials,
//! binomials, their products); partial fractions with digamma and
//! polygamma functions for the rest of the rational functions (complex
//! roots of quadratic factors included); `sin(a k + b)` and
//! `cos(a k + b)`; `ln k` (`lgamma`); `k^s` for symbolic `s` (Hurwitz
//! zeta).
//!
//! Indefinite products are multiplicative: a constant `c` gives `c^k`;
//! `exp(g)` and `c^g` give `exp(Σ g)`; powers with constant exponents
//! carry over; a rational function is split into its leading coefficient
//! and the roots of its irreducible linear and quadratic factors over the
//! rationals, `Π (k - r) ↦ Γ(k - r)`.
//!
//! Every answer is checked numerically against its defining relation, and
//! the definite `sum(f, k, a, b)` and `product(f, k, a, b)` with symbolic
//! bounds fall back on `T(b + 1) - T(a)` and `P(b + 1) / P(a)`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::poly::best;
use crate::rules::poly::univariate;

fn depends(
    graph: &Graph,
    e: NodeId,
    k: NodeId,
) -> bool {
    graph.symbol_of(k).is_some_and(|s| graph.depends_on(graph.find(e), s))
}

fn call(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

fn rat(
    graph: &mut Graph,
    r: BigRational,
) -> NodeId {
    graph.num(Number::rat(r))
}

fn shifted(
    graph: &mut Graph,
    k: NodeId,
    by: NodeId,
) -> NodeId {
    graph.node(core::ADD, &[k, by])
}

/// Samples the relation `lhs(k) = rhs(k)` at a few points (integers and a
/// non-integer, with every other symbol bound to a fixed value).
fn holds(
    graph: &Graph,
    difference: NodeId,
    k: NodeId,
) -> bool {
    let Some(ks) = graph.symbol_of(k) else {
        return false;
    };
    let symbols = graph.free_symbols(graph.find(difference)).to_vec();
    let mut checked = 0;
    for point in [3.0, 4.0, 7.0, 5.37] {
        let mut env = Env::numeric(0.0);
        for (j, &s) in symbols.iter().enumerate() {
            let value = if s == ks { point } else { 0.43 + 0.29 * f64::from(u32::try_from(j % 7).unwrap_or(0)) };
            env.bind(s, value);
        }
        if let Some(v) = graph.eval(difference, &env) {
            if !v.is_finite() {
                continue;
            }
            if v.abs() > 1e-7 {
                return false;
            }
            checked += 1;
        }
    }
    checked > 0
}

/// Gosper's antidifference, first for any integer `k`, then for `k >= 0`
/// (factorials and binomials need the latter).
fn gosper(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    super::gosper::antidifference(cx, f, k, false).or_else(|| super::gosper::antidifference(cx, f, k, true))
}

/// `T(k)` with `T(k + 1) - T(k) = f(k)`, checked.
pub fn indefinite_sum(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let f = cx.simplify(f);
    let t = sum_terms(cx, f, k)?;
    let t = cx.simplify(t);
    let one = cx.graph.int(1);
    let next = shifted(cx.graph, k, one);
    let ahead = cx.graph.substitute(t, k, next);
    let minus = cx.graph.int(-1);
    let back = cx.graph.node(core::MUL, &[minus, t]);
    let minus_f = cx.graph.node(core::MUL, &[minus, f]);
    let difference = cx.graph.node(core::ADD, &[ahead, back, minus_f]);
    holds(cx.graph, difference, k).then_some(t)
}

fn sum_terms(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    if !depends(cx.graph, f, k) {
        return Some(cx.graph.node(core::MUL, &[f, k]));
    }
    // Gosper first, on the whole term (it handles sums of similar terms).
    if let Some(t) = gosper(cx, f, k) {
        return Some(t);
    }
    let expanded = crate::rules::poly::expand_form(cx.graph, f).unwrap_or(f);
    if cx.graph.op(expanded) == core::ADD {
        let terms = cx.graph.children(expanded).to_vec();
        if terms.len() > 1 {
            let mut parts = Vec::with_capacity(terms.len());
            for term in terms {
                parts.push(single_term(cx, term, k)?);
            }
            return Some(cx.graph.node(core::ADD, &parts));
        }
    }
    single_term(cx, expanded, k)
}

fn single_term(
    cx: &mut Cx<'_>,
    term: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    if !depends(cx.graph, term, k) {
        return Some(cx.graph.node(core::MUL, &[term, k]));
    }
    if let Some(t) = gosper(cx, term, k) {
        return Some(t);
    }
    // Constant factors out.
    if cx.graph.op(term) == core::MUL {
        let factors = cx.graph.children(term).to_vec();
        let (constant, varying): (Vec<NodeId>, Vec<NodeId>) = factors.iter().partition(|&&n| !depends(cx.graph, n, k));
        if !constant.is_empty() {
            let rest = if varying.len() == 1 { varying[0] } else { cx.graph.node(core::MUL, &varying) };
            let inner = single_term(cx, rest, k)?;
            let mut all = constant;
            all.push(inner);
            return Some(cx.graph.node(core::MUL, &all));
        }
    }
    if let Some(t) = rational(cx, term, k) {
        return Some(t);
    }
    special(cx, term, k)
}

/// `Σ c / (k + α)^m` by digamma/polygamma, `Σ` polynomial by Gosper.
fn rational(
    cx: &mut Cx<'_>,
    term: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let (numer, denom) = crate::rules::poly::rational_function_in(cx.graph, term, k)?;
    if denom.len() < 2 {
        return None;
    }
    let parts = crate::rules::poly::apart::apart(&numer, &denom)?;
    let mut pieces = Vec::new();
    if !parts.quotient.is_empty() {
        let q = poly_term(cx.graph, &parts.quotient, k);
        pieces.push(gosper(cx, q, k)?);
    }
    for piece in &parts.pieces {
        match (piece.factor.as_slice(), piece.numerator.as_slice()) {
            | ([p0, p1], [c]) => {
                // c / (p1 k + p0)^m = c p1^-m / (k + α)^m
                let alpha = p0 / p1;
                let mut scale = c.clone();
                for _ in 0..piece.power {
                    scale /= p1;
                }
                let alpha = rat(cx.graph, alpha);
                let z = shifted(cx.graph, k, alpha);
                let value = polygamma_tail(cx.graph, z, piece.power)?;
                let scale = rat(cx.graph, scale);
                pieces.push(cx.graph.node(core::MUL, &[scale, value]));
            },
            | ([c0, b0, a0], numerator) if piece.power == 1 => {
                // (B k + C)/(a k² + b k + c) over the complex roots r1, r2:
                // A1/(k - r1) + A2/(k - r2), A_i = (B r_i + C)/(a (r_i - r_j)).
                let big_c = numerator.first().cloned().unwrap_or_else(BigRational::zero);
                let big_b = numerator.get(1).cloned().unwrap_or_else(BigRational::zero);
                let g = &mut *cx.graph;
                let disc = b0 * b0 - BigRational::from_integer(BigInt::from(4)) * a0 * c0;
                let disc = rat(g, disc);
                let half = g.num(Number::fraction(1, 2)?);
                let root = g.node(core::POW, &[disc, half]);
                let two_a = rat(g, BigRational::from_integer(BigInt::from(2)) * a0);
                let minus_b = rat(g, -b0.clone());
                let inv_two_a = {
                    let m = g.int(-1);
                    g.node(core::POW, &[two_a, m])
                };
                let minus = g.int(-1);
                let neg_root = g.node(core::MUL, &[minus, root]);
                let r1_num = g.node(core::ADD, &[minus_b, root]);
                let r2_num = g.node(core::ADD, &[minus_b, neg_root]);
                let r1 = g.node(core::MUL, &[r1_num, inv_two_a]);
                let r2 = g.node(core::MUL, &[r2_num, inv_two_a]);
                let (bb, cc, aa) = (rat(g, big_b), rat(g, big_c), rat(g, a0.clone()));
                for (ri, rj) in [(r1, r2), (r2, r1)] {
                    let br = g.node(core::MUL, &[bb, ri]);
                    let numerator = g.node(core::ADD, &[br, cc]);
                    let neg_rj = g.node(core::MUL, &[minus, rj]);
                    let gap = g.node(core::ADD, &[ri, neg_rj]);
                    let denominator = g.node(core::MUL, &[aa, gap]);
                    let inv = g.node(core::POW, &[denominator, minus]);
                    let coefficient = g.node(core::MUL, &[numerator, inv]);
                    let neg_ri = g.node(core::MUL, &[minus, ri]);
                    let z = g.node(core::ADD, &[k, neg_ri]);
                    let digamma = call(g, "digamma", &[z])?;
                    pieces.push(g.node(core::MUL, &[coefficient, digamma]));
                }
            },
            | _ => return None,
        }
    }
    Some(cx.graph.node(core::ADD, &pieces))
}

/// `Σ 1/(z)^m` in `k`: `ψ(z)` for `m = 1`, `(-1)^(m-1) ψ^(m-1)(z)/(m-1)!`.
fn polygamma_tail(
    graph: &mut Graph,
    z: NodeId,
    m: u32,
) -> Option<NodeId> {
    if m == 1 {
        return call(graph, "digamma", &[z]);
    }
    let order = graph.int(i64::from(m) - 1);
    let value = call(graph, "polygamma", &[order, z])?;
    let factorial: BigInt = (1..m).fold(BigInt::one(), |a, j| a * BigInt::from(j));
    let sign = if m.is_multiple_of(2) { -BigInt::one() } else { BigInt::one() };
    let c = rat(graph, BigRational::new(sign, factorial));
    Some(graph.node(core::MUL, &[c, value]))
}

fn poly_term(
    graph: &mut Graph,
    coefficients: &[BigRational],
    k: NodeId,
) -> NodeId {
    let mut terms = Vec::new();
    for (p, c) in coefficients.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let c = rat(graph, c.clone());
        let e = graph.int(i64::try_from(p).unwrap_or(0));
        let power = graph.node(core::POW, &[k, e]);
        terms.push(graph.node(core::MUL, &[c, power]));
    }
    if terms.is_empty() { graph.int(0) } else { graph.node(core::ADD, &terms) }
}

/// `(a, b)` with `u = a k + b`, `a` free of `k` and non-zero.
fn linear(
    cx: &mut Cx<'_>,
    u: NodeId,
    k: NodeId,
) -> Option<(NodeId, NodeId)> {
    let d = super::derivative(cx.graph, u, k)?;
    let a = cx.simplify(d);
    if depends(cx.graph, a, k) || cx.graph.number_of(a).is_some_and(Number::is_zero) {
        return None;
    }
    let zero = cx.graph.int(0);
    let b = cx.graph.substitute(u, k, zero);
    let b = cx.simplify(b);
    Some((a, b))
}

/// `sin(a k + b)`, `cos(a k + b)`, `ln k`, `k^s`.
fn special(
    cx: &mut Cx<'_>,
    term: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let (sin, cos, ln) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?, cx.graph.ops().lookup("ln")?);
    let op = cx.graph.op(term);
    let children = cx.graph.children(term).to_vec();
    match children.as_slice() {
        | &[u] if op == sin || op == cos => {
            // Σ sin(a k + b) = -cos(a k + b - a/2) / (2 sin(a/2)),
            // Σ cos(a k + b) =  sin(a k + b - a/2) / (2 sin(a/2)).
            let (a, _) = linear(cx, u, k)?;
            let g = &mut *cx.graph;
            let half = g.num(Number::fraction(-1, 2)?);
            let shift = g.node(core::MUL, &[half, a]);
            let argument = g.node(core::ADD, &[u, shift]);
            let minus_half = g.num(Number::fraction(1, 2)?);
            let half_a = g.node(core::MUL, &[minus_half, a]);
            let s = call(g, "sin", &[half_a])?;
            let two = g.int(2);
            let minus = g.int(-1);
            let denominator = g.node(core::MUL, &[two, s]);
            let inv = g.node(core::POW, &[denominator, minus]);
            if op == sin {
                let c = call(g, "cos", &[argument])?;
                Some(g.node(core::MUL, &[minus, c, inv]))
            } else {
                let s2 = call(g, "sin", &[argument])?;
                Some(g.node(core::MUL, &[s2, inv]))
            }
        },
        | &[u] if op == ln && u == k => call(cx.graph, "lgamma", &[k]),
        | &[base, e] if op == core::POW && base == k && !depends(cx.graph, e, k) => {
            // Σ k^s = -ζ(-s, k)
            let minus = cx.graph.int(-1);
            let s = cx.graph.node(core::MUL, &[minus, e]);
            let zeta = call(cx.graph, "hurwitz_zeta", &[s, k])?;
            Some(cx.graph.node(core::MUL, &[minus, zeta]))
        },
        | _ => None,
    }
}

/// `P(k)` with `P(k + 1) / P(k) = f(k)`, checked.
pub fn indefinite_product(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let f = cx.simplify(f);
    let p = product_terms(cx, f, k, 0)?;
    let p = cx.simplify(p);
    let one = cx.graph.int(1);
    let next = shifted(cx.graph, k, one);
    let ahead = cx.graph.substitute(p, k, next);
    let minus = cx.graph.int(-1);
    let inv = cx.graph.node(core::POW, &[p, minus]);
    let ratio = cx.graph.node(core::MUL, &[ahead, inv]);
    let minus_f = cx.graph.node(core::MUL, &[minus, f]);
    let difference = cx.graph.node(core::ADD, &[ratio, minus_f]);
    holds(cx.graph, difference, k).then_some(p)
}

fn product_terms(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    depth: usize,
) -> Option<NodeId> {
    if depth > 8 {
        return None;
    }
    if !depends(cx.graph, f, k) {
        return Some(cx.graph.node(core::POW, &[f, k]));
    }
    if f == k {
        return call(cx.graph, "gamma", &[k]);
    }
    let op = cx.graph.op(f);
    let children = cx.graph.children(f).to_vec();
    let exp = cx.graph.ops().lookup("exp");
    if op == core::MUL {
        // A rational function as a whole (so that common roots combine),
        // else factor by factor.
        if let Some(p) = rational_product(cx, f, k) {
            return Some(p);
        }
        let mut parts = Vec::with_capacity(children.len());
        for c in children {
            parts.push(product_terms(cx, c, k, depth + 1)?);
        }
        return Some(cx.graph.node(core::MUL, &parts));
    }
    if let (true, &[base, e]) = (op == core::POW, children.as_slice()) {
        if !depends(cx.graph, e, k) {
            let inner = product_terms(cx, base, k, depth + 1)?;
            return Some(cx.graph.node(core::POW, &[inner, e]));
        }
        if !depends(cx.graph, base, k) {
            // c^g(k) = exp(ln c · g)
            let s = indefinite_sum(cx, e, k)?;
            return Some(cx.graph.node(core::POW, &[base, s]));
        }
        return None;
    }
    if Some(op) == exp {
        let &[g] = children.as_slice() else {
            return None;
        };
        let s = indefinite_sum(cx, g, k)?;
        return call(cx.graph, "exp", &[s]);
    }
    rational_product(cx, f, k)
}

/// `Π f` for a rational function `f` of `k` over `Q`: `lc^k Π Γ(k - r)^m`
/// over the roots of the irreducible linear and quadratic factors.
fn rational_product(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let (numer, denom) = crate::rules::poly::rational_function_in(cx.graph, f, k)?;
    let mut factors = Vec::new();
    let mut lead = BigRational::one();
    for (poly, sign) in [(&numer, 1_i64), (&denom, -1_i64)] {
        if poly.is_empty() {
            return None;
        }
        let (unit, parts) = univariate::factor(poly);
        lead = if sign > 0 { lead * unit } else { lead / unit };
        for (z, m) in parts {
            let coefficients: Vec<BigRational> = z.iter().cloned().map(BigRational::from_integer).collect();
            let leading = coefficients.last()?.clone();
            let multiplicity = i64::from(m) * sign;
            for _ in 0..m {
                lead = if sign > 0 { lead * &leading } else { lead / &leading };
            }
            factors.push((coefficients, multiplicity));
        }
    }
    let g = &mut *cx.graph;
    let lead_node = rat(g, lead);
    let mut pieces = vec![g.node(core::POW, &[lead_node, k])];
    for (c, m) in factors {
        let roots: Vec<NodeId> = match c.as_slice() {
            | [c0, c1] => vec![rat(g, -(c0 / c1))],
            | [c0, b0, a0] => {
                let disc = b0 * b0 - BigRational::from_integer(BigInt::from(4)) * a0 * c0;
                let disc = rat(g, disc);
                let half = g.num(Number::fraction(1, 2)?);
                let root = g.node(core::POW, &[disc, half]);
                let minus = g.int(-1);
                let neg_root = g.node(core::MUL, &[minus, root]);
                let minus_b = rat(g, -b0.clone());
                let two_a = rat(g, BigRational::from_integer(BigInt::from(2)) * a0);
                let inv = g.node(core::POW, &[two_a, minus]);
                let plus = g.node(core::ADD, &[minus_b, root]);
                let less = g.node(core::ADD, &[minus_b, neg_root]);
                vec![g.node(core::MUL, &[plus, inv]), g.node(core::MUL, &[less, inv])]
            },
            | _ => return None,
        };
        for r in roots {
            let minus = g.int(-1);
            let neg_r = g.node(core::MUL, &[minus, r]);
            let z = g.node(core::ADD, &[k, neg_r]);
            let gamma = call(g, "gamma", &[z])?;
            let e = g.int(m);
            pieces.push(g.node(core::POW, &[gamma, e]));
        }
    }
    Some(g.node(core::MUL, &pieces))
}

/// `Σ_(k=a)^b f = T(b + 1) - T(a)` from the indefinite sum.
pub fn definite_sum(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let t = indefinite_sum(cx, f, k)?;
    let g = &mut *cx.graph;
    let one = g.int(1);
    let after = g.node(core::ADD, &[upper, one]);
    let end = g.substitute(t, k, after);
    let start = g.substitute(t, k, lower);
    let minus = g.int(-1);
    let neg = g.node(core::MUL, &[minus, start]);
    let value = g.node(core::ADD, &[end, neg]);
    Some(cx.simplify(value))
}

/// `Π_(k=a)^b f = P(b + 1) / P(a)` from the indefinite product.
pub fn definite_product(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let p = indefinite_product(cx, f, k)?;
    let g = &mut *cx.graph;
    let one = g.int(1);
    let after = g.node(core::ADD, &[upper, one]);
    let end = g.substitute(p, k, after);
    let start = g.substitute(p, k, lower);
    let minus = g.int(-1);
    let inv = g.node(core::POW, &[start, minus]);
    let value = g.node(core::MUL, &[end, inv]);
    Some(cx.simplify(value))
}

/// Kernel for `indefinite_sum` and `indefinite_product`.
pub(super) struct Indefinite {
    pub(super) sum: crate::graph::OpId,
    pub(super) product: crate::graph::OpId,
}

impl crate::graph::Kernel for Indefinite {
    fn ops(&self) -> Vec<crate::graph::OpId> {
        vec![self.sum, self.product]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> crate::graph::Outcome {
        let &[f, k] = cx.graph.children(node) else {
            return crate::graph::Outcome::Pass;
        };
        let Some(f) = best(cx.graph, f) else {
            return crate::graph::Outcome::Pass;
        };
        let found = if cx.graph.op(node) == self.sum { indefinite_sum(cx, f, k) } else { indefinite_product(cx, f, k) };
        found.map_or(crate::graph::Outcome::Pass, crate::graph::Outcome::Equal)
    }
}
