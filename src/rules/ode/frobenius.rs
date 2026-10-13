//! Linear second-order equations with polynomial coefficients: the
//! hypergeometric families and Frobenius series.
//!
//! * **Hypergeometric matching.** `α x (1 - x) y'' + (β + γ x) y' + δ y = 0`
//!   is Gauss' equation with `c = β/α`, `a + b + 1 = -γ/α`, `a b = -δ/α`:
//!   `y = C1 2F1(a, b; c; x) + C2 x^{1-c} 2F1(a-c+1, b-c+1; 2-c; x)`.
//!   `α x y'' + (β + γ x) y' + δ y = 0` is Kummer's equation in `z = s x`
//!   with `s = -γ/α`, `b = β/α`, `a = δ/γ`:
//!   `y = C1 1F1(a; b; s x) + C2 x^{1-b} 1F1(a-b+1; 2-b; s x)`. (Bessel and
//!   Airy equations are handled by the special-function method; Whittaker
//!   equations reduce to Kummer's by `y = exp(-z/2) z^{b/2} w`.) The
//!   non-integer `c`, `b` restriction is what makes the second solution
//!   a pure power times a hypergeometric function; integer values fall to
//!   the Frobenius method below.
//! * **Frobenius.** At an ordinary point or a regular singular point `x0`
//!   (`x0 = 0` first, then the rational roots of the leading coefficient)
//!   the exponents are the roots `r1 >= r2` of the indicial polynomial
//!   `F(r) = r(r - 1) + P0 r + Q0`, and the solutions are
//!   `t^r Σ a_k t^k` with `F(r + k) a_k = -Σ_j ((k - j + r) P_j + Q_j) a_{k-j}`
//!   (`t = x - x0`, `P`, `Q` the analytic coefficients of the normal form).
//!   When `r1 - r2` is a non-negative integer the second solution is
//!   `c y1 ln t + t^{r2} Σ b_k t^k` with `b_0 = 1` (`b_0 = 0`, `c = 1` for a
//!   double root), `c` fixed at `k = r1 - r2` and the rest by the same
//!   recurrence with the logarithmic forcing; `c = 0` at ordinary points.
//!   The answer is truncated, with the order written as `big_o(t^n)`.
//!   Exponents must be rational (the discriminant of the indicial
//!   polynomial a rational square).
//! * `ode_frobenius(eq, y(x), x0, n)` returns
//!   `list(y(x) = C1 y1 + C2 y2, list(r1, r2), recurrence)` where the
//!   recurrence is the explicit finite-order relation
//!   `Σ_d c_d(k) a(k - d) = 0` of the larger exponent.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::Problem;
use super::add;
use super::explicit;
use super::mul;
use super::neg;
use super::reduce::fresh_function;
use super::reduce::linear_second_order;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::poly::univariate;

type Q = BigRational;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

/// Coefficient lists of `p2`, `p1`, `p0` (ascending powers of `x`) cleared
/// of denominators, from the coefficients of a linear equation.
fn polynomials(
    cx: &mut Cx<'_>,
    a: &[NodeId; 3],
    x: NodeId,
) -> Option<[Vec<Q>; 3]> {
    let mut fractions = Vec::new();
    for &c in a {
        let (n, d) = crate::rules::poly::rational_function_in(cx.graph, c, x)?;
        fractions.push((n, d));
    }
    let mut product_of_denominators = vec![q(1)];
    for (_, d) in &fractions {
        product_of_denominators = multiply(&product_of_denominators, d);
    }
    let mut out: [Vec<Q>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (k, (n, _)) in fractions.iter().enumerate() {
        let mut numerator = n.clone();
        for (j, (_, d)) in fractions.iter().enumerate() {
            if j != k {
                numerator = multiply(&numerator, d);
            }
        }
        trim(&mut numerator);
        *out.get_mut(k)? = numerator;
    }
    let _ = product_of_denominators;
    Some(out)
}

fn trim(p: &mut Vec<Q>) {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
}

fn multiply(
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

/// `p(t + h)`.
fn shift(
    p: &[Q],
    h: &Q,
) -> Vec<Q> {
    let mut work = p.to_vec();
    let n = work.len();
    for i in 0..n {
        for j in (i..n.saturating_sub(1)).rev() {
            let next = work.get(j + 1).cloned().unwrap_or_else(Q::zero);
            if let Some(slot) = work.get_mut(j) {
                *slot = &*slot + h * next;
            }
        }
    }
    work
}

fn node_of(
    graph: &mut Graph,
    value: &Q,
) -> NodeId {
    graph.num(Number::rat(value.clone()))
}

fn rational_sqrt(value: &Q) -> Option<Q> {
    if value.is_negative() {
        return None;
    }
    let (n, d) = (value.numer().sqrt(), value.denom().sqrt());
    (&n * &n == *value.numer() && &d * &d == *value.denom()).then(|| Q::new(n, d))
}

fn is_integer(value: &Q) -> bool {
    value.denom().is_one()
}

// ---------------------------------------------------------------------------
// Hypergeometric matching

/// Gauss and Kummer equations with polynomial coefficients.
pub(super) fn hypergeometric(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    _depth: u32,
) -> Option<NodeId> {
    let (a, forcing) = linear_second_order(cx, problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    let [p0, p1, p2] = polynomials(cx, &a, x)?;
    let coefficient = |p: &[Q], i: usize| p.get(i).cloned().unwrap_or_else(Q::zero);
    if p0.len() > 1 || p1.len() > 2 || p2.len() < 2 || !coefficient(&p2, 0).is_zero() {
        return None;
    }
    let alpha = coefficient(&p2, 1);
    let (beta, gamma, delta) = (coefficient(&p1, 0), coefficient(&p1, 1), coefficient(&p0, 0));
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let (first, second) = if p2.len() == 3 && coefficient(&p2, 2) == -alpha.clone() {
        // Gauss.
        let c = &beta / &alpha;
        if is_integer(&c) {
            return None;
        }
        let s = -&gamma / &alpha - q(1);
        let disc = rational_sqrt(&(&s * &s + q(4) * &delta / &alpha))?;
        let (pa, pb) = ((&s + &disc) / q(2), (&s - &disc) / q(2));
        let f = cx.graph.ops().lookup("hyper2f1")?;
        let (na, nb, nc) = (node_of(cx.graph, &pa), node_of(cx.graph, &pb), node_of(cx.graph, &c));
        let y1 = cx.graph.node(f, &[na, nb, nc, x]);
        let (a2, b2, c2n) = (&pa - &c + q(1), &pb - &c + q(1), q(2) - &c);
        let (na2, nb2, nc2) = (node_of(cx.graph, &a2), node_of(cx.graph, &b2), node_of(cx.graph, &c2n));
        let h = cx.graph.node(f, &[na2, nb2, nc2, x]);
        let exponent = node_of(cx.graph, &(q(1) - &c));
        let power = cx.graph.node(core::POW, &[x, exponent]);
        (y1, mul(cx.graph, &[power, h]))
    } else if p2.len() == 2 && !gamma.is_zero() {
        // Kummer in z = s x.
        let s = -&gamma / &alpha;
        let b = &beta / &alpha;
        if is_integer(&b) {
            return None;
        }
        let pa = &delta / &gamma;
        let f = cx.graph.ops().lookup("hyper1f1")?;
        let sn = node_of(cx.graph, &s);
        let z = mul(cx.graph, &[sn, x]);
        let (na, nb) = (node_of(cx.graph, &pa), node_of(cx.graph, &b));
        let y1 = cx.graph.node(f, &[na, nb, z]);
        let (a2, b2) = (&pa - &b + q(1), q(2) - &b);
        let (na2, nb2) = (node_of(cx.graph, &a2), node_of(cx.graph, &b2));
        let h = cx.graph.node(f, &[na2, nb2, z]);
        let exponent = node_of(cx.graph, &(q(1) - &b));
        let power = cx.graph.node(core::POW, &[x, exponent]);
        (y1, mul(cx.graph, &[power, h]))
    } else {
        return None;
    };
    let terms = [mul(cx.graph, &[c1, first]), mul(cx.graph, &[c2, second])];
    let solution = add(cx.graph, &terms);
    Some(explicit(cx, problem, solution))
}

// ---------------------------------------------------------------------------
// Frobenius

/// The solutions at a point: exponents, the recurrence, series.
pub(super) struct Frobenius {
    /// Exponents `r1 >= r2`.
    pub(super) roots: (Q, Q),
    /// `a_k` for `r1` (starting at `a_0 = 1`).
    pub(super) first: Vec<Q>,
    /// Second solution `c y1 ln t + t^{r2} Σ b_k t^k`.
    pub(super) log_coefficient: Q,
    /// `b_k`.
    pub(super) second: Vec<Q>,
    /// `c_d(k)` of the recurrence for `r1` as polynomials in `k`
    /// (ascending coefficients), `d = 0, 1, ...`.
    pub(super) recurrence: Vec<Vec<Q>>,
}

/// Power-series quotient `numerator / denominator` to `n` terms.
fn series_quotient(
    numerator: &[Q],
    denominator: &[Q],
    n: usize,
) -> Option<Vec<Q>> {
    let lead = denominator.first()?.clone();
    if lead.is_zero() {
        return None;
    }
    let mut out: Vec<Q> = Vec::with_capacity(n);
    for k in 0..n {
        let mut value = numerator.get(k).cloned().unwrap_or_else(Q::zero);
        for j in 1..=k {
            value -= denominator.get(j).cloned().unwrap_or_else(Q::zero) * out.get(k - j)?;
        }
        out.push(value / &lead);
    }
    Some(out)
}

/// Frobenius data at `t = 0` of `p2 y'' + p1 y' + p0 y = 0` (polynomials in
/// `t`), with `n + 1` terms.
pub(super) fn frobenius_at_zero(
    p: &[Vec<Q>; 3],
    n: usize,
) -> Option<Frobenius> {
    let [p0, p1, p2] = p;
    let order = |v: &[Q]| v.iter().position(|c| !c.is_zero());
    let m = order(p2)?;
    // Regular singular (or ordinary): ord p1 >= m - 1, ord p0 >= m - 2.
    if order(p1).is_some_and(|o| o + 1 < m) || order(p0).is_some_and(|o| o + 2 < m) {
        return None;
    }
    let s: Vec<Q> = p2.iter().skip(m).cloned().collect();
    // P = t p1 / (t^m s), Q = t^2 p0 / (t^m s): shift exponents by 1 - m, 2 - m.
    let shifted = |v: &[Q], by: i64| -> Vec<Q> {
        // coefficients of t^(by) * v
        let by = usize::try_from(by).ok();
        match by {
            | Some(b) => std::iter::repeat_n(Q::zero(), b).chain(v.iter().cloned()).collect(),
            | None => Vec::new(),
        }
    };
    let drop_low = |v: &[Q], count: usize| -> Vec<Q> { v.iter().skip(count).cloned().collect() };
    let m_i = i64::try_from(m).ok()?;
    let p_num = if m_i >= 1 { drop_low(p1, usize::try_from(m_i - 1).ok()?) } else { shifted(p1, 1 - m_i) };
    let q_num = if m_i >= 2 { drop_low(p0, usize::try_from(m_i - 2).ok()?) } else { shifted(p0, 2 - m_i) };
    let terms = n + 2;
    let pp = series_quotient(&p_num, &s, terms)?;
    let qq = series_quotient(&q_num, &s, terms)?;
    let (p_0, q_0) = (pp.first()?.clone(), qq.first()?.clone());
    // F(r) = r^2 + (P0 - 1) r + Q0.
    let b = &p_0 - q(1);
    let disc = rational_sqrt(&(&b * &b - q(4) * &q_0))?;
    let (r1, r2) = ((-&b + &disc) / q(2), (-&b - &disc) / q(2));
    let f = |r: &Q| -> Q { r * r + &b * r + &q_0 };
    let step = |r: &Q, k: usize, series: &[Q]| -> Option<Q> {
        // -Σ_{j=1..k} ((k - j + r) P_j + Q_j) a_{k-j}
        let mut sum = Q::zero();
        for j in 1..=k {
            let kj = Q::from_integer(BigInt::from(k - j));
            let weight = (&kj + r) * pp.get(j)? + qq.get(j)?;
            sum += weight * series.get(k - j)?;
        }
        Some(-sum)
    };
    let kq = |k: usize| Q::from_integer(BigInt::from(k));
    // First solution.
    let mut first = vec![q(1)];
    for k in 1..=n {
        let denominator = f(&(&r1 + kq(k)));
        if denominator.is_zero() {
            return None;
        }
        let value = step(&r1, k, &first)? / denominator;
        first.push(value);
    }
    // Second solution.
    let difference = &r1 - &r2;
    let (log_coefficient, second);
    if !is_integer(&difference) {
        let mut series = vec![q(1)];
        for k in 1..=n {
            let value = step(&r2, k, &series)? / f(&(&r2 + kq(k)));
            series.push(value);
        }
        log_coefficient = Q::zero();
        second = series;
    } else {
        let big_n = usize::try_from(difference.to_integer().to_i64()?).ok()?;
        // G_k = Σ_{m + j = k - N} a_m (P_j + δ_{j0} (2 (m + r1) - 1)).
        let g = |k: usize| -> Option<Q> {
            if k < big_n {
                return Some(Q::zero());
            }
            let rest = k - big_n;
            let mut sum = Q::zero();
            for m_index in 0..=rest {
                let j = rest - m_index;
                let mut weight = pp.get(j)?.clone();
                if j == 0 {
                    weight += q(2) * (kq(m_index) + &r1) - q(1);
                }
                sum += weight * first.get(m_index)?;
            }
            Some(sum)
        };
        let mut series = vec![if big_n == 0 { Q::zero() } else { q(1) }];
        let mut c = if big_n == 0 { q(1) } else { Q::zero() };
        for k in 1..=n {
            if big_n > 0 && k == big_n {
                // F(r1) = 0: the logarithmic coefficient absorbs the residue.
                let forcing = step(&r2, k, &series)?;
                c = -forcing / g(k)?;
                series.push(Q::zero());
                continue;
            }
            let value = (step(&r2, k, &series)? - &c * g(k)?) / f(&(&r2 + kq(k)));
            series.push(value);
        }
        log_coefficient = c;
        second = series;
    }
    // The explicit recurrence for r1: Σ_d c_d(k) a_{k-d} = 0.
    let lowest = [order(p2).unwrap_or(usize::MAX), order(p1).map_or(usize::MAX, |o| o + 1), order(p0).map_or(usize::MAX, |o| o + 2)]
        .into_iter()
        .min()?;
    let depth = p2.len().max(p1.len() + 1).max(p0.len() + 2).saturating_sub(lowest);
    let mut recurrence = Vec::new();
    let at = |v: &[Q], i: i64| -> Q { usize::try_from(i).ok().and_then(|i| v.get(i).cloned()).unwrap_or_else(Q::zero) };
    for d in 0..depth {
        let d_i = i64::try_from(d).ok()?;
        let mu = i64::try_from(lowest).ok()?;
        // (k - d + r)(k - d + r - 1) A + (k - d + r) B + C as polynomial in k.
        let shift_by = &r1 - Q::from_integer(BigInt::from(d_i));
        let (a, b_c, c_c) = (at(p2, mu + d_i), at(p1, mu + d_i - 1), at(p0, mu + d_i - 2));
        let u = &shift_by;
        let poly = vec![
            &a * u * (u - q(1)) + &b_c * u + &c_c,
            &a * (q(2) * u - q(1)) + &b_c,
            a,
        ];
        let mut poly = poly;
        trim(&mut poly);
        recurrence.push(poly);
    }
    while recurrence.last().is_some_and(Vec::is_empty) {
        recurrence.pop();
    }
    Some(Frobenius { roots: (r1, r2), first, log_coefficient, second, recurrence })
}

/// The point to expand about and the shifted polynomials.
fn expansion(
    p: &[Vec<Q>; 3],
    n: usize,
    preferred: Option<&Q>,
) -> Option<(Q, Frobenius)> {
    let mut points: Vec<Q> = Vec::new();
    if let Some(x0) = preferred {
        points.push(x0.clone());
    } else {
        points.push(Q::zero());
        let leading: Vec<Q> = p[2].clone();
        for r in univariate::rational_roots(&leading) {
            if !points.contains(&r) {
                points.push(r);
            }
        }
    }
    for x0 in points {
        let shifted = [shift(&p[0], &x0), shift(&p[1], &x0), shift(&p[2], &x0)];
        if let Some(found) = frobenius_at_zero(&shifted, n) {
            return Some((x0, found));
        }
    }
    None
}

fn series_node(
    cx: &mut Cx<'_>,
    t: NodeId,
    coefficients: &[Q],
    exponent: &Q,
    degree: usize,
) -> Option<NodeId> {
    let mut terms = Vec::new();
    for (k, c) in coefficients.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let cn = node_of(cx.graph, c);
        let power = if k == 0 {
            cx.graph.int(1)
        } else {
            let e = cx.graph.int(i64::try_from(k).ok()?);
            cx.graph.node(core::POW, &[t, e])
        };
        terms.push(mul(cx.graph, &[cn, power]));
    }
    let order_symbol = cx.graph.ops().lookup("big_o")?;
    let e = cx.graph.int(i64::try_from(degree + 1).ok()?);
    let tail = cx.graph.node(core::POW, &[t, e]);
    terms.push(cx.graph.node(order_symbol, &[tail]));
    let polynomial = add(cx.graph, &terms);
    if exponent.is_zero() {
        return Some(polynomial);
    }
    let en = node_of(cx.graph, exponent);
    let power = cx.graph.node(core::POW, &[t, en]);
    Some(mul(cx.graph, &[power, polynomial]))
}

/// The two solutions as terms in `t = x - x0`.
fn solutions(
    cx: &mut Cx<'_>,
    x: NodeId,
    x0: &Q,
    data: &Frobenius,
    n: usize,
) -> Option<(NodeId, NodeId)> {
    let t = if x0.is_zero() {
        x
    } else {
        let h = node_of(cx.graph, &-x0.clone());
        cx.graph.node(core::ADD, &[x, h])
    };
    let y1 = series_node(cx, t, &data.first, &data.roots.0, n)?;
    let mut y2 = series_node(cx, t, &data.second, &data.roots.1, n)?;
    if !data.log_coefficient.is_zero() {
        let ln = cx.graph.ops().lookup("ln")?;
        let log = cx.graph.node(ln, &[t]);
        let c = node_of(cx.graph, &data.log_coefficient);
        let logged = mul(cx.graph, &[c, y1, log]);
        y2 = add(cx.graph, &[logged, y2]);
    }
    Some((y1, y2))
}

/// `dsolve` fallback: the truncated Frobenius series (eight terms).
pub(super) fn series_solution(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
) -> Option<NodeId> {
    const TERMS: usize = 8;
    let (a, forcing) = linear_second_order(cx, problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    let p = polynomials(cx, &a, x)?;
    let (x0, data) = expansion(&p, TERMS, None)?;
    let (y1, y2) = solutions(cx, x, &x0, &data, TERMS)?;
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
    let solution = add(cx.graph, &terms);
    Some(explicit(cx, problem, solution))
}

/// `ode_frobenius(equation, y(x), x0, n)`.
pub(super) fn report(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    x0: NodeId,
    n: NodeId,
) -> Option<NodeId> {
    let n = usize::try_from(cx.graph.number_of(n)?.to_i64()?).ok().filter(|&n| n <= 60)?;
    let center = cx.graph.number_of(x0)?.to_rational()?;
    let mut problem = super::parse(cx.graph, equation, unknown)?;
    let (a, forcing) = linear_second_order(cx, &problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    let p = polynomials(cx, &a, x)?;
    let (x0, data) = expansion(&p, n, Some(&center))?;
    let (y1, y2) = solutions(cx, x, &x0, &data, n)?;
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
    let solution = add(cx.graph, &terms);
    let answer = explicit(cx, &problem, solution);
    let roots = {
        let (r1, r2) = (node_of(cx.graph, &data.roots.0), node_of(cx.graph, &data.roots.1));
        cx.graph.node(core::LIST, &[r1, r2])
    };
    // Σ_d c_d(k) a(k - d) = 0.
    let k_symbol = cx.graph.interner_mut().fresh_symbol("k");
    let k = cx.graph.symbol_node(k_symbol);
    let a_symbol = cx.graph.interner_mut().fresh_symbol("a");
    let a_fn = cx.graph.symbol_node(a_symbol);
    let mut parts = Vec::new();
    for (d, coefficients) in data.recurrence.iter().enumerate() {
        let mut terms = Vec::new();
        for (i, c) in coefficients.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let cn = node_of(cx.graph, c);
            let e = cx.graph.int(i64::try_from(i).ok()?);
            let power = cx.graph.node(core::POW, &[k, e]);
            terms.push(mul(cx.graph, &[cn, power]));
        }
        let weight = add(cx.graph, &terms);
        let index = if d == 0 {
            k
        } else {
            let back = cx.graph.int(-i64::try_from(d).ok()?);
            cx.graph.node(core::ADD, &[k, back])
        };
        let term = cx.graph.node(core::APPLY, &[a_fn, index]);
        parts.push(mul(cx.graph, &[weight, term]));
    }
    let recurrence_lhs = add(cx.graph, &parts);
    let zero = cx.graph.int(0);
    let recurrence = cx.graph.node(core::EQ, &[recurrence_lhs, zero]);
    let _ = (neg, fresh_function);
    Some(cx.graph.node(core::LIST, &[answer, roots, recurrence]))
}

// ---------------------------------------------------------------------------
// Numeric semantics of the hypergeometric functions

/// `1F1(a; b; z)` by its series.
pub(super) fn hyper1f1(args: &[f64]) -> f64 {
    let (Some(&a), Some(&b), Some(&z)) = (args.first(), args.get(1), args.get(2)) else { return f64::NAN };
    let (mut term, mut sum) = (1.0_f64, 1.0_f64);
    for n in 0..5000_u32 {
        let n = f64::from(n);
        term *= (a + n) / ((b + n) * (n + 1.0)) * z;
        sum += term;
        if !term.is_finite() {
            return f64::NAN;
        }
        if term.abs() < 1e-17 * sum.abs().max(1e-300) && n > 2.0 + z.abs() {
            break;
        }
    }
    sum
}

/// `2F1(a, b; c; z)` for `z < 1` (series, with Pfaff's transformation for
/// `z < -1/2`).
pub(super) fn hyper2f1(args: &[f64]) -> f64 {
    let (Some(&a), Some(&b), Some(&c), Some(&z)) = (args.first(), args.get(1), args.get(2), args.get(3)) else {
        return f64::NAN;
    };
    if z >= 1.0 {
        return f64::NAN;
    }
    if z < -0.5 {
        let w = z / (z - 1.0);
        return (1.0 - z).powf(-a) * hyper2f1(&[a, c - b, c, w]);
    }
    let (mut term, mut sum) = (1.0_f64, 1.0_f64);
    for n in 0..20_000_u32 {
        let n = f64::from(n);
        term *= (a + n) * (b + n) / ((c + n) * (n + 1.0)) * z;
        sum += term;
        if !term.is_finite() {
            return f64::NAN;
        }
        if term.abs() < 1e-17 * sum.abs().max(1e-300) {
            break;
        }
    }
    sum
}
