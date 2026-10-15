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

/// Gauss, Kummer and Whittaker equations with polynomial coefficients.
///
/// When the parameter `c` (or `b`, or `2μ`) is an integer the two power
/// solutions coincide or one is undefined: the solution that exists is
/// returned together with a second one by reduction of order
/// (`y1 ∫ exp(-∫ a1/a2) / y1²`), which is the degenerate case of the
/// formulas (a logarithm in the hypergeometric form).
pub(super) fn hypergeometric(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    _depth: u32,
) -> Option<NodeId> {
    hypergeometric_with(cx, problem, false)
}

/// The degenerate (integer-parameter) cases, tried after the other
/// methods because the second solution is a quadrature.
pub(super) fn hypergeometric_degenerate(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    _depth: u32,
) -> Option<NodeId> {
    hypergeometric_with(cx, problem, true)
}

fn hypergeometric_with(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    degenerate: bool,
) -> Option<NodeId> {
    let (a, forcing) = linear_second_order(cx, problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    let [p0, p1, p2] = polynomials(cx, &a, x)?;
    let coefficient = |p: &[Q], i: usize| p.get(i).cloned().unwrap_or_else(Q::zero);
    let power_of = |cx: &mut Cx<'_>, e: &Q| {
        let exponent = node_of(cx.graph, e);
        cx.graph.node(core::POW, &[x, exponent])
    };
    let (first, second): (NodeId, Option<NodeId>) = if p2.len() == 3 && coefficient(&p2, 0).is_zero() && coefficient(&p2, 1).is_zero() && p1.is_empty() && p0.len() <= 3 {
        // Whittaker: y'' + (-s²/4 + κ s/x + (1/4 - μ²)/x²) y = 0, z = s x.
        let alpha = coefficient(&p2, 2);
        let (big_a, big_b, big_c) = (coefficient(&p0, 2) / &alpha, coefficient(&p0, 1) / &alpha, coefficient(&p0, 0) / &alpha);
        if !big_a.is_negative() || big_a.is_zero() {
            return None;
        }
        let s = q(2) * rational_sqrt(&-big_a)?;
        let kappa = &big_b / &s;
        let mu = rational_sqrt(&(Q::new(BigInt::one(), BigInt::from(4)) - big_c))?;
        let f = cx.graph.ops().lookup("hyper1f1")?;
        let exp = cx.graph.ops().lookup("exp")?;
        let sn = node_of(cx.graph, &s);
        let z = mul(cx.graph, &[sn, x]);
        let m_function = |cx: &mut Cx<'_>, mu: &Q| -> Option<NodeId> {
            // exp(-z/2) z^(μ+1/2) 1F1(μ - κ + 1/2; 1 + 2μ; z).
            let minus_half = cx.graph.num(Number::fraction(-1, 2)?);
            let exponent = mul(cx.graph, &[minus_half, z]);
            let envelope = cx.graph.node(exp, &[exponent]);
            let e = node_of(cx.graph, &(mu + Q::new(BigInt::one(), BigInt::from(2))));
            let power = cx.graph.node(core::POW, &[z, e]);
            let (na, nb) = (node_of(cx.graph, &(mu - &kappa + Q::new(BigInt::one(), BigInt::from(2)))), node_of(cx.graph, &(q(1) + q(2) * mu)));
            let h = cx.graph.node(f, &[na, nb, z]);
            Some(mul(cx.graph, &[envelope, power, h]))
        };
        let two_mu = q(2) * &mu;
        if is_integer(&two_mu) != degenerate {
            return None;
        }
        if is_integer(&two_mu) {
            (m_function(cx, &mu)?, None)
        } else {
            (m_function(cx, &mu)?, Some(m_function(cx, &-mu.clone())?))
        }
    } else {
        if p0.len() > 1 || p1.len() > 2 || p2.len() < 2 || !coefficient(&p2, 0).is_zero() {
            return None;
        }
        let alpha = coefficient(&p2, 1);
        let (beta, gamma, delta) = (coefficient(&p1, 0), coefficient(&p1, 1), coefficient(&p0, 0));
        if p2.len() == 3 && coefficient(&p2, 2) == -alpha.clone() {
            // Gauss.
            let c = &beta / &alpha;
            if is_integer(&c) != degenerate {
                return None;
            }
            let s = -&gamma / &alpha - q(1);
            let disc = rational_sqrt(&(&s * &s + q(4) * &delta / &alpha))?;
            let (pa, pb) = ((&s + &disc) / q(2), (&s - &disc) / q(2));
            let f = cx.graph.ops().lookup("hyper2f1")?;
            let regular = |cx: &mut Cx<'_>| {
                let (na, nb, nc) = (node_of(cx.graph, &pa), node_of(cx.graph, &pb), node_of(cx.graph, &c));
                cx.graph.node(f, &[na, nb, nc, x])
            };
            let singular = |cx: &mut Cx<'_>| {
                let (a2, b2, c2n) = (&pa - &c + q(1), &pb - &c + q(1), q(2) - &c);
                let (na2, nb2, nc2) = (node_of(cx.graph, &a2), node_of(cx.graph, &b2), node_of(cx.graph, &c2n));
                let h = cx.graph.node(f, &[na2, nb2, nc2, x]);
                let power = power_of(cx, &(q(1) - &c));
                mul(cx.graph, &[power, h])
            };
            if is_integer(&c) {
                // 2F1(a, b; c; x) exists for c >= 1, the other form for c <= 0.
                (if c >= q(1) { regular(cx) } else { singular(cx) }, None)
            } else {
                (regular(cx), Some(singular(cx)))
            }
        } else if p2.len() == 2 && !gamma.is_zero() {
            // Kummer in z = s x.
            let s = -&gamma / &alpha;
            let b = &beta / &alpha;
            if is_integer(&b) != degenerate {
                return None;
            }
            let pa = &delta / &gamma;
            let f = cx.graph.ops().lookup("hyper1f1")?;
            let sn = node_of(cx.graph, &s);
            let z = mul(cx.graph, &[sn, x]);
            let regular = |cx: &mut Cx<'_>| {
                let (na, nb) = (node_of(cx.graph, &pa), node_of(cx.graph, &b));
                cx.graph.node(f, &[na, nb, z])
            };
            let singular = |cx: &mut Cx<'_>| {
                let (a2, b2) = (&pa - &b + q(1), q(2) - &b);
                let (na2, nb2) = (node_of(cx.graph, &a2), node_of(cx.graph, &b2));
                let h = cx.graph.node(f, &[na2, nb2, z]);
                let power = power_of(cx, &(q(1) - &b));
                mul(cx.graph, &[power, h])
            };
            if is_integer(&b) {
                (if b >= q(1) { regular(cx) } else { singular(cx) }, None)
            } else {
                (regular(cx), Some(singular(cx)))
            }
        } else {
            return None;
        }
    };
    let second = match second {
        | Some(y) => y,
        | None => super::reduce::reduction_of_order(cx, &a, first, x)?,
    };
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
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

/// The power series `P`, `Q` (`n + 2` terms) of the normal form
/// `t² y'' + t P y' + Q y = 0` at `t = 0`; `None` unless the point is
/// ordinary or regular singular.
fn normal_form_series(
    p: &[Vec<Q>; 3],
    n: usize,
) -> Option<(Vec<Q>, Vec<Q>)> {
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
    Some((series_quotient(&p_num, &s, terms)?, series_quotient(&q_num, &s, terms)?))
}

/// An element `a + b √D` of a quadratic field (the field's `D` is kept by
/// the caller).
#[derive(Clone)]
struct Quad {
    a: Q,
    b: Q,
}

impl Quad {
    fn rational(a: Q) -> Self {
        Self { a, b: Q::zero() }
    }

    fn add(&self, o: &Self) -> Self {
        Self { a: &self.a + &o.a, b: &self.b + &o.b }
    }

    fn mul(&self, o: &Self, d: &Q) -> Self {
        Self { a: &self.a * &o.a + &self.b * &o.b * d, b: &self.a * &o.b + &self.b * &o.a }
    }

    fn neg(&self) -> Self {
        Self { a: -self.a.clone(), b: -self.b.clone() }
    }

    fn div(&self, o: &Self, d: &Q) -> Option<Self> {
        let norm = &o.a * &o.a - &o.b * &o.b * d;
        if norm.is_zero() {
            return None;
        }
        let conjugate = Self { a: o.a.clone(), b: -o.b.clone() };
        let numerator = self.mul(&conjugate, d);
        Some(Self { a: numerator.a / &norm, b: numerator.b / &norm })
    }
}

/// A regular singular point whose indicial roots `σ ± τ √D` are irrational
/// (or complex, `D < 0`): the two solutions are conjugate, and the series of
/// the first has coefficients `u_k + v_k √D`.
pub(super) struct QuadraticData {
    d: Q,
    sigma: Q,
    tau: Q,
    series: Vec<(Q, Q)>,
}

fn frobenius_quadratic(
    p: &[Vec<Q>; 3],
    n: usize,
) -> Option<QuadraticData> {
    let (pp, qq) = normal_form_series(p, n)?;
    let (p_0, q_0) = (pp.first()?.clone(), qq.first()?.clone());
    let b = &p_0 - q(1);
    let disc = &b * &b - q(4) * &q_0;
    if disc.is_zero() || rational_sqrt(&disc).is_some() {
        return None;
    }
    let (sigma, tau) = (-&b / q(2), Q::new(BigInt::one(), BigInt::from(2)));
    let root = Quad { a: sigma.clone(), b: tau.clone() };
    let f = |r: &Quad| -> Quad {
        r.mul(r, &disc).add(&r.mul(&Quad::rational(b.clone()), &disc)).add(&Quad::rational(q_0.clone()))
    };
    let mut series: Vec<Quad> = vec![Quad::rational(q(1))];
    for k in 1..=n {
        let mut sum = Quad::rational(Q::zero());
        for j in 1..=k {
            let shift_by = Quad { a: &sigma + Q::from_integer(BigInt::from(k - j)), b: tau.clone() };
            let weight = shift_by.mul(&Quad::rational(pp.get(j)?.clone()), &disc).add(&Quad::rational(qq.get(j)?.clone()));
            sum = sum.add(&weight.mul(series.get(k - j)?, &disc));
        }
        let at = Quad { a: &root.a + Q::from_integer(BigInt::from(k)), b: root.b.clone() };
        series.push(sum.neg().div(&f(&at), &disc)?);
    }
    Some(QuadraticData { d: disc, sigma, tau, series: series.into_iter().map(|c| (c.a, c.b)).collect() })
}

/// Frobenius data at `t = 0` of `p2 y'' + p1 y' + p0 y = 0` (polynomials in
/// `t`), with `n + 1` terms.
pub(super) fn frobenius_at_zero(
    p: &[Vec<Q>; 3],
    n: usize,
) -> Option<Frobenius> {
    let [p0, p1, p2] = p;
    let order = |v: &[Q]| v.iter().position(|c| !c.is_zero());
    let (pp, qq) = normal_form_series(p, n)?;
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
    let non_integer = !is_integer(&difference);
    if non_integer {
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

/// Where to expand.
#[derive(Clone)]
enum Point {
    /// A finite point; `None` lets the search try `0` and the rational roots
    /// of the leading coefficient.
    Finite(Option<Q>),
    /// `x = ∞`, through `t = 1/x`.
    Infinity,
}

/// A formal solution `exp(λ/t) t^ρ Σ a_k t^k` at an irregular singular
/// point of rank one.
pub(super) struct Formal {
    lambda: Q,
    rho: Q,
    coefficients: Vec<Q>,
    /// `c_d(k)` of `Σ_d c_d(k) a(k - d) = 0`, ascending in `k`.
    recurrence: Vec<Vec<Q>>,
}

/// What is known at the chosen point.
pub(super) enum Local {
    /// Ordinary or regular singular: Frobenius.
    Regular { centre: Option<Q>, data: Frobenius },
    /// Regular singular with irrational or complex exponents.
    Quadratic { centre: Option<Q>, data: QuadraticData },
    /// Irregular singular: two formal solutions with exponential factors.
    Formal { centre: Option<Q>, items: Vec<Formal> },
}

/// The equation in `t = 1/x`: `P2 (t^4 Y'' + 2 t^3 Y') - P1 t^2 Y' + P0 Y`
/// with `P_k(t) = t^D p_k(1/t)`.
fn at_infinity(p: &[Vec<Q>; 3]) -> [Vec<Q>; 3] {
    let degree = p.iter().map(Vec::len).max().unwrap_or(1).saturating_sub(1);
    let reversed = |v: &[Q]| -> Vec<Q> {
        let mut out = vec![Q::zero(); degree + 1];
        for (i, c) in v.iter().enumerate() {
            if let Some(slot) = out.get_mut(degree - i) {
                *slot = c.clone();
            }
        }
        out
    };
    let times = |v: &[Q], k: usize| -> Vec<Q> { std::iter::repeat_n(Q::zero(), k).chain(v.iter().cloned()).collect() };
    let (r0, r1, r2) = (reversed(&p[0]), reversed(&p[1]), reversed(&p[2]));
    let b2 = times(&r2, 4);
    let b1 = univariate::sub(&univariate::mul(&[q(2)], &times(&r2, 3)), &times(&r1, 2));
    let b0 = r0;
    [univariate::sub(&b0, &[]), b1, b2]
}

/// Formal power-series solutions `t^ρ Σ a_k t^k` of `B2 w'' + B1 w' + B0 w = 0`
/// for each admissible rational `ρ`: the lowest power of `t` in
/// `Σ_d B_d D^d t^σ` fixes `ρ` as a root of the leading polynomial `Q(σ)`,
/// and `a_k Q(k + ρ) = -Σ_{j >= 1} Q_j(k + ρ - j) a_{k-j}`.
fn formal_series(
    b: &[Vec<Q>; 3],
    n: usize,
) -> Vec<(Q, Vec<Q>, Vec<Vec<Q>>)> {
    let coefficient = |d: usize, j: i64| -> Q { usize::try_from(j).ok().and_then(|j| b[d].get(j).cloned()).unwrap_or_else(Q::zero) };
    // Q_e(σ) = B0[e] + B1[e + 1] σ + B2[e + 2] (σ² - σ), ascending in σ.
    let q_e = |e: i64| -> Vec<Q> {
        let (c0, c1, c2) = (coefficient(0, e), coefficient(1, e + 1), coefficient(2, e + 2));
        trim_poly(vec![c0, &c1 - &c2, c2])
    };
    let top = b.iter().enumerate().map(|(d, v)| v.len() as i64 - 1 - d as i64).max().unwrap_or(0);
    let Some(e_min) = (-2..=top).find(|&e| !q_e(e).is_empty()) else { return Vec::new() };
    let leading = q_e(e_min);
    // Rational roots of the leading polynomial.
    let roots: Vec<Q> = match leading.as_slice() {
        | [c0, c1] => vec![-c0 / c1],
        | [c0, c1, c2] => match rational_sqrt(&(c1 * c1 - q(4) * c0 * c2)) {
            | Some(s) => {
                let (r1, r2) = ((-c1 + &s) / (q(2) * c2), (-c1 - &s) / (q(2) * c2));
                if r1 == r2 { vec![r1] } else { vec![r1, r2] }
            },
            | None => Vec::new(),
        },
        | _ => Vec::new(),
    };
    let depth = usize::try_from(top - e_min).unwrap_or(0);
    let mut out = Vec::new();
    for rho in roots {
        let mut a = vec![Q::one()];
        let mut ok = true;
        for k in 1..=n {
            let sigma = &rho + Q::from_integer(BigInt::from(k));
            let diagonal = univariate::eval(&leading, &sigma);
            if diagonal.is_zero() {
                ok = false;
                break;
            }
            let mut sum = Q::zero();
            for j in 1..=depth.min(k) {
                let e = e_min + i64::try_from(j).unwrap_or(0);
                let shifted = &sigma - Q::from_integer(BigInt::from(j));
                if let Some(ak) = a.get(k - j) {
                    sum += univariate::eval(&q_e(e), &shifted) * ak;
                }
            }
            a.push(-sum / diagonal);
        }
        if !ok {
            continue;
        }
        // Recurrence: coefficient of a(k - j) is Q_{e_min + j}(k + ρ - j), as a polynomial in k.
        let mut recurrence = Vec::new();
        for j in 0..=depth {
            let e = e_min + i64::try_from(j).unwrap_or(0);
            let offset = &rho - Q::from_integer(BigInt::from(j));
            recurrence.push(trim_poly(shift(&q_e(e), &offset)));
        }
        out.push((rho, a, recurrence));
    }
    out
}

fn trim_poly(mut p: Vec<Q>) -> Vec<Q> {
    trim(&mut p);
    p
}

/// Irregular singular point of rank one at `t = 0`: for each root `λ` of
/// `s0 λ² - p1[m-2] λ + p0[m-4]` the transformed equation for
/// `y = exp(λ/t) w` has a formal series solution.
fn formal_solutions(
    p: &[Vec<Q>; 3],
    n: usize,
) -> Option<Vec<Formal>> {
    let [p0, p1, p2] = p;
    let m = p2.iter().position(|c| !c.is_zero())?;
    let s0 = p2.get(m)?.clone();
    let at = |v: &[Q], index: i64| -> Q { usize::try_from(index).ok().and_then(|i| v.get(i).cloned()).unwrap_or_else(Q::zero) };
    let m_i = i64::try_from(m).ok()?;
    // Rank one: ord p1 >= m - 2 and ord p0 >= m - 4.
    if p1.iter().position(|c| !c.is_zero()).is_some_and(|o| (o as i64) < m_i - 2)
        || p0.iter().position(|c| !c.is_zero()).is_some_and(|o| (o as i64) < m_i - 4)
    {
        return None;
    }
    let (b1c, b0c) = (at(p1, m_i - 2), at(p0, m_i - 4));
    let disc = rational_sqrt(&(&b1c * &b1c - q(4) * &s0 * &b0c))?;
    let (l1, l2) = ((&b1c + &disc) / (q(2) * &s0), (&b1c - &disc) / (q(2) * &s0));
    if l1 == l2 {
        return None;
    }
    let times = |v: &[Q], k: usize| -> Vec<Q> { std::iter::repeat_n(Q::zero(), k).chain(v.iter().cloned()).collect() };
    let scaled = |v: &[Q], c: &Q| -> Vec<Q> { v.iter().map(|x| x * c).collect() };
    let mut out = Vec::new();
    for lambda in [l1, l2] {
        // B2 = t^4 p2; B1 = t^4 p1 - 2 λ t^2 p2; B0 = t^4 p0 - λ t^2 p1 + (2 λ t + λ²) p2.
        let b2 = times(p2, 4);
        let b1 = univariate::sub(&times(p1, 4), &scaled(&times(p2, 2), &(q(2) * &lambda)));
        let b0 = univariate::add(
            &univariate::sub(&times(p0, 4), &scaled(&times(p1, 2), &lambda)),
            &univariate::add(&scaled(&times(p2, 1), &(q(2) * &lambda)), &scaled(p2, &(&lambda * &lambda))),
        );
        let (rho, coefficients, recurrence) = formal_series(&[b0, b1, b2], n).into_iter().next()?;
        out.push(Formal { lambda, rho, coefficients, recurrence });
    }
    Some(out)
}

/// Frobenius data, or failing that the formal solutions, at the chosen point.
fn local(
    p: &[Vec<Q>; 3],
    n: usize,
    point: &Point,
) -> Option<Local> {
    match point {
        | Point::Infinity => {
            let transformed = at_infinity(p);
            if let Some(data) = frobenius_at_zero(&transformed, n) {
                return Some(Local::Regular { centre: None, data });
            }
            if let Some(data) = frobenius_quadratic(&transformed, n) {
                return Some(Local::Quadratic { centre: None, data });
            }
            formal_solutions(&transformed, n).map(|items| Local::Formal { centre: None, items })
        },
        | Point::Finite(preferred) => {
            let mut points: Vec<Q> = Vec::new();
            if let Some(x0) = preferred {
                points.push(x0.clone());
            } else {
                points.push(Q::zero());
                for r in univariate::rational_roots(&p[2]) {
                    if !points.contains(&r) {
                        points.push(r);
                    }
                }
            }
            let shifted_at = |x0: &Q| [shift(&p[0], x0), shift(&p[1], x0), shift(&p[2], x0)];
            for x0 in &points {
                if let Some(data) = frobenius_at_zero(&shifted_at(x0), n) {
                    return Some(Local::Regular { centre: Some(x0.clone()), data });
                }
            }
            for x0 in &points {
                if let Some(data) = frobenius_quadratic(&shifted_at(x0), n) {
                    return Some(Local::Quadratic { centre: Some(x0.clone()), data });
                }
            }
            for x0 in &points {
                if let Some(items) = formal_solutions(&shifted_at(x0), n) {
                    return Some(Local::Formal { centre: Some(x0.clone()), items });
                }
            }
            None
        },
    }
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

/// The local variable `t` as a term.
fn local_variable(
    cx: &mut Cx<'_>,
    x: NodeId,
    centre: Option<&Q>,
) -> NodeId {
    match centre {
        | None => {
            let minus_one = cx.graph.int(-1);
            cx.graph.node(core::POW, &[x, minus_one])
        },
        | Some(x0) if x0.is_zero() => x,
        | Some(x0) => {
            let h = node_of(cx.graph, &-x0.clone());
            cx.graph.node(core::ADD, &[x, h])
        },
    }
}

/// The two solutions as terms.
fn solutions(
    cx: &mut Cx<'_>,
    x: NodeId,
    local: &Local,
    n: usize,
) -> Option<(NodeId, NodeId)> {
    match local {
        | Local::Regular { centre, data } => {
            let t = local_variable(cx, x, centre.as_ref());
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
        },
        | Local::Quadratic { centre, data } => {
            let t = local_variable(cx, x, centre.as_ref());
            let negative = data.d.is_negative();
            let magnitude = if negative { -data.d.clone() } else { data.d.clone() };
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let magnitude_node = node_of(cx.graph, &magnitude);
            let root = cx.graph.node(core::POW, &[magnitude_node, half]);
            let order_symbol = cx.graph.ops().lookup("big_o")?;
            let sigma_node = node_of(cx.graph, &data.sigma);
            let tau_node = node_of(cx.graph, &data.tau);
            let terms_of = |cx: &mut Cx<'_>, coefficients: Vec<NodeId>| -> Option<NodeId> {
                let mut terms = Vec::new();
                for (k, c) in coefficients.into_iter().enumerate() {
                    let power = if k == 0 {
                        cx.graph.int(1)
                    } else {
                        let e = cx.graph.int(i64::try_from(k).ok()?);
                        cx.graph.node(core::POW, &[t, e])
                    };
                    terms.push(mul(cx.graph, &[c, power]));
                }
                let e = cx.graph.int(i64::try_from(n + 1).ok()?);
                let tail = cx.graph.node(core::POW, &[t, e]);
                terms.push(cx.graph.node(order_symbol, &[tail]));
                Some(add(cx.graph, &terms))
            };
            if negative {
                // r = σ + i τ √m; y = t^σ e^{iθ} (A + i B), θ = τ √m ln t.
                let ln = cx.graph.ops().lookup("ln")?;
                let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
                let a_coefficients: Vec<NodeId> = data.series.iter().map(|(u, _)| node_of(cx.graph, u)).collect();
                let b_coefficients: Vec<NodeId> = data
                    .series
                    .iter()
                    .map(|(_, v)| {
                        let vn = node_of(cx.graph, v);
                        mul(cx.graph, &[vn, root])
                    })
                    .collect();
                let (big_a, big_b) = (terms_of(cx, a_coefficients)?, terms_of(cx, b_coefficients)?);
                let log = cx.graph.node(ln, &[t]);
                let theta = mul(cx.graph, &[tau_node, root, log]);
                let (s_node, c_node) = (cx.graph.node(sin, &[theta]), cx.graph.node(cos, &[theta]));
                let sigma_power = cx.graph.node(core::POW, &[t, sigma_node]);
                let neg_b = neg(cx.graph, big_b);
                let neg_a = neg(cx.graph, big_a);
                let _ = neg_a;
                let real_inner = {
                    let first = mul(cx.graph, &[c_node, big_a]);
                    let second = mul(cx.graph, &[s_node, neg_b]);
                    add(cx.graph, &[first, second])
                };
                let imaginary_inner = {
                    let first = mul(cx.graph, &[s_node, big_a]);
                    let second = mul(cx.graph, &[c_node, big_b]);
                    add(cx.graph, &[first, second])
                };
                return Some((mul(cx.graph, &[sigma_power, real_inner]), mul(cx.graph, &[sigma_power, imaginary_inner])));
            }
            // Real irrational exponents σ ± τ √D.
            let mut pair = Vec::new();
            for sign in [1_i64, -1] {
                let sg = cx.graph.int(sign);
                let offset = mul(cx.graph, &[sg, tau_node, root]);
                let exponent = add(cx.graph, &[sigma_node, offset]);
                let coefficients: Vec<NodeId> = data
                    .series
                    .iter()
                    .map(|(u, v)| {
                        let (un, vn) = (node_of(cx.graph, u), node_of(cx.graph, v));
                        let signed = mul(cx.graph, &[sg, vn, root]);
                        add(cx.graph, &[un, signed])
                    })
                    .collect();
                let polynomial = terms_of(cx, coefficients)?;
                let power = cx.graph.node(core::POW, &[t, exponent]);
                pair.push(mul(cx.graph, &[power, polynomial]));
            }
            let (&first, &second) = (pair.first()?, pair.get(1)?);
            Some((first, second))
        },
        | Local::Formal { centre, items } => {
            let t = local_variable(cx, x, centre.as_ref());
            let exp = cx.graph.ops().lookup("exp")?;
            let mut nodes = Vec::new();
            for item in items {
                let series = series_node(cx, t, &item.coefficients, &item.rho, n)?;
                if item.lambda.is_zero() {
                    nodes.push(series);
                    continue;
                }
                let minus_one = cx.graph.int(-1);
                let inverse = cx.graph.node(core::POW, &[t, minus_one]);
                let lambda = node_of(cx.graph, &item.lambda);
                let exponent = mul(cx.graph, &[lambda, inverse]);
                let factor = cx.graph.node(exp, &[exponent]);
                nodes.push(mul(cx.graph, &[factor, series]));
            }
            let (&first, &second) = (nodes.first()?, nodes.get(1)?);
            Some((first, second))
        },
    }
}

/// `dsolve` fallback: the truncated Frobenius or formal series (eight terms).
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
    let found = local(&p, TERMS, &Point::Finite(None))?;
    let (y1, y2) = solutions(cx, x, &found, TERMS)?;
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
    let solution = add(cx.graph, &terms);
    Some(explicit(cx, problem, solution))
}

/// The recurrence `Σ_d c_d(k) a(k - d) = 0` as an equation.
fn recurrence_node(
    cx: &mut Cx<'_>,
    recurrence: &[Vec<Q>],
) -> Option<NodeId> {
    let k_symbol = cx.graph.interner_mut().fresh_symbol("k");
    let k = cx.graph.symbol_node(k_symbol);
    let a_symbol = cx.graph.interner_mut().fresh_symbol("a");
    let a_fn = cx.graph.symbol_node(a_symbol);
    let mut parts = Vec::new();
    for (d, coefficients) in recurrence.iter().enumerate() {
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
        if terms.is_empty() {
            continue;
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
    let lhs = add(cx.graph, &parts);
    let zero = cx.graph.int(0);
    Some(cx.graph.node(core::EQ, &[lhs, zero]))
}

/// `ode_frobenius(equation, y(x), x0, n)`; `x0` may be `oo`.
pub(super) fn report(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    x0: NodeId,
    n: NodeId,
) -> Option<NodeId> {
    let n = usize::try_from(cx.graph.number_of(n)?.to_i64()?).ok().filter(|&n| n <= 60)?;
    let infinity = cx.graph.ops().lookup("oo").is_some_and(|op| cx.graph.op(x0) == op);
    let point = if infinity { Point::Infinity } else { Point::Finite(Some(cx.graph.number_of(x0)?.to_rational()?)) };
    let mut problem = super::parse(cx.graph, equation, unknown)?;
    let (a, forcing) = linear_second_order(cx, &problem)?;
    if !cx.is_zero(forcing) {
        return None;
    }
    let x = problem.x;
    let p = polynomials(cx, &a, x)?;
    let found = local(&p, n, &point)?;
    let (y1, y2) = solutions(cx, x, &found, n)?;
    let (c1, c2) = (problem.constant(cx.graph), problem.constant(cx.graph));
    let terms = [mul(cx.graph, &[c1, y1]), mul(cx.graph, &[c2, y2])];
    let solution = add(cx.graph, &terms);
    let answer = explicit(cx, &problem, solution);
    let (exponents, recurrence) = match &found {
        | Local::Regular { data, .. } => {
            let (r1, r2) = (node_of(cx.graph, &data.roots.0), node_of(cx.graph, &data.roots.1));
            (cx.graph.node(core::LIST, &[r1, r2]), recurrence_node(cx, &data.recurrence)?)
        },
        | Local::Quadratic { data, .. } => {
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let negative = data.d.is_negative();
            let magnitude = if negative { -data.d.clone() } else { data.d.clone() };
            let magnitude_node = node_of(cx.graph, &magnitude);
            let root = cx.graph.node(core::POW, &[magnitude_node, half]);
            let mut factors = vec![node_of(cx.graph, &data.tau), root];
            if negative {
                factors.push(cx.graph.node(cx.graph.ops().lookup("I")?, &[]));
            }
            let offset = mul(cx.graph, &factors);
            let sigma = node_of(cx.graph, &data.sigma);
            let minus_offset = neg(cx.graph, offset);
            let (r1, r2) = (add(cx.graph, &[sigma, offset]), add(cx.graph, &[sigma, minus_offset]));
            (cx.graph.node(core::LIST, &[r1, r2]), cx.graph.node(core::LIST, &[]))
        },
        | Local::Formal { items, .. } => {
            let mut pairs = Vec::new();
            for item in items {
                let (l, r) = (node_of(cx.graph, &item.lambda), node_of(cx.graph, &item.rho));
                pairs.push(cx.graph.node(core::LIST, &[l, r]));
            }
            let first = items.first()?;
            (cx.graph.node(core::LIST, &pairs), recurrence_node(cx, &first.recurrence)?)
        },
    };
    let _ = (neg, fresh_function);
    Some(cx.graph.node(core::LIST, &[answer, exponents, recurrence]))
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
