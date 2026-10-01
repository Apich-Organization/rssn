//! Series, sums and products.
//!
//! * `taylor(f, x, a, n)` and `laurent(f, x, a, n)`: expansions about `a`
//!   up to `(x - a)^n`, returned in that form (pinned). `laurent(f, x, oo,
//!   n)` (alias `asymptotic(f, x, n)`) expands in powers of `1/x`.
//! * `fourier_series(f, x, L, n)`: the first `n` harmonics of `f` on
//!   `[-L, L]`.
//! * `sum(f, k, a, b)` and `product(f, k, a, b)`: finite sums are written
//!   out; polynomial, geometric and a few classical infinite summands have
//!   closed forms; anything else is summed numerically, with convergence
//!   acceleration, when a number is wanted.
//! * `converges(f, k)`: divergence, ratio, root, alternating, p-series and
//!   condensation tests for the series with general term `f`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::backend::Backend;
use crate::backend::Interpreter;
use crate::graph::op::core;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::kernels::series::sum_range;
use crate::kernels::series::sum_to_infinity;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;

use super::diff::Differentiate;
use super::integrate::Functions;
use super::limits::limit_at;
use super::limits::Side;
use super::powerseries::laurent_expansion;

/// Operators shared by the kernels of this module.
#[derive(Copy, Clone)]
pub(super) struct SeriesOps {
    pub(super) functions: Functions,
    pub(super) infinity: OpId,
    pub(super) defint: OpId,
    pub(super) taylor: OpId,
    pub(super) laurent: OpId,
    pub(super) fourier: OpId,
    pub(super) sum: OpId,
    pub(super) product: OpId,
    pub(super) converges: OpId,
}

fn factorial(n: usize) -> BigInt {
    (1..=n).fold(BigInt::one(), |acc, k| acc * BigInt::from(k))
}

fn mul(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn add(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn sub(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    let negated = graph.node(core::MUL, &[minus_one, b]);
    graph.node(core::ADD, &[a, negated])
}

fn is_zero(
    graph: &Graph,
    node: NodeId,
) -> bool {
    graph.number_of(node).is_some_and(Number::is_zero)
}

/// Coefficients `c_0..=c_n` of the Taylor expansion of `f` about `a`.
fn taylor_coefficients(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    f: NodeId,
    x: NodeId,
    a: NodeId,
    n: usize,
) -> Option<Vec<NodeId>> {
    let symbol = cx.graph.symbol_of(x)?;
    let mut derivative = cx.simplify(f);
    let mut out = Vec::with_capacity(n + 1);
    for k in 0..=n {
        // The value of the k-th derivative at the point, as a limit so
        // that removable singularities (sin(x)/x at 0) are handled.
        let value = limit_at(cx, ops.functions, ops.infinity, derivative, x, a, Side::Both)?;
        if cx.graph.eval(value, &Env::numeric(0.0)).is_some_and(f64::is_infinite) {
            return None;
        }
        let scale = cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial(k))));
        let coefficient = mul(cx.graph, &[scale, value]);
        out.push(cx.simplify(coefficient));
        if k < n {
            let raw = Differentiate { diff: ops.functions.diff }.derive(cx.graph, derivative, symbol, x);
            derivative = cx.simplify(raw);
        }
    }
    Some(out)
}

/// `sum_k coefficients[k] * (x - a)^(k + shift)`.
fn power_series(
    graph: &mut Graph,
    coefficients: &[NodeId],
    x: NodeId,
    a: NodeId,
    shift: i64,
) -> NodeId {
    let offset = if is_zero(graph, a) { x } else { sub(graph, x, a) };
    let mut terms = Vec::new();
    for (k, &c) in coefficients.iter().enumerate() {
        if is_zero(graph, c) {
            continue;
        }
        let power = i64::try_from(k).unwrap_or(i64::MAX).saturating_add(shift);
        let unit = graph.number_of(c).is_some_and(Number::is_one);
        let raised = match power {
            | 0 => None,
            | 1 => Some(offset),
            | p => {
                let e = graph.int(p);
                Some(graph.node(core::POW, &[offset, e]))
            },
        };
        let term = match raised {
            | None => c,
            | Some(r) if unit => r,
            | Some(r) => mul(graph, &[c, r]),
        };
        terms.push(term);
    }
    add(graph, &terms)
}

fn order(
    graph: &Graph,
    n: NodeId,
) -> Option<usize> {
    usize::try_from(graph.number_of(n)?.to_i64()?).ok().filter(|&n| n <= 40)
}

fn taylor(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, x, a, n] = args else {
        return None;
    };
    let n = order(cx.graph, n)?;
    let a = best(cx.graph, a)?;
    let order_i = i64::try_from(n).ok()?;
    if let Some((valuation, coefficients)) = laurent_expansion(cx, f, x, a, order_i) {
        // A Taylor series has no negative powers.
        if valuation < 0 {
            return None;
        }
        return Some(power_series(cx.graph, &coefficients, x, a, valuation));
    }
    let coefficients = taylor_coefficients(cx, ops, f, x, a, n)?;
    Some(power_series(cx.graph, &coefficients, x, a, 0))
}

/// The expansion of `f` in powers of `1/x` as `x → ∞`, through
/// `x^(-n)`: `f(1/t)` is expanded about `t = 0` (a positive `t`, so
/// `sqrt(t^2) = t`) and `t = 1/x` substituted back.
fn at_infinity(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    n: usize,
) -> Option<NodeId> {
    let fresh = cx.graph.interner_mut().fresh_symbol("t");
    cx.graph.assume(fresh, crate::graph::Facts::POSITIVE);
    let t = cx.graph.symbol_node(fresh);
    let minus_one = cx.graph.int(-1);
    let inverse_t = cx.graph.node(core::POW, &[t, minus_one]);
    let substituted = cx.graph.substitute(f, x, inverse_t);
    let g = cx.simplify(substituted);
    let zero = cx.graph.int(0);
    let (valuation, coefficients) = laurent_expansion(cx, g, t, zero, i64::try_from(n).ok()?)?;
    let series = power_series(cx.graph, &coefficients, t, zero, valuation);
    let inverse_x = cx.graph.node(core::POW, &[x, minus_one]);
    let back = cx.graph.substitute(series, t, inverse_x);
    Some(cx.simplify(back))
}

fn laurent(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, x, a, n] = args else {
        return None;
    };
    let n = order(cx.graph, n)?;
    let a = best(cx.graph, a)?;
    if cx.graph.op(a) == ops.infinity {
        return at_infinity(cx, f, x, n);
    }
    if let Some((valuation, coefficients)) = laurent_expansion(cx, f, x, a, i64::try_from(n).ok()?) {
        return Some(power_series(cx.graph, &coefficients, x, a, valuation));
    }
    // Fallback: the order of the pole by limits, then a Taylor series of
    // the regular part by differentiation.
    let offset = if is_zero(cx.graph, a) { x } else { sub(cx.graph, x, a) };
    for m in 0..=8_usize {
        let regular = if m == 0 {
            f
        } else {
            let e = cx.graph.int(i64::try_from(m).ok()?);
            let factor = cx.graph.node(core::POW, &[offset, e]);
            mul(cx.graph, &[factor, f])
        };
        let Some(value) = limit_at(cx, ops.functions, ops.infinity, regular, x, a, Side::Both) else {
            continue;
        };
        if cx.graph.eval(value, &Env::numeric(0.0)).is_some_and(f64::is_infinite) {
            continue;
        }
        let coefficients = taylor_coefficients(cx, ops, regular, x, a, n + m)?;
        let shift = -i64::try_from(m).ok()?;
        return Some(power_series(cx.graph, &coefficients, x, a, shift));
    }
    None
}

/// `a_0/2 + sum_{k=1}^{n} a_k cos(k pi x / L) + b_k sin(k pi x / L)` with
/// the coefficients left as definite-integral requests for the engine.
fn fourier(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, x, half_period, n] = args else {
        return None;
    };
    let n = order(cx.graph, n)?;
    let graph = &mut *cx.graph;
    let pi = graph.ops().lookup("pi")?;
    let pi = graph.node(pi, &[]);
    let minus_one = graph.int(-1);
    let lower = graph.node(core::MUL, &[minus_one, half_period]);
    let inverse_l = graph.node(core::POW, &[half_period, minus_one]);
    let mut terms = Vec::new();
    let mean = graph.node(ops.defint, &[f, x, lower, half_period]);
    let half = graph.num(Number::fraction(1, 2)?);
    terms.push(mul(graph, &[half, inverse_l, mean]));
    for k in 1..=n {
        let k_node = graph.int(i64::try_from(k).ok()?);
        let angle = mul(graph, &[k_node, pi, x, inverse_l]);
        for wave in [ops.functions.cos, ops.functions.sin] {
            let basis = graph.node(wave, &[angle]);
            let integrand = mul(graph, &[f, basis]);
            let coefficient = graph.node(ops.defint, &[integrand, x, lower, half_period]);
            terms.push(mul(graph, &[inverse_l, coefficient, basis]));
        }
    }
    Some(add(graph, &terms))
}

/// Stirling numbers of the second kind `S(p, 0..=p)`.
fn stirling2_row(p: usize) -> Vec<BigInt> {
    let mut row = vec![BigInt::one()];
    for n in 1..=p {
        let mut next = vec![BigInt::zero(); n + 1];
        for k in 1..=n {
            let stay = row.get(k).cloned().unwrap_or_else(BigInt::zero) * BigInt::from(k);
            next[k] = stay + row[k - 1].clone();
        }
        row = next;
    }
    row
}

/// `sum_{k=0}^{n} k^p` as a term in `n`, through falling factorials:
/// `k^p = sum_j S(p, j) k^(j)` and `sum_{k=0}^{n} k^(j) = (n+1)^(j+1)/(j+1)`.
fn power_sum(
    graph: &mut Graph,
    p: usize,
    n: NodeId,
) -> NodeId {
    let row = stirling2_row(p);
    let one = graph.int(1);
    let n_plus_one = graph.node(core::ADD, &[n, one]);
    let mut terms = Vec::new();
    for (j, s) in row.iter().enumerate() {
        if s.is_zero() {
            continue;
        }
        // (n + 1)(n)(n - 1)...(n - j + 1): j + 1 factors.
        let mut factors = Vec::with_capacity(j + 2);
        factors.push(graph.num(Number::rat(BigRational::new(s.clone(), BigInt::from(j + 1)))));
        for i in 0..=j {
            let shift = graph.int(-i64::try_from(i).unwrap_or(0));
            factors.push(graph.node(core::ADD, &[n_plus_one, shift]));
        }
        terms.push(mul(graph, &factors));
    }
    add(graph, &terms)
}

/// Closed form of `sum_{k=lower}^{upper} f`, if `f` is a polynomial in `k`
/// or geometric.
fn symbolic_sum(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(k)?;
    let term = cx.simplify(f);
    let graph = &mut *cx.graph;
    let depends = |graph: &Graph, node: NodeId| graph.depends_on(graph.find(node), symbol);
    let infinite_upper = graph.eval(upper, &Env::numeric(0.0)).is_some_and(|v| v == f64::INFINITY);

    // Polynomial in k: Faulhaber.
    if !infinite_upper {
        let mut gens = Gens::default();
        let gk = gens.index(graph, k);
        if let Some(poly) = from_term(graph, &mut gens, term, Limits::default()) {
            let pure = poly.support().iter().all(|&g| g == gk || gens.node(g).is_some_and(|n| !depends(graph, n)));
            if pure {
                let one = graph.int(1);
                let before = sub(graph, lower, one);
                let mut pieces = Vec::new();
                for (p, coefficient) in poly.coefficients_in(gk).iter().enumerate() {
                    if coefficient.is_zero() {
                        continue;
                    }
                    let c = to_term(graph, &gens, coefficient);
                    let up_to_upper = power_sum(graph, p, upper);
                    let up_to_before = power_sum(graph, p, before);
                    let span = sub(graph, up_to_upper, up_to_before);
                    pieces.push(mul(graph, &[c, span]));
                }
                return Some(add(graph, &pieces));
            }
        }
    }

    // c * r^k with c, r free of k.
    let factors: Vec<NodeId> = if graph.op(term) == core::MUL { graph.children(term).to_vec() } else { vec![term] };
    let (constant, varying): (Vec<NodeId>, Vec<NodeId>) = factors.iter().partition(|&&n| !depends(graph, n));
    if let &[power] = varying.as_slice() {
        if let (true, &[ratio, exponent]) = (graph.op(power) == core::POW, graph.children(power)) {
            if exponent == k && !depends(graph, ratio) {
                let one = graph.int(1);
                let start = graph.node(core::POW, &[ratio, lower]);
                let denominator = sub(graph, one, ratio);
                let minus_one = graph.int(-1);
                let inverse = graph.node(core::POW, &[denominator, minus_one]);
                let numerator = if infinite_upper {
                    // Converges only for |r| < 1, which must be known.
                    let r = graph.eval(ratio, &Env::numeric(0.0))?;
                    if r.abs() >= 1.0 {
                        return None;
                    }
                    start
                } else {
                    let next = graph.node(core::ADD, &[upper, one]);
                    let end = graph.node(core::POW, &[ratio, next]);
                    sub(graph, start, end)
                };
                let mut all = constant;
                all.push(numerator);
                all.push(inverse);
                return Some(mul(graph, &all));
            }
            // 1/k^s from 1 to infinity is zeta(s).
            let from_one = graph.number_of(lower).is_some_and(Number::is_one);
            if infinite_upper && from_one && ratio == k && !depends(graph, exponent) {
                let zeta = graph.ops().lookup("zeta")?;
                let s = graph.number_of(exponent)?.neg();
                if s.to_f64() <= 1.0 {
                    return None;
                }
                let s = graph.num(s);
                let value = graph.node(zeta, &[s]);
                let mut all = constant;
                all.push(value);
                return Some(mul(graph, &all));
            }
        }
    }
    // Hypergeometric terms: Gosper's antidifference; rational terms that
    // are not Gosper-summable: partial fractions and polygamma.
    if infinite_upper {
        return rational_sum(cx, term, k, lower, None);
    }
    let from_nonnegative = cx.graph.number_of(lower).is_some_and(|n| !n.is_negative());
    let Some(big_t) = super::gosper::antidifference(cx, term, k, from_nonnegative) else {
        if let Some(found) = rational_sum(cx, term, k, lower, Some(upper)) {
            return Some(found);
        }
        // Definite sums such as Σ binomial(n, k) x^k: a fitted recurrence.
        return super::recurrence::sum_by_recurrence(cx, term, k, lower, upper);
    };
    let graph = &mut *cx.graph;
    let one = graph.int(1);
    let after = graph.node(core::ADD, &[upper, one]);
    let at_end = graph.substitute(big_t, k, after);
    let at_start = graph.substitute(big_t, k, lower);
    let difference = sub(graph, at_end, at_start);
    Some(cx.simplify(difference))
}

/// `sum_{k=lower}^{upper} r(k)` for a rational function `r` over `Q` whose
/// denominator splits into linear factors, by partial fractions:
/// `sum 1/(k + a)^m` is a difference of polygamma values (harmonic numbers
/// for `m = 1` and integer `a`). `upper = None` sums to infinity, which
/// requires the residues of the simple poles to cancel.
fn rational_sum(
    cx: &mut Cx<'_>,
    term: NodeId,
    k: NodeId,
    lower: NodeId,
    upper: Option<NodeId>,
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let (numer, denom) = crate::rules::poly::rational_function_in(graph, term, k)?;
    if denom.len() < 2 {
        return None;
    }
    let parts = crate::rules::poly::apart::apart(&numer, &denom)?;
    if !parts.quotient.is_empty() && upper.is_none() {
        return None;
    }
    let polygamma = graph.ops().lookup("polygamma")?;
    let digamma = graph.ops().lookup("digamma")?;
    let harmonic = graph.ops().lookup("harmonic");
    let one = graph.int(1);
    // Simple poles must have residues summing to zero for convergence.
    let mut residue_sum = BigRational::zero();
    let mut pieces = Vec::new();
    for piece in &parts.pieces {
        // factor = p1 k + p0 = p1 (k + alpha); numerator c' constant.
        let ([p0, p1], [c]) = (piece.factor.as_slice(), piece.numerator.as_slice()) else {
            return None;
        };
        let alpha = p0 / p1;
        let mut scale = BigRational::one();
        for _ in 0..piece.power {
            scale /= p1;
        }
        let c = c * scale;
        let m = i64::from(piece.power);
        if m == 1 {
            residue_sum += &c;
        }
        let alpha_term = graph.num(Number::rat(alpha.clone()));
        let start = graph.node(core::ADD, &[lower, alpha_term]);
        let integer_shift = alpha.is_integer();
        // S(z) = sum_{j >= 0} 1/(z + j)^m, up to a constant for m = 1.
        let tail = |graph: &mut Graph, z: NodeId| -> NodeId {
            if m == 1 {
                // -digamma(z), or -harmonic(z - 1) for an integer shift.
                let value = match harmonic {
                    | Some(h) if integer_shift => {
                        let minus_one = graph.int(-1);
                        let before = graph.node(core::ADD, &[z, minus_one]);
                        graph.node(h, &[before])
                    },
                    | _ => graph.node(digamma, &[z]),
                };
                let minus_one = graph.int(-1);
                graph.node(core::MUL, &[minus_one, value])
            } else {
                let order = graph.int(m - 1);
                let value = graph.node(polygamma, &[order, z]);
                let mut factorial = BigInt::one();
                for j in 2..m {
                    factorial *= BigInt::from(j);
                }
                let sign = if m % 2 == 0 { BigInt::one() } else { -BigInt::one() };
                let coefficient = graph.num(Number::rat(BigRational::new(sign, factorial)));
                graph.node(core::MUL, &[coefficient, value])
            }
        };
        let from = tail(graph, start);
        let span = match upper {
            | Some(b) => {
                let after = graph.node(core::ADD, &[b, alpha_term, one]);
                let to = tail(graph, after);
                sub(graph, from, to)
            },
            | None => from,
        };
        let c = graph.num(Number::rat(c));
        pieces.push(graph.node(core::MUL, &[c, span]));
    }
    if upper.is_none() && !residue_sum.is_zero() {
        return None;
    }
    if let Some(b) = upper {
        if !parts.quotient.is_empty() {
            let before = sub(graph, lower, one);
            for (p, coefficient) in parts.quotient.iter().enumerate() {
                if coefficient.is_zero() {
                    continue;
                }
                let c = graph.num(Number::rat(coefficient.clone()));
                let up_to_upper = power_sum(graph, p, b);
                let up_to_before = power_sum(graph, p, before);
                let span = sub(graph, up_to_upper, up_to_before);
                pieces.push(graph.node(core::MUL, &[c, span]));
            }
        }
    }
    let total = add(graph, &pieces);
    Some(cx.simplify(total))
}

/// The integer bounds of a sum or product if both are literal and the
/// range is small enough to write out.
fn literal_range(
    graph: &Graph,
    lower: NodeId,
    upper: NodeId,
) -> Option<(i64, i64)> {
    let (a, b) = (graph.number_of(lower)?.to_i64()?, graph.number_of(upper)?.to_i64()?);
    (b.saturating_sub(a) <= 10_000).then_some((a, b))
}

fn sum(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, k, lower, upper] = args else {
        return None;
    };
    let (lower, upper) = (best(cx.graph, lower)?, best(cx.graph, upper)?);
    if let Some((a, b)) = literal_range(cx.graph, lower, upper) {
        let body = best(cx.graph, f)?;
        let terms: Vec<NodeId> = (a..=b)
            .map(|value| {
                let v = cx.graph.int(value);
                cx.graph.substitute(body, k, v)
            })
            .collect();
        return Some(add(cx.graph, &terms));
    }
    symbolic_sum(cx, f, k, lower, upper)
}

fn product(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, k, lower, upper] = args else {
        return None;
    };
    let (lower, upper) = (best(cx.graph, lower)?, best(cx.graph, upper)?);
    let body = best(cx.graph, f)?;
    let graph = &mut *cx.graph;
    if let Some((a, b)) = literal_range(graph, lower, upper) {
        let factors: Vec<NodeId> = (a..=b)
            .map(|value| {
                let v = graph.int(value);
                graph.substitute(body, k, v)
            })
            .collect();
        return Some(mul(graph, &factors));
    }
    let symbol = graph.symbol_of(k)?;
    let one = graph.int(1);
    if !graph.depends_on(graph.find(body), symbol) {
        // c^(upper - lower + 1)
        let span = sub(graph, upper, lower);
        let count = graph.node(core::ADD, &[span, one]);
        return Some(graph.node(core::POW, &[body, count]));
    }
    if body == k {
        // upper! / (lower - 1)!
        let factorial = graph.ops().lookup("factorial")?;
        let top = graph.node(factorial, &[upper]);
        let before = sub(graph, lower, one);
        let bottom = graph.node(factorial, &[before]);
        let minus_one = graph.int(-1);
        let inverse = graph.node(core::POW, &[bottom, minus_one]);
        return Some(mul(graph, &[top, inverse]));
    }
    None
}

/// The numeric value of `lim_{k→∞} e`, if the limit engine finds one.
fn limit_value(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    e: NodeId,
    k: NodeId,
) -> Option<f64> {
    let infinity = cx.graph.node(ops.infinity, &[]);
    let limit = limit_at(cx, ops.functions, ops.infinity, e, k, infinity, Side::Both)?;
    cx.graph.eval(limit, &Env::numeric(0.0))
}

/// Whether `sum_k a_k` converges, by a battery of tests in order of cost:
///
/// 1. divergence test: `a_k ↛ 0` diverges;
/// 2. ratio test and 3. root test (decisive away from 1);
/// 4. alternating series (Leibniz): `(-1)^k b_k` with `b_k ↓ 0`;
/// 5. limit comparison with `k^-p`, `p = lim -ln|a_k| / ln k`;
/// 6. Cauchy condensation for terms with logarithms
///    (`sum a_k` ~ `sum 2^k a_{2^k}` for decreasing positive terms),
///    which decides the borderline `p = 1` cases such as `1/(k ln(k)^2)`.
fn converges(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, k] = args else {
        return None;
    };
    // The index is a positive integer, which lets `abs` and `ln` simplify.
    let fresh = cx.graph.interner_mut().fresh_symbol("k");
    cx.graph.assume(fresh, crate::graph::Facts::POSITIVE | crate::graph::Facts::INTEGER);
    let index = cx.graph.symbol_node(fresh);
    let f = cx.graph.substitute(f, k, index);
    let verdict = convergence(cx, ops, f, index, 0)?;
    Some(cx.graph.lit(Payload::Bool(verdict)))
}

fn convergence(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    f: NodeId,
    k: NodeId,
    depth: u32,
) -> Option<bool> {
    let term = cx.simplify(f);
    let abs = cx.graph.ops().lookup("abs")?;
    let ln = cx.graph.ops().lookup("ln")?;
    // |a_k|: the term itself (or its negative) when its sign is eventually
    // fixed, which keeps it simplifiable; `abs` otherwise.
    let same_sign = eventually_one_sign(cx, term, k);
    let magnitude = match sample_sign(cx, term, k) {
        | Some(true) if same_sign => term,
        | Some(false) if same_sign => {
            let minus_one = cx.graph.int(-1);
            let negated = cx.graph.node(core::MUL, &[minus_one, term]);
            cx.simplify(negated)
        },
        | _ => cx.graph.node(abs, &[term]),
    };

    // 1. Divergence test.
    if let Some(v) = limit_value(cx, ops, magnitude, k) {
        if v.is_nan() || v > 1e-12 {
            return (!v.is_nan()).then_some(false);
        }
    }

    // 2./3. Ratio and root tests. Both quantities are first estimated
    // numerically far out, which decides clear cases cheaply; the symbolic
    // limit is taken only when that is inconclusive. A term whose base and
    // exponent both vary with k (`(k/(2k+1))^k`) goes to the root test
    // first, where it simplifies.
    let graph = &mut *cx.graph;
    let one = graph.int(1);
    let minus_one = graph.int(-1);
    let next_index = graph.node(core::ADD, &[k, one]);
    let next = graph.substitute(magnitude, k, next_index);
    let inverse = graph.node(core::POW, &[magnitude, minus_one]);
    let quotient = graph.node(core::MUL, &[next, inverse]);
    let ratio = cx.simplify(quotient);
    let graph = &mut *cx.graph;
    let reciprocal = graph.node(core::POW, &[k, minus_one]);
    let raw_root = graph.node(core::POW, &[magnitude, reciprocal]);
    let root = cx.simplify(raw_root);
    for test in [ratio, root] {
        if let Some(v) = settled_value(cx, test, k) {
            if (v - 1.0).abs() > 0.05 {
                return Some(v < 1.0);
            }
        }
    }
    let tests = if varying_power(cx.graph, term, k) { [root, ratio] } else { [ratio, root] };
    for test in tests {
        if let Some(v) = limit_value(cx, ops, test, k) {
            if (v - 1.0).abs() > 1e-9 {
                return Some(v < 1.0);
            }
        }
    }

    // 4. Alternating series test: the sign alternates and |a_k| decreases
    // to zero (checked on a long stretch of indices).
    if alternates_and_decreases(cx, term, k) {
        return Some(true);
    }

    // 5. Limit comparison with a p-series.
    let graph = &mut *cx.graph;
    let minus_one = graph.int(-1);
    let log_magnitude = graph.node(ln, &[magnitude]);
    let log_k = graph.node(ln, &[k]);
    let inverse_log_k = graph.node(core::POW, &[log_k, minus_one]);
    let neg = graph.node(core::MUL, &[minus_one, log_magnitude, inverse_log_k]);
    let exponent = limit_value(cx, ops, neg, k);
    if let Some(p) = exponent.filter(|p| p.is_finite()) {
        if p > 1.0 + 1e-9 {
            return Some(true);
        }
        if p < 1.0 - 1e-9 && same_sign {
            return Some(false);
        }
    }

    // p = 1 exactly: compare with 1/k directly, then condense.
    if same_sign {
        let weighted = cx.graph.node(core::MUL, &[k, magnitude]);
        if let Some(v) = limit_value(cx, ops, weighted, k) {
            if v.is_finite() && v > 1e-12 {
                return Some(false);
            }
        }
        if depth < 2 {
            let graph = &mut *cx.graph;
            let two = graph.int(2);
            let power = graph.node(core::POW, &[two, k]);
            let at_power = graph.substitute(magnitude, k, power);
            let condensed = graph.node(core::MUL, &[power, at_power]);
            return convergence(cx, ops, condensed, k, depth.saturating_add(1));
        }
    }
    None
}

/// Samples of `e` at integer `k` from `from` on.
fn samples(
    cx: &Cx<'_>,
    e: NodeId,
    k: NodeId,
    from: i32,
    count: i32,
) -> Option<Vec<f64>> {
    let symbol = cx.graph.symbol_of(k)?;
    (from..from.saturating_add(count))
        .map(|i| {
            let mut env = Env::numeric(0.0);
            env.bind(symbol, f64::from(i));
            cx.graph.eval(e, &env).filter(|v| v.is_finite())
        })
        .collect()
}

/// The value `e` settles to at large `k`, judged from samples at
/// `k = 10^3 .. 10^6`: `None` unless they agree to within 1%.
fn settled_value(
    cx: &Cx<'_>,
    e: NodeId,
    k: NodeId,
) -> Option<f64> {
    let symbol = cx.graph.symbol_of(k)?;
    let values: Option<Vec<f64>> = [1e3, 1e4, 1e5, 1e6]
        .iter()
        .map(|&at| {
            let mut env = Env::numeric(0.0);
            env.bind(symbol, at);
            cx.graph.eval(e, &env).filter(|v| v.is_finite())
        })
        .collect();
    let values = values?;
    let last = *values.last()?;
    let spread = values.iter().map(|v| (v - last).abs()).fold(0.0, f64::max);
    (spread <= 0.01 * last.abs().max(1e-3)).then_some(last)
}

/// Whether `e` contains a power whose base and exponent both depend on `k`.
fn varying_power(
    graph: &Graph,
    e: NodeId,
    k: NodeId,
) -> bool {
    let Some(symbol) = graph.symbol_of(k) else {
        return false;
    };
    let depends = |n: NodeId| graph.depends_on(graph.find(n), symbol);
    let mut stack = vec![e];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if let (true, &[base, exponent]) = (graph.op(n) == core::POW, graph.children(n)) {
            if depends(base) && depends(exponent) {
                return true;
            }
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// Whether `e` is positive at a large index.
fn sample_sign(
    cx: &Cx<'_>,
    e: NodeId,
    k: NodeId,
) -> Option<bool> {
    samples(cx, e, k, 100, 1).and_then(|v| v.first().map(|x| *x > 0.0))
}

fn eventually_one_sign(
    cx: &Cx<'_>,
    term: NodeId,
    k: NodeId,
) -> bool {
    samples(cx, term, k, 50, 200).is_some_and(|v| v.iter().all(|x| *x > 0.0) || v.iter().all(|x| *x < 0.0))
}

fn alternates_and_decreases(
    cx: &Cx<'_>,
    term: NodeId,
    k: NodeId,
) -> bool {
    let Some(values) = samples(cx, term, k, 20, 400) else {
        return false;
    };
    values.windows(2).all(|w| {
        let (a, b) = (w[0], w[1]);
        a * b < 0.0 && b.abs() <= a.abs()
    }) && values.last().is_some_and(|v| v.abs() < values[0].abs())
}

/// Symbolic kernel for every request of this module.
pub(super) struct SeriesKernel {
    pub(super) ops: SeriesOps,
}

impl Kernel for SeriesKernel {
    fn ops(&self) -> Vec<OpId> {
        let o = self.ops;
        vec![o.taylor, o.laurent, o.fourier, o.sum, o.product, o.converges]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let args = cx.graph.children(node).to_vec();
        let o = self.ops;
        if op == o.taylor {
            taylor(cx, o, &args).map_or(Outcome::Pass, Outcome::Pinned)
        } else if op == o.laurent {
            laurent(cx, o, &args).map_or(Outcome::Pass, Outcome::Pinned)
        } else if op == o.fourier {
            fourier(cx, o, &args).map_or(Outcome::Pass, Outcome::Equal)
        } else if op == o.sum {
            sum(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
        } else if op == o.product {
            product(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
        } else {
            converges(cx, o, &args).map_or(Outcome::Pass, Outcome::Equal)
        }
    }
}

/// Numeric kernel for `sum` and `product`.
pub(super) struct NumericSum {
    pub(super) sum: OpId,
    pub(super) product: OpId,
}

impl Kernel for NumericSum {
    fn ops(&self) -> Vec<OpId> {
        vec![self.sum, self.product]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        // Already evaluated in this run.
        if cx.graph.approx(cx.graph.find(node)).is_some() {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let heavy = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).flags.has(OpFlags::HEAVY);
        if graph.enodes(graph.find(node)).any(|n| !heavy(graph, n)) {
            return Outcome::Pass;
        }
        let &[f, k, lower, upper] = graph.children(node) else {
            return Outcome::Pass;
        };
        let is_sum = graph.op(node) == self.sum;
        let (Some(symbol), Some(term)) = (graph.symbol_of(k), best(graph, f)) else {
            return Outcome::Pass;
        };
        let bound = |graph: &mut Graph, n: NodeId| best(graph, n).and_then(|t| graph.eval(t, cx.env));
        let (Some(a), Some(b)) = (bound(graph, lower), bound(graph, upper)) else {
            return Outcome::Pass;
        };
        let mut inputs = vec![symbol];
        let mut values = vec![0.0];
        for &(s, v) in cx.env.bindings() {
            if s != symbol {
                inputs.push(s);
                values.push(v);
            }
        }
        let Ok(compiled) = Interpreter.compile(graph, term, &inputs) else {
            return Outcome::Pass;
        };
        let at = |index: f64| {
            let mut args = values.clone();
            if let Some(slot) = args.first_mut() {
                *slot = index;
            }
            compiled.call(&args)
        };
        if a.fract() != 0.0 || !a.is_finite() {
            return Outcome::Pass;
        }
        #[allow(clippy::cast_possible_truncation)]
        let from = a as i64;
        let tolerance = cx.env.tolerance.max(1e-14);
        if b == f64::INFINITY && is_sum {
            let result = sum_to_infinity(at, from, tolerance, 300_000);
            return if result.value.is_finite() && result.error.is_finite() {
                Outcome::Approx(Ball { mid: result.value, rad: result.error })
            } else {
                Outcome::Pass
            };
        }
        if !b.is_finite() || b.fract() != 0.0 || b - a > 5e7 {
            return Outcome::Pass;
        }
        #[allow(clippy::cast_possible_truncation)]
        let to = b as i64;
        let value = if is_sum {
            sum_range(at, from, to)
        } else {
            #[allow(clippy::cast_precision_loss)]
            (from..=to).fold(1.0, |acc, index| acc * at(index as f64))
        };
        if value.is_finite() {
            Outcome::Approx(Ball { mid: value, rad: f64::EPSILON * value.abs() * (b - a + 1.0) })
        } else {
            Outcome::Pass
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}
