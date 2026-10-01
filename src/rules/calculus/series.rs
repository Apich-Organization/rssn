//! Series, sums and products.
//!
//! * `taylor(f, x, a, n)` and `laurent(f, x, a, n)`: expansions about `a`
//!   up to `(x - a)^n`, returned in that form (pinned).
//! * `fourier_series(f, x, L, n)`: the first `n` harmonics of `f` on
//!   `[-L, L]`.
//! * `sum(f, k, a, b)` and `product(f, k, a, b)`: finite sums are written
//!   out; polynomial, geometric and a few classical infinite summands have
//!   closed forms; anything else is summed numerically, with convergence
//!   acceleration, when a number is wanted.
//! * `converges(f, k)`: the ratio test for the series with general term
//!   `f`.

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
    None
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

/// Ratio test: the limit of `|a_{k+1} / a_k|` decides convergence unless
/// it equals one.
fn converges(
    cx: &mut Cx<'_>,
    ops: SeriesOps,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, k] = args else {
        return None;
    };
    let term = cx.simplify(f);
    let graph = &mut *cx.graph;
    let abs = graph.ops().lookup("abs")?;
    let one = graph.int(1);
    let next_index = graph.node(core::ADD, &[k, one]);
    let next = graph.substitute(term, k, next_index);
    let minus_one = graph.int(-1);
    let inverse = graph.node(core::POW, &[term, minus_one]);
    let quotient = graph.node(core::MUL, &[next, inverse]);
    let ratio = graph.node(abs, &[quotient]);
    let infinity = graph.node(ops.infinity, &[]);
    let limit = limit_at(cx, ops.functions, ops.infinity, ratio, k, infinity, Side::Both)?;
    let value = cx.graph.eval(limit, &Env::numeric(0.0))?;
    if (value - 1.0).abs() < 1e-9 {
        // Inconclusive: compare with a p-series through k^2 * a_k.
        let two = cx.graph.int(2);
        let square = cx.graph.node(core::POW, &[k, two]);
        let scaled = cx.graph.node(core::MUL, &[square, term]);
        let bounded = limit_at(cx, ops.functions, ops.infinity, scaled, k, infinity, Side::Both)
            .and_then(|l| cx.graph.eval(l, &Env::numeric(0.0)))
            .is_some_and(f64::is_finite);
        if bounded {
            return Some(cx.graph.lit(Payload::Bool(true)));
        }
        let weighted = cx.graph.node(core::MUL, &[k, term]);
        let harmonic_like = limit_at(cx, ops.functions, ops.infinity, weighted, k, infinity, Side::Both)
            .and_then(|l| cx.graph.eval(l, &Env::numeric(0.0)))
            .is_some_and(|v| v != 0.0);
        return harmonic_like.then(|| cx.graph.lit(Payload::Bool(false)));
    }
    Some(cx.graph.lit(Payload::Bool(value < 1.0)))
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
