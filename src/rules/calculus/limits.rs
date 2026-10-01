//! Limits.
//!
//! `limit(f, x, a)` is the two-sided limit of `f` as `x` tends to `a`
//! (`a` may be `oo` or `-oo`); a fourth argument, the symbol `plus` or
//! `minus`, asks for a one-sided limit.
//!
//! The symbolic kernel substitutes where the function is continuous and
//! otherwise resolves the indeterminate form: leading terms for rational
//! functions at infinity, l'Hôpital's rule for `0/0` and `∞/∞`, products
//! `0·∞` rewritten as quotients, and powers through `exp(g·ln f)`. Every
//! candidate is compared with a numeric probe of the function near the
//! point before it is accepted, so a misapplied rule cannot produce a
//! wrong limit. The numeric kernel returns that probe, extrapolated, when
//! a number is what was asked for.

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
use crate::graph::SymbolId;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;

use super::diff::Differentiate;
use super::integrate::Functions;

/// From which side the point is approached.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(super) enum Side {
    Both,
    Above,
    Below,
}

/// What a numeric probe of a function near a point suggests.
#[derive(Copy, Clone, Debug, PartialEq)]
pub(super) enum Probe {
    Finite(f64),
    PlusInfinity,
    MinusInfinity,
    Unknown,
}

/// Bindings giving every free symbol of `term` other than `x` a fixed
/// generic value, on top of `base`. The value depends only on the symbol,
/// so different terms of one problem see the same point.
fn generic_env(
    graph: &Graph,
    term: NodeId,
    x: SymbolId,
    base: &Env,
) -> Env {
    let mut env = base.clone();
    for &symbol in graph.free_symbols(graph.find(term)) {
        if symbol != x && env.value(symbol).is_none() {
            let id = symbol.raw();
            env.bind(symbol, 0.7 + 0.13 * f64::from(id % 7) + 0.017 * f64::from(id % 13));
        }
    }
    env
}

/// Evaluates `term` along a sequence approaching `target` from one side.
fn probe_one_side(
    graph: &Graph,
    term: NodeId,
    x: SymbolId,
    target: f64,
    from_above: bool,
    env: &Env,
) -> Probe {
    let mut env = env.clone();
    let mut values: Vec<f64> = Vec::new();
    for k in 6..=26 {
        let step = 0.5_f64.powi(k);
        let point = if target.is_infinite() {
            // Large arguments of alternating sign size: 2^k.
            target.signum() * 2.0_f64.powi(k)
        } else if from_above {
            target + step * target.abs().max(1.0)
        } else {
            target - step * target.abs().max(1.0)
        };
        env.bind(x, point);
        match graph.eval(term, &env) {
            | Some(v) if !v.is_nan() => values.push(v),
            | _ => {},
        }
    }
    let n = values.len();
    if n < 6 {
        return Probe::Unknown;
    }
    let tail = &values[n - 6..];
    let last = tail[5];
    if tail.iter().all(|v| v.is_infinite()) || (last.abs() > 1e6 && tail.windows(2).all(|w| w[1].abs() >= w[0].abs())) {
        return if last > 0.0 { Probe::PlusInfinity } else { Probe::MinusInfinity };
    }
    if !last.is_finite() {
        return Probe::Unknown;
    }
    // Slow divergence (logarithms): the values move steadily in one
    // direction and the steps do not shrink.
    let steps: Vec<f64> = values.windows(2).map(|w| w[1] - w[0]).collect();
    let recent = &steps[steps.len().saturating_sub(10)..];
    let one_way = recent.iter().all(|d| *d > 0.0) || recent.iter().all(|d| *d < 0.0);
    if one_way && recent.len() == 10 && recent[9].abs() >= 0.5 * recent[0].abs() && last.abs() > 5.0 {
        return if recent[9] > 0.0 { Probe::PlusInfinity } else { Probe::MinusInfinity };
    }
    // Converging: successive differences must shrink.
    let d1 = (tail[5] - tail[4]).abs();
    let d2 = (tail[4] - tail[3]).abs();
    let d3 = (tail[3] - tail[2]).abs();
    let settled = d1 <= 1e-9 * (1.0 + last.abs());
    if settled || (d1 <= d2 && d2 <= d3 * 1.0001 && d1 <= 1e-3 * (1.0 + last.abs())) {
        // One Richardson step on the last three values (first-order error).
        let extrapolated = if settled { last } else { 2.0 * tail[5] - tail[4] };
        Probe::Finite(extrapolated)
    } else {
        Probe::Unknown
    }
}

/// Numeric evidence for the limit of `term` as `x` tends to `target`.
pub(super) fn probe(
    graph: &Graph,
    term: NodeId,
    x: SymbolId,
    target: f64,
    side: Side,
    env: &Env,
) -> Probe {
    if target.is_infinite() {
        return probe_one_side(graph, term, x, target, true, env);
    }
    let above = probe_one_side(graph, term, x, target, true, env);
    let below = probe_one_side(graph, term, x, target, false, env);
    match side {
        | Side::Above => above,
        | Side::Below => below,
        | Side::Both => match (above, below) {
            | (Probe::Finite(a), Probe::Finite(b)) if (a - b).abs() <= 1e-5 * (1.0 + a.abs()) => {
                Probe::Finite(f64::midpoint(a, b))
            },
            | (a, b) if a == b => a,
            // Defined on one side only (sqrt, ln): that side decides.
            | (a, Probe::Unknown) => a,
            | (Probe::Unknown, b) => b,
            | _ => Probe::Unknown,
        },
    }
}

struct Limiter<'c, 'a> {
    cx: &'c mut Cx<'a>,
    f: Functions,
    infinity: OpId,
    x: NodeId,
    symbol: SymbolId,
    point: NodeId,
    target: f64,
    side: Side,
}

const MAX_DEPTH: usize = 6;

impl Limiter<'_, '_> {
    fn depends(
        &self,
        node: NodeId,
    ) -> bool {
        self.cx.graph.depends_on(self.cx.graph.find(node), self.symbol)
    }

    fn probe(
        &self,
        term: NodeId,
    ) -> Probe {
        let env = generic_env(self.cx.graph, term, self.symbol, &Env::numeric(0.0));
        probe(self.cx.graph, term, self.symbol, self.target, self.side, &env)
    }

    fn infinite(
        &mut self,
        positive: bool,
    ) -> NodeId {
        let oo = self.cx.graph.node(self.infinity, &[]);
        if positive {
            oo
        } else {
            let minus_one = self.cx.graph.int(-1);
            self.cx.graph.node(core::MUL, &[minus_one, oo])
        }
    }

    fn derivative(
        &mut self,
        term: NodeId,
    ) -> NodeId {
        let raw = Differentiate { diff: self.f.diff }.derive(self.cx.graph, term, self.symbol, self.x);
        self.cx.simplify(raw)
    }

    /// The limit of `f`, unverified.
    fn limit(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        if depth > MAX_DEPTH {
            return None;
        }
        if !self.depends(f) {
            return Some(f);
        }
        // Where the function is continuous, substitute. Whether it is
        // defined at the point is decided on the raw substituted term:
        // simplifying first would turn 0/0 into 0.
        if self.target.is_finite() {
            let raw = self.cx.graph.substitute(f, self.x, self.point);
            let env = generic_env(self.cx.graph, raw, self.symbol, &Env::numeric(0.0));
            if self.cx.graph.eval(raw, &env).is_some_and(f64::is_finite) {
                return Some(self.cx.simplify(raw));
            }
            // Not evaluable for lack of semantics (an undetermined
            // function), rather than undefined: nothing to be done.
            self.cx.graph.eval(raw, &env)?;
        }
        match self.probe(f) {
            | Probe::PlusInfinity if self.simple_divergence(f) => return Some(self.infinite(true)),
            | Probe::MinusInfinity if self.simple_divergence(f) => return Some(self.infinite(false)),
            | _ => {},
        }
        if let Some(found) = self.rational_at_infinity(f) {
            return Some(found);
        }
        let op = self.cx.graph.op(f);
        if op == core::ADD || op == core::MUL {
            // Sum or product of limits when every part has a finite one.
            let parts = self.cx.graph.children(f).to_vec();
            let mut limits = Vec::with_capacity(parts.len());
            for part in parts {
                match self.limit(part, depth + 1) {
                    | Some(l) if !self.is_infinite(l) => limits.push(l),
                    | _ => break,
                }
            }
            if limits.len() == self.cx.graph.children(f).len() {
                let combined = self.cx.graph.node(op, &limits);
                return Some(self.cx.simplify(combined));
            }
        }
        if let Some(found) = self.composition(f, depth) {
            return Some(found);
        }
        if let Some(found) = self.power_form(f, depth) {
            return Some(found);
        }
        self.hospital(f, depth)
    }

    fn is_infinite(
        &self,
        term: NodeId,
    ) -> bool {
        self.cx.graph.eval(term, &Env::numeric(0.0)).is_some_and(f64::is_infinite)
    }

    /// Whether `f` is a single elementary building block (a power, a
    /// function application) whose probe can be trusted to mean divergence
    /// rather than a large finite value.
    fn simple_divergence(
        &self,
        f: NodeId,
    ) -> bool {
        let op = self.cx.graph.op(f);
        op != core::ADD
    }

    /// Ratio of polynomials in `x` as `x` tends to infinity: compare
    /// degrees.
    fn rational_at_infinity(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        if self.target.is_finite() {
            return None;
        }
        let mut gens = Gens::default();
        let x = self.x;
        let gx = gens.index(self.cx.graph, x);
        let fraction = ratio(self.cx.graph, &mut gens, f, Limits::default())?;
        for g in fraction.numer.support().into_iter().chain(fraction.denom.support()) {
            if g != gx && gens.node(g).is_some_and(|n| self.depends(n)) {
                return None;
            }
        }
        let (dn, dd) = (fraction.numer.degree_in(gx), fraction.denom.degree_in(gx));
        let lead_n = fraction.numer.coefficients_in(gx).get(dn as usize)?.clone();
        let lead_d = fraction.denom.coefficients_in(gx).get(dd as usize)?.clone();
        if dn < dd {
            return Some(self.cx.graph.int(0));
        }
        let top = to_term(self.cx.graph, &gens, &lead_n);
        let bottom = to_term(self.cx.graph, &gens, &lead_d);
        let minus_one = self.cx.graph.int(-1);
        let inverse = self.cx.graph.node(core::POW, &[bottom, minus_one]);
        let quotient = self.cx.graph.node(core::MUL, &[top, inverse]);
        let quotient = self.cx.simplify(quotient);
        if dn == dd {
            return Some(quotient);
        }
        // Diverges: the sign needs a number.
        let value = self.cx.graph.eval(quotient, &Env::numeric(0.0))?;
        let odd = (dn - dd) % 2 == 1;
        let positive = (value > 0.0) == (self.target > 0.0 || !odd);
        Some(self.infinite(positive))
    }

    /// A function of something that has a limit: continuity where the
    /// inner limit is finite, the known behaviour at infinity otherwise.
    fn composition(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let op = self.cx.graph.op(f);
        let children = self.cx.graph.children(f).to_vec();
        let (inner_term, exponent) = match children.as_slice() {
            | &[u] if op != core::ADD && op != core::MUL => (u, None),
            | &[base, exp] if op == core::POW && !self.depends(exp) => (base, Some(exp)),
            | _ => return None,
        };
        let inner = self.limit(inner_term, depth + 1)?;
        if !self.is_infinite(inner) {
            let rebuilt = match exponent {
                | Some(exp) => self.cx.graph.node(core::POW, &[inner, exp]),
                | None => self.cx.graph.node(op, &[inner]),
            };
            let env = generic_env(self.cx.graph, rebuilt, self.symbol, &Env::numeric(0.0));
            return match self.cx.graph.eval(rebuilt, &env) {
                | Some(v) if v.is_finite() => Some(self.cx.simplify(rebuilt)),
                // The outer function has a pole or a logarithmic
                // singularity exactly there.
                | Some(v) if v.is_infinite() => Some(self.infinite(v > 0.0)),
                | _ => None,
            };
        }
        let up = self.cx.graph.eval(inner, &Env::numeric(0.0))? > 0.0;
        if let Some(exp) = exponent {
            let negative = self.cx.graph.facts(exp).has(crate::graph::Facts::NEGATIVE);
            let positive = self.cx.graph.facts(exp).has(crate::graph::Facts::POSITIVE);
            return match (negative, positive, up) {
                | (true, _, _) => Some(self.cx.graph.int(0)),
                | (_, true, true) => Some(self.infinite(true)),
                | _ => None,
            };
        }
        let graph = &mut *self.cx.graph;
        let pi = graph.ops().lookup("pi")?;
        let name = graph.ops().get(op).name.to_string();
        let sign = graph.int(if up { 1 } else { -1 });
        match name.as_str() {
            | "atan" => {
                let pi = graph.node(pi, &[]);
                let half = graph.num(Number::fraction(1, 2)?);
                Some(graph.node(core::MUL, &[sign, half, pi]))
            },
            | "tanh" => Some(sign),
            | "acot" | "sech" | "csch" => Some(graph.int(0)),
            | "coth" => Some(sign),
            | "exp" if up => Some(self.infinite(true)),
            | "exp" => Some(graph.int(0)),
            | "ln" | "sqrt" | "cosh" | "abs" | "asinh" | "acosh" if up || name == "cosh" || name == "abs" => {
                Some(self.infinite(true))
            },
            | "sinh" | "asinh" => Some(self.infinite(up)),
            | _ => {
                let known = graph.ops().attr::<super::AtInfinity>(op)?;
                let pattern = if up { known.plus.clone() } else { known.minus.clone() }?;
                pattern.instantiate(graph, &[])
            },
        }
    }

    /// `b^e` with a varying exponent: `exp(lim e·ln b)`.
    fn power_form(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let &[base, exponent] = self.cx.graph.children(f) else {
            return None;
        };
        if self.cx.graph.op(f) != core::POW || !self.depends(exponent) {
            return None;
        }
        let log = self.cx.graph.node(self.f.ln, &[base]);
        let product = self.cx.graph.node(core::MUL, &[exponent, log]);
        let inner = self.limit(product, depth + 1)?;
        if self.is_infinite(inner) {
            let positive = self.cx.graph.eval(inner, &Env::numeric(0.0))? > 0.0;
            return Some(if positive { self.infinite(true) } else { self.cx.graph.int(0) });
        }
        let result = self.cx.graph.node(self.f.exp, &[inner]);
        Some(self.cx.simplify(result))
    }

    /// Splits `f` into numerator and denominator for l'Hôpital's rule.
    fn quotient(
        &mut self,
        f: NodeId,
    ) -> Option<(NodeId, NodeId)> {
        let graph = &mut *self.cx.graph;
        let factors: Vec<NodeId> = if graph.op(f) == core::MUL { graph.children(f).to_vec() } else { vec![f] };
        let mut numer = Vec::new();
        let mut denom = Vec::new();
        for factor in &factors {
            let reciprocal = match *graph.children(*factor) {
                | [base, exp] if graph.op(*factor) == core::POW => {
                    graph.number_of(exp).filter(|n| n.is_negative()).map(|n| (base, n.neg()))
                },
                | _ => None,
            };
            match reciprocal {
                | Some((base, positive)) => {
                    let e = graph.num(positive);
                    denom.push(if graph.number_of(e).is_some_and(Number::is_one) {
                        base
                    } else {
                        graph.node(core::POW, &[base, e])
                    });
                },
                | None => numer.push(*factor),
            }
        }
        let build = |graph: &mut Graph, parts: &[NodeId]| match parts {
            | [] => graph.int(1),
            | [only] => *only,
            | _ => graph.node(core::MUL, parts),
        };
        if !denom.is_empty() {
            return Some((build(graph, &numer), build(graph, &denom)));
        }
        if numer.len() < 2 {
            return None;
        }
        // A product 0·∞ or ∞·0: the simplest factor goes underneath as its
        // reciprocal (x·ln(1 + 1/x) becomes ln(1 + 1/x) / (1/x)), which
        // keeps the derivatives manageable.
        let size = |graph: &Graph, node: NodeId| {
            let mut count = 0_usize;
            let mut stack = vec![node];
            while let Some(n) = stack.pop() {
                count += 1;
                stack.extend_from_slice(graph.children(n));
            }
            count
        };
        // An exponential that vanishes goes underneath instead
        // (x·exp(-x²) becomes x / exp(x²)): differentiating the exponential
        // never makes it simpler, so it must end up where it only grows.
        let exp = self.f.exp;
        let decaying = numer
            .iter()
            .position(|&n| self.cx.graph.op(n) == exp && matches!(self.probe(n), Probe::Finite(v) if v.abs() < 1e-6));
        let graph = &mut *self.cx.graph;
        let (index, _) = match decaying {
            | Some(k) => (k, 0),
            | None => numer.iter().enumerate().map(|(k, &n)| (k, size(graph, n))).min_by_key(|&(_, s)| s)?,
        };
        let simplest = *numer.get(index)?;
        let rest: Vec<NodeId> = numer.iter().enumerate().filter(|&(i, _)| i != index).map(|(_, &n)| n).collect();
        let minus_one = graph.int(-1);
        let reciprocal = graph.node(core::POW, &[simplest, minus_one]);
        let top = build(graph, &rest);
        let reciprocal = self.cx.simplify(reciprocal);
        Some((top, reciprocal))
    }

    /// L'Hôpital's rule for `0/0` and `∞/∞`.
    fn hospital(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let (numer, denom) = self.quotient(f)?;
        let (pn, pd) = (self.probe(numer), self.probe(denom));
        let vanishing = |p: Probe| matches!(p, Probe::Finite(v) if v.abs() < 1e-6);
        let diverging = |p: Probe| matches!(p, Probe::PlusInfinity | Probe::MinusInfinity);
        if !((vanishing(pn) && vanishing(pd)) || (diverging(pn) && diverging(pd))) {
            // Not indeterminate: finite over non-zero, or a genuine pole.
            return match (pn, pd) {
                | (Probe::Finite(_), Probe::Finite(d)) if d.abs() > 1e-6 => {
                    let top = self.limit(numer, depth + 1)?;
                    let bottom = self.limit(denom, depth + 1)?;
                    let minus_one = self.cx.graph.int(-1);
                    let inverse = self.cx.graph.node(core::POW, &[bottom, minus_one]);
                    let quotient = self.cx.graph.node(core::MUL, &[top, inverse]);
                    Some(self.cx.simplify(quotient))
                },
                | (Probe::Finite(_), p) if diverging(p) => Some(self.cx.graph.int(0)),
                | _ => match self.probe(f) {
                    | Probe::PlusInfinity => Some(self.infinite(true)),
                    | Probe::MinusInfinity => Some(self.infinite(false)),
                    | _ => None,
                },
            };
        }
        let (dn, dd) = (self.derivative(numer), self.derivative(denom));
        let minus_one = self.cx.graph.int(-1);
        let inverse = self.cx.graph.node(core::POW, &[dd, minus_one]);
        let next = self.cx.graph.node(core::MUL, &[dn, inverse]);
        let next = self.cx.simplify(next);
        if next == f {
            // The rule reproduced the problem (abs(x)/x): no progress.
            return self.recognised(f);
        }
        self.limit(next, depth + 1).or_else(|| self.recognised(f))
    }

    /// Last resort for a parameter-free function: if the numeric probe
    /// settles on a small rational, that is the limit.
    fn recognised(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let graph = &*self.cx.graph;
        if graph.free_symbols(graph.find(f)).iter().any(|&s| s != self.symbol) {
            return None;
        }
        let Probe::Finite(value) = self.probe(f) else {
            return None;
        };
        (1..=12_i64).find_map(|denominator| {
            #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
            let numerator = (value * denominator as f64).round() as i64;
            #[allow(clippy::cast_precision_loss)]
            let exact = numerator as f64 / denominator as f64;
            ((exact - value).abs() <= 1e-9 * (1.0 + value.abs()))
                .then(|| Number::fraction(numerator, denominator))
                .flatten()
                .map(|n| self.cx.graph.num(n))
        })
    }

    /// Whether `candidate` agrees with the numeric behaviour of `f`.
    fn verified(
        &self,
        f: NodeId,
        candidate: NodeId,
    ) -> bool {
        let env = generic_env(self.cx.graph, f, self.symbol, &Env::numeric(0.0));
        let env = generic_env(self.cx.graph, candidate, self.symbol, &env);
        let observed = probe(self.cx.graph, f, self.symbol, self.target, self.side, &env);
        let Some(claimed) = self.cx.graph.eval(candidate, &env) else {
            // Nothing to compare with: no evidence against.
            return true;
        };
        match observed {
            | Probe::Finite(v) => claimed.is_finite() && (claimed - v).abs() <= 1e-4 * (1.0 + v.abs()),
            | Probe::PlusInfinity => claimed == f64::INFINITY,
            | Probe::MinusInfinity => claimed == f64::NEG_INFINITY,
            | Probe::Unknown => claimed.is_finite(),
        }
    }
}

/// Reads the arguments of a `limit` node.
fn arguments(
    graph: &Graph,
    node: NodeId,
) -> Option<(NodeId, NodeId, NodeId, Side)> {
    match *graph.children(node) {
        | [f, x, a] => Some((f, x, a, Side::Both)),
        | [f, x, a, direction] => {
            let name = graph.symbol_of(direction).map(|s| graph.interner().symbol_name(s).to_owned())?;
            let side = match name.as_str() {
                | "plus" => Side::Above,
                | "minus" => Side::Below,
                | _ => return None,
            };
            Some((f, x, a, side))
        },
        | _ => None,
    }
}

/// Symbolic kernel for `limit`.
pub(super) struct SymbolicLimit {
    pub(super) limit: OpId,
    pub(super) infinity: OpId,
    pub(super) functions: Functions,
}

impl Kernel for SymbolicLimit {
    fn ops(&self) -> Vec<OpId> {
        vec![self.limit]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let Some((f, x, a, side)) = arguments(cx.graph, node) else {
            return Outcome::Pass;
        };
        limit_at(cx, self.functions, self.infinity, f, x, a, side).map_or(Outcome::Pass, Outcome::Equal)
    }
}

/// The verified limit of `f` as `x` tends to `a`, if one can be found.
pub(super) fn limit_at(
    cx: &mut Cx<'_>,
    functions: Functions,
    infinity: OpId,
    f: NodeId,
    x: NodeId,
    a: NodeId,
    side: Side,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(x)?;
    let point = best(cx.graph, a)?;
    // A symbolic point is given a generic value: enough to decide
    // continuity there, which is all that can be done symbolically.
    let env = generic_env(cx.graph, point, symbol, &Env::numeric(0.0));
    let target = cx.graph.eval(point, &env)?;
    let term = cx.simplify(f);
    let mut limiter = Limiter { cx, f: functions, infinity, x, symbol, point, target, side };
    let candidate = limiter.limit(term, 0)?;
    limiter.verified(term, candidate).then_some(candidate)
}

/// Numeric kernel for `limit`: the extrapolated probe.
pub(super) struct NumericLimit {
    pub(super) limit: OpId,
}

impl Kernel for NumericLimit {
    fn ops(&self) -> Vec<OpId> {
        vec![self.limit]
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
        let Some((f, x, a, side)) = arguments(graph, node) else {
            return Outcome::Pass;
        };
        let (Some(symbol), Some(term), Some(point)) = (graph.symbol_of(x), best(graph, f), best(graph, a)) else {
            return Outcome::Pass;
        };
        let Some(target) = graph.eval(point, cx.env) else {
            return Outcome::Pass;
        };
        match probe(graph, term, symbol, target, side, cx.env) {
            | Probe::Finite(value) => Outcome::Approx(Ball { mid: value, rad: 1e-6 * (1.0 + value.abs()) }),
            | _ => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}
