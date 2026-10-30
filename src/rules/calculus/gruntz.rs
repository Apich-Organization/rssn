//! A Gruntz-style limit kernel for exp-log expressions.
//!
//! The heuristics in [`super::limits`] decide limits by continuity,
//! l'Hôpital's rule and leading terms of rational functions. They cannot see
//! the relative growth of nested exponentials and logarithms: the limit of
//! `(exp(exp(x)) + 1) / (exp(exp(x)) - 1)` or of `ln(ln(x)) / ln(x)` needs a
//! comparison of growth rates.
//!
//! This module implements the core of Gruntz's algorithm. The variable is
//! first moved to `t → +∞` (a finite point `a` becomes `a ± 1/t`). The
//! most rapidly varying exponential `E = exp(h)` of the expression is
//! identified, every other maximal exponential `exp(g)` is rewritten as
//! `E^(c) exp(g - c h)` with a rational `c`, and the expression is expanded
//! in the small variable `q = E` (when `E → 0`) or `q = 1/E` (when
//! `E → ∞`) with the series kernel of [`super::powerseries`]. The valuation
//! of that series decides the limit; a zero valuation recurses on the
//! constant coefficient, which involves only slower functions. Expressions
//! the algorithm does not understand (trigonometric functions, algebraic
//! cancellation at the chosen order) make it return `None`, and the caller
//! falls back on the existing heuristics.

#![allow(
    clippy::single_call_fn,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::rules::poly::best;

use super::laurent_expansion;
use super::limits::Side;

/// Deepest recursion of the growth comparison.
const MAX_DEPTH: usize = 8;
/// Order of the series in the mrv variable.
const SERIES_ORDER: i64 = 6;
/// Largest denominator accepted for a ratio of exponential arguments.
const MAX_RATIO_DENOMINATOR: i64 = 12;

/// The symbolic value `±∞` of an operator.
fn infinity(
    graph: &mut Graph,
    op: OpId,
    positive: bool,
) -> NodeId {
    let node = graph.node(op, &[]);
    if positive {
        node
    } else {
        let minus_one = graph.int(-1);
        graph.node(core::MUL, &[minus_one, node])
    }
}

/// `Some(true)` for `+∞`, `Some(false)` for `-∞`.
fn infinite_sign(
    graph: &Graph,
    op: OpId,
    node: NodeId,
) -> Option<bool> {
    if graph.op(node) == op {
        return Some(true);
    }
    if graph.op(node) == core::MUL {
        let children = graph.children(node);
        if children.len() == 2
            && children
                .iter()
                .any(|&c| graph.number_of(c).is_some_and(Number::is_negative))
            && children.iter().any(|&c| graph.op(c) == op)
        {
            return Some(false);
        }
    }
    None
}

fn depends(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> bool {
    graph.depends_on(graph.find(node), symbol)
}

fn mul_of(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn pow_of(
    graph: &mut Graph,
    base: NodeId,
    exponent: NodeId,
) -> NodeId {
    graph.node(core::POW, &[base, exponent])
}

fn sub_of(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    let negated = graph.node(core::MUL, &[minus_one, b]);
    graph.node(core::ADD, &[a, negated])
}

/// Bindings for the free symbols other than the limit variable, giving each
/// a fixed positive integer.
fn generic_env(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> Env {
    let mut env = Env::numeric(0.0);
    for &other in graph.free_symbols(graph.find(node)) {
        if other != symbol {
            env.bind(other, 2.0 + f64::from(other.raw() % 5));
        }
    }
    env
}

/// Whether the expression is built from the limit variable, numbers,
/// `+`, `*`, `^`, `exp` and `ln` only.
fn exp_log(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> bool {
    if !depends(graph, node, symbol) {
        return true;
    }
    let op = graph.op(node);
    if op == core::SYM {
        return true;
    }
    if op == core::ADD || op == core::MUL || op == core::POW {
        return graph
            .children(node)
            .iter()
            .all(|&child| exp_log(graph, child, symbol));
    }
    let name = graph.ops().get(op).name.as_ref();
    matches!(name, "exp" | "ln")
        && graph
            .children(node)
            .iter()
            .all(|&child| exp_log(graph, child, symbol))
}

/// The level of an exp-log expression in the tower: `t` is level 0, `exp`
/// raises it by one and `ln` lowers it by one.
fn rank(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> Option<i64> {
    if !depends(graph, node, symbol) {
        return None;
    }
    let op = graph.op(node);
    if op == core::SYM {
        return (graph.symbol_of(node) == Some(symbol)).then_some(0);
    }
    if op == core::ADD || op == core::MUL {
        return graph
            .children(node)
            .iter()
            .filter_map(|&child| rank(graph, child, symbol))
            .max();
    }
    if op == core::POW {
        let &[base, exponent] = graph.children(node) else {
            return None;
        };
        let base_rank = rank(graph, base, symbol);
        let exponent_rank = rank(graph, exponent, symbol);
        return match exponent_rank {
            | None => base_rank,
            | Some(high) => Some(1 + high.max(base_rank.unwrap_or(0).saturating_sub(1))),
        };
    }
    let name = graph.ops().get(op).name.as_ref();
    let &[argument] = graph.children(node) else {
        return None;
    };
    let inner = rank(graph, argument, symbol)?;
    match name {
        | "exp" => Some(inner.saturating_add(1)),
        | "ln" => Some(inner.saturating_sub(1)),
        | _ => None,
    }
}

/// Rewrites powers with a varying exponent as exponentials.
fn rewrite_powers(
    cx: &mut Cx<'_>,
    node: NodeId,
    symbol: SymbolId,
    order: u32,
) -> Option<NodeId> {
    if order > 64 {
        return None;
    }
    if !depends(cx.graph, node, symbol) {
        return Some(node);
    }
    let op = cx.graph.op(node);
    if op == core::POW {
        let &[base, exponent] = cx.graph.children(node) else {
            return None;
        };
        let base = rewrite_powers(cx, base, symbol, order + 1)?;
        let exponent = rewrite_powers(cx, exponent, symbol, order + 1)?;
        if !depends(cx.graph, exponent, symbol) {
            return cx.graph.try_node(core::POW, &[base, exponent]);
        }
        if !positive_at_infinity(cx, base, symbol) {
            return None;
        }
        let ln = cx.graph.ops().lookup("ln")?;
        let exp = cx.graph.ops().lookup("exp")?;
        let log = cx.graph.node(ln, &[base]);
        let product = cx.graph.node(core::MUL, &[exponent, log]);
        return Some(cx.graph.node(exp, &[product]));
    }
    let children = cx.graph.children(node).to_vec();
    let mut rebuilt = Vec::with_capacity(children.len());
    let mut changed = false;
    for child in children {
        let new = rewrite_powers(cx, child, symbol, order + 1)?;
        changed |= new != child;
        rebuilt.push(new);
    }
    if changed {
        cx.graph.try_node(op, &rebuilt)
    } else {
        Some(node)
    }
}

/// Whether `base` is positive for large values of the limit variable.
fn positive_at_infinity(
    cx: &Cx<'_>,
    base: NodeId,
    symbol: SymbolId,
) -> bool {
    let env = generic_env(cx.graph, base, symbol);
    let mut checked = 0;
    for exponent in 2..=5 {
        let point = 10.0_f64.powi(exponent);
        let mut env = env.clone();
        env.bind(symbol, point);
        let Some(value) = cx.graph.eval(base, &env) else {
            continue;
        };
        if !value.is_finite() || value <= 0.0 {
            return false;
        }
        checked += 1;
    }
    checked >= 2
}

/// The exponentials of maximal rank, with their arguments.
fn maximal_exponentials(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> Option<(i64, Vec<(NodeId, NodeId)>)> {
    let exp = graph.ops().lookup("exp")?;
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    let mut found: Vec<(i64, NodeId, NodeId)> = Vec::new();
    while let Some(current) = stack.pop() {
        if !seen.insert(current) {
            continue;
        }
        if graph.op(current) == exp
            && let Some(level) = rank(graph, current, symbol)
            && let &[argument] = graph.children(current)
        {
            found.push((level, current, argument));
        }
        stack.extend_from_slice(graph.children(current));
    }
    let highest = found.iter().map(|&(level, _, _)| level).max()?;
    if highest <= 0 {
        return None;
    }
    Some((
        highest,
        found
            .into_iter()
            .filter(|&(level, _, _)| level == highest)
            .map(|(_, e, arg)| (e, arg))
            .collect(),
    ))
}

/// The numeric value of `node` at a large point, with generic parameters.
fn sample(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
    exponent: i32,
) -> Option<f64> {
    let env = generic_env(graph, node, symbol);
    let mut env = env;
    env.bind(symbol, 10.0_f64.powi(exponent));
    graph.eval(node, &env)
}

/// A rational close to `value`, or `None`.
fn small_rational(value: f64) -> Option<(i64, i64)> {
    if !value.is_finite() || value.abs() > 1e6 {
        return None;
    }
    for denominator in 1..=MAX_RATIO_DENOMINATOR {
        #[allow(clippy::cast_precision_loss)]
        let numerator = (value * denominator as f64).round();
        #[allow(clippy::cast_possible_truncation)]
        let numerator = numerator as i64;
        #[allow(clippy::cast_precision_loss)]
        let exact = numerator as f64 / denominator as f64;
        if (exact - value).abs() <= 1e-6 * value.abs().max(1.0) {
            return Some((numerator, denominator));
        }
    }
    None
}

fn lcm(
    a: i64,
    b: i64,
) -> Option<i64> {
    if a == 0 || b == 0 {
        return Some(0);
    }
    let (mut x, mut y) = (a.abs(), b.abs());
    while y != 0 {
        let r = x % y;
        x = y;
        y = r;
    }
    a.checked_div(x)?.checked_mul(b)
}

/// The rational `c` with `g / h → c`, verified to leave a subdominant rest.
fn rational_ratio(
    cx: &mut Cx<'_>,
    g: NodeId,
    h: NodeId,
    symbol: SymbolId,
) -> Option<(i64, i64)> {
    let mut ratio = 0.0;
    let mut sampled = 0;
    for exponent in 3..=5 {
        let Some(gv) = sample(cx.graph, g, symbol, exponent) else {
            continue;
        };
        let Some(hv) = sample(cx.graph, h, symbol, exponent) else {
            continue;
        };
        if !gv.is_finite() || !hv.is_finite() || hv.abs() < 1e-12 {
            continue;
        }
        ratio = gv / hv;
        sampled += 1;
    }
    if sampled < 2 {
        return None;
    }
    let (numerator, denominator) = small_rational(ratio)?;
    // exp(g) = exp(c h) exp(r): r must be subdominant to h.
    let scale = cx.graph.num(Number::fraction(numerator, denominator)?);
    let product = cx.graph.node(core::MUL, &[scale, h]);
    let rest = sub_of(cx.graph, g, product);
    let rest = best(cx.graph, rest)?;
    if !depends(cx.graph, rest, symbol) {
        return Some((numerator, denominator));
    }
    for exponent in 3..=5 {
        let Some(rv) = sample(cx.graph, rest, symbol, exponent) else {
            continue;
        };
        let Some(hv) = sample(cx.graph, h, symbol, exponent) else {
            continue;
        };
        if hv.abs() < 1e-12 {
            continue;
        }
        let relative = (rv / hv).abs();
        if relative > 0.05 {
            return None;
        }
    }
    Some((numerator, denominator))
}

/// The limit of an exp-log expression as the variable tends to `+∞`.
fn infinity_limit(
    cx: &mut Cx<'_>,
    infinity_op: OpId,
    f: NodeId,
    t: NodeId,
    symbol: SymbolId,
    depth: usize,
) -> Option<NodeId> {
    if depth > MAX_DEPTH {
        return None;
    }
    let simplified = cx.simplify(f);
    if !depends(cx.graph, simplified, symbol) {
        return Some(simplified);
    }
    // The variable itself tends to +∞ by construction.
    if cx.graph.same(simplified, t) {
        return Some(infinity(cx.graph, infinity_op, true));
    }
    if !exp_log(cx.graph, simplified, symbol) {
        return None;
    }
    let rewritten = rewrite_powers(cx, simplified, symbol, 0)?;
    let rewritten = cx.simplify(rewritten);
    let rewritten = reduce_logs(cx, rewritten, symbol, 0)?;
    let rewritten = cx.simplify(rewritten);
    if !depends(cx.graph, rewritten, symbol) {
        return Some(rewritten);
    }
    let Some((_, candidates)) = maximal_exponentials(cx.graph, rewritten, symbol) else {
        // No growing exponential: t = exp(s) turns logarithms into powers.
        let fresh = cx.graph.interner_mut().fresh_symbol("s");
        cx.graph.assume(fresh, Facts::POSITIVE);
        let s = cx.graph.symbol_node(fresh);
        let exp = cx.graph.ops().lookup("exp")?;
        let lifted = cx.graph.node(exp, &[s]);
        let substituted = cx.graph.substitute(rewritten, t, lifted);
        let substituted = cx.simplify(substituted);
        return infinity_limit(cx, infinity_op, substituted, s, fresh, depth + 1);
    };
    let mut scored: Vec<(f64, NodeId, NodeId)> = Vec::new();
    for (node, argument) in candidates {
        let score = sample(cx.graph, argument, symbol, 4).unwrap_or(f64::NEG_INFINITY);
        scored.push((score, node, argument));
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (_, _, h) in &scored {
        let Some(result) = limit_around(cx, infinity_op, rewritten, t, symbol, *h, &scored, depth)
        else {
            continue;
        };
        return Some(result);
    }
    None
}

/// Rewrites every maximal exponential around `h` and takes the limit in the
/// resulting single exponential.
fn limit_around(
    cx: &mut Cx<'_>,
    infinity_op: OpId,
    f: NodeId,
    t: NodeId,
    symbol: SymbolId,
    h: NodeId,
    scored: &[(f64, NodeId, NodeId)],
    depth: usize,
) -> Option<NodeId> {
    let h_limit = infinity_limit(cx, infinity_op, h, t, symbol, depth + 1)?;
    if depends(cx.graph, h_limit, symbol) {
        return None;
    }
    let direction = infinite_sign(cx.graph, infinity_op, h_limit);
    let fresh = cx.graph.interner_mut().fresh_symbol("q");
    cx.graph.assume(fresh, Facts::POSITIVE);
    let q = cx.graph.symbol_node(fresh);
    let exp = cx.graph.ops().lookup("exp")?;
    let in_q = match direction {
        | Some(growing) => {
            let mut denominator = 1_i64;
            let mut parts: Vec<(NodeId, i64, NodeId)> = Vec::new();
            for &(_, node, argument) in scored {
                let (numerator, den) = rational_ratio(cx, argument, h, symbol)?;
                denominator = lcm(denominator, den)?;
                if denominator > 64 {
                    return None;
                }
                let exponent = numerator.checked_mul(denominator.checked_div(den)?)?;
                let scale = cx.graph.num(Number::fraction(numerator, den)?);
                let product = cx.graph.node(core::MUL, &[scale, h]);
                let rest = sub_of(cx.graph, argument, product);
                let rest = cx.simplify(rest);
                parts.push((node, exponent, rest));
            }
            let denominator_node = cx.graph.num(Number::fraction(1, denominator)?);
            let scaled_h = if denominator == 1 {
                h
            } else {
                cx.graph.node(core::MUL, &[denominator_node, h])
            };
            let prime = cx.graph.node(exp, &[scaled_h]);
            let mut rewritten = f;
            for (node, exponent, rest) in parts {
                let mut factors = Vec::new();
                if exponent != 0 {
                    let e = cx.graph.int(exponent);
                    factors.push(pow_of(cx.graph, prime, e));
                }
                if !cx.graph.number_of(rest).is_some_and(Number::is_zero) {
                    factors.push(cx.graph.node(exp, &[rest]));
                }
                let replacement = mul_of(cx.graph, &factors);
                rewritten = cx.graph.replace_subterm(rewritten, node, replacement);
            }
            let small = if growing {
                let minus_one = cx.graph.int(-1);
                pow_of(cx.graph, q, minus_one)
            } else {
                q
            };
            cx.graph.replace_subterm(rewritten, prime, small)
        },
        | None => {
            let mut rewritten = f;
            for &(_, node, argument) in scored {
                let (numerator, den) = rational_ratio(cx, argument, h, symbol)?;
                let scale = cx.graph.num(Number::fraction(numerator, den)?);
                let product = if numerator == den {
                    q
                } else {
                    cx.graph.node(core::MUL, &[scale, q])
                };
                let exponential = cx.graph.node(exp, &[product]);
                let scaled_h = cx.graph.node(core::MUL, &[scale, h]);
                let rest = sub_of(cx.graph, argument, scaled_h);
                let rest = cx.simplify(rest);
                let replacement = if cx.graph.number_of(rest).is_some_and(Number::is_zero) {
                    exponential
                } else {
                    let tail = cx.graph.node(exp, &[rest]);
                    mul_of(cx.graph, &[exponential, tail])
                };
                rewritten = cx.graph.replace_subterm(rewritten, node, replacement);
            }
            let rewritten = cx.graph.substitute(rewritten, h, q);
            cx.simplify(rewritten)
        },
    };
    let zero = cx.graph.int(0);
    let (valuation, coefficients) = laurent_expansion(cx, in_q, q, zero, SERIES_ORDER)?;
    if valuation > 0 {
        return Some(cx.graph.int(0));
    }
    if valuation < 0 {
        return divergence(cx, infinity_op, f, symbol);
    }
    let leading = *coefficients.first()?;
    let leading = cx.simplify(leading);
    if !depends(cx.graph, leading, symbol) {
        return Some(leading);
    }
    infinity_limit(cx, infinity_op, leading, t, symbol, depth + 1)
}

/// Removes `ln(exp(u))` (the inner exponential disappears from the growth
/// comparison): without this the series in the mrv variable would have to
/// expand `ln(q)`, which the power-series kernel does not do.
fn reduce_logs(
    cx: &mut Cx<'_>,
    node: NodeId,
    symbol: SymbolId,
    depth: usize,
) -> Option<NodeId> {
    if depth > 64 {
        return None;
    }
    if !depends(cx.graph, node, symbol) {
        return Some(node);
    }
    let op = cx.graph.op(node);
    let children = cx.graph.children(node).to_vec();
    let name = cx.graph.ops().get(op).name.to_string();
    if name == "ln"
        && let &[argument] = children.as_slice()
        && cx.graph.ops().get(cx.graph.op(argument)).name.as_ref() == "exp"
        && let &[inner] = cx.graph.children(argument)
    {
        return reduce_logs(cx, inner, symbol, depth + 1);
    }
    let mut rebuilt = Vec::with_capacity(children.len());
    let mut changed = false;
    for child in children {
        let new = reduce_logs(cx, child, symbol, depth + 1)?;
        changed |= new != child;
        rebuilt.push(new);
    }
    if changed {
        cx.graph.try_node(op, &rebuilt)
    } else {
        Some(node)
    }
}

/// `±∞` from the sign of `f` far out.
fn divergence(
    cx: &mut Cx<'_>,
    infinity_op: OpId,
    f: NodeId,
    symbol: SymbolId,
) -> Option<NodeId> {
    let mut positive = 0;
    let mut negative = 0;
    for exponent in 3..=6 {
        let Some(value) = sample(cx.graph, f, symbol, exponent) else {
            continue;
        };
        if value.is_nan() {
            continue;
        }
        if value > 0.0 {
            positive += 1;
        } else if value < 0.0 {
            negative += 1;
        }
    }
    if positive >= 2 && negative == 0 {
        Some(infinity(cx.graph, infinity_op, true))
    } else if negative >= 2 && positive == 0 {
        Some(infinity(cx.graph, infinity_op, false))
    } else {
        None
    }
}

/// Whether the expression contains an exp-log atom or a power with a
/// varying exponent — the cases this kernel can decide.
fn interesting(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> bool {
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    while let Some(current) = stack.pop() {
        if !seen.insert(current) {
            continue;
        }
        let op = graph.op(current);
        if graph.depends_on(graph.find(current), symbol) {
            let name = graph.ops().get(op).name.as_ref();
            if name == "exp" || name == "ln" {
                return true;
            }
            if op == core::POW
                && let &[base, exponent] = graph.children(current)
                && graph.depends_on(graph.find(base), symbol)
                && graph.depends_on(graph.find(exponent), symbol)
            {
                return true;
            }
        }
        stack.extend_from_slice(graph.children(current));
    }
    false
}

/// The limit of `f` as `x → a` from `side`, by the Gruntz core.
pub(super) fn limit(
    cx: &mut Cx<'_>,
    infinity_op: OpId,
    f: NodeId,
    x: NodeId,
    point: NodeId,
    target: f64,
    side: Side,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(x)?;
    if !interesting(cx.graph, f, symbol) {
        return None;
    }
    let fresh = cx.graph.interner_mut().fresh_symbol("t");
    cx.graph.assume(fresh, Facts::POSITIVE);
    let t = cx.graph.symbol_node(fresh);
    let term = cx.simplify(f);
    let inverse = {
        let minus_one = cx.graph.int(-1);
        cx.graph.node(core::POW, &[t, minus_one])
    };
    let run = |cx: &mut Cx<'_>, replacement: NodeId| -> Option<NodeId> {
        let substituted = cx.graph.substitute(term, x, replacement);
        let substituted = cx.simplify(substituted);
        infinity_limit(cx, infinity_op, substituted, t, fresh, 0)
    };
    if target == f64::INFINITY {
        return run(&mut *cx, t);
    }
    if target == f64::NEG_INFINITY {
        let minus_one = cx.graph.int(-1);
        let negated = cx.graph.node(core::MUL, &[minus_one, t]);
        return run(&mut *cx, negated);
    }
    let point = best(cx.graph, point)?;
    let above = cx.graph.node(core::ADD, &[point, inverse]);
    let below = sub_of(cx.graph, point, inverse);
    match side {
        | Side::Above => run(&mut *cx, above),
        | Side::Below => run(&mut *cx, below),
        | Side::Both => {
            let above = run(&mut *cx, above);
            let below = run(&mut *cx, below);
            let (Some(above), Some(below)) = (above, below) else {
                return None;
            };
            if cx.graph.same(above, below) {
                return Some(above);
            }
            let sign = infinite_sign(cx.graph, infinity_op, above);
            if sign.is_some() && sign == infinite_sign(cx.graph, infinity_op, below) {
                return Some(above);
            }
            let env = generic_env(cx.graph, above, fresh);
            let (Some(a), Some(b)) = (cx.graph.eval(above, &env), cx.graph.eval(below, &env))
            else {
                return None;
            };
            if a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-8 * (1.0 + a.abs()) {
                Some(above)
            } else {
                None
            }
        },
    }
}
