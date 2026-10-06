//! Rewrites that put an equation in a form where the polynomial solver
//! sees its common generators.
//!
//! * [`exponentials`] makes every exponential in one unknown a power of a
//!   single generator: `exp(x + 1)` becomes `E exp(x)`, `exp(-x)` becomes
//!   `exp(x)^-1`, `exp(3 x/2)` and `exp(x)` become `exp(x/2)^3` and
//!   `exp(x/2)^2`, `4^x` becomes `(2^x)^2`, `(1/2)^x` becomes `(2^x)^-1`.
//! * [`hyperbolic_to_exponentials`] writes the hyperbolic functions through
//!   `exp`.
//! * [`powers_to_exp`] writes `b^u` as `exp(u ln b)`.

use std::collections::HashMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;

use super::product;
use crate::graph::op::core;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::SymbolId;

/// Rebuilds `term` bottom-up; `f` may replace a node given its rebuilt
/// children.
pub(super) fn map_term(
    graph: &mut Graph,
    term: NodeId,
    f: &mut dyn FnMut(&mut Graph, NodeId, &[NodeId]) -> Option<NodeId>,
) -> NodeId {
    let mut done: HashMap<NodeId, NodeId> = HashMap::new();
    let mut stack = vec![(term, false)];
    while let Some((node, expanded)) = stack.pop() {
        if done.contains_key(&node) {
            continue;
        }
        let children = graph.children(node).to_vec();
        if children.is_empty() {
            done.insert(node, node);
            continue;
        }
        if !expanded {
            stack.push((node, true));
            stack.extend(children.iter().filter(|c| !done.contains_key(c)).map(|&c| (c, false)));
            continue;
        }
        let rebuilt: Vec<NodeId> = children.iter().map(|c| done.get(c).copied().unwrap_or(*c)).collect();
        let new = match f(graph, node, &rebuilt) {
            | Some(replacement) => replacement,
            | None if rebuilt == children => node,
            | None => graph.try_node(graph.op(node), &rebuilt).unwrap_or(node),
        };
        done.insert(node, new);
    }
    done.get(&term).copied().unwrap_or(term)
}

/// What `node` is as an exponential of the unknown: `(op, base, argument)`
/// for `exp(argument)` (base `NONE`) and `base^argument` with a constant
/// base and an exponent that involves the unknown.
fn exponential_parts(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
    exp: Option<OpId>,
) -> Option<(OpId, NodeId, NodeId)> {
    let op = graph.op(node);
    let children = graph.children(node);
    if Some(op) == exp {
        let &[argument] = children else {
            return None;
        };
        return graph.depends_on(graph.find(argument), symbol).then_some((op, NodeId::NONE, argument));
    }
    if op == core::POW {
        let &[base, argument] = children else {
            return None;
        };
        let constant_base = !graph.depends_on(graph.find(base), symbol);
        return (constant_base && graph.depends_on(graph.find(argument), symbol)).then_some((op, base, argument));
    }
    None
}

/// `argument = c * rest` with a rational numeric factor `c`.
fn split_multiple(
    graph: &mut Graph,
    argument: NodeId,
) -> (BigRational, NodeId) {
    if graph.op(argument) == core::MUL {
        let factors = graph.children(argument).to_vec();
        let mut c = BigRational::one();
        let mut rest = Vec::new();
        for &f in &factors {
            match graph.number_of(f).and_then(Number::to_rational) {
                | Some(n) => c *= n,
                | None => rest.push(f),
            }
        }
        if !rest.is_empty() && c != BigRational::one() {
            let rest = product(graph, &rest);
            return (c, rest);
        }
    }
    (BigRational::one(), argument)
}

fn rebuild_exponential(
    graph: &mut Graph,
    op: OpId,
    base: NodeId,
    argument: NodeId,
) -> Option<NodeId> {
    if base == NodeId::NONE {
        graph.try_node(op, &[argument])
    } else {
        graph.try_node(op, &[base, argument])
    }
}

/// The integer `r >= 2` and `k >= 2` with `b = r^k`, for the largest `k`.
fn perfect_power(b: &BigInt) -> Option<(BigInt, u32)> {
    if *b < BigInt::from(4) {
        return None;
    }
    let bits = u32::try_from(b.bits()).ok()?;
    for k in (2..=bits.min(40)).rev() {
        let r = b.nth_root(k);
        if r.pow(k) == *b && r >= BigInt::from(2) {
            return Some((r, k));
        }
    }
    None
}

/// See the module documentation.
pub(super) fn exponentials(
    graph: &mut Graph,
    term: NodeId,
    symbol: SymbolId,
) -> NodeId {
    let exp = graph.ops().lookup("exp");
    // Pass 1: split sums in exponents, normalise integer bases.
    let mut sums = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
        let rebuilt_node = if children == graph.children(node) { node } else { graph.try_node(graph.op(node), children)? };
        let (op, base, argument) = exponential_parts(graph, rebuilt_node, symbol, exp)?;
        // Perfect-power base: 4^u = 2^(2u).
        if base != NodeId::NONE {
            if let Some(n) = graph.number_of(base).and_then(Number::to_rational) {
                if n.is_integer() {
                    if let Some((r, k)) = perfect_power(&n.to_integer()) {
                        let r_node = graph.num(Number::Int(r));
                        let k_node = graph.int(i64::from(k));
                        let scaled = graph.node(core::MUL, &[k_node, argument]);
                        return graph.try_node(op, &[r_node, scaled]);
                    }
                } else if n.numer().is_one() && n.is_positive() {
                    if let Some(inverse) = n.recip().to_integer().to_i64() {
                        let r_node = graph.int(inverse);
                        let minus_one = graph.int(-1);
                        let negated = graph.node(core::MUL, &[minus_one, argument]);
                        return graph.try_node(op, &[r_node, negated]);
                    }
                }
            }
        }
        if graph.op(argument) == core::ADD {
            let (mut dependent, mut free) = (Vec::new(), Vec::new());
            for &t in graph.children(argument) {
                if graph.depends_on(graph.find(t), symbol) { dependent.push(t) } else { free.push(t) }
            }
            if !dependent.is_empty() && !free.is_empty() {
                let dep_sum = if dependent.len() == 1 { dependent[0] } else { graph.node(core::ADD, &dependent) };
                let free_sum = if free.len() == 1 { free[0] } else { graph.node(core::ADD, &free) };
                let a = rebuild_exponential(graph, op, base, dep_sum)?;
                let b = rebuild_exponential(graph, op, base, free_sum)?;
                return Some(graph.node(core::MUL, &[a, b]));
            }
        }
        (rebuilt_node != node).then_some(rebuilt_node)
    };
    let mut current = map_term(graph, term, &mut sums);
    // Pass 2: the least common denominator of the numeric multiples per
    // (kind, base, argument without the multiple).
    let mut groups: HashMap<(OpId, NodeId, NodeId), BigInt> = HashMap::new();
    let mut stack = vec![current];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        stack.extend_from_slice(graph.children(n));
        if let Some((op, base, argument)) = exponential_parts(graph, n, symbol, exp) {
            let (c, rest) = split_multiple(graph, argument);
            let entry = groups.entry((op, base, rest)).or_insert_with(BigInt::one);
            *entry = entry.lcm(c.denom());
        }
    }
    let mut multiples = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
        let rebuilt_node = if children == graph.children(node) { node } else { graph.try_node(graph.op(node), children)? };
        let (op, base, argument) = exponential_parts(graph, rebuilt_node, symbol, exp)?;
        let (c, rest) = split_multiple(graph, argument);
        let d = groups.get(&(op, base, rest)).cloned().unwrap_or_else(BigInt::one);
        let multiplier = &c * BigRational::from_integer(d.clone());
        if !multiplier.is_integer() {
            return (rebuilt_node != node).then_some(rebuilt_node);
        }
        let m = multiplier.to_integer();
        if m.is_one() && d.is_one() {
            return (rebuilt_node != node).then_some(rebuilt_node);
        }
        let unit = if d.is_one() {
            rest
        } else {
            let scale = graph.num(Number::rat(BigRational::new(BigInt::one(), d)));
            graph.node(core::MUL, &[scale, rest])
        };
        let inner = rebuild_exponential(graph, op, base, unit)?;
        if m.is_one() {
            return Some(inner);
        }
        let m_node = graph.num(Number::Int(m));
        Some(graph.node(core::POW, &[inner, m_node]))
    };
    current = map_term(graph, current, &mut multiples);
    current
}

/// `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch` of the unknown written
/// through `exp`. `None` when none occurs.
pub(super) fn hyperbolic_to_exponentials(
    graph: &mut Graph,
    term: NodeId,
    symbol: SymbolId,
) -> Option<NodeId> {
    let names = ["sinh", "cosh", "tanh", "coth", "sech", "csch"];
    let ops: Vec<Option<OpId>> = names.iter().map(|n| graph.ops().lookup(n)).collect();
    let exp = graph.ops().lookup("exp")?;
    let mut changed = false;
    let mut f = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
        let which = ops.iter().position(|&o| o == Some(graph.op(node)))?;
        let &[a] = children else {
            return None;
        };
        if !graph.depends_on(graph.find(a), symbol) {
            return None;
        }
        changed = true;
        let minus_one = graph.int(-1);
        let neg = graph.node(core::MUL, &[minus_one, a]);
        let (p, m) = (graph.node(exp, &[a]), graph.node(exp, &[neg]));
        let half = graph.num(Number::fraction(1, 2)?);
        let minus_m = graph.node(core::MUL, &[minus_one, m]);
        let sinh = {
            let diff = graph.node(core::ADD, &[p, minus_m]);
            graph.node(core::MUL, &[half, diff])
        };
        let cosh = {
            let sum = graph.node(core::ADD, &[p, m]);
            graph.node(core::MUL, &[half, sum])
        };
        let inverse = |graph: &mut Graph, v: NodeId| graph.node(core::POW, &[v, minus_one]);
        Some(match which {
            | 0 => sinh,
            | 1 => cosh,
            | 2 => {
                let c = inverse(graph, cosh);
                graph.node(core::MUL, &[sinh, c])
            },
            | 3 => {
                let s = inverse(graph, sinh);
                graph.node(core::MUL, &[cosh, s])
            },
            | 4 => inverse(graph, cosh),
            | _ => inverse(graph, sinh),
        })
    };
    let out = map_term(graph, term, &mut f);
    changed.then_some(out)
}

/// `b^u` (constant base, `u` involving the unknown) as `exp(u ln b)`.
pub(super) fn powers_to_exp(
    graph: &mut Graph,
    term: NodeId,
    symbol: SymbolId,
) -> Option<NodeId> {
    let (exp, ln) = (graph.ops().lookup("exp")?, graph.ops().lookup("ln")?);
    let mut changed = false;
    let mut f = |graph: &mut Graph, node: NodeId, children: &[NodeId]| -> Option<NodeId> {
        if graph.op(node) != core::POW {
            return None;
        }
        let &[base, argument] = children else {
            return None;
        };
        if graph.depends_on(graph.find(base), symbol) || !graph.depends_on(graph.find(argument), symbol) {
            return None;
        }
        changed = true;
        let log = graph.node(ln, &[base]);
        let scaled = graph.node(core::MUL, &[argument, log]);
        Some(graph.node(exp, &[scaled]))
    };
    let out = map_term(graph, term, &mut f);
    changed.then_some(out)
}
