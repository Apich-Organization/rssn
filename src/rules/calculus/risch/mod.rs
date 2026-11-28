//! Risch–Bronstein integration over a tower of monomial extensions.
//!
//! The integrator in `integrate` is a search over rewriting heuristics. This
//! module adds the decision procedure for integrands that live in a tower
//!
//! ```text
//! Q(x) = K_1 < K_2 < ... < K_n,   theta_i' = a_i*theta_i + b_i,
//! ```
//!
//! where every `theta_i` is a monomial over the field below it: a logarithm
//! (`theta = log u`, so `a_i = 0` and `b_i = u'/u`) or an exponential
//! (`theta = exp u`, so `a_i = u'` and `b_i = 0`). The variable `x` is the
//! first level of the tower, which makes rational functions the base case of
//! the same code.
//!
//! Integration at the top level follows Bronstein's algorithm:
//!
//! 1. **Hermite reduction** (`hermite`): `f = g' + A/D` with `D` squarefree.
//!    The linear system `A*D_2 = (B'*D_2 - B*D_2')*D_sqf + C*D_2^2` is solved
//!    over the coefficient field, so a repeated factor of the denominator
//!    costs one elimination instead of a search.
//! 2. **The logarithmic part** (`log_part`): the residues are the roots of
//!    `res_theta(D, A - z*D'*theta')`. The resultant is interpolated from
//!    `deg D + 1` evaluations, and its rational roots `c` give the logarithms
//!    `c*log(gcd(D, A - c*D'*theta'))`. When a residue is not rational the
//!    integral needs an algebraic extension and the module declines.
//! 3. **The polynomial part** (`poly_part`): what remains is a polynomial in
//!    `theta` (a Laurent polynomial for an exponential), and its
//!    antiderivative is a bounded ansatz over the monomials of the tower
//!    solved by one linear system over `Q`. The ansatz carries the constant
//!    freedom of the top coefficient, which a top-down recursion would have
//!    to search for. `∫ log(x)/x = log(x)^2/2` needs exactly that freedom.
//!
//! For an exponential monomial the place `theta = 0` is not a pole of the
//! differential (it lies above `x = -infinity`), so the principal part there
//! is integrated as negative powers instead of as a logarithm; that is what
//! `∫ dx/e^x = -e^{-x}` needs.
//!
//! Every candidate is checked by differentiating it back (`Integral` is only
//! returned when `p' == h` holds exactly), and the caller verifies the answer
//! numerically like every other antiderivative.

mod alg;
mod elem;
mod flat;
mod poly;
mod solve;
mod tower;

#[cfg(test)]
mod tests;

use std::collections::HashMap;

use num_rational::BigRational;

use super::diff::Differentiate;
use super::integrate::Functions;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::SymbolId;
use crate::graph::op::core;
use elem::Coeff;
use elem::Elem;
use tower::Kind;
use tower::Tower;

/// Rational numbers, the constant field of the tower.
pub(super) type Q = BigRational;

/// Deepest tower the module builds. Each extra level multiplies the cost of
/// the Hermite system and the residue interpolation, so the bound is a budget
/// rather than a mathematical limit.
const MAX_LEVELS: usize = 5;

/// A generator of the tower together with the node that spells it.
struct Generator {
    node: NodeId,
    kind: Kind,
}

/// An antiderivative of `f` with respect to `x`, or `None` when the tower
/// algorithm does not decide it.
pub(super) fn integrate(
    cx: &mut Cx<'_>,
    functions: Functions,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let symbol = cx.graph.symbol_of(x)?;
    let mut generators = Vec::new();
    collect(cx.graph, f, symbol, functions, &mut generators);
    if generators.len() > MAX_LEVELS {
        return None;
    }
    let (tower, index) = build_tower(cx, functions, &generators, x, symbol)?;
    let mut memo = HashMap::new();
    let depth = tower.depth();
    let element = elem_of(cx.graph, f, depth, x, &index, &mut memo)?;
    let found = alg::integrate(&element, &tower)?;
    let nodes: Vec<NodeId> = generators.iter().map(|g| g.node).collect();
    let mut pieces = Vec::new();
    if !found.rational.is_zero() {
        pieces.push(term_of(cx.graph, &found.rational, x, &nodes)?);
    }
    for (coefficient, argument) in &found.logs {
        let inner = term_of(cx.graph, argument, x, &nodes)?;
        let log = cx.graph.node(functions.ln, &[inner]);
        let scale = cx.graph.num(Number::rat(coefficient.clone()));
        pieces.push(product(cx.graph, &[scale, log]));
    }
    match pieces.len() {
        | 0 => Some(cx.graph.int(0)),
        | 1 => pieces.first().copied(),
        | _ => Some(sum(cx.graph, &pieces)),
    }
}

/// Collects the exponential and logarithmic generators of `f`, ordered by
/// nesting depth.
fn collect(
    graph: &Graph,
    f: NodeId,
    symbol: SymbolId,
    functions: Functions,
    out: &mut Vec<Generator>,
) {
    let mut found: Vec<(NodeId, Kind, usize)> = Vec::new();
    let mut memo: HashMap<NodeId, usize> = HashMap::new();
    height(graph, f, symbol, functions, &mut found, &mut memo);
    // Shallower generators come first: each one is a monomial over the
    // field its predecessors generate.
    found.sort_by_key(|&(node, _, level)| (level, node));
    for (node, kind, _) in found {
        if out.iter().any(|g| g.node == node) {
            continue;
        }
        out.push(Generator { node, kind });
    }
}

/// Height of `node`: one more than the deepest generator inside it.
fn height(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
    functions: Functions,
    found: &mut Vec<(NodeId, Kind, usize)>,
    memo: &mut HashMap<NodeId, usize>,
) -> usize {
    if let Some(&value) = memo.get(&node) {
        return value;
    }
    let op = graph.op(node);
    let mut inner = 0;
    for &child in graph.children(node) {
        inner = inner.max(height(graph, child, symbol, functions, found, memo));
    }
    let value = if (op == functions.exp || op == functions.ln)
        && graph.depends_on(graph.find(node), symbol)
    {
        let kind = if op == functions.exp {
            Kind::Exp
        } else {
            Kind::Log
        };
        let level = inner + 1;
        found.push((node, kind, level));
        level
    } else {
        inner
    };
    memo.insert(node, value);
    value
}

/// Builds the tower and the map from generator node to level.
fn build_tower(
    cx: &mut Cx<'_>,
    functions: Functions,
    generators: &[Generator],
    x: NodeId,
    symbol: SymbolId,
) -> Option<(Tower, HashMap<NodeId, usize>)> {
    let mut tower = Tower::with_variable();
    let mut index: HashMap<NodeId, usize> = HashMap::new();
    for (position, generator) in generators.iter().enumerate() {
        let level = position + 2;
        let argument = *cx.graph.children(generator.node).first()?;
        let derivative =
            Differentiate { diff: functions.diff }.derive(cx.graph, argument, symbol, x);
        let mut memo = HashMap::new();
        let slope = elem_of(cx.graph, derivative, level - 1, x, &index, &mut memo)?;
        let monomial = match generator.kind {
            | Kind::Exp => (slope, Elem::zero_at(level - 1)),
            | Kind::Log => {
                let mut memo = HashMap::new();
                let value = elem_of(cx.graph, argument, level - 1, x, &index, &mut memo)?;
                (Elem::zero_at(level - 1), slope.mul(&value.inv()?))
            },
            | Kind::Primitive => return None,
        };
        tower.push(generator.kind, monomial.0, monomial.1);
        index.insert(generator.node, level);
    }
    Some((tower, index))
}

/// Converts a term to an element of level `depth`.
fn elem_of(
    graph: &Graph,
    node: NodeId,
    depth: usize,
    x: NodeId,
    index: &HashMap<NodeId, usize>,
    memo: &mut HashMap<NodeId, Elem>,
) -> Option<Elem> {
    if let Some(value) = memo.get(&node) {
        return Some(value.clone());
    }
    let value = if let Some(number) = graph.number_of(node) {
        Elem::Const(number.to_rational()?)
    } else if node == x {
        Elem::gen_at(1)
    } else if let Some(&level) = index.get(&node) {
        if level > depth {
            return None;
        }
        Elem::gen_at(level)
    } else {
        let children = graph.children(node);
        match graph.op(node) {
            | op if op == core::ADD => {
                let mut acc = Elem::zero_at(depth);
                for &child in children {
                    acc = acc.add(&elem_of(graph, child, depth, x, index, memo)?);
                }
                acc
            },
            | op if op == core::MUL => {
                let mut acc = Elem::one_at(depth);
                for &child in children {
                    acc = acc.mul(&elem_of(graph, child, depth, x, index, memo)?);
                }
                acc
            },
            | op if op == core::POW => {
                let [base, exponent] = children else {
                    return None;
                };
                let power = graph.number_of(*exponent)?.to_i64()?;
                elem_of(graph, *base, depth, x, index, memo)?.pow_i(power)?
            },
            | _ => return None,
        }
    };
    memo.insert(node, value.clone());
    Some(value)
}

/// Converts an element of the tower back to a term.
fn term_of(
    graph: &mut Graph,
    element: &Elem,
    x: NodeId,
    generators: &[NodeId],
) -> Option<NodeId> {
    match element {
        | Elem::Const(value) => Some(graph.num(Number::rat(value.clone()))),
        | Elem::Ext(ext) => {
            let theta = if ext.depth == 1 {
                x
            } else {
                *generators.get(ext.depth - 2)?
            };
            let numerator = polynomial_of(graph, &ext.num, theta, x, generators)?;
            let denominator = polynomial_of(graph, &ext.den, theta, x, generators)?;
            if denominator == graph.int(1) {
                return Some(numerator);
            }
            let inverse = {
                let minus_one = graph.int(-1);
                graph.node(core::POW, &[denominator, minus_one])
            };
            Some(product(graph, &[numerator, inverse]))
        },
    }
}

/// `sum c_i theta^i` by Horner's rule.
fn polynomial_of(
    graph: &mut Graph,
    coefficients: &[Elem],
    theta: NodeId,
    x: NodeId,
    generators: &[NodeId],
) -> Option<NodeId> {
    let mut acc = graph.int(0);
    for coefficient in coefficients.iter().rev() {
        let shifted = product(graph, &[acc, theta]);
        acc = if coefficient.is_zero() {
            shifted
        } else {
            let value = term_of(graph, coefficient, x, generators)?;
            sum(graph, &[shifted, value])
        };
    }
    Some(acc)
}

fn sum(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn product(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}
