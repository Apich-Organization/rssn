//! Building concrete rssn terms from the exact coefficient data used by the
//! Hermite reduction.

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;

use super::field::Frac;

/// The literal node for the rational number `q`.
pub(super) fn rat_node(
    graph: &mut Graph,
    q: &BigRational,
) -> NodeId {
    graph.num(Number::rat(q.clone()))
}

/// `-1`.
fn minus_one(graph: &mut Graph) -> NodeId {
    graph.int(-1)
}

/// Sum of `terms`, folding the empty sum to zero and a single term to
/// itself.
pub(super) fn add_all(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

/// Additive inverse.
pub(super) fn neg_node(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let m = minus_one(graph);
    graph.node(core::MUL, &[m, node])
}

/// Difference.
pub(super) fn sub_node(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let nb = neg_node(graph, b);
    graph.node(core::ADD, &[a, nb])
}

/// `a / b`.
pub(super) fn div_node(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    if graph.number_of(b).is_some_and(Number::is_one) {
        return a;
    }
    let m = minus_one(graph);
    let inv = graph.node(core::POW, &[b, m]);
    graph.node(core::MUL, &[a, inv])
}

/// `value^(1/2)` for an arbitrary term.
pub(super) fn sqrt_node(
    graph: &mut Graph,
    value: NodeId,
) -> NodeId {
    let half = Number::fraction(1, 2).unwrap_or_else(|| Number::from(0));
    let half_node = graph.num(half);
    graph.node(core::POW, &[value, half_node])
}

/// The polynomial `p` evaluated at the term `x`.
pub(super) fn poly_node(
    graph: &mut Graph,
    x: NodeId,
    p: &[BigRational],
) -> NodeId {
    let mut terms = Vec::new();
    for (k, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let coeff = rat_node(graph, c);
        if k == 0 {
            terms.push(coeff);
            continue;
        }
        let power = u32::try_from(k).unwrap_or(u32::MAX);
        let base = if power == 1 {
            x
        } else {
            let e = graph.int(i64::from(power));
            graph.node(core::POW, &[x, e])
        };
        terms.push(graph.node(core::MUL, &[coeff, base]));
    }
    add_all(graph, &terms)
}

/// The rational function `f` of the term `x`.
pub(super) fn frac_node(
    graph: &mut Graph,
    x: NodeId,
    f: &Frac,
) -> NodeId {
    let num = poly_node(graph, x, &f.num);
    if f.den.len() == 1 && f.den.first().is_some_and(One::is_one) {
        return num;
    }
    let den = poly_node(graph, x, &f.den);
    div_node(graph, num, den)
}

/// `ln(argument)`.
pub(super) fn log_node(
    graph: &mut Graph,
    argument: NodeId,
) -> Option<NodeId> {
    let ln = graph.ops().lookup("ln")?;
    Some(graph.node(ln, &[argument]))
}

/// `atan(argument)`.
pub(super) fn atan_node(
    graph: &mut Graph,
    argument: NodeId,
) -> Option<NodeId> {
    let atan = graph.ops().lookup("atan")?;
    Some(graph.node(atan, &[argument]))
}

/// `asin(argument)`.
pub(super) fn asin_node(
    graph: &mut Graph,
    argument: NodeId,
) -> Option<NodeId> {
    let asin = graph.ops().lookup("asin")?;
    Some(graph.node(asin, &[argument]))
}

/// `integral(integrand, x)`.
pub(super) fn integral_node(
    graph: &mut Graph,
    integrand: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let integral = graph.ops().lookup("integral")?;
    Some(graph.node(integral, &[integrand, x]))
}
