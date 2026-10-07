//! Small builders shared by the PDE solvers.

use crate::graph::Cx;
use crate::graph::SymbolId;
use num_complex::Complex64;
use std::collections::HashMap;
use crate::graph::Env;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;

/// `a / b`.
pub(super) fn div(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let inverse = powi(cx.graph, b, -1);
    mul(cx.graph, &[a, inverse])
}

/// The application of the registered operator `name` to `args`.
pub(super) fn call(
    cx: &mut Cx<'_>,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = cx.graph.ops().lookup(name)?;
    Some(cx.graph.node(op, args))
}

/// `pi`.
pub(super) fn pi(cx: &mut Cx<'_>) -> Option<NodeId> {
    call(cx, "pi", &[])
}

/// `oo`.
pub(super) fn infinity(cx: &mut Cx<'_>) -> Option<NodeId> {
    call(cx, "oo", &[])
}

/// The imaginary unit.
pub(super) fn imaginary(cx: &mut Cx<'_>) -> Option<NodeId> {
    call(cx, "I", &[])
}

/// The exact fraction `n / d`.
pub(super) fn fraction(
    cx: &mut Cx<'_>,
    n: i64,
    d: i64,
) -> Option<NodeId> {
    Some(cx.graph.num(Number::fraction(n, d)?))
}

/// `x^(1/2)`.
pub(super) fn sqrt(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    let half = fraction(cx, 1, 2)?;
    Some(pow(cx.graph, x, half))
}

/// `defint(body, var, lo, hi)`.
pub(super) fn defint(
    cx: &mut Cx<'_>,
    body: NodeId,
    var: NodeId,
    lo: NodeId,
    hi: NodeId,
) -> Option<NodeId> {
    call(cx, "defint", &[body, var, lo, hi])
}

/// `sum(body, var, lo, hi)`.
pub(super) fn series(
    cx: &mut Cx<'_>,
    body: NodeId,
    var: NodeId,
    lo: NodeId,
    hi: NodeId,
) -> Option<NodeId> {
    call(cx, "sum", &[body, var, lo, hi])
}

/// Whether `n` is the number zero.
pub(super) fn is_zero_number(
    graph: &Graph,
    n: NodeId,
) -> bool {
    graph.number_of(n).is_some_and(Number::is_zero)
}

/// Whether some node of `term` has the operator named `name`.
pub(super) fn contains_op(
    graph: &Graph,
    term: NodeId,
    name: &str,
) -> bool {
    let Some(op) = graph.ops().lookup(name) else {
        return false;
    };
    let mut stack = vec![term];
    while let Some(n) = stack.pop() {
        if graph.op(n) == op {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// An environment binding every free symbol of `node` to a distinct
/// moderate number: integers for symbols known to be integers.
pub(super) fn sample_env(
    graph: &Graph,
    node: NodeId,
    k: u32,
) -> Env {
    let symbols = graph.free_symbols(graph.find(node)).to_vec();
    let mut env = Env::numeric(0.0);
    for &s in &symbols {
        let value = if graph.assumption(s).has(Facts::INTEGER) {
            f64::from(2 + (s.raw() + k) % 5)
        } else {
            0.31 + 0.17 * f64::from(k) + 0.13 * f64::from(s.raw() % 7)
        };
        env.bind(s, value);
    }
    env
}

/// The value of `node` at the sample point `k`.
pub(super) fn sample(
    graph: &Graph,
    node: NodeId,
    k: u32,
) -> Option<f64> {
    let env = sample_env(graph, node, k);
    graph.eval(node, &env).filter(|v| v.is_finite())
}

/// The derivative `∂^α f` as unevaluated `diff` nodes.
pub(super) fn jet(
    cx: &mut Cx<'_>,
    f: NodeId,
    vars: &[NodeId],
    index: &[u32],
) -> Option<NodeId> {
    let diff = cx.graph.ops().lookup("diff")?;
    let mut out = f;
    for (k, &order) in index.iter().enumerate() {
        for _ in 0..order {
            out = cx.graph.node(diff, &[out, *vars.get(k)?]);
        }
    }
    Some(out)
}

/// The magnitude of `node` at the sample point `k`: the real value, or
/// the complex one when the term contains complex numbers.
pub(super) fn sample_abs(
    graph: &Graph,
    node: NodeId,
    k: u32,
) -> Option<f64> {
    if let Some(v) = sample(graph, node, k) {
        return Some(v.abs());
    }
    let symbols = graph.free_symbols(graph.find(node)).to_vec();
    let mut bindings: HashMap<SymbolId, Complex64> = HashMap::new();
    for &s in &symbols {
        let value = if graph.assumption(s).has(Facts::INTEGER) {
            f64::from(2 + (s.raw() + k) % 5)
        } else {
            0.31 + 0.17 * f64::from(k) + 0.13 * f64::from(s.raw() % 7)
        };
        bindings.insert(s, Complex64::new(value, 0.0));
    }
    graph.eval_complex(node, &bindings).filter(|v| v.re.is_finite() && v.im.is_finite()).map(|v| v.norm())
}
