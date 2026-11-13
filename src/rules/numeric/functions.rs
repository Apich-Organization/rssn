//! Numeric kernels that take a function, exposed as operators over
//! expressions.
//!
//! The kernels of [`crate::kernels`] want closures. Here the closure is an
//! expression in named variables, compiled once with the current backend
//! (the route `odeint` takes), and, where the kernel needs derivatives, the
//! derivatives are taken symbolically.
//!
//! | operator | value |
//! |---|---|
//! | `oscillatory_integral(f, g, x, a, b, w, n)` | `list(re, im)` of `int_a^b f(x) exp(i w g(x)) dx` by Levin collocation at `n` Chebyshev points |
//! | `anderson_fixed_point(g, vars, x0, tol)` | the fixed point `x = g(x)` by Anderson acceleration; `g`, `vars`, `x0` are lists (or single items) |
//! | `newton_cg_minimize(f, vars, x0)` | the minimiser of `f` by the trust-region Newton-CG method with the exact symbolic gradient and Hessian |
//!
//! Every parameter other than the variables must be a number; anything else,
//! a failed convergence, or an expression the backend cannot compile leaves
//! the call as a request.

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::kernels::lowrank;
use crate::kernels::newton_cg;
use crate::kernels::oscillatory;
use crate::rules::calculus::derivative;
use crate::rules::poly::best;

/// Most variables of a vector problem.
pub(super) const MAX_VARIABLES: usize = 12;

/// The operators of this module, in the order of [`Functions::ops`].
#[derive(Clone, Copy)]
pub(super) struct Functions {
    pub(super) oscillatory: OpId,
    pub(super) anderson: OpId,
    pub(super) newton_cg: OpId,
}

/// The nodes a list holds, or the single node.
pub(super) fn items(
    graph: &Graph,
    node: NodeId,
) -> Vec<NodeId> {
    graph
        .enodes(graph.find(node))
        .find(|&n| graph.op(n) == core::LIST)
        .map_or_else(|| vec![node], |list| graph.children(list).to_vec())
}

pub(super) fn symbols(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<SymbolId>> {
    let list = items(graph, node);
    (!list.is_empty() && list.len() <= MAX_VARIABLES)
        .then(|| list.iter().map(|&n| graph.symbol_of(n)).collect())
        .flatten()
}

pub(super) fn number(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<f64> {
    let env = cx.env;
    let value = best(cx.graph, node).and_then(|t| cx.graph.eval(t, env))?;
    value.is_finite().then_some(value)
}

pub(super) fn numbers(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<Vec<f64>> {
    items(cx.graph, node)
        .into_iter()
        .map(|n| number(cx, n))
        .collect()
}

/// A compiled function of the variables, writing into an output slice.
pub(super) type Compiled = Box<dyn Fn(&[f64], &mut [f64])>;

/// The terms compiled into one function of `vars`.
pub(super) fn compile(
    graph: &mut Graph,
    terms: &[NodeId],
    vars: &[SymbolId],
) -> Option<Compiled> {
    let terms: Vec<NodeId> = terms
        .iter()
        .map(|&t| best(graph, t))
        .collect::<Option<_>>()?;
    let compiled = crate::backend::current()
        .compile_multi(graph, &terms, vars, &[])
        .ok()?;
    Some(Box::new(move |args: &[f64], out: &mut [f64]| {
        compiled.eval(args, &[], out);
    }))
}

pub(super) fn list(
    graph: &mut Graph,
    values: &[f64],
) -> NodeId {
    let nodes: Vec<NodeId> = values.iter().map(|&v| graph.float(v)).collect();
    graph.node(core::LIST, &nodes)
}

impl Functions {
    fn oscillatory_integral(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[f, g, x, a, b, w, n] = args else {
            return None;
        };
        let var = cx.graph.symbol_of(x)?;
        let (a, b, w, n) = (
            number(cx, a)?,
            number(cx, b)?,
            number(cx, w)?,
            number(cx, n)?,
        );
        ((2.0..=512.0).contains(&n) && n.fract() == 0.0).then_some(())?;
        let dg = derivative(cx.graph, g, x)?;
        let value = compile(cx.graph, &[f, g, dg], &[var])?;
        let at = |i: usize, t: f64| {
            let mut out = [0.0; 3];
            value(&[t], &mut out);
            out[i]
        };
        let (fv, gv, dgv) = (|t| at(0, t), |t| at(1, t), |t| at(2, t));
        let result = oscillatory::levin(&fv, &gv, &dgv, a, b, w, n as usize).ok()?;
        (result.re.is_finite() && result.im.is_finite()).then_some(())?;
        Some(list(cx.graph, &[result.re, result.im]))
    }

    fn anderson_fixed_point(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[g, vars, x0, tol] = args else {
            return None;
        };
        let vars = symbols(cx.graph, vars)?;
        let rhs = items(cx.graph, g);
        let start = numbers(cx, x0)?;
        let tol = number(cx, tol)?;
        (rhs.len() == vars.len() && start.len() == vars.len() && tol > 0.0).then_some(())?;
        let map = compile(cx.graph, &rhs, &vars)?;
        let n = vars.len();
        let step = |x: &[f64]| {
            let mut out = vec![0.0; n];
            map(x, &mut out);
            out
        };
        let found = lowrank::anderson(&step, &start, 5, 1.0, tol, 10_000).ok()?;
        Some(list(cx.graph, &found.x))
    }

    fn newton_cg_minimize(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[f, vars, x0] = args else {
            return None;
        };
        let var_nodes = items(cx.graph, vars);
        let vars = symbols(cx.graph, vars)?;
        let start = numbers(cx, x0)?;
        let n = vars.len();
        (start.len() == n).then_some(())?;
        // The gradient, then the Hessian row by row.
        let mut terms = vec![f];
        let mut gradient = Vec::with_capacity(n);
        for &v in &var_nodes {
            gradient.push(derivative(cx.graph, f, v)?);
        }
        terms.extend(&gradient);
        for &gi in &gradient {
            for &v in &var_nodes {
                terms.push(derivative(cx.graph, gi, v)?);
            }
        }
        let all = compile(cx.graph, &terms, &vars)?;
        let table = |x: &[f64]| {
            let mut out = vec![0.0; 1 + n + n * n];
            all(x, &mut out);
            out
        };
        let result = newton_cg::minimize(
            |x| table(x)[0],
            |x, out| out.copy_from_slice(&table(x)[1..=n]),
            |x, d, out| {
                let t = table(x);
                for (i, o) in out.iter_mut().enumerate() {
                    *o = (0..n).map(|j| t[1 + n + i * n + j] * d[j]).sum();
                }
            },
            &start,
            &newton_cg::NewtonCgOptions::default(),
        )
        .ok()?;
        matches!(
            result.status,
            newton_cg::NewtonCgStatus::FirstOrderConverged
        )
        .then_some(())?;
        Some(list(cx.graph, &result.x))
    }
}

impl Kernel for Functions {
    fn ops(&self) -> Vec<OpId> {
        vec![self.oscillatory, self.anderson, self.newton_cg]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let args = cx.graph.children(node).to_vec();
        let found = if op == self.oscillatory {
            Self::oscillatory_integral(cx, &args)
        } else if op == self.anderson {
            Self::anderson_fixed_point(cx, &args)
        } else if op == self.newton_cg {
            Self::newton_cg_minimize(cx, &args)
        } else {
            None
        };
        found.map_or(Outcome::Pass, Outcome::Equal)
    }
}
