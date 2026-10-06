//! Numeric optimisation on the graph.
//!
//! | operator | value |
//! |---|---|
//! | `nminimize(f, list(x, ...), list(x0, ...))` | `list(list(x*, ...), f(x*))`: a local minimum by BFGS from the start point |
//! | `nmaximize(f, list(x, ...), list(x0, ...))` | the same for a local maximum |
//! | `nminimize_global(f, list(x, ...), list(list(lo, hi), ...))` | a global minimum over the box by differential evolution, polished by BFGS |
//! | `nmaximize_global(f, list(x, ...), list(list(lo, hi), ...))` | the same for the maximum |
//!
//! The objective and its symbolic gradient are compiled together with the
//! [current](crate::backend::current) backend (machine code under the JIT),
//! so the minimiser never evaluates finite differences. Symbols other than
//! the variables take their values from the bindings of the run.

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::kernels::optim;
use crate::rules::poly::best;

use super::super::calculus::derivative;
use super::list_items;

/// Gradient norm at which BFGS stops.
const TOLERANCE: f64 = 1e-12;
const MAX_ITERATIONS: usize = 2_000;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let heavy = |name: &str| OpDescriptor::new(name, Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100);
    for (name, sign, global) in [
        ("nminimize", 1.0, false),
        ("nmaximize", -1.0, false),
        ("nminimize_global", 1.0, true),
        ("nmaximize_global", -1.0, true),
    ] {
        let op = i.op(heavy(name))?;
        i.kernel(&format!("optimize/{name}"), Tier::Reduce, Minimize { op, sign, global });
    }
    Ok(())
}

struct Minimize {
    op: OpId,
    /// `1` to minimise, `-1` to maximise.
    sign: f64,
    global: bool,
}

impl Kernel for Minimize {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, vars, start] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        self.minimize(cx, f, vars, start).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        // The objective or the start point may simplify to something
        // evaluable later.
        true
    }
}

/// Bounds of a box `list(list(lo, hi), ...)`.
fn bounds(
    cx: &mut Cx<'_>,
    node: NodeId,
    dimension: usize,
) -> Option<Vec<(f64, f64)>> {
    let items = list_items(cx.graph, node)?;
    if items.len() != dimension {
        return None;
    }
    let mut out = Vec::with_capacity(items.len());
    for item in items {
        let pair = list_items(cx.graph, item)?;
        let &[lo, hi] = pair.as_slice() else {
            return None;
        };
        let (lo, hi) = (value(cx, lo)?, value(cx, hi)?);
        if lo >= hi {
            return None;
        }
        out.push((lo, hi));
    }
    Some(out)
}

fn value(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<f64> {
    let term = best(cx.graph, node)?;
    cx.graph.eval(term, cx.env).filter(|v| v.is_finite())
}

impl Minimize {
    fn minimize(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        vars: NodeId,
        start: NodeId,
    ) -> Option<NodeId> {
        let vars = list_items(cx.graph, vars)?;
        let symbols: Vec<SymbolId> = vars.iter().map(|&v| cx.graph.symbol_of(v)).collect::<Option<_>>()?;
        let objective = best(cx.graph, f)?;
        // Every other free symbol must be bound.
        let bindings: Vec<(SymbolId, f64)> =
            cx.env.bindings().iter().copied().filter(|(s, _)| !symbols.contains(s)).collect();
        let free = cx.graph.free_symbols(cx.graph.find(objective)).to_vec();
        if free.iter().any(|s| !symbols.contains(s) && !bindings.iter().any(|(b, _)| b == s)) {
            return None;
        }
        let mut roots = vec![objective];
        for &x in &vars {
            let partial = derivative(cx.graph, objective, x)?;
            let partial = cx.simplify(partial);
            roots.push(best(cx.graph, partial)?);
        }
        let params: Vec<SymbolId> = bindings.iter().map(|(s, _)| *s).collect();
        let values: Vec<f64> = bindings.iter().map(|(_, v)| *v).collect();
        let compiled = crate::backend::current().compile_multi(cx.graph, &roots, &symbols, &params).ok()?;
        let n = symbols.len();
        let sign = self.sign;
        let run = |x: &[f64]| {
            let mut out = vec![0.0; n + 1];
            compiled.eval(x, &values, &mut out);
            out
        };
        let objective = |x: &[f64]| run(x).first().map_or(f64::NAN, |v| sign * v);
        let gradient = |x: &[f64]| run(x).iter().skip(1).map(|g| sign * g).collect::<Vec<f64>>();
        let x0 = if self.global {
            let bounds = bounds(cx, start, n)?;
            let guarded = |x: &[f64]| {
                let v = objective(x);
                if v.is_finite() { v } else { f64::INFINITY }
            };
            let population = (15 * n).max(20);
            optim::differential_evolution(guarded, &bounds, population, 300, 0x5eed).x
        } else {
            let items = list_items(cx.graph, start)?;
            if items.len() != n {
                return None;
            }
            items.into_iter().map(|e| value(cx, e)).collect::<Option<Vec<f64>>>()?
        };
        let result = optim::bfgs(objective, Some(gradient), &x0, TOLERANCE, MAX_ITERATIONS);
        let fx = sign * result.fx;
        if !fx.is_finite() || result.x.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let point: Vec<NodeId> = result.x.iter().map(|&v| cx.graph.float(v)).collect();
        let point = cx.graph.node(core::LIST, &point);
        let fx = cx.graph.float(fx);
        Some(cx.graph.node(core::LIST, &[point, fx]))
    }
}

#[cfg(test)]
mod tests {
    use crate::rules::optimize::optimize;
    use crate::rules::standard;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    /// The numbers in a result string.
    fn numbers(s: &str) -> Vec<f64> {
        s.split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == 'e'))
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    #[test]
    fn local_minimum_of_rosenbrock() {
        let out = simplify(&[optimize()], "nminimize((1 - x)^2 + 100*(y - x^2)^2, list(x, y), list(-1.2, 1))");
        let v = numbers(&out);
        assert!(v.len() == 3, "{out}");
        assert!((v[0] - 1.0).abs() < 1e-6 && (v[1] - 1.0).abs() < 1e-6 && v[2].abs() < 1e-10, "{out}");
    }

    #[test]
    fn maximum_and_special_functions() {
        // x*exp(-x^2) peaks at 1/sqrt(2).
        let out = simplify(&standard(), "nmaximize(x*exp(-x^2), list(x), list(0.3))");
        let v = numbers(&out);
        assert!((v[0] - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-8, "{out}");
        // lgamma has its minimum at x = 1.4616321449683623.
        let out = simplify(&standard(), "nminimize(lgamma(x), list(x), list(2))");
        let v = numbers(&out);
        assert!((v[0] - 1.461_632_144_968_362_3).abs() < 1e-7, "{out}");
    }

    #[test]
    fn global_minimum_of_a_multimodal_function() {
        // Rastrigin in two dimensions: global minimum 0 at the origin among
        // many local minima.
        let out = simplify(
            &standard(),
            "nminimize_global(20 + x^2 + y^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y), list(x, y), list(list(-5, 5), list(-5, 5)))",
        );
        let v = numbers(&out);
        assert!(v.len() == 3 && v[0].abs() < 1e-6 && v[1].abs() < 1e-6 && v[2].abs() < 1e-9, "{out}");
    }

    #[test]
    fn unbound_parameters_leave_the_request() {
        let (out, reduced) = reduce_with(&[optimize()], "nminimize((x - a)^2, list(x), list(0))", &[]);
        assert!(!reduced && out.starts_with("nminimize("), "{out}");
    }
}
