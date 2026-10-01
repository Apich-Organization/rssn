//! Symbolic optimisation: critical points and their classification.
//!
//! | operator | value |
//! |---|---|
//! | `find_extrema(f, list(x, y, ...))` | `list(list(point, kind, value), ...)`: every real critical point (where the gradient vanishes, solved exactly), its kind — `local_min`, `local_max`, `saddle` or `degenerate` — and `f` there |
//! | `find_constrained_extrema(f, list(g1, ...), list(x, ...))` | `list(list(point, multipliers, kind, value), ...)` for the extrema of `f` subject to `g_i = 0`, by Lagrange multipliers |
//!
//! A point is `list(x0, y0, ...)` in the order of the variables. With one
//! variable the kind comes from the first non-vanishing higher
//! derivative; with several from Sylvester's criterion on the Hessian, and
//! for constrained problems from the bordered Hessian (two variables, one
//! constraint) when it is not degenerate (`critical` otherwise).

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;

use super::calculus::derivative;
use super::linalg::linalg;

/// The optimisation rule set.
#[must_use]
pub fn optimize() -> RuleSet {
    RuleSet::new("optimize", install).needs(linalg())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let heavy = |name: &str, arity: u8| OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
    let op = i.op(heavy("find_extrema", 2))?;
    i.kernel("optimize/find_extrema", Tier::Reduce, Extrema { op, constrained: false });
    let op = i.op(heavy("find_constrained_extrema", 3))?;
    i.kernel("optimize/find_constrained_extrema", Tier::Reduce, Extrema { op, constrained: true });
    Ok(())
}

struct Extrema {
    op: OpId,
    constrained: bool,
}

impl Kernel for Extrema {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = if self.constrained { constrained(cx, &args) } else { unconstrained(cx, &args) };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn list_items(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let node = best(graph, node)?;
    Some(if graph.op(node) == core::LIST { graph.children(node).to_vec() } else { vec![node] })
}

/// The solutions of `equations = 0` for `unknowns`, each a tuple in the
/// order of the unknowns.
fn solve_system(
    cx: &mut Cx<'_>,
    equations: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let solve = cx.graph.ops().lookup("solve")?;
    let zero = cx.graph.int(0);
    let eqs: Vec<NodeId> = equations.iter().map(|&e| cx.graph.node(core::EQ, &[e, zero])).collect();
    let eqs = cx.graph.node(core::LIST, &eqs);
    let vars = cx.graph.node(core::LIST, unknowns);
    let request = cx.graph.node(solve, &[eqs, vars]);
    let solved = cx.simplify(request);
    if cx.graph.op(solved) != core::LIST {
        return None;
    }
    let mut out = Vec::new();
    for tuple in cx.graph.children(solved).to_vec() {
        let values = if cx.graph.op(tuple) == core::LIST { cx.graph.children(tuple).to_vec() } else { vec![tuple] };
        if values.len() != unknowns.len() {
            return None;
        }
        out.push(values);
    }
    Some(out)
}

fn substitute_point(
    graph: &mut Graph,
    term: NodeId,
    vars: &[NodeId],
    point: &[NodeId],
) -> NodeId {
    vars.iter().zip(point).fold(term, |acc, (&v, &p)| graph.substitute(acc, v, p))
}

/// The value of a closed numeric term (parameters get generic values).
fn number(
    cx: &mut Cx<'_>,
    term: NodeId,
) -> Option<f64> {
    let term = cx.simplify(term);
    let mut env = Env::numeric(0.0);
    for &s in cx.graph.free_symbols(cx.graph.find(term)) {
        env.bind(s, 1.0 + 0.1 * f64::from(s.raw() % 7));
    }
    cx.graph.eval(term, &env).filter(|v| v.is_finite())
}

fn symbol(
    graph: &mut Graph,
    name: &str,
) -> NodeId {
    graph.sym(name)
}

/// Leading principal minors of a numeric symmetric matrix.
fn leading_minors(h: &[Vec<f64>]) -> Vec<f64> {
    (1..=h.len())
        .map(|k| {
            let sub: Vec<Vec<f64>> = h[..k].iter().map(|r| r[..k].to_vec()).collect();
            determinant(sub)
        })
        .collect()
}

/// Determinant by Gaussian elimination with partial pivoting.
fn determinant(mut m: Vec<Vec<f64>>) -> f64 {
    let n = m.len();
    let mut det = 1.0;
    for c in 0..n {
        let Some(p) = (c..n).max_by(|&a, &b| m[a][c].abs().total_cmp(&m[b][c].abs())) else {
            return 0.0;
        };
        if m[p][c].abs() < 1e-300 {
            return 0.0;
        }
        if p != c {
            m.swap(p, c);
            det = -det;
        }
        det *= m[c][c];
        for r in c + 1..n {
            let f = m[r][c] / m[c][c];
            for k in c..n {
                m[r][k] -= f * m[c][k];
            }
        }
    }
    det
}

fn unconstrained(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, vars] = args else {
        return None;
    };
    let f = best(cx.graph, f)?;
    let vars = list_items(cx.graph, vars)?;
    if vars.iter().any(|&v| cx.graph.symbol_of(v).is_none()) {
        return None;
    }
    let mut gradient = Vec::with_capacity(vars.len());
    for &v in &vars {
        let d = derivative(cx.graph, f, v)?;
        gradient.push(cx.simplify(d));
    }
    let points = solve_system(cx, &gradient, &vars)?;
    // Second derivatives.
    let mut hessian = vec![vec![NodeId::NONE; vars.len()]; vars.len()];
    for (i, &gi) in gradient.iter().enumerate() {
        for (j, &v) in vars.iter().enumerate() {
            let d = derivative(cx.graph, gi, v)?;
            hessian[i][j] = cx.simplify(d);
        }
    }
    let mut out = Vec::with_capacity(points.len());
    for point in points {
        let kind = if vars.len() == 1 {
            one_variable_kind(cx, f, vars[0], point[0])
        } else {
            let mut h = Vec::with_capacity(vars.len());
            for row in &hessian {
                let mut values = Vec::with_capacity(vars.len());
                for &entry in row {
                    let at = substitute_point(cx.graph, entry, &vars, &point);
                    values.push(number(cx, at)?);
                }
                h.push(values);
            }
            let minors = leading_minors(&h);
            let tiny = |v: f64| v.abs() < 1e-12;
            let last = *minors.last()?;
            if minors.iter().all(|&m| m > 0.0 && !tiny(m)) {
                "local_min"
            } else if minors.iter().enumerate().all(|(k, &m)| !tiny(m) && (m < 0.0) == (k % 2 == 0)) {
                "local_max"
            } else if !tiny(last) {
                "saddle"
            } else {
                "degenerate"
            }
        };
        let tuple = cx.graph.node(core::LIST, &point);
        let value = substitute_point(cx.graph, f, &vars, &point);
        let value = cx.simplify(value);
        let kind = symbol(cx.graph, kind);
        out.push(cx.graph.node(core::LIST, &[tuple, kind, value]));
    }
    Some(cx.graph.node(core::LIST, &out))
}

/// The kind of a critical point of `f(x)` from the first non-vanishing
/// derivative of order two or more.
fn one_variable_kind(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    point: NodeId,
) -> &'static str {
    let Some(mut d) = derivative(cx.graph, f, x) else {
        return "degenerate";
    };
    for order in 2..=8 {
        let Some(next) = derivative(cx.graph, d, x) else {
            return "degenerate";
        };
        d = cx.simplify(next);
        let at = cx.graph.substitute(d, x, point);
        let Some(value) = number(cx, at) else {
            return "degenerate";
        };
        if value.abs() > 1e-12 {
            return match (order % 2 == 0, value > 0.0) {
                | (true, true) => "local_min",
                | (true, false) => "local_max",
                | (false, _) => "saddle",
            };
        }
    }
    "degenerate"
}

fn constrained(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, constraints, vars] = args else {
        return None;
    };
    let f = best(cx.graph, f)?;
    let vars = list_items(cx.graph, vars)?;
    let constraints: Vec<NodeId> = list_items(cx.graph, constraints)?
        .into_iter()
        .map(|g| match *cx.graph.children(g) {
            | [lhs, rhs] if cx.graph.op(g) == core::EQ => sub(cx.graph, lhs, rhs),
            | _ => g,
        })
        .collect();
    if vars.iter().any(|&v| cx.graph.symbol_of(v).is_none()) || constraints.is_empty() {
        return None;
    }
    let multipliers: Vec<NodeId> = (0..constraints.len())
        .map(|k| {
            let s = cx.graph.interner_mut().fresh_symbol(&format!("lambda{}", k + 1));
            cx.graph.symbol_node(s)
        })
        .collect();
    // L = f - Σ λ_i g_i
    let mut terms = vec![f];
    for (&l, &g) in multipliers.iter().zip(&constraints) {
        let lg = mul(cx.graph, &[l, g]);
        let minus_one = cx.graph.int(-1);
        terms.push(mul(cx.graph, &[minus_one, lg]));
    }
    let lagrangian = cx.graph.node(core::ADD, &terms);
    let mut all = vars.clone();
    all.extend(&multipliers);
    let mut equations = Vec::with_capacity(all.len());
    for &v in &vars {
        let d = derivative(cx.graph, lagrangian, v)?;
        equations.push(cx.simplify(d));
    }
    equations.extend(constraints.iter().copied());
    let solutions = solve_system(cx, &equations, &all)?;
    let mut out = Vec::with_capacity(solutions.len());
    for solution in solutions {
        let (point, lambdas) = solution.split_at(vars.len());
        let kind = bordered_kind(cx, lagrangian, &constraints, &all, &solution).unwrap_or("critical");
        let value = substitute_point(cx.graph, f, &vars, point);
        let value = cx.simplify(value);
        let tuple = cx.graph.node(core::LIST, point);
        let lambdas = cx.graph.node(core::LIST, lambdas);
        let kind = symbol(cx.graph, kind);
        out.push(cx.graph.node(core::LIST, &[tuple, lambdas, kind, value]));
    }
    Some(cx.graph.node(core::LIST, &out))
}

/// Two variables, one constraint: the sign of the bordered Hessian.
fn bordered_kind(
    cx: &mut Cx<'_>,
    lagrangian: NodeId,
    constraints: &[NodeId],
    all: &[NodeId],
    solution: &[NodeId],
) -> Option<&'static str> {
    let &[g] = constraints else {
        return None;
    };
    let &[x, y, _] = all else {
        return None;
    };
    let at = |cx: &mut Cx<'_>, term: NodeId| -> Option<f64> {
        let point = substitute_point(cx.graph, term, all, solution);
        number(cx, point)
    };
    let gx = derivative(cx.graph, g, x)?;
    let gy = derivative(cx.graph, g, y)?;
    let lx = derivative(cx.graph, lagrangian, x)?;
    let ly = derivative(cx.graph, lagrangian, y)?;
    let lxx = derivative(cx.graph, lx, x)?;
    let lxy = derivative(cx.graph, lx, y)?;
    let lyy = derivative(cx.graph, ly, y)?;
    let m = vec![
        vec![0.0, at(cx, gx)?, at(cx, gy)?],
        vec![at(cx, gx)?, at(cx, lxx)?, at(cx, lxy)?],
        vec![at(cx, gy)?, at(cx, lxy)?, at(cx, lyy)?],
    ];
    let det = determinant(m);
    if det.abs() < 1e-12 {
        return None;
    }
    Some(if det > 0.0 { "local_max" } else { "local_min" })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[optimize()], src)
    }

    #[test]
    fn critical_points() {
        assert_eq!(run("find_extrema(x^2 - 4*x + 1, list(x))"), "list(list(list(2), local_min, -3))");
        assert_eq!(run("find_extrema(x^3, list(x))"), "list(list(list(0), saddle, 0))");
        assert_eq!(run("find_extrema(-x^4, list(x))"), "list(list(list(0), local_max, 0))");
        assert_eq!(run("find_extrema(x^2 + y^2 - 2*x, list(x, y))"), "list(list(list(1, 0), local_min, -1))");
        assert_eq!(run("find_extrema(x^2 - y^2, list(x, y))"), "list(list(list(0, 0), saddle, 0))");
        let two = run("find_extrema(x^3 - 3*x, list(x))");
        assert!(two.contains("list(list(1), local_min, -2)") && two.contains("list(list(-1), local_max, 2)"), "{two}");
    }

    #[test]
    fn lagrange_multipliers() {
        // Extremes of x + y on the unit circle.
        let text = run("find_constrained_extrema(x + y, list(x^2 + y^2 - 1), list(x, y))");
        assert!(text.contains("local_max, 2^(1/2)") && text.contains("local_min, -2^(1/2)"), "{text}");
        // The closest point of the line x + y = 2 to the origin.
        let text = run("find_constrained_extrema(x^2 + y^2, list(x + y = 2), list(x, y))");
        assert_eq!(text, "list(list(list(1, 1), list(2), local_min, 2))");
    }
}
