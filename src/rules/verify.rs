//! Verification: checking a claimed result instead of computing one.
//!
//! Every operator answers `true` or `false`. A claim is accepted when the
//! defining identity simplifies to zero, or — when simplification cannot
//! decide — when it holds numerically at sample points of every free
//! symbol; a single sample that refutes it gives `false`. A claim that can
//! be neither simplified nor evaluated anywhere is not accepted.
//!
//! | operator | checks |
//! |---|---|
//! | `verify_solution(equations, list(x = a, ...))` | the substitution satisfies the equation (or list of equations) |
//! | `verify_integral(f, x, F)` | `F' = f` |
//! | `verify_definite_integral(f, x, a, b, value)` | `∫_a^b f dx = value`, by quadrature |
//! | `verify_ode_solution(equation, y(x), solution)` | the solution satisfies the differential equation |
//! | `verify_inverse(A, B)` | `A B = I` |
//! | `verify_derivative(f, x, g)` | `f' = g` |
//! | `verify_limit(f, x, a, L)` | `f → L` as `x → a`, by probing near `a` |

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
use crate::kernels::integrate::gauss_kronrod_any;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;

use super::calculus::derivative;
use super::linalg::linalg;
use super::logic::logic;

/// The verification rule set.
#[must_use]
pub fn verify() -> RuleSet {
    RuleSet::new("verify", install).needs(linalg()).needs(logic())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Solution,
    Integral,
    DefiniteIntegral,
    Ode,
    Inverse,
    Derivative,
    Limit,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, request) in [
        ("verify_solution", 2, Request::Solution),
        ("verify_integral", 3, Request::Integral),
        ("verify_definite_integral", 5, Request::DefiniteIntegral),
        ("verify_ode_solution", 3, Request::Ode),
        ("verify_inverse", 2, Request::Inverse),
        ("verify_derivative", 3, Request::Derivative),
        ("verify_limit", 4, Request::Limit),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("verify/{name}"), Tier::Reduce, Verify { op, request });
    }
    Ok(())
}

struct Verify {
    op: OpId,
    request: Request,
}

impl Kernel for Verify {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let verdict = match self.request {
            | Request::Solution => solution(cx, &args),
            | Request::Integral => integral(cx, &args),
            | Request::DefiniteIntegral => definite(cx, &args),
            | Request::Ode => ode(cx, &args),
            | Request::Inverse => inverse(cx, &args),
            | Request::Derivative => derivative_claim(cx, &args),
            | Request::Limit => limit(cx, &args),
        };
        let Some(verdict) = verdict else {
            return Outcome::Pass;
        };
        let name = if verdict { "true" } else { "false" };
        match cx.graph.ops().lookup(name) {
            | Some(op) => Outcome::Equal(cx.graph.node(op, &[])),
            | None => Outcome::Pass,
        }
    }
}

/// `lhs - rhs` of an equation, or the term itself.
fn residual(
    graph: &mut Graph,
    term: NodeId,
) -> NodeId {
    match *graph.children(term) {
        | [lhs, rhs] if graph.op(term) == core::EQ => sub(graph, lhs, rhs),
        | _ => term,
    }
}

/// Whether `term` is identically zero: `Some(true)` when it simplifies to
/// zero or vanishes at every sample, `Some(false)` when a sample refutes
/// it, `None` when there is no evidence either way.
fn vanishes(
    cx: &mut Cx<'_>,
    term: NodeId,
) -> Option<bool> {
    let simplified = cx.simplify(term);
    if let Some(n) = cx.graph.number_of(simplified) {
        return Some(n.is_zero() || n.to_f64().abs() < 1e-12);
    }
    let symbols = cx.graph.free_symbols(cx.graph.find(simplified)).to_vec();
    let mut evidence = false;
    for sample in 0..12_u32 {
        let mut env = Env::numeric(0.0);
        for &s in &symbols {
            let k = f64::from(s.raw() % 13);
            env.bind(s, 0.31 + 0.47 * f64::from(sample) - 0.09 * k + 0.013 * k * k);
        }
        let Some(value) = cx.graph.eval(simplified, &env) else {
            continue;
        };
        if !value.is_finite() {
            continue;
        }
        if value.abs() > 1e-8 {
            return Some(false);
        }
        evidence = true;
    }
    evidence.then_some(true)
}

fn solution(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[equations, assignment] = args else {
        return None;
    };
    let equations = best(cx.graph, equations)?;
    let equations = if cx.graph.op(equations) == core::LIST { cx.graph.children(equations).to_vec() } else { vec![equations] };
    let assignment = best(cx.graph, assignment)?;
    let pairs = if cx.graph.op(assignment) == core::LIST { cx.graph.children(assignment).to_vec() } else { vec![assignment] };
    let mut substitutions = Vec::new();
    for pair in pairs {
        let &[x, value] = cx.graph.children(pair) else {
            return None;
        };
        if cx.graph.op(pair) != core::EQ {
            return None;
        }
        cx.graph.symbol_of(x)?;
        substitutions.push((x, value));
    }
    for equation in equations {
        let mut r = residual(cx.graph, equation);
        for &(x, value) in &substitutions {
            r = cx.graph.substitute(r, x, value);
        }
        if vanishes(cx, r) != Some(true) {
            return Some(false);
        }
    }
    Some(true)
}

fn integral(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[f, x, big_f] = args else {
        return None;
    };
    cx.graph.symbol_of(x)?;
    let big_f = best(cx.graph, big_f)?;
    let d = derivative(cx.graph, big_f, x)?;
    let difference = sub(cx.graph, d, f);
    Some(vanishes(cx, difference) == Some(true))
}

fn derivative_claim(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[f, x, g] = args else {
        return None;
    };
    cx.graph.symbol_of(x)?;
    let f = best(cx.graph, f)?;
    let d = derivative(cx.graph, f, x)?;
    let difference = sub(cx.graph, d, g);
    Some(vanishes(cx, difference) == Some(true))
}

fn definite(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[f, x, a, b, value] = args else {
        return None;
    };
    let symbol = cx.graph.symbol_of(x)?;
    let f = best(cx.graph, f)?;
    let env = Env::numeric(0.0);
    let (a, b) = (cx.graph.eval(a, &env)?, cx.graph.eval(b, &env)?);
    let claimed = cx.graph.eval(value, &env)?;
    let graph = &*cx.graph;
    let q = gauss_kronrod_any(
        |t| {
            let mut env = Env::numeric(0.0);
            env.bind(symbol, t);
            graph.eval(f, &env).unwrap_or(f64::NAN)
        },
        a,
        b,
        1e-11,
        2_000,
    );
    if !q.value.is_finite() {
        return None;
    }
    Some((q.value - claimed).abs() <= 1e-7 * claimed.abs().max(1.0) + q.error)
}

fn ode(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[equation, unknown, solution] = args else {
        return None;
    };
    let equation = best(cx.graph, equation)?;
    let unknown = best(cx.graph, unknown)?;
    if cx.graph.op(unknown) != core::APPLY {
        return None;
    }
    let mut solution = best(cx.graph, solution)?;
    // `y(x) = expression` is accepted as well as the bare expression.
    if let [lhs, rhs] = *cx.graph.children(solution)
        && cx.graph.op(solution) == core::EQ && cx.graph.same(lhs, unknown) {
            solution = rhs;
        }
    let r = residual(cx.graph, equation);
    let replaced = cx.graph.replace_subterm(r, unknown, solution);
    Some(vanishes(cx, replaced) == Some(true))
}

fn inverse(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[a, b] = args else {
        return None;
    };
    let matmul = cx.graph.ops().lookup("matmul")?;
    let product = cx.graph.node(matmul, &[a, b]);
    let product = cx.simplify(product);
    if cx.graph.op(product) != core::LIST {
        return None;
    }
    let rows = cx.graph.children(product).to_vec();
    for (i, row) in rows.iter().enumerate() {
        let entries = cx.graph.children(*row).to_vec();
        if entries.len() != rows.len() {
            return Some(false);
        }
        for (j, &e) in entries.iter().enumerate() {
            let expected = cx.graph.int(i64::from(i == j));
            let difference = sub(cx.graph, e, expected);
            if vanishes(cx, difference) != Some(true) {
                return Some(false);
            }
        }
    }
    Some(true)
}

fn limit(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[f, x, a, l] = args else {
        return None;
    };
    let symbol = cx.graph.symbol_of(x)?;
    let f = best(cx.graph, f)?;
    let env = Env::numeric(0.0);
    let point = cx.graph.eval(a, &env)?;
    let claimed = cx.graph.eval(l, &env)?;
    let probe = |graph: &Graph, t: f64| {
        let mut env = Env::numeric(0.0);
        env.bind(symbol, t);
        graph.eval(f, &env)
    };
    let mut evidence = false;
    for k in 3..=7 {
        let eps = 10f64.powi(-k);
        let points: Vec<f64> = if point.is_infinite() {
            vec![point.signum() / eps]
        } else {
            vec![point + eps, point - eps]
        };
        for t in points {
            let Some(v) = probe(cx.graph, t) else {
                continue;
            };
            if !v.is_finite() {
                if claimed.is_infinite() && v.signum() == claimed.signum() {
                    evidence = true;
                    continue;
                }
                continue;
            }
            if claimed.is_infinite() {
                if k == 7 && v.abs() < 1e6 {
                    return Some(false);
                }
                continue;
            }
            let tolerance = 1e-6_f64.max(100.0 * eps.sqrt()) * claimed.abs().max(1.0);
            if k >= 6 && (v - claimed).abs() > tolerance {
                return Some(false);
            }
            evidence = true;
        }
    }
    Some(evidence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[verify(), crate::rules::ode()], src)
    }

    #[test]
    fn solutions_and_inverses() {
        assert_eq!(run("verify_solution(x^2 - 5*x + 6 = 0, list(x = 2))"), "true");
        assert_eq!(run("verify_solution(x^2 - 5*x + 6 = 0, list(x = 4))"), "false");
        assert_eq!(run("verify_solution(list(x + y = 3, x - y = 1), list(x = 2, y = 1))"), "true");
        assert_eq!(run("verify_inverse(list(list(2, 1), list(1, 1)), list(list(1, -1), list(-1, 2)))"), "true");
        assert_eq!(run("verify_inverse(list(list(2, 1), list(1, 1)), list(list(1, 1), list(-1, 2)))"), "false");
    }

    #[test]
    fn calculus_claims() {
        assert_eq!(run("verify_integral(x*cos(x), x, x*sin(x) + cos(x))"), "true");
        assert_eq!(run("verify_integral(x*cos(x), x, x*sin(x))"), "false");
        assert_eq!(run("verify_derivative(sin(x)^2, x, sin(2*x))"), "true");
        assert_eq!(run("verify_derivative(sin(x)^2, x, cos(x)^2)"), "false");
        assert_eq!(run("verify_definite_integral(exp(-x^2), x, -oo, oo, pi^(1/2))"), "true");
        assert_eq!(run("verify_definite_integral(x, x, 0, 1, 1)"), "false");
        assert_eq!(run("verify_limit(sin(x)/x, x, 0, 1)"), "true");
        assert_eq!(run("verify_limit(sin(x)/x, x, 0, 2)"), "false");
        assert_eq!(run("verify_limit((1 + 1/x)^x, x, oo, E)"), "true");
        assert_eq!(run("verify_ode_solution(diff(diff(y(x), x), x) + y(x) = 0, y(x), 3*cos(x) - sin(x))"), "true");
        assert_eq!(run("verify_ode_solution(diff(y(x), x) = y(x), y(x), y(x) = 2*exp(x))"), "true");
        assert_eq!(run("verify_ode_solution(diff(y(x), x) = y(x), y(x), exp(2*x))"), "false");
    }
}
