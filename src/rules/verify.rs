//! Verification: checking a claimed result instead of computing one.
//!
//! Verification operators return `true` or `false` for supported claims.
//! A residual that simplifies to an exact number must equal zero. Approximate
//! residuals use a numerical tolerance. For an unresolved residual, numerical
//! sampling rejects the claim if a finite sample refutes it or if no finite
//! sample is available. Definite integrals use quadrature. Limit checks call
//! the shared calculus limit reducer.
//!
//! | operator | checks |
//! |---|---|
//! | `verify_solution(equations, list(x = a, ...))` | the substitution satisfies the equation (or list of equations) |
//! | `verify_integral(f, x, F)` | `F' = f` |
//! | `verify_definite_integral(f, x, a, b, value)` | `∫_a^b f dx = value`, by quadrature |
//! | `verify_ode_solution(equation, y(x), solution)` | the solution satisfies the differential equation |
//! | `verify_inverse(A, B)` | `A B = I` |
//! | `verify_derivative(f, x, g)` | `f' = g` |
//! | `verify_limit(f, x, a, L, [dir])` | the shared calculus reducer finds `f → L` as `x → a`, optionally from `plus` or `minus` |

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernel::Signature;
use crate::kernel::TypeExpr;
use crate::kernel::infer::term;
use crate::kernels::integrate::gauss_kronrod_any;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;

use super::calculus::calculus;
use super::calculus::derivative;
use super::calculus::limit_for_claim;
use super::linalg::linalg;
use super::logic::logic;

/// The verification rule set.
#[must_use]
pub fn verify() -> RuleSet {
    RuleSet::new("verify", install)
        .needs(calculus())
        .needs(linalg())
        .needs(logic())
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
    for (name, arity, request, doc) in [
        (
            "verify_solution",
            2,
            Request::Solution,
            "`verify_solution(eq_or_expr, subst)`: checks whether substituting `subst` into `eq_or_expr` satisfies the equation or expression (returns `true` or `false`).",
        ),
        (
            "verify_integral",
            3,
            Request::Integral,
            "`verify_integral(f, x, F)` checks whether `diff(F, x) == f`.",
        ),
        (
            "verify_definite_integral",
            5,
            Request::DefiniteIntegral,
            "`verify_definite_integral(integrand, var, a, b, value)`: numerically and symbolically checks whether the definite integral of `integrand` over `[a, b]` equals `value`.",
        ),
        (
            "verify_ode_solution",
            3,
            Request::Ode,
            "`verify_ode_solution(equation, y(x), solution)` checks whether the proposed `solution` for `y(x)` satisfies the differential equation.",
        ),
        (
            "verify_inverse",
            2,
            Request::Inverse,
            "`verify_inverse(A, B)` checks whether the matrix product `A*B` equals the identity.",
        ),
        (
            "verify_derivative",
            3,
            Request::Derivative,
            "`verify_derivative(f, x, g)` checks whether `diff(f, x) == g`. An exact residual must be zero. Approximate residuals use the numerical tolerance.",
        ),
        (
            "verify_limit",
            4,
            Request::Limit,
            "`verify_limit(f, x, a, L, [dir])` checks the limit of `f` as `x` approaches `a` against `L`, which may be finite or signed infinity. The optional `plus` or `minus` selects one side. Uses the same calculus reducer as `limit`. A claim with insufficient evidence returns `false`.",
        ),
    ] {
        let descriptor = if request == Request::Limit {
            let mut parameters = vec![term(), TypeExpr::Symbol, term(), term()];
            let two_sided = Signature::new(parameters.clone(), TypeExpr::Bool);
            parameters.push(TypeExpr::Symbol);
            OpDescriptor::new(name, Arity::Variadic)
                .flags(OpFlags::HEAVY.with(OpFlags::PREDICATE))
                .binder(1, 0b1)
                .signature(two_sided)
                .signature(Signature::new(parameters, TypeExpr::Bool))
        } else {
            OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY)
        };
        let op = i.op(descriptor.cost(100).doc(doc))?;
        i.kernel(
            &format!("verify/{name}"),
            Tier::Reduce,
            Verify { op, request },
        );
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
        Outcome::Equal(cx.graph.lit(Payload::Bool(verdict)))
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

/// Reachable classes of a term from which [`vanishes`] evaluates its sample points concurrently.
const VERIFY_PAR_NODES: usize = 40;

/// Classifies finite numerical samples and ignores unavailable evaluations.
fn sample_verdict(values: impl IntoIterator<Item = Option<f64>>) -> Option<bool> {
    let mut evidence = false;
    for value in values.into_iter().flatten() {
        if value.abs() > 1e-8 {
            return Some(false);
        }
        evidence = true;
    }
    evidence.then_some(true)
}

/// Checks an exact residual directly and samples an unresolved expression.
fn vanishes(
    cx: &mut Cx<'_>,
    term: NodeId,
) -> Option<bool> {
    let simplified = cx.simplify(term);
    if let Some(n) = cx.graph.number_of(simplified) {
        return Some(n.is_zero() || (!n.is_exact() && n.to_f64().abs() < 1e-12));
    }
    let symbols = cx.graph.free_symbols(cx.graph.find(simplified)).to_vec();
    // The twelve sample points are independent evaluations of an immutable graph. The verdict
    // (some sample refutes / some sample is finite) does not depend on their order, so a large
    // term is evaluated at the samples concurrently and a small one sequentially.
    let graph: &Graph = cx.graph;
    let sample_value = |sample: u32| -> Option<f64> {
        let mut env = Env::numeric(0.0);
        for &s in &symbols {
            let k = f64::from(s.raw() % 13);
            env.bind(
                s,
                0.31 + 0.47 * f64::from(sample) - 0.09 * k + 0.013 * k * k,
            );
        }
        graph.eval(simplified, &env).filter(|v| v.is_finite())
    };
    if crate::graph::extract::reachable(graph, &[simplified]).len() >= VERIFY_PAR_NODES
        && crate::execution::should_split(12, 2)
    {
        sample_verdict(crate::kernels::par::map_indexed(12, 2, |s| {
            sample_value(s as u32)
        }))
    } else {
        sample_verdict((0..12_u32).map(sample_value))
    }
}

fn solution(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<bool> {
    let &[equations, assignment] = args else {
        return None;
    };
    let equations = best(cx.graph, equations)?;
    let equations = if cx.graph.op(equations) == core::LIST {
        cx.graph.children(equations).to_vec()
    } else {
        vec![equations]
    };
    let assignment = best(cx.graph, assignment)?;
    let pairs = if cx.graph.op(assignment) == core::LIST {
        cx.graph.children(assignment).to_vec()
    } else {
        vec![assignment]
    };
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
        && cx.graph.op(solution) == core::EQ
        && cx.graph.same(lhs, unknown)
    {
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
    let (f, x, a, claim, direction) = match *args {
        | [f, x, a, claim] => (f, x, a, claim, None),
        | [f, x, a, claim, direction] => (f, x, a, claim, Some(direction)),
        | _ => return None,
    };
    let Some(found) = limit_for_claim(cx, f, x, a, direction) else {
        return Some(false);
    };
    let claim = cx.simplify(claim);
    if cx.graph.same(found, claim) {
        return Some(true);
    }
    let infinity = cx.graph.ops().lookup("oo")?;
    let positive = cx.graph.node(infinity, &[]);
    let minus_one = cx.graph.int(-1);
    let negative = cx.graph.node(core::MUL, &[minus_one, positive]);
    let sign = |value| {
        if cx.graph.same(value, positive) {
            Some(true)
        } else if cx.graph.same(value, negative) {
            Some(false)
        } else {
            match cx.graph.number_of(value) {
                | Some(Number::Float(n)) if n.is_infinite() => Some(n.is_sign_positive()),
                | _ => None,
            }
        }
    };
    let signs = (sign(found), sign(claim));
    Some(match signs {
        | (None, None) => {
            let difference = sub(cx.graph, found, claim);
            vanishes(cx, difference) == Some(true)
        },
        | (found, claim) => found == claim,
    })
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
        assert_eq!(
            run("verify_solution(x^2 - 5*x + 6 = 0, list(x = 2))"),
            "true"
        );
        assert_eq!(
            run("verify_solution(x^2 - 5*x + 6 = 0, list(x = 4))"),
            "false"
        );
        assert_eq!(
            run("verify_solution(list(x + y = 3, x - y = 1), list(x = 2, y = 1))"),
            "true"
        );
        assert_eq!(
            run("verify_inverse(list(list(2, 1), list(1, 1)), list(list(1, -1), list(-1, 2)))"),
            "true"
        );
        assert_eq!(
            run("verify_inverse(list(list(2, 1), list(1, 1)), list(list(1, 1), list(-1, 2)))"),
            "false"
        );
    }

    #[test]
    fn calculus_claims() {
        assert_eq!(
            run("verify_integral(x*cos(x), x, x*sin(x) + cos(x))"),
            "true"
        );
        assert_eq!(run("verify_integral(x*cos(x), x, x*sin(x))"), "false");
        assert_eq!(run("verify_derivative(sin(x)^2, x, sin(2*x))"), "true");
        assert_eq!(run("verify_derivative(sin(x)^2, x, cos(x)^2)"), "false");
        assert_eq!(
            run("verify_definite_integral(exp(-x^2), x, -oo, oo, pi^(1/2))"),
            "true"
        );
        assert_eq!(run("verify_definite_integral(x, x, 0, 1, 1)"), "false");
        assert_eq!(run("verify_limit(sin(x)/x, x, 0, 1)"), "true");
        assert_eq!(run("verify_limit(sin(x)/x, x, 0, 2)"), "false");
        assert_eq!(run("verify_limit((1 + 1/x)^x, x, oo, E)"), "true");
        assert_eq!(
            run("verify_ode_solution(diff(diff(y(x), x), x) + y(x) = 0, y(x), 3*cos(x) - sin(x))"),
            "true"
        );
        assert_eq!(
            run("verify_ode_solution(diff(y(x), x) = y(x), y(x), y(x) = 2*exp(x))"),
            "true"
        );
        assert_eq!(
            run("verify_ode_solution(diff(y(x), x) = y(x), y(x), exp(2*x))"),
            "false"
        );
    }
}
