//! D-finite functions: differential equations, recurrences and a decision
//! procedure for identities.
//!
//! The rewriting engine proves an identity by finding a chain of rules from
//! one side to the other. For the large class of *D-finite* expressions
//! (rational functions, `exp`, `ln`, `sqrt` and rational powers, the circular
//! and hyperbolic functions and their inverses, closed under `+`, `*`, powers
//! and derivatives) there is something better than a search: an algorithm
//! that always terminates with a proof of "equal" or "different"
//! ([`crate::kernels::holonomic`]). This rule set exposes it.
//!
//! | operator | value |
//! |---|---|
//! | `holonomic_ode(f, x)` | `list(p0, p1, ..., pr)`: polynomials with `p0 f + p1 f' + ... + pr f^(r) = 0` |
//! | `holonomic_recurrence(f, x, n)` | `list(list(s, c), ...)`: the Taylor coefficients `a_n` of `f` satisfy `sum c(n) a_(n+s) = 0` |
//! | `holonomic_series(f, x, k)` | the Taylor polynomial of `f` at `0` with the terms below `x^k` |
//! | `holonomic_zero(f, x)` | `true` or `false`: whether `f` is identically zero, decided |
//! | `holonomic_equal(f, g, x)` | whether `f` and `g` are the same function, decided |
//!
//! The expression must be univariate in `x` with rational constants, and each
//! function must be analytic at `0` with rational Taylor coefficients there:
//! the argument of `exp`, `sin`, `cos`, `sinh`, `cosh`, `atan`, `atanh` and
//! `asin` vanishes at `0` (and is a polynomial for the four circular and
//! hyperbolic functions), the argument of `ln` and the base of a fractional
//! power equal `1` at `0`. Anything else is left unreduced, as is a
//! computation that exceeds [`crate::kernels::holonomic::Limits`].

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::holonomic as hol;
use crate::kernels::holonomic::Atom;
use crate::kernels::holonomic::Expr;
use crate::kernels::holonomic::Limits;
use crate::kernels::holonomic::Poly;
use crate::kernels::holonomic::RatFn;
use crate::rules::linalg;
use crate::rules::poly::best;

use super::elementary::elementary;
use super::logic::logic;

/// The D-finite rule set.
#[must_use]
pub fn holonomic() -> RuleSet {
    RuleSet::new("holonomic", install)
        .needs(elementary())
        .needs(logic())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Kind {
    Ode,
    Recurrence,
    Series,
    Zero,
    Equal,
}

struct Holonomic {
    op: OpId,
    kind: Kind,
}

/// A purely rational sub-term as a rational function of `x`.
fn to_rat(
    graph: &Graph,
    node: NodeId,
    x: SymbolId,
) -> Option<RatFn> {
    if let Some(number) = graph.as_number(node) {
        return number.to_rational().map(RatFn::constant);
    }
    if let Some(symbol) = graph.as_symbol(node) {
        return (symbol == x).then(|| RatFn::poly(vec![BigRational::zero(), BigRational::one()]));
    }
    let op = graph.op(node);
    let children = graph.children(node);
    if op == core::ADD {
        return children.iter().try_fold(RatFn::zero(), |acc, &c| {
            Some(acc.add(&to_rat(graph, c, x)?))
        });
    }
    if op == core::MUL {
        return children
            .iter()
            .try_fold(RatFn::constant(BigRational::one()), |acc, &c| {
                Some(acc.mul(&to_rat(graph, c, x)?))
            });
    }
    if op == core::POW {
        let [base, exponent] = children else {
            return None;
        };
        let e = integer(graph, *exponent)?;
        let base = to_rat(graph, *base, x)?;
        let one = RatFn::constant(BigRational::one());
        let positive = (0..e.unsigned_abs().min(4096)).fold(one.clone(), |acc, _| acc.mul(&base));
        if e.unsigned_abs() > 4096 {
            return None;
        }
        return if e < 0 {
            one.div(&positive)
        } else {
            Some(positive)
        };
    }
    None
}

fn integer(
    graph: &Graph,
    node: NodeId,
) -> Option<i64> {
    let r = graph.as_number(node)?.to_rational()?;
    r.is_integer().then(|| r.to_integer().to_i64()).flatten()
}

fn polynomial(r: RatFn) -> Option<Poly> {
    (r.den.len() == 1).then_some(r.num)
}

/// The term as a D-finite expression, or `None` when it is outside the
/// class this rule set handles.
fn to_expr(
    graph: &Graph,
    node: NodeId,
    x: SymbolId,
) -> Option<Expr> {
    if let Some(r) = to_rat(graph, node, x) {
        return Some(Expr::Atom(Atom::Rational(r)));
    }
    let op = graph.op(node);
    let children = graph.children(node);
    if op == core::ADD {
        let mut parts = children.iter().map(|&c| to_expr(graph, c, x));
        let first = parts.next()??;
        return parts.try_fold(first, |acc, part| {
            Some(Expr::Add(Box::new(acc), Box::new(part?)))
        });
    }
    if op == core::MUL {
        let mut parts = children.iter().map(|&c| to_expr(graph, c, x));
        let first = parts.next()??;
        return parts.try_fold(first, |acc, part| {
            Some(Expr::Mul(Box::new(acc), Box::new(part?)))
        });
    }
    if op == core::POW {
        let [base, exponent] = children else {
            return None;
        };
        let e = graph.as_number(*exponent)?.to_rational()?;
        if e.is_integer()
            && let Some(n) = e.to_integer().to_u32()
        {
            return Some(Expr::Pow(Box::new(to_expr(graph, *base, x)?), n));
        }
        // A fractional or negative power of a rational function.
        return Some(Expr::Atom(Atom::Power(to_rat(graph, *base, x)?, e)));
    }
    let [argument] = children else {
        return None;
    };
    let g = to_rat(graph, *argument, x)?;
    let half = BigRational::new(BigInt::one(), BigInt::from(2));
    let atom = match &*graph.ops().get(op).name {
        | "exp" => Atom::Exp(g),
        | "ln" => Atom::Log(g),
        | "sqrt" => Atom::Power(g, half),
        | "sin" => Atom::Sin(polynomial(g)?),
        | "cos" => Atom::Cos(polynomial(g)?),
        | "sinh" => Atom::Sinh(polynomial(g)?),
        | "cosh" => Atom::Cosh(polynomial(g)?),
        | "atan" => Atom::Atan(g),
        | "atanh" => Atom::Atanh(g),
        | "asin" => Atom::Asin(g),
        | _ => return None,
    };
    atom.is_valid().then_some(Expr::Atom(atom))
}

/// `sum c_k v^k` as a term in the variable node `v`.
fn poly_term(
    graph: &mut Graph,
    p: &[BigRational],
    v: NodeId,
) -> NodeId {
    let mut terms = Vec::new();
    for (k, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let c = linalg::rational_node(graph, c);
        if k == 0 {
            terms.push(c);
            continue;
        }
        let e = graph.num(Number::Int(BigInt::from(k)));
        let power = graph.node(core::POW, &[v, e]);
        terms.push(linalg::mul(graph, &[c, power]));
    }
    linalg::add(graph, &terms)
}

/// A truth value as the kernel's boolean literal (what the dynamic ABI
/// returns as `Value::Bool`).
#[allow(clippy::unnecessary_wraps)]
fn truth(
    graph: &mut Graph,
    value: bool,
) -> Option<NodeId> {
    Some(graph.lit(crate::graph::Payload::Bool(value)))
}

impl Holonomic {
    /// The expression of argument `f` in the variable argument `x`, with the
    /// variable's node.
    fn read(
        graph: &mut Graph,
        f: NodeId,
        x: NodeId,
    ) -> Option<(Expr, NodeId)> {
        let x = best(graph, x)?;
        let symbol = graph.as_symbol(x)?;
        let f = best(graph, f)?;
        Some((to_expr(graph, f, symbol)?, x))
    }

    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let limits = Limits::default();
        let graph = &mut *cx.graph;
        match (self.kind, args) {
            | (Kind::Ode, [f, x]) => {
                let (expr, x) = Self::read(graph, *f, *x)?;
                let op = hol::annihilator(&expr, &limits)?;
                let items: Vec<NodeId> = op.iter().map(|p| poly_term(graph, p, x)).collect();
                Some(linalg::list(graph, &items))
            },
            | (Kind::Recurrence, [f, x, n]) => {
                let (expr, _) = Self::read(graph, *f, *x)?;
                let n = best(graph, *n)?;
                graph.as_symbol(n)?;
                let op = hol::annihilator(&expr, &limits)?;
                let items: Vec<NodeId> = hol::recurrence(&op)
                    .iter()
                    .map(|(shift, c)| {
                        let s = graph.num(Number::Int(BigInt::from(*shift)));
                        let c = poly_term(graph, c, n);
                        linalg::list(graph, &[s, c])
                    })
                    .collect();
                Some(linalg::list(graph, &items))
            },
            | (Kind::Series, [f, x, k]) => {
                let (expr, x) = Self::read(graph, *f, *x)?;
                let k = best(graph, *k)?;
                let k = usize::try_from(integer(graph, k)?).ok()?;
                if k > limits.max_terms {
                    return None;
                }
                let series = expr.series(k)?;
                Some(poly_term(graph, &series, x))
            },
            | (Kind::Zero, [f, x]) => {
                let (expr, _) = Self::read(graph, *f, *x)?;
                truth(graph, hol::is_zero(&expr, &limits)?)
            },
            | (Kind::Equal, [f, g, x]) => {
                let (left, _) = Self::read(graph, *f, *x)?;
                let (right, _) = Self::read(graph, *g, *x)?;
                let minus = BigRational::from(BigInt::from(-1));
                let difference = Expr::Add(
                    Box::new(left),
                    Box::new(Expr::Scale(minus, Box::new(right))),
                );
                truth(graph, hol::is_zero(&difference, &limits)?)
            },
            | _ => None,
        }
    }
}

impl Kernel for Holonomic {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, &args)
            .map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, kind, rule) in [
        ("holonomic_ode", 2, Kind::Ode, DimRule::Container),
        (
            "holonomic_recurrence",
            3,
            Kind::Recurrence,
            DimRule::Container,
        ),
        ("holonomic_series", 3, Kind::Series, DimRule::Same(0)),
        ("holonomic_zero", 2, Kind::Zero, DimRule::NonMetric),
        ("holonomic_equal", 3, Kind::Equal, DimRule::NonMetric),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        i.dim(op, rule);
        i.kernel(
            &format!("holonomic/{name}"),
            Tier::Reduce,
            Holonomic { op, kind },
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::holonomic;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[holonomic()], src)
    }

    #[test]
    fn identities_are_decided_not_searched_for() {
        for (src, want) in [
            ("holonomic_zero(sin(x)^2 + cos(x)^2 - 1, x)", "true"),
            ("holonomic_zero(sin(x)^2 + cos(x)^2 - 2, x)", "false"),
            ("holonomic_equal(sin(2*x), 2*sin(x)*cos(x), x)", "true"),
            ("holonomic_equal(cosh(x)^2 - sinh(x)^2, 1, x)", "true"),
            ("holonomic_equal(exp(x)*exp(-x), 1, x)", "true"),
            ("holonomic_equal(sqrt(1 + x)^2, 1 + x, x)", "true"),
            (
                "holonomic_equal((1 - x)^(-1/2)*(1 - x)^(-1/2), 1/(1 - x), x)",
                "true",
            ),
            (
                "holonomic_equal(cos(x)^3, (3*cos(x) + cos(3*x))/4, x)",
                "true",
            ),
            ("holonomic_equal(atan(x), x - x^3/3, x)", "false"),
            (
                "holonomic_equal(ln(1 + x) + ln(1 - x), ln(1 - x^2), x)",
                "true",
            ),
            (
                "holonomic_equal(2*atanh(x), ln(1 + x) - ln(1 - x), x)",
                "true",
            ),
            (
                "holonomic_equal(sin(x + x^2), sin(x)*cos(x^2) + cos(x)*sin(x^2), x)",
                "true",
            ),
        ] {
            assert_eq!(run(src), want, "{src}");
        }
    }

    #[test]
    fn equations_series_and_recurrences() {
        assert_eq!(run("holonomic_ode(sin(x), x)"), "list(1, 0, 1)");
        // The equation every product of two solutions of f'' + f = 0 satisfies.
        assert_eq!(run("holonomic_ode(sin(x)*cos(x), x)"), "list(0, 4, 0, 1)");
        assert_eq!(
            run("holonomic_series(exp(x), x, 4)"),
            run("1 + x + x^2/2 + x^3/6")
        );
        assert_eq!(
            run("holonomic_series(1/(1 - x - x^2), x, 6)"),
            run("1 + x + 2*x^2 + 3*x^3 + 5*x^4 + 8*x^5")
        );
        assert_eq!(
            run("holonomic_series(asin(x), x, 6)"),
            run("x + x^3/6 + 3*x^5/40")
        );
        assert_eq!(
            run("holonomic_recurrence(exp(x), x, n)"),
            run("list(list(0, -1), list(1, n + 1))")
        );
    }

    #[test]
    fn expressions_outside_the_class_are_left_alone() {
        for src in [
            // A pole at the origin.
            "holonomic_zero(1/x, x)",
            "holonomic_series(1/x, x, 3)",
            // Irrational Taylor coefficients.
            "holonomic_zero(exp(x + 1), x)",
            "holonomic_zero(ln(2 + x), x)",
            // A second symbol.
            "holonomic_zero(sin(a*x), x)",
            // Not D-finite.
            "holonomic_zero(tan(x), x)",
            // The variable must be a symbol.
            "holonomic_zero(sin(x), 2)",
        ] {
            let (_, reduced) = reduce_with(&[holonomic()], src, &[]);
            assert!(!reduced, "{src}");
        }
    }
}
