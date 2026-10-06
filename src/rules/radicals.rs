//! Radicals: denesting and rationalising denominators.
//!
//! * `√(a + b√c)` with rational `a, b, c` denests when `d = a² - b²c` is
//!   the square of a rational: `√(a + b√c) = √((a + √d)/2) + sign(b)
//!   √((a - √d)/2)` (a kernel, so it happens automatically; also the
//!   operator `denest_sqrt(e)`).
//! * `simplify_radicals(e)`: denesting everywhere, and denominators
//!   `p + q√c` multiplied by their conjugate.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
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
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

/// The radicals rule set.
#[must_use]
pub fn radicals() -> RuleSet {
    RuleSet::new("radicals", install)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let identity = |a: &[f64]| a.first().copied().unwrap_or(f64::NAN);
    let denest = i.op(OpDescriptor::new("denest_sqrt", Arity::Fixed(1)).flags(OpFlags::HEAVY).cost(100).eval(identity))?;
    let simplify = i.op(OpDescriptor::new("simplify_radicals", Arity::Fixed(1)).flags(OpFlags::HEAVY).cost(100).eval(identity))?;
    i.kernel("radicals/requests", Tier::Reduce, Requests { denest, simplify });
    i.kernel("radicals/denest", Tier::Normalize, Denest);
    Ok(())
}

fn rational_sqrt(r: &BigRational) -> Option<BigRational> {
    if r.is_negative() {
        return None;
    }
    let (n, d) = (r.numer().sqrt(), r.denom().sqrt());
    (&n * &n == *r.numer() && &d * &d == *r.denom()).then(|| BigRational::new(n, d))
}

/// `(p, q, c)` with `e = p + q √c`, all rational, `c` positive and not a
/// rational square.
fn split_quadratic(
    graph: &Graph,
    e: NodeId,
) -> Option<(BigRational, BigRational, BigRational)> {
    let half = BigRational::new(BigInt::from(1), BigInt::from(2));
    let surd = |n: NodeId| -> Option<(BigRational, BigRational)> {
        // q √c  or  √c
        let factors: Vec<NodeId> = if graph.op(n) == core::MUL { graph.children(n).to_vec() } else { vec![n] };
        let mut q = BigRational::from_integer(BigInt::from(1));
        let mut c = None;
        for f in factors {
            if let Some(v) = graph.number_of(f).and_then(Number::to_rational) {
                q *= v;
            } else if let (true, &[base, exp]) = (graph.op(f) == core::POW, graph.children(f)) {
                let e = graph.number_of(exp)?.to_rational()?;
                if e != half || c.is_some() {
                    return None;
                }
                c = Some(graph.number_of(base)?.to_rational()?);
            } else {
                return None;
            }
        }
        Some((q, c?))
    };
    let terms: Vec<NodeId> = if graph.op(e) == core::ADD { graph.children(e).to_vec() } else { vec![e] };
    let mut p = BigRational::zero();
    let mut found = None;
    for t in terms {
        if let Some(v) = graph.number_of(t).and_then(Number::to_rational) {
            p += v;
        } else {
            let (q, c) = surd(t)?;
            match &mut found {
                | None => found = Some((q, c)),
                | Some((q0, c0)) if *c0 == c => *q0 += q,
                | Some(_) => return None,
            }
        }
    }
    let (q, c) = found?;
    (c.is_positive() && rational_sqrt(&c).is_none()).then_some((p, q, c))
}

fn rat(
    graph: &mut Graph,
    r: BigRational,
) -> NodeId {
    graph.num(Number::rat(r))
}

fn sqrt_of(
    graph: &mut Graph,
    r: BigRational,
) -> Option<NodeId> {
    if let Some(s) = rational_sqrt(&r) {
        return Some(rat(graph, s));
    }
    let base = rat(graph, r);
    let half = graph.num(Number::fraction(1, 2)?);
    Some(graph.node(core::POW, &[base, half]))
}

/// `√(p + q √c)` denested, if possible.
fn denest(
    graph: &mut Graph,
    radicand: NodeId,
) -> Option<NodeId> {
    let (p, q, c) = split_quadratic(graph, radicand)?;
    let d = &p * &p - &q * &q * &c;
    let root_d = rational_sqrt(&d)?;
    let two = BigRational::from_integer(BigInt::from(2));
    let (u, v) = ((&p + &root_d) / &two, (&p - &root_d) / &two);
    if u.is_negative() || v.is_negative() {
        return None;
    }
    let (su, sv) = (sqrt_of(graph, u)?, sqrt_of(graph, v)?);
    let sv = if q.is_negative() {
        let minus = graph.int(-1);
        graph.node(core::MUL, &[minus, sv])
    } else {
        sv
    };
    Some(graph.node(core::ADD, &[su, sv]))
}

/// Applies `f` bottom-up to every node of `e`.
fn rewrite(
    graph: &mut Graph,
    e: NodeId,
    f: &dyn Fn(&mut Graph, NodeId) -> Option<NodeId>,
) -> NodeId {
    let children = graph.children(e).to_vec();
    let op = graph.op(e);
    let new: Vec<NodeId> = children.iter().map(|&c| rewrite(graph, c, f)).collect();
    let rebuilt = if new == children { e } else { graph.try_node(op, &new).unwrap_or(e) };
    f(graph, rebuilt).unwrap_or(rebuilt)
}

fn as_sqrt(
    graph: &Graph,
    n: NodeId,
) -> Option<NodeId> {
    let sqrt = graph.ops().lookup("sqrt");
    match graph.children(n) {
        | &[base] if Some(graph.op(n)) == sqrt => Some(base),
        | &[base, exp]
            if graph.op(n) == core::POW
                && graph.number_of(exp).and_then(Number::to_rational) == Some(BigRational::new(BigInt::from(1), BigInt::from(2))) =>
        {
            Some(base)
        },
        | _ => None,
    }
}

fn denest_node(
    graph: &mut Graph,
    n: NodeId,
) -> Option<NodeId> {
    let base = as_sqrt(graph, n)?;
    denest(graph, base)
}

/// `1/(p + q√c) = (p - q√c)/(p² - q²c)`.
fn rationalise(
    graph: &mut Graph,
    n: NodeId,
) -> Option<NodeId> {
    let &[base, exp] = graph.children(n) else {
        return None;
    };
    if graph.op(n) != core::POW || graph.number_of(exp).and_then(Number::to_i64) != Some(-1) {
        return None;
    }
    let (p, q, c) = split_quadratic(graph, base)?;
    let norm = &p * &p - &q * &q * &c;
    if norm.is_zero() {
        return None;
    }
    let p_n = rat(graph, &p / &norm);
    let q_n = rat(graph, -(&q / &norm));
    let root = sqrt_of(graph, c)?;
    let surd = graph.node(core::MUL, &[q_n, root]);
    Some(graph.node(core::ADD, &[p_n, surd]))
}

struct Requests {
    denest: OpId,
    simplify: OpId,
}

impl Kernel for Requests {
    fn ops(&self) -> Vec<OpId> {
        vec![self.denest, self.simplify]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[e] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(e) = crate::rules::poly::best(cx.graph, e) else {
            return Outcome::Pass;
        };
        let op = cx.graph.op(node);
        let rewritten = if op == self.denest {
            rewrite(cx.graph, e, &denest_node)
        } else {
            let once = rewrite(cx.graph, e, &denest_node);
            rewrite(cx.graph, once, &rationalise)
        };
        let simplified = cx.simplify(rewritten);
        Outcome::Pinned(simplified)
    }
}

/// Automatic denesting of `√(a + b√c)`.
struct Denest;

impl Kernel for Denest {
    fn ops(&self) -> Vec<OpId> {
        vec![core::POW]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        denest_node(cx.graph, node).map_or(Outcome::Pass, Outcome::Equal)
    }
}

#[cfg(test)]
mod tests {
    use crate::rules::testing::eval;
    use crate::rules::testing::simplify;

    #[test]
    fn denesting_and_rationalising() {
        let rules = crate::rules::standard();
        let run = |src: &str| simplify(&rules, src);
        // The results are free of nested radicals and equal in value.
        let nested = |text: &str| text.matches("^(1/2)").count() > 0 && text.contains(")^(1/2)") && text.contains(" + ") && text.starts_with('(');
        for (src, want) in [
            ("denest_sqrt((3 + 2*2^(1/2))^(1/2))", 2.0_f64.sqrt() + 1.0),
            ("denest_sqrt((5 - 2*6^(1/2))^(1/2))", 3.0_f64.sqrt() - 2.0_f64.sqrt()),
        ] {
            let got = run(src);
            assert!(!nested(&got), "{got}");
            assert!((eval(&rules, &got, &[]) - want).abs() < 1e-12, "{src} = {got}");
        }
        // No denesting possible: unchanged value.
        let stay = run("denest_sqrt((2 + 3^(1/2)*5)^(1/2))");
        assert!((eval(&rules, &stay, &[]) - (2.0 + 5.0 * 3.0_f64.sqrt()).sqrt()).abs() < 1e-12);
        let rational = run("simplify_radicals(1/(1 + 2^(1/2)))");
        assert!(!rational.contains("/(") && (eval(&rules, &rational, &[]) - (2.0_f64.sqrt() - 1.0)).abs() < 1e-12, "{rational}");
        let both = run("simplify_radicals(1/(3 + 2*2^(1/2))^(1/2))");
        assert!((eval(&rules, &both, &[]) - (2.0_f64.sqrt() - 1.0)).abs() < 1e-12 && !both.contains("(3 +"), "{both}");
    }
}
