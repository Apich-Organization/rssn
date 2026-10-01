//! Calculus of variations.
//!
//! A functional `S[u] = ∫ L(x, u, u', u'', ...) dx` is given by its
//! Lagrangian `L`, written with the undetermined function itself —
//! `y(x)`, `diff(y(x), x)`, `diff(diff(y(x), x), x)` — exactly as in
//! [`dsolve`](super::ode). Functions of several variables (`u(x, y)`,
//! `diff(u(x, y), y)`) give field equations.
//!
//! | operator | value |
//! |---|---|
//! | `euler_lagrange(L, y(x), x)` | the Euler–Lagrange (Euler–Poisson) expression `Σ (-1)^k d^k/dx^k ∂L/∂y^(k)`; its vanishing is the condition for an extremal |
//! | `euler_lagrange(L, list(q1(t), ...), t)` | one expression per function |
//! | `euler_lagrange(L, u(x, y), list(x, y))` | the field equation `∂L/∂u - Σ ∂_i ∂L/∂u_i + ...` |
//! | `solve_euler_lagrange(L, y(x), x)` | the extremals: `dsolve` applied to the Euler–Lagrange equation |
//! | `hamiltons_principle(L, q(t), t)` | the equations of motion, `euler_lagrange(L, q(t), t)` |
//! | `action(L, y(x), path, x, a, b)` | `∫_a^b L dx` along the given path |
//! | `first_integral(L, y(x), x)` | the Beltrami identity `L - y' ∂L/∂y'`, constant along extremals of an `L` free of `x` |
//!
//! The derivatives are taken by replacing every derivative of the
//! function by a fresh symbol, differentiating with respect to that
//! symbol and putting the derivative back, then applying the *total*
//! derivative — the structural differentiator treats `y(x)` as an
//! undetermined function, so the chain rule produces `y''` and higher
//! derivatives on its own.

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
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
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;

use super::calculus::derivative;
use super::ode::ode;

/// The calculus-of-variations rule set.
#[must_use]
pub fn variational() -> RuleSet {
    RuleSet::new("variational", install).needs(ode())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    EulerLagrange,
    Action,
    FirstIntegral,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let heavy = |name: &str, arity: u8| OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
    let diff = i.graph().ops().lookup("diff").ok_or(RuleError::Invalid {
        rule: "variational".to_owned(),
        reason: "needs the calculus rule set",
    })?;
    let defint = i.graph().ops().lookup("defint").ok_or(RuleError::Invalid {
        rule: "variational".to_owned(),
        reason: "needs the calculus rule set",
    })?;
    for (name, arity, request) in [
        ("euler_lagrange", 3, Request::EulerLagrange),
        ("action", 6, Request::Action),
        ("first_integral", 3, Request::FirstIntegral),
    ] {
        let op = i.op(heavy(name, arity))?;
        i.kernel(&format!("variational/{name}"), Tier::Reduce, Variational { op, request, diff, defint });
    }
    i.define(&[
        "solve_euler_lagrange(L, y, x) := dsolve(euler_lagrange(L, y, x) = 0, y)",
        "hamiltons_principle(L, q, t) := euler_lagrange(L, q, t)",
    ])
}

struct Variational {
    op: OpId,
    request: Request,
    diff: OpId,
    defint: OpId,
}

impl Kernel for Variational {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = match self.request {
            | Request::EulerLagrange => self.euler_lagrange(cx, &args),
            | Request::Action => self.action(cx, &args),
            | Request::FirstIntegral => self.first_integral(cx, &args),
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

/// A derivative `∂^α u` of the function: the term and the variables it
/// was differentiated by, in order.
struct Jet {
    term: NodeId,
    by: Vec<NodeId>,
}

impl Variational {
    /// Every derivative of `f` (including `f` itself) occurring in `term`,
    /// highest order first.
    fn jets(
        &self,
        graph: &Graph,
        term: NodeId,
        f: NodeId,
    ) -> Vec<Jet> {
        let mut found: Vec<Jet> = vec![Jet { term: f, by: Vec::new() }];
        let mut stack = vec![term];
        let mut seen = std::collections::HashSet::new();
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            if let Some(by) = self.chain(graph, n, f) {
                if !by.is_empty() && !found.iter().any(|j| j.term == n) {
                    found.push(Jet { term: n, by });
                }
                continue;
            }
            stack.extend(graph.children(n).iter().copied());
        }
        found.sort_by_key(|j| std::cmp::Reverse(j.by.len()));
        found
    }

    /// If `n` is `diff(...diff(f, x1)..., xk)`, the variables `x1..xk`.
    fn chain(
        &self,
        graph: &Graph,
        n: NodeId,
        f: NodeId,
    ) -> Option<Vec<NodeId>> {
        if n == f {
            return Some(Vec::new());
        }
        if graph.op(n) != self.diff {
            return None;
        }
        let &[inner, x] = graph.children(n) else {
            return None;
        };
        let mut by = self.chain(graph, inner, f)?;
        by.push(x);
        Some(by)
    }

    /// `∂L/∂J` for the jet `J`, all jets frozen as symbols.
    fn partial(
        cx: &mut Cx<'_>,
        l: NodeId,
        jets: &[Jet],
        which: usize,
    ) -> Option<NodeId> {
        let mut frozen = l;
        let mut symbols = Vec::with_capacity(jets.len());
        for jet in jets {
            let fresh = cx.graph.interner_mut().fresh_symbol("jet");
            let s = cx.graph.symbol_node(fresh);
            frozen = cx.graph.replace_subterm(frozen, jet.term, s);
            symbols.push(s);
        }
        let mut d = derivative(cx.graph, frozen, symbols[which])?;
        d = cx.simplify(d);
        for (k, jet) in jets.iter().enumerate().rev() {
            d = cx.graph.substitute(d, symbols[k], jet.term);
        }
        Some(d)
    }

    fn equation(
        &self,
        cx: &mut Cx<'_>,
        l: NodeId,
        f: NodeId,
    ) -> Option<NodeId> {
        if cx.graph.op(f) != core::APPLY && cx.graph.symbol_of(f).is_none() {
            return None;
        }
        let jets = self.jets(cx.graph, l, f);
        let mut terms = Vec::with_capacity(jets.len());
        for k in 0..jets.len() {
            let mut term = Self::partial(cx, l, &jets, k)?;
            for &x in jets[k].by.iter().rev() {
                term = derivative(cx.graph, term, x)?;
            }
            if jets[k].by.len() % 2 == 1 {
                term = neg(cx.graph, term);
            }
            terms.push(term);
        }
        // Conventionally d/dt ∂L/∂q' - ∂L/∂q: the negative of the sum.
        let sum = add(cx.graph, &terms);
        let value = neg(cx.graph, sum);
        Some(cx.simplify(value))
    }

    fn euler_lagrange(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[l, f, _] = args else {
            return None;
        };
        let l = best(cx.graph, l)?;
        let f = best(cx.graph, f)?;
        if cx.graph.op(f) == core::LIST {
            let functions = cx.graph.children(f).to_vec();
            let mut out = Vec::with_capacity(functions.len());
            for g in functions {
                out.push(self.equation(cx, l, g)?);
            }
            return Some(cx.graph.node(core::LIST, &out));
        }
        self.equation(cx, l, f)
    }

    fn action(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[l, f, path, x, a, b] = args else {
            return None;
        };
        cx.graph.symbol_of(x)?;
        let l = best(cx.graph, l)?;
        let f = best(cx.graph, f)?;
        let along = cx.graph.replace_subterm(l, f, path);
        Some(cx.graph.node(self.defint, &[along, x, a, b]))
    }

    fn first_integral(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[l, f, x] = args else {
            return None;
        };
        let l = best(cx.graph, l)?;
        let f = best(cx.graph, f)?;
        let prime = cx.graph.node(self.diff, &[f, x]);
        let jets = vec![Jet { term: prime, by: vec![x] }, Jet { term: f, by: Vec::new() }];
        let p = Self::partial(cx, l, &jets, 0)?;
        let term = mul(cx.graph, &[prime, p]);
        let value = sub(cx.graph, l, term);
        Some(cx.simplify(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[variational()], src)
    }

    #[test]
    fn harmonic_oscillator() {
        assert_eq!(
            run("euler_lagrange(m*diff(q(t), t)^2/2 - k*q(t)^2/2, q(t), t)"),
            run("k*q(t) + m*diff(diff(q(t), t), t)")
        );
    }

    #[test]
    fn shortest_path_is_a_line() {
        // L = sqrt(1 + y'^2): the extremals have y'' = 0.
        let el = run("euler_lagrange((1 + diff(y(x), x)^2)^(1/2), y(x), x)");
        assert!(el.contains("diff(diff(y(x), x), x)"), "{el}");
        assert!(!el.contains("diff(diff(diff"), "{el}");
    }

    #[test]
    fn several_coordinates_and_higher_derivatives() {
        assert_eq!(
            run("euler_lagrange(diff(x(t), t)^2/2 + diff(y(t), t)^2/2 - x(t)*y(t), list(x(t), y(t)), t)"),
            run("list(y(t) + diff(diff(x(t), t), t), x(t) + diff(diff(y(t), t), t))")
        );
        // Euler–Poisson: L = y''^2 gives 2 y'''' = 0.
        assert_eq!(
            run("euler_lagrange(diff(diff(y(x), x), x)^2, y(x), x)"),
            run("-2*diff(diff(diff(diff(y(x), x), x), x), x)")
        );
    }

    #[test]
    fn field_equation() {
        // L = (u_x^2 + u_y^2)/2 gives Laplace's equation (with a sign).
        assert_eq!(
            run("euler_lagrange((diff(u(x, y), x)^2 + diff(u(x, y), y)^2)/2, u(x, y), list(x, y))"),
            run("diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y)")
        );
    }

    #[test]
    fn extremals_and_action() {
        let solution = run("solve_euler_lagrange(diff(y(x), x)^2, y(x), x)");
        assert_eq!(solution, run("y(x) = C1 + C2*x"));
        assert_eq!(run("action(diff(y(t), t)^2/2, y(t), t, t, 0, 1)"), "1/2");
        assert_eq!(run("first_integral(diff(y(x), x)^2 - y(x)^2, y(x), x)"), run("-diff(y(x), x)^2 - y(x)^2"));
    }
}
