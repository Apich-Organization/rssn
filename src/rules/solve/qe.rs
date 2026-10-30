//! Real quantifier elimination.
//!
//! Quantified formulas over the real numbers are decided and eliminated by
//! a cylindrical algebraic decomposition: a projection to fewer variables,
//! isolation of the real roots of the projection polynomials, lifting over
//! every cell, and a truth evaluation on a sample point of each cell. The
//! answer is a quantifier-free formula in the free parameters, built from
//! sign conditions of the parameter projection polynomials. Algebraic
//! sample points are exact (a number field `Q(theta)` plus a real embedding;
//! see [`field`]), so every sign and every comparison is decided without
//! floating point.
//!
//! Operators:
//!
//! * `real_qe(formula, list(x, ...))` eliminates the quantifiers of
//!   `formula` (and the listed variables additionally as existentials when
//!   they are free) and returns an equivalent quantifier-free formula.
//!   Quantifiers are eliminated from the inside out, so nested and mixed
//!   `exists`/`forall` blocks work.
//! * `cad_cells(polys, list(x1, ..., xn))` returns one
//!   `list(sample_point, open_sections, signs)` per cell of a
//!   sign-invariant cylindrical algebraic decomposition of `R^n`:
//!   `sample_point` is a `list` of exact coordinates (a number or a
//!   `root_of`-expression), `open_sections` is a `list` of booleans that are
//!   `true` for a sector (an open interval) and `false` for a section (a
//!   root), and `signs` is a `list` of `-1`, `0`, `1`, the signs of the
//!   input polynomials at the sample point. The cells are the leaves of the
//!   decomposition, in lifting order.
//! * The logic rule set's `exists` and `forall` are routed here when their
//!   bodies are polynomial comparisons over the reals; genuine Boolean
//!   atoms keep the Boolean decision-diagram semantics.
//!
//! # Input syntax
//!
//! Atoms are comparisons `=`, `!=`, `<`, `<=`, `>`, `>=` (`eq`, `ne`, `lt`,
//! `le`, `gt`, `ge`) whose two sides are polynomials over the rationals in
//! the occurring symbols. Boolean structure is `and`, `or`, `not` (also
//! `implies`, `iff` and `xor`, which are expanded), `true` and `false`.
//! Quantifiers are `exists(prefix, body)` and `forall(prefix, body)` where
//! `prefix` is a symbol or a `list` of symbols. Terms containing a
//! non-polynomial operation are rejected, which leaves the request
//! unevaluated; they are never decided with Boolean semantics.
//!
//! # Exactness and limits
//!
//! The algorithm is exact: projection uses all coefficients, the resultant
//! with the derivative, and the pairwise resultants (Collins' set, which
//! also handles the nullification of leading coefficients), and lifting
//! refines algebraic intervals until signs are decided, using gcd and Sturm
//! counts for equality. Input atoms and projection polynomials are replaced
//! by their distinct irreducible factors over `Q`, so repeated factors are
//! handled on a square-free basis. Every block result is verified
//! independently on a deterministic grid of rational parameter assignments;
//! if the verification disagrees, or any resource bound is reached, the
//! request is left unevaluated rather than answered wrongly.
//!
//! Documented limits (all returning `None`/`Pass`, never a wrong answer):
//! at most four variables, at most 256 projection polynomials per level, at
//! most 20 000 cells, degree at most 48 for a fibre polynomial, at most
//! 20 000 terms in a projection polynomial, a number-field degree of at
//! most 64, and at most 256 interpolation points when a primitive element is
//! built. The dimension limit is the practical one: the cost grows with
//! the product of the projection sizes.

#![allow(
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

pub(crate) mod alg;
pub(crate) mod cad;
pub(crate) mod field;
pub(crate) mod formula;
pub(crate) mod poly;
mod qe_impl;

/// The outcome of routing one `exists`/`forall` request.
pub(crate) enum Routed {
    /// The body is not a real polynomial formula; the Boolean path owns it.
    NotReal,
    /// The eliminated formula.
    Result(NodeId),
    /// It is a real formula, but elimination failed or was not verified.
    Unevaluated,
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let qe_op = i.op(OpDescriptor::new("real_qe", Arity::Fixed(2))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    let cells_op = i.op(OpDescriptor::new("cad_cells", Arity::Fixed(2))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    // `real_qe` returns an equivalent formula; `cad_cells` describes cells of
    // the space of the listed variables, so it carries their dimension.
    i.dim(qe_op, DimRule::NonMetric);
    i.dim(cells_op, DimRule::Same(1));
    i.kernel("solve/real-qe", Tier::Reduce, RealQe { op: qe_op });
    i.kernel("solve/cad-cells", Tier::Reduce, CadCells { op: cells_op });
    Ok(())
}

/// `real_qe(formula, list(x, ...))`.
struct RealQe {
    op: OpId,
}

impl Kernel for RealQe {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let args = graph.children(node).to_vec();
        let &[formula, listed] = args.as_slice() else {
            return Outcome::Pass;
        };
        let extras: Vec<NodeId> = if graph.op(listed) == core::LIST {
            graph.children(listed).to_vec()
        } else {
            vec![listed]
        };
        if extras.iter().any(|&n| graph.as_symbol(n).is_none()) {
            return Outcome::Pass;
        }
        qe_impl::run(graph, formula, &extras).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `cad_cells(polys, list(x1, ..., xn))`.
struct CadCells {
    op: OpId,
}

impl Kernel for CadCells {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let args = graph.children(node).to_vec();
        let &[polys, variables] = args.as_slice() else {
            return Outcome::Pass;
        };
        if graph.op(variables) != core::LIST {
            return Outcome::Pass;
        }
        let variable_nodes = graph.children(variables).to_vec();
        if variable_nodes.is_empty() || variable_nodes.iter().any(|&v| graph.as_symbol(v).is_none())
        {
            return Outcome::Pass;
        }
        let limits = cad::CadLimits::default();
        if variable_nodes.len() > limits.variables {
            return Outcome::Pass;
        }
        let Some(polynomials) = qe_impl::read_polynomials(graph, polys, &variable_nodes) else {
            return Outcome::Pass;
        };
        let Some(cells) = cad::cells(&polynomials, variable_nodes.len(), &limits) else {
            return Outcome::Pass;
        };
        let mut items = Vec::with_capacity(cells.len());
        for cell in &cells {
            let Some(term) = qe_impl::cell_term(graph, cell) else {
                return Outcome::Pass;
            };
            items.push(term);
        }
        Outcome::Pinned(graph.node(core::LIST, &items))
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Routes an `exists`/`forall` request to real QE when its body is a
/// polynomial formula over the reals.
#[must_use]
pub(crate) fn try_quantifier(
    graph: &mut crate::graph::Graph,
    prefix: NodeId,
    body: NodeId,
    universal: bool,
) -> Routed {
    if !contains_comparison(graph, body, 0) {
        return Routed::NotReal;
    }
    qe_impl::try_quantifier(graph, prefix, body, universal)
        .map_or(Routed::Unevaluated, Routed::Result)
}

/// Whether the term contains a comparison, i.e. whether it claims to be a
/// real formula.
fn contains_comparison(
    graph: &crate::graph::Graph,
    node: NodeId,
    depth: usize,
) -> bool {
    if depth > 512 {
        return false;
    }
    let op = graph.op(node);
    if formula::FormulaOps::lookup(graph).is_some_and(|ops| ops.comparison(op).is_some()) {
        return true;
    }
    graph
        .children(node)
        .iter()
        .any(|&child| contains_comparison(graph, child, depth.saturating_add(1)))
}

#[cfg(test)]
mod tests {
    use crate::graph::Env;
    use crate::graph::Graph;
    use crate::graph::NodeId;
    use crate::graph::Payload;
    use crate::graph::op::core;
    use crate::rules::standard;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(source: &str) -> String {
        simplify(&standard(), source)
    }

    /// Truth of a concrete formula under the bindings, with the comparisons
    /// evaluated directly (the answer may contain `=`, which has no scalar
    /// evaluator of its own).
    fn truth(
        graph: &Graph,
        node: NodeId,
        env: &Env,
        depth: usize,
    ) -> Option<bool> {
        if depth > 256 {
            return None;
        }
        if let Some(Payload::Bool(value)) = graph.payload(node) {
            return Some(*value);
        }
        let op = graph.op(node);
        let children = graph.children(node);
        let name = graph.ops().get(op).name.as_ref();
        if name == "true" {
            return Some(true);
        }
        if name == "false" {
            return Some(false);
        }
        if name == "and" {
            let mut result = true;
            for &child in children {
                result &= truth(graph, child, env, depth.saturating_add(1))?;
            }
            return Some(result);
        }
        if name == "or" {
            let mut result = false;
            for &child in children {
                result |= truth(graph, child, env, depth.saturating_add(1))?;
            }
            return Some(result);
        }
        if name == "not" {
            return Some(!truth(
                graph,
                *children.first()?,
                env,
                depth.saturating_add(1),
            )?);
        }
        if name == "implies" {
            let &[a, b] = children else { return None };
            return Some(
                !truth(graph, a, env, depth.saturating_add(1))?
                    || truth(graph, b, env, depth.saturating_add(1))?,
            );
        }
        if name == "iff" {
            let &[a, b] = children else { return None };
            return Some(
                truth(graph, a, env, depth.saturating_add(1))?
                    == truth(graph, b, env, depth.saturating_add(1))?,
            );
        }
        if children.len() == 2
            && (op == core::EQ || matches!(name, "lt" | "le" | "gt" | "ge" | "ne"))
        {
            let a = graph.eval(*children.first()?, env)?;
            let b = graph.eval(*children.get(1)?, env)?;
            let tolerance = 1e-9 * (1.0 + a.abs() + b.abs());
            return Some(match name {
                | "lt" => a < b - tolerance,
                | "le" => a <= b + tolerance,
                | "gt" => a > b + tolerance,
                | "ge" => a >= b - tolerance,
                | "ne" => (a - b).abs() > tolerance,
                | _ => (a - b).abs() <= tolerance,
            });
        }
        None
    }

    fn holds(
        source: &str,
        bindings: &[(&str, f64)],
    ) -> bool {
        let mut graph = Graph::new();
        assert!(crate::graph::Engine::install(&mut graph, &standard()).is_ok());
        let root = graph
            .parse(source)
            .unwrap_or_else(|e| panic!("cannot parse `{source}`: {e}"));
        let mut env = Env::numeric(0.0);
        for (name, value) in bindings {
            env.bind(graph.interner_mut().symbol(name), *value);
        }
        truth(&graph, root, &env, 0).unwrap_or_else(|| panic!("cannot evaluate `{source}`"))
    }

    #[test]
    fn one_variable_quantifiers() {
        assert_eq!(s("real_qe(exists(x, x^2 - 2 = 0), list(x))"), "true");
        assert_eq!(s("real_qe(forall(x, gt(x^2 + 1, 0)), list(x))"), "true");
        assert_eq!(s("real_qe(exists(x, lt(x^2, 0)), list(x))"), "false");
        assert_eq!(s("real_qe(forall(x, ge(x^2, 0)), list(x))"), "true");
    }

    #[test]
    fn parameter_quantifiers() {
        let ball = s("real_qe(exists(y, lt(x^2 + y^2, 1)), list(y))");
        assert!(holds(&ball, &[("x", 0.0)]));
        assert!(holds(&ball, &[("x", 0.5)]));
        assert!(holds(&ball, &[("x", -0.5)]));
        assert!(!holds(&ball, &[("x", 1.0)]));
        assert!(!holds(&ball, &[("x", 2.0)]));
        let square = s("real_qe(forall(y, ge(y^2, x)), list(y))");
        assert!(holds(&square, &[("x", -2.0)]));
        assert!(holds(&square, &[("x", 0.0)]));
        assert!(!holds(&square, &[("x", 1.0)]));
        assert!(!holds(&square, &[("x", 0.5)]));
        let root = s("real_qe(exists(y, y^2 = x), list(y))");
        assert!(holds(&root, &[("x", 0.0)]));
        assert!(holds(&root, &[("x", 3.0)]));
        assert!(!holds(&root, &[("x", -1.0)]));
        let inverse = s("real_qe(exists(y, x*y = 1), list(y))");
        assert!(!holds(&inverse, &[("x", 0.0)]));
        assert!(holds(&inverse, &[("x", 2.0)]));
        assert!(holds(&inverse, &[("x", -2.0)]));
    }

    #[test]
    fn sections_and_degenerate_leading_coefficients() {
        let positive = s("real_qe(exists(y, and(y^2 = x, gt(y, 0))), list(y))");
        assert!(holds(&positive, &[("x", 1.0)]));
        assert!(!holds(&positive, &[("x", 0.0)]));
        assert!(!holds(&positive, &[("x", -1.0)]));
        let linear = s("real_qe(exists(y, a*y = 1), list(y))");
        assert!(holds(&linear, &[("a", 2.0)]));
        assert!(!holds(&linear, &[("a", 0.0)]));
    }

    #[test]
    fn quadratic_case_split() {
        let exists = s("real_qe(exists(x, a*x^2 + b*x + c = 0), list(x))");
        for (a, b, c, want) in [
            (1.0, 0.0, -1.0, true),
            (1.0, 0.0, 1.0, false),
            (0.0, 1.0, 0.0, true),
            (0.0, 1.0, 1.0, true),
            (0.0, 0.0, 0.0, true),
            (0.0, 0.0, 1.0, false),
            (1.0, 2.0, 1.0, true),
            (-1.0, 0.0, 1.0, true),
        ] {
            assert_eq!(
                holds(&exists, &[("a", a), ("b", b), ("c", c)]),
                want,
                "exists a={a} b={b} c={c}"
            );
        }
        let positive = s("real_qe(forall(x, gt(a*x^2 + b*x + c, 0)), list(x))");
        for (a, b, c, want) in [
            (1.0, 0.0, 1.0, true),
            (1.0, 0.0, -1.0, false),
            (0.0, 0.0, 1.0, true),
            (0.0, 0.0, 0.0, false),
            (0.0, 1.0, 0.0, false),
            (-1.0, 0.0, -1.0, false),
            (1.0, 2.0, 1.0, false),
            (1.0, 2.0, 2.0, true),
        ] {
            assert_eq!(
                holds(&positive, &[("a", a), ("b", b), ("c", c)]),
                want,
                "forall a={a} b={b} c={c}"
            );
        }
    }

    #[test]
    fn nested_quantifiers() {
        assert_eq!(
            s("real_qe(forall(a, exists(x, gt(x^2 + a, 0))), list(a))"),
            "true"
        );
        // A closed nested formula: there is an `a` with `x^2 >= a` for all x.
        assert_eq!(
            s("real_qe(exists(a, forall(x, ge(x^2, a))), list(a))"),
            "true"
        );
        // An open nested formula: for every `a` some `x` has `a*x^2 >= b`,
        // exactly when `b <= 0`.
        let bounded = s("real_qe(forall(a, exists(x, ge(a*x^2, b))), list(a))");
        assert!(holds(&bounded, &[("b", -1.0)]));
        assert!(holds(&bounded, &[("b", 0.0)]));
        assert!(!holds(&bounded, &[("b", 1.0)]));
    }

    #[test]
    fn boolean_atoms_are_untouched() {
        // A genuine Boolean atom keeps the decision-diagram semantics.
        assert_eq!(s("exists(y, and(x, y))"), "x");
        let (text, _) = reduce_with(&standard(), "exists(y, sin(y) = 0)", &[]);
        // Non-polynomial real atoms are not silently decided as Booleans.
        assert!(text.contains("exists"), "{text}");
    }

    #[test]
    fn cad_cells_shape() {
        let text = s("cad_cells(list(x^2 - 1), list(x))");
        assert!(text.starts_with("list(list(list("), "{text}");
        // One `list(list(` per cell, plus the outer list.
        assert!(text.matches("list(list(").count() >= 5, "{text}");
        assert_eq!(text.matches("list(false)").count(), 2, "{text}");
    }

    #[test]
    fn repeated_factors_are_decided() {
        // `(y^2 + a*y + b)^2 = 0` is `y^2 + a*y + b = 0`: the discriminant
        // has to enter through the square-free factorisation.
        let repeated = s("real_qe(exists(y, (y^2 + a*y + b)^2 = 0), list(y))");
        for (a, b, want) in [
            (0.0, -1.0, true),
            (0.0, 1.0, false),
            (2.0, 1.0, true),
            (1.0, 1.0, false),
            (0.0, 0.0, true),
        ] {
            assert_eq!(
                holds(&repeated, &[("a", a), ("b", b)]),
                want,
                "repeated at a={a} b={b}: `{repeated}`"
            );
        }
        // A square is never negative.
        assert_eq!(s("real_qe(forall(x, ge((x - 1)^2, 0)), list(x))"), "true");
    }

    #[test]
    fn cad_cells_uses_root_of_for_irrational_sections() {
        // x^3 - 2 has one real, non-quadratic root: the sample point must be
        // an exact `root_of` object, not a float.
        let text = s("cad_cells(list(x^3 - 2), list(x))");
        assert!(text.contains("root_of("), "{text}");
        assert!(!text.contains("0."), "{text}");
    }
}
