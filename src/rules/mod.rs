//! # Domain rule sets
//!
//! Every mathematical domain is a [`RuleSet`]: the operators it introduces,
//! the identities it knows as declarative rewrites, and its algorithms as
//! reduction kernels — symbolic, exact and numeric side by side.
//!
//! # Writing a rule set
//!
//! A rule set is a function from an [`Installer`](crate::graph::Installer)
//! to `Result<(), RuleError>`, wrapped with [`RuleSet::new`] and given its
//! dependencies with [`RuleSet::needs`]. Inside it:
//!
//! **Operators.** `i.op(OpDescriptor::new("name", Arity::Fixed(n)))`
//! registers an operator and returns its `OpId`. Give it
//! `.eval(|args| ...)` if it has a value on floats — that alone makes it
//! work in the numeric phase, in [`Term::eval`](crate::api::Term::eval) and
//! in compiled functions. Mark *requests* — operators that stand for a
//! computation still to be done (`diff`, `integral`, `solve`) — with
//! `.flags(OpFlags::HEAVY).cost(100)`; an answer containing one is reported
//! as not reduced. Operators that bind a variable declare it with
//! `.binder(var_index, scope_mask)`.
//!
//! **Rewrites.** `i.rewrites(tier, &["name: lhs => rhs if guard, ..."])`.
//! Terms use `+ - * / ^`, function-call syntax for every other operator,
//! `?x` for pattern variables and bare identifiers for symbols and nullary
//! operators. `<=>` adds both directions. Subtraction, negation and
//! division are sugar for `add`/`mul`/`pow`, so rules never mention them as
//! operators. Guards: `free_of(?a, ?x)`, `number`, `integer`, `symbol`,
//! `nonzero`, `positive`, `negative`, `nonnegative`, `real`, `same(?a, ?b)`,
//! each negatable with `!`.
//!
//! Matching is modulo associativity and commutativity of `+` and `*`. At
//! the root of a pattern extra operands are carried along automatically
//! (`sin(?x)^2 + cos(?x)^2 => 1` fires inside a longer sum). Below the
//! root, the **last** operand of a `+` or `*` pattern, if it is a variable,
//! absorbs whatever the other operands leave: write
//! `integral(?c * ?f, ?x) => ?c * integral(?f, ?x) if free_of(?c, ?x)` and
//! `?f` will be the product of all remaining factors.
//!
//! **Tiers.** [`Tier::Reduce`](crate::graph::Tier::Reduce) for rules that eliminate a heavy operator;
//! [`Tier::Normalize`](crate::graph::Tier::Normalize) for simplifications that never make a term larger
//! (these also run destructively inside tree windows);
//! [`Tier::Explore`](crate::graph::Tier::Explore) for identities that change structure without
//! simplifying (expansion, angle addition) — budgeted and backed off.
//!
//! **Kernels.** Anything procedural implements
//! [`Kernel`](crate::graph::Kernel): look at a node, return
//! [`Outcome::Equal`](crate::graph::Outcome::Equal) with an identical term,
//! [`Outcome::Pinned`](crate::graph::Outcome::Pinned) when the result is a
//! requested *form* (`expand`, `factor`), or
//! [`Outcome::Approx`](crate::graph::Outcome::Approx) with a numeric
//! enclosure valid under the bindings of `cx.env`. Numeric kernels should
//! return [`Outcome::Pass`](crate::graph::Outcome::Pass) unless
//! `cx.env.numeric`. Differentiate/integrate/solve the *best known form* of
//! an argument — `Extractor::new(graph, &[arg], &SizeCost).build(graph, arg)`
//! — rather than the spelling the node was built with.
//!
//! **Attributes.** Per-operator data that another layer consumes is
//! attached with `i.graph().ops_mut().set_attr(op, value)`: derivative
//! rules ([`calculus::Partials`]), parity ([`elementary::Parity`]), facts
//! on real arguments ([`OnReals`](crate::graph::OnReals)).
//!
//! **Soundness.** Every rewrite of a standard set is tested against random
//! numeric samples by this module's test suite; a rule whose two sides
//! disagree anywhere fails the build. Requests are tested in both phases:
//! the symbolic answer, evaluated, must agree with the numeric answer.

pub mod arith;
pub mod calculus;
pub mod combinatorics;
pub mod complex;
pub mod discrete;
pub mod elementary;
pub mod functional;
pub mod geometry;
pub mod linalg;
pub mod logic;
pub mod number_theory;
pub mod ode;
pub mod optimize;
pub mod pde;
pub mod physics;
pub mod poly;
pub mod solve;
pub mod special;
pub mod stats;
pub mod transforms;
pub mod units;
pub mod variational;
pub mod verify;
#[cfg(test)]
pub(crate) mod testing;

pub use arith::arith;
pub use calculus::calculus;
pub use combinatorics::combinatorics;
pub use complex::complex;
pub use discrete::discrete;
pub use elementary::elementary;
pub use functional::functional;
pub use geometry::geometry;
pub use linalg::linalg;
pub use logic::logic;
pub use number_theory::number_theory;
pub use ode::ode;
pub use optimize::optimize;
pub use pde::pde;
pub use physics::physics;
pub use poly::poly;
pub use solve::solve;
pub use special::special;
pub use stats::stats;
pub use transforms::transforms;
pub use units::units;
pub use variational::variational;
pub use verify::verify;

use crate::graph::RuleSet;

/// Every rule set shipped with rssn.
#[must_use]
pub fn standard() -> Vec<RuleSet> {
    vec![arith(), elementary(), calculus(), poly(), solve(), ode(), linalg(), geometry(), complex(), number_theory(), combinatorics(), logic(), special(), stats(), transforms(), variational(), pde(), functional(), physics(), optimize(), verify(), discrete(), units()]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::soundness::check_program;
    use crate::graph::Engine;
    use crate::graph::Graph;

    /// Every declarative rule shipped with rssn is tested against random
    /// numeric samples. Adding a rule to any standard set puts it under
    /// this test automatically.
    #[test]
    fn no_standard_rewrite_is_refuted() {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &standard()).unwrap_or_else(|e| panic!("{e}"));
        let verdicts = check_program(&graph, engine.program(), 60).unwrap_or_else(|e| panic!("{e}"));
        // Rules about requests (derivatives, integrals) have no numeric
        // value to compare; they are covered by the dual-phase tests of
        // their rule sets and listed here explicitly so that the list
        // cannot grow unnoticed.
        let exempt = ["calculus/ftc", "calculus/defint-empty"];
        let untested: Vec<&str> = verdicts
            .iter()
            .filter(|(name, v)| v.agreed == 0 && !exempt.contains(&name.as_str()))
            .map(|(name, _)| name.as_str())
            .collect();
        assert!(untested.is_empty(), "rules that were never actually tested: {untested:?}");
    }
}
