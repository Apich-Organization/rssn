use super::egraph::EGraph;

pub mod algebra;
pub mod calculus;
pub mod trig;
pub mod oracles;

/// Trait implemented by all E-Graph rewrite rules and algorithm oracles.
pub trait Rule: Send + Sync + std::fmt::Debug {
    /// Human-readable identifier of the rule.
    fn name(&self) -> &str;

    /// Priority Tier for heuristic scheduling:
    /// - Tier 0: High-weight operator de-cocooning (Derivatives, Integrals, ODEs, Solves)
    /// - Tier 1: Constant folding & Identity reductions (e.g. x+0 -> x, x*0 -> 0)
    /// - Tier 2: Structural algebraic exploration (e.g. Commutativity, Associativity, Distributivity)
    fn tier(&self) -> u8;

    /// Applies the rule to the E-Graph.
    /// Returns the number of new nodes or unions introduced by this rule.
    fn apply(&self, egraph: &mut EGraph) -> usize;
}
