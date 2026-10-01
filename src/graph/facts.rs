//! Sign and reality facts about classes.
//!
//! Many identities hold only under side conditions: `ln(a*b) = ln a + ln b`
//! for positive `a`, `b`; `sqrt(x^2) = x` for non-negative `x`; `x/x = 1`
//! strictly only for non-zero `x`. [`Facts`] is the small lattice those
//! conditions are checked against. Facts come from three sources: literal
//! values, assumptions the user attaches to symbols, and structure
//! (a sum of positives is positive, `exp` of a real is positive).

use std::collections::HashMap;
use std::ops::BitOr;

use super::id::ClassId;
use super::id::NodeId;
use super::id::SymbolId;
use super::number::Number;
use super::op::core;
use super::store::Graph;

/// A set of properties known to hold for a value.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Facts(u8);

impl Facts {
    /// Nothing is known.
    pub const NONE: Self = Self(0);
    /// The value is a real number.
    pub const REAL: Self = Self(1);
    /// The value is not zero.
    pub const NONZERO: Self = Self(1 << 1);
    /// The value is real and `>= 0`.
    pub const NONNEGATIVE: Self = Self(1 | 1 << 2);
    /// The value is real and `<= 0`.
    pub const NONPOSITIVE: Self = Self(1 | 1 << 3);
    /// The value is real and `> 0`.
    pub const POSITIVE: Self = Self(1 | 1 << 1 | 1 << 2);
    /// The value is real and `< 0`.
    pub const NEGATIVE: Self = Self(1 | 1 << 1 | 1 << 3);
    /// The value is an integer.
    pub const INTEGER: Self = Self(1 | 1 << 4);

    /// Whether every fact in `other` is known.
    #[must_use]
    pub const fn has(
        self,
        other: Self,
    ) -> bool {
        self.0 & other.0 == other.0
    }

    /// Facts common to both.
    #[must_use]
    pub const fn and(
        self,
        other: Self,
    ) -> Self {
        Self(self.0 & other.0)
    }

    fn of_f64(value: f64) -> Self {
        if value.is_nan() {
            return Self::NONE;
        }
        let mut facts = Self::REAL;
        if value > 0.0 {
            facts = facts | Self::POSITIVE;
        } else if value < 0.0 {
            facts = facts | Self::NEGATIVE;
        } else {
            facts = facts | Self::NONNEGATIVE | Self::NONPOSITIVE;
        }
        facts
    }
}

impl BitOr for Facts {
    type Output = Self;

    fn bitor(
        self,
        rhs: Self,
    ) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// Operator attribute: facts that hold for the operator's value whenever
/// all its arguments are real (`exp` is positive, `abs` is non-negative,
/// `sin` is real).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct OnReals(pub Facts);

/// Operator attribute: facts that hold for the operator's value whatever
/// its arguments (`re`, `im` and `abs` are real for complex arguments too;
/// the imaginary unit is non-zero).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Always(pub Facts);

impl Graph {
    /// Declares facts about a symbol. Repeated calls accumulate.
    pub fn assume(
        &mut self,
        symbol: SymbolId,
        facts: Facts,
    ) {
        let old = self.assumption(symbol);
        self.set_assumption(symbol, old | facts);
    }

    /// Everything that can be inferred about the value of `node`'s class.
    #[must_use]
    pub fn facts(
        &self,
        node: NodeId,
    ) -> Facts {
        let mut memo = HashMap::new();
        self.class_facts(self.find(node), &mut memo)
    }

    /// Facts of a class, computed once per query. A class that is reached
    /// again while it is being computed (a cycle such as `x = x^1`)
    /// contributes nothing, which is sound: facts are only ever added.
    fn class_facts(
        &self,
        class: ClassId,
        memo: &mut HashMap<ClassId, Facts>,
    ) -> Facts {
        if let Some(&known) = memo.get(&class) {
            return known;
        }
        memo.insert(class, Facts::NONE);
        let facts = self.derive_facts(class, memo);
        memo.insert(class, facts);
        facts
    }

    fn derive_facts(
        &self,
        class: ClassId,
        memo: &mut HashMap<ClassId, Facts>,
    ) -> Facts {
        if let Some(n) = self.class_number(class) {
            let mut facts = Facts::of_f64(n.to_f64());
            if n.is_integer() {
                facts = facts | Facts::INTEGER;
            }
            return facts;
        }
        // A witness on a class without free symbols holds for every binding.
        let mut facts = Facts::NONE;
        if self.free_symbols(class).is_empty() {
            if let Some(ball) = self.approx(class) {
                if ball.mid - ball.rad > 0.0 {
                    facts = facts | Facts::POSITIVE;
                } else if ball.mid + ball.rad < 0.0 {
                    facts = facts | Facts::NEGATIVE;
                }
            }
        }
        // All members denote the same value: what any of them proves
        // holds. Members hidden from matching by leaf collapse still count.
        for member in self.members(class) {
            if !self.is_redundant(member) {
                facts = facts | self.node_facts(member, memo);
            }
        }
        facts
    }

    fn node_facts(
        &self,
        node: NodeId,
        memo: &mut HashMap<ClassId, Facts>,
    ) -> Facts {
        if let Some(symbol) = self.as_symbol(node) {
            return self.assumption(symbol);
        }
        let children = self.children(node);
        let of: Vec<Facts> = children
            .iter()
            .map(|&c| self.class_facts(self.find(c), memo))
            .collect();
        let all = |wanted: Facts| of.iter().all(|f| f.has(wanted));
        let any = |wanted: Facts| of.iter().any(|f| f.has(wanted));
        let mut facts = Facts::NONE;
        match self.op(node) {
            | core::ADD => {
                if all(Facts::REAL) {
                    facts = facts | Facts::REAL;
                }
                if all(Facts::INTEGER) {
                    facts = facts | Facts::INTEGER;
                }
                if all(Facts::NONNEGATIVE) {
                    facts = facts | Facts::NONNEGATIVE;
                    if any(Facts::POSITIVE) {
                        facts = facts | Facts::POSITIVE;
                    }
                }
                if all(Facts::NONPOSITIVE) {
                    facts = facts | Facts::NONPOSITIVE;
                    if any(Facts::NEGATIVE) {
                        facts = facts | Facts::NEGATIVE;
                    }
                }
            },
            | core::MUL => {
                if all(Facts::REAL) {
                    facts = facts | Facts::REAL;
                }
                if all(Facts::INTEGER) {
                    facts = facts | Facts::INTEGER;
                }
                if all(Facts::NONZERO) {
                    facts = facts | Facts::NONZERO;
                }
                if of
                    .iter()
                    .all(|f| f.has(Facts::NONNEGATIVE) || f.has(Facts::NONPOSITIVE))
                {
                    let flips = of.iter().filter(|f| !f.has(Facts::NONNEGATIVE)).count();
                    facts = facts
                        | if flips % 2 == 0 {
                            Facts::NONNEGATIVE
                        } else {
                            Facts::NONPOSITIVE
                        };
                }
            },
            | core::POW => {
                let base = of.first().copied().unwrap_or(Facts::NONE);
                let exp = of.get(1).copied().unwrap_or(Facts::NONE);
                if base.has(Facts::POSITIVE) && exp.has(Facts::REAL) {
                    facts = facts | Facts::POSITIVE;
                }
                if base.has(Facts::REAL) && exp.has(Facts::INTEGER) {
                    facts = facts | Facts::REAL;
                    let exponent = children
                        .get(1)
                        .and_then(|&e| self.number_of(e))
                        .and_then(Number::to_i64);
                    if exponent.is_some_and(|n| n % 2 == 0) {
                        facts = facts | Facts::NONNEGATIVE;
                    }
                    if base.has(Facts::NONZERO) {
                        facts = facts | Facts::NONZERO;
                    }
                }
                if base.has(Facts::INTEGER)
                    && exp.has(Facts::INTEGER)
                    && exp.has(Facts::NONNEGATIVE)
                {
                    facts = facts | Facts::INTEGER;
                }
            },
            | op => {
                if let Some(on_reals) = self.ops().attr::<OnReals>(op) {
                    if all(Facts::REAL) {
                        facts = on_reals.0;
                    }
                }
                if let Some(always) = self.ops().attr::<Always>(op) {
                    facts = facts | always.0;
                }
            },
        }
        facts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;

    fn facts_of(
        g: &mut Graph,
        src: &str,
    ) -> Facts {
        let node = g.parse(src).unwrap_or(NodeId::NONE);
        g.facts(node)
    }

    fn setup() -> Graph {
        let mut g = Graph::new();
        let exp = g
            .ops_mut()
            .register(OpDescriptor::new("exp", Arity::Fixed(1)))
            .unwrap_or(crate::graph::OpId::NONE);
        g.ops_mut().set_attr(exp, OnReals(Facts::POSITIVE));
        let (p, r) = (g.interner_mut().symbol("p"), g.interner_mut().symbol("r"));
        g.assume(p, Facts::POSITIVE);
        g.assume(r, Facts::REAL);
        g
    }

    #[test]
    fn literals() {
        let mut g = setup();
        assert!(facts_of(&mut g, "3").has(Facts::POSITIVE | Facts::INTEGER));
        assert!(facts_of(&mut g, "-1/2").has(Facts::NEGATIVE));
        assert!(!facts_of(&mut g, "-1/2").has(Facts::INTEGER));
        let zero = facts_of(&mut g, "0");
        assert!(
            zero.has(Facts::NONNEGATIVE)
                && zero.has(Facts::NONPOSITIVE)
                && !zero.has(Facts::NONZERO)
        );
    }

    #[test]
    fn assumptions_and_structure() {
        let mut g = setup();
        assert_eq!(facts_of(&mut g, "x"), Facts::NONE);
        assert!(facts_of(&mut g, "p + 1").has(Facts::POSITIVE));
        assert!(facts_of(&mut g, "p * 2 * p").has(Facts::POSITIVE));
        assert!(facts_of(&mut g, "-p").has(Facts::NEGATIVE));
        assert!(facts_of(&mut g, "(-2) * (-p)").has(Facts::POSITIVE));
        assert!(facts_of(&mut g, "r^2").has(Facts::NONNEGATIVE));
        assert!(!facts_of(&mut g, "r^2").has(Facts::NONZERO));
        assert!(facts_of(&mut g, "r^2 + 1").has(Facts::POSITIVE));
        assert!(!facts_of(&mut g, "r^3").has(Facts::NONNEGATIVE));
        assert!(facts_of(&mut g, "p^x").has(Facts::NONE));
        assert!(facts_of(&mut g, "p^r").has(Facts::POSITIVE));
        assert!(facts_of(&mut g, "exp(r)").has(Facts::POSITIVE));
        assert_eq!(
            facts_of(&mut g, "exp(x)"),
            Facts::NONE,
            "exp of a complex number has no sign"
        );
        assert!(!facts_of(&mut g, "p - 1").has(Facts::NONNEGATIVE));
    }

    #[test]
    fn facts_travel_through_equalities() {
        let mut g = setup();
        let x = g.sym("x");
        let known = g.parse("r^2 + 1").unwrap_or(NodeId::NONE);
        assert_eq!(g.facts(x), Facts::NONE);
        g.union(x, known);
        g.rebuild();
        assert!(g.facts(x).has(Facts::POSITIVE));
    }
}
