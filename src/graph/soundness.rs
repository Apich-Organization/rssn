//! Randomised refutation of rewrite rules.
//!
//! A rewrite claims `lhs = rhs` for every value of its variables that
//! passes its guards. That claim can be tested without any mathematics:
//! substitute random numbers, evaluate both sides, compare. A single
//! disagreement refutes the rule; agreement on many samples does not prove
//! it, but it catches the sign errors, swapped arguments and missing
//! factors that make up nearly all real-world rule bugs.
//!
//! Rules over operators without scalar semantics (unevaluated derivatives,
//! integrals, ...) cannot be tested this way and are reported as
//! inconclusive; they are covered by checking the reduced results of whole
//! requests instead.

use std::collections::HashMap;

use num_complex::Complex64;

use super::facts::Facts;
use super::id::SymbolId;
use super::id::NodeId;
use super::number::Number;
use super::rule::Action;
use super::rule::Env;
use super::rule::Guard;
use super::rule::Program;
use super::rule::Rewrite;
use super::store::Graph;

/// Outcome of checking one rewrite.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Verdict {
    /// Samples on which both sides evaluated and agreed.
    pub agreed: usize,
    /// Samples that could not be evaluated or did not pass the guards.
    pub inconclusive: usize,
}

/// A refuted rule: a sample on which the two sides differ.
#[derive(Clone, Debug, PartialEq)]
pub struct Refutation {
    /// Name of the rule.
    pub rule: String,
    /// The instantiated left-hand side.
    pub lhs: String,
    /// The instantiated right-hand side.
    pub rhs: String,
    /// Value of the left-hand side.
    pub lhs_value: f64,
    /// Value of the right-hand side.
    pub rhs_value: f64,
    /// The bindings of the sample, as `symbol = value`.
    pub sample: Vec<(String, f64)>,
}

impl std::fmt::Display for Refutation {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(
            f,
            "rule `{}` is unsound: {} = {} but {} = {} at {:?}",
            self.rule, self.lhs, self.lhs_value, self.rhs, self.rhs_value, self.sample
        )
    }
}

impl std::error::Error for Refutation {}

/// Deterministic generator so that failures reproduce.
struct Lcg(u64);

impl Lcg {
    const fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    /// A value in `[lo, hi)`.
    fn range(
        &mut self,
        lo: f64,
        hi: f64,
    ) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let unit = (self.next() % 1_000_003) as f64 / 1_000_003.0;
        lo + unit * (hi - lo)
    }
}

/// What the guards of a rule require of one variable.
#[derive(Copy, Clone, Default)]
#[allow(clippy::struct_excessive_bools)] // independent guard flags, not a state machine
struct Need {
    number: bool,
    integer: bool,
    positive: bool,
    negative: bool,
    nonnegative: bool,
}

fn needs(
    rewrite: &Rewrite,
    var: u32,
) -> Need {
    let mut need = Need::default();
    for guard in &rewrite.guards {
        match guard {
            | Guard::IsNumber(v) if *v == var => need.number = true,
            | Guard::IsInteger(v) if *v == var => {
                need.number = true;
                need.integer = true;
            },
            | Guard::Positive(v) if *v == var => need.positive = true,
            | Guard::Negative(v) if *v == var => need.negative = true,
            | Guard::NonNegative(v) if *v == var => need.nonnegative = true,
            | _ => {},
        }
    }
    need
}

/// Tests `rewrite` on `samples` random substitutions.
///
/// Works on a private copy of `graph`, so the caller's graph is untouched.
///
/// # Errors
/// Returns the first sample on which the two sides disagree.
pub fn check_rewrite(
    graph: &Graph,
    name: &str,
    rewrite: &Rewrite,
    samples: usize,
) -> Result<Verdict, Refutation> {
    let mut graph = graph.clone();
    let mut rng = Lcg(0x5eed_0000 ^ name.bytes().fold(0_u64, |h, b| h.wrapping_mul(131).wrapping_add(u64::from(b))));
    let mut verdict = Verdict { agreed: 0, inconclusive: 0 };
    for _ in 0..samples {
        let mut env = Env::numeric(0.0);
        let mut subst = Vec::with_capacity(rewrite.nvars);
        let mut sample = Vec::new();
        for var in 0..rewrite.nvars {
            let need = needs(rewrite, u32::try_from(var).unwrap_or(u32::MAX));
            let (lo, hi) = if need.positive || need.nonnegative {
                (0.25, 3.0)
            } else if need.negative {
                (-3.0, -0.25)
            } else {
                (-3.0, 3.0)
            };
            let node = if need.integer {
                #[allow(clippy::cast_possible_truncation)]
                let value = rng.range(lo, hi + 1.0).floor() as i64;
                let value = if value == 0 && !need.nonnegative { 2 } else { value };
                graph.int(value)
            } else if need.number {
                #[allow(clippy::cast_possible_truncation)]
                let numer = (rng.range(lo, hi) * 4.0).round() as i64;
                let numer = if numer == 0 { 1 } else { numer };
                graph.num(Number::fraction(numer, 4).unwrap_or_else(|| Number::from(1)))
            } else {
                let symbol = graph.interner_mut().fresh_symbol("sample");
                let value = rng.range(lo, hi);
                env.bind(symbol, value);
                sample.push((graph.interner().symbol_name(symbol).to_owned(), value));
                let mut facts = Facts::REAL | Facts::NONZERO;
                if need.positive || need.nonnegative {
                    facts = facts | Facts::POSITIVE;
                }
                if need.negative {
                    facts = facts | Facts::NEGATIVE;
                }
                graph.assume(symbol, facts);
                graph.symbol_node(symbol)
            };
            subst.push(node);
        }
        if !rewrite.admits(&graph, &subst) {
            verdict.inconclusive = verdict.inconclusive.saturating_add(1);
            continue;
        }
        let (Some(lhs), Some(rhs)) =
            (rewrite.lhs.instantiate(&mut graph, &subst), rewrite.rhs.instantiate(&mut graph, &subst))
        else {
            verdict.inconclusive = verdict.inconclusive.saturating_add(1);
            continue;
        };
        match (graph.eval(lhs, &env), graph.eval(rhs, &env)) {
            | (Some(a), Some(b)) if a.is_finite() && b.is_finite() => {
                if (a - b).abs() <= 1e-7 * (1.0 + a.abs().max(b.abs())) {
                    verdict.agreed = verdict.agreed.saturating_add(1);
                } else {
                    return Err(Refutation {
                        rule: name.to_owned(),
                        lhs: graph.display(lhs),
                        rhs: graph.display(rhs),
                        lhs_value: a,
                        rhs_value: b,
                        sample,
                    });
                }
            },
            | _ => {
                // Terms with complex values (the imaginary unit, branches
                // of multi-valued functions) are compared over C.
                let bindings: HashMap<SymbolId, Complex64> =
                    env.bindings().iter().map(|&(s, v)| (s, Complex64::new(v, 0.0))).collect();
                match (graph.eval_complex(lhs, &bindings), graph.eval_complex(rhs, &bindings)) {
                    | (Some(a), Some(b)) if a.is_finite() && b.is_finite() => {
                        if (a - b).norm() <= 1e-7 * (1.0 + a.norm().max(b.norm())) {
                            verdict.agreed = verdict.agreed.saturating_add(1);
                        } else {
                            return Err(Refutation {
                                rule: name.to_owned(),
                                lhs: graph.display(lhs),
                                rhs: graph.display(rhs),
                                lhs_value: a.re,
                                rhs_value: b.re,
                                sample,
                            });
                        }
                    },
                    | _ => verdict.inconclusive = verdict.inconclusive.saturating_add(1),
                }
            },
        }
    }
    Ok(verdict)
}

/// Checks every rewrite of `program`. Returns, per rule, its verdict.
///
/// # Errors
/// Returns the first refutation found.
pub fn check_program(
    graph: &Graph,
    program: &Program,
    samples: usize,
) -> Result<Vec<(String, Verdict)>, Refutation> {
    let mut out = Vec::new();
    for rule in &program.rules {
        if let Action::Rewrite(rewrite) = &rule.action {
            out.push((rule.name.to_string(), check_rewrite(graph, &rule.name, rewrite, samples)?));
        }
    }
    Ok(out)
}

/// Builds the instantiated sides of a rewrite for inspection in tests.
#[must_use]
pub fn sides(
    graph: &mut Graph,
    rewrite: &Rewrite,
    subst: &[NodeId],
) -> Option<(NodeId, NodeId)> {
    Some((rewrite.lhs.instantiate(graph, subst)?, rewrite.rhs.instantiate(graph, subst)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::rule::RuleSet;
    use crate::graph::rule::Tier;
    use crate::graph::schedule::Engine;

    fn program(texts: &'static [&'static str]) -> (Graph, Program) {
        let mut g = Graph::new();
        let set = RuleSet::new("t", move |i| {
            i.op(crate::graph::OpDescriptor::new("ln", crate::graph::Arity::Fixed(1))
                .eval(|a| a.first().map_or(f64::NAN, |x| x.ln())))?;
            i.rewrites(Tier::Normalize, texts)
        });
        let engine = Engine::install(&mut g, &[set]).unwrap_or_else(|e| panic!("{e}"));
        (g, engine.program().clone())
    }

    #[test]
    fn sound_rules_pass() {
        let (g, p) = program(&[
            "distribute: ?a * (?b + ?c) => ?a*?b + ?a*?c",
            "square-of-sum: (?a + ?b)^2 => ?a^2 + 2*?a*?b + ?b^2",
            "ln-prod: ln(?a * ?b) => ln(?a) + ln(?b) if positive(?a), positive(?b)",
            "pow-add: ?a^?m * ?a^?n => ?a^(?m + ?n) if integer(?m), integer(?n)",
        ]);
        let verdicts = check_program(&g, &p, 40).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(verdicts.len(), 4);
        for (name, verdict) in verdicts {
            assert!(verdict.agreed >= 20, "{name}: {verdict:?}");
        }
    }

    #[test]
    fn unsound_rules_are_refuted() {
        for bad in [
            "sign: ?a - ?b => ?b - ?a",
            "freshman: (?a + ?b)^2 => ?a^2 + ?b^2",
            "missing-factor: ?a * (?b + ?c) => ?a*?b + ?c",
            "ln-sum: ln(?a + ?b) => ln(?a) + ln(?b) if positive(?a), positive(?b)",
        ] {
            let texts: &'static [&'static str] = Box::leak(vec![bad].into_boxed_slice());
            let (g, p) = program(texts);
            let result = check_program(&g, &p, 40);
            assert!(result.is_err(), "`{bad}` should have been refuted");
        }
    }

    #[test]
    fn domain_restrictions_are_inconclusive_not_failures() {
        // Without the guard the identity is tested where ln is undefined;
        // NaN samples are skipped, the defined ones agree.
        let (g, p) = program(&["ln-square: ln(?a ^ 2) => 2 * ln(?a) if positive(?a)"]);
        let verdicts = check_program(&g, &p, 30).unwrap_or_else(|e| panic!("{e}"));
        assert!(verdicts.first().is_some_and(|(_, v)| v.agreed == 30));
    }
}
