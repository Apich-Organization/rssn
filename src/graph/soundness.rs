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
//!
//! Parallelism: [`check_program`] checks the rewrites concurrently (each on its own graph copy and
//! with a name-seeded generator) and returns the first refutation in rule order.

use std::collections::HashMap;

use num_complex::Complex64;

use super::facts::Facts;
use super::id::NodeId;
use super::id::SymbolId;
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
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
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

/// Values sitting on or next to the singular loci of the elementary
/// functions: zero and the log branch cut, `±1` for `asin`/`acos`, `±π/2` for
/// `tan`.
const SINGULAR_VALUES: [f64; 14] = [
    f64::MIN_POSITIVE,
    -f64::MIN_POSITIVE,
    f64::EPSILON,
    -f64::EPSILON,
    1e-8,
    -1e-8,
    1.0 - 1e-8,
    1.0 + 1e-8,
    -1.0 + 1e-8,
    -1.0 - 1e-8,
    std::f64::consts::FRAC_PI_2 - 1e-9,
    std::f64::consts::FRAC_PI_2 + 1e-9,
    -std::f64::consts::FRAC_PI_2 + 1e-9,
    -std::f64::consts::FRAC_PI_2 - 1e-9,
];

/// Large magnitudes used to build near-cancellation samples: two operands of
/// opposite sign leave a small, ill-conditioned result even though every
/// intermediate is large. The magnitudes stay modest because several special
/// functions evaluate very slowly (or lose the real branch) for huge
/// arguments, and they are deliberately *not* integer-valued: a large
/// negative integer sends integer-indexed special functions (Hurwitz zeta,
/// Bernoulli polynomials) into their recursive non-closed-form paths.
const LARGE_VALUES: [f64; 4] = [1e3 + 0.37, -1e3 + 0.37, 1e3 + 0.11, -1e3 + 0.11];

/// A sample value either near a singular locus or half of a near-cancellation
/// pair, adapted to the variable's guard requirements.
///
/// `first_symbol` is the value chosen for variable 0 when it is a symbol: an
/// unconstrained later variable is then set to its negation plus a small
/// offset, which is where `(a + b)` cancels.
fn adversarial_value(
    round: usize,
    var: usize,
    need: Need,
    first_symbol: Option<f64>,
) -> f64 {
    let constrained = need.positive || need.nonnegative || need.negative;
    let raw = if round % 2 == 1 {
        if var == 0 {
            LARGE_VALUES[(round / 2) % LARGE_VALUES.len()]
        } else if !constrained {
            let base = first_symbol.unwrap_or(LARGE_VALUES[0]);
            // The partner `-base + offset` is deliberately not
            // integer-valued: a large negative integer sends integer-indexed
            // special functions into their slow recursion paths.
            let offset = if base.abs() > 1.0 {
                1e-2
            } else {
                1e-9
            };
            -base + offset
        } else {
            SINGULAR_VALUES[(round.saturating_add(var)) % SINGULAR_VALUES.len()]
        }
    } else {
        SINGULAR_VALUES[(round / 2 + var) % SINGULAR_VALUES.len()]
    };
    if need.positive || need.nonnegative {
        // Keep the sign the guard asserted while staying near the locus.
        raw.abs().max(f64::MIN_POSITIVE)
    } else if need.negative {
        -raw.abs().max(f64::MIN_POSITIVE)
    } else {
        raw
    }
}

/// Evaluates `node` like [`Graph::eval`] while also summing the absolute
/// values of every subexpression.
///
/// That sum is the scale at which floating-point rounding enters: an identity
/// evaluated through cancelling intermediates can differ by `eps * scale`
/// even when its final value is tiny. Returns `None` when the term does not
/// evaluate to a finite number.
fn eval_with_scale(
    graph: &Graph,
    node: NodeId,
    env: &Env,
) -> Option<(f64, f64)> {
    let mut values: HashMap<NodeId, (f64, f64)> = HashMap::new();
    let mut stack = vec![(node, false)];
    let mut args: Vec<f64> = Vec::new();
    while let Some((cur, expanded)) = stack.pop() {
        if values.contains_key(&cur) {
            continue;
        }
        let children = graph.children(cur);
        if children.is_empty() {
            let value = match (graph.as_number(cur), graph.as_symbol(cur)) {
                | (Some(n), _) => n.to_f64(),
                | (None, Some(s)) => env.value(s)?,
                | (None, None) => (graph.ops().get(graph.op(cur)).eval?)(&[]),
            };
            values.insert(cur, (value, value.abs()));
        } else if expanded {
            let eval = graph.ops().get(graph.op(cur)).eval?;
            args.clear();
            let mut scale = 0.0;
            let mut complete = true;
            for child in children {
                let Some(&(value, below)) = values.get(child) else {
                    complete = false;
                    break;
                };
                args.push(value);
                scale += below;
            }
            if !complete {
                return None;
            }
            let value = eval(&args);
            if !value.is_finite() {
                return None;
            }
            values.insert(cur, (value, scale + value.abs()));
        } else {
            stack.push((cur, true));
            stack.extend(
                children
                    .iter()
                    .filter(|c| !values.contains_key(c))
                    .map(|&c| (c, false)),
            );
        }
    }
    values.get(&node).copied()
}

/// The largest change of any of `terms` caused by nudging one bound symbol by
/// a relative `1e-6`, over both directions: the local Lipschitz scale of the
/// identity at the sample.
fn finite_difference_scale(
    graph: &Graph,
    terms: &[NodeId],
    env: &Env,
) -> f64 {
    let mut scale = 0.0_f64;
    for &(symbol, value) in env.bindings() {
        let step = 1e-6 * (1.0 + value.abs());
        for direction in [1.0, -1.0] {
            let mut nudged = env.clone();
            nudged.bind(symbol, value + direction * step);
            for &term in terms {
                if let (Some(base), Some(shifted)) =
                    (graph.eval(term, env), graph.eval(term, &nudged))
                    && base.is_finite()
                    && shifted.is_finite()
                {
                    scale = scale.max((shifted - base).abs());
                }
            }
        }
    }
    scale
}

/// Tolerance for a singular-adjacent sample.
///
/// Near a pole, a branch cut or a cancellation, two mathematically equal
/// expressions can disagree by much more than `1e-7` relative: the rounding
/// error is proportional to the magnitude of the *intermediate* results, and
/// the derivative says how far the neighbourhood is swept by one rounding
/// step. The tolerance is the largest of the plain relative bound, one part in
/// `1e9` of the intermediate magnitude, and one part in `1e6` of the local
/// finite-difference scale.
fn singular_tolerance(
    graph: &Graph,
    lhs: NodeId,
    rhs: NodeId,
    env: &Env,
    a: f64,
    b: f64,
) -> f64 {
    let base = 1e-7 * (1.0 + a.abs().max(b.abs()));
    let magnitude = eval_with_scale(graph, lhs, env)
        .zip(eval_with_scale(graph, rhs, env))
        .map_or(1.0, |((_, lhs_scale), (_, rhs_scale))| {
            lhs_scale.max(rhs_scale).max(1.0)
        });
    let derivative = finite_difference_scale(graph, &[lhs, rhs], env);
    base.max(1e-9 * magnitude).max(1e-6 * derivative)
}

/// Whether two evaluations disagree *robustly*, or only at the sampled point.
///
/// At a discontinuity the two sides of a sound identity can land on different
/// roundings of a boundary: `floor(x + 2)` and `floor(x) + 2` differ for a
/// negative `x` so small that `x + 2` rounds to exactly `2.0`. Such a
/// disagreement disappears when the sample is nudged by an absolute
/// `1e-9`-scale amount, whereas a genuine counterexample survives the nudge.
/// Requires the disagreement on at least two of the three nudges.
fn refutation_survives_jitter(
    graph: &Graph,
    lhs: NodeId,
    rhs: NodeId,
    env: &Env,
) -> bool {
    let mut disagreements = 0_usize;
    for &jitter in &[1e-9, -1e-9, 1e-6] {
        let mut nudged = env.clone();
        for &(symbol, value) in env.bindings() {
            nudged.bind(symbol, value + jitter);
        }
        let (Some(a), Some(b)) = (graph.eval(lhs, &nudged), graph.eval(rhs, &nudged)) else {
            continue;
        };
        if !a.is_finite() || !b.is_finite() {
            continue;
        }
        if (a - b).abs() > singular_tolerance(graph, lhs, rhs, &nudged, a, b) {
            disagreements = disagreements.saturating_add(1);
        }
    }
    disagreements >= 2
}

/// Tests `rewrite` on `samples` random substitutions.
///
/// Works on a private copy of `graph`, so the caller's graph is untouched.
/// Every fourth sample (when at least eight are requested) is drawn near a
/// singular locus or a cancellation instead of uniformly, and is compared
/// with the conditioning-aware tolerance of [`singular_tolerance`]. A
/// refuted sample is still reported as a refutation: sampling is a smoke
/// test, not a proof, and the singular-adjacent samples deliberately trade
/// refutation power for robustness against rounding at the loci. A
/// singular-adjacent sample whose sides do not evaluate over the reals is
/// inconclusive, not a refutation, because the rewrite's claim is about real
/// values; integer and rational sample magnitudes are clamped so that a
/// sample cannot send a special function into a pathological recursion.
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
    let mut rng = Lcg(0x5eed_0000
        ^ name
            .bytes()
            .fold(0_u64, |h, b| h.wrapping_mul(131).wrapping_add(u64::from(b))));
    let mut verdict = Verdict {
        agreed: 0,
        inconclusive: 0,
    };
    for sample_index in 0..samples {
        let adversarial = samples >= 8 && sample_index % 4 == 3;
        let round = sample_index / 4;
        let mut env = Env::numeric(0.0);
        let mut subst = Vec::with_capacity(rewrite.nvars);
        let mut sample = Vec::new();
        let mut first_symbol = None;
        for var in 0..rewrite.nvars {
            let need = needs(rewrite, u32::try_from(var).unwrap_or(u32::MAX));
            let (lo, hi) = if need.positive || need.nonnegative {
                (0.25, 3.0)
            } else if need.negative {
                (-3.0, -0.25)
            } else {
                (-3.0, 3.0)
            };
            // The adversarial value is cheap and depends only on the guard
            // requirements and the first symbol chosen; the uniform path uses
            // the generator instead.
            let near = adversarial_value(round, var, need, first_symbol);
            let node = if need.integer {
                #[allow(clippy::cast_possible_truncation)]
                let value = if adversarial {
                    // Integer arguments index recursions and powers: keep them
                    // small so an adversarial sample cannot ask a special
                    // function for its millionth Bernoulli polynomial.
                    near.round().clamp(-12.0, 12.0) as i64
                } else {
                    rng.range(lo, hi + 1.0).floor() as i64
                };
                let value = if value == 0 && !need.nonnegative {
                    if need.negative { -1 } else { 2 }
                } else {
                    value
                };
                graph.int(value)
            } else if need.number {
                #[allow(clippy::cast_possible_truncation)]
                let numer = (if adversarial {
                    near
                } else {
                    rng.range(lo, hi)
                } * 4.0)
                    .round()
                    .clamp(-64.0, 64.0) as i64;
                let numer = if numer == 0 {
                    if need.negative { -1 } else { 1 }
                } else {
                    numer
                };
                graph.num(Number::fraction(numer, 4).unwrap_or_else(|| Number::from(1)))
            } else {
                let symbol = graph.interner_mut().fresh_symbol("sample");
                let value = if adversarial {
                    near
                } else {
                    rng.range(lo, hi)
                };
                if first_symbol.is_none() {
                    first_symbol = Some(value);
                }
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
        let (Some(lhs), Some(rhs)) = (
            rewrite.lhs.instantiate(&mut graph, &subst),
            rewrite.rhs.instantiate(&mut graph, &subst),
        ) else {
            verdict.inconclusive = verdict.inconclusive.saturating_add(1);
            continue;
        };
        match (graph.eval(lhs, &env), graph.eval(rhs, &env)) {
            | (Some(a), Some(b)) if a.is_finite() && b.is_finite() => {
                let residual = (a - b).abs();
                let plain = 1e-7 * (1.0 + a.abs().max(b.abs()));
                // The conditioning estimate costs several extra evaluations,
                // so it only runs when the plain tolerance already failed.
                let agrees = if adversarial {
                    residual <= plain
                        || residual <= singular_tolerance(&graph, lhs, rhs, &env, a, b)
                } else {
                    residual <= plain
                };
                if agrees {
                    verdict.agreed = verdict.agreed.saturating_add(1);
                } else if adversarial && !refutation_survives_jitter(&graph, lhs, rhs, &env) {
                    // The disagreement may be a rounding artefact at a
                    // discontinuity; only a robust one refutes.
                    verdict.inconclusive = verdict.inconclusive.saturating_add(1);
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
                if adversarial {
                    // The claim is about real values. If a singular-adjacent
                    // sample left the real axis, the principal complex values
                    // can differ by branch and rounding artifacts that say
                    // nothing about the rewrite: report it as inconclusive
                    // rather than a refutation.
                    verdict.inconclusive = verdict.inconclusive.saturating_add(1);
                    continue;
                }
                // Terms with complex values (the imaginary unit, branches
                // of multi-valued functions) are compared over C.
                let bindings: HashMap<SymbolId, Complex64> = env
                    .bindings()
                    .iter()
                    .map(|&(s, v)| (s, Complex64::new(v, 0.0)))
                    .collect();
                match (
                    graph.eval_complex(lhs, &bindings),
                    graph.eval_complex(rhs, &bindings),
                ) {
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
    // Every rewrite is checked on its own private copy of the graph with a
    // generator seeded from the rule name, so the rules are independent: they
    // run in parallel and the first refutation in rule order is returned,
    // exactly as in the sequential scan.
    let rewrites: Vec<(&str, &Rewrite)> = program
        .rules
        .iter()
        .filter_map(|rule| {
            if let Action::Rewrite(rewrite) = &rule.action {
                Some((&*rule.name, rewrite))
            } else {
                None
            }
        })
        .collect();
    let results = crate::kernels::par::map_slice(&rewrites, 8, |&(name, rewrite)| {
        check_rewrite(graph, name, rewrite, samples).map(|verdict| (name.to_owned(), verdict))
    });
    results.into_iter().collect()
}

/// Builds the instantiated sides of a rewrite for inspection in tests.
#[must_use]
pub fn sides(
    graph: &mut Graph,
    rewrite: &Rewrite,
    subst: &[NodeId],
) -> Option<(NodeId, NodeId)> {
    Some((
        rewrite.lhs.instantiate(graph, subst)?,
        rewrite.rhs.instantiate(graph, subst)?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::core;
    use crate::graph::rule::RuleSet;
    use crate::graph::rule::Tier;
    use crate::graph::schedule::Engine;

    fn program(texts: &'static [&'static str]) -> (Graph, Program) {
        let mut g = Graph::new();
        let set = RuleSet::new("t", move |i| {
            i.op(
                crate::graph::OpDescriptor::new("ln", crate::graph::Arity::Fixed(1))
                    .eval(|a| a.first().map_or(f64::NAN, |x| x.ln())),
            )?;
            i.op(
                crate::graph::OpDescriptor::new("floor", crate::graph::Arity::Fixed(1))
                    .eval(|a| a.first().map_or(f64::NAN, |x| x.floor())),
            )?;
            i.op(
                crate::graph::OpDescriptor::new("ceil", crate::graph::Arity::Fixed(1))
                    .eval(|a| a.first().map_or(f64::NAN, |x| x.ceil())),
            )?;
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
        // NaN samples are skipped, the defined ones agree. The adversarial
        // samples deliberately include values so small that `a^2` underflows
        // to zero, which is one more inconclusive domain case.
        let (g, p) = program(&["ln-square: ln(?a ^ 2) => 2 * ln(?a) if positive(?a)"]);
        let verdicts = check_program(&g, &p, 30).unwrap_or_else(|e| panic!("{e}"));
        assert!(
            verdicts
                .first()
                .is_some_and(|(_, v)| v.agreed >= 20 && v.inconclusive >= 1),
            "{verdicts:?}"
        );
    }

    #[test]
    fn adversarial_values_pair_for_cancellation() {
        let need = Need::default();
        let base = adversarial_value(1, 0, need, None);
        assert_eq!(base, 1e3 + 0.37);
        let partner = adversarial_value(1, 1, need, Some(base));
        assert!((partner + base).abs() < 0.02, "{partner}");
        // Neither value is integer-valued: a large negative integer sends
        // integer-indexed special functions into a very slow recursion.
        assert!((base - base.round()).abs() > 0.01, "{base}");
        assert!((partner - partner.round()).abs() > 0.01, "{partner}");
        // A guarded variable keeps the sign its guard asserted.
        let positive = Need {
            positive: true,
            ..Need::default()
        };
        assert!(adversarial_value(0, 0, positive, None) > 0.0);
        let negative = Need {
            negative: true,
            ..Need::default()
        };
        assert!(adversarial_value(0, 0, negative, None) < 0.0);
        // The singular vector carries both signs of zero, values next to ±1
        // and the tan poles.
        let values: Vec<f64> = (0..SINGULAR_VALUES.len())
            .map(|i| adversarial_value(2 * i, 0, need, None))
            .collect();
        assert!(values.iter().any(|v| v.abs() < 1e-300 && *v > 0.0));
        assert!(values.iter().any(|v| v.abs() < 1e-300 && *v < 0.0));
        assert!(values.iter().any(|v| (v.abs() - 1.0).abs() < 1e-6));
        assert!(
            values
                .iter()
                .any(|v| (v.abs() - std::f64::consts::FRAC_PI_2).abs() < 1e-6)
        );
    }

    #[test]
    fn near_cancellation_needs_the_conditioning_tolerance() {
        // `(a + b)^2` and its expansion through cancelling large terms agree
        // mathematically but not to the plain relative tolerance when the
        // result is tiny.
        let mut g = Graph::new();
        let a = g.sym("a");
        let b = g.sym("b");
        let two = g.int(2);
        let sum = g.node(core::ADD, &[a, b]);
        let lhs = g.node(core::POW, &[sum, two]);
        let a_sq = g.node(core::POW, &[a, two]);
        let b_sq = g.node(core::POW, &[b, two]);
        let cross = g.node(core::MUL, &[two, a, b]);
        let rhs = g.node(core::ADD, &[a_sq, cross, b_sq]);
        let (sa, sb) = (
            g.interner().find_symbol("a").unwrap_or(SymbolId::NONE),
            g.interner().find_symbol("b").unwrap_or(SymbolId::NONE),
        );
        let mut env = Env::numeric(0.0);
        env.bind(sa, 1e10);
        env.bind(sb, -1e10 + 1e-2);
        let left = g.eval(lhs, &env).unwrap_or(f64::NAN);
        let right = g.eval(rhs, &env).unwrap_or(f64::NAN);
        let residual = (left - right).abs();
        let plain = 1e-7 * (1.0 + left.abs().max(right.abs()));
        assert!(
            residual > plain,
            "residual {residual} should exceed the plain tolerance {plain}"
        );
        let conditioning = singular_tolerance(&g, lhs, rhs, &env, left, right);
        assert!(
            residual <= conditioning,
            "conditioning tolerance {conditioning} must cover {residual}"
        );
    }

    #[test]
    fn well_conditioned_identity_near_singular_loci_is_not_refuted() {
        let (g, p) =
            program(&["ln-prod: ln(?a * ?b) => ln(?a) + ln(?b) if positive(?a), positive(?b)"]);
        let verdicts = check_program(&g, &p, 40).unwrap_or_else(|e| panic!("{e}"));
        let (_, verdict) = verdicts.first().unwrap_or_else(|| panic!("no verdicts"));
        assert!(verdict.agreed >= 20, "{verdict:?}");
        assert!(
            verdict.inconclusive >= 1,
            "a sample next to zero underflows and is skipped, {verdict:?}"
        );
    }

    #[test]
    fn rounding_at_a_discontinuity_is_not_a_refutation() {
        // Mathematically exact, but for a negative `x` within one ulp of zero
        // `x + 2` rounds to exactly 2.0, so the two sides evaluate to 2 and 1.
        // The disagreement is an artefact of the sample, not of the rewrite.
        let (g, p) = program(&[
            "floor-shift: floor(?a + 2) => floor(?a) + 2",
            "ceil-shift: ceil(?a + 2) => ceil(?a) + 2",
        ]);
        let verdicts = check_program(&g, &p, 40).unwrap_or_else(|e| panic!("{e}"));
        for (name, verdict) in verdicts {
            assert!(verdict.agreed >= 20, "{name}: {verdict:?}");
        }
    }
}
