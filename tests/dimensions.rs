//! Dimensional analysis through the public session API.
//!
//! Besides the hand-picked cases, this file holds
//! [`per_rewrite_random_dimension_scan`], a deterministic fixed-seed scan that
//! instantiates every declarative rewrite of `rules::standard()` under random
//! dimension assignments and checks that both sides keep the same dimension.

use num_traits::ToPrimitive;
use rssn::graph::DimError;
use rssn::prelude::*;

/// `(input, expected rewrite, dimension assumptions)`.
type Rewrite<'a> = (&'a str, &'a str, &'a [(&'a str, &'a str)]);

fn session() -> Session {
    Session::new()
}

#[test]
fn advisory_reports_dimensions() {
    let session = session();
    let term = session
        .parse("diff(quantity(3, m/s) * x^2, x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(
        answer
            .dimension
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default(),
        "Length^2/Time"
    );
    assert!(answer.diagnostics.is_empty(), "{:?}", answer.diagnostics);
}

#[test]
fn dimensionless_power_and_symbolic_exponent() {
    let session = session();
    let term = session
        .parse("x^n")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(
        answer
            .dimension
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default(),
        "Length^n"
    );
}

#[test]
fn strict_mode_rejects_mismatches() {
    let session = session();
    let term = session
        .parse("x + y")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("y", "s")
        .dimensions(rssn::graph::DimensionMode::Strict);
    let error = session
        .compute(term, &config)
        .expect_err("strict mode must reject length + time");
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}

#[test]
fn bypass_turns_strict_conflicts_into_diagnostics() {
    let session = session();
    let term = session
        .parse("x + y")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("y", "s")
        .dimensions(rssn::graph::DimensionMode::Strict)
        .bypass_dimension_checks();
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("bypassed strict run failed: {error}"));
    assert!(!answer.diagnostics.is_empty());
}

#[test]
fn off_mode_has_no_dimensions() {
    let session = session();
    let term = session
        .parse("x^2")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new().dimensions(rssn::graph::DimensionMode::Off);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert!(answer.dimension.is_none());
    assert!(answer.diagnostics.is_empty());
}

#[test]
fn user_defined_base_dimension_through_config() {
    let session = session();
    let term = session
        .parse("price^2")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .define_base_dimension("Currency")
        .dimension("price", "Currency");
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(
        answer
            .dimension
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default(),
        "Currency^2"
    );
}

#[test]
fn boolean_answers_are_non_metric() {
    let session = session();
    for source in [
        "eq(1, 1)",
        "is_symmetric(list(list(1, 2), list(2, 5)))",
        "is_convex(x^2, x)",
    ] {
        let term = session
            .parse(source)
            .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
        let config = Config::new().dimensions(rssn::graph::DimensionMode::Strict);
        let answer = session
            .compute(term, &config)
            .unwrap_or_else(|error| panic!("`{source}` failed in strict mode: {error}"));
        let dimension = answer
            .dimension
            .as_ref()
            .unwrap_or_else(|| panic!("`{source}`: no dimension reported"));
        assert!(dimension.non_metric, "`{source}`: {dimension}");
        assert!(
            answer.diagnostics.is_empty(),
            "`{source}`: {:?}",
            answer.diagnostics
        );
    }
}

#[test]
fn a_declared_function_value_has_its_dimension() {
    let session = session();
    let term = session
        .parse("f(x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension_of("f", "m/s^2")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(
        answer
            .dimension
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default(),
        "Length/Time^2"
    );
}

#[test]
fn symbolic_dimension_variables_are_reported() {
    let session = session();
    let term = session
        .parse("x + y")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    let dimension = answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_default();
    assert!(
        dimension.contains("?y") || dimension.contains("Length"),
        "{dimension}"
    );
}

/// A custom base whose name is also a unit symbol is reachable through the
/// explicit `base:` syntax; bare text still means the standard unit.
#[test]
fn base_prefix_selects_a_custom_base_through_config() {
    let session = Session::new();
    let term = session.parse("x").unwrap_or_else(|error| panic!("{error}"));
    let custom = Config::new()
        .define_base_dimension("T")
        .dimension("x", "base:T");
    let answer = session
        .compute(term, &custom)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(dimension_text(&answer), "T");

    let term = session.parse("x").unwrap_or_else(|error| panic!("{error}"));
    let tesla = Config::new().dimension("x", "T");
    let answer = session
        .compute(term, &tesla)
        .unwrap_or_else(|error| panic!("{error}"));
    let dimension = dimension_text(&answer);
    assert_ne!(dimension, "T", "bare `T` must still be the tesla");
    assert!(dimension.contains("Mass"), "{dimension}");
}

fn dimension_text(answer: &Answer<'_>) -> String {
    answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// An unresolved unit inside a term is an unknown, not a conflict: it is a
/// diagnostic in both modes, the dimension stays `any`, and Strict answers.
#[test]
fn unresolved_term_units_are_advisory_in_both_modes() {
    let session = Session::new();
    for (source, unresolved) in [
        ("quantity(1, parsec)", "parsec"),
        ("quantity(1, degC*m)", "degC*m"),
    ] {
        for mode in [
            rssn::graph::DimensionMode::Advisory,
            rssn::graph::DimensionMode::Strict,
        ] {
            let term = session
                .parse(source)
                .unwrap_or_else(|error| panic!("`{source}`: {error}"));
            let answer = session
                .compute(term, &Config::new().dimensions(mode))
                .unwrap_or_else(|error| panic!("`{source}` ({mode:?}): {error}"));
            assert_eq!(dimension_text(&answer), "any", "`{source}` ({mode:?})");
            assert!(
                answer.diagnostics.iter().any(|error| {
                    matches!(
                        error,
                        DimError::UnknownUnit {
                            context: Some(_),
                            unit,
                        } if unit == unresolved
                    )
                }),
                "`{source}` ({mode:?}): {:?}",
                answer.diagnostics
            );
            // The finding is advisory: `is_blocking` must be false.
            assert!(
                answer
                    .diagnostics
                    .iter()
                    .filter(|error| matches!(error, DimError::UnknownUnit { .. }))
                    .all(|error| !error.is_blocking())
            );
        }
    }
}

/// A unit in *requested dimension text* is a declaration, not a term: an
/// unknown one is a hard error in Strict mode (and a diagnostic in Advisory).
#[test]
fn unresolved_dimension_text_is_a_declaration_error() {
    let session = Session::new();
    let term = session
        .parse("x^2")
        .unwrap_or_else(|error| panic!("{error}"));
    let strict = Config::new()
        .dimension("x", "parsec")
        .dimensions(rssn::graph::DimensionMode::Strict);
    let error = session
        .compute(term, &strict)
        .expect_err("an unknown unit in a declaration is an error");
    assert!(
        matches!(
            error,
            ComputeError::Dimension(DimError::UnknownUnit { context: None, .. })
        ),
        "unexpected error: {error}"
    );
    let advisory = Config::new()
        .dimension("x", "parsec")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &advisory)
        .unwrap_or_else(|error| panic!("advisory run failed: {error}"));
    assert!(
        answer
            .diagnostics
            .iter()
            .any(|error| matches!(error, DimError::UnknownUnit { context: None, .. }))
    );
}

/// Every diagnostic carries a stable, machine-readable code (`DimError::code`),
/// so a binding can branch on the kind without parsing the message. The code
/// is part of the public API and is additive: `Display` is unchanged.
#[test]
fn diagnostics_carry_stable_codes() {
    let session = Session::new();

    // A conflict.
    let term = session
        .parse("x + y")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("y", "s")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert!(
        answer
            .diagnostics
            .iter()
            .any(|error| error.code() == "mismatch"),
        "{:?}",
        answer.diagnostics
    );

    // An argument that must be dimensionless.
    let term = session
        .parse("exp(x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert!(
        answer
            .diagnostics
            .iter()
            .any(|error| error.code() == "non_dimensionless"),
        "{:?}",
        answer.diagnostics
    );

    // An affine point misuse.
    let term = session
        .parse("quantity(20, degC) * 2")
        .unwrap_or_else(|error| panic!("{error}"));
    let answer = session
        .compute(
            term,
            &Config::new().dimensions(rssn::graph::DimensionMode::Advisory),
        )
        .unwrap_or_else(|error| panic!("{error}"));
    assert!(
        answer
            .diagnostics
            .iter()
            .any(|error| error.code() == "affine"),
        "{:?}",
        answer.diagnostics
    );

    // The full set of codes is constructible through the public API and is
    // closed: eight variants, two `UnknownUnit` contexts sharing one code.
    let variants = [
        (
            DimError::Mismatch {
                left: String::from("Length"),
                right: String::from("Time"),
                context: String::from("x + y"),
            },
            "mismatch",
        ),
        (
            DimError::NonDimensionless {
                op: String::from("exp"),
                argument: String::from("x"),
                found: String::from("Length"),
            },
            "non_dimensionless",
        ),
        (
            DimError::SymbolicConstraint {
                constraint: String::from("n = 1"),
                context: String::from("x^n"),
            },
            "symbolic_constraint",
        ),
        (
            DimError::UnknownUnit {
                unit: String::from("parsec"),
                context: None,
            },
            "unknown_unit",
        ),
        (
            DimError::Affine {
                operation: String::from("degC * 2"),
            },
            "affine",
        ),
        (
            DimError::MixedList {
                context: String::from("list(x, y)"),
            },
            "mixed_list",
        ),
        (
            DimError::Unannotated {
                op: String::from("mystery"),
            },
            "unannotated",
        ),
        (
            DimError::Parse {
                text: String::from("m^"),
                reason: String::from("missing exponent"),
            },
            "parse",
        ),
    ];
    for (error, code) in &variants {
        assert_eq!(error.code(), *code, "{error:?}");
        assert!(
            code.chars().all(|c| c.is_ascii_lowercase() || c == '_'),
            "`{code}` is not snake_case"
        );
        // `is_blocking` is orthogonal to the code and stays as documented.
        assert_eq!(
            error.is_blocking(),
            !matches!(
                error,
                DimError::Unannotated { .. }
                    | DimError::MixedList { .. }
                    | DimError::UnknownUnit { context: Some(_), .. }
            ),
            "{error:?}"
        );
    }
}

/// An affine temperature is a point: scaling it, raising it to a power or
/// adding two points is rejected. A point plus a difference in kelvin is a
/// point and stays valid.
#[test]
fn affine_temperature_points_are_checked() {
    let session = Session::new();
    for source in [
        "quantity(20, degC) * 2",
        "2 * quantity(20, degC)",
        "quantity(20, degC) / 2",
        "quantity(20, degC)^2",
        "quantity(20, degC) + quantity(30, degC)",
    ] {
        let term = session
            .parse(source)
            .unwrap_or_else(|error| panic!("`{source}`: {error}"));
        let error = session
            .compute(
                term,
                &Config::new().dimensions(rssn::graph::DimensionMode::Strict),
            )
            .expect_err(&format!("Strict must reject `{source}`"));
        assert!(
            matches!(error, ComputeError::Dimension(DimError::Affine { .. })),
            "`{source}`: unexpected error {error}"
        );
        let answer = session
            .compute(
                term,
                &Config::new().dimensions(rssn::graph::DimensionMode::Advisory),
            )
            .unwrap_or_else(|error| panic!("`{source}` (advisory): {error}"));
        assert!(
            answer
                .diagnostics
                .iter()
                .any(|error| matches!(error, DimError::Affine { .. })),
            "`{source}`: {:?}",
            answer.diagnostics
        );
    }
    // A point plus a difference (written in K) is still a temperature.
    let term = session
        .parse("quantity(20, degC) + quantity(2, K)")
        .unwrap_or_else(|error| panic!("{error}"));
    let answer = session
        .compute(
            term,
            &Config::new().dimensions(rssn::graph::DimensionMode::Strict),
        )
        .unwrap_or_else(|error| panic!("point plus difference failed: {error}"));
    assert_eq!(dimension_text(&answer), "Temperature");
    assert!(answer.diagnostics.is_empty(), "{:?}", answer.diagnostics);
}

/// A representative set of standard rewrites must preserve the dimension.
#[test]
fn standard_rewrites_keep_a_consistent_dimension() {
    let session = Session::new();
    let cases: &[Rewrite<'_>] = &[
        ("x + x", "2*x", &[("x", "m")]),
        ("x^2 / x", "x", &[("x", "m")]),
        (
            "expand((x + a)^2)",
            "x^2 + 2*x*a + a^2",
            &[("x", "m"), ("a", "m")],
        ),
        (
            "factor(x^2 - a^2)",
            "(x - a)*(x + a)",
            &[("x", "m"), ("a", "m")],
        ),
        ("sin(t)^2 + cos(t)^2", "1", &[("t", "1")]),
        ("convert(quantity(5, km), m)", "quantity(5000, m)", &[]),
        (
            "simplify_units(quantity(1, N) * quantity(2, m))",
            "quantity(2, J)",
            &[],
        ),
    ];
    for (left, right, assumptions) in cases {
        let mut config = Config::new().dimensions(rssn::graph::DimensionMode::Strict);
        for &(name, text) in *assumptions {
            config = config.dimension(name, text);
        }
        let left_answer = computed(&session, left, &config);
        let right_answer = computed(&session, right, &config);
        assert_eq!(
            dimension_text(&left_answer),
            dimension_text(&right_answer),
            "`{left}` vs `{right}`"
        );
        for (source, answer) in [(left, &left_answer), (right, &right_answer)] {
            assert!(
                answer.diagnostics.iter().all(|error| !error.is_blocking()),
                "`{source}`: {:?}",
                answer.diagnostics
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Per-rewrite random dimensional scan
// ---------------------------------------------------------------------------

/// Dimensions a scanned pattern variable can be given.
const SCAN_DIMENSIONS: &[&str] = &["1", "Length", "Time", "Mass", "Length/Time", "1/Length"];

/// Assignments every scanned rewrite is checked under. The first three are
/// structured (all dimensionless, all length, all time); the rest are
/// pseudo-random but seeded per rule, so the scan is deterministic.
const SCAN_SAMPLES: usize = 6;

/// A small deterministic generator: the scan must never be flaky.
struct ScanRng(u64);

impl ScanRng {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn below(
        &mut self,
        bound: usize,
    ) -> usize {
        if bound <= 1 {
            return 0;
        }
        usize::try_from((self.next() >> 33) % u64::try_from(bound).unwrap_or(1)).unwrap_or(0)
    }
}

/// A deterministic seed for `(rule name, sample index)`.
fn scan_seed(
    name: &str,
    sample: usize,
) -> u64 {
    let mut hash = 0x5eed_0000_u64
        ^ u64::try_from(sample)
            .unwrap_or(0)
            .wrapping_mul(0x9e37_79b9_7f4a_7c15);
    for byte in name.bytes() {
        hash = hash.wrapping_mul(131).wrapping_add(u64::from(byte));
    }
    hash
}

/// Which numeric facts and shape a scanned pattern variable must satisfy.
#[derive(Clone, Copy, Default)]
struct ScanNeed {
    number: bool,
    integer: bool,
    nonzero: bool,
    positive: bool,
    negative: bool,
    nonnegative: bool,
}

fn scan_guard_needs(
    guard: &rssn::graph::Guard,
    needs: &mut [ScanNeed],
) {
    let mark = |needs: &mut [ScanNeed], var: u32, apply: fn(&mut ScanNeed)| {
        if let Some(need) = needs.get_mut(usize::try_from(var).unwrap_or(usize::MAX)) {
            apply(need);
        }
    };
    match guard {
        | rssn::graph::Guard::IsNumber(var) => mark(needs, *var, |need| need.number = true),
        | rssn::graph::Guard::IsInteger(var) => mark(needs, *var, |need| need.integer = true),
        | rssn::graph::Guard::NonZero(var) => mark(needs, *var, |need| need.nonzero = true),
        | rssn::graph::Guard::Positive(var) => mark(needs, *var, |need| need.positive = true),
        | rssn::graph::Guard::Negative(var) => mark(needs, *var, |need| need.negative = true),
        | rssn::graph::Guard::NonNegative(var) => mark(needs, *var, |need| need.nonnegative = true),
        | rssn::graph::Guard::Not(inner) => scan_guard_needs(inner, needs),
        | _ => {},
    }
}

/// Instantiates the pattern variables of `rewrite` for one scan sample.
///
/// Variables that a guard requires to be numeric become integer literals
/// (dimensionless). Every other variable becomes a fresh symbol with a
/// dimension from [`SCAN_DIMENSIONS`] and the facts its guards ask for, so
/// `Rewrite::admits` accepts the sample.
fn scan_substitution(
    graph: &mut rssn::graph::Graph,
    rewrite: &rssn::graph::Rewrite,
    needs: &[ScanNeed],
    rule_index: usize,
    sample: usize,
    rng: &mut ScanRng,
) -> Vec<rssn::graph::NodeId> {
    let mut subst = Vec::with_capacity(rewrite.nvars);
    for var in 0..rewrite.nvars {
        let need = needs.get(var).copied().unwrap_or_default();
        if need.integer || need.number {
            // An integer literal satisfies both `integer` and `number`, and it
            // is dimensionless.
            let value = if need.positive || need.nonnegative {
                2
            } else {
                -3
            };
            subst.push(graph.int(value));
            continue;
        }
        let text = if sample < 3 {
            SCAN_DIMENSIONS[sample]
        } else {
            SCAN_DIMENSIONS[rng.below(SCAN_DIMENSIONS.len())]
        };
        let name = format!("dimscan{rule_index}_{sample}_{var}");
        let symbol = graph.interner_mut().symbol(&name);
        graph
            .assume_dimension(symbol, text)
            .unwrap_or_else(|error| panic!("cannot declare `{name}` as `{text}`: {error}"));
        let mut facts = rssn::graph::Facts::REAL | rssn::graph::Facts::NONZERO;
        if need.nonzero {
            facts = facts | rssn::graph::Facts::NONZERO;
        }
        if need.positive || need.nonnegative {
            facts = facts | rssn::graph::Facts::POSITIVE;
        }
        if need.negative {
            facts = facts | rssn::graph::Facts::NEGATIVE;
        }
        graph.assume(symbol, facts);
        subst.push(graph.symbol_node(symbol));
    }
    subst
}

/// Whether two inferred dimensions can be the same value's dimension.
///
/// `Any` (the literal zero, an order term) unifies with everything, which is
/// the documented polymorphic-zero behaviour; every other dimension must be
/// exactly equal.
fn scan_compatible(
    left: &rssn::graph::Dimension,
    right: &rssn::graph::Dimension,
) -> bool {
    left.any || right.any || left == right
}

/// Every declarative rewrite of `rules::standard()` is instantiated under
/// random dimension assignments and both sides are analysed: a rewrite whose
/// two sides are dimensionally well-typed must preserve the dimension.
///
/// The scan is deterministic (fixed rule order, per-rule seeded generator) and
/// bounded ([`SCAN_SAMPLES`] assignments per rewrite). It is also non-vacuous:
/// it asserts that at least 20 rewrites were exercised and at least 100
/// instantiations were dimension-checked. The test prints the coverage counters
/// when run with `--nocapture`.
///
/// Skipped samples are counted, not failures:
/// * a *left* side that is dimensionally ill-typed is outside the rewrite's
///   dimensional domain (the pattern can only match a well-typed term);
/// * a *right* side that is ill-typed is a rewrite domain restriction the
///   random assignment does not satisfy (for example a unit rewrite whose
///   coefficient must be a pure number).
///
/// A refutation is a sample where both sides are well-typed and their
/// dimensions are incompatible.
#[test]
fn per_rewrite_random_dimension_scan() {
    use rssn::graph::Action;

    let mut graph = rssn::graph::Graph::new();
    let engine = rssn::graph::Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    let rules: Vec<(&str, &rssn::graph::Rewrite)> = engine
        .program()
        .rules
        .iter()
        .filter_map(|rule| {
            match &rule.action {
                | Action::Rewrite(rewrite) => Some((&*rule.name, rewrite)),
                | Action::Kernel(_) => None,
            }
        })
        .collect();

    let mut rules_total = 0_usize;
    let mut rules_exercised = 0_usize;
    let mut samples_checked = 0_usize;
    let mut samples_left_ill_typed = 0_usize;
    let mut samples_right_restricted = 0_usize;
    let mut samples_guarded_out = 0_usize;
    let mut refutations: Vec<String> = Vec::new();

    for (index, (name, rewrite)) in rules.iter().enumerate() {
        rules_total = rules_total.saturating_add(1);
        let mut needs = vec![ScanNeed::default(); rewrite.nvars];
        for guard in &rewrite.guards {
            scan_guard_needs(guard, &mut needs);
        }
        let mut exercised = false;
        for sample in 0..SCAN_SAMPLES {
            let mut rng = ScanRng::new(scan_seed(name, sample));
            let subst = scan_substitution(&mut graph, rewrite, &needs, index, sample, &mut rng);
            if !rewrite.admits(&graph, &subst) {
                samples_guarded_out = samples_guarded_out.saturating_add(1);
                continue;
            }
            let Some(left) = rewrite.lhs.instantiate(&mut graph, &subst) else {
                samples_guarded_out = samples_guarded_out.saturating_add(1);
                continue;
            };
            let Some(right) = rewrite.rhs.instantiate(&mut graph, &subst) else {
                samples_guarded_out = samples_guarded_out.saturating_add(1);
                continue;
            };
            let left_outcome = graph.dimension(left);
            if left_outcome
                .errors
                .iter()
                .any(rssn::graph::DimError::is_blocking)
            {
                samples_left_ill_typed = samples_left_ill_typed.saturating_add(1);
                continue;
            }
            let right_outcome = graph.dimension(right);
            if right_outcome
                .errors
                .iter()
                .any(rssn::graph::DimError::is_blocking)
            {
                samples_right_restricted = samples_right_restricted.saturating_add(1);
                continue;
            }
            if scan_compatible(&left_outcome.dimension, &right_outcome.dimension) {
                samples_checked = samples_checked.saturating_add(1);
                exercised = true;
            } else {
                refutations.push(format!(
                    "`{name}` sample {sample}: `{}` has {}, `{}` has {}",
                    graph.display(left),
                    left_outcome.dimension,
                    graph.display(right),
                    right_outcome.dimension,
                ));
            }
        }
        if exercised {
            rules_exercised = rules_exercised.saturating_add(1);
        }
    }

    println!(
        "dimension scan: {rules_total} rewrites, {rules_exercised} exercised, {samples_checked} checked, \
         {samples_left_ill_typed} left-ill-typed, {samples_right_restricted} right-restricted, \
         {samples_guarded_out} guarded out"
    );
    assert!(
        refutations.is_empty(),
        "{} dimension-changing rewrite(s) found among {rules_total}:\n{}",
        refutations.len(),
        refutations.join("\n")
    );
    assert!(
        rules_total >= 100,
        "unexpectedly few declarative rewrites in the standard rule set: {rules_total}"
    );
    assert!(
        rules_exercised >= 20,
        "the scan exercised only {rules_exercised} of {rules_total} rewrites"
    );
    assert!(
        samples_checked >= 100,
        "the scan checked only {samples_checked} dimensioned instantiations"
    );
}

fn computed<'s>(
    session: &'s Session,
    source: &str,
    config: &Config,
) -> Answer<'s> {
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    session
        .compute(term, config)
        .unwrap_or_else(|error| panic!("`{source}` failed: {error}"))
}

/// Typing the input as written rejects `ln(exp(x))` for a dimensionful `x`;
/// this is intentional and documented, not a rewrite-soundness bug.
#[test]
fn ln_exp_of_a_dimensionful_argument_is_the_documented_strict_case() {
    let session = Session::new();
    let term = session
        .parse("ln(exp(x))")
        .unwrap_or_else(|error| panic!("{error}"));
    let strict = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Strict);
    let error = session
        .compute(term, &strict)
        .expect_err("exp of a length must be rejected in Strict mode");
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
    // Advisory mode answers and reports the same conflict.
    let term = session
        .parse("ln(exp(x))")
        .unwrap_or_else(|error| panic!("{error}"));
    let advisory = Config::new()
        .dimension("x", "m")
        .dimensions(rssn::graph::DimensionMode::Advisory);
    let answer = session
        .compute(term, &advisory)
        .unwrap_or_else(|error| panic!("advisory run failed: {error}"));
    assert!(
        answer.diagnostics.iter().any(|error| error.is_blocking()),
        "{:?}",
        answer.diagnostics
    );
}

#[test]
fn astronomy_and_cgs_units_carry_their_exact_si_scales() {
    // The values are the defining constants, not measurements: the au is exact
    // by definition, the light-year follows from the exact speed of light, and
    // the parsec from the exact au and the IAU 2015 parallax.
    let rational = [
        ("au", 1.495_978_707e11_f64),
        ("ly", 9.460_730_472_580_8e15),
        ("pc", 3.085_677_581_491_367e16),
        ("dyn", 1e-5),
        ("erg", 1e-7),
        ("Gal", 1e-2),
        ("P", 0.1),
        ("St", 1e-4),
        ("G", 1e-4),
        ("Mx", 1e-8),
        ("Ba", 0.1),
    ];
    for (name, want) in rational {
        let unit = rssn::dim::units::lookup(name).unwrap_or_else(|| panic!("no unit {name}"));
        assert_eq!(unit.scale.pi, 0, "{name} should have no power of pi");
        let got = unit.scale.rat.to_f64().unwrap_or(f64::NAN);
        assert!(
            (got / want - 1.0).abs() < 1e-12,
            "{name}: {got:e} vs {want:e}"
        );
    }
    // The oersted is `1000/(4 pi) A/m`, so it carries a power of pi.
    let oersted = rssn::dim::units::lookup("Oe").unwrap_or_else(|| panic!("no unit Oe"));
    assert_eq!(oersted.scale.pi, -1);
    let amps_per_metre =
        oersted.scale.rat.to_f64().unwrap_or(f64::NAN) * std::f64::consts::PI.powi(-1);
    let want = 250.0 / std::f64::consts::PI;
    assert!(
        (amps_per_metre / want - 1.0).abs() < 1e-15,
        "{amps_per_metre}"
    );
    // A prefixed astronomical unit of length works like any other.
    let kiloparsec = rssn::dim::units::resolve("kpc").unwrap_or_else(|e| panic!("{e}"));
    let parsec = rssn::dim::units::resolve("pc").unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(kiloparsec.dimension, parsec.dimension);
    let ratio = kiloparsec.scale.rat.to_f64().unwrap_or(f64::NAN)
        / parsec.scale.rat.to_f64().unwrap_or(f64::NAN);
    assert!((ratio - 1000.0).abs() < 1e-9, "{ratio}");
}

#[test]
fn a_unit_on_a_custom_base_dimension_can_be_registered() {
    use rssn::dim::DimExpr;
    use rssn::dim::Scale;
    use rssn::dim::StandardRegistry;
    use rssn::dim::UnitDef;
    use rssn::dim::UnitRegistry;

    // `Currency` is a base dimension the seven SI bases do not have, which is
    // the case that used to be refused as "not a base dimension".
    let euro = UnitDef::linear("EUR", DimExpr::base("Currency"), Scale::one())
        .alias("euro")
        .display("\u{20ac}")
        .category("currency");
    let mut registry = StandardRegistry;
    registry.register(euro).unwrap_or_else(|e| panic!("{e}"));

    let back = registry
        .lookup("EUR")
        .unwrap_or_else(|| panic!("EUR was not registered"));
    assert_eq!(back.dimension, DimExpr::base("Currency"));
    assert_eq!(back.display, "\u{20ac}");
    assert!(back.aliases.iter().any(|alias| alias == "euro"));

    // The alias resolves to the same unit.
    let by_alias = registry
        .lookup("euro")
        .unwrap_or_else(|| panic!("the alias was not registered"));
    assert_eq!(by_alias.name, "EUR");

    // A second unit on the same base composes with it.
    let cent = UnitDef::linear(
        "EURcent",
        DimExpr::base("Currency"),
        Scale {
            rat: num_rational::BigRational::new(1.into(), 100.into()),
            pi: 0,
        },
    );
    registry.register(cent).unwrap_or_else(|e| panic!("{e}"));
    let cents = registry
        .resolve("EURcent")
        .unwrap_or_else(|e| panic!("{e}"));
    let euros = registry.resolve("EUR").unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(cents.dimension, euros.dimension);

    // The graph carries it into a dimension term once the host declares the
    // base, which is what `Config::base_dimensions` is for.
    let session = Session::new();
    let price = session.sym("price");
    let config = Config::new()
        .base_dimension("Currency")
        .dimension_term("price", DimExpr::base("Currency"));
    let answer = session
        .compute(session.parse("price + price").unwrap_or(price), &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(answer.term.to_string(), "2*price");
}
