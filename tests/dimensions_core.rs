//! Dimension rules of the core rule sets.
//!
//! `tests/dimensions.rs` covers the dimension typing layer through the
//! session API. This file exercises the semantic rules that the core rule
//! sets attach to their own operators: arithmetic and elementary functions,
//! calculus requests, polynomials, solving, logic, units, radicals,
//! rewriting and common-subexpression elimination.

use rssn::graph::DimError;
use rssn::graph::DimOutcome;
use rssn::graph::DimRule;
use rssn::graph::DimensionMode;
use rssn::graph::Engine;
use rssn::graph::Graph;
use rssn::graph::Payload;
use rssn::prelude::*;

/// The dimension of `source` in a fresh graph with the given symbol
/// assumptions.
fn dimension_of(
    source: &str,
    assumptions: &[(&str, &str)],
) -> DimOutcome {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    for (name, text) in assumptions {
        let symbol = graph.interner_mut().symbol(name);
        graph
            .assume_dimension(symbol, text)
            .unwrap_or_else(|error| panic!("cannot declare `{name}` as `{text}`: {error}"));
    }
    let root = graph
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    graph.dimension(root)
}

/// The printed dimension of `source`, requiring a conflict-free inference.
fn dim(
    source: &str,
    assumptions: &[(&str, &str)],
) -> String {
    let outcome = dimension_of(source, assumptions);
    assert!(
        outcome.errors.is_empty(),
        "`{source}`: unexpected dimension errors {:?}",
        outcome.errors
    );
    outcome.dimension.to_string()
}

/// Whether the request is a structural, non-scalar value.
fn is_non_metric(
    source: &str,
    assumptions: &[(&str, &str)],
) -> bool {
    let outcome = dimension_of(source, assumptions);
    assert!(
        outcome.errors.is_empty(),
        "`{source}`: unexpected dimension errors {:?}",
        outcome.errors
    );
    outcome.dimension.non_metric
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

fn answer_dimension(answer: &Answer<'_>) -> String {
    answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// Unit text wins over base-dimension text; `base:Name` forces a base and
/// resolves a custom base whose name is also a standard unit (`T`).
#[test]
fn dimension_text_uses_units_first_and_supports_base_prefix() {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    graph.define_base_dimension("T");
    graph.define_base_dimension("Currency");
    let arena = graph.dim_arena();

    let tesla = graph
        .dimension_from_text("Wb/m^2")
        .unwrap_or_else(|error| panic!("{error}"));
    // A bare `T` is still the tesla, even with a custom base of that name.
    assert_eq!(
        graph
            .dimension_from_text("T")
            .unwrap_or_else(|error| panic!("{error}"))
            .describe(&arena),
        tesla.describe(&arena)
    );
    // The explicit prefix selects the base dimension.
    assert_eq!(
        graph
            .dimension_from_text("base:T")
            .unwrap_or_else(|error| panic!("{error}"))
            .describe(&arena),
        "T"
    );
    // The prefix also works for the SI bases and for exponents.
    assert_eq!(
        graph
            .dimension_from_text("base:Length")
            .unwrap_or_else(|error| panic!("{error}"))
            .describe(&arena),
        "Length"
    );
    assert_eq!(
        graph
            .dimension_from_text("base:Time^2")
            .unwrap_or_else(|error| panic!("{error}"))
            .describe(&arena),
        "Time^2"
    );
    // A custom base without a colliding unit needs no prefix.
    assert_eq!(
        graph
            .dimension_from_text("Currency")
            .unwrap_or_else(|error| panic!("{error}"))
            .describe(&arena),
        "Currency"
    );
    assert!(graph.dimension_from_text("base:Nope").is_err());
    // The registered name is exact: the `t` alias for `Time` is not applied.
    assert!(graph.dimension_from_text("base:t").is_err());
}

/// A class that contains the literal zero is dimensionally polymorphic:
/// kernels simplify intermediate expressions to zero, and a zero-valued
/// expression must not impose its syntactic dimension on the class.
#[test]
fn a_zero_class_is_dimensionally_polymorphic() {
    // A product with the literal zero is already polymorphic.
    assert_eq!(dim("0 * quantity(2, m)", &[]), "any");
    // Once the engine simplifies a zero-valued difference to the literal
    // zero, the class is polymorphic and reports no conflict.
    let session = Session::new();
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(DimensionMode::Advisory);
    for source in ["x - x", "quantity(2, m) - quantity(2, m)"] {
        let answer = computed(&session, source, &config);
        assert_eq!(answer_dimension(&answer), "any", "`{source}`");
        assert!(
            answer.diagnostics.is_empty(),
            "`{source}`: {:?}",
            answer.diagnostics
        );
    }
}

#[test]
fn a_dimensionful_exponent_is_rejected() {
    let outcome = dimension_of("x^t", &[("x", "m"), ("t", "s")]);
    assert!(
        outcome
            .errors
            .iter()
            .any(|error| matches!(error, DimError::NonDimensionless { .. })),
        "{:?}",
        outcome.errors
    );
    // The dimension still records the symbolic exponent rather than guessing.
    assert_eq!(outcome.dimension.to_string(), "Length^t");
    // `atan2` needs equal dimensions on both arguments.
    let mixed = dimension_of("atan2(x, t)", &[("x", "m"), ("t", "s")]);
    assert!(
        mixed
            .errors
            .iter()
            .any(|error| matches!(error, DimError::Mismatch { .. })),
        "{:?}",
        mixed.errors
    );
}

#[test]
fn non_numeric_literals_are_non_metric() {
    // A truth value, a name or opaque data is not a scalar. Before this was
    // handled, a boolean literal had no rule (`lit`), so every predicate
    // answer was an unknown dimension plus a coverage-gap diagnostic, and
    // Strict mode rejected it as a conflict.
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    for payload in [
        Payload::Bool(true),
        Payload::Bool(false),
        Payload::Str(std::sync::Arc::from("m")),
    ] {
        let node = graph.lit(payload);
        let outcome = graph.dimension(node);
        assert_eq!(outcome.dimension.to_string(), "non-metric");
        assert!(
            outcome.errors.is_empty(),
            "unexpected dimension errors {:?}",
            outcome.errors
        );
    }
}

#[test]
fn dimension_text_rejects_a_malformed_exponent() {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    // A missing or non-numeric exponent used to be read as zero, silently
    // making the unit dimensionless.
    for text in ["m^", "m^x", "m^(1/0)", "m^1/0"] {
        assert!(
            graph.dimension_from_text(text).is_err(),
            "`{text}` should not parse"
        );
    }
    let arena = graph.dim_arena();
    let length_per_time = graph
        .dimension_from_text("m/s")
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(length_per_time.describe(&arena), "Length/Time");
    let root = graph
        .dimension_from_text("m^(1/2)")
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(root.describe(&arena), "sqrt(Length)");
}

#[test]
fn arithmetic_and_elementary_rules() {
    assert_eq!(
        dim("x * y / t^2", &[("x", "m"), ("y", "kg"), ("t", "s")]),
        "Length*Mass/Time^2"
    );
    assert_eq!(dim("sqrt(x^2)", &[("x", "m")]), "Length");
    assert_eq!(dim("E", &[]), "1");
    assert_eq!(dim("atan2(x, x)", &[("x", "m")]), "1");
    assert_eq!(dim("exp(t) + sin(t)", &[("t", "1")]), "1");
}

#[test]
fn advisory_mode_reports_core_dimensions() {
    let session = Session::new();
    let config = Config::new()
        .dimension("x", "m")
        .dimension("y", "s")
        .dimensions(DimensionMode::Advisory);
    let answer = computed(&session, "x * y", &config);
    assert_eq!(answer_dimension(&answer), "Length*Time");
    assert!(answer.diagnostics.is_empty(), "{:?}", answer.diagnostics);
}

#[test]
fn strict_mode_rejects_non_dimensionless_arguments() {
    let session = Session::new();
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(DimensionMode::Strict);
    let term = session
        .parse("exp(x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let error = session
        .compute(term, &config)
        .expect_err("exp of a length must be rejected in strict mode");
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}

#[test]
fn strict_mode_rejects_unit_mismatches() {
    let session = Session::new();
    let config = Config::new().dimensions(DimensionMode::Strict);
    let term = session
        .parse("quantity(1, m) + quantity(1, s)")
        .unwrap_or_else(|error| panic!("{error}"));
    let error = session
        .compute(term, &config)
        .expect_err("adding metres and seconds must be rejected");
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}

#[test]
fn calculus_requests_keep_their_dimension() {
    assert_eq!(
        dim("diff(quantity(3, m/s) * x^2, x)", &[("x", "m")]),
        "Length^2/Time"
    );
    assert_eq!(
        dim("integral(quantity(3, m/s) * x, x)", &[("x", "m")]),
        "Length^3/Time"
    );
    assert_eq!(
        dim(
            "defint(quantity(9.8, m/s^2), t, 0, quantity(3, s))",
            &[("t", "s")]
        ),
        "Length/Time"
    );
    assert_eq!(
        dim("sum(x^2, k, 1, n)", &[("x", "m"), ("k", "1"), ("n", "1")]),
        "Length^2"
    );
    assert_eq!(
        dim("product(x, k, 1, 3)", &[("x", "m"), ("k", "1")]),
        "Length^3"
    );
    assert_eq!(
        dim("antidifference(x, k)", &[("x", "m"), ("k", "1")]),
        "Length"
    );
    assert_eq!(
        dim("indefinite_sum(x^2, k)", &[("x", "m"), ("k", "1")]),
        "Length^2"
    );
    assert_eq!(
        dim("indefinite_product(quantity(2, m), k)", &[("k", "1")]),
        "Length"
    );
    assert_eq!(
        dim(
            "mseries(quantity(2, m) * x^2, list(x), list(0), 2)",
            &[("x", "1")]
        ),
        "Length"
    );
    assert_eq!(
        dim("asymptotic(quantity(2, m) * exp(1/x), x, 2)", &[("x", "1")]),
        "Length"
    );
    assert_eq!(
        dim("mtaylor(sin(x), list(x), list(0), 3)", &[("x", "1")]),
        "1"
    );
}

#[test]
fn strict_mode_accepts_consistent_calculus() {
    let session = Session::new();
    let config = Config::new()
        .dimension("x", "m")
        .dimensions(DimensionMode::Strict);
    let answer = computed(&session, "diff(quantity(3, m/s) * x^2, x)", &config);
    assert_eq!(answer_dimension(&answer), "Length^2/Time");
    assert!(answer.diagnostics.is_empty(), "{:?}", answer.diagnostics);
}

#[test]
fn polynomial_form_requests_keep_their_dimension() {
    assert_eq!(
        dim("expand((x + y)^2)", &[("x", "m"), ("y", "m")]),
        "Length^2"
    );
    assert_eq!(
        dim("factor(x^2 - y^2)", &[("x", "m"), ("y", "m")]),
        "Length^2"
    );
    assert_eq!(dim("degree(x^3, x)", &[("x", "m")]), "1");
    assert_eq!(
        dim("coeff(a * x^2, x, 2)", &[("a", "m"), ("x", "1")]),
        "Length"
    );
    assert_eq!(dim("quo(x^2, x, x)", &[("x", "m")]), "Length");
    assert_eq!(
        dim("expand_trig(sin(x + y))", &[("x", "1"), ("y", "1")]),
        "1"
    );
}

#[test]
fn groebner_bases_keep_the_generator_dimension() {
    let assumptions = [("x", "m"), ("y", "m")];
    // A basis element is a combination of the generators.
    assert_eq!(
        dim("groebner(list(x - y, x + y), list(x, y))", &assumptions),
        "Length"
    );
    assert_eq!(
        dim("groebner_f4(list(x - y, x + y), list(x, y))", &assumptions),
        "Length"
    );
    assert_eq!(
        dim("groebner_f5(list(x - y, x + y), list(x, y))", &assumptions),
        "Length"
    );
}

#[test]
fn strict_mode_accepts_reduced_solutions() {
    let session = Session::new();
    let dimensionless = Config::new()
        .dimension("x", "1")
        .dimensions(DimensionMode::Strict);
    let answer = computed(&session, "solve(x^2 = 4, x)", &dimensionless);
    assert_eq!(answer_dimension(&answer), "1");

    let dimensionful = Config::new()
        .dimension("x", "m")
        .dimension("a", "m")
        .dimensions(DimensionMode::Strict);
    let answer = computed(&session, "solve(x - a, x)", &dimensionful);
    assert_eq!(answer_dimension(&answer), "Length");
}

#[test]
fn strict_mode_accepts_a_dimensionful_groebner_basis() {
    let session = Session::new();
    let config = Config::new()
        .dimension("x", "m")
        .dimension("y", "m")
        .dimensions(DimensionMode::Strict);
    let answer = computed(
        &session,
        "groebner(list(x - y, x + y), list(x, y))",
        &config,
    );
    assert_eq!(answer_dimension(&answer), "Length");
}

#[test]
fn root_of_takes_the_root_of_the_homogeneous_polynomial() {
    // -4 m^2 + x^2 has roots +/- 2 m.
    assert_eq!(dim("root_of(1, quantity(-4, m^2), 0, 1)", &[]), "Length");
}

#[test]
fn strict_mode_accepts_root_of_with_units() {
    let session = Session::new();
    let config = Config::new().dimensions(DimensionMode::Strict);
    let answer = computed(&session, "root_of(1, quantity(-4, m^2), 0, 1)", &config);
    assert_eq!(answer_dimension(&answer), "Length");
}

#[test]
fn algebraic_polynomial_invariants_follow_their_degrees() {
    // p has degree 2 in x and dimension Length^2; q has degree 1 and
    // dimension Length. Res_x(p, q) = Length^(1 * 2 + 2 * 1) = Length^4.
    assert_eq!(
        dim(
            "resultant(quantity(1, m^2) * x^2 + quantity(1, m^2), \
             quantity(1, m) * x + quantity(1, m), x)",
            &[("x", "1")]
        ),
        "Length^4"
    );
    // disc of a quadratic of dimension Length^2 is Length^4.
    assert_eq!(
        dim(
            "discriminant(quantity(1, m^2) * x^2 + quantity(3, m^2) * x \
             + quantity(2, m^2), x)",
            &[("x", "1")]
        ),
        "Length^4"
    );
    assert_eq!(
        dim(
            "leading_coeff(quantity(2, m^3) * x^2 + quantity(3, m^3) * x, x)",
            &[("x", "1")]
        ),
        "Length^3"
    );
    assert_eq!(
        dim(
            "leading_term(quantity(2, m^3) * x^2 + quantity(3, m^3) * x, x)",
            &[("x", "1")]
        ),
        "Length^3"
    );
    assert_eq!(
        dim(
            "leading_monomial(quantity(2, m^3) * x^2 + quantity(3, m^3) * x, x)",
            &[("x", "1")]
        ),
        "1"
    );
}

#[test]
fn solving_returns_values_of_the_unknown() {
    let x = [("x", "m")];
    // Solution sets, boxes and complex roots are values of `x`.
    assert_eq!(dim("solve(x^2 - quantity(4, m^2), x)", &x), "Length");
    assert_eq!(
        dim("solve_general(x^2 - quantity(4, m^2), x)", &x),
        "Length"
    );
    assert_eq!(
        dim("solve_complex(x^2 - quantity(4, m^2), x)", &x),
        "Length"
    );
    assert_eq!(dim("solve_box(x^2 - quantity(4, m^2), x, 1)", &x), "Length");
    assert_eq!(dim("nsolve(x^2 - quantity(4, m^2), x, 1)", &x), "Length");
    // Quantifier elimination returns an equivalent formula.
    assert!(is_non_metric(
        "real_qe(exists(x, x^2 = quantity(4, m^2)), list(x))",
        &x
    ));
    // A decomposition lives in the space of the listed variables.
    assert_eq!(
        dim("cad_cells(list(x^2 - quantity(1, m^2)), list(x))", &x),
        "Length"
    );
}

#[test]
fn logic_requests_are_non_metric_and_counts_are_numbers() {
    let booleans = [("a", "1"), ("b", "1")];
    let x = [("x", "m")];
    assert!(is_non_metric("exists(x, gt(x^2, 0))", &x));
    assert!(is_non_metric(
        "forall(x, gt(x^2 + quantity(1, m^2), 0))",
        &x
    ));
    assert!(is_non_metric("cnf(or(and(a, b), a))", &booleans));
    assert!(is_non_metric("satisfiable(and(a, not(b)))", &booleans));
    assert!(is_non_metric("tautology(or(a, not(a)))", &booleans));
    assert!(is_non_metric("truth_table(or(a, b))", &booleans));
    assert!(is_non_metric("entails(and(a, b), a)", &booleans));
    assert!(is_non_metric(
        "equivalent(implies(a, b), or(not(a), b))",
        &booleans
    ));
    assert_eq!(dim("count_models(or(a, b))", &booleans), "1");
    assert_eq!(dim("bdd_size(or(a, b))", &booleans), "1");
}

#[test]
fn units_requests_report_their_dimension() {
    assert_eq!(dim("quantity(3, m/s)", &[]), "Length/Time");
    assert_eq!(dim("convert(quantity(1, km), m)", &[]), "Length");
    assert_eq!(
        dim("simplify_units(quantity(1, km) + quantity(500, m))", &[]),
        "Length"
    );
    assert_eq!(
        dim("unify_expression(quantity(1, km) + quantity(500, m))", &[]),
        "Length"
    );
    assert!(is_non_metric("dimension(quantity(3, m))", &[]));
}

#[test]
fn radicals_and_rewriting_keep_the_value() {
    assert_eq!(dim("denest_sqrt(sqrt(x^2))", &[("x", "m")]), "Length");
    assert_eq!(dim("simplify_radicals(sqrt(x^2) / x)", &[("x", "m")]), "1");
    assert_eq!(
        dim(
            "rewrite_with(x + y, list(x + y = z), list(x, y, z))",
            &[("x", "m"), ("y", "m"), ("z", "m")]
        ),
        "Length"
    );
    assert!(is_non_metric(
        "knuth_bendix(list(m(x, y) = m(y, x)), list(x, y))",
        &[]
    ));
    assert!(is_non_metric(
        "critical_pairs(list(f(f(x)) = x), list(x))",
        &[]
    ));
}

#[test]
fn cse_keeps_the_expression_dimension_and_reports_counts() {
    let assumptions = [("x", "m"), ("y", "m")];
    // The rewritten expressions keep the dimension of the inputs.
    assert_eq!(dim("cse(list(x + y, x - y))", &assumptions), "Length");
    // The cost is a pair of operation counts.
    assert_eq!(dim("cse_cost(list(x + y))", &assumptions), "1");
}

/// Every operator registered by the core rule sets carries a semantic rule.
/// This is the subset of `graph::dim::tests::every_standard_operator_has_a_dimension_rule`
/// owned by this workstream; other rule sets are annotated separately.
#[test]
fn every_core_operator_has_a_semantic_rule() {
    const OPERATORS: &[&str] = &[
        // elementary
        "E",
        // calculus
        "antidifference",
        "indefinite_sum",
        "indefinite_product",
        "zeilberger",
        "definite_sum",
        "mseries",
        "mtaylor",
        "asymptotic",
        // poly
        "expand",
        "collect",
        "factor",
        "together",
        "cancel",
        "apart",
        "degree",
        "coeff",
        "quo",
        "rem",
        "pgcd",
        "groebner",
        "count_real_roots",
        "real_root_intervals",
        "sturm_sequence",
        "resultant",
        "discriminant",
        "leading_term",
        "leading_coeff",
        "leading_monomial",
        "normal_form",
        "ideal_member",
        "simplify_with_relations",
        "cad",
        "expand_trig",
        "root_of",
        "groebner_f4",
        "groebner_f5",
        // solve
        "solve",
        "solve_general",
        "nsolve",
        "nsolve_certified",
        "solve_box",
        "solve_complex",
        "real_qe",
        "cad_cells",
        // logic
        "cnf",
        "dnf",
        "nnf",
        "simplify_logic",
        "satisfiable",
        "tautology",
        "truth_table",
        "count_models",
        "sat_model",
        "all_models",
        "bdd",
        "bdd_size",
        "exists",
        "forall",
        "restrict",
        "entails",
        "equivalent",
        "prime_implicants",
        // radicals
        "denest_sqrt",
        "simplify_radicals",
        // rewriting
        "knuth_bendix",
        "rewrite_with",
        "critical_pairs",
        // cse
        "cse",
        "cse_cost",
    ];
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    for name in OPERATORS {
        let op = graph
            .ops()
            .lookup(name)
            .unwrap_or_else(|| panic!("operator `{name}` is not registered"));
        let rule = graph
            .dim_rule(op)
            .unwrap_or_else(|| panic!("operator `{name}` has no dimension rule"));
        assert!(
            !matches!(rule, DimRule::Unspecified),
            "operator `{name}` still has an unspecified dimension rule"
        );
    }
}
