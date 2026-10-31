//! Dimensional analysis through the public session API.

use rssn::graph::DimError;
use rssn::prelude::*;

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
    let cases: &[(&str, &str, &[(&str, &str)])] = &[
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
