//! Dimensional analysis through the public session API.

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
