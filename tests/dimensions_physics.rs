//! Dimension rules of the physics and optimisation rule sets.
//!
//! The checks drive the public `Session` API with explicit symbol
//! dimensions. In `Strict` mode a missing or contradictory rule is a hard
//! error, so a green test means the operator's semantic rule is installed
//! and agrees with the dimension of the computed answer.
//!
//! The `I` used by the quantum normalisation rules is registered by the
//! `complex` rule set, which has no dimension rule yet; the resulting
//! `Unannotated` finding is explicitly non-blocking (`DimError::is_blocking`),
//! so the tests below require the absence of *blocking* findings rather
//! than an empty diagnostics list.

use rssn::graph::DimensionMode;
use rssn::prelude::*;

/// The printed dimension of an answer.
fn dimension(answer: &Answer<'_>) -> String {
    answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// Whether the answer has a blocking dimensional finding.
fn blocked(answer: &Answer<'_>) -> bool {
    answer.diagnostics.iter().any(|error| error.is_blocking())
}

#[test]
fn nminimize_follows_the_objective_dimension() {
    let session = Session::new();
    // The start point and the objective are symbolic (no numeric literal is
    // mixed into a dimensionful expression), so the request stays a
    // minimisation request and its rule is what is checked.
    let term = session
        .parse("nminimize((x - a)^2 + b, list(x), list(0))")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("a", "m")
        .dimension("b", "m^2")
        .dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(dimension(&answer), "Length^2");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn euler_lagrange_of_a_lagrangian_is_a_force() {
    let session = Session::new();
    // `L = m y'^2/2 - k y^2/2 + c y`. The linear term pins the dimension of
    // the undetermined function `y(x)` through the Lagrangian's own
    // unification.
    let source = "euler_lagrange(m*diff(y(x), x)^2/2 - k*y(x)^2/2 + c*y(x), y(x), x)";
    let config = Config::new()
        .dimension("m", "kg")
        .dimension("x", "s")
        .dimension("k", "kg/s^2")
        .dimension("c", "kg*m/s^2")
        .dimensions(DimensionMode::Strict);
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("{error}"));
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    // `d/dt ∂L/∂y' - ∂L/∂y` is `m y'' + k y - c`, a force.
    assert_eq!(dimension(&answer), "Length*Mass/Time^2");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn action_value_scales_with_the_lagrangian_and_time() {
    let session = Session::new();
    let term = session
        .parse("action(m*diff(y(x), x)^2/2 - k*y(x)^2/2 + c*y(x), y(x), y(x), x, 0, t1)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("m", "kg")
        .dimension("x", "s")
        .dimension("k", "kg/s^2")
        .dimension("c", "kg*m/s^2")
        .dimension("t1", "s")
        .dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    // `∫ L dx` is an energy times a time: an action.
    assert_eq!(dimension(&answer), "Length^2*Mass/Time");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn quantum_expectation_has_the_observable_dimension() {
    let session = Session::new();
    // A one-dimensional wave packet: the amplitude `A` has dimension
    // `Length^(-1/2)`, the packet width `s` and centre `a` are lengths.
    let term = session
        .parse("expectation_value(position_operator(x), A*exp(-(x - a)^2/(2*s^2)), x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("a", "m")
        .dimension("s", "m")
        .dimension("A", "m^(-0.5)")
        .dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    // `⟨psi| x |psi⟩ / ⟨psi|psi⟩` is a length, the observable's dimension.
    assert_eq!(dimension(&answer), "Length");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn probability_density_of_a_normalised_state_is_an_inverse_length() {
    let session = Session::new();
    let term = session
        .parse("probability_density(psi)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("psi", "m^(-0.5)")
        .dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(dimension(&answer), "1/Length");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn normalisation_adds_the_inverse_square_root_length() {
    let session = Session::new();
    let term = session
        .parse("normalize_wavefunction(A*exp(-x^2/(2*s^2)), x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("x", "m")
        .dimension("s", "m")
        .dimension("A", "m^(-0.5)")
        .dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(dimension(&answer), "1/sqrt(Length)");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn a_physics_constant_has_its_si_dimension() {
    let session = Session::new();
    let term = session
        .parse("hbar")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new().dimensions(DimensionMode::Strict);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    // Action: `J s = Length^2*Mass/Time`.
    assert_eq!(dimension(&answer), "Length^2*Mass/Time");
    assert!(!blocked(&answer), "{:?}", answer.diagnostics);
}

#[test]
fn variational_equation_object_is_non_metric() {
    let session = Session::new();
    // `solve_euler_lagrange` returns the solution as an equation object
    // (`y(x) = C1 + C2 x`), which is `NonMetric`. Advisory mode is used
    // because the inner `dsolve` is registered by the `ode` rule set,
    // which does not carry a dimension rule yet; its finding is an
    // `Unannotated` coverage gap in another set.
    let term = session
        .parse("solve_euler_lagrange(diff(y(x), x)^2, y(x), x)")
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new().dimensions(DimensionMode::Advisory);
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    let dimension = answer
        .dimension
        .as_ref()
        .unwrap_or_else(|| panic!("no dimension inferred"));
    assert!(dimension.non_metric, "{dimension:?}");
}

#[test]
fn advisory_reports_a_physics_dimension_conflict() {
    let source = "first_law_thermodynamics(dU, Q, W)";
    let advisory = Session::new();
    let term = advisory
        .parse(source)
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("dU", "J")
        .dimension("Q", "m")
        .dimension("W", "m")
        .dimensions(DimensionMode::Advisory);
    let answer = advisory
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("{error}"));
    assert!(blocked(&answer), "{:?}", answer.diagnostics);

    // The same request is a hard error in strict mode.
    let strict = Session::new();
    let term = strict
        .parse(source)
        .unwrap_or_else(|error| panic!("{error}"));
    let config = Config::new()
        .dimension("dU", "J")
        .dimension("Q", "m")
        .dimension("W", "m")
        .dimensions(DimensionMode::Strict);
    let error = strict
        .compute(term, &config)
        .expect_err("strict mode must reject energy - length");
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}
