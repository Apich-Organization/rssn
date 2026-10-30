//! Dimension rules for the special, transform, linear-algebra, statistics,
//! number-theory, combinatorics, discrete, Lie, geometric-algebra and
//! geometry rule sets.
//!
//! Every test names a representative operator and asserts the dimension of
//! the answer in Strict mode (a conflict would be a hard error) or, for the
//! negative tests, that Strict mode does reject an ill-typed request.

use rssn::graph::DimensionMode;
use rssn::prelude::*;

fn config(assumptions: &[(&str, &str)]) -> Config {
    let mut config = Config::new().dimensions(DimensionMode::Strict);
    for &(name, text) in assumptions {
        config = config.dimension(name, text);
    }
    config
}

fn dimension(
    source: &str,
    assumptions: &[(&str, &str)],
) -> String {
    let session = Session::new();
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    let answer = session
        .compute(term, &config(assumptions))
        .unwrap_or_else(|error| panic!("`{source}`: {error}"));
    answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_else(|| panic!("`{source}`: no dimension reported"))
}

fn dimension_advisory(
    source: &str,
    assumptions: &[(&str, &str)],
) -> (String, usize) {
    let session = Session::new();
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    let mut config = Config::new().dimensions(DimensionMode::Advisory);
    for &(name, text) in assumptions {
        config = config.dimension(name, text);
    }
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("`{source}`: {error}"));
    (
        answer
            .dimension
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default(),
        answer.diagnostics.len(),
    )
}

fn strict_error(
    source: &str,
    assumptions: &[(&str, &str)],
) -> ComputeError {
    let session = Session::new();
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    session
        .compute(term, &config(assumptions))
        .expect_err(&format!("strict mode must reject `{source}`"))
}

#[test]
fn special_functions_are_dimensionless() {
    assert_eq!(dimension("gamma(2)", &[]), "1");
    assert_eq!(dimension("besselj(0, 2)", &[]), "1");
    assert_eq!(dimension("meijerg(1, 0, list(1), list(0), 2)", &[]), "1");
}

#[test]
fn special_function_rejects_a_dimensionful_argument() {
    let error = strict_error("gamma(x)", &[("x", "m")]);
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}

#[test]
fn integral_transforms_carry_the_dual_dimension() {
    // `L{t^3}(s) = 6/s^4`; with `s` the inverse of `t` this is a length^4.
    assert_eq!(
        dimension("laplace(t^3, t, s)", &[("t", "m"), ("s", "1/m")]),
        "Length^4"
    );
    // A Fourier transform of a dimensionless function stays dimensionless.
    assert_eq!(
        dimension("fourier(1/(t^2 + 4), t, w)", &[("t", "1"), ("w", "1")]),
        "1"
    );
    // The z transform sums over its index and keeps the dimension of `f`.
    assert_eq!(dimension("ztransform(1, n, z)", &[]), "1");
}

#[test]
fn linear_algebra_dimensions() {
    assert_eq!(
        dimension("det(list(list(x, 0), list(0, x)))", &[("x", "m")]),
        "Length^2"
    );
    assert_eq!(
        dimension("inverse(list(list(x)))", &[("x", "m")]),
        "1/Length"
    );
    assert_eq!(
        dimension(
            "matmul(list(list(x)), list(list(y)))",
            &[("x", "m"), ("y", "s")]
        ),
        "Length*Time"
    );
    assert_eq!(
        dimension("eigenvals(list(list(1, 0), list(0, 2)))", &[]),
        "1"
    );
}

#[test]
fn determinant_of_a_mixed_matrix_is_diagnosed() {
    // A matrix whose entries have different dimensions is not a scalar
    // conflict: the determinant still has the product dimension, and the
    // mixed list is reported as a diagnostic.
    let (dimension, diagnostics) = dimension_advisory(
        "det(list(list(x, 0), list(0, y)))",
        &[("x", "m"), ("y", "s")],
    );
    assert_eq!(dimension, "Length*Time");
    assert!(diagnostics >= 1);
}

#[test]
fn tensor_addition_rejects_mixed_dimensions() {
    let error = strict_error("tensor_add(list(x), list(y))", &[("x", "m"), ("y", "s")]);
    assert!(
        matches!(error, ComputeError::Dimension(_)),
        "unexpected error: {error}"
    );
}

#[test]
fn descriptive_statistics_compose_with_dimensions() {
    assert_eq!(
        dimension("mean(list(x, y))", &[("x", "m"), ("y", "m")]),
        "Length"
    );
    assert_eq!(
        dimension("variance(list(x, y))", &[("x", "m"), ("y", "m")]),
        "Length^2"
    );
    assert_eq!(
        dimension("std(list(x, y))", &[("x", "m"), ("y", "m")]),
        "Length"
    );
    // A density is a value per unit of its variable.
    assert_eq!(
        dimension(
            "pdf(normal(mu, sigma), t)",
            &[("mu", "m"), ("sigma", "m"), ("t", "m")]
        ),
        "1/Length"
    );
}

#[test]
fn statistics_advisory_reports_without_failing() {
    let (dimension, diagnostics) =
        dimension_advisory("mean(list(x, y))", &[("x", "m"), ("y", "m")]);
    assert_eq!(dimension, "Length");
    assert_eq!(diagnostics, 0);
}

#[test]
fn number_theory_and_combinatorics_are_dimensionless() {
    assert_eq!(dimension("factorint(12)", &[]), "1");
    assert_eq!(dimension("divisor_sum(12)", &[]), "1");
    assert_eq!(dimension("binomial(5, 2)", &[]), "1");
    assert_eq!(dimension("catalan(5)", &[]), "1");
}

#[test]
fn discrete_invariants_are_dimensionless() {
    assert_eq!(dimension("graph_node_count(graph_complete(4))", &[]), "1");
    assert_eq!(
        dimension("graph_chromatic_number(graph_cycle(5))", &[]),
        "1"
    );
    assert_eq!(dimension("group_order(dihedral_group(4))", &[]), "1");
    // A distance keeps the dimension of the points.
    assert_eq!(
        dimension(
            "euclidean_distance(list(x, y), list(0, 0))",
            &[("x", "m"), ("y", "m")]
        ),
        "Length"
    );
}

#[test]
fn lie_and_geometry_values() {
    assert_eq!(dimension("dim_irrep(A2, list(1, 0))", &[]), "1");
    assert_eq!(
        dimension("tensor_rank(list(list(1, 2), list(3, 4)))", &[]),
        "1"
    );
    // A rank is a count, even when the tensor entries carry a dimension.
    assert_eq!(
        dimension("tensor_rank(list(list(x, x), list(x, x)))", &[("x", "m")]),
        "1"
    );
}

#[test]
fn lie_bracket_multiplies_dimensions() {
    // The bracket is bilinear, so the dimension rule is `Product` and the
    // request is dimensionally coherent for the declared generators.
    assert_eq!(
        dimension(
            "lie_bracket(list(a, 0), list(0, b))",
            &[("a", "1"), ("b", "1")]
        ),
        "1"
    );
}
