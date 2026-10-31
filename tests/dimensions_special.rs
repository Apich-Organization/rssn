//! Dimension rules for the special, transform, linear-algebra, statistics,
//! number-theory, combinatorics, discrete, Lie, geometric-algebra and
//! geometry rule sets.
//!
//! Every test names a representative operator and asserts the dimension of
//! the answer in Strict mode (a conflict would be a hard error) or, for the
//! negative tests, that Strict mode does reject an ill-typed request.

use rssn::graph::DimensionMode;
use rssn::graph::Engine;
use rssn::graph::Graph;
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
fn matrix_functions_carry_the_degree_of_an_entry() {
    // `adj(A)` and `cof(A)` of a 3x3 matrix have 2x2 minors as entries.
    assert_eq!(
        dimension(
            "adjugate(list(list(x, 0, 0), list(0, x, 0), list(0, 0, x)))",
            &[("x", "m")]
        ),
        "Length^2"
    );
    assert_eq!(
        dimension(
            "cofactors(list(list(x, 0, 0), list(0, x, 0), list(0, 0, x)))",
            &[("x", "m")]
        ),
        "Length^2"
    );
    // A square root of a matrix of squared lengths has length entries.
    assert_eq!(
        dimension("matsqrt(list(list(quantity(4, m^2))))", &[]),
        "Length"
    );
    assert_eq!(
        dimension("cholesky(list(list(quantity(4, m^2))))", &[]),
        "Length"
    );
    // `A^n` has entries `entry^n`, for positive and negative powers.
    assert_eq!(
        dimension("mpow(list(list(quantity(2, m))), 3)", &[]),
        "Length^3"
    );
    assert_eq!(
        dimension("mpow(list(list(x)), -1)", &[("x", "m")]),
        "1/Length"
    );
    // Unit vectors are dimensionless.
    assert_eq!(
        dimension(
            "orthonormalize(list(list(x, 0), list(0, x)))",
            &[("x", "m")]
        ),
        "1"
    );
    // A non-homogeneous matrix function, or a matrix exponential, of a
    // dimensionful matrix is rejected rather than guessed.
    for source in [
        "matfun(exp(x), x, list(list(quantity(1, m))))",
        "matexp(list(list(quantity(1, m))), 4)",
    ] {
        let error = strict_error(source, &[]);
        assert!(
            matches!(error, ComputeError::Dimension(_)),
            "`{source}`: unexpected error {error}"
        );
    }
}

#[test]
fn row_reduction_returns_dimensionless_entries() {
    // RREF divides pivot rows by pivots, and a null-space basis normalises a
    // free variable to 1: their entries are ratios, hence dimensionless.
    // A full-rank matrix has the empty basis, whose element dimension is
    // unconstrained (like a literal zero).
    for source in [
        "rref(list(list(x, 0), list(0, x)))",
        "nullspace(list(list(x, 0), list(0, x)))",
        "nullspace(list(list(x, 0), list(0, 0)))",
        "rowspace(list(list(x, 0), list(0, x)))",
        "left_nullspace(list(list(x, 0), list(0, x)))",
    ] {
        assert_eq!(dimension(source, &[("x", "m")]), "1", "{source}");
    }
    // The column space keeps the entry dimension: the columns are not scaled.
    assert_eq!(
        dimension("colspace(list(list(x, 0), list(0, x)))", &[("x", "m")]),
        "Length"
    );
}

#[test]
fn distribution_parameters_are_consistent() {
    // A location-scale family: both parameters have the distribution's
    // dimension, and `expectation` agrees with the reduced closed form.
    assert_eq!(
        dimension("normal(quantity(1, m), quantity(2, m))", &[]),
        "Length"
    );
    assert_eq!(dimension("gamma_dist(2, quantity(1, 1/s))", &[]), "Time");
    assert_eq!(
        dimension("expectation(gamma_dist(2, quantity(1, 1/s)))", &[]),
        "Time"
    );
    assert_eq!(dimension("exponential(quantity(2, 1/s))", &[]), "Time");
    assert_eq!(
        dimension("weibull(quantity(2, 1), quantity(3, m))", &[]),
        "Length"
    );
    // A location in metres with a scale in seconds is ill-typed.
    for source in [
        "normal(quantity(1, m), quantity(2, s))",
        "cauchy(quantity(1, m), quantity(2, s))",
        "variance_of(normal(quantity(1, m), quantity(2, s)))",
    ] {
        let error = strict_error(source, &[]);
        assert!(
            matches!(error, ComputeError::Dimension(_)),
            "`{source}`: unexpected error {error}"
        );
    }
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

/// The printed parts of the reduced answer, each typed on its own.
fn expanded_part_dimensions(source: &str) -> Vec<String> {
    let session = Session::new();
    let term = session
        .parse(source)
        .unwrap_or_else(|error| panic!("cannot parse `{source}`: {error}"));
    let answer = session
        .compute(term, &Config::new().dimensions(DimensionMode::Advisory))
        .unwrap_or_else(|error| panic!("`{source}` failed: {error}"));
    let text = answer.term.to_string();
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard())
        .unwrap_or_else(|error| panic!("cannot install the standard rules: {error}"));
    let root = graph
        .parse(&text)
        .unwrap_or_else(|error| panic!("cannot re-parse the answer `{text}`: {error}"));
    graph
        .children(root)
        .iter()
        .map(|&child| graph.dimension(child).dimension.to_string())
        .collect()
}

/// Structured decompositions are non-scalar containers, not mixed scalar
/// lists: Advisory mode reports their dimension without a `MixedList`
/// finding, and the parts keep the documented dimensions.
#[test]
fn structured_decompositions_are_typed_containers() {
    let length2 = "list(list(quantity(2, m), 0), list(0, quantity(3, m)))";
    let time3 = "list(list(quantity(2, s), 0, 0), \
                 list(0, quantity(3, s), 0), list(0, 0, quantity(4, s)))";
    let unit2 = "list(list(1, 0), list(0, 2))";
    for source in [
        format!("lu({length2})"),
        format!("lu({time3})"),
        format!("lu({unit2})"),
        format!("qr({length2})"),
        format!("qr({time3})"),
        format!("svd({length2})"),
        format!("eigenvects({length2})"),
        format!("jordan({length2})"),
        format!("jordan_blocks({time3})"),
        format!("spectral({length2})"),
        format!("ldl({length2})"),
        format!("ldl({unit2})"),
    ] {
        let (dimension, diagnostics) = dimension_advisory(&source, &[]);
        assert_eq!(dimension, "non-metric", "`{source}`");
        assert_eq!(diagnostics, 0, "`{source}` produced {diagnostics} findings");
    }

    // Eigenvalues are entries of the matrix: they keep its dimension.
    assert_eq!(dimension(&format!("eigenvals({length2})"), &[]), "Length");
    assert_eq!(dimension(&format!("eigenvals({time3})"), &[]), "Time");

    // `lu` = [P, L, U]: the permutation and unit-triangular factor are
    // dimensionless, U carries the entry dimension.
    assert_eq!(
        expanded_part_dimensions(&format!("lu({length2})")),
        ["1", "1", "Length"]
    );
    assert_eq!(
        expanded_part_dimensions(&format!("lu({time3})")),
        ["1", "1", "Time"]
    );
    // `qr` = [Q, R]: the orthonormal factor is dimensionless, R keeps it.
    assert_eq!(
        expanded_part_dimensions(&format!("qr({length2})")),
        ["1", "Length"]
    );
    // `ldl` = [L, d]: L is dimensionless, the pivots keep the entry
    // dimension.
    assert_eq!(
        expanded_part_dimensions(&format!("ldl({length2})")),
        ["1", "Length"]
    );
    // `jordan` = [P, J]: the chain basis is normalised and dimensionless,
    // the diagonal of J holds the eigenvalues.
    assert_eq!(
        expanded_part_dimensions(&format!("jordan({length2})")),
        ["1", "Length"]
    );
}

/// Shape parameters are pure numbers; location/scale parameters agree with
/// the random variable; a density point lives in the variable's dimension.
#[test]
fn distribution_parameters_are_typed() {
    assert_eq!(dimension("gamma_dist(2, quantity(1, 1/s))", &[]), "Time");
    assert_eq!(
        dimension("weibull(quantity(2, 1), quantity(3, m))", &[]),
        "Length"
    );
    assert_eq!(
        dimension("pareto(quantity(2, m), quantity(3, 1))", &[]),
        "Length"
    );
    assert_eq!(
        dimension("normal(quantity(1, m), quantity(2, m))", &[]),
        "Length"
    );
    assert_eq!(
        dimension(
            "pdf(normal(quantity(1, m), quantity(2, m)), quantity(3, m))",
            &[]
        ),
        "1/Length"
    );
    assert_eq!(
        dimension(
            "cdf(normal(quantity(1, m), quantity(2, m)), quantity(3, m))",
            &[]
        ),
        "1"
    );
    for source in [
        // shape parameters must be dimensionless
        "gamma_dist(quantity(2, m), quantity(3, 1/s))",
        "weibull(quantity(2, m), quantity(3, s))",
        "pareto(quantity(2, m), quantity(3, s))",
        "beta_dist(quantity(2, m), quantity(3, 1))",
        "student_t(quantity(2, m))",
        "lognormal(quantity(1, m), quantity(2, 1))",
        // location and scale must agree
        "normal(quantity(1, m), quantity(2, s))",
        // the evaluation point must be the random variable
        "pdf(normal(quantity(1, m), quantity(2, m)), quantity(3, s))",
        "cdf(normal(quantity(1, m), quantity(2, m)), quantity(3, s))",
    ] {
        let error = strict_error(source, &[]);
        assert!(
            matches!(error, ComputeError::Dimension(_)),
            "`{source}`: unexpected error {error}"
        );
    }
}

/// A homogeneous matrix polynomial carries its degree; a non-homogeneous
/// one still requires a dimensionless matrix.
#[test]
fn homogeneous_matrix_functions_carry_their_degree() {
    let length2 = "list(list(quantity(2, m), 0), list(0, quantity(3, m)))";
    let time3 = "list(list(quantity(2, s), 0, 0), \
                 list(0, quantity(3, s), 0), list(0, 0, quantity(4, s)))";
    let unit2 = "list(list(1, 0), list(0, 2))";
    assert_eq!(
        dimension(&format!("matfun(x^2, x, {length2})"), &[]),
        "Length^2"
    );
    assert_eq!(
        dimension(&format!("matfun(x^2, x, {time3})"), &[]),
        "Time^2"
    );
    assert_eq!(
        dimension(&format!("matfun(x^3, x, {length2})"), &[]),
        "Length^3"
    );
    assert_eq!(
        dimension(&format!("matfun(quantity(2, m) * x^2, x, {length2})"), &[]),
        "Length^3"
    );
    // A constant function is independent of the matrix.
    assert_eq!(
        dimension(&format!("matfun(quantity(5, m), x, {length2})"), &[]),
        "Length"
    );
    // The dimensionless matrix stays dimensionless.
    assert_eq!(dimension(&format!("matfun(x^2, x, {unit2})"), &[]), "1");
    // Degrees 2 and 1 do not agree: a dimensionful matrix is refused.
    for source in [
        format!("matfun(x^2 + x, x, {length2})"),
        format!("matfun(exp(x), x, {length2})"),
    ] {
        let error = strict_error(&source, &[]);
        assert!(
            matches!(error, ComputeError::Dimension(_)),
            "`{source}`: unexpected error {error}"
        );
    }
}
