//! Multivariate series through the public `Session` API.
//!
//! The tests check the shape of the expansions, their agreement with the
//! univariate `taylor`/`laurent` operators and, most importantly, that every
//! returned series approximates the function at small points with the error
//! scaling of the first omitted degree.

use rssn::prelude::*;

/// Reduces `src` symbolically and returns the printed answer and whether the
/// request was fully reduced.
fn run(src: &str) -> (String, bool) {
    let session = Session::new();
    let term = session
        .parse(src)
        .unwrap_or_else(|e| panic!("parse `{src}`: {e}"));
    let answer = session
        .compute(term, &Config::new())
        .unwrap_or_else(|e| panic!("compute `{src}`: {e}"));
    (answer.term.to_string(), answer.reduced)
}

/// The printed, fully reduced answer; panics when the request stayed
/// unevaluated.
fn series(src: &str) -> String {
    let (text, reduced) = run(src);
    assert!(reduced, "`{src}` was not reduced: {text}");
    text
}

/// Computes `request` symbolically and compares the value of the answer with
/// the value of `function` at `at`.
fn check(
    request: &str,
    function: &str,
    at: &[(&str, f64)],
    tolerance: f64,
) {
    let session = Session::new();
    let node = session
        .parse(request)
        .unwrap_or_else(|e| panic!("parse `{request}`: {e}"));
    let answer = session
        .compute(node, &Config::new())
        .unwrap_or_else(|e| panic!("compute `{request}`: {e}"));
    assert!(answer.reduced, "`{request}` was not reduced");
    let exact = session
        .parse(function)
        .unwrap_or_else(|e| panic!("parse `{function}`: {e}"));
    let got = answer
        .term
        .eval(at)
        .unwrap_or_else(|| panic!("cannot evaluate `{request}` at {at:?}"));
    let want = exact
        .eval(at)
        .unwrap_or_else(|| panic!("cannot evaluate `{function}` at {at:?}"));
    assert!(
        (got - want).abs() <= tolerance * (1.0 + want.abs()),
        "`{request}` = {} = {got}, but `{function}` = {want} at {at:?}",
        answer.term
    );
}

/// Two-dimensional points at which a truncated expansion is expected to be
/// accurate.
const SMALL: &[(&str, f64)] = &[("x", 0.03), ("y", -0.02)];

#[test]
fn taylor_expansions_are_exact() {
    assert_eq!(
        series("mseries(exp(x + y), list(x, y), list(0, 0), 3)"),
        "1/6*x^3 + 1/2*x^2*y + 1/2*x*y^2 + 1/6*y^3 + 1/2*x^2 + x*y + 1/2*y^2 + x + y + 1"
    );
    assert_eq!(
        series("mseries(1/(1 - x - y), list(x, y), list(0, 0), 3)"),
        "x^3 + 3*x^2*y + 3*x*y^2 + y^3 + x^2 + 2*x*y + y^2 + x + y + 1"
    );
    assert_eq!(
        series("mseries((1 + x + y)^(1/2), list(x, y), list(0, 0), 3)"),
        "1/16*x^3 + 3/16*x^2*y + 3/16*x*y^2 + 1/16*y^3 - 1/8*x^2 - 1/4*x*y - 1/8*y^2 + 1/2*x + 1/2*y + 1"
    );
    assert_eq!(
        series("mseries(ln(1 + x + y), list(x, y), list(0, 0), 3)"),
        "1/3*x^3 + x^2*y + x*y^2 + 1/3*y^3 - 1/2*x^2 - x*y - 1/2*y^2 + x + y"
    );
    assert_eq!(
        series("mseries(sin(x*y), list(x, y), list(0, 0), 6)"),
        "x*y - 1/6*x^3*y^3"
    );
    assert_eq!(
        series("mseries(cos(x + y), list(x, y), list(0, 0), 4)"),
        "1/24*x^4 + 1/6*x^3*y + 1/4*x^2*y^2 + 1/6*x*y^3 + 1/24*y^4 - 1/2*x^2 - x*y - 1/2*y^2 + 1"
    );
    assert_eq!(
        series("mseries(exp(x)/(1 - y), list(x, y), list(0, 0), 3)"),
        "1/6*x^3 + 1/2*x^2*y + x*y^2 + y^3 + 1/2*x^2 + x*y + y^2 + x + y + 1"
    );
}

#[test]
fn elementary_functions_compose() {
    assert_eq!(
        series("mseries(tan(x + y), list(x, y), list(0, 0), 3)"),
        "1/3*x^3 + x^2*y + x*y^2 + 1/3*y^3 + x + y"
    );
    assert_eq!(
        series("mseries(atan(x + y), list(x, y), list(0, 0), 3)"),
        "x - 1/3*x^3 - x^2*y - x*y^2 - 1/3*y^3 + y"
    );
    assert_eq!(
        series("mseries(sinh(x + y), list(x, y), list(0, 0), 3)"),
        "1/6*x^3 + 1/2*x^2*y + 1/2*x*y^2 + 1/6*y^3 + x + y"
    );
    assert_eq!(
        series("mseries(cosh(x + y), list(x, y), list(0, 0), 4)"),
        "1/24*x^4 + 1/6*x^3*y + 1/4*x^2*y^2 + 1/6*x*y^3 + 1/24*y^4 + 1/2*x^2 + x*y + 1/2*y^2 + 1"
    );
    assert_eq!(
        series("mseries(sin(x + y), list(x, y), list(0, 0), 5)"),
        "1/120*x^5 + 1/24*x^4*y + 1/12*x^3*y^2 + 1/12*x^2*y^3 + 1/24*x*y^4 + 1/120*y^5 - 1/6*x^3 - 1/2*x^2*y - 1/2*x*y^2 - 1/6*y^3 + x + y"
    );
    assert_eq!(
        series("mseries(exp(x)*cos(y), list(x, y), list(0, 0), 3)"),
        "1/6*x^3 - 1/2*x*y^2 + 1/2*x^2 - 1/2*y^2 + x + 1"
    );
    assert_eq!(
        series("mseries(1/(2 + x + y), list(x, y), list(0, 0), 2)"),
        "1/8*x^2 + 1/4*x*y + 1/8*y^2 - 1/4*x - 1/4*y + 1/2"
    );
    // `atan`, `asin` and `acos` are recovered by differentiating and
    // integrating; check them numerically instead of by shape.
    for (request, function) in [
        (
            "mseries(atan(x + y), list(x, y), list(0, 0), 6)",
            "atan(x + y)",
        ),
        (
            "mseries(asin(x + y), list(x, y), list(0, 0), 6)",
            "asin(x + y)",
        ),
        (
            "mseries(acos(x + y), list(x, y), list(0, 0), 6)",
            "acos(x + y)",
        ),
    ] {
        check(request, function, SMALL, 1e-8);
    }
}

#[test]
fn rational_functions_and_cancellation() {
    assert_eq!(
        series("mseries((x^2 - y^2)/(x - y), list(x, y), list(0, 0), 2)"),
        "x + y",
        "the removable singularity cancels"
    );
    assert_eq!(
        series("mseries(exp(x)*y^2, list(x, y), list(0, 0), 3)"),
        "x*y^2 + y^2"
    );
    check(
        "mseries(1/(1 - x - y), list(x, y), list(0, 0), 5)",
        "1/(1 - x - y)",
        &[("x", 0.05), ("y", 0.03)],
        1e-6,
    );
    check(
        "mseries((1 + x + y)^(1/2), list(x, y), list(0, 0), 5)",
        "sqrt(1 + x + y)",
        SMALL,
        1e-8,
    );
}

#[test]
fn total_degree_truncation() {
    // No term above the order: the order-1 truncation keeps only the linear
    // part, and a purely higher-order function truncates to zero.
    assert_eq!(
        series("mseries(exp(x + y), list(x, y), list(0, 0), 1)"),
        "x + y + 1"
    );
    assert_eq!(
        series("mseries(exp(x + y), list(x, y), list(0, 0), 0)"),
        "1"
    );
    assert_eq!(series("mseries(x^3 + x*y, list(x, y), list(0, 0), 1)"), "0");
    // The error of the geometric series scales like the first omitted degree.
    for order in 1..=4 {
        let request = format!("mseries(1/(1 - x - y), list(x, y), list(0, 0), {order})");
        let session = Session::new();
        let term = session
            .parse(&request)
            .unwrap_or_else(|e| panic!("parse `{request}`: {e}"));
        let answer = session
            .compute(term, &Config::new())
            .unwrap_or_else(|e| panic!("compute `{request}`: {e}"));
        let at = [("x", 0.04), ("y", 0.03)];
        let got = answer
            .term
            .eval(&at)
            .unwrap_or_else(|| panic!("cannot evaluate the series"));
        let want = 1.0 / (1.0 - 0.04 - 0.03);
        let bound = 2.0 * 0.1_f64.powi(i32::try_from(order + 1).unwrap_or(2));
        assert!(
            (got - want).abs() <= bound,
            "order {order}: {got} vs {want} (error {})",
            (got - want).abs()
        );
    }
    // Mixed terms of the total-degree boundary are present.
    let text = series("mseries(exp(x + y), list(x, y), list(0, 0), 3)");
    for mixed in ["x^2*y", "x*y^2"] {
        assert!(
            text.contains(mixed),
            "missing the mixed term {mixed}: {text}"
        );
    }
}

#[test]
fn laurent_series_with_monomial_factor() {
    assert_eq!(
        series("mseries(1/(x*y), list(x, y), list(0, 0), 2)"),
        "1/(x*y)"
    );
    assert_eq!(
        series("mseries(exp(x)/y^2, list(x, y), list(0, 0), 3)"),
        "1/120*x^5/y^2 + 1/24*x^4/y^2 + 1/6*x^3/y^2 + 1/2*x^2/y^2 + x/y^2 + 1/y^2"
    );
    assert_eq!(
        series("mseries(x^(-1)*sin(y), list(x, y), list(0, 0), 3)"),
        "y/x - 1/6*y^3/x"
    );
    // Multiplication back recovers the truncated numerator: y^2 times the
    // series above is the Taylor polynomial of exp(x) through x^5.
    let session = Session::new();
    let request = session
        .parse("mseries(exp(x)/y^2, list(x, y), list(0, 0), 3)")
        .unwrap_or_else(|e| panic!("parse: {e}"));
    let answer = session
        .compute(request, &Config::new())
        .unwrap_or_else(|e| panic!("compute: {e}"));
    let y = session.sym("y");
    let product = y.pow(2) * answer.term;
    let taylor = session
        .parse("taylor(exp(x), x, 0, 5)")
        .unwrap_or_else(|e| panic!("parse: {e}"));
    let taylor = session
        .compute(taylor, &Config::new())
        .unwrap_or_else(|e| panic!("compute: {e}"));
    for at in [0.2_f64, -0.3] {
        let got = product
            .eval(&[("x", at), ("y", 0.7)])
            .unwrap_or_else(|| panic!("cannot evaluate the product"));
        let want = taylor
            .term
            .eval(&[("x", at)])
            .unwrap_or_else(|| panic!("cannot evaluate the Taylor polynomial"));
        assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    }
}

#[test]
fn directional_laurent_series() {
    // `1/(x^2 - y^2)` has no Laurent expansion at the origin. The result is
    // the expansion in `x/y` through the requested order:
    //   -sum_{k=0}^{floor(N/2)} x^(2k) y^(-2-2k).
    assert_eq!(
        series("mseries(1/(x^2 - y^2), list(x, y), list(0, 0), 6)"),
        "-x^6/y^8 - x^4/y^6 - x^2/y^4 - 1/y^2"
    );
    // `(x + y)^(-2) = sum_p (-1)^p (p+1) x^p y^(-2-p)`.
    assert_eq!(
        series("mseries((x + y)^(-2), list(x, y), list(0, 0), 4)"),
        "5*x^4/y^6 - 4*x^3/y^5 + 3*x^2/y^4 - 2*x/y^3 + 1/y^2"
    );
    // Multiplying the truncated expansion back telescopes: with N = 6 the
    // remainder is exactly (x/y)^8.
    let session = Session::new();
    let request = session
        .parse("mseries(1/(x^2 - y^2), list(x, y), list(0, 0), 6)")
        .unwrap_or_else(|e| panic!("parse: {e}"));
    let answer = session
        .compute(request, &Config::new())
        .unwrap_or_else(|e| panic!("compute: {e}"));
    let x = session.sym("x");
    let y = session.sym("y");
    let product = (x.pow(2) - y.pow(2)) * answer.term;
    for (ratio, tolerance) in [(0.01_f64, 1e-12), (0.05, 1e-9)] {
        let at = [("x", ratio * 0.1), ("y", 0.1)];
        let got = product
            .eval(&at)
            .unwrap_or_else(|| panic!("cannot evaluate the product"));
        let want = 1.0 - ratio.powi(8);
        assert!(
            (got - want).abs() <= tolerance,
            "ratio {ratio}: {got} vs {want}"
        );
    }
    check(
        "mseries(1/(x^2 - y^2), list(x, y), list(0, 0), 8)",
        "1/(x^2 - y^2)",
        &[("x", 0.002), ("y", 0.1)],
        1e-6,
    );
    check(
        "mseries((x + y)^(-2), list(x, y), list(0, 0), 8)",
        "1/(x + y)^2",
        &[("x", 0.002), ("y", 0.1)],
        1e-6,
    );
}

#[test]
fn points_finite_and_infinite() {
    assert_eq!(
        series("mseries(1/(x*y), list(x, y), list(1, 2), 2)"),
        "1/8*(y - 2)^2 + 1/4*(y - 2)*(x - 1) + 1/2*(x - 1)^2 - 1/4*(y - 2) - 1/2*(x - 1) + 1/2"
    );
    assert_eq!(
        series("mseries(x/(x - 1), list(x), list(oo), 3)"),
        "1 + 1/x + 1/x^2 + 1/x^3"
    );
    assert_eq!(
        series("mseries(1/(x + 1), list(x), list(-oo), 2)"),
        "1/x - 1/x^2"
    );
    assert_eq!(
        series("mseries(exp(x)/y, list(x, y), list(0, oo), 2)"),
        "x/y + 1/y"
    );
    check(
        "mseries(1/(x*y), list(x, y), list(1, 2), 4)",
        "1/(x*y)",
        &[("x", 1.05), ("y", 1.95)],
        1e-6,
    );
    check(
        "mseries(1/(x^2 - y^2), list(x, y), list(1, 2), 2)",
        "1/(x^2 - y^2)",
        &[("x", 0.95), ("y", 1.9)],
        1e-3,
    );
    check(
        "mseries(exp(x)/(1 - y), list(x, y), list(1, 2), 3)",
        "exp(x)/(1 - y)",
        &[("x", 0.98), ("y", 2.03)],
        1e-5,
    );
    check(
        "mseries(sqrt(x), list(x), list(1), 4)",
        "sqrt(x)",
        &[("x", 1.1)],
        1e-5,
    );
}

#[test]
fn univariate_agreement_and_low_orders() {
    for (request, univariate) in [
        (
            "mseries(exp(x), list(x), list(0), 4)",
            "taylor(exp(x), x, 0, 4)",
        ),
        (
            "mseries(sin(x), list(x), list(0), 5)",
            "taylor(sin(x), x, 0, 5)",
        ),
        (
            "mseries(ln(x), list(x), list(1), 3)",
            "taylor(ln(x), x, 1, 3)",
        ),
        (
            "mseries(1/(x*(1 - x)), list(x), list(0), 2)",
            "laurent(1/(x*(1 - x)), x, 0, 2)",
        ),
        (
            "mseries(exp(x)/x^2, list(x), list(0), 1)",
            "laurent(exp(x)/x^2, x, 0, 1)",
        ),
        (
            "mseries(x/(x - 1), list(x), list(oo), 3)",
            "laurent(x/(x - 1), x, oo, 3)",
        ),
        (
            "mseries(sin(x)/x, list(x), list(0), 4)",
            "taylor(sin(x)/x, x, 0, 4)",
        ),
    ] {
        assert_eq!(series(request), series(univariate), "{request}");
    }
    assert_eq!(series("mseries(exp(x), list(x), list(0), 0)"), "1");
    assert_eq!(series("mseries(exp(x), list(x), list(0), 1)"), "x + 1");
}

#[test]
fn symbolic_parameters() {
    assert_eq!(
        series("mseries(exp(x), list(x), list(a), 2)"),
        "1/2*(x - a)^2*exp(a) + (x - a)*exp(a) + exp(a)"
    );
    assert_eq!(
        series("mseries(a*exp(x + y), list(x, y), list(0, 0), 2)"),
        "1/2*a*x^2 + a*x*y + 1/2*a*y^2 + a*x + a*y + a"
    );
    check(
        "mseries(exp(a*x + y), list(x, y), list(0, 0), 4)",
        "exp(a*x + y)",
        &[("x", 0.02), ("y", -0.03), ("a", 1.7)],
        1e-8,
    );
    check(
        "mseries(sin(a*x + y), list(x, y), list(0, 0), 5)",
        "sin(a*x + y)",
        &[("x", 0.02), ("y", -0.03), ("a", 0.6)],
        1e-8,
    );
}

#[test]
fn mtaylor_is_the_taylor_special_case() {
    assert_eq!(
        series("mtaylor(exp(x + y), list(x, y), list(0, 0), 2)"),
        series("mseries(exp(x + y), list(x, y), list(0, 0), 2)")
    );
    for src in [
        "mtaylor(1/(x*y), list(x, y), list(0, 0), 2)",
        "mtaylor(exp(x)/y^2, list(x, y), list(0, 0), 2)",
        "mtaylor(exp(x), list(x), list(oo), 2)",
        "mtaylor(sqrt(x), list(x), list(0), 2)",
    ] {
        let (text, reduced) = run(src);
        assert!(!reduced, "`{src}` should stay unevaluated: {text}");
    }
    assert_eq!(
        series("mtaylor(1/(1 - x - y), list(x, y), list(0, 0), 2)"),
        "x^2 + 2*x*y + y^2 + x + y + 1"
    );
}

#[test]
fn unevaluated_cases_stay_unevaluated() {
    for src in [
        // Essential singularities.
        "mseries(exp(1/x), list(x), list(0), 2)",
        "mseries(sin(1/x), list(x), list(0), 2)",
        "mseries(exp(x/y), list(x, y), list(0, 0), 2)",
        "mseries(sin(1/y), list(x, y), list(0, 0), 2)",
        // A logarithmic singularity.
        "mseries(ln(x), list(x), list(0), 2)",
        // A branch locus through the origin: no Puiseux expansion exists.
        "mseries((x^2 + y^2)^(1/2), list(x, y), list(0, 0), 4)",
        // Malformed requests.
        "mseries(exp(x), x, list(0), 2)",
        "mseries(exp(x), list(x), list(0), -1)",
        "mseries(exp(x), list(x, y), list(0), 2)",
    ] {
        let (text, reduced) = run(src);
        assert!(!reduced, "`{src}` should stay unevaluated: {text}");
        assert!(
            text.contains("mseries("),
            "`{src}` must not be mis-expanded: {text}"
        );
    }
}

#[test]
fn puiseux_monomial_factors_are_exact() {
    assert_eq!(
        series("mseries(sqrt(x*y), list(x, y), list(0, 0), 4)"),
        "x^(1/2)*y^(1/2)"
    );
    assert_eq!(series("mseries(x^(4/3), list(x), list(0), 3)"), "x^(4/3)");
    assert_eq!(
        series("mseries((x*y)^(5/2), list(x, y), list(0, 0), 5)"),
        "x^(5/2)*y^(5/2)"
    );
    assert_eq!(
        series("mseries(x^(1/2)*exp(y), list(x, y), list(0, 0), 3)"),
        "1/2*x^(1/2)*y^2 + x^(1/2)*y + x^(1/2)"
    );
    assert_eq!(
        series("mseries(sqrt(x^2*(1 + y)), list(x, y), list(0, 0), 3)"),
        "1/2*x*y - 1/8*x*y^2 + x"
    );
    assert_eq!(
        series("mseries(x^(1/2)/(1 - x - y), list(x, y), list(0, 0), 2)"),
        "x^(3/2) + x^(1/2)*y + x^(1/2)"
    );
}

/// Simplifies `src` in a fresh session and returns the printed result.
fn simplify(src: &str) -> String {
    let session = Session::new();
    let term = session
        .parse(src)
        .unwrap_or_else(|e| panic!("parse `{src}`: {e}"));
    session
        .simplify(term)
        .unwrap_or_else(|e| panic!("simplify `{src}`: {e}"))
        .to_string()
}

/// The numeric value of `src` at `at`, without running any rules.
fn value(
    src: &str,
    at: &[(&str, f64)],
) -> f64 {
    let session = Session::new();
    session
        .parse(src)
        .unwrap_or_else(|e| panic!("parse `{src}`: {e}"))
        .eval(at)
        .unwrap_or_else(|| panic!("cannot evaluate `{src}` at {at:?}"))
}

#[test]
fn puiseux_expansions_raise_back() {
    // Repeated multiplication and powers recover the monomial part exactly.
    assert_eq!(simplify("x^(1/2)*y^(1/2)*x^(1/2)*y^(1/2)"), "x*y");
    assert_eq!(simplify("(x^(4/3))^3"), "x^4");
    assert_eq!(simplify("x^(5/2)*y^(5/2)*x^(5/2)*y^(5/2)"), "x^5*y^5");
    // Dividing the series by its monomial factor recovers the unit.
    let at = &[("x", 0.04), ("y", 0.03)][..];
    let got = value("(1/2*x^(1/2)*y^2 + x^(1/2)*y + x^(1/2))*x^(-1/2)", at);
    let want = value("exp(y)", at);
    assert!((got - want).abs() < 1e-4, "{got} vs {want}");
    // Multiplying by the denominator recovers the truncation of x^(1/2):
    // the remainder has total degree above the requested order.
    let got = value("(1 - x - y)*(x^(3/2) + x^(1/2)*y + x^(1/2))", at);
    let want = value("x^(1/2)", at);
    assert!((got - want).abs() < 1e-2, "{got} vs {want}");
}

#[test]
fn puiseux_expansions_match_numerically() {
    // Positive small points: the principal branch of every fractional power.
    let at = &[("x", 0.04), ("y", 0.09)][..];
    for (request, function, tolerance) in [
        (
            "mseries(sqrt(x*y), list(x, y), list(0, 0), 4)",
            "sqrt(x*y)",
            1e-15,
        ),
        (
            "mseries((x*y)^(5/2), list(x, y), list(0, 0), 5)",
            "(x*y)^(5/2)",
            1e-15,
        ),
        (
            "mseries(x^(1/2)*exp(y), list(x, y), list(0, 0), 3)",
            "x^(1/2)*exp(y)",
            1e-4,
        ),
        (
            "mseries(sqrt(x^2*(1 + y)), list(x, y), list(0, 0), 3)",
            "sqrt(x^2*(1 + y))",
            1e-4,
        ),
        (
            "mseries(x^(1/2)/(1 - x - y), list(x, y), list(0, 0), 4)",
            "x^(1/2)/(1 - x - y)",
            1e-3,
        ),
    ] {
        check(request, function, at, tolerance);
    }
    // Univariate fractional powers.
    let session = Session::new();
    let request = session
        .parse("mseries(x^(4/3), list(x), list(0), 3)")
        .unwrap_or_else(|e| panic!("parse: {e}"));
    let answer = session
        .compute(request, &Config::new())
        .unwrap_or_else(|e| panic!("compute: {e}"));
    let got = answer
        .term
        .eval(&[("x", 0.008)])
        .unwrap_or_else(|| panic!("cannot evaluate the series"));
    let want = 0.008_f64.powf(4.0 / 3.0);
    assert!((got - want).abs() < 1e-15, "{got} vs {want}");
}

#[test]
fn numeric_tolerance_profile() {
    // (request, function, order): the error must be bounded by the size of
    // the first omitted terms, not merely "small". At `|x|,|y| <= 0.04` the
    // radius is at most 0.08, so an order-N error of about radius^(N+1),
    // with a generous constant, is a meaningful profile.
    let analytic = [
        (
            "mseries(exp(x + y), list(x, y), list(0, 0), 5)",
            "exp(x + y)",
            5,
        ),
        (
            "mseries(sin(x + y), list(x, y), list(0, 0), 7)",
            "sin(x + y)",
            7,
        ),
        (
            "mseries(cos(x*y), list(x, y), list(0, 0), 6)",
            "cos(x*y)",
            6,
        ),
        (
            "mseries(1/(1 - x - y), list(x, y), list(0, 0), 5)",
            "1/(1 - x - y)",
            5,
        ),
        (
            "mseries(ln(1 + x + y), list(x, y), list(0, 0), 5)",
            "ln(1 + x + y)",
            5,
        ),
        (
            "mseries((1 + x + y)^(1/2), list(x, y), list(0, 0), 5)",
            "sqrt(1 + x + y)",
            5,
        ),
        (
            "mseries(exp(x)/(1 - y), list(x, y), list(0, 0), 5)",
            "exp(x)/(1 - y)",
            5,
        ),
        (
            "mseries(tan(x + y), list(x, y), list(0, 0), 5)",
            "tan(x + y)",
            5,
        ),
        (
            "mseries(atan(x + y), list(x, y), list(0, 0), 7)",
            "atan(x + y)",
            7,
        ),
        (
            "mseries(sinh(x + y), list(x, y), list(0, 0), 5)",
            "sinh(x + y)",
            5,
        ),
        (
            "mseries(cosh(x + y), list(x, y), list(0, 0), 6)",
            "cosh(x + y)",
            6,
        ),
        (
            "mseries(exp(x)*cos(y), list(x, y), list(0, 0), 5)",
            "exp(x)*cos(y)",
            5,
        ),
    ];
    for (request, function, order) in analytic {
        for at in [
            &[("x", 0.04), ("y", -0.03)][..],
            &[("x", -0.02), ("y", 0.04)][..],
            &[("x", 0.01), ("y", 0.01)][..],
            &[("x", 0.035), ("y", 0.02)][..],
        ] {
            let radius: f64 = at.iter().map(|&(_, value)| f64::abs(value)).sum();
            #[allow(clippy::cast_possible_wrap)]
            let tolerance = 20.0 * radius.powi(order + 1);
            check(request, function, at, tolerance);
        }
    }
    // A truncated Taylor series with the same order is a strictly worse
    // approximation than the next order: the error profile is monotone.
    for (low, high) in [(2_u32, 4_u32), (4, 6)] {
        let at = [("x", 0.05), ("y", -0.04)];
        let value = |order: u32| {
            let session = Session::new();
            let request = format!("mseries(exp(x + y), list(x, y), list(0, 0), {order})");
            let term = session
                .parse(&request)
                .unwrap_or_else(|e| panic!("parse `{request}`: {e}"));
            let answer = session
                .compute(term, &Config::new())
                .unwrap_or_else(|e| panic!("compute `{request}`: {e}"));
            answer
                .term
                .eval(&at)
                .unwrap_or_else(|| panic!("cannot evaluate the series"))
        };
        let want = (0.05_f64 - 0.04).exp();
        let coarse = (value(low) - want).abs();
        let fine = (value(high) - want).abs();
        assert!(
            fine < coarse * 0.05,
            "order {high} error {fine} is not much smaller than order {low} error {coarse}"
        );
    }
    // Laurent expansions with a monomial factor are far more accurate at
    // small points.
    for (request, function, at, tolerance) in [
        (
            "mseries(exp(x)/y^2, list(x, y), list(0, 0), 4)",
            "exp(x)/y^2",
            &[("x", 0.01), ("y", 0.5)][..],
            1e-9,
        ),
        (
            "mseries(x^(-1)*sin(y), list(x, y), list(0, 0), 5)",
            "sin(y)/x",
            &[("x", 0.4), ("y", 0.02)][..],
            1e-9,
        ),
        (
            "mseries(1/(x*y), list(x, y), list(0, 0), 2)",
            "1/(x*y)",
            &[("x", 0.3), ("y", -0.7)][..],
            0.0,
        ),
    ] {
        check(request, function, at, tolerance);
    }
}

#[test]
fn three_or_more_variables() {
    // A three-variable Taylor expansion: every monomial of total degree at
    // most 3 appears exactly once.
    let text = series("mseries(exp(x + y + z), list(x, y, z), list(0, 0, 0), 3)");
    for monomial in [
        "x^3", "x^2*y", "x^2*z", "x*y^2", "x*y*z", "x*z^2", "y^3", "y^2*z", "y*z^2", "z^3", "x^2",
        "x*y", "x*z", "y^2", "y*z", "z^2", "x", "y", "z",
    ] {
        assert!(text.contains(monomial), "missing {monomial}: {text}");
    }
    check(
        "mseries(exp(x + y + z), list(x, y, z), list(0, 0, 0), 5)",
        "exp(x + y + z)",
        &[("x", 0.02), ("y", -0.03), ("z", 0.01)],
        1e-9,
    );
    // Laurent monomials in three variables, and a three-variable directional
    // expansion.
    assert_eq!(
        series("mseries(1/(x*y*z), list(x, y, z), list(0, 0, 0), 2)"),
        "1/(x*y*z)"
    );
    assert_eq!(
        series("mseries(1/(x*y^2), list(x, y), list(0, 0), 2)"),
        "1/(x*y^2)"
    );
    check(
        "mseries(1/(x*y*z), list(x, y, z), list(0, 0, 0), 2)",
        "1/(x*y*z)",
        &[("x", 0.3), ("y", -0.2), ("z", 0.5)],
        0.0,
    );
    // Three variables at a finite non-zero point.
    check(
        "mseries(exp(x + 2*y - z), list(x, y, z), list(1, 1, 1), 4)",
        "exp(x + 2*y - z)",
        &[("x", 0.99), ("y", 1.01), ("z", 0.98)],
        1e-8,
    );
    // Mixed finite and infinite points.
    check(
        "mseries(exp(x)/y + z^2, list(x, y, z), list(0, oo, 0), 3)",
        "exp(x)/y + z^2",
        &[("x", 0.01), ("y", 40.0), ("z", 0.02)],
        1e-6,
    );
}
