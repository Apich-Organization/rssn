//! Exact inputs give exact answers: no float in the symbolic phase, closed
//! forms at known special values, odd/even normalisation, and printing.
mod common;
use common::*;
use rssn::api::Config;
use rssn::api::Session;

fn exact(
    src: &str,
    want: &str,
) {
    let o = sym(src);
    assert!(o.reduced, "{src} not reduced: {}", o.text);
    assert!(!has_float(&o.text), "{src} leaked a float: {}", o.text);
    assert_eq!(o.text, want, "{src}");
}

// ---- 1. exact limits

#[test]
fn limit_one_minus_cos() {
    exact("limit((1 - cos(t))/t^2, t, 0)", "1/2");
}

#[test]
fn limit_rational_infinity() {
    exact("limit((t^2+1)/(2*t^2 - t), t, oo)", "1/2");
}

/// A numeric target is allowed to answer with a float, a symbolic one is not.
#[test]
fn numeric_target_is_a_float_symbolic_is_not() {
    let n = num("limit((1 - cos(t))/t^2, t, 0)", 1e-9);
    assert!((n.value.unwrap_or(f64::NAN) - 0.5).abs() < 1e-6, "{n:?}");
    let s = sym("limit((1 - cos(t))/t^2, t, 0)");
    assert!(!has_float(&s.text), "{s:?}");
}

// ---- 2. battery of growth comparisons

#[test]
fn limit_battery() {
    exact("limit(t*ln(t), t, 0)", "0");
    exact("limit(t*ln(t), t, 0, plus)", "0");
    exact("limit(x^x, x, 0, plus)", "1");
    exact("limit((1 + 1/x)^x, x, oo)", "E");
    exact("limit(x*exp(-x), x, oo)", "0");
    exact("limit(ln(x)/x, x, oo)", "0");
    exact("limit(sin(x)/x, x, oo)", "0");
    exact("limit((exp(x) - 1)/x, x, 0)", "1");
    exact("limit(x^2*ln(x), x, 0)", "0");
    exact("limit(x^2*ln(x), x, 0, plus)", "0");
}

// ---- 3. special values

#[test]
fn special_values_are_exact() {
    exact("si(0)", "0");
    exact("erf(0)", "0");
    exact("gamma(1)", "1");
    exact("gamma(5)", "24");
    exact("gamma(1/2)", "pi^(1/2)");
    exact("zeta(2)", "1/6*pi^2");
    exact("zeta(0)", "-1/2");
    exact("li(0)", "0");
}

#[test]
fn si_integral_is_closed() {
    let o = sym("defint(sin(t)/t, t, 0, 1)");
    assert!(o.reduced, "{o:?}");
    assert!(!o.text.contains("si(0)"), "{}", o.text);
    assert_eq!(o.text, "si(1)");
}

// ---- 4. parity

#[test]
fn odd_even_functions_normalise() {
    exact("asin(-x)", "-asin(x)");
    exact("sin(-x)", "-sin(x)");
    exact("atan(-x)", "-atan(x)");
    exact("sinh(-x)", "-sinh(x)");
    exact("cos(-x)", "cos(x)");
    exact("cosh(-x)", "cosh(x)");
    exact("tan(-x)", "-tan(x)");
    exact("erf(-x)", "-erf(x)");
}

#[test]
fn antiderivative_of_semicircle() {
    exact(
        "integral(sqrt(1 - t^2), t)",
        "1/2*t*sqrt(1 - t^2) + 1/2*asin(t)",
    );
}

// ---- 5. printing

#[test]
fn derivative_prints_as_quotient() {
    let o = sym("diff(ln(t^2 + 1), t)");
    assert_eq!(o.text, "2*t/(t^2 + 1)");
    assert!(o.latex.contains("\\frac"), "{}", o.latex);
}

#[test]
fn negative_powers_in_every_printer() {
    let s = Session::new();
    let t = s
        .parse("2*x*(x^2 + 1)^(-1)")
        .unwrap_or_else(|e| panic!("{e}"));
    let a = s
        .compute(t, &Config::new())
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(a.term.to_string(), "2*x/(x^2 + 1)");
    assert!(
        a.term.to_latex().contains("\\frac"),
        "{}",
        a.term.to_latex()
    );
    assert!(a.term.to_typst().contains("frac("), "{}", a.term.to_typst());
}

#[test]
fn rational_coefficients_print_consistently() {
    // Documented rule: a rational coefficient is written in front, `1/120*x^5`.
    let o = sym("x^5/120");
    assert_eq!(o.text, "1/120*x^5");
}

// ---- 6. series order

#[test]
fn series_prints_in_ascending_order() {
    let o = sym("taylor(sin(t)/t, t, 0, 6)");
    assert!(o.reduced, "{o:?}");
    assert_eq!(o.text, "1 - 1/6*t^2 + 1/120*t^4 - 1/5040*t^6");
}

#[test]
fn polynomials_print_in_descending_order() {
    assert_eq!(sym("1 + x + x^2").text, "x^2 + x + 1");
}

// ---- assumptions: logs and powers

fn assuming(
    src: &str,
    facts: &[(&str, rssn::graph::Facts)],
) -> String {
    let mut config = Config::new();
    for (name, f) in facts {
        config = config.assume(name, *f);
    }
    run(src, &config).text
}

#[test]
fn assumptions_unlock_log_and_power_identities() {
    use rssn::graph::Facts;
    let positive = [("x", Facts::POSITIVE), ("y", Facts::POSITIVE)];
    assert_eq!(assuming("sqrt(x^2)", &positive), "x");
    assert_eq!(assuming("ln(x^2)", &positive), "2*ln(x)");
    assert_eq!(assuming("abs(x)", &positive), "x");
    assert_eq!(assuming("exp(ln(x))", &[]), "x");
    assert_eq!(assuming("sqrt(x^2)", &[("x", Facts::REAL)]), "abs(x)");
    assert_eq!(assuming("abs(x)", &[("x", Facts::NEGATIVE)]), "-x");
    // Without an assumption the complex-valued identities are not applied.
    assert_eq!(assuming("sqrt(x^2)", &[]), "sqrt(x^2)");
    assert_eq!(assuming("ln(x^2)", &[]), "ln(x^2)");
}

#[test]
fn polynomial_identities_expand() {
    assert_eq!(sym("(x + y)^2 - x^2 - 2*x*y - y^2").text, "0");
    assert_eq!(sym("(x + 1)^3 - x^3 - 3*x^2 - 3*x - 1").text, "0");
}

#[test]
fn polynomial_roots_never_leak_floats() {
    for src in [
        "solve(x^5 - x - 1, x)",
        "solve(x^4 - x - 1, x)",
        "solve(x^6 - x - 1, x)",
        "solve(x^7 - 3*x + 1, x)",
        "solve(x - cos(x), x)",
    ] {
        let o = sym(src);
        assert!(!has_float(&o.text), "{src}: {}", o.text);
    }
}

// ---- 7. multiplicities

#[test]
fn multiplicities_match_roots() {
    // A polynomial over the rationals reports one multiplicity per root; a
    // transcendental equation may report none ("not reported" is not the
    // claim "simple"), but never a list of a different length.
    let s = Session::new();
    for (src, polynomial) in [
        ("exp(t) - 2", false),
        ("t^2 - 2*t + 1", true),
        ("t^2 - 1", true),
        ("ln(t) - 1", false),
        ("sin(t) - 1/2", false),
        ("2*t + 3", true),
        ("t^3 - 2", true),
    ] {
        let eq = s.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let v = s.parse("t").unwrap_or_else(|e| panic!("{e}"));
        let sol = s
            .solve_structured(eq, v)
            .unwrap_or_else(|e| panic!("{src}: {e}"));
        if polynomial {
            assert_eq!(sol.roots.len(), sol.multiplicities.len(), "{src}: {sol:?}");
        } else {
            assert!(
                sol.multiplicities.is_empty() || sol.multiplicities.len() == sol.roots.len(),
                "{src}: {sol:?}"
            );
        }
    }
}
