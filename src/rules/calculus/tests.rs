//! Tests of the calculus rule set as a whole.

use crate::rules::calculus;
use crate::rules::testing::eval;
use crate::rules::testing::numeric;
use crate::rules::testing::reduce_with;
use crate::rules::testing::simplify;

fn run(src: &str) -> String {
    simplify(&[calculus()], src)
}

/// The antiderivative the engine finds for `f` must differentiate back to
/// `f`: checked numerically at several points with a central difference
/// of the returned closed form.
fn check_antiderivative(f: &str) {
    let (primitive, reduced) = reduce_with(&[calculus()], &format!("integral({f}, x)"), &[]);
    assert!(reduced, "∫ {f} dx was not found: {primitive}");
    for at in [0.4, 0.9, 1.7] {
        let h = 1e-5;
        let up = eval(&[calculus()], &primitive, &[("x", at + h), ("a", 1.3), ("b", 0.6)]);
        let down = eval(&[calculus()], &primitive, &[("x", at - h), ("a", 1.3), ("b", 0.6)]);
        let derivative = (up - down) / (2.0 * h);
        let want = eval(&[calculus()], f, &[("x", at), ("a", 1.3), ("b", 0.6)]);
        assert!(
            (derivative - want).abs() < 1e-5 * (1.0 + want.abs()),
            "∫ {f} dx = {primitive}, but its derivative at {at} is {derivative}, not {want}"
        );
    }
}

/// A definite integral in both phases: the closed form, evaluated, must
/// agree with the numeric quadrature of the same request.
fn check_definite(
    f: &str,
    a: &str,
    b: &str,
    expected: f64,
) {
    let request = format!("defint({f}, x, {a}, {b})");
    let (closed, reduced) = reduce_with(&[calculus()], &request, &[]);
    assert!(reduced, "{request} has no closed form: {closed}");
    let from_closed = eval(&[calculus()], &closed, &[]);
    let (value, error) = numeric(&[calculus()], &request, &[], 1e-10);
    assert!((from_closed - expected).abs() < 1e-9 * (1.0 + expected.abs()), "{request} = {closed} = {from_closed}, expected {expected}");
    assert!((value - expected).abs() < 1e-8 * (1.0 + expected.abs()), "{request} numerically {value} ± {error}, expected {expected}");
}

#[test]
fn table_integrals() {
    assert_eq!(run("integral(x^2, x)"), "1/3*x^3");
    assert_eq!(run("integral(1/x, x)"), "ln(x)");
    assert_eq!(run("integral(cos(x), x)"), "sin(x)");
    assert_eq!(run("integral(exp(2*x), x)"), "1/2*exp(2*x)");
    assert_eq!(run("integral(a, x)"), "a*x");
    assert_eq!(run("integral(3*x^2 + 2*x + 1, x)"), "x^3 + x^2 + x");
    assert_eq!(run("integral(1/(2*x + 1), x)"), "1/2*ln(2*x + 1)");
}

#[test]
fn antiderivatives_differentiate_back() {
    for f in [
        "x^5 - 3*x + 2",
        "(3*x + 1)^4",
        "sqrt(x)",
        "1/x^3",
        "2^x",
        "sin(3*x + 1)",
        "tan(x)",
        "ln(x)",
        "atan(x)",
        "asin(x/2)",
        "sinh(2*x) + cosh(x/3) + tanh(x)",
        "sin(x)^2",
        "cos(2*x)^2",
        "1/cos(x)^2",
        "exp(a*x)*sin(b*x)",
        "exp(-x)*cos(2*x)",
    ] {
        check_antiderivative(f);
    }
}

#[test]
fn rational_functions() {
    for f in [
        "1/(x^2 + 1)",
        "1/(x^2 - 1)",
        "x/(x^2 + 1)",
        "(x^3 + 1)/(x^2 + 3*x + 2)",
        "1/(x*(x + 1)^2)",
        "(2*x + 3)/(x^2 + 2*x + 5)",
        "1/(x^3 + x)",
        "x^4/(x^2 + 4)",
    ] {
        check_antiderivative(f);
    }
    assert_eq!(run("integral(1/(x^2 + 1), x)"), "atan(x)");
}

#[test]
fn substitution() {
    for f in [
        "x*exp(x^2)",
        "cos(x)*sin(x)^3",
        "ln(x)/x",
        "x/sqrt(x^2 + 1)",
        "exp(x)/(1 + exp(x))",
        "sin(x)*exp(cos(x))",
        "x^2*cos(x^3)",
        "1/(x*ln(x))",
    ] {
        check_antiderivative(f);
    }
}

#[test]
fn by_parts() {
    for f in ["x*exp(x)", "x^2*exp(x)", "x*sin(x)", "x^2*cos(2*x)", "x*ln(x)", "x^3*ln(x)", "x*atan(x)"] {
        check_antiderivative(f);
    }
}

#[test]
fn powers_of_sines_and_cosines() {
    for f in ["sin(x)^3", "cos(x)^5", "sin(x)^2*cos(x)^3", "sin(x)^4", "sin(x)^2*cos(x)^2", "sin(2*x)^3*cos(2*x)^2"] {
        check_antiderivative(f);
    }
}

#[test]
fn what_cannot_be_integrated_stays_a_request() {
    for src in ["integral(exp(x^2), x)", "integral(sin(x)/x, x)", "integral(f(x), x)", "integral(x^x, x)"] {
        let (text, reduced) = reduce_with(&[calculus()], src, &[]);
        assert!(!reduced, "{src} unexpectedly gave {text}");
    }
}

#[test]
fn fundamental_theorem() {
    assert_eq!(run("diff(integral(f(x), x), x)"), "apply(f, x)");
    assert_eq!(run("diff(integral(x*exp(x), x), x)"), "x*exp(x)");
    assert_eq!(run("defint(f(x), x, a, a)"), "0");
}

#[test]
fn definite_integrals_agree_across_phases() {
    check_definite("x^2", "0", "1", 1.0 / 3.0);
    check_definite("sin(x)", "0", "pi", 2.0);
    check_definite("exp(x)", "0", "1", std::f64::consts::E - 1.0);
    check_definite("1/x", "1", "2", std::f64::consts::LN_2);
    check_definite("1/(1 + x^2)", "0", "1", std::f64::consts::FRAC_PI_4);
    check_definite("x*exp(-x)", "0", "2", 1.0 - 3.0 * (-2.0_f64).exp());
    check_definite("cos(x)^2", "0", "pi", std::f64::consts::FRAC_PI_2);
    assert_eq!(run("defint(x^2, x, a, b)"), "1/3*b^3 - 1/3*a^3");
    assert_eq!(run("defint(2*x, x, 0, 3)"), "9");
}

#[test]
fn numeric_quadrature_without_a_closed_form() {
    // Gaussian integral over the whole line.
    let (value, error) = numeric(&[calculus()], "defint(exp(-x^2), x, -oo, oo)", &[], 1e-10);
    assert!((value - std::f64::consts::PI.sqrt()).abs() < 1e-9, "{value} ± {error}");
    assert!(error < 1e-8);
    // Sine integral Si(1).
    let (value, _) = numeric(&[calculus()], "defint(sin(x)/x, x, 0, 1)", &[], 1e-12);
    assert!((value - 0.946_083_070_367_183).abs() < 1e-10, "{value}");
    // Parameters come from the bindings.
    let (value, _) = numeric(&[calculus()], "defint(exp(-a*x^2), x, 0, oo)", &[("a", 2.0)], 1e-10);
    assert!((value - 0.5 * (std::f64::consts::PI / 2.0).sqrt()).abs() < 1e-9, "{value}");
    // Without a binding for the parameter there is no number to give.
    let (value, _) = numeric(&[calculus()], "defint(exp(-a*x^2), x, 0, 1)", &[], 1e-10);
    assert!(value.is_nan());
}

#[test]
fn nested_requests() {
    assert_eq!(run("diff(defint(t^2, t, 0, x), x)"), "x^2");
    assert_eq!(run("integral(diff(x^3, x), x)"), "x^3");
    let (value, _) = numeric(&[calculus()], "defint(defint(x*y, y, 0, x), x, 0, 1)", &[], 1e-10);
    assert!((value - 0.125).abs() < 1e-9, "{value}");
}

