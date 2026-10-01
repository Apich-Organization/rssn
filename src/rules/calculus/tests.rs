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
    let mut checked = 0;
    for at in [0.4, 0.9, 1.7, 2.6] {
        let h = 1e-5;
        let up = eval(&[calculus()], &primitive, &[("x", at + h), ("a", 1.3), ("b", 0.6)]);
        let down = eval(&[calculus()], &primitive, &[("x", at - h), ("a", 1.3), ("b", 0.6)]);
        let derivative = (up - down) / (2.0 * h);
        let want = eval(&[calculus()], f, &[("x", at), ("a", 1.3), ("b", 0.6)]);
        // A real logarithm of a negative number is NaN here although the
        // antiderivative is fine as a complex function: skip such points.
        if !derivative.is_finite() || !want.is_finite() {
            continue;
        }
        checked += 1;
        assert!(
            (derivative - want).abs() < 1e-5 * (1.0 + want.abs()),
            "∫ {f} dx = {primitive}, but its derivative at {at} is {derivative}, not {want}"
        );
    }
    assert!(checked > 0, "∫ {f} dx = {primitive} could not be evaluated anywhere");
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
    for f in ["sin(2*x)*sin(x)", "cos(3*x)*cos(x)", "x*sin(x)*cos(2*x)", "sin(x)^3", "cos(x)^5", "sin(x)^2*cos(x)^3", "sin(x)^4", "sin(x)^2*cos(x)^2", "sin(2*x)^3*cos(2*x)^2"] {
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
    assert_eq!(run("diff(integral(f(x), x), x)"), "f(x)");
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



/// A limit in both phases.
fn check_limit(
    request: &str,
    expected: f64,
) {
    let (closed, reduced) = reduce_with(&[calculus()], request, &[]);
    assert!(reduced, "{request} was not resolved: {closed}");
    let value = eval(&[calculus()], &closed, &[]);
    assert!(
        (value - expected).abs() < 1e-9 * (1.0 + expected.abs()) || value.total_cmp(&expected).is_eq(),
        "{request} = {closed} = {value}, expected {expected}"
    );
    if expected.is_finite() {
        let (numeric_value, _) = numeric(&[calculus()], request, &[], 1e-6);
        assert!((numeric_value - expected).abs() < 1e-5 * (1.0 + expected.abs()), "{request} numerically {numeric_value}");
    }
}

#[test]
fn limits_by_continuity_and_cancellation() {
    assert_eq!(run("limit(x^2 + 1, x, 2)"), "5");
    assert_eq!(run("limit((x^2 - 1)/(x - 1), x, 1)"), "2");
    assert_eq!(run("limit(sin(x)/x, x, 0)"), "1");
    assert_eq!(run("limit(a*x + b, x, c)"), "a*c + b");
    check_limit("limit((1 - cos(x))/x^2, x, 0)", 0.5);
    check_limit("limit((exp(x) - 1 - x)/x^2, x, 0)", 0.5);
    check_limit("limit(tan(x)/x, x, 0)", 1.0);
    check_limit("limit((x^3 - 8)/(x - 2), x, 2)", 12.0);
    check_limit("limit(ln(1 + x)/x, x, 0)", 1.0);
    check_limit("limit(x*ln(x), x, 0, plus)", 0.0);
}

#[test]
fn limits_at_infinity() {
    assert_eq!(run("limit(1/x, x, oo)"), "0");
    assert_eq!(run("limit((2*x^2 + 1)/(x^2 - 3), x, oo)"), "2");
    assert_eq!(run("limit((3*x + 1)/(x^2 + 1), x, -oo)"), "0");
    assert_eq!(run("limit(x^2, x, -oo)"), "oo");
    assert_eq!(run("limit(x^3/(x + 1), x, -oo)"), "oo");
    assert_eq!(run("limit(exp(-x), x, oo)"), "0");
    assert_eq!(run("limit(atan(x), x, oo)"), "1/2*pi");
    assert_eq!(run("limit(atan(x), x, -oo)"), "-1/2*pi");
    check_limit("limit(x*exp(-x), x, oo)", 0.0);
    check_limit("limit(ln(x)/x, x, oo)", 0.0);
    check_limit("limit((1 + 1/x)^x, x, oo)", std::f64::consts::E);
    check_limit("limit(x^(1/x), x, oo)", 1.0);
    check_limit("limit(tanh(x), x, oo)", 1.0);
}

#[test]
fn infinite_and_one_sided_limits() {
    assert_eq!(run("limit(1/x^2, x, 0)"), "oo");
    assert_eq!(run("limit(1/x, x, 0, plus)"), "oo");
    assert_eq!(run("limit(1/x, x, 0, minus)"), "-oo");
    assert_eq!(run("limit(ln(x), x, 0, plus)"), "-oo");
    assert_eq!(run("limit(abs(x)/x, x, 0, minus)"), "-1");
    // No two-sided limit: the request stays.
    for src in ["limit(1/x, x, 0)", "limit(abs(x)/x, x, 0)", "limit(sin(1/x), x, 0)", "limit(sin(x), x, oo)"] {
        let (text, reduced) = reduce_with(&[calculus()], src, &[]);
        assert!(!reduced, "{src} unexpectedly gave {text}");
    }
}

#[test]
fn taylor_and_laurent() {
    assert_eq!(run("taylor(exp(x), x, 0, 4)"), "1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1");
    assert_eq!(run("taylor(sin(x), x, 0, 5)"), "1/120*x^5 - 1/6*x^3 + x");
    assert_eq!(run("taylor(1/(1 - x), x, 0, 3)"), "x^3 + x^2 + x + 1");
    assert_eq!(run("taylor(ln(x), x, 1, 3)"), "1/3*(x - 1)^3 - 1/2*(x - 1)^2 + x - 1");
    assert_eq!(run("taylor(sin(x)/x, x, 0, 4)"), "1/120*x^4 - 1/6*x^2 + 1", "removable singularity");
    assert_eq!(run("laurent(1/(x*(1 - x)), x, 0, 2)"), "x^2 + x + 1 + 1/x");
    assert_eq!(run("laurent(exp(x)/x^2, x, 0, 1)"), "1/6*x + 1/2 + 1/x + 1/x^2");
    // The series approximates the function near the point.
    let series = run("taylor(cos(x)*exp(x), x, 0, 8)");
    let (got, want) = (eval(&[calculus()], &series, &[("x", 0.3)]), 0.3_f64.cos() * 0.3_f64.exp());
    assert!((got - want).abs() < 1e-8, "{series}: {got} vs {want}");
}

#[test]
fn sums_and_products() {
    assert_eq!(run("sum(k^2, k, 1, 4)"), "30");
    assert_eq!(run("sum(a^k, k, 0, 2)"), "a^2 + a + 1");
    assert_eq!(run("product(k, k, 1, 5)"), "120");
    assert_eq!(run("product(c, k, 1, n)"), "c^n");
    // Closed forms, checked against the written-out sum.
    for (summand, bindings) in [("k", 1.0), ("k^2", 1.0), ("k^3 - 2*k + 1", 1.0), ("3*r^k", 0.7), ("k*a + 2", 1.3)] {
        let (closed, reduced) = reduce_with(&[calculus()], &format!("sum({summand}, k, 2, n)"), &[]);
        assert!(reduced, "sum of {summand}: {closed}");
        let at = [("n", 9.0), ("r", bindings), ("a", bindings)];
        let got = eval(&[calculus()], &closed, &at);
        let want: f64 = (2..=9)
            .map(|k| eval(&[calculus()], summand, &[("k", f64::from(k)), ("r", bindings), ("a", bindings)]))
            .sum();
        assert!((got - want).abs() < 1e-9 * (1.0 + want.abs()), "sum of {summand} = {closed}: {got} vs {want}");
    }
    assert_eq!(run("sum((1/2)^k, k, 0, oo)"), "2");
    assert_eq!(run("sum(3*(1/3)^k, k, 1, oo)"), "3/2");
    let (text, reduced) = reduce_with(&[calculus()], "sum(2^k, k, 0, oo)", &[]);
    assert!(!reduced, "a divergent geometric series has no sum: {text}");
}

#[test]
fn numeric_sums() {
    let (value, error) = numeric(&[calculus()], "sum(1/k^2, k, 1, oo)", &[], 1e-10);
    assert!((value - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-5, "{value} ± {error}");
    let (value, _) = numeric(&[calculus()], "sum((-1)^(k + 1)/k, k, 1, oo)", &[], 1e-12);
    assert!((value - std::f64::consts::LN_2).abs() < 1e-9, "{value}");
    let (value, _) = numeric(&[calculus()], "sum(sin(k)/k^2, k, 1, 1000)", &[], 1e-12);
    let want: f64 = (1..=1000).map(|k| f64::from(k).sin() / f64::from(k * k)).sum();
    assert!((value - want).abs() < 1e-12, "{value}");
    let (value, _) = numeric(&[calculus()], "product(1 + 1/k^2, k, 1, 50)", &[], 1e-12);
    let want: f64 = (1..=50).map(|k| 1.0 + 1.0 / f64::from(k * k)).product();
    assert!((value - want).abs() < 1e-12, "{value}");
}

#[test]
fn convergence_tests() {
    assert_eq!(run("converges(1/2^k, k)"), "true");
    assert_eq!(run("converges(2^k/k^5, k)"), "false");
    assert_eq!(run("converges(1/k^2, k)"), "true");
    assert_eq!(run("converges(1/k, k)"), "false");
}

#[test]
fn fourier_series_of_simple_functions() {
    // x on [-pi, pi]: 2 sin x - sin 2x + (2/3) sin 3x - ...
    let series = run("fourier_series(x, x, pi, 3)");
    for at in [0.5, 1.5, 2.5] {
        let got = eval(&[calculus()], &series, &[("x", at)]);
        let want = 2.0 * (at.sin() - (2.0 * at).sin() / 2.0 + (3.0 * at).sin() / 3.0);
        assert!((got - want).abs() < 1e-9, "{series} at {at}: {got} vs {want}");
    }
    // An even function has no sine terms.
    let series = run("fourier_series(x^2, x, 1, 2)");
    let at = 0.4_f64;
    let pi = std::f64::consts::PI;
    let want = 1.0 / 3.0 - 4.0 / (pi * pi) * (pi * at).cos() + 1.0 / (pi * pi) * (2.0 * pi * at).cos();
    let got = eval(&[calculus()], &series, &[("x", at)]);
    assert!((got - want).abs() < 1e-9, "{series}: {got} vs {want}");
}


#[test]
fn improper_integrals() {
    assert_eq!(run("defint(exp(-x), x, 0, oo)"), "1");
    assert_eq!(run("defint(x*exp(-x), x, 0, oo)"), "1");
    assert_eq!(run("defint(1/(1 + x^2), x, -oo, oo)"), "pi");
    assert_eq!(run("defint(1/x^2, x, 1, oo)"), "1");
}

#[test]
fn derivatives_map_over_lists_and_equations() {
    assert_eq!(run("diff(list(t^2, 3, sin(t)), t)"), "list(2*t, 0, cos(t))");
    assert_eq!(run("diff(x^2 = y, x)"), "2*x = 0");
}

#[test]
fn polynomials_given_as_powers_of_sums() {
    assert_eq!(run("integral((x^2 - 1)^2, x)"), "1/5*x^5 - 2/3*x^3 + x");
    assert_eq!(run("defint((x^2 + x)^2, x, -1, 1)"), "16/15");
    // The variable of a definite integral knows its range.
    assert_eq!(run("defint(abs(x), x, 0, 2)"), "2");
    assert_eq!(run("defint((x^2)^(1/2), x, -3, 0)"), "9/2");
}

/// Gosper's algorithm closes hypergeometric sums that are not polynomial or
/// geometric: each closed form is checked against the written-out sum.
#[test]
fn gosper_sums() {
    let rules = crate::rules::standard();
    for (summand, from) in [("1/(k*(k+1))", 1), ("k*2^k", 0), ("factorial(k)*k", 0), ("k^3*3^k", 1), ("1/(4*k^2-1)", 1)] {
        let (closed, reduced) = reduce_with(&rules, &format!("sum({summand}, k, {from}, n)"), &[]);
        assert!(reduced, "sum of {summand}: {closed}");
        assert!(!closed.contains("sum("), "sum of {summand}: {closed}");
        for n in [5, 8] {
            let got = eval(&rules, &closed, &[("n", f64::from(n))]);
            let want: f64 = (from..=n).map(|k| eval(&rules, summand, &[("k", f64::from(k))])).sum();
            assert!((got - want).abs() < 1e-9 * (1.0 + want.abs()), "sum of {summand} = {closed}: {got} vs {want}");
        }
    }
    assert_eq!(simplify(&rules, "sum(1/(k*(k+1)), k, 1, n)"), "1 - 1/(n + 1)");
    assert_eq!(simplify(&rules, "sum(factorial(k)*k, k, 0, n)"), "factorial(n + 1) - 1");
}
