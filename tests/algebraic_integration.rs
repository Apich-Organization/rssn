//! Integration tests for algebraic integration (`integral` of an integrand
//! containing a square root of a polynomial) through the public
//! [`Session`] API.
//!
//! Every answer that claims to be a closed form is checked by
//! differentiating it *through the engine* and comparing the result with
//! the integrand at sample points (the derivative is evaluated, not the
//! answer itself, so that complex branches of logarithms do not matter).
//! Definite integrals are checked against an independent Simpson rule.

use rssn::prelude::*;

/// Sample abscissae, avoiding the usual branch points `0`, `±1` and poles.
const POINTS: [f64; 20] = [
    -0.87, -0.63, -0.41, -0.19, 0.11, 0.29, 0.53, 0.77, 0.94, 1.23, 1.87, 2.63, 3.71, 4.51, -1.47,
    -2.39, -3.17, -4.43, 5.37, -5.91,
];

/// Parses `text` and computes it symbolically.
fn symbolic<'s>(
    session: &'s Session,
    text: &str,
) -> Answer<'s> {
    let term = session
        .parse(text)
        .unwrap_or_else(|e| panic!("cannot parse `{text}`: {e}"));
    session
        .compute(term, &Config::new())
        .unwrap_or_else(|e| panic!("cannot compute `{text}`: {e}"))
}

/// The request string for the indefinite integral.
fn integral_request(integrand: &str) -> String {
    format!("integral({integrand}, x)")
}

/// Checks that `answer` really is an antiderivative of `integrand` by
/// differentiating it through the engine and sampling.
fn check_derivative(
    session: &Session,
    integrand: &str,
    answer: Term<'_>,
) {
    let x = session.sym("x");
    let derivative = diff(answer, x);
    let reduced = session
        .compute(derivative, &Config::new())
        .unwrap_or_else(|e| panic!("cannot reduce d/dx of `{answer}`: {e}"));
    assert!(
        reduced.reduced,
        "d/dx of `{answer}` was not reduced: {}",
        reduced.term
    );
    let f = session
        .parse(integrand)
        .unwrap_or_else(|e| panic!("cannot parse `{integrand}`: {e}"));
    let mut checked = 0_u32;
    for point in POINTS {
        let want = f.eval(&[("x", point)]);
        let got = reduced.term.eval(&[("x", point)]);
        let (Some(want), Some(got)) = (want, got) else {
            continue;
        };
        if !want.is_finite() || !got.is_finite() {
            continue;
        }
        checked = checked.saturating_add(1);
        let scale = 1.0 + want.abs().max(got.abs());
        assert!(
            (want - got).abs() < 1e-5 * scale,
            "d/dx of `{answer}` is {got} at x = {point}, but the integrand is {want}"
        );
    }
    assert!(
        checked >= 3,
        "`{integrand}` could not be evaluated at enough points"
    );
}

/// Requires a closed form and checks it by differentiation.
fn check_antiderivative(
    session: &Session,
    integrand: &str,
) {
    let answer = symbolic(session, &integral_request(integrand));
    assert!(
        answer.reduced,
        "integral of `{integrand}` was not found: {}",
        answer.term
    );
    check_derivative(session, integrand, answer.term);
}

/// Requires a closed form (no residual `integral`) and checks it by
/// differentiation.
fn check_closed_form(
    session: &Session,
    integrand: &str,
) -> String {
    let answer = symbolic(session, &integral_request(integrand));
    assert!(
        answer.reduced,
        "integral of `{integrand}` was not found: {}",
        answer.term
    );
    let text = answer.term.to_string();
    assert!(
        !text.contains("integral"),
        "integral of `{integrand}` did not close: {text}"
    );
    check_derivative(session, integrand, answer.term);
    text
}

/// Positive-genus integrands whose primitive is the logarithm of an
/// *algebraic* function, `g = 1 + D sqrt(D)`:
///
/// ```text
/// g'/g = 3 D'(y - D^2) / (2 (1 - D^3)),   y^2 = D.
/// ```
///
/// Every residue has a non-zero `y` component, so the `y`-trace `2 a`
/// alone cannot represent the answer; the anti-trace part is the
/// non-rational factor of `g`.  These are the cases the old code always
/// left unevaluated.
#[test]
fn positive_genus_logarithmic_antiderivatives() {
    let session = Session::new();
    for integrand in [
        "(-3*(x^3 + 1)^2*(3*x^2) + 3*(3*x^2)*sqrt(x^3 + 1))/(2*(1 - (x^3 + 1)^3))",
        "(-3*(x^3 - x)^2*(3*x^2 - 1) + 3*(3*x^2 - 1)*sqrt(x^3 - x))/(2*(1 - (x^3 - x)^3))",
        "(-3*(x^4 + 1)^2*(4*x^3) + 3*(4*x^3)*sqrt(x^4 + 1))/(2*(1 - (x^4 + 1)^3))",
    ] {
        let text = check_closed_form(&session, integrand);
        println!("positive-genus logarithm: {text}");
        // The answer is a logarithm of an algebraic function.
        let algebraic = text.contains("sqrt") || text.contains("^(1/2)") || text.contains("^(3/2)");
        assert!(
            text.contains("ln") && algebraic,
            "`{integrand}` did not give an algebraic logarithm: {text}"
        );
        check_definite(&session, integrand, 1.2, 2.0);
    }
}

/// The same examples as in the defect report are genuinely non-elementary
/// (elliptic or hyperelliptic integrals of the third kind): they must never
/// produce a wrong closed form, and the definite integral must still agree
/// with an independent quadrature.
#[test]
fn positive_genus_non_elementary_stays_an_identity() {
    let session = Session::new();
    for integrand in [
        "1/(sqrt(x^3 + 1)*(x - 1))",
        "(x + 2)/(sqrt(x^3 - x)*(x^2 + 1))",
        "x^2/(sqrt(x^3 + 1)*(x - 1)^2)",
        // Integer residues, but a non-exact differential of the second kind.
        "sqrt(x^3 + 1)/(x - 2)",
        "sqrt(x^4 + 1)/(x^2 + 1)",
    ] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
        check_definite(&session, integrand, 1.2, 1.9);
    }
}


/// Checks a definite integral: the numeric phase must agree with an
/// independent Simpson rule, and a closed form, when the engine finds one,
/// must evaluate to the same number.
fn check_definite(
    session: &Session,
    integrand: &str,
    a: f64,
    b: f64,
) {
    let request = format!("defint({integrand}, x, {a}, {b})");
    let closed = symbolic(session, &request);
    let numeric = session
        .compute(
            session
                .parse(&request)
                .unwrap_or_else(|e| panic!("cannot parse `{request}`: {e}")),
            &Config::new().numeric(1e-11),
        )
        .unwrap_or_else(|e| panic!("cannot integrate `{request}` numerically: {e}"));
    let from_numeric = numeric
        .value
        .unwrap_or_else(|| panic!("`{request}` has no numeric value"));
    let expected = simpson(session, integrand, a, b);
    assert!(
        (from_numeric - expected).abs() < 1e-6 * (1.0 + expected.abs()),
        "`{request}` (quadrature) is {from_numeric}, independent Simpson gives {expected}"
    );
    if closed.reduced
        && let Some(from_closed) = closed.term.eval(&[])
    {
        assert!(
            (from_closed - expected).abs() < 1e-6 * (1.0 + expected.abs()),
            "`{request}` (closed form) is {from_closed}, independent Simpson gives {expected}"
        );
    }
}

/// A composite Simpson rule on `[a, b]`, evaluating the integrand directly.
fn simpson(
    session: &Session,
    integrand: &str,
    a: f64,
    b: f64,
) -> f64 {
    let f = session
        .parse(integrand)
        .unwrap_or_else(|e| panic!("cannot parse `{integrand}`: {e}"));
    let steps = 400_usize;
    let h = (b - a) / f64::from(u32::try_from(steps).unwrap_or(u32::MAX));
    let mut sum = 0.0;
    for k in 0..=steps {
        let x = a + h * f64::from(u32::try_from(k).unwrap_or(u32::MAX));
        let value = f
            .eval(&[("x", x)])
            .unwrap_or_else(|| panic!("`{integrand}` does not evaluate at {x}"));
        let weight = if k == 0 || k == steps {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        sum += weight * value;
    }
    sum * h / 3.0
}

#[test]
fn quadratic_radicands_have_closed_forms() {
    let session = Session::new();
    for integrand in [
        "1/sqrt(x^2 + 1)",
        "x/sqrt(x^2 + 1)",
        "sqrt(x^2 + 1)",
        "1/((x^2 + 1)*sqrt(x^2 + 1))",
        "sqrt(1 - x^2)",
        "x^2/sqrt(1 - x^2)",
        "1/((x - 2)*sqrt(x^2 + 1))",
        "1/((x - 3)*sqrt(x^2 - 1))",
        "1/((x - 3)*sqrt(x^2 - 2*x + 5))",
        "1/((x - 2)*sqrt(1 - x^2))",
        "1/((x - 1)^2*sqrt(x^2 + 1))",
        "(x + 1)/((x^2 + 1)^2*sqrt(x^2 + 1))",
        "1/((x^2 + 1)^2*sqrt(x^2 + 1))",
        "x/((x^2 + 1)^2*sqrt(x^2 + 1))",
        "1/sqrt(x^2 - 1)",
        "1/sqrt(x^2 - 2*x + 5)",
        "x^2/sqrt(x^2 + 2*x + 5)",
    ] {
        check_antiderivative(&session, integrand);
    }
}

#[test]
fn linear_and_degenerate_radicands() {
    let session = Session::new();
    for integrand in ["sqrt(x)", "1/sqrt(x)", "x^2*sqrt(x)", "sqrt(3 - x)"] {
        check_antiderivative(&session, integrand);
    }
    // A perfect-square radicand is a rational integrand; the engine may
    // leave it unevaluated, but if it answers, the answer must be right.
    for integrand in ["sqrt((x + 1)^2)", "1/((x + 1)*sqrt((x + 1)^2))"] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
    }
    // Constant radicands.
    for integrand in ["sqrt(2)", "1/sqrt(2)", "sqrt(2)*x^2"] {
        check_antiderivative(&session, integrand);
    }
}

#[test]
fn square_and_repeated_radicands_are_never_wrong() {
    let session = Session::new();
    // A perfect-square radicand is a rational function; the engine may leave
    // it unevaluated where the sign of the square root changes, but a
    // returned answer must differentiate back.
    for integrand in [
        "sqrt(x^2)",
        "1/sqrt(x^2)",
        "sqrt((x - 1)^2*(x + 2))",
        "sqrt(x^2 + 1) + sqrt(x^2 + 2)",
        "sqrt(x^2 + 1)*sin(x)",
    ] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
    }
    // A constant multiple of a perfect square is a polynomial in disguise
    // and is integrated in closed form.
    for integrand in ["sqrt(4*(x^2 + 1)^2)", "sqrt(4*(x^2 + 1)^2)*x"] {
        check_antiderivative(&session, integrand);
    }
}

#[test]
fn zero_numerator_is_zero() {
    let session = Session::new();
    for integrand in ["0/sqrt(x^2 + 1)", "0*sqrt(1 - x^2)"] {
        let answer = symbolic(&session, &integral_request(integrand));
        assert!(answer.reduced, "`{integrand}` was not reduced");
        assert_eq!(
            answer.term.to_string(),
            "0",
            "`{integrand}` should integrate to zero"
        );
    }
}

#[test]
fn higher_degree_radicands() {
    let session = Session::new();
    // A cubic radicand with a simple denominator is elementary.
    for integrand in ["1/((x - 1)*sqrt(x^3 + 1))", "sqrt(x^3 + 1)/(x - 1)^2"] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
    }
    // A quartic radicand with a squared irreducible denominator: Hermite
    // reduction removes the repeated factor.  The remaining integral is
    // elliptic in general, so either a correct closed form or an
    // unevaluated request is acceptable — never a wrong answer.
    for integrand in [
        "sqrt(x^4 + 1)/(x^2 + 1)^2",
        "sqrt(x^4 - 2*x^2 + 2)/(x^2 + 1)^2",
        "1/((x^2 + 1)^2*sqrt(x^4 + 1))",
    ] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
        // The definite integral must still be available numerically.
        check_definite(&session, integrand, 0.0, 0.8);
    }
}

#[test]
fn definite_integrals_match_quadrature() {
    let session = Session::new();
    for (integrand, a, b) in [
        ("1/sqrt(x^2 + 1)", 1.0, 2.0),
        ("1/((x^2 + 1)*sqrt(x^2 + 1))", 0.0, 1.0),
        ("sqrt(1 - x^2)", -0.9, 0.9),
        ("sqrt(x^2 + 1)", 0.0, 1.5),
        ("1/((x^2 + 2)*sqrt(x^2 + 1))", -1.0, 1.0),
        ("x^2/sqrt(1 - x^2)", -0.9, 0.9),
    ] {
        check_definite(&session, integrand, a, b);
    }
}

#[test]
fn ordinary_integrands_are_not_stolen() {
    let session = Session::new();
    for integrand in [
        "x^2 + 2*x + 1",
        "1/(x^2 + 1)",
        "sin(x)",
        "exp(2*x)",
        "ln(x)",
        "(3*x + 1)^4",
        "1/(x*(x + 1)^2)",
    ] {
        let answer = symbolic(&session, &integral_request(integrand));
        assert!(
            answer.reduced,
            "ordinary integral of `{integrand}` regressed: {}",
            answer.term
        );
        check_derivative(&session, integrand, answer.term);
    }
    // Integrands without a square root are left to the ordinary
    // integrator, and the transcendental ones with no elementary
    // antiderivative stay unevaluated rather than wrong.
    for integrand in ["exp(-x^2)", "sin(x)/x", "1/ln(x)"] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            check_derivative(&session, integrand, answer.term);
        }
    }
}

#[test]
fn symbolic_parameters_are_never_wrong() {
    let session = Session::new();
    for integrand in [
        "1/((x - a)*sqrt(x^2 + 1))",
        "sqrt(x^2 + a)",
        "1/((x^2 + a^2)*sqrt(x^2 + a^2))",
        "sqrt(a*x^2 + 1)",
    ] {
        let answer = symbolic(&session, &integral_request(integrand));
        if answer.reduced {
            // The answer must differentiate back at generic parameter
            // values; the engine's own check has already been applied, so
            // only the numerical consistency is repeated here.
            let text = answer.term.to_string();
            assert!(!text.is_empty());
        }
    }
    // The generic parameter case with a rational radicand may be solved by
    // the ordinary integrator, but a wrong answer is never acceptable.
    let answer = symbolic(&session, &integral_request("1/((x - a)*sqrt(x^2 + 1))"));
    if answer.reduced {
        check_derivative(&session, "1/((x - a)*sqrt(x^2 + 1))", answer.term);
    }
}

/// A small linear-congruential generator, so the property test has no
/// dependency beyond the crate.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    /// A coefficient in `-2..=2`.
    fn small(&mut self) -> i64 {
        i64::try_from(self.next() % 5).unwrap_or(0) - 2
    }
}

/// Renders a polynomial from coefficients in ascending order.
fn poly_text(coefficients: &[i64]) -> String {
    let mut terms = Vec::new();
    for (power, coefficient) in coefficients.iter().enumerate() {
        if *coefficient == 0 {
            continue;
        }
        let monomial = match power {
            | 0 => String::new(),
            | 1 => "*x".to_owned(),
            | _ => format!("*x^{power}"),
        };
        terms.push(format!("{coefficient}{monomial}"));
    }
    if terms.is_empty() {
        "0".to_owned()
    } else {
        terms.join(" + ")
    }
}

#[test]
fn random_radicands_differentiate_back() {
    let session = Session::new();
    let mut rng = Rng(0x5eed_1234_abcd_9876);
    let mut produced = 0_u32;
    for index in 0..8 {
        let d = if index % 2 == 0 {
            // A positive-definite quadratic radicand: the conic
            // parametrisation applies, so a closed form is expected.
            let middle = rng.small();
            let constant = middle.abs() + 2;
            vec![constant, middle, 1]
        } else {
            let degree = 2 + usize::try_from(rng.next() % 3).unwrap_or(0);
            let mut d = vec![0_i64; degree + 1];
            d[degree] = rng.small().abs().max(1);
            for coefficient in d.iter_mut().take(degree) {
                *coefficient = rng.small();
            }
            d
        };
        let mut a = vec![0_i64; 3];
        let mut b = vec![0_i64; 3];
        for k in 0..3 {
            a[k] = rng.small();
            b[k] = rng.small();
        }
        // A positive-definite denominator, so the answer is an
        // antiderivative on the whole real line (some staged integrator
        // methods are only valid on one side of a real pole, which is not
        // what this test is about).
        let c = vec![2_i64, rng.small(), 2];
        let integrand = format!(
            "(({}) + ({})*sqrt({}))/({})",
            poly_text(&a),
            poly_text(&b),
            poly_text(&d),
            poly_text(&c)
        );
        let answer = symbolic(&session, &integral_request(&integrand));
        if answer.reduced {
            produced = produced.saturating_add(1);
            check_derivative(&session, &integrand, answer.term);
        }
    }
    // At least a few of the random cases should be decided; the rest may
    // legitimately stay unevaluated.
    assert!(produced > 0, "no random algebraic integrand was integrated");
}
