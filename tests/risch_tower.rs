//! The Risch–Bronstein tower integrator, exercised through the public API.
//!
//! Each antiderivative must come back reduced and differentiate back to the
//! integrand. The nested towers here are the ones the heuristic search does
//! not reach.

use rssn::api::Config;
use rssn::api::Session;
use rssn::api::diff;

/// `∫ f dx` is reduced and differentiates back to `f`.
fn check(
    session: &Session,
    src: &str,
) {
    let integrand = session
        .parse(src)
        .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    let request = session
        .parse(&format!("integral({src}, x)"))
        .unwrap_or_else(|e| panic!("cannot parse the request for `{src}`: {e}"));
    let answer = session
        .compute(request, &Config::new())
        .unwrap_or_else(|e| panic!("`{src}` failed: {e}"));
    assert!(answer.reduced, "`{src}` stayed a request: {}", answer.term);
    let back = session
        .compute(diff(answer.term, session.sym("x")), &Config::new())
        .unwrap_or_else(|e| panic!("cannot differentiate the answer for `{src}`: {e}"));
    let mut checked = 0;
    for point in [0.37, 0.91, 1.63, 2.3] {
        let (Some(want), Some(got)) = (
            integrand.eval(&[("x", point)]),
            back.term.eval(&[("x", point)]),
        ) else {
            continue;
        };
        if want.is_finite() && got.is_finite() {
            assert!(
                (want - got).abs() <= 1e-6 * (1.0 + want.abs().max(got.abs())),
                "`{src}` at {point}: the answer differentiates to {got}, not {want}"
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "`{src}`: no sample point was usable");
}

#[test]
fn rational_functions() {
    let session = Session::new();
    check(&session, "x^3");
    check(&session, "1/(x + 1)");
    check(&session, "(2*x + 1)/(x^2 + x + 1)");
    check(&session, "x/(x^2 + 1)^2");
    check(&session, "1/(x^2 - 1)");
}

#[test]
fn exponential_extensions() {
    let session = Session::new();
    check(&session, "exp(x)");
    check(&session, "1/exp(x)");
    check(&session, "x^2*exp(x)");
    check(&session, "1/(exp(x) + 1)");
    check(&session, "1/(exp(x) + 1)^2");
    check(&session, "exp(-2*x)");
}

#[test]
fn logarithmic_extensions() {
    let session = Session::new();
    check(&session, "ln(x)");
    check(&session, "ln(x)^3");
    check(&session, "ln(x)/x");
    check(&session, "1/(x*ln(x))");
}

#[test]
fn nested_towers() {
    let session = Session::new();
    // Q(x) < Q(x)(ln x) < Q(x)(ln x)(ln ln x)
    check(&session, "ln(ln(x))/(x*ln(x))");
    check(&session, "1/(x*ln(x)*ln(ln(x)))");
    // Q(x) < Q(x)(exp x) < Q(x)(exp x)(ln(exp x + 1))
    check(&session, "exp(x)*ln(exp(x) + 1)");
}
