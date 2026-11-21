//! Symbolic calculus end to end: differentiate, integrate, take a limit,
//! expand a series — then ask the same requests for numbers.
//!
//! Run with `cargo run --release --example symbolic_calculus`.

use rssn::api::*;

/// The numeric value of a numeric-target answer, or a readable failure.
fn number(answer: &Answer<'_>) -> Result<f64, String> {
    answer
        .value
        .ok_or_else(|| format!("no numeric value for `{}`", answer.term))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();
    let x = s.sym("x");

    // A derivative is a request term; the closed form is the answer.
    let f = sin(x) * cos(x);
    let d = s.compute(diff(f, x), &Config::new())?;
    // sin(x)cos(x) = sin(2x)/2, so the derivative is cos(2x).
    println!("d/dx sin(x)cos(x) = {}", d.term);

    // The *same* request with a numeric target: one request, one answer,
    // differing only in the phase that was asked for.
    let at = s.compute(diff(f, x), &Config::new().numeric(1e-12).bind("x", 0.3))?;
    println!(
        "  at x = 0.3: {:.12}  (cos 0.6 = {:.12})",
        number(&at)?,
        0.6_f64.cos()
    );

    // An antiderivative is checked by differentiating it again.
    let g = s.parse("x^2*exp(x)")?;
    let capital = s.compute(s.parse("integral(x^2*exp(x), x)")?, &Config::new())?;
    println!("\nintegral(x^2 e^x) = {}", capital.term);
    let back = s.compute(diff(capital.term, x), &Config::new())?;
    let (lhs, rhs) = (back.term.eval(&[("x", 0.7)]), g.eval(&[("x", 0.7)]));
    println!("  d/dx of it at 0.7: {lhs:?}, integrand: {rhs:?}");

    // A definite and an improper integral, symbolically and as a float.
    let gauss = s.parse("defint(exp(-x^2), x, -oo, oo)")?;
    let exact = s.compute(gauss, &Config::new())?;
    let approx = s.compute(gauss, &Config::new().numeric(1e-12))?;
    println!(
        "\ndefint(e^-x^2, -oo, oo) = {} = {:.12}",
        exact.term,
        number(&approx)?
    );

    // A limit that the rules resolve without a series.
    let l = s.compute(s.parse("limit(sin(x)/x, x, 0)")?, &Config::new())?;
    println!("limit sin(x)/x at 0 = {}", l.term);

    // A Taylor series: the answer is a pinned polynomial form, so its
    // numeric value is only an approximation of the function.
    let series = s.compute(s.parse("taylor(sin(x), x, 0, 5)")?, &Config::new())?;
    let s5 = series.term.eval(&[("x", 0.3)]).unwrap_or(f64::NAN);
    println!("\ntaylor(sin x, 0, 5) = {}", series.term);
    println!(
        "  at x = 0.3: {s5:.12}  (sin 0.3 = {:.12}, error {:.1e})",
        0.3_f64.sin(),
        (s5 - 0.3_f64.sin()).abs()
    );

    // LaTeX and Unicode renderings are one call away.
    println!("\nLaTeX: {}", d.term.to_latex());
    Ok(())
}
