//! A tour of the compute API: one entry point, two target phases.
//!
//! Run with `cargo run --release --example quickstart`.

use rssn::api::*;
use rssn::graph::Facts;

fn main() -> Result<(), ComputeError> {
    let s = Session::new();
    let x = s.sym("x");

    // Symbolic: a closed form.
    let d = s.compute(diff(sin(x) * cos(x), x), &Config::new())?;
    println!("d/dx sin(x)cos(x) = {}    LaTeX: {}", d.term, d.term.to_latex());

    // Numeric: the same request evaluated at a point, with an error bound.
    let n = s.compute(diff(sin(x) * cos(x), x), &Config::new().numeric(1e-12).bind("x", 0.3))?;
    println!("  at x = 0.3: {:?} ± {:?}", n.value, n.error);

    // Everything else is a request term too, parsed or built.
    for request in [
        "integral(x^2*exp(x), x)",
        "defint(exp(-x^2), x, -oo, oo)",
        "sum(1/(k*(k+1)), k, 1, n)",
        "limit(sin(x)/x, x, 0)",
        "solve(x^2 - 5*x + 6, x)",
        "dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))",
        "factor(x^4 - 1)",
        "det(list(list(a, b), list(c, d)))",
    ] {
        let answer = s.compute(s.parse(request)?, &Config::new())?;
        println!("{request:45} = {}", answer.term);
    }

    // Assumptions enable conditional identities.
    let a = s.parse("sqrt(a^2)")?;
    let plain = s.compute(a, &Config::new())?;
    let positive = s.compute(a, &Config::new().assume("a", Facts::POSITIVE))?;
    println!("sqrt(a^2) = {}   (a > 0: {})", plain.term, positive.term);

    // A pretty drawing, and a compiled function.
    let f = s.parse("x^2/(1 + x)")?;
    println!("{}", f.to_pretty());
    if let Ok(compiled) = f.compile(&["x"]) {
        println!("f(2) = {}", compiled.call(&[2.0]));
    }
    Ok(())
}
