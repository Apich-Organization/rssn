//! Solving: a polynomial, a system, and an inequality. The answer is a term
//! like any other, so it can be read back as data or evaluated.
//!
//! Run with `cargo run --release --example solve_equations`.

use rssn::api::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();

    // --- one polynomial -------------------------------------------------
    let p = s.parse("x^2 - 5*x + 6")?;
    let roots = s.compute(s.parse("solve(x^2 - 5*x + 6, x)")?, &Config::new())?;
    println!("solve x^2 - 5x + 6 = 0  ->  {}", roots.term);

    // A `list(...)` answer is a value: read it as a vector and check it.
    let (shape, data) = roots.term.as_tensor().ok_or("the answer is not a list")?;
    for &r in &data {
        let value = p.eval(&[("x", r)]).unwrap_or(f64::NAN);
        println!("  p({r}) = {value:.1e}");
    }
    println!("  shape {shape:?}, {} roots", data.len());

    // An equation with no closed-form roots: the numeric solver reports one
    // root near the guess, and the numeric phase carries an error bound.
    let fixed = s.compute(
        s.parse("nsolve(cos(x) = x, x, 1)")?,
        &Config::new().numeric(1e-12),
    )?;
    let r = fixed.value.ok_or("nsolve returned no value")?;
    println!(
        "\nnsolve cos(x) = x       ->  {r:.12}  (cos r - r = {:.1e})",
        r.cos() - r
    );

    // --- a system -------------------------------------------------------
    let system = s.compute(
        s.parse("solve(list(x + y = 3, x - y = 1), list(x, y))")?,
        &Config::new(),
    )?;
    println!("\nsolve {{x + y = 3, x - y = 1}}  ->  {}", system.term);
    if let Some((_, tuple)) = system.term.as_tensor() {
        let (x, y) = (tuple[0], tuple[1]);
        println!(
            "  x + y = {:.1} (want 3), x - y = {:.1} (want 1)",
            x + y,
            x - y
        );
    }

    // An underdetermined system answers with a free parameter rather than
    // pretending the solution is unique.
    let line = s.compute(
        s.parse("solve(list(x + y = 1, 2*x + 2*y = 2), list(x, y))")?,
        &Config::new(),
    )?;
    println!("underdetermined        ->  {}", line.term);

    // --- an inequality --------------------------------------------------
    let set = s.compute(s.parse("solve(lt(x^2 - 4, 0), x)")?, &Config::new())?;
    println!("\nsolve x^2 - 4 < 0      ->  {}", set.term);

    // The solution set is the open interval (-2, 2): a point inside
    // satisfies the inequality, a point outside does not.
    for probe in [0.0, 3.0] {
        let inside = probe * probe - 4.0 < 0.0;
        println!(
            "  x = {probe}: x^2 - 4 = {:.0}, satisfies: {inside}",
            probe * probe - 4.0
        );
    }
    Ok(())
}
