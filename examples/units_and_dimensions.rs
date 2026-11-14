//! Units and dimensions: exact unit arithmetic, inference over an ordinary
//! expression, and the strict mode that turns a mismatch into an error.
//!
//! Run with `cargo run --release --example units_and_dimensions`.

use rssn::api::*;
use rssn::graph::DimensionMode;

fn ask(
    s: &Session,
    request: &str,
) -> Result<String, Box<dyn std::error::Error>> {
    let term = s.parse(request)?;
    Ok(s.compute(term, &Config::new())?.term.to_string())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();

    // A quantity keeps its unit: nothing is converted behind your back.
    println!(
        "36 km/h                 = {}",
        ask(&s, "quantity(36, km/h)")?
    );
    // Conversion is explicit and exact: 36 km/h is 10 m/s.
    println!(
        "  in m/s                = {}",
        ask(&s, "convert(quantity(36, km/h), m/s)")?
    );
    // A derived unit keeps the units it was built from.
    println!(
        "  2 kg * 3 m/s^2        = {}",
        ask(&s, "unify_expression(quantity(2, kg) * quantity(3, m/s^2))")?
    );
    // The dimension of a unit is a base-dimension product, not a unit name.
    println!(
        "  dimension of N*m      = {}",
        ask(&s, "dimension(quantity(1, N*m))")?
    );

    // Dimensions are also inferred over an ordinary expression, from the
    // declarations in the config. This is advisory: it never changes the
    // answer, it only reports what the answer's dimension is.
    let x = s.sym("x");
    let t = s.sym("t");
    let speed = x / t;
    let config = Config::new().dimension("x", "m").dimension("t", "s");
    let answer = s.compute(speed, &config)?;
    let dimension = answer
        .dimension
        .map(|d| d.to_string())
        .unwrap_or_else(|| "none".to_owned());
    println!("\ndimension of x/t        = {dimension}  (x in m, t in s)");

    // Strict mode makes an ill-typed expression a hard error instead of a
    // diagnostic, which is what a generated front end wants.
    let strict = Config::new()
        .dimensions(DimensionMode::Strict)
        .dimension("x", "m")
        .dimension("y", "s");
    match s.compute(s.parse("x + y")?, &strict) {
        | Err(e) => println!("strict: x + y rejected: {e}"),
        | Ok(a) => println!("strict: x + y accepted as {}", a.term),
    }

    // The same request in the default advisory mode answers and reports.
    let advisory = s.compute(
        s.parse("x + y")?,
        &Config::new().dimension("x", "m").dimension("y", "s"),
    )?;
    println!(
        "advisory: x + y = {} with {} diagnostic(s)",
        advisory.term,
        advisory.diagnostics.len()
    );
    Ok(())
}
