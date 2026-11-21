//! Differential equations: one ODE with initial conditions and one PDE with
//! boundary conditions, each checked at a point against its closed form.
//!
//! Run with `cargo run --release --example ode_pde`.

use rssn::api::*;

/// The right-hand side of an `lhs = rhs` answer.
fn rhs<'s>(answer: &Answer<'s>) -> Option<Term<'s>> {
    answer.term.children().get(1).copied()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();

    // --- an ODE, general and with an initial condition -------------------
    let general = s.compute(
        s.parse("dsolve(diff(y(x), x) = y(x), y(x))")?,
        &Config::new(),
    )?;
    // The integration constant appears as C1; the equation is the answer.
    println!("y' = y                 ->  {}", general.term);

    let ivp = s.compute(
        s.parse("dsolve(diff(y(x), x) = y(x), y(x), list(y(0) = 3))")?,
        &Config::new(),
    )?;
    println!("y' = y, y(0) = 3       ->  {}", ivp.term);

    // The right-hand side is an ordinary term: evaluate it and compare.
    if let Some(solution) = rhs(&ivp) {
        let at = solution.eval(&[("x", 0.4)]).unwrap_or(f64::NAN);
        println!("  y(0.4) = {at:.9}  (3 e^0.4 = {:.9})", 3.0 * 0.4_f64.exp());
    }

    // Second order, two constants.
    let oscillator = s.compute(
        s.parse("dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))")?,
        &Config::new(),
    )?;
    println!("y'' + y = 0            ->  {}", oscillator.term);

    // --- a PDE with boundary and initial data ---------------------------
    // Heat equation on [0, pi] with fixed ends and a single sine mode:
    // the mode decays as exp(-4kt) without changing shape.
    let heat = s.compute(
        s.parse(
            "pdsolve(diff(u(x, t), t) = k*diff(diff(u(x, t), x), x), u(x, t), \
             list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = 3*sin(2*x)))",
        )?,
        &Config::new(),
    )?;
    println!("\nheat, u(x,0) = 3 sin 2x ->  {}", heat.term);

    if let Some(solution) = rhs(&heat) {
        // At t = 0 the answer must reproduce the initial data.
        let (u0, later) = (
            solution.eval(&[("k", 1.0), ("x", 0.5), ("t", 0.0)]),
            solution.eval(&[("k", 1.0), ("x", 0.5), ("t", 0.5)]),
        );
        println!(
            "  u(0.5, 0)   = {:.9}  (3 sin 1 = {:.9})",
            u0.unwrap_or(f64::NAN),
            3.0 * 1.0_f64.sin()
        );
        println!(
            "  u(0.5, 0.5) = {:.9}  (3 e^-2 sin 1 = {:.9})",
            later.unwrap_or(f64::NAN),
            3.0 * (-2.0_f64).exp() * 1.0_f64.sin()
        );
    }
    Ok(())
}
