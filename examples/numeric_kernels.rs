//! The numeric side on its own: quadrature, root finding, an ODE
//! integration and a minimisation, with the closed-form value next to each.
//!
//! Run with `cargo run --release --example numeric_kernels`.

use std::f64::consts::PI;

use rssn::kernels::ode;
use rssn::kernels::optim;
use rssn::kernels::quadrature;
use rssn::kernels::rootfind;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- quadrature -----------------------------------------------------
    // int_0^pi sin x dx = 2, and the rule reports its own error estimate.
    let q = quadrature::integrate(|x: f64| x.sin(), 0.0, PI, 1e-12);
    println!(
        "int_0^pi sin x        = {:.15} +- {:.1e} ({} evals, converged {})",
        q.value, q.error, q.evaluations, q.converged
    );

    // An improper integral by the tanh-sinh rule: int e^-x^2 = sqrt(pi).
    let gauss = quadrature::integrate_infinite(|x| (-x * x).exp(), 1e-12);
    println!(
        "int_-oo^oo e^-x^2     = {:.15} +- {:.1e} (sqrt(pi) = {:.15})",
        gauss.value,
        gauss.error,
        PI.sqrt()
    );

    // --- root finding ---------------------------------------------------
    // The fixed point of cos is the root of cos(x) - x.
    let root = rootfind::brent(|x| x.cos() - x, 0.0, 1.0, 1e-14, 100)
        .map_err(|e| format!("brent: {e:?}"))?;
    println!(
        "\ncos(x) - x = 0        = {root:.15}  (residual {:.1e})",
        root.cos() - root
    );

    // --- an initial value problem ---------------------------------------
    // y' = -2y, y(0) = 1 has y(2) = e^-4; adaptive Dormand-Prince 5(4).
    let rhs = |_t: f64, y: &[f64], d: &mut [f64]| d[0] = -2.0 * y[0];
    let traj = ode::solve_adaptive(rhs, &[1.0], (0.0, 2.0), 1e-10, 1e-12, 10_000)?;
    let y2 = traj.last()[0];
    println!(
        "\ny' = -2y, y(2)        = {y2:.15}  (e^-4 = {:.15}, {} accepted points)",
        (-4.0_f64).exp(),
        traj.t.len()
    );

    // --- minimisation ---------------------------------------------------
    // Rosenbrock: a curved valley with its minimum at (1, 1).
    let rosenbrock = |p: &[f64]| {
        let (x, y) = (p[0], p[1]);
        (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2)
    };
    let best = optim::nelder_mead(rosenbrock, &[-1.2, 1.0], 0.1, 1e-12, 10_000);
    println!(
        "\nRosenbrock minimum    = ({:.9}, {:.9}), f = {:.2e} in {} iterations",
        best.x[0], best.x[1], best.fx, best.iterations
    );
    Ok(())
}
