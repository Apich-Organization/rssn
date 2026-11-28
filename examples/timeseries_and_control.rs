//! Time-series analysis, state-space control, and generalised linear models.
//!
//! Run with `cargo run --locked --example timeseries_and_control`.

use rssn::kernels::control::StateSpace;
use rssn::kernels::control::mpc::MpcOptions;
use rssn::kernels::control::mpc::mpc_step;
use rssn::kernels::glm::Family;
use rssn::kernels::glm::Link;
use rssn::kernels::glm::fit;
use rssn::kernels::glm::predict;
use rssn::kernels::timeseries::acf;
use rssn::kernels::timeseries::difference;
use rssn::kernels::timeseries::pacf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Time series: AR(1)-like series, ACF, PACF, and differencing.
    let series: Vec<f64> = (0..64)
        .map(|t| 0.7_f64.powi(t % 8) + 0.05 * (t as f64 * 0.25).sin())
        .collect();
    let r = acf(&series, 4);
    let p = pacf(&series, 4);
    let d1 = difference(&series, 1);
    println!(
        "ACF(0..=3)  = [{:.3}, {:.3}, {:.3}, {:.3}]",
        r[0], r[1], r[2], r[3]
    );
    println!(
        "PACF(0..=3) = [{:.3}, {:.3}, {:.3}, {:.3}]",
        p[0], p[1], p[2], p[3]
    );
    println!("First difference length: {} -> {}", series.len(), d1.len());

    // 2. Control: discrete-time double integrator with Model Predictive Control.
    let sys = StateSpace::discrete(
        vec![vec![1.0, 0.1], vec![0.0, 1.0]],
        vec![vec![0.005], vec![0.1]],
        vec![vec![1.0, 0.0]],
        vec![vec![0.0]],
        0.1,
    )?;
    let step = mpc_step(
        &sys,
        &[1.0, 0.0],
        &[0.0, 0.0],
        &MpcOptions {
            horizon: 10,
            control_horizon: 4,
            q: vec![10.0, 1.0],
            r: vec![0.1],
            u_min: Some(vec![-2.0]),
            u_max: Some(vec![2.0]),
            x_min: None,
            x_max: None,
        },
    )?;
    println!(
        "MPC optimal first move u[0] = {:.4} (objective = {:.4}, iters = {})",
        step.u[0], step.objective, step.iterations
    );

    // 3. Generalised linear model: Poisson regression with log link.
    let x = vec![
        vec![0.0],
        vec![1.0],
        vec![2.0],
        vec![3.0],
        vec![4.0],
        vec![5.0],
    ];
    let y = vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
    let model = fit(&x, &y, Family::Poisson, Link::Log, None, 50, 1e-10)?;
    let preds = predict(&model, &x, Link::Log)?;
    println!(
        "Poisson GLM beta = [{:.4}, {:.4}], converged = {}, pred(5) = {:.2}",
        model.beta[0], model.beta[1], model.converged, preds[5]
    );

    Ok(())
}
