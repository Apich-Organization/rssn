//! ODE integrators (ported from the old `numerical_ode_test.rs` and
//! `numerical/ode.rs`; the `Expr`-based right-hand sides are now closures).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::ode::{OdeSolverMethod, solve_adaptive, solve_fixed};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn grow(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = y[0];
}

fn oscillator(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = y[1];
    dy[1] = -y[0];
}

fn last(f: impl Fn(f64, &[f64], &mut [f64]), y0: &[f64], span: (f64, f64), steps: usize, m: OdeSolverMethod) -> Vec<f64> {
    solve_fixed(f, y0, span, steps, m)
        .unwrap_or_else(|e| panic!("{e}"))
        .last()
        .to_vec()
}

#[test]
fn rk4_exponential() {
    let y = last(grow, &[1.0], (0.0, 1.0), 100, OdeSolverMethod::RungeKutta4);
    assert_approx_eq!(y[0], std::f64::consts::E, 1e-8);
}

#[test]
fn euler_exponential() {
    let y = last(grow, &[1.0], (0.0, 1.0), 1000, OdeSolverMethod::Euler);
    // Euler with h = 1e-3: (1 + h)^1000 ~ e - e*h/2
    assert_approx_eq!(y[0], std::f64::consts::E, 2e-3);
    assert_approx_eq!(y[0], 1.001_f64.powi(1000), 1e-9);
}

#[test]
fn heun_exponential() {
    let y = last(grow, &[1.0], (0.0, 1.0), 100, OdeSolverMethod::Heun);
    assert_approx_eq!(y[0], std::f64::consts::E, 1e-3);
}

#[test]
fn rk4_harmonic_oscillator_half_period() {
    // y0 = cos t, y1 = -sin t; at t = pi: (-1, 0)
    let y = last(oscillator, &[1.0, 0.0], (0.0, std::f64::consts::PI), 100, OdeSolverMethod::RungeKutta4);
    assert_approx_eq!(y[0], -1.0, 1e-5);
    assert_approx_eq!(y[1], 0.0, 1e-5);
}

#[test]
fn rk4_conserves_oscillator_energy() {
    let tr = solve_fixed(oscillator, &[1.0, 0.0], (0.0, 10.0), 1000, OdeSolverMethod::RungeKutta4)
        .unwrap_or_else(|e| panic!("{e}"));
    for y in &tr.y {
        let energy = 0.5 * (y[0] * y[0] + y[1] * y[1]);
        assert!((energy - 0.5).abs() < 1e-8, "energy drifted to {energy}");
    }
}

#[test]
fn trajectory_grid_is_uniform_and_hits_endpoints() {
    let tr = solve_fixed(grow, &[1.0], (1.0, 3.0), 20, OdeSolverMethod::Heun)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(tr.t.len(), 21);
    assert_eq!(tr.y.len(), 21);
    assert_approx_eq!(tr.t[0], 1.0, 1e-15);
    assert_approx_eq!(tr.t[20], 3.0, 1e-12);
    assert_approx_eq!(tr.y[0][0], 1.0, 1e-15);
}

#[test]
fn rk4_reports_overflow_of_finite_time_blowup() {
    // y' = y^2, y(0) = 1 blows up at t = 1; integrating to t = 2 must fail.
    let r = solve_fixed(|_t, y, dy| dy[0] = y[0] * y[0], &[1.0], (0.0, 2.0), 100, OdeSolverMethod::RungeKutta4);
    assert_eq!(
        r.err().as_deref(),
        Some("Overflow or invalid value encountered during ODE solving.")
    );
}

#[test]
fn zero_steps_is_an_error() {
    assert!(solve_fixed(grow, &[1.0], (0.0, 1.0), 0, OdeSolverMethod::Euler).is_err());
}

#[test]
fn adaptive_matches_exponential() {
    let tr = solve_adaptive(grow, &[1.0], (0.0, 1.0), 1e-10, 1e-12, 10_000)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(tr.last()[0], std::f64::consts::E, 1e-8);
    assert_approx_eq!(*tr.t.last().unwrap_or(&f64::NAN), 1.0, 1e-12);
}

#[test]
fn adaptive_oscillator_backwards_in_time() {
    // Integrate from pi to 0 starting at (-1, 0): ends at (1, 0).
    let tr = solve_adaptive(oscillator, &[-1.0, 0.0], (std::f64::consts::PI, 0.0), 1e-10, 1e-12, 10_000)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(tr.last()[0], 1.0, 1e-7);
    assert_approx_eq!(tr.last()[1], 0.0, 1e-7);
}

#[test]
fn adaptive_step_budget_is_enforced() {
    let r = solve_adaptive(oscillator, &[1.0, 0.0], (0.0, 100.0), 1e-12, 1e-14, 3);
    assert!(r.is_err());
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_ode_constant_rhs_is_exact(a in -2.0..2.0f64, x_end in 0.1..2.0f64) {
        let y = last(move |_t, _y, dy| dy[0] = a, &[0.0], (0.0, x_end), 50, OdeSolverMethod::RungeKutta4);
        prop_assert!((y[0] - a * x_end).abs() < 1e-9);
    }

    #[test]
    fn prop_rk4_matches_exp_for_linear_decay_or_growth(a in -3.0..3.0f64, y0 in 0.5..2.0f64) {
        let y = last(move |_t, y, dy| dy[0] = a * y[0], &[y0], (0.0, 1.0), 200, OdeSolverMethod::RungeKutta4);
        let exact = y0 * a.exp();
        prop_assert!((y[0] - exact).abs() < 1e-7 * exact.abs().max(1.0), "{} vs {}", y[0], exact);
    }

    #[test]
    fn prop_adaptive_matches_exp(a in -3.0..3.0f64) {
        let tr = solve_adaptive(move |_t, y, dy| dy[0] = a * y[0], &[1.0], (0.0, 1.0), 1e-9, 1e-12, 100_000)
            .map_err(|e| TestCaseError::fail(e))?;
        prop_assert!((tr.last()[0] - a.exp()).abs() < 1e-6 * a.exp().max(1.0));
    }

    #[test]
    fn prop_heun_error_shrinks_when_steps_double(a in 0.5..2.0f64) {
        let exact = a.exp();
        let err = |n| (last(move |_t, y, dy| dy[0] = a * y[0], &[1.0], (0.0, 1.0), n, OdeSolverMethod::Heun)[0] - exact).abs();
        prop_assert!(err(200) < err(100));
    }
}
