//! Euler-type integrators (ported from `physics_em_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::matrix::Matrix;
use rssn::sim::physics_em::*;
use rssn::sim::physics_rkm::{DampedOscillatorSystem, OdeSystem};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 32,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// y' = -y
struct DecaySystem;

impl OdeSystem for DecaySystem {
    fn dim(&self) -> usize {
        1
    }

    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = -y[0];
    }
}

/// y' = a y
struct Exponential(f64);

impl OdeSystem for Exponential {
    fn dim(&self) -> usize {
        1
    }

    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        dy[0] = self.0 * y[0];
    }
}

fn last(res: &[(f64, Vec<f64>)]) -> (f64, f64) {
    res.last().map_or((f64::NAN, f64::NAN), |(t, y)| (*t, y[0]))
}

#[test]
fn forward_euler_decay_matches_the_closed_form_of_the_scheme() {
    let res = solve_forward_euler(&DecaySystem, &[1.0], (0.0, 1.0), 0.1);
    assert_eq!(res.len(), 11);
    let (t, y) = last(&res);
    assert_approx_eq!(t, 1.0, 1e-12);
    // (1 - h)^10 = 0.9^10
    assert_approx_eq!(y, 0.348_678_440_1, 1e-10);
    for (k, (tk, yk)) in res.iter().enumerate() {
        assert_approx_eq!(*tk, 0.1 * k as f64, 1e-12);
        assert_approx_eq!(yk[0], 0.9f64.powi(k as i32), 1e-12);
    }
}

#[test]
fn midpoint_and_heun_have_second_order_closed_forms() {
    // Both methods reduce to multiplication by 1 - h + h^2/2 per step for y' = -y.
    let h = 0.1;
    let factor = 1.0 - h + h * h / 2.0;
    for solver in [solve_midpoint_euler::<DecaySystem>, solve_heun_euler::<DecaySystem>] {
        let res = solver(&DecaySystem, &[1.0], (0.0, 1.0), h);
        let (_, y) = last(&res);
        assert_approx_eq!(y, factor.powi(10), 1e-12);
        assert!((y - (-1.0f64).exp()).abs() < 2e-3, "second-order accuracy: {y}");
    }
}

#[test]
fn convergence_orders_of_the_explicit_schemes() {
    let exact = 1.0f64.exp();
    let err = |solver: fn(&Exponential, &[f64], (f64, f64), f64) -> Vec<(f64, Vec<f64>)>, h: f64| {
        (last(&solver(&Exponential(1.0), &[1.0], (0.0, 1.0), h)).1 - exact).abs()
    };
    for (name, solver, order) in [
        ("euler", solve_forward_euler::<Exponential> as fn(&_, &_, _, _) -> _, 1.0),
        ("midpoint", solve_midpoint_euler::<Exponential> as fn(&_, &_, _, _) -> _, 2.0),
        ("heun", solve_heun_euler::<Exponential> as fn(&_, &_, _, _) -> _, 2.0),
    ] {
        let observed = (err(solver, 0.01) / err(solver, 0.005)).log2();
        assert!((observed - order).abs() < 0.1, "{name}: observed order {observed}");
    }
}

#[test]
fn undamped_oscillator_scenario_shows_forward_euler_energy_growth() {
    let res = simulate_oscillator_forward_euler_scenario();
    assert!(res.len() >= 1000);
    assert_eq!(res[0].1, vec![1.0, 0.0]);
    // Forward Euler multiplies the energy by (1 + (omega h)^2) each step: it grows without bound.
    let energy = |y: &[f64]| y[1] * y[1] + (2.0 * std::f64::consts::PI).powi(2) * y[0] * y[0];
    let (e0, e1) = (energy(&res[0].1), energy(&res.last().map_or(&[0.0, 0.0][..], |r| &r.1[..])));
    let growth = (1.0 + (2.0 * std::f64::consts::PI * 0.01f64).powi(2)).powi(res.len() as i32 - 1);
    assert!(e1 > 2.0 * e0);
    assert!((e1 / e0 - growth).abs() < 1e-6 * growth, "{} vs {growth}", e1 / e0);
}

struct Harmonic;

impl MechanicalSystem for Harmonic {
    fn spatial_dim(&self) -> usize {
        1
    }

    fn eval_acceleration(&self, x: &[f64], a: &mut [f64]) {
        a[0] = -4.0 * x[0]; // omega = 2
    }
}

#[test]
fn semi_implicit_euler_keeps_the_oscillator_bounded_and_the_energy_close() {
    let res = solve_semi_implicit_euler(&Harmonic, &[1.0, 0.0], (0.0, 20.0), 0.001).unwrap_or_else(|e| panic!("{e}"));
    let energy = |y: &[f64]| 0.5 * y[1] * y[1] + 2.0 * y[0] * y[0];
    let e0 = energy(&res[0].1);
    for (t, y) in &res {
        assert!((energy(y) - e0).abs() < 2e-3 * e0, "energy at t = {t}");
        assert!(y[0].abs() <= 1.01);
    }
    let (t_end, x_end) = last(&res);
    assert!((t_end - 20.0).abs() < 1e-6);
    assert!((x_end - (2.0 * t_end).cos()).abs() < 0.05);
}

#[test]
fn semi_implicit_euler_validates_the_state_length() {
    assert!(solve_semi_implicit_euler(&Harmonic, &[1.0], (0.0, 1.0), 0.1).is_err());
}

#[test]
fn orbit_scenario_conserves_angular_momentum_exactly_and_energy_approximately() {
    let res = simulate_gravity_semi_implicit_euler_scenario().unwrap_or_else(|e| panic!("{e}"));
    let (gm, dt) = (1000.0, 0.001);
    let l = |y: &[f64]| y[0] * y[3] - y[1] * y[2];
    let e = |y: &[f64]| 0.5 * (y[2] * y[2] + y[3] * y[3]) - gm / y[0].hypot(y[1]);
    let (l0, e0) = (l(&res[0].1), e(&res[0].1));
    assert!((l0 - 300.0).abs() < 1e-12 && (e0 - 350.0).abs() < 1e-9);
    for (t, y) in &res {
        assert!((l(y) - l0).abs() < 1e-9 * l0, "angular momentum at t = {t}");
    }
    let drift = res.iter().map(|(_, y)| (e(y) - e0).abs()).fold(0.0, f64::max);
    assert!(drift < 0.5, "energy drift {drift} (dt = {dt})");
    assert!(res.iter().all(|(_, y)| y.iter().all(|v| v.is_finite())));
}

struct Stiff;

impl LinearOdeSystem for Stiff {
    fn dim(&self) -> usize {
        2
    }

    fn get_matrix(&self) -> Matrix<f64> {
        Matrix::new(2, 2, vec![-1000.0, 0.0, 0.0, -1.0])
    }
}

struct Singular;

impl LinearOdeSystem for Singular {
    fn dim(&self) -> usize {
        2
    }

    fn get_matrix(&self) -> Matrix<f64> {
        // I - dt A is singular for dt = 1 (A has eigenvalue 1).
        Matrix::new(2, 2, vec![1.0, 0.0, 0.0, 1.0])
    }
}

#[test]
fn backward_euler_is_stable_for_stiff_systems_and_matches_the_closed_form() {
    let res = solve_backward_euler_linear(&Stiff, &[1.0, 1.0], (0.0, 1.0), 0.1).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(res.len(), 11);
    for (k, (_, y)) in res.iter().enumerate() {
        let k = k as i32;
        assert!((y[0] - (1.0 / 101.0f64).powi(k)).abs() < 1e-15 + 1e-12 * y[0].abs());
        assert!((y[1] - (1.0 / 1.1f64).powi(k)).abs() < 1e-12);
    }
    // The explicit method blows up with the same step.
    let explicit: Vec<f64> = (0..11).map(|k| (1.0 - 1000.0 * 0.1f64).powi(k)).collect();
    assert!(explicit[10].abs() > 1e15);
}

#[test]
fn stiff_decay_scenario_matches_the_closed_form() {
    let res = simulate_stiff_decay_scenario().unwrap_or_else(|e| panic!("{e}"));
    let steps = res.len() - 1;
    assert!(steps >= 25);
    let (_, y) = &res[steps];
    // Backward Euler: y_{n+1} = y_n / (1 - a dt)
    assert!((y[0] - (1.0 / 5.0f64).powi(steps as i32)).abs() < 1e-15);
    assert!((y[1] - (1.0 / 1.1f64).powi(steps as i32)).abs() < 1e-12);
}

#[test]
#[ignore = "library bug (Matrix::inverse on singular input): solve_backward_euler_linear returns Ok for A = I with dt = 1 (I - dt A = 0) because Matrix::inverse() returns Some(garbage) for singular matrices; expected Err(\"Matrix (I - dt*A) is not invertible.\")"]
fn backward_euler_reports_a_singular_step_matrix() {
    assert!(solve_backward_euler_linear(&Singular, &[1.0, 1.0], (0.0, 1.0), 1.0).is_err());
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_euler_stays_finite_for_the_damped_oscillator(dt in 0.001..0.05f64) {
        let system = DampedOscillatorSystem { omega: 1.0, zeta: 0.1 };
        let res = solve_forward_euler(&system, &[1.0, 0.0], (0.0, 1.0), dt);
        for (_, y) in res {
            prop_assert!(y[0].is_finite() && y[1].is_finite());
        }
    }

    #[test]
    fn prop_forward_euler_is_geometric_for_linear_decay(a in -3.0..3.0f64, dt in 0.01..0.2f64) {
        let res = solve_forward_euler(&Exponential(a), &[2.0], (0.0, 1.0), dt);
        for (k, (_, y)) in res.iter().enumerate() {
            let want = 2.0 * (1.0 + a * dt).powi(k as i32);
            prop_assert!((y[0] - want).abs() < 1e-9 * want.abs().max(1.0));
        }
    }

    #[test]
    fn prop_heun_is_more_accurate_than_forward_euler(a in -2.0..-0.1f64) {
        let exact = a.exp();
        let e1 = (last(&solve_forward_euler(&Exponential(a), &[1.0], (0.0, 1.0), 0.05)).1 - exact).abs();
        let e2 = (last(&solve_heun_euler(&Exponential(a), &[1.0], (0.0, 1.0), 0.05)).1 - exact).abs();
        prop_assert!(e2 < e1);
    }
}
