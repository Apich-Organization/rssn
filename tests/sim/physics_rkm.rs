//! Runge-Kutta ODE solvers and classic systems (ported from `physics_rkm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_rkm::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// dy/dt = a y
struct Linear(f64);

impl OdeSystem for Linear {
    fn dim(&self) -> usize {
        1
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = self.0 * y[0];
    }
}

fn last(res: &[(f64, Vec<f64>)]) -> &(f64, Vec<f64>) {
    res.last().unwrap_or_else(|| panic!("empty trajectory"))
}

#[test]
fn rk4_reproduces_the_exponential() {
    let res = solve_rk4(&Linear(-0.5), &[2.0], (0.0, 1.0), 0.01);
    assert_eq!(res[0].0, 0.0);
    assert_eq!(res[0].1, vec![2.0]);
    let (t, y) = last(&res);
    assert_approx_eq!(*t, 1.0, 1e-9);
    assert_approx_eq!(y[0], 2.0 * (-0.5f64).exp(), 1e-9);
}

#[test]
fn rk4_is_fourth_order() {
    let err = |dt: f64| {
        let res = solve_rk4(&Linear(1.0), &[1.0], (0.0, 1.0), dt);
        (last(&res).1[0] - 1f64.exp()).abs()
    };
    let ratio = err(0.1) / err(0.05);
    assert!((12.0..20.0).contains(&ratio), "ratio = {ratio}");
}

#[test]
fn rk4_damped_oscillator_matches_the_closed_form() {
    let (omega, zeta) = (1.0f64, 0.15f64);
    let sys = DampedOscillatorSystem { omega, zeta };
    let res = solve_rk4(&sys, &[1.0, 0.0], (0.0, 10.0), 0.01);
    let wd = omega * (1.0 - zeta * zeta).sqrt();
    let (t, y) = last(&res);
    let exact = (-zeta * omega * t).exp() * ((wd * t).cos() + zeta * omega / wd * (wd * t).sin());
    assert_approx_eq!(y[0], exact, 1e-7);
}

#[test]
fn adaptive_solvers_reach_the_exponential() {
    let sys = Linear(1.0);
    let tol = (1e-8, 1e-8);
    let e = 1f64.exp();
    let dp = DormandPrince54::new().solve(&sys, &[1.0], (0.0, 1.0), 0.01, tol);
    let ck = CashKarp45::default().solve(&sys, &[1.0], (0.0, 1.0), 0.01, tol);
    let bs = BogackiShampine23::default().solve(&sys, &[1.0], (0.0, 1.0), 0.01, tol);
    assert!((last(&dp).1[0] - e).abs() < 1e-6);
    assert!((last(&ck).1[0] - e).abs() < 1e-6);
    assert!((last(&bs).1[0] - e).abs() < 1e-4);
    for r in [&dp, &ck, &bs] {
        assert_approx_eq!(last(r).0, 1.0, 1e-9);
        // Time stamps strictly increase.
        assert!(r.windows(2).all(|w| w[1].0 > w[0].0));
    }
}

#[test]
fn tighter_tolerance_gives_smaller_error_and_more_steps() {
    let sys = VanDerPolSystem { mu: 1.0 };
    let dp = DormandPrince54::new();
    let loose = dp.solve(&sys, &[2.0, 0.0], (0.0, 5.0), 0.1, (1e-3, 1e-3));
    let tight = dp.solve(&sys, &[2.0, 0.0], (0.0, 5.0), 0.1, (1e-9, 1e-9));
    let reference = solve_rk4(&sys, &[2.0, 0.0], (0.0, 5.0), 1e-3);
    let r = &last(&reference).1;
    let d = |a: &[f64]| ((a[0] - r[0]).powi(2) + (a[1] - r[1]).powi(2)).sqrt();
    assert!(tight.len() > loose.len());
    assert!(
        d(&last(&tight).1) < 1e-6,
        "tight error {}",
        d(&last(&tight).1)
    );
    assert!(d(&last(&tight).1) < d(&last(&loose).1));
}

#[test]
fn pendulum_conserves_energy_with_dormand_prince() {
    let sys = PendulumSystem { g: 9.81, l: 1.0 };
    let energy = |y: &[f64]| 0.5 * y[1] * y[1] - 9.81 * y[0].cos();
    let res = DormandPrince54::new().solve(&sys, &[0.5, 0.0], (0.0, 5.0), 0.01, (1e-10, 1e-10));
    let e0 = energy(&res[0].1);
    for (_, y) in &res {
        assert!((energy(y) - e0).abs() < 1e-6);
    }
}

#[test]
fn pendulum_small_angle_period_matches_two_pi_sqrt_l_over_g() {
    let sys = PendulumSystem { g: 9.81, l: 1.0 };
    let period = 2.0 * std::f64::consts::PI / 9.81f64.sqrt();
    let res = solve_rk4(&sys, &[0.01, 0.0], (0.0, period), 1e-3);
    let y = &last(&res).1;
    // Back near the start after one small-angle period (nonlinear correction ~ theta^2 / 16).
    assert!((y[0] - 0.01).abs() < 1e-5, "theta = {}", y[0]);
}

#[test]
fn lorenz_rhs_at_a_known_point() {
    let sys = LorenzSystem {
        sigma: 10.0,
        rho: 28.0,
        beta: 8.0 / 3.0,
    };
    assert_eq!(sys.dim(), 3);
    let mut dy = [0.0; 3];
    sys.eval(0.0, &[1.0, 1.0, 1.0], &mut dy);
    assert_approx_eq!(dy[0], 0.0);
    assert_approx_eq!(dy[1], 26.0);
    assert_approx_eq!(dy[2], 1.0 - 8.0 / 3.0);
}

#[test]
fn lotka_volterra_conserves_its_invariant() {
    let sys = LotkaVolterraSystem {
        alpha: 1.5,
        beta: 1.0,
        delta: 1.0,
        gamma: 3.0,
    };
    let v = |y: &[f64]| 1.0 * y[0] - 3.0 * y[0].ln() + 1.0 * y[1] - 1.5 * y[1].ln();
    let res = DormandPrince54::new().solve(&sys, &[10.0, 5.0], (0.0, 2.0), 0.001, (1e-10, 1e-10));
    let v0 = v(&res[0].1);
    for (_, y) in &res {
        assert!((v(y) - v0).abs() < 1e-6);
    }
}

#[test]
fn vanderpol_rhs_and_limit_cycle_amplitude() {
    let sys = VanDerPolSystem { mu: 1.0 };
    let mut dy = [0.0; 2];
    sys.eval(0.0, &[2.0, 1.0], &mut dy);
    assert_approx_eq!(dy[0], 1.0);
    assert_approx_eq!(dy[1], (1.0 - 4.0) * 1.0 - 2.0);
    let res = CashKarp45::default().solve(&sys, &[0.1, 0.0], (0.0, 30.0), 0.1, (1e-8, 1e-8));
    let amp = res
        .iter()
        .skip(res.len() / 2)
        .map(|(_, y)| y[0].abs())
        .fold(0.0, f64::max);
    assert!((amp - 2.0).abs() < 0.05, "amplitude {amp}");
}

#[test]
fn scenarios_start_at_their_initial_conditions_and_reach_the_end() {
    let lz = simulate_lorenz_attractor_scenario();
    assert_eq!(lz[0], (0.0, vec![1.0, 1.0, 1.0]));
    assert_approx_eq!(last(&lz).0, 50.0, 1e-6);
    // The Lorenz attractor is bounded.
    assert!(lz.iter().all(|(_, y)| y.iter().all(|v| v.abs() < 100.0)));

    let dm = simulate_damped_oscillator_scenario();
    assert_eq!(dm[0], (0.0, vec![1.0, 0.0]));
    assert_approx_eq!(last(&dm).0, 40.0, 1e-6);
    assert!(last(&dm).1[0].abs() < 0.01);

    let vp = simulate_vanderpol_scenario();
    assert_eq!(vp[0], (0.0, vec![2.0, 0.0]));
    assert_approx_eq!(last(&vp).0, 20.0, 1e-6);

    let lv = simulate_lotka_volterra_scenario();
    assert_eq!(lv[0], (0.0, vec![10.0, 5.0]));
    assert!(lv.iter().all(|(_, y)| y[0] > 0.0 && y[1] > 0.0));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_rk4_matches_exponential(y0 in -100.0f64..100.0, a in -1.0f64..1.0, dt in 0.001f64..0.05) {
        let res = solve_rk4(&Linear(a), &[y0], (0.0, 1.0), dt);
        let exact = y0 * a.exp();
        prop_assert!((last(&res).1[0] - exact).abs() < 1e-6 * (1.0 + y0.abs()));
    }

    #[test]
    fn prop_dormand_prince_meets_tolerance(y0 in -100.0f64..100.0, tol in 1e-9f64..1e-4) {
        let res = DormandPrince54::new().solve(&Linear(-0.5), &[y0], (0.0, 1.0), 0.1, (tol, tol));
        let exact = y0 * (-0.5f64).exp();
        prop_assert!((last(&res).1[0] - exact).abs() < tol * 100.0 * (1.0 + y0.abs()));
    }
}
