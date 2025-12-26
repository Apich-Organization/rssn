//! Comprehensive tests for physics RKM module.

use rssn::physics::physics_rkm::*;

#[test]

fn test_lorenz_simulate() {

    let results = simulate_lorenz_attractor_scenario();

    assert!(results.len() > 0);

    // Initial condition was [1, 1, 1] at t=0
    assert_eq!(results[0].0, 0.0);

    assert_eq!(results[0].1, vec![1.0, 1.0, 1.0]);
}

#[test]

fn test_damped_oscillator_simulate() {

    let results = simulate_damped_oscillator_scenario();

    assert!(results.len() > 0);

    assert_eq!(results[0].0, 0.0);

    assert_eq!(results[0].1, vec![1.0, 0.0]);
}

#[test]

fn test_vanderpol_simulate() {

    let results = simulate_vanderpol_scenario();

    assert!(results.len() > 0);

    assert_eq!(results[0].0, 0.0);

    assert_eq!(results[0].1, vec![2.0, 0.0]);
}

#[test]

fn test_lotka_volterra_simulate() {

    let results = simulate_lotka_volterra_scenario();

    assert!(results.len() > 0);

    assert_eq!(results[0].0, 0.0);

    assert_eq!(results[0].1, vec![10.0, 5.0]);
}

#[test]

fn test_pendulum_rk4() {

    let system = PendulumSystem { g: 9.81, l: 1.0 };

    let y0 = vec![0.1, 0.0];

    let t_span = (0.0, 1.0);

    let dt = 0.01;

    let results = solve_rk4(&system, &y0, t_span, dt);

    assert!(results.len() > 0);
}

#[test]

fn test_adaptive_solvers_accuracy() {

    // Simple linear ODE: dy/dt = y, y(0) = 1 => y(t) = e^t
    struct SimpleSystem;

    impl OdeSystem for SimpleSystem {
        fn dim(&self) -> usize {

            1
        }

        fn eval(
            &self,
            _t: f64,
            y: &[f64],
            dy: &mut [f64],
        ) {

            dy[0] = y[0];
        }
    }

    let system = SimpleSystem;

    let y0 = vec![1.0];

    let t_span = (0.0, 1.0);

    let dt_initial = 0.01;

    let tol = (1e-8, 1e-8);

    // Dormand-Prince
    let dp54 = DormandPrince54::new();

    let res_dp = dp54.solve(&system, &y0, t_span, dt_initial, tol);

    let last_dp = res_dp
        .last()
        .unwrap();

    assert!((last_dp.1[0] - 1.0f64.exp()).abs() < 1e-6);

    // Cash-Karp
    let ck45 = CashKarp45::default();

    let res_ck = ck45.solve(&system, &y0, t_span, dt_initial, tol);

    let last_ck = res_ck
        .last()
        .unwrap();

    assert!((last_ck.1[0] - 1.0f64.exp()).abs() < 1e-6);

    // Bogacki-Shampine
    let bs23 = BogackiShampine23::default();

    let res_bs = bs23.solve(&system, &y0, t_span, dt_initial, tol);

    let last_bs = res_bs
        .last()
        .unwrap();

    assert!((last_bs.1[0] - 1.0f64.exp()).abs() < 1e-4);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_rk4_simple_linear(y0_val in -100.0..100.0f64, dt in 0.001..0.1f64) {
            struct LinearSystem;
            impl OdeSystem for LinearSystem {
                fn dim(&self) -> usize { 1 }
                fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
                    dy[0] = -0.5 * y[0]; // dy/dt = -0.5y => y(t) = y0 * e^(-0.5t)
                }
            }

            let system = LinearSystem;
            let y0 = vec![y0_val];
            let t_span = (0.0, 1.0);
            let results = solve_rk4(&system, &y0, t_span, dt);

            if !results.is_empty() {
                let last = results.last().unwrap();
                let expected = y0_val * (-0.5f64).exp();
                println!("Expected: {}, Actual: {}", expected, last.1[0]);
                println!("Diff: {}", (last.1[0] - expected).abs());
                prop_assert!((last.1[0] - expected).abs() < 1e-4);
            }
        }

        #[test]
        fn prop_dp54_simple_linear(y0_val in -100.0..100.0f64, tol_val in 1e-9..1e-3f64) {
             struct LinearSystem;
            impl OdeSystem for LinearSystem {
                fn dim(&self) -> usize { 1 }
                fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
                    dy[0] = -0.5 * y[0];
                }
            }

            let system = LinearSystem;
            let y0 = vec![y0_val];
            let t_span = (0.0, 1.0);
            let solver = DormandPrince54::new();
            let results = solver.solve(&system, &y0, t_span, 0.1, (tol_val, tol_val));

            if !results.is_empty() {
                let last = results.last().unwrap();
                let expected = y0_val * (-0.5f64).exp();
                prop_assert!((last.1[0] - expected).abs() < tol_val * 100.0);
            }
        }
    }
}
