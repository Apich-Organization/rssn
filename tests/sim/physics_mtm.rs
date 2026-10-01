//! Multigrid Poisson solvers (ported from `physics_mtm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_mtm::*;
use std::f64::consts::PI;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn poisson_1d_multigrid_matches_the_parabola() {
    // -u'' = 2, u(0) = u(1) = 0  =>  u = x (1 - x).
    let n = 31;
    let u = solve_poisson_1d_multigrid(n, &vec![2.0; n], 12).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), n + 2);
    assert_eq!(u[0], 0.0);
    assert_eq!(u[n + 1], 0.0);
    let h = 1.0 / (n + 1) as f64;
    for (i, &v) in u.iter().enumerate() {
        let x = i as f64 * h;
        assert_approx_eq!(v, x * (1.0 - x), 1e-6);
    }
}

#[test]
fn poisson_1d_scenario_is_accurate() {
    let u = simulate_1d_poisson_multigrid_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), 129);
    for (i, &v) in u.iter().enumerate() {
        let x = i as f64 / 128.0;
        assert_approx_eq!(v, x * (1.0 - x), 1e-3);
    }
}

#[test]
fn poisson_1d_rejects_sizes_that_are_not_two_to_the_k_minus_one() {
    for n in [2usize, 6, 10, 20, 100] {
        assert!(
            solve_poisson_1d_multigrid(n, &vec![1.0; n], 3).is_err(),
            "n = {n}"
        );
    }
    assert!(solve_poisson_1d_multigrid(15, &[1.0; 15], 1).is_ok());
}

#[test]
fn more_v_cycles_reduce_the_error() {
    let n = 63;
    let f = vec![2.0; n];
    let err = |cycles| {
        let u = solve_poisson_1d_multigrid(n, &f, cycles).unwrap_or_else(|e| panic!("{e}"));
        let h = 1.0 / (n + 1) as f64;
        u.iter()
            .enumerate()
            .map(|(i, &v)| (v - (i as f64 * h) * (1.0 - i as f64 * h)).abs())
            .fold(0.0, f64::max)
    };
    let (e1, e3, e8) = (err(1), err(3), err(8));
    assert!(e3 < e1 && e8 < e3, "{e1} {e3} {e8}");
    assert!(e8 < 1e-6);
}

#[test]
fn poisson_2d_multigrid_matches_the_sine_solution() {
    let u = simulate_2d_poisson_multigrid_scenario().unwrap_or_else(|e| panic!("{e}"));
    let n = 33;
    assert_eq!(u.len(), n * n);
    let h = 1.0 / (n - 1) as f64;
    for i in 0..n {
        for j in 0..n {
            let exact = (PI * i as f64 * h).sin() * (PI * j as f64 * h).sin();
            assert!(
                (u[i * n + j] - exact).abs() < 5e-3,
                "({i},{j}): {} vs {exact}",
                u[i * n + j]
            );
        }
    }
}

#[test]
fn poisson_2d_rejects_bad_sizes() {
    assert!(solve_poisson_2d_multigrid(10, &vec![0.0; 100], 2).is_err());
    assert!(solve_poisson_2d_multigrid(9, &vec![0.0; 81], 2).is_ok());
}

#[test]
fn grid_helpers_are_not_public_but_zero_forcing_gives_zero() {
    let u = solve_poisson_1d_multigrid(7, &[0.0; 7], 4).unwrap_or_else(|e| panic!("{e}"));
    assert!(u.iter().all(|&v| v == 0.0));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_poisson_1d_is_linear_in_the_force(scale in 0.1f64..10.0) {
        let n = 31;
        let f1 = vec![1.0; n];
        let f2: Vec<f64> = f1.iter().map(|x| x * scale).collect();
        let u1 = solve_poisson_1d_multigrid(n, &f1, 5).unwrap_or_else(|e| panic!("{e}"));
        let u2 = solve_poisson_1d_multigrid(n, &f2, 5).unwrap_or_else(|e| panic!("{e}"));
        for i in 0..u1.len() {
            prop_assert!((u1[i] * scale - u2[i]).abs() < 1e-9 * (1.0 + scale));
        }
    }

    #[test]
    fn prop_positive_force_gives_nonnegative_solution(k in 2u32..6, c in 0.1f64..5.0) {
        let n = 2usize.pow(k) - 1;
        let u = solve_poisson_1d_multigrid(n, &vec![c; n], 10).unwrap_or_else(|e| panic!("{e}"));
        prop_assert!(u.iter().all(|&v| v >= -1e-12));
    }
}
