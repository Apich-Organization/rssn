//! Finite element Poisson solvers (ported from `physics_fem_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_fem::*;
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
fn poisson_1d_constant_force_is_nodally_exact() {
    // -u'' = 2 on (0, 1), u(0) = u(1) = 0  =>  u = x (1 - x).
    let n = 10;
    let u = solve_poisson_1d(n, 1.0, |_| 2.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), n + 1);
    assert_eq!(u[0], 0.0);
    assert_eq!(u[n], 0.0);
    for (i, &v) in u.iter().enumerate() {
        let x = i as f64 / n as f64;
        assert_approx_eq!(v, x * (1.0 - x), 1e-6);
    }
}

#[test]
fn poisson_1d_scales_with_domain_length() {
    // -u'' = 1 on (0, L)  =>  u = x (L - x) / 2, so the midpoint is L^2 / 8.
    let u = solve_poisson_1d(20, 3.0, |_| 1.0).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(u[10], 9.0 / 8.0, 1e-6);
}

#[test]
fn poisson_1d_scenario_matches_the_parabola() {
    let u = simulate_1d_poisson_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), 51);
    assert_approx_eq!(u[25], 0.25, 1e-6);
}

#[test]
fn poisson_2d_sine_force_reproduces_sine_solution() {
    // -lap u = 2 pi^2 sin(pi x) sin(pi y)  =>  u = sin(pi x) sin(pi y).
    let n = 16;
    let u = solve_poisson_2d(n, n, |x, y| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin())
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), (n + 1) * (n + 1));
    let centre = u[(n / 2) * (n + 1) + n / 2];
    assert!((centre - 1.0).abs() < 0.02, "centre = {centre}");
    // Boundary nodes are pinned to zero.
    for i in 0..=n {
        assert_eq!(u[i], 0.0);
        assert_eq!(u[n * (n + 1) + i], 0.0);
    }
}

#[test]
fn poisson_2d_basic_shape_and_positivity() {
    let u = solve_poisson_2d(5, 5, |_, _| 2.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), 36);
    assert!(u.iter().all(|&v| v >= -1e-12));
    assert!(u[2 * 6 + 2] > 0.0);
}

#[test]
fn poisson_2d_scenario_peak_is_near_one() {
    let u = simulate_2d_poisson_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), 21 * 21);
    let max = u.iter().cloned().fold(f64::MIN, f64::max);
    assert!((max - 1.0).abs() < 0.02, "max = {max}");
}

#[test]
fn poisson_3d_sine_force_reproduces_sine_solution() {
    let n = 6;
    let u = solve_poisson_3d(n, |x, y, z| {
        3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin()
    })
    .unwrap_or_else(|e| panic!("{e}"));
    let m = n + 1;
    assert_eq!(u.len(), m * m * m);
    let max = u.iter().cloned().fold(f64::MIN, f64::max);
    assert!((max - 1.0).abs() < 0.15, "max = {max}");
    // The maximum sits at the centre node.
    let centre = u[(n / 2) * m * m + (n / 2) * m + n / 2];
    assert_approx_eq!(centre, max, 1e-9);
}

#[test]
fn poisson_3d_zero_force_gives_zero_solution() {
    let u = solve_poisson_3d(4, |_, _, _| 0.0).unwrap_or_else(|e| panic!("{e}"));
    assert!(u.iter().all(|&v| v.abs() < 1e-12));
}

#[test]
fn poisson_3d_symmetric_force_gives_symmetric_solution() {
    let n = 4;
    let m = n + 1;
    let v = solve_poisson_3d(n, |x, y, z| x + y + z).unwrap_or_else(|e| panic!("{e}"));
    for k in 0..m {
        for j in 0..m {
            for i in 0..m {
                assert_approx_eq!(v[(k * m + j) * m + i], v[(k * m + i) * m + j], 1e-8);
            }
        }
    }
}

#[test]
fn poisson_3d_is_second_order_accurate_against_the_exact_sine_solution() {
    // -lap u = 3 pi^2 sin sin sin has u = sin sin sin; the nodal error must fall ~4x per halving.
    let err = |n: usize| {
        let u = solve_poisson_3d(n, |x, y, z| {
            3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin()
        })
        .unwrap();
        let m = n + 1;
        let mut e = 0.0f64;
        for k in 0..m {
            for j in 0..m {
                for i in 0..m {
                    let h = 1.0 / n as f64;
                    let ex = (PI * i as f64 * h).sin()
                        * (PI * j as f64 * h).sin()
                        * (PI * k as f64 * h).sin();
                    e = e.max((u[(k * m + j) * m + i] - ex).abs());
                }
            }
        }
        e
    };
    let (e4, e8) = (err(4), err(8));
    assert!(e8 < 0.05, "e8 = {e8}");
    assert!(e4 / e8 > 3.0, "ratio {}", e4 / e8);
}

#[test]
fn poisson_3d_uniform_force_is_symmetric_under_all_axis_swaps() {
    let n = 4;
    let m = n + 1;
    let v = solve_poisson_3d(n, |_, _, _| 1.0).unwrap();
    for k in 0..m {
        for j in 0..m {
            for i in 0..m {
                let a = v[(k * m + j) * m + i];
                assert_approx_eq!(a, v[(i * m + j) * m + k], 1e-8);
                assert_approx_eq!(a, v[(j * m + k) * m + i], 1e-8);
            }
        }
    }
}

#[test]
fn poisson_3d_scenario_runs() {
    let u = simulate_3d_poisson_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.len(), 11 * 11 * 11);
    assert!(u.iter().all(|v| v.is_finite()));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_poisson_1d_symmetry(n in 4usize..40) {
        let u = solve_poisson_1d(n, 1.0, |_| 1.0).unwrap_or_else(|e| panic!("{e}"));
        for i in 0..=n {
            prop_assert!((u[i] - u[n - i]).abs() < 1e-8);
        }
    }

    #[test]
    fn prop_poisson_1d_is_linear_in_the_force(c in 0.1f64..10.0) {
        let u1 = solve_poisson_1d(16, 1.0, |_| 1.0).unwrap_or_else(|e| panic!("{e}"));
        let uc = solve_poisson_1d(16, 1.0, |_| c).unwrap_or_else(|e| panic!("{e}"));
        for i in 0..u1.len() {
            prop_assert!((c * u1[i] - uc[i]).abs() < 1e-6 * (1.0 + c));
        }
    }
}
