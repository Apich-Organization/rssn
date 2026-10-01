//! Crank-Nicolson time stepping (ported from `physics_cnm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use num_complex::Complex;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_cnm::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 32,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn heat_1d_dirichlet_ends_leak_heat() {
    let n = 20;
    let res = solve_heat_equation_1d_cn(&vec![1.0; n], 0.1, 0.01, 0.1, 10);
    // The implementation pins u_0 = u_{n-1} = 0, so a uniform field loses heat through the ends.
    assert_eq!(res.len(), n);
    assert!(res[10] > 0.0);
    assert_eq!(res[0], 0.0);
    assert_eq!(res[n - 1], 0.0);
    assert!(res.iter().sum::<f64>() < n as f64 - 2.0);
    // Interior stays below the initial value and is symmetric about the middle.
    assert!(res.iter().all(|&v| (0.0..=1.0).contains(&v)));
    for i in 0..n / 2 {
        assert!((res[i] - res[n - 1 - i]).abs() < 1e-12);
    }
}

#[test]
fn heat_1d_sine_mode_decays_with_the_exact_crank_nicolson_factor() {
    let (n, d_coeff, dt, steps) = (51usize, 0.3, 0.01, 20usize);
    let dx = 1.0 / (n as f64 - 1.0);
    let u0: Vec<f64> = (0..n)
        .map(|i| (std::f64::consts::PI * i as f64 * dx).sin())
        .collect();
    let res = solve_heat_equation_1d_cn(&u0, dx, dt, d_coeff, steps);
    // Eigenvalue of the discrete Laplacian for the first sine mode.
    let lambda = -d_coeff * 2.0 * (1.0 - (std::f64::consts::PI * dx).cos()) / (dx * dx);
    let g = (1.0 + lambda * dt / 2.0) / (1.0 - lambda * dt / 2.0);
    let amp = g.powi(steps as i32);
    for i in 0..n {
        assert!((res[i] - amp * u0[i]).abs() < 1e-12, "node {i}");
    }
    // ...which approximates exp(-D pi^2 t).
    let cont = (-d_coeff * std::f64::consts::PI.powi(2) * steps as f64 * dt).exp();
    assert!((amp - cont).abs() < 1e-3);
}

#[test]
fn heat_1d_scenario_matches_the_analytic_decay() {
    let res = simulate_1d_heat_conduction_cn_scenario();
    assert_eq!(res.len(), 100);
    // D = 0.01, t = 50 * 0.001 = 0.05
    let decay = (-0.01 * std::f64::consts::PI.powi(2) * 0.05).exp();
    let dx = 1.0 / 99.0;
    for (i, &v) in res.iter().enumerate() {
        let want = decay * (std::f64::consts::PI * i as f64 * dx).sin();
        assert!((v - want).abs() < 1e-5, "node {i}: {v} vs {want}");
    }
}

#[test]
fn schrodinger_1d_conserves_the_norm() {
    let n = 100;
    let dx = 0.1;
    let mut psi0: Vec<Complex<f64>> = (0..n)
        .map(|i| {
            let x = (i as f64 - n as f64 / 2.0) * dx;
            Complex::new((-x * x).exp(), 0.0)
        })
        .collect();
    let norm = (psi0.iter().map(|p| p.norm_sqr()).sum::<f64>() * dx).sqrt();
    for p in psi0.iter_mut() {
        *p /= norm;
    }
    let res = solve_schrodinger_1d_cn(&psi0, &vec![0.0; n], dx, 0.01, 10);
    let final_norm = (res.iter().map(|p| p.norm_sqr()).sum::<f64>() * dx).sqrt();
    assert_approx_eq!(final_norm, 1.0, 1e-10);
}

#[test]
fn schrodinger_1d_eigenstate_only_acquires_the_exact_crank_nicolson_phase() {
    // Box eigenstate j = 2 with V = 0 (hbar = m = 1): psi -> ((1 - i E dt / 2) / (1 + i E dt / 2))^steps psi.
    let (n, j, dt, steps) = (41usize, 2usize, 0.02, 25usize);
    let dx = 0.05;
    let psi0: Vec<Complex<f64>> = (0..n)
        .map(|i| {
            Complex::new(
                (std::f64::consts::PI * (j * i) as f64 / (n as f64 - 1.0)).sin(),
                0.0,
            )
        })
        .collect();
    let energy = (1.0 - (j as f64 * std::f64::consts::PI / (n as f64 - 1.0)).cos()) / (dx * dx);
    let res = solve_schrodinger_1d_cn(&psi0, &vec![0.0; n], dx, dt, steps);
    let e = Complex::new(0.0, energy * dt / 2.0);
    let factor = ((Complex::new(1.0, 0.0) - e) / (Complex::new(1.0, 0.0) + e)).powu(steps as u32);
    assert!((factor.norm() - 1.0).abs() < 1e-12);
    for i in 0..n {
        assert!((res[i] - factor * psi0[i]).norm() < 1e-10, "node {i}");
    }
}

#[test]
fn schrodinger_1d_with_a_constant_potential_only_shifts_the_phase() {
    // A constant potential V0 multiplies the state by the same unimodular factor everywhere (interior).
    let n = 31;
    let psi0: Vec<Complex<f64>> = (0..n)
        .map(|i| {
            Complex::new(
                (std::f64::consts::PI * i as f64 / (n as f64 - 1.0)).sin(),
                0.0,
            )
        })
        .collect();
    let free = solve_schrodinger_1d_cn(&psi0, &vec![0.0; n], 0.1, 0.01, 5);
    let shifted = solve_schrodinger_1d_cn(&psi0, &vec![3.0; n], 0.1, 0.01, 5);
    let ratio = shifted[n / 2] / free[n / 2];
    assert!((ratio.norm() - 1.0).abs() < 1e-9);
    for i in 1..n - 1 {
        assert!((shifted[i] / free[i] - ratio).norm() < 1e-9, "node {i}");
    }
}

#[test]
fn heat_2d_adi_decays_a_product_mode_with_the_exact_factor() {
    let (nx, ny, steps) = (17usize, 17usize, 6usize);
    let (dx, dy, dt, d_coeff) = (1.0 / 16.0, 1.0 / 16.0, 0.01, 0.4);
    let pi = std::f64::consts::PI;
    let mut u0 = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            u0[j * nx + i] = (pi * i as f64 * dx).sin() * (pi * j as f64 * dy).sin();
        }
    }
    let cfg = HeatEquationSolverConfig {
        nx,
        ny,
        dx,
        dy,
        dt,
        d_coeff,
        steps,
    };
    let res = solve_heat_equation_2d_cn_adi(&u0, &cfg);
    let lx = -d_coeff * 2.0 * (1.0 - (pi * dx).cos()) / (dx * dx);
    let ly = -d_coeff * 2.0 * (1.0 - (pi * dy).cos()) / (dy * dy);
    let (ax, ay) = (lx * dt / 2.0, ly * dt / 2.0);
    // Peaceman-Rachford amplification factor of the mode.
    let g = (1.0 + ax) * (1.0 + ay) / ((1.0 - ax) * (1.0 - ay));
    let amp = g.powi(steps as i32);
    for k in 0..nx * ny {
        assert!((res[k] - amp * u0[k]).abs() < 1e-11, "cell {k}");
    }
    assert!((amp - (-d_coeff * 2.0 * pi * pi * steps as f64 * dt).exp()).abs() < 5e-3);
}

#[test]
fn heat_2d_adi_scenario_is_symmetric_bounded_and_dissipative() {
    let res = simulate_2d_heat_conduction_cn_adi_scenario();
    assert_eq!(res.len(), 50 * 50);
    assert!(res.iter().all(|v| v.is_finite()));
    for j in 0..50 {
        for i in 0..50 {
            assert!(
                (res[j * 50 + i] - res[i * 50 + j]).abs() < 1e-9,
                "x/y symmetry at ({i}, {j})"
            );
        }
    }
    let max = res.iter().cloned().fold(f64::MIN, f64::max);
    assert!(
        max < 100.0 && max > 1.0,
        "peak {max} must lie below the initial 100"
    );
    // Zero Dirichlet boundary.
    for k in 0..50 {
        assert_eq!(
            (res[k], res[49 * 50 + k], res[k * 50], res[k * 50 + 49]),
            (0.0, 0.0, 0.0, 0.0)
        );
    }
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_heat_1d_stays_in_range(d_coeff in 0.01..1.0f64) {
        let res = solve_heat_equation_1d_cn(&[1.0; 10], 0.1, 0.001, d_coeff, 5);
        for &val in &res {
            prop_assert!(val.is_finite());
            prop_assert!((-1e-12..=1.0 + 1e-12).contains(&val));
        }
    }

    #[test]
    fn prop_heat_1d_is_linear_in_the_initial_data(scale in -5.0..5.0f64, seed in 0u64..100) {
        let base: Vec<f64> = (0..15).map(|i| (((i as u64 * 37 + seed) % 11) as f64) / 10.0).collect();
        let scaled: Vec<f64> = base.iter().map(|v| v * scale).collect();
        let a = solve_heat_equation_1d_cn(&base, 0.1, 0.01, 0.2, 8);
        let b = solve_heat_equation_1d_cn(&scaled, 0.1, 0.01, 0.2, 8);
        for (x, y) in a.iter().zip(&b) {
            prop_assert!((x * scale - y).abs() < 1e-10);
        }
    }

    #[test]
    fn prop_schrodinger_norm_is_conserved_for_any_potential(v0 in -5.0..5.0f64, dt in 0.001..0.05f64) {
        let n = 40;
        let psi0: Vec<Complex<f64>> = (0..n).map(|i| {
            let x = i as f64 - 20.0;
            Complex::new((-x * x / 20.0).exp(), (0.3 * x).sin() * (-x * x / 20.0).exp())
        }).collect();
        let v: Vec<f64> = (0..n).map(|i| v0 * (i as f64 / n as f64)).collect();
        let res = solve_schrodinger_1d_cn(&psi0, &v, 0.5, dt, 10);
        let (n0, n1): (f64, f64) = (psi0.iter().map(|p| p.norm_sqr()).sum(), res.iter().map(|p| p.norm_sqr()).sum());
        // Boundary values are pinned to zero, so the norm can only shrink by the (tiny) boundary mass.
        prop_assert!((n1 - n0).abs() < 1e-6 * n0, "{n0} -> {n1}");
    }
}
