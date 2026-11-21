//! Finite-difference CFD kernels (ported from `numerical_physics_cfd_test.rs`).
//!
//! Tests for fluid properties, dimensionless numbers, and PDE solvers.

use rssn::kernels::matrix::Matrix;
use rssn::sim::physics_cfd::*;

// ============================================================================
// Fluid Properties Tests
// ============================================================================

#[test]

fn test_fluid_properties_new() {
    let fluid = FluidProperties::new(1000.0, 0.001, 0.6, 4200.0);

    assert_eq!(fluid.density, 1000.0);

    assert_eq!(fluid.dynamic_viscosity, 0.001);
}

#[test]

fn test_fluid_properties_air() {
    let air = FluidProperties::air();

    assert!((air.density - 1.204).abs() < 0.01);

    assert!(air.dynamic_viscosity > 1e-5 && air.dynamic_viscosity < 2e-5);
}

#[test]

fn test_fluid_properties_water() {
    let water = FluidProperties::water();

    assert!((water.density - 998.2).abs() < 1.0);

    assert!(water.dynamic_viscosity > 9e-4 && water.dynamic_viscosity < 1.1e-3);
}

#[test]

fn test_kinematic_viscosity() {
    let water = FluidProperties::water();

    let nu = water.kinematic_viscosity();

    // ν ≈ 1e-6 m²/s for water
    assert!(nu > 9e-7 && nu < 1.1e-6);
}

#[test]

fn test_thermal_diffusivity() {
    let water = FluidProperties::water();

    let alpha = water.thermal_diffusivity();

    // α ≈ 1.4e-7 m²/s for water
    assert!(alpha > 1e-7 && alpha < 2e-7);
}

#[test]

fn test_prandtl_number() {
    let water = FluidProperties::water();

    let pr = water.prandtl_number();

    // Pr ≈ 7 for water at 20°C
    assert!(pr > 6.0 && pr < 8.0);
}

// ============================================================================
// Dimensionless Numbers Tests
// ============================================================================

#[test]

fn test_reynolds_number() {
    // Flow at 1 m/s, length 1 m, water
    let re = reynolds_number(1.0, 1.0, 1e-6);

    assert!((re - 1e6).abs() < 1.0);
}

#[test]

fn test_mach_number() {
    // Aircraft at 250 m/s in air (c ≈ 340 m/s)
    let ma = mach_number(250.0, 340.0);

    assert!((ma - 0.735).abs() < 0.01);
}

#[test]

fn test_froude_number() {
    // Ship at 5 m/s, length 10 m
    let fr = froude_number(5.0, 10.0, 9.81);

    assert!(fr > 0.4 && fr < 0.6);
}

#[test]

fn test_cfl_number() {
    let cfl = cfl_number(1.0, 0.01, 0.1);

    assert!((cfl - 0.1).abs() < 1e-10);
}

#[test]

fn test_check_cfl_stability() {
    assert!(check_cfl_stability(1.0, 0.01, 0.1, 1.0));

    assert!(!check_cfl_stability(10.0, 0.1, 0.1, 1.0));
}

#[test]

fn test_diffusion_number() {
    let r = diffusion_number(0.01, 0.001, 0.01);

    // r = 0.01 * 0.001 / 0.0001 = 0.1
    assert!((r - 0.1).abs() < 1e-10);
}

// ============================================================================
// 1D Solver Tests
// ============================================================================

#[test]

fn test_solve_advection_1d() {
    // Initial step function
    let n = 20;

    let mut u0 = vec![0.0; n];

    for i in 5..10 {
        u0[i] = 1.0;
    }

    let results = solve_advection_1d(&u0, 1.0, 0.1, 0.05, 10);

    assert_eq!(results.len(), 11);

    // Check that the solution evolved
    assert_ne!(results[10], results[0]);
}

#[test]

fn test_solve_diffusion_1d() {
    // Initial peak
    let n = 20;

    let mut u0 = vec![0.0; n];

    u0[10] = 1.0;

    let results = solve_diffusion_1d(&u0, 0.1, 0.1, 0.01, 10);

    assert_eq!(results.len(), 11);

    // Peak should spread out
    assert!(results[10][10] < results[0][10]);

    assert!(results[10][9] > results[0][9]);
}

#[test]

fn test_solve_advection_diffusion_1d() {
    let n = 20;

    let mut u0 = vec![0.0; n];

    for i in 5..10 {
        u0[i] = 1.0;
    }

    let results = solve_advection_diffusion_1d(&u0, 0.5, 0.01, 0.1, 0.01, 10);

    assert_eq!(results.len(), 11);
}

#[test]

fn test_solve_burgers_1d() {
    // Initial smooth profile
    let n = 50;

    let mut u0 = vec![0.0; n];

    for i in 10..30 {
        u0[i] = 1.0 - ((i as f64 - 20.0).abs() / 10.0);
    }

    let results = solve_burgers_1d(&u0, 0.01, 0.02, 0.001, 50);

    assert_eq!(results.len(), 51);
}

// ============================================================================
// 2D Solver Tests
// ============================================================================

#[test]

fn test_solve_poisson_2d_jacobi() {
    let n = 10;

    let f = Matrix::zeros(n, n);

    let u0 = Matrix::zeros(n, n);

    let result = solve_poisson_2d_jacobi(&f, &u0, 0.1, 0.1, 100, 1e-6);

    assert_eq!(result.rows(), n);

    assert_eq!(result.cols(), n);
}

#[test]

fn test_solve_poisson_2d_gauss_seidel() {
    let n = 10;

    let f = Matrix::zeros(n, n);

    let u0 = Matrix::zeros(n, n);

    let result = solve_poisson_2d_gauss_seidel(&f, &u0, 0.1, 0.1, 100, 1e-6);

    assert_eq!(result.rows(), n);

    assert_eq!(result.cols(), n);
}

#[test]

fn test_solve_poisson_2d_sor() {
    let n = 10;

    let f = Matrix::zeros(n, n);

    let u0 = Matrix::zeros(n, n);

    let result = solve_poisson_2d_sor(&f, &u0, 0.1, 0.1, 1.5, 100, 1e-6);

    assert_eq!(result.rows(), n);

    assert_eq!(result.cols(), n);
}

// ============================================================================
// Vorticity and Stream Function Tests
// ============================================================================

#[test]

fn test_compute_vorticity() {
    let n = 10;

    let u = Matrix::zeros(n, n);

    let v = Matrix::zeros(n, n);

    let omega = compute_vorticity(&u, &v, 0.1, 0.1);

    assert_eq!(omega.rows(), n);

    assert_eq!(omega.cols(), n);
}

#[test]

fn test_compute_stream_function() {
    let n = 10;

    let omega = Matrix::zeros(n, n);

    let psi = compute_stream_function(&omega, 0.1, 0.1, 50, 1e-6);

    assert_eq!(psi.rows(), n);

    assert_eq!(psi.cols(), n);
}

#[test]

fn test_velocity_from_stream_function() {
    let n = 10;

    let psi = Matrix::zeros(n, n);

    let (u, v) = velocity_from_stream_function(&psi, 0.1, 0.1);

    assert_eq!(u.rows(), n);

    assert_eq!(v.cols(), n);
}

// ============================================================================
// Field Operations Tests
// ============================================================================

#[test]

fn test_compute_divergence() {
    let n = 10;

    let u = Matrix::zeros(n, n);

    let v = Matrix::zeros(n, n);

    let div = compute_divergence(&u, &v, 0.1, 0.1);

    assert_eq!(div.rows(), n);

    // Divergence of zero field is zero
    for i in 0..n {
        for j in 0..n {
            assert!(div.get(i, j).abs() < 1e-10);
        }
    }
}

#[test]

fn test_compute_gradient() {
    let n = 10;

    let p = Matrix::zeros(n, n);

    let (dp_dx, dp_dy) = compute_gradient(&p, 0.1, 0.1);

    assert_eq!(dp_dx.rows(), n);

    assert_eq!(dp_dy.cols(), n);
}

#[test]

fn test_compute_laplacian() {
    let n = 10;

    let f = Matrix::zeros(n, n);

    let lap = compute_laplacian(&f, 0.1, 0.1);

    assert_eq!(lap.rows(), n);

    assert_eq!(lap.cols(), n);
}

// ============================================================================
// Boundary Condition Tests
// ============================================================================

#[test]

fn test_apply_dirichlet_bc() {
    let n = 5;

    let mut field = Matrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            *field.get_mut(i, j) = 1.0;
        }
    }

    apply_dirichlet_bc(&mut field, 0.0);

    // Check boundaries are zero
    for i in 0..n {
        assert!(field.get(i, 0).abs() < 1e-10);

        assert!(field.get(i, n - 1).abs() < 1e-10);
    }

    for j in 0..n {
        assert!(field.get(0, j).abs() < 1e-10);

        assert!(field.get(n - 1, j).abs() < 1e-10);
    }

    // Interior should still be 1
    assert!((*field.get(2, 2) - 1.0).abs() < 1e-10);
}

#[test]

fn test_apply_neumann_bc() {
    let n = 5;

    let mut field = Matrix::zeros(n, n);

    *field.get_mut(2, 2) = 1.0;

    apply_neumann_bc(&mut field);

    // Boundaries should copy from neighbors
    // This is a simple test
    assert!(field.get(0, 2).abs() < 1e-10); // Copies from interior
}

// ============================================================================
// Utility Function Tests
// ============================================================================

#[test]

fn test_max_velocity_magnitude() {
    let n = 5;

    let mut u = Matrix::zeros(n, n);

    let mut v = Matrix::zeros(n, n);

    *u.get_mut(2, 2) = 3.0;

    *v.get_mut(2, 2) = 4.0;

    let max_vel = max_velocity_magnitude(&u, &v);

    assert!((max_vel - 5.0).abs() < 1e-10);
}

#[test]

fn test_l2_norm() {
    let n = 4;

    let mut field = Matrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            *field.get_mut(i, j) = 1.0;
        }
    }

    let norm = l2_norm(&field);

    assert!((norm - 1.0).abs() < 1e-10);
}

#[test]

fn test_max_abs() {
    let n = 5;

    let mut field = Matrix::zeros(n, n);

    *field.get_mut(2, 2) = -3.0;

    *field.get_mut(1, 1) = 2.0;

    let max_val = max_abs(&field);

    assert!((max_val - 3.0).abs() < 1e-10);
}

#[test]

fn test_lid_driven_cavity_simple() {
    let (psi, omega) = lid_driven_cavity_simple(10, 10, 100.0, 1.0, 10, 0.001);

    assert_eq!(psi.rows(), 10);

    assert_eq!(omega.cols(), 10);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use proptest::prelude::*;

    use super::*;

    proptest! {
        #[test]
        fn prop_reynolds_positive(v in 0.1..100.0f64, l in 0.1..10.0f64, nu in 1e-7..1e-3f64) {
            let re = reynolds_number(v, l, nu);
            prop_assert!(re > 0.0);
        }

        #[test]
        fn prop_cfl_positive(v in 0.1..10.0f64, dt in 1e-4..0.1f64, dx in 0.01..1.0f64) {
            let cfl = cfl_number(v, dt, dx);
            prop_assert!(cfl >= 0.0);
        }

        #[test]
        fn prop_diffusion_number_positive(alpha in 1e-6..1.0f64, dt in 1e-4..0.1f64, dx in 0.01..1.0f64) {
            let r = diffusion_number(alpha, dt, dx);
            prop_assert!(r >= 0.0);
        }

        #[test]
        fn prop_prandtl_positive(mu in 1e-5..1.0f64, cp in 100.0..5000.0f64, k in 0.01..100.0f64) {
            let fluid = FluidProperties::new(1000.0, mu, k, cp);
            prop_assert!(fluid.prandtl_number() > 0.0);
        }
    }
}

// ============================================================================
// Added: analytic and conservation checks
// ============================================================================

mod strengthened {
    use std::f64::consts::PI;

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::matrix::Matrix;
    use rssn::sim::physics_cfd::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            cases: 32,
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    /// Matrix sampled on a uniform grid over the unit square (row index = x, column index = y).
    fn grid(
        nx: usize,
        ny: usize,
        f: impl Fn(f64, f64) -> f64,
    ) -> Matrix<f64> {
        let mut m = Matrix::zeros(nx, ny);
        for i in 0..nx {
            for j in 0..ny {
                *m.get_mut(i, j) = f(i as f64 / (nx as f64 - 1.0), j as f64 / (ny as f64 - 1.0));
            }
        }
        m
    }

    fn interior_max_err(
        a: &Matrix<f64>,
        b: &Matrix<f64>,
    ) -> f64 {
        let mut e = 0.0f64;
        for i in 1..a.rows() - 1 {
            for j in 1..a.cols() - 1 {
                e = e.max((a.get(i, j) - b.get(i, j)).abs());
            }
        }
        e
    }

    #[test]
    fn fluid_property_formulas() {
        let f = FluidProperties::new(1000.0, 0.002, 0.5, 4000.0);
        assert_eq!(f.kinematic_viscosity(), 2e-6);
        assert!((f.thermal_diffusivity() - 0.5 / (1000.0 * 4000.0)).abs() < 1e-20);
        assert!((f.prandtl_number() - 16.0).abs() < 1e-12); // mu cp / k
        // Air at 20 C: nu = 1.516e-5, Pr = 0.71
        let air = FluidProperties::air();
        assert!((air.kinematic_viscosity() - 1.516e-5).abs() < 2e-8);
        assert!((air.prandtl_number() - 0.714).abs() < 0.01);
        let water = FluidProperties::water();
        assert!((water.prandtl_number() - 7.0).abs() < 0.1);
    }

    #[test]
    fn dimensionless_numbers_exact_values() {
        assert_eq!(reynolds_number(2.0, 3.0, 1e-6), 6e6);
        assert_eq!(mach_number(340.0, 340.0), 1.0);
        assert!(
            (froude_number(3.132_091_952_673_165, 10.0, 9.81) - 0.316_227_766_016_837_9).abs()
                < 1e-12
        );
        assert_eq!(cfl_number(-2.0, 0.1, 0.5), 0.4); // uses |v|
        assert!(check_cfl_stability(1.0, 0.1, 0.1, 1.0)); // boundary: cfl == max
        assert!(!check_cfl_stability(1.0, 0.11, 0.1, 1.0));
    }

    #[test]
    fn advection_at_cfl_one_translates_the_profile_exactly() {
        let n = 30;
        let mut u0 = vec![0.0; n];
        u0[5] = 1.0;
        u0[6] = 2.0;
        u0[7] = 3.0;
        let res = solve_advection_1d(&u0, 1.0, 0.1, 0.1, 10);
        assert_eq!(res.len(), 11);
        for k in 0..=10usize {
            for i in 0..n {
                let want = if i > k && i - k < n {
                    u0[i - k]
                } else {
                    0.0
                };
                if i >= 1 && i < n - 1 {
                    assert!((res[k][i] - want).abs() < 1e-12, "step {k}, cell {i}");
                }
            }
        }
    }

    #[test]
    fn advection_upwinds_in_the_direction_of_the_wind() {
        let n = 20;
        let mut u0 = vec![0.0; n];
        u0[10] = 1.0;
        let right = solve_advection_1d(&u0, 1.0, 0.1, 0.05, 1);
        let left = solve_advection_1d(&u0, -1.0, 0.1, 0.05, 1);
        assert!(right[1][11] > 0.0 && right[1][9] == 0.0);
        assert!(left[1][9] > 0.0 && left[1][11] == 0.0);
    }

    #[test]
    fn advection_conserves_the_periodic_ring_sum_and_keeps_bounds() {
        let n = 40;
        let mut u0 = vec![0.0; n];
        for i in 10..20 {
            u0[i] = (i as f64 - 9.0) * 0.1;
        }
        let total = |u: &[f64]| u[1..n - 1].iter().sum::<f64>();
        let res = solve_advection_1d(&u0, 1.0, 0.1, 0.05, 200);
        let t0 = total(&res[0]);
        for u in &res {
            assert!((total(u) - t0).abs() < 1e-10);
            assert!(
                u.iter().all(|&v| (-1e-12..=1.0 + 1e-12).contains(&v)),
                "upwind must stay monotone for CFL <= 1"
            );
        }
    }

    #[test]
    fn diffusion_of_a_sine_mode_matches_the_ftcs_amplification_factor() {
        let (n, alpha, dt) = (41usize, 0.1, 0.002);
        let dx = 1.0 / (n as f64 - 1.0);
        let u0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let steps = 50;
        let res = solve_diffusion_1d(&u0, alpha, dx, dt, steps);
        let r = alpha * dt / (dx * dx);
        assert!(r <= 0.5);
        let amp = (1.0 - 4.0 * r * (PI * dx / 2.0).sin().powi(2)).powi(steps as i32);
        for i in 1..n - 1 {
            assert!((res[steps][i] - amp * u0[i]).abs() < 1e-12, "cell {i}");
        }
        // And it approximates the continuous decay exp(-alpha pi^2 t).
        let cont = (-alpha * PI * PI * steps as f64 * dt).exp();
        assert!((amp - cont).abs() < 1e-3);
    }

    #[test]
    fn diffusion_conserves_mass_away_from_the_boundary_and_respects_the_maximum_principle() {
        let n = 41;
        let mut u0 = vec![0.0; n];
        for i in 15..25 {
            u0[i] = 1.0;
        }
        let res = solve_diffusion_1d(&u0, 0.01, 0.05, 0.05, 10); // r = 0.2; the spread stays far from the fixed ends
        let mass0: f64 = u0.iter().sum();
        for u in &res {
            assert!((u.iter().sum::<f64>() - mass0).abs() < 1e-12);
            assert!(u.iter().all(|&v| (-1e-15..=1.0 + 1e-15).contains(&v)));
        }
        // Fixed boundary values are kept.
        assert!(res.iter().all(|u| u[0] == 0.0 && u[n - 1] == 0.0));
    }

    #[test]
    fn advection_diffusion_reduces_to_its_limits() {
        let n = 30;
        let mut u0 = vec![0.0; n];
        for i in 8..14 {
            u0[i] = 1.0;
        }
        let pure_d = solve_diffusion_1d(&u0, 0.02, 0.1, 0.05, 15);
        let ad_no_wind = solve_advection_diffusion_1d(&u0, 0.0, 0.02, 0.1, 0.05, 15);
        for (a, b) in pure_d.iter().zip(&ad_no_wind) {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() < 1e-14);
            }
        }
        let pure_a = solve_advection_1d(&u0, 1.0, 0.1, 0.05, 6);
        let ad_no_diff = solve_advection_diffusion_1d(&u0, 1.0, 0.0, 0.1, 0.05, 6);
        for (a, b) in pure_a.iter().zip(&ad_no_diff) {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn burgers_keeps_constants_and_stays_bounded() {
        let n = 50;
        let res = solve_burgers_1d(&vec![0.7; n], 0.01, 0.02, 0.001, 20);
        assert!(
            res.iter()
                .all(|u| u.iter().all(|&v| (v - 0.7).abs() < 1e-12))
        );
        let mut u0 = vec![0.0; n];
        for i in 10..30 {
            u0[i] = 1.0 - ((i as f64 - 20.0).abs() / 10.0);
        }
        let res = solve_burgers_1d(&u0, 0.01, 0.02, 0.001, 200);
        for u in &res {
            assert!(
                u.iter().all(|&v| (-1e-12..=1.0 + 1e-12).contains(&v)),
                "max principle violated"
            );
        }
        // Viscosity smooths the pulse: its L2 norm decreases.
        let l2 = |u: &[f64]| u.iter().map(|v| v * v).sum::<f64>();
        assert!(l2(&res[200]) < l2(&res[0]));
    }

    #[test]
    fn poisson_solvers_agree_with_the_exact_discrete_solution() {
        // laplace(u) = -2 pi^2 sin(pi x) sin(pi y), u = 0 on the boundary.
        let n = 21;
        let d = 1.0 / (n as f64 - 1.0);
        let f = grid(n, n, |x, y| {
            -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin()
        });
        let u0 = Matrix::zeros(n, n);
        // The 5-point Laplacian has eigenvalue -4 (1 - cos(pi d)) / d^2 on this mode.
        let amp = 2.0 * PI * PI * d * d / (4.0 * (1.0 - (PI * d).cos()));
        let exact = grid(n, n, |x, y| amp * (PI * x).sin() * (PI * y).sin());

        let jac = solve_poisson_2d_jacobi(&f, &u0, d, d, 6000, 1e-13);
        let gs = solve_poisson_2d_gauss_seidel(&f, &u0, d, d, 6000, 1e-13);
        let sor = solve_poisson_2d_sor(&f, &u0, d, d, 1.7, 3000, 1e-13);
        for (name, sol) in [("jacobi", &jac), ("gauss-seidel", &gs), ("sor", &sor)] {
            let e = interior_max_err(sol, &exact);
            assert!(e < 1e-8, "{name}: max error vs exact discrete solution {e}");
            // And within second-order truncation error of the continuous solution.
            let cont = grid(n, n, |x, y| (PI * x).sin() * (PI * y).sin());
            assert!(interior_max_err(sol, &cont) < 5e-3, "{name}");
        }
    }

    #[test]
    fn poisson_solvers_leave_boundary_values_untouched() {
        let n = 8;
        let u0 = grid(n, n, |x, _| if x == 0.0 { 5.0 } else { 0.0 });
        let f = Matrix::zeros(n, n);
        let sol = solve_poisson_2d_sor(&f, &u0, 0.1, 0.1, 1.5, 500, 1e-12);
        for j in 0..n {
            assert_eq!(*sol.get(0, j), 5.0);
        }
        // Harmonic function with a maximum on the boundary: interior values between 0 and 5.
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let v = *sol.get(i, j);
                assert!(v > 0.0 && v < 5.0);
            }
        }
    }

    #[test]
    fn finite_difference_operators_are_exact_on_polynomials() {
        let n = 9;
        let d = 1.0 / (n as f64 - 1.0);
        let p = grid(n, n, |x, y| x * x + 3.0 * y * y + x * y);
        let lap = compute_laplacian(&p, d, d);
        let (gx, gy) = compute_gradient(&p, d, d);
        let u = grid(n, n, |x, _| x);
        let v = grid(n, n, |_, y| y);
        let div = compute_divergence(&u, &v, d, d);
        let mut omega_src = (grid(n, n, |_, y| -y), grid(n, n, |x, _| x));
        let vort = compute_vorticity(&omega_src.0, &omega_src.1, d, d);
        omega_src.0 = Matrix::zeros(1, 1);
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let (x, y) = (i as f64 * d, j as f64 * d);
                assert!((lap.get(i, j) - 8.0).abs() < 1e-9, "laplacian");
                assert!((gx.get(i, j) - (2.0 * x + y)).abs() < 1e-9, "d/dx");
                assert!((gy.get(i, j) - (6.0 * y + x)).abs() < 1e-9, "d/dy");
                assert!((div.get(i, j) - 2.0).abs() < 1e-9, "div");
                assert!(
                    (vort.get(i, j) - 2.0).abs() < 1e-9,
                    "vorticity of solid-body rotation"
                );
            }
        }
        // Boundary entries are left at zero.
        assert_eq!(*lap.get(0, 3), 0.0);
        assert_eq!(*gx.get(n - 1, 3), 0.0);
    }

    #[test]
    fn stream_function_velocities_are_divergence_free_and_match_the_gradient() {
        let n = 21;
        let d = 1.0 / (n as f64 - 1.0);
        let psi = grid(n, n, |x, y| (PI * x).sin() * (PI * y).sin());
        let (u, v) = velocity_from_stream_function(&psi, d, d);
        // u = dpsi/dy, v = -dpsi/dx
        for i in 2..n - 2 {
            for j in 2..n - 2 {
                let (x, y) = (i as f64 * d, j as f64 * d);
                // central differences: error ~ h^2 pi^3 / 6 = 1.3e-2
                assert!((u.get(i, j) - PI * (PI * x).sin() * (PI * y).cos()).abs() < 2e-2);
                assert!((v.get(i, j) + PI * (PI * x).cos() * (PI * y).sin()).abs() < 2e-2);
            }
        }
        let div = compute_divergence(&u, &v, d, d);
        for i in 2..n - 2 {
            for j in 2..n - 2 {
                assert!(
                    div.get(i, j).abs() < 1e-9,
                    "central differences of a stream function commute"
                );
            }
        }
    }

    #[test]
    fn stream_function_solves_the_vorticity_poisson_equation() {
        let n = 21;
        let d = 1.0 / (n as f64 - 1.0);
        let omega = grid(n, n, |x, y| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin());
        let psi = compute_stream_function(&omega, d, d, 5000, 1e-13);
        let amp = 2.0 * PI * PI * d * d / (4.0 * (1.0 - (PI * d).cos()));
        let exact = grid(n, n, |x, y| amp * (PI * x).sin() * (PI * y).sin());
        assert!(interior_max_err(&psi, &exact) < 1e-8);
    }

    #[test]
    fn boundary_condition_helpers() {
        let mut m = grid(5, 6, |x, y| 10.0 * x + y);
        apply_dirichlet_bc(&mut m, -1.0);
        for i in 0..5 {
            assert_eq!((*m.get(i, 0), *m.get(i, 5)), (-1.0, -1.0));
        }
        for j in 0..6 {
            assert_eq!((*m.get(0, j), *m.get(4, j)), (-1.0, -1.0));
        }
        assert!(
            (m.get(2, 2) - (10.0 * 0.5 + 0.4)).abs() < 1e-12,
            "interior untouched"
        );

        let mut f = grid(5, 6, |x, y| 100.0 * x + 7.0 * y + 3.0 * x * y);
        let interior = f.clone();
        apply_neumann_bc(&mut f);
        for j in 1..5 {
            assert_eq!(*f.get(0, j), *interior.get(1, j));
            assert_eq!(*f.get(4, j), *interior.get(3, j));
        }
        for i in 1..4 {
            assert_eq!(*f.get(i, 0), *interior.get(i, 1));
            assert_eq!(*f.get(i, 5), *interior.get(i, 4));
        }
    }

    #[test]
    fn norms_and_maxima() {
        let mut u = Matrix::zeros(4, 4);
        let mut v = Matrix::zeros(4, 4);
        *u.get_mut(1, 1) = 3.0;
        *v.get_mut(1, 1) = 4.0;
        *u.get_mut(2, 2) = -1.0;
        assert_eq!(max_velocity_magnitude(&u, &v), 5.0);
        // l2_norm is the root-mean-square: sqrt(sum / (nx ny))
        assert!((l2_norm(&u) - (10.0f64 / 16.0).sqrt()).abs() < 1e-12);
        assert_eq!(max_abs(&u), 3.0);
        assert_eq!(max_abs(&Matrix::zeros(3, 3)), 0.0);
    }

    #[test]
    fn lid_driven_cavity_develops_a_primary_vortex() {
        let n = 17;
        let (psi, omega) = lid_driven_cavity_simple(n, n, 100.0, 1.0, 40, 1e-3);
        assert_eq!(
            (psi.rows(), psi.cols(), omega.rows(), omega.cols()),
            (n, n, n, n)
        );
        // No penetration through the walls.
        for k in 0..n {
            assert_eq!(
                (
                    *psi.get(k, 0),
                    *psi.get(k, n - 1),
                    *psi.get(0, k),
                    *psi.get(n - 1, k)
                ),
                (0.0, 0.0, 0.0, 0.0)
            );
        }
        // The lid drives the fluid in +x, so the interior stream function is negative (u = psi_y).
        let interior_min = (1..n - 1)
            .flat_map(|i| (1..n - 1).map(move |j| (i, j)))
            .map(|(i, j)| *psi.get(i, j))
            .fold(f64::MAX, f64::min);
        let interior_max = (1..n - 1)
            .flat_map(|i| (1..n - 1).map(move |j| (i, j)))
            .map(|(i, j)| *psi.get(i, j))
            .fold(f64::MIN, f64::max);
        assert!(
            interior_min < -1e-3 && interior_max < 1e-9,
            "psi range [{interior_min}, {interior_max}]"
        );
        // Lid vorticity is strongly negative: about -2 U / dy.
        let d = 1.0 / (n as f64 - 1.0);
        assert!(*omega.get(n / 2, n - 1) < -1.0 / d);
        assert!(psi.data().iter().chain(omega.data()).all(|v| v.is_finite()));
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_laplacian_of_quadratic_is_constant(a in -3.0..3.0f64, b in -3.0..3.0f64, c in -3.0..3.0f64) {
            let n = 7;
            let d = 1.0 / (n as f64 - 1.0);
            let f = grid(n, n, |x, y| a * x * x + b * y * y + c * x * y);
            let lap = compute_laplacian(&f, d, d);
            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    prop_assert!((lap.get(i, j) - 2.0 * (a + b)).abs() < 1e-8);
                }
            }
        }

        #[test]
        fn prop_diffusion_with_stable_r_never_amplifies(
            vals in proptest::collection::vec(0.0..1.0f64, 12..30), r in 0.01..0.5f64,
        ) {
            let dx = 0.1;
            let alpha = 1.0;
            let dt = r * dx * dx / alpha;
            let res = solve_diffusion_1d(&vals, alpha, dx, dt, 30);
            let (lo, hi) = (vals.iter().cloned().fold(f64::MAX, f64::min), vals.iter().cloned().fold(f64::MIN, f64::max));
            for u in &res {
                prop_assert!(u.iter().all(|&v| v >= lo - 1e-12 && v <= hi + 1e-12));
            }
        }

        #[test]
        fn prop_upwind_advection_preserves_the_max_principle(
            vals in proptest::collection::vec(0.0..1.0f64, 12..30), cfl in 0.05..1.0f64, wind in prop::bool::ANY,
        ) {
            let dx = 0.1;
            let c = if wind { 1.0 } else { -1.0 };
            let dt = cfl * dx;
            let mut u0 = vals.clone();
            let n = u0.len();
            // Make the ghost cells consistent with the periodic ring the solver uses.
            u0[0] = u0[n - 2];
            u0[n - 1] = u0[1];
            let res = solve_advection_1d(&u0, c, dx, dt, 25);
            for u in &res {
                prop_assert!(u.iter().all(|&v| (-1e-12..=1.0 + 1e-12).contains(&v)));
            }
        }
    }
}
