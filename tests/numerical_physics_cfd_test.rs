//! Comprehensive tests for numerical CFD module.
//!
//! Tests for fluid properties, dimensionless numbers, and PDE solvers.

use rssn::numerical::matrix::Matrix;
use rssn::numerical::physics_cfd::*;

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
    use super::*;
    use proptest::prelude::*;

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
