//! Comprehensive tests for physics FDM module.

use proptest::prelude::*;
use rssn::physics::physics_fdm::*;

#[test]
fn test_grid_indexing_1d() {
    let mut grid = FdmGrid::new(Dimensions::D1(10));
    grid[5] = 42.0;
    assert_eq!(grid[5], 42.0);
}

#[test]
fn test_grid_indexing_2d() {
    let mut grid = FdmGrid::new(Dimensions::D2(10, 10));
    grid[(5, 5)] = 42.0;
    assert_eq!(grid[(5, 5)], 42.0);
}

#[test]
fn test_heat_equation_2d_stability() {
    // Small simulation to check stability and convergence
    let grid = solve_heat_equation_2d(20, 20, 0.01, 1.0, 1.0, 0.1, 100, |x, y| {
        if x == 10 && y == 10 {
            100.0
        } else {
            0.0
        }
    });

    // Total energy should be conserved (roughly, for zero boundaries it leaks)
    // Here we just check it doesn't blow up
    for &val in grid.as_slice() {
        assert!(val.is_finite());
        assert!(val >= -1e-10); // Temperature shouldn't go negative in this setup
    }
}

#[test]
fn test_wave_equation_2d_basic() {
    let grid = solve_wave_equation_2d(30, 30, 1.0, 1.0, 1.0, 0.1, 50, |x, y| {
        if x == 15 && y == 15 {
            1.0
        } else {
            0.0
        }
    });

    for &val in grid.as_slice() {
        assert!(val.is_finite());
    }
}

#[test]
fn test_poisson_solver() {
    let width = 20;
    let height = 20;
    let mut source = FdmGrid::new(Dimensions::D2(width, height));
    source[(10, 10)] = 10.0; // Positive source => Concave up => Minimum at source

    let u = solve_poisson_2d(width, height, &source, 1.0, 1.0, 1.5, 1000, 1e-6);

    // Potential should be minimum (most negative) at the negative source
    let min_val = u.as_slice().iter().fold(f64::INFINITY, |a, &b| a.min(b));
    println!("min_val: {}", min_val);
    println!("u[(10, 10)]: {}", u[(10, 10)]);
    println!("min_val + 1e-10: {}", min_val + 1e-10);
    assert!(u[(10, 10)] <= min_val + 1e-10);
}

#[test]
fn test_burgers_1d_shocks() {
    let mut initial_u = vec![0.0; 100];
    for i in 0..50 {
        initial_u[i] = 1.0;
    } // Step function

    let result = solve_burgers_1d(&initial_u, 1.0, 0.1, 0.1, 100);

    // Step should smooth out and move to the right
    assert!(result[40] < 1.0);
    assert!(result[60] > 0.0);
}

// ============================================================================
// Property Tests
// ============================================================================

proptest! {
    #[test]
    fn prop_heat_2d_not_blow_up(
        alpha in 0.001..0.05f64,
        dt in 0.01..0.1f64,
        steps in 1usize..20usize
    ) {
        // Condition for stability: dt <= dx^2 / (4 * alpha)
        // With DX=1, dt <= 1 / (4 * alpha)
        // If alpha = 0.05, dt <= 5.0. Our range 0.01..0.1 is safe.
        let grid = solve_heat_equation_2d(10, 10, alpha, 1.0, 1.0, dt, steps, |_, _| 1.0);
        for &val in grid.as_slice() {
            prop_assert!(val.is_finite());
            prop_assert!(val <= 1.0 + 1e-10); // Heat shouldn't increase beyond initial
        }
    }

    #[test]
    fn prop_advection_diffusion_constant(
        c in -1.0..1.0f64,
        d in 0.01..0.2f64,
        val_in in -10.0..10.0f64
    ) {
        let initial = vec![val_in; 20];
        let res = solve_advection_diffusion_1d(&initial, 1.0, c, d, 0.01, 10);
        for &v in &res {
            prop_assert!((v - val_in).abs() < 1e-10);
        }
    }
}
