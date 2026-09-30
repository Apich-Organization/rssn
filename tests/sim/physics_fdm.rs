//! Finite-difference PDE solvers on regular grids (ported from `physics_fdm_test.rs`).

use proptest::prelude::*;
use rssn::sim::physics_fdm::*;

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
    let config = FdmSolverConfig2D {
        width: 20,
        height: 20,
        dx: 1.0,
        dy: 1.0,
        dt: 0.1,
        steps: 100,
    };

    let grid = solve_heat_equation_2d(&config, 0.01, |x, y| {
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
    let config = FdmSolverConfig2D {
        width: 30,
        height: 30,
        dx: 1.0,
        dy: 1.0,
        dt: 0.1,
        steps: 50,
    };

    let grid = solve_wave_equation_2d(&config, 1.0, |x, y| {
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

    let config = PoissonSolverConfig2D {
        width,
        height,
        dx: 1.0,
        dy: 1.0,
        omega: 1.5,
        max_iter: 1000,
        tolerance: 1e-6,
    };

    let u = solve_poisson_2d(&config, &source);

    // Potential should be minimum (most negative) at the negative source
    let min_val = u.as_slice().iter().fold(f64::INFINITY, |a, &b| a.min(b));




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
    #![proptest_config(ProptestConfig { rng_seed: proptest::test_runner::RngSeed::Fixed(0x5EED), failure_persistence: None, ..ProptestConfig::default() })]

    #[test]
    fn prop_heat_2d_not_blow_up(
        alpha in 0.001..0.05f64,
        dt in 0.01..0.1f64,
        steps in 1usize..20usize
    ) {
        // Condition for stability: dt <= dx^2 / (4 * alpha)
        // With DX=1, dt <= 1 / (4 * alpha)
        // If alpha = 0.05, dt <= 5.0. Our range 0.01..0.1 is safe.
        let config = FdmSolverConfig2D {
            width: 10,
            height: 10,
            dx: 1.0,
            dy: 1.0,
            dt,
            steps,
        };
        let grid = solve_heat_equation_2d(&config, alpha, |_, _| 1.0);
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

// ============================================================================
// Added: exact discrete solutions and structural checks
// ============================================================================

mod strengthened {
    use std::f64::consts::PI;

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::sim::physics_fdm::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            cases: 32,
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    #[test]
    fn grid_construction_and_row_major_indexing() {
        let g: FdmGrid<f64> = FdmGrid::new(Dimensions::D2(4, 3));
        assert_eq!((g.len(), g.is_empty()), (12, false));
        assert_eq!(g.dimensions(), &Dimensions::D2(4, 3));
        assert!(g.as_slice().iter().all(|&v| v == 0.0));
        let v = FdmGrid::with_value(Dimensions::D3(2, 3, 4), 7.0);
        assert_eq!(v.len(), 24);
        assert!(v.as_slice().iter().all(|&x| x == 7.0));
        let mut d2 = FdmGrid::from_data((0..12).map(f64::from).collect(), Dimensions::D2(4, 3));
        assert_eq!(d2[(1, 2)], 9.0, "index = y * width + x");
        d2[(3, 0)] = -1.0;
        assert_eq!(d2[3], -1.0);
        let d3 = FdmGrid::from_data((0..24).map(f64::from).collect(), Dimensions::D3(2, 3, 4));
        assert_eq!(d3[(1, 2, 3)], 23.0, "index = z * w * h + y * w + x");
        assert_eq!(d3[(0, 1, 2)], 14.0);
        let mut d1: FdmGrid<f64> = FdmGrid::new(Dimensions::D1(0));
        assert!(d1.is_empty());
        d1.as_mut_slice();
    }

    #[test]
    #[should_panic(expected = "2D")]
    fn two_d_indexing_on_a_1d_grid_panics() {
        let g: FdmGrid<f64> = FdmGrid::new(Dimensions::D1(5));
        let _ = g[(1, 1)];
    }

    #[test]
    #[should_panic(expected = "3D")]
    fn three_d_indexing_on_a_2d_grid_panics() {
        let g: FdmGrid<f64> = FdmGrid::new(Dimensions::D2(5, 5));
        let _ = g[(1, 1, 1)];
    }

    #[test]
    fn heat_2d_product_mode_decays_with_the_exact_ftcs_factor() {
        let (w, h, steps) = (21usize, 21usize, 30usize);
        let (alpha, dt) = (0.4, 0.1); // r_x + r_y = 0.08
        let cfg = FdmSolverConfig2D { width: w, height: h, dx: 1.0, dy: 1.0, dt, steps };
        let mode = |x: usize, y: usize| (PI * x as f64 / (w as f64 - 1.0)).sin() * (PI * y as f64 / (h as f64 - 1.0)).sin();
        let grid = solve_heat_equation_2d(&cfg, alpha, mode);
        assert_eq!(grid.dimensions(), &Dimensions::D2(w, h));
        let lam = -2.0 * (1.0 - (PI / (w as f64 - 1.0)).cos());
        let g = 1.0 + alpha * dt * (lam + lam);
        let amp = g.powi(steps as i32);
        for y in 0..h {
            for x in 0..w {
                assert!((grid[(x, y)] - amp * mode(x, y)).abs() < 1e-12, "({x}, {y})");
            }
        }
    }

    #[test]
    fn heat_2d_point_source_spreads_symmetrically_and_conserves_heat_far_from_the_boundary() {
        let cfg = FdmSolverConfig2D { width: 31, height: 31, dx: 1.0, dy: 1.0, dt: 0.2, steps: 15 };
        let grid = solve_heat_equation_2d(&cfg, 1.0, |x, y| if (x, y) == (15, 15) { 1000.0 } else { 0.0 });
        let total: f64 = grid.as_slice().iter().sum();
        assert!((total - 1000.0).abs() < 1e-4, "heat {total}");
        for k in 1..10 {
            assert!((grid[(15 + k, 15)] - grid[(15 - k, 15)]).abs() < 1e-10);
            assert!((grid[(15, 15 + k)] - grid[(15, 15 - k)]).abs() < 1e-10);
            assert!((grid[(15 + k, 15)] - grid[(15, 15 + k)]).abs() < 1e-10, "x/y symmetry");
        }
        // Profile decreases away from the source.
        assert!(grid[(15, 15)] > grid[(16, 15)] && grid[(16, 15)] > grid[(17, 15)]);
        // Second moment grows by 2 r per step and axis: sum x^2 u / sum u = 2 * alpha * dt * steps = 6.
        let m2: f64 = (0..31).map(|x| (0..31).map(|y| grid[(x, y)]).sum::<f64>() * ((x as f64 - 15.0).powi(2))).sum::<f64>() / total;
        assert!((m2 - 2.0 * 1.0 * 0.2 * 15.0).abs() < 0.01, "variance {m2}");
    }

    #[test]
    fn heat_2d_keeps_boundary_values_fixed() {
        let cfg = FdmSolverConfig2D { width: 12, height: 9, dx: 1.0, dy: 1.0, dt: 0.1, steps: 25 };
        let grid = solve_heat_equation_2d(&cfg, 0.5, |x, y| if x == 0 { 10.0 } else if x == 11 || y == 0 || y == 8 { 0.0 } else { 3.0 });
        for y in 0..9 {
            assert_eq!((grid[(0, y)], grid[(11, y)]), (10.0, 0.0));
        }
        // Interior relaxes toward the boundary data (values between the extremes).
        assert!(grid.as_slice().iter().all(|&v| (0.0..=10.0).contains(&v)));
    }

    #[test]
    fn heat_scenario_runs_and_cools_the_hot_disc() {
        let g = simulate_2d_heat_conduction_scenario();
        assert_eq!(g.dimensions(), &Dimensions::D2(100, 100));
        let peak = g.as_slice().iter().cloned().fold(f64::MIN, f64::max);
        assert!(peak < 100.0 && peak > 0.0);
        assert!(g.as_slice().iter().all(|v| v.is_finite() && *v >= -1e-9));
    }

    #[test]
    fn wave_2d_is_bounded_and_symmetric_for_a_symmetric_pulse() {
        let cfg = FdmSolverConfig2D { width: 41, height: 41, dx: 1.0, dy: 1.0, dt: 0.5, steps: 30 };
        let g = solve_wave_equation_2d(&cfg, 1.0, |x, y| (-((x as f64 - 20.0).powi(2) + (y as f64 - 20.0).powi(2)) / 10.0).exp());
        assert!(g.as_slice().iter().all(|v| v.abs() <= 1.05));
        for k in 1..15 {
            assert!((g[(20 + k, 20)] - g[(20 - k, 20)]).abs() < 1e-12);
            assert!((g[(20, 20 + k)] - g[(20 + k, 20)]).abs() < 1e-12);
        }
        // Boundary rows stay zero.
        for k in 0..41 {
            assert_eq!((g[(0, k)], g[(40, k)], g[(k, 0)], g[(k, 40)]), (0.0, 0.0, 0.0, 0.0));
        }
    }

    /// Discrete dispersion of the standing mode sin(pi x / (w-1)) sin(pi y / (h-1)): cos(theta) = 1 + s lambda / 2.
    fn theta(w: usize, dt: f64, c: f64) -> f64 {
        let lam = -2.0 * (1.0 - (PI / (w as f64 - 1.0)).cos());
        (1.0 + (c * dt).powi(2) * (lam + lam) / 2.0).acos()
    }

    #[test]
    #[ignore = "library bug: solve_wave_equation_2d starts leapfrog with u^{-1} = u^0 (\"u_t = 0\") instead of u^1 = u^0 + s/2 Lap u^0, so the solution is advanced by half a step: for the standing mode (w = 41, c = 1, dt = 0.5, 40 steps) the centre value is cos(40.5 theta)/cos(theta/2) = -0.0184 instead of cos(40 theta) = 0.0012"]
    fn wave_2d_standing_mode_follows_the_discrete_dispersion_relation() {
        let w = 41usize;
        let (c, dt, steps) = (1.0, 0.5, 40usize);
        let cfg = FdmSolverConfig2D { width: w, height: w, dx: 1.0, dy: 1.0, dt, steps };
        let mode = |x: usize, y: usize| (PI * x as f64 / (w as f64 - 1.0)).sin() * (PI * y as f64 / (w as f64 - 1.0)).sin();
        let g = solve_wave_equation_2d(&cfg, c, mode);
        let th = theta(w, dt, c);
        let expected = (steps as f64 * th).cos();
        assert!((g[(20, 20)] - expected).abs() < 1e-6, "centre {} vs {expected}", g[(20, 20)]);
    }

    #[test]
    fn wave_2d_standing_mode_oscillates_with_the_right_period_to_first_order() {
        // The half-step start-up error above is O(dt): the solution still tracks cos(omega t) within ~2 theta.
        let w = 41usize;
        let (c, dt) = (1.0, 0.5);
        let mode = |x: usize, y: usize| (PI * x as f64 / (w as f64 - 1.0)).sin() * (PI * y as f64 / (w as f64 - 1.0)).sin();
        let th = theta(w, dt, c);
        for steps in [10usize, 40, 80] {
            let cfg = FdmSolverConfig2D { width: w, height: w, dx: 1.0, dy: 1.0, dt, steps };
            let g = solve_wave_equation_2d(&cfg, c, mode);
            assert!((g[(20, 20)] - (steps as f64 * th).cos()).abs() < 2.0 * th, "steps {steps}: {}", g[(20, 20)]);
        }
    }

    #[test]
    fn wave_3d_symmetric_pulse_stays_symmetric_and_bounded() {
        let n = 15;
        let cfg = FdmSolverConfig3D { width: n, height: n, depth: n, dx: 1.0, dy: 1.0, dz: 1.0, dt: 0.4, steps: 8 };
        let g = solve_wave_equation_3d(&cfg, 1.0, |x, y, z| if (x, y, z) == (7, 7, 7) { 1.0 } else { 0.0 });
        assert_eq!(g.dimensions(), &Dimensions::D3(n, n, n));
        assert!(g.as_slice().iter().all(|v| v.is_finite() && v.abs() <= 1.0 + 1e-12));
        for k in 1..6 {
            let a = g[(7 + k, 7, 7)];
            assert!((a - g[(7 - k, 7, 7)]).abs() < 1e-12);
            assert!((a - g[(7, 7 + k, 7)]).abs() < 1e-12);
            assert!((a - g[(7, 7, 7 + k)]).abs() < 1e-12);
        }
        // The disturbance travels at most one cell per step: nothing beyond radius 8 has moved.
        assert_eq!(g[(0, 0, 0)], 0.0);
    }

    #[test]
    fn poisson_2d_matches_the_exact_discrete_solution() {
        // laplace(u) = -2 pi^2 sin(pi x) sin(pi y) on the unit square with u = 0 on the boundary.
        let n = 21usize;
        let d = 1.0 / (n as f64 - 1.0);
        let mut src = FdmGrid::new(Dimensions::D2(n, n));
        let mode = |x: usize, y: usize| (PI * x as f64 * d).sin() * (PI * y as f64 * d).sin();
        for y in 0..n {
            for x in 0..n {
                src[(x, y)] = -2.0 * PI * PI * mode(x, y);
            }
        }
        let cfg = PoissonSolverConfig2D { width: n, height: n, dx: d, dy: d, omega: 1.7, max_iter: 3000, tolerance: 1e-13 };
        let u = solve_poisson_2d(&cfg, &src);
        let amp = 2.0 * PI * PI * d * d / (4.0 * (1.0 - (PI * d).cos()));
        for y in 0..n {
            for x in 0..n {
                assert!((u[(x, y)] - amp * mode(x, y)).abs() < 1e-8, "({x}, {y}): {}", u[(x, y)]);
            }
        }
    }

    #[test]
    fn poisson_2d_point_source_is_symmetric_and_convex() {
        let n = 21usize;
        let mut src = FdmGrid::new(Dimensions::D2(n, n));
        src[(10, 10)] = 5.0;
        let cfg = PoissonSolverConfig2D { width: n, height: n, dx: 1.0, dy: 1.0, omega: 1.6, max_iter: 2000, tolerance: 1e-12 };
        let u = solve_poisson_2d(&cfg, &src);
        for k in 1..10 {
            assert!((u[(10 + k, 10)] - u[(10 - k, 10)]).abs() < 1e-9);
            assert!((u[(10, 10 + k)] - u[(10 + k, 10)]).abs() < 1e-9);
        }
        // Positive source => minimum at the source, boundary stays zero, all values <= 0.
        assert!(u.as_slice().iter().all(|&v| v <= 1e-12));
        assert!(u[(10, 10)] < u[(11, 10)] && u[(11, 10)] < u[(12, 10)]);
        // Discrete equation at the source: 4 u0 - sum of neighbours = -5.
        let lap = u[(11, 10)] + u[(9, 10)] + u[(10, 11)] + u[(10, 9)] - 4.0 * u[(10, 10)];
        assert!((lap - 5.0).abs() < 1e-8);
    }

    #[test]
    fn burgers_1d_basic_properties() {
        // Constants and zero stay constant.
        assert!(solve_burgers_1d(&vec![2.5; 30], 0.1, 0.05, 0.01, 50).iter().all(|&v| (v - 2.5).abs() < 1e-12));
        assert!(solve_burgers_1d(&vec![0.0; 30], 0.1, 0.05, 0.01, 50).iter().all(|&v| v == 0.0));
        // Very short inputs are returned unchanged.
        assert_eq!(solve_burgers_1d(&[1.0], 0.1, 0.05, 0.01, 5), vec![1.0]);
        // Fixed end values.
        let mut u0 = vec![1.0; 40];
        u0[20..].iter_mut().for_each(|v| *v = 0.0);
        let r = solve_burgers_1d(&u0, 1.0, 0.1, 0.1, 100);
        assert_eq!((r[0], r[39]), (1.0, 0.0));
        assert!(r[15] < 1.0 + 1e-9 && r[30] > -1e-9);
    }

    #[test]
    fn advection_diffusion_reduces_to_translation_and_pure_diffusion() {
        // c dt / dx = 1, d = 0: upwind translates the profile exactly one cell per step.
        let mut u0 = vec![0.0; 30];
        u0[5] = 1.0;
        u0[6] = 2.0;
        let r = solve_advection_diffusion_1d(&u0, 0.1, 1.0, 0.0, 0.1, 7);
        assert_eq!((r[12], r[13]), (1.0, 2.0));
        // Right-to-left wind translates the other way.
        let l = solve_advection_diffusion_1d(&u0, 0.1, -1.0, 0.0, 0.1, 3);
        assert_eq!((l[2], l[3]), (1.0, 2.0));
        // c = 0: FTCS diffusion of a sine mode with the exact amplification factor.
        let n = 41usize;
        let dx = 1.0 / (n as f64 - 1.0);
        let s: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let (d, dt, steps) = (0.1, 0.001, 20usize);
        let res = solve_advection_diffusion_1d(&s, dx, 0.0, d, dt, steps);
        let g = 1.0 - 4.0 * d * dt / (dx * dx) * (PI * dx / 2.0).sin().powi(2);
        for i in 1..n - 1 {
            assert!((res[i] - g.powi(steps as i32) * s[i]).abs() < 1e-12);
        }
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_heat_2d_is_linear_in_the_initial_data(a in -3.0..3.0f64, b in -3.0..3.0f64) {
            let cfg = FdmSolverConfig2D { width: 10, height: 10, dx: 1.0, dy: 1.0, dt: 0.1, steps: 5 };
            let f = |x: usize, y: usize| ((x * 3 + y * 7) % 5) as f64;
            let g = |x: usize, y: usize| ((x * 5 + y) % 4) as f64;
            let lhs = solve_heat_equation_2d(&cfg, 0.2, |x, y| a * f(x, y) + b * g(x, y));
            let fa = solve_heat_equation_2d(&cfg, 0.2, f);
            let gb = solve_heat_equation_2d(&cfg, 0.2, g);
            for i in 0..100 {
                prop_assert!((lhs[i] - (a * fa[i] + b * gb[i])).abs() < 1e-9);
            }
        }

        #[test]
        fn prop_upwind_advection_diffusion_obeys_the_maximum_principle(
            vals in proptest::collection::vec(0.0..1.0f64, 10..30), c in -1.0..1.0f64, d in 0.0..0.2f64,
        ) {
            // dt chosen so that |c| dt / dx + 2 d dt / dx^2 <= 1.
            let dx = 1.0;
            let dt = 0.9 / (c.abs() / dx + 2.0 * d / (dx * dx) + 1e-12);
            let dt = dt.min(0.4);
            let res = solve_advection_diffusion_1d(&vals, dx, c, d, dt, 20);
            prop_assert!(res.iter().all(|&v| v >= -1e-12 && v <= 1.0 + 1e-12));
        }
    }
}
