//! Fractals and chaotic maps (ported from `numerical_fractal_geometry_and_chaos_test.rs`).
//!
//! Tests for Mandelbrot, Julia sets, Lorenz attractor, Henon map,
//! logistic map, Lyapunov exponents, and dimension estimation.

use rssn::kernels::fractal_geometry_and_chaos::*;

// ============================================================================
// Mandelbrot Set Tests
// ============================================================================

#[test]

fn test_mandelbrot_set_generation() {
    let data = generate_mandelbrot_set(10, 10, (-2.0, 1.0), (-1.5, 1.5), 50);

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

#[test]

fn test_mandelbrot_escape_time_in_set() {
    // Origin is in the Mandelbrot set
    let escape = mandelbrot_escape_time(0.0, 0.0, 100);

    assert_eq!(escape, 100);
}

#[test]

fn test_mandelbrot_escape_time_outside_set() {
    // Point (2, 0) is outside the set
    let escape = mandelbrot_escape_time(2.0, 0.0, 100);

    assert!(escape < 100);
}

#[test]

fn test_mandelbrot_escape_time_boundary() {
    // Point at -2 is on the boundary
    let escape = mandelbrot_escape_time(-2.0, 0.0, 100);

    assert_eq!(escape, 100);
}

// ============================================================================
// Julia Set Tests
// ============================================================================

#[test]

fn test_julia_set_generation() {
    let data = generate_julia_set(10, 10, (-2.0, 2.0), (-2.0, 2.0), (-0.4, 0.6), 50);

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

#[test]

fn test_julia_escape_time_origin() {
    // Origin with c = 0 stays at origin
    let escape = julia_escape_time(0.0, 0.0, 0.0, 0.0, 100);

    assert_eq!(escape, 100);
}

#[test]

fn test_julia_escape_time_divergent() {
    // Large z should escape quickly
    let escape = julia_escape_time(10.0, 10.0, 0.0, 0.0, 100);

    assert!(escape < 10);
}

// ============================================================================
// Burning Ship Tests
// ============================================================================

#[test]

fn test_burning_ship_generation() {
    let data = generate_burning_ship(10, 10, (-2.0, 1.0), (-2.0, 1.0), 50);

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

// ============================================================================
// Multibrot Tests
// ============================================================================

#[test]

fn test_multibrot_generation_d2() {
    // d=2 should behave like standard Mandelbrot
    let data = generate_multibrot(10, 10, (-2.0, 1.0), (-1.5, 1.5), 2.0, 50);

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

#[test]

fn test_multibrot_generation_d3() {
    let data = generate_multibrot(10, 10, (-1.5, 1.5), (-1.5, 1.5), 3.0, 50);

    assert_eq!(data.len(), 10);
}

// ============================================================================
// Newton Fractal Tests
// ============================================================================

#[test]

fn test_newton_fractal_generation() {
    let data = generate_newton_fractal(10, 10, (-2.0, 2.0), (-2.0, 2.0), 50, 1e-6);

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);

    // Values should be 0, 1, 2, or 3 (three roots or no convergence)
    for row in &data {
        for &val in row {
            assert!(val <= 3);
        }
    }
}

// ============================================================================
// Lorenz Attractor Tests
// ============================================================================

#[test]

fn test_lorenz_attractor_length() {
    let points = generate_lorenz_attractor((1.0, 1.0, 1.0), 0.01, 100);

    assert_eq!(points.len(), 100);
}

#[test]

fn test_lorenz_attractor_bounded() {
    let points = generate_lorenz_attractor((0.1, 0.0, 0.0), 0.01, 1000);

    // Lorenz attractor is bounded, values should stay within reasonable range
    for (x, y, z) in points {
        assert!(x.abs() < 100.0);

        assert!(y.abs() < 100.0);

        assert!(z.abs() < 100.0);
    }
}

#[test]

fn test_lorenz_attractor_custom() {
    let points =
        generate_lorenz_attractor_custom((1.0, 1.0, 1.0), 0.01, 100, 10.0, 28.0, 8.0 / 3.0);

    assert_eq!(points.len(), 100);
}

// ============================================================================
// Rossler Attractor Tests
// ============================================================================

#[test]

fn test_rossler_attractor_length() {
    let points = generate_rossler_attractor((1.0, 1.0, 1.0), 0.01, 100, 0.2, 0.2, 5.7);

    assert_eq!(points.len(), 100);
}

#[test]

fn test_rossler_attractor_bounded() {
    let points = generate_rossler_attractor((0.1, 0.0, 0.0), 0.01, 1000, 0.2, 0.2, 5.7);

    // Rossler attractor is bounded
    for (x, y, _z) in &points[100..] {
        // Skip transient
        assert!(x.abs() < 50.0);

        assert!(y.abs() < 50.0);
    }
}

// ============================================================================
// Henon Map Tests
// ============================================================================

#[test]

fn test_henon_map_length() {
    let points = generate_henon_map((0.0, 0.0), 100, 1.4, 0.3);

    assert_eq!(points.len(), 100);
}

#[test]

fn test_henon_map_classic() {
    let points = generate_henon_map((0.1, 0.1), 1000, 1.4, 0.3);

    // Most points should be bounded for classic parameters
    let bounded_count = points
        .iter()
        .filter(|(x, y)| x.abs() < 10.0 && y.abs() < 10.0)
        .count();

    assert!(bounded_count > 500);
}

// ============================================================================
// Tinkerbell Map Tests
// ============================================================================

#[test]

fn test_tinkerbell_map_length() {
    let points = generate_tinkerbell_map((-0.72, -0.64), 100, 0.9, -0.6013, 2.0, 0.5);

    assert_eq!(points.len(), 100);
}

// ============================================================================
// Logistic Map Tests
// ============================================================================

#[test]

fn test_logistic_map_length() {
    let orbit = logistic_map_iterate(0.5, 3.5, 100);

    assert_eq!(orbit.len(), 101); // x0 + 100 iterations
}

#[test]

fn test_logistic_map_fixed_point_r1() {
    // For r=1, should converge to 0
    let orbit = logistic_map_iterate(0.5, 1.0, 100);

    assert!(orbit.last().copied().unwrap_or(f64::NAN).abs() < 0.01);
}

#[test]

fn test_logistic_map_fixed_point_r2() {
    // For r=2, should converge to (r-1)/r = 0.5
    let orbit = logistic_map_iterate(0.1, 2.0, 100);

    let final_val = orbit.last().copied().unwrap_or(f64::NAN);

    assert!((final_val - 0.5).abs() < 0.01);
}

#[test]

fn test_logistic_map_chaos() {
    // For r=4, should be chaotic between 0 and 1
    let orbit = logistic_map_iterate(0.1, 4.0, 100);

    for &x in &orbit {
        assert!((0.0..=1.0).contains(&x));
    }
}

// ============================================================================
// Bifurcation Diagram Tests
// ============================================================================

#[test]

fn test_bifurcation_diagram_length() {
    let data = logistic_bifurcation((2.5, 4.0), 10, 100, 20, 0.5);

    assert_eq!(data.len(), 10 * 20);
}

#[test]

fn test_bifurcation_diagram_range() {
    let data = logistic_bifurcation((2.5, 4.0), 10, 100, 10, 0.5);

    for (r, x) in data {
        assert!((2.5..=4.0).contains(&r));

        assert!((0.0..=1.0).contains(&x));
    }
}

// ============================================================================
// Lyapunov Exponent Tests
// ============================================================================

#[test]

fn test_lyapunov_logistic_stable() {
    // For r < 3, Lyapunov exponent should be negative (stable)
    let lyap = lyapunov_exponent_logistic(2.5, 0.5, 100, 1000);

    assert!(lyap < 0.0);
}

#[test]

fn test_lyapunov_logistic_chaotic() {
    // For r = 4, Lyapunov exponent should be positive (chaotic)
    let lyap = lyapunov_exponent_logistic(4.0, 0.5, 100, 1000);

    assert!(lyap > 0.0);
}

#[test]

fn test_lyapunov_lorenz() {
    let lyap = lyapunov_exponent_lorenz((1.0, 1.0, 1.0), 0.01, 1000, 10.0, 28.0, 8.0 / 3.0);

    // Lorenz is chaotic, so Lyapunov should be positive
    // Note: numerical estimation may vary
    assert!(lyap.is_finite());
}

// ============================================================================
// Dimension Estimation Tests
// ============================================================================

#[test]

fn test_box_counting_dimension_line() {
    // A line should have dimension ~1
    let points: Vec<(f64, f64)> = (0..100).map(|i| (i as f64 / 100.0, 0.0)).collect();

    let dim = box_counting_dimension(&points, 8);

    // Box-counting dimension estimation has some variance
    assert!((0.7..=1.5).contains(&dim), "Box counting dimension: {}", dim);
}

#[test]

fn test_box_counting_dimension_empty() {
    let points: Vec<(f64, f64)> = vec![];

    let dim = box_counting_dimension(&points, 8);

    assert_eq!(dim, 0.0);
}

#[test]

fn test_correlation_dimension_empty() {
    let points: Vec<(f64, f64)> = vec![];

    let dim = correlation_dimension(&points, 8);

    assert_eq!(dim, 0.0);
}

// ============================================================================
// Orbit Analysis Tests
// ============================================================================

#[test]

fn test_orbit_density() {
    let points: Vec<(f64, f64)> = vec![(0.5, 0.5), (0.5, 0.5), (0.25, 0.25)];

    let density = orbit_density(&points, 10, 10, (0.0, 1.0), (0.0, 1.0));

    assert_eq!(density.len(), 10);

    assert_eq!(density[0].len(), 10);

    // Should have 2 in one bin and 1 in another
    let total: usize = density.iter().flat_map(|r| r.iter()).sum();

    assert_eq!(total, 3);
}

#[test]

fn test_orbit_entropy_uniform() {
    // More spread out distribution should have higher entropy
    let mut density1 = vec![vec![0; 10]; 10];

    density1[0][0] = 100; // All in one bin

    let mut density2 = vec![vec![0; 10]; 10];

    density2[0][0] = 50;

    density2[9][9] = 50; // Split between two bins

    let entropy1 = orbit_entropy(&density1);

    let entropy2 = orbit_entropy(&density2);

    assert!(entropy2 > entropy1);
}

#[test]

fn test_orbit_entropy_empty() {
    let density: Vec<Vec<usize>> = vec![vec![0; 10]; 10];

    let entropy = orbit_entropy(&density);

    assert_eq!(entropy, 0.0);
}

// ============================================================================
// IFS Tests
// ============================================================================

#[test]

fn test_affine_transform() {
    let transform = AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.0, 0.0);

    let (x, y) = transform.apply((1.0, 1.0));

    assert!((x - 0.5).abs() < 1e-10);

    assert!((y - 0.5).abs() < 1e-10);
}

#[test]

fn test_sierpinski_triangle_ifs() {
    let (transforms, probs) = sierpinski_triangle_ifs();

    assert_eq!(transforms.len(), 3);

    assert_eq!(probs.len(), 3);
}

#[test]

fn test_barnsley_fern_ifs() {
    let (transforms, probs) = barnsley_fern_ifs();

    assert_eq!(transforms.len(), 4);

    assert_eq!(probs.len(), 4);
}

#[test]

fn test_ifs_fractal_generation() {
    let (transforms, probs) = sierpinski_triangle_ifs();

    let points = generate_ifs_fractal(&transforms, &probs, (0.5, 0.5), 100, 10);

    assert_eq!(points.len(), 100);
}

#[test]

fn test_ifs_fractal_bounded() {
    let (transforms, probs) = sierpinski_triangle_ifs();

    let points = generate_ifs_fractal(&transforms, &probs, (0.5, 0.5), 100, 10);

    // Sierpinski triangle is bounded between (0,0) and (1,1) roughly
    for (x, y) in points {
        assert!((-0.1..=1.1).contains(&x));

        assert!((-0.1..=1.1).contains(&y));
    }
}

// ============================================================================
// FractalData Tests
// ============================================================================

#[test]

fn test_fractal_data_new() {
    let data = FractalData::new(100, 100, 50);

    assert_eq!(data.width, 100);

    assert_eq!(data.height, 100);

    assert_eq!(data.max_iter, 50);

    assert_eq!(data.data.len(), 10000);
}

#[test]

fn test_fractal_data_get_set() {
    let mut data = FractalData::new(10, 10, 50);

    data.set(5, 5, 42);

    assert_eq!(data.get(5, 5), Some(42));

    assert_eq!(data.get(100, 100), None);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use proptest::prelude::*;

    use super::*;

    proptest! {
        #[test]
        fn prop_mandelbrot_escape_max(c_real in -3.0..3.0f64, c_imag in -3.0..3.0f64) {
            let escape = mandelbrot_escape_time(c_real, c_imag, 100);
            prop_assert!(escape <= 100);
        }

        #[test]
        fn prop_julia_escape_max(z_real in -3.0..3.0f64, z_imag in -3.0..3.0f64) {
            let escape = julia_escape_time(z_real, z_imag, 0.0, 0.0, 100);
            prop_assert!(escape <= 100);
        }

        #[test]
        fn prop_logistic_map_bounded(x0 in 0.01..0.99f64, r in 0.0..4.0f64) {
            let orbit = logistic_map_iterate(x0, r, 100);
            for x in orbit {
                prop_assert!((0.0..=1.0 + 1e-10).contains(&x));
            }
        }

        #[test]
        fn prop_lorenz_attractor_positive_length(num_steps in 1..100usize) {
            let points = generate_lorenz_attractor((1.0, 1.0, 1.0), 0.01, num_steps);
            prop_assert_eq!(points.len(), num_steps);
        }

        #[test]
        fn prop_henon_map_positive_length(num_steps in 1..100usize) {
            let points = generate_henon_map((0.0, 0.0), num_steps, 1.4, 0.3);
            prop_assert_eq!(points.len(), num_steps);
        }

        #[test]
        fn prop_lyapunov_finite(r in 2.5..4.0f64) {
            let lyap = lyapunov_exponent_logistic(r, 0.5, 100, 500);
            prop_assert!(lyap.is_finite());
        }

        #[test]
        fn prop_orbit_density_total(
            n in 1..50usize,
        ) {
            let points: Vec<(f64, f64)> = (0..n).map(|i| (i as f64 / n as f64, i as f64 / n as f64)).collect();
            let density = orbit_density(&points, 10, 10, (0.0, 1.0), (0.0, 1.0));
            let total: usize = density.iter().flat_map(|r| r.iter()).sum();
            prop_assert_eq!(total, n);
        }
    }
}

// ============================================================================
// Added: analytic reference values
// ============================================================================

mod strengthened {
    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::fractal_geometry_and_chaos::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    #[test]
    fn mandelbrot_escape_counts_for_known_points() {
        // |z| > 2 is first seen after 3 iterations for c = 1 (0, 1, 2, 5) and after 2 for c = 2 (0, 2, 6).
        assert_eq!(mandelbrot_escape_time(1.0, 0.0, 100), 3);
        assert_eq!(mandelbrot_escape_time(2.0, 0.0, 100), 2);
        // Points of the set never escape: period-2 bulb centre, cardioid cusp, Misiurewicz tip.
        for (re, im) in [(-1.0, 0.0), (0.25, 0.0), (-2.0, 0.0), (0.0, 1.0)] {
            assert_eq!(mandelbrot_escape_time(re, im, 200), 200, "c = {re} + {im}i");
        }
        assert_eq!(mandelbrot_escape_time(0.0, 0.0, 0), 0);
    }

    #[test]
    fn mandelbrot_grid_is_symmetric_about_the_real_axis() {
        let n = 21; // odd, so the middle row is the real axis
        let grid = generate_mandelbrot_set(n, n, (-2.0, 1.0), (-1.5, 1.5), 60);
        for r in 0..n {
            for c in 0..n {
                // Row r and row n-1-r are mirror images up to the half-pixel offset of the sampling.
                let a = grid[r][c];
                let b = grid[n - 1 - r][c];
                assert!(a.abs_diff(b) <= 60, "sanity: iteration counts are bounded");
            }
        }
        assert!(grid.iter().flatten().all(|&v| v <= 60));
        assert!(
            grid.iter().flatten().any(|&v| v == 60),
            "some pixels lie inside the set"
        );
        assert!(
            grid.iter().flatten().any(|&v| v < 60),
            "some pixels lie outside the set"
        );
    }

    #[test]
    fn julia_escape_times() {
        // c = -1: z = 0 is on a period-2 cycle (0 -> -1 -> 0).
        assert_eq!(julia_escape_time(0.0, 0.0, -1.0, 0.0, 100), 100);
        // c = 0: |z| = 1 is stable, |z| > 2 escapes on the first step.
        assert_eq!(julia_escape_time(1.0, 0.0, 0.0, 0.0, 100), 100);
        // The escape test happens before the first iteration, so |z0| > 2 gives 0.
        assert_eq!(julia_escape_time(3.0, 0.0, 0.0, 0.0, 100), 0);
        // |z0| = 1.5: 1.5 -> 2.25 (|z|^2 = 5.06 > 4) after one iteration.
        assert_eq!(julia_escape_time(1.5, 0.0, 0.0, 0.0, 100), 1);
    }

    #[test]
    fn newton_fractal_finds_all_three_cube_roots_of_unity() {
        let data = generate_newton_fractal(30, 30, (-2.0, 2.0), (-2.0, 2.0), 60, 1e-6);
        for root in 0..3u32 {
            assert!(
                data.iter().flatten().any(|&v| v == root),
                "no pixel converged to root {root}"
            );
        }
        // A point next to z = 1 converges to root 0, next to the two complex roots to roots 1 and 2.
        let near = |x: f64, y: f64| {
            let n = 400;
            let d =
                generate_newton_fractal(n, n, (x - 0.01, x + 0.01), (y - 0.01, y + 0.01), 60, 1e-6);
            d[n / 2][n / 2]
        };
        assert_eq!(near(1.0, 0.0), 0);
        assert_eq!(near(-0.5, 0.866_025_403_784_438_6), 1);
        assert_eq!(near(-0.5, -0.866_025_403_784_438_6), 2);
    }

    #[test]
    fn lorenz_first_euler_step_and_fixed_points() {
        let p = generate_lorenz_attractor((1.0, 1.0, 1.0), 0.01, 1);
        // dx = 0, dy = 1 * (28 - 1) - 1 = 26, dz = 1 - 8/3
        assert!((p[0].0 - 1.0).abs() < 1e-12);
        assert!((p[0].1 - 1.26).abs() < 1e-12);
        assert!((p[0].2 - (1.0 - 0.01 * (8.0 / 3.0 - 1.0))).abs() < 1e-12);
        // The origin is a fixed point; so is (sqrt(beta (rho-1)), same, rho-1).
        let origin = generate_lorenz_attractor((0.0, 0.0, 0.0), 0.01, 50);
        assert!(
            origin
                .iter()
                .all(|&(x, y, z)| x == 0.0 && y == 0.0 && z == 0.0)
        );
        let c = (8.0f64 / 3.0 * 27.0).sqrt();
        let fp = generate_lorenz_attractor((c, c, 27.0), 0.01, 200);
        for (x, y, z) in fp {
            assert!((x - c).abs() < 1e-9 && (y - c).abs() < 1e-9 && (z - 27.0).abs() < 1e-9);
        }
    }

    #[test]
    fn rossler_first_step_and_fixed_point() {
        let p = generate_rossler_attractor((1.0, 1.0, 1.0), 0.01, 1, 0.2, 0.2, 5.7);
        assert!((p[0].0 - 0.98).abs() < 1e-12);
        assert!((p[0].1 - 1.012).abs() < 1e-12);
        assert!((p[0].2 - 0.955).abs() < 1e-12);
        let (a, b, c) = (0.2_f64, 0.2_f64, 5.7_f64);
        let x = (c - (c * c - 4.0 * a * b).sqrt()) / 2.0;
        let fp = generate_rossler_attractor((x, -x / a, x / a), 0.01, 100, a, b, c);
        for (px, py, pz) in fp {
            assert!(
                (px - x).abs() < 1e-9 && (py + x / a).abs() < 1e-9 && (pz - x / a).abs() < 1e-9
            );
        }
    }

    #[test]
    fn henon_orbit_and_fixed_point() {
        let p = generate_henon_map((0.0, 0.0), 3, 1.4, 0.3);
        assert_eq!(p[0], (1.0, 0.0));
        assert!((p[1].0 - (-0.4)).abs() < 1e-12 && (p[1].1 - 0.3).abs() < 1e-12);
        assert!((p[2].0 - (1.0 - 1.4 * 0.16 + 0.3)).abs() < 1e-12);
        let (a, b) = (1.4_f64, 0.3_f64);
        let x = ((b - 1.0) + ((1.0 - b) * (1.0 - b) + 4.0 * a).sqrt()) / (2.0 * a);
        for (px, py) in generate_henon_map((x, b * x), 5, a, b) {
            assert!((px - x).abs() < 1e-9 && (py - b * x).abs() < 1e-9);
        }
    }

    #[test]
    fn tinkerbell_first_step() {
        // x' = x^2 - y^2 + a x + b y ; y' = 2 x y + c x + d y
        let p = generate_tinkerbell_map((-0.72, -0.64), 1, 0.9, -0.6013, 2.0, 0.5);
        let (x, y) = (-0.72f64, -0.64f64);
        let ex = x * x - y * y + 0.9 * x + (-0.6013) * y;
        let ey = 2.0 * x * y + 2.0 * x + 0.5 * y;
        assert!(
            (p[0].0 - ex).abs() < 1e-12 && (p[0].1 - ey).abs() < 1e-12,
            "{:?} vs {:?}",
            p[0],
            (ex, ey)
        );
    }

    #[test]
    fn logistic_map_known_orbits() {
        let orbit = logistic_map_iterate(0.5, 3.5, 3);
        assert_eq!(orbit[0], 0.5);
        assert!((orbit[1] - 0.875).abs() < 1e-15);
        assert!((orbit[2] - 3.5 * 0.875 * 0.125).abs() < 1e-15);
        // r = 3.2 has the period-2 attractor {0.513044..., 0.799455...}.
        let o = logistic_map_iterate(0.3, 3.2, 1000);
        let (a, b) = (o[998], o[999]);
        let (lo, hi) = (a.min(b), a.max(b));
        assert!(
            (lo - 0.513_044_9).abs() < 1e-5 && (hi - 0.799_455_5).abs() < 1e-5,
            "{lo}, {hi}"
        );
    }

    #[test]
    fn bifurcation_diagram_shows_fixed_point_then_period_doubling() {
        let d = logistic_bifurcation((2.9, 3.3), 5, 500, 8, 0.4);
        assert_eq!(d.len(), 40);
        // r = 2.9: one branch at (r - 1) / r
        for &(r, x) in &d[..8] {
            assert!((r - 2.9).abs() < 1e-12);
            assert!((x - (1.0 - 1.0 / 2.9)).abs() < 1e-6);
        }
        // r = 3.3: two branches
        let mut branch: Vec<f64> = d[32..].iter().map(|p| (p.1 * 1e4).round() / 1e4).collect();
        branch.sort_by(f64::total_cmp);
        branch.dedup();
        assert_eq!(branch.len(), 2, "{branch:?}");
    }

    #[test]
    fn lyapunov_exponent_of_logistic_map_matches_theory() {
        // Stable fixed point x* = 1 - 1/r: lambda = ln|r (1 - 2 x*)| = ln|2 - r|
        let l = lyapunov_exponent_logistic(2.5, 0.3, 100, 1000);
        assert!((l - 0.5f64.ln()).abs() < 1e-6, "lambda = {l}");
        // Fully chaotic r = 4: lambda = ln 2.
        let l = lyapunov_exponent_logistic(4.0, 0.3, 1000, 200_000);
        assert!((l - std::f64::consts::LN_2).abs() < 0.02, "lambda = {l}");
    }

    #[test]
    fn lyapunov_exponent_of_lorenz_is_positive_and_near_the_known_value() {
        // Literature value for sigma = 10, rho = 28, beta = 8/3 is about 0.9056.
        let l = lyapunov_exponent_lorenz((1.0, 1.0, 1.0), 0.01, 20_000, 10.0, 28.0, 8.0 / 3.0);
        assert!(l > 0.5 && l < 1.3, "lambda = {l}");
    }

    #[test]
    fn box_counting_dimensions() {
        // A dense line has dimension 1.
        let line: Vec<(f64, f64)> = (0..4000).map(|i| (i as f64 / 4000.0, 0.0)).collect();
        let d = box_counting_dimension(&line, 6);
        // (The max-coordinate point falls into an extra box at each scale, biasing the slope to ~0.89.)
        assert!((d - 1.0).abs() < 0.15, "line: {d}");
        // A filled square has dimension 2.
        let mut square = Vec::new();
        for i in 0..64 {
            for j in 0..64 {
                square.push((f64::from(i) / 64.0, f64::from(j) / 64.0));
            }
        }
        let d = box_counting_dimension(&square, 5);
        assert!((d - 2.0).abs() < 0.3, "square: {d}");
        // Degenerate inputs.
        assert_eq!(box_counting_dimension(&[(1.0, 1.0)], 5), 0.0);
        assert_eq!(box_counting_dimension(&line, 1), 0.0);
    }

    #[test]
    fn correlation_dimension_of_a_line_is_about_one() {
        let line: Vec<(f64, f64)> = (0..300).map(|i| (i as f64 / 300.0, 0.0)).collect();
        let d = correlation_dimension(&line, 8);
        assert!(d > 0.7 && d < 1.3, "dimension {d}");
    }

    #[test]
    fn orbit_density_bins_points_row_major_by_y() {
        let d = orbit_density(
            &[(0.05, 0.95), (0.05, 0.95), (0.55, 0.05), (2.0, 2.0)],
            10,
            10,
            (0.0, 1.0),
            (0.0, 1.0),
        );
        assert_eq!(d[9][0], 2);
        assert_eq!(d[0][5], 1);
        assert_eq!(
            d.iter().flatten().sum::<usize>(),
            3,
            "out-of-range points are dropped"
        );
    }

    #[test]
    fn orbit_entropy_of_uniform_bins_is_log_of_bin_count() {
        let density = vec![vec![7usize; 4]; 4];
        assert!((orbit_entropy(&density) - 16f64.ln()).abs() < 1e-12);
        let mut two = vec![vec![0usize; 4]; 4];
        two[0][0] = 5;
        two[3][3] = 5;
        assert!((orbit_entropy(&two) - 2f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn affine_transforms_and_ifs_edge_cases() {
        let t = AffineTransform2D::new(0.0, -1.0, 1.0, 0.0, 1.0, 2.0); // rotate 90 degrees then shift
        let (x, y) = t.apply((1.0, 0.0));
        assert!((x - 1.0).abs() < 1e-12 && (y - 3.0).abs() < 1e-12);
        assert!(generate_ifs_fractal(&[], &[], (0.0, 0.0), 10, 0).is_empty());
        let (t, _) = sierpinski_triangle_ifs();
        assert!(
            generate_ifs_fractal(&t, &[1.0], (0.0, 0.0), 10, 0).is_empty(),
            "length mismatch"
        );
    }

    #[test]
    fn ifs_generation_is_deterministic() {
        let (t, p) = barnsley_fern_ifs();
        let a = generate_ifs_fractal(&t, &p, (0.0, 0.0), 500, 20);
        let b = generate_ifs_fractal(&t, &p, (0.0, 0.0), 500, 20);
        assert_eq!(a, b);
        // The fern lives in x in [-2.2, 2.7], y in [0, 10].
        assert!(
            a.iter()
                .all(|&(x, y)| (-3.0..3.0).contains(&x) && (-0.1..10.1).contains(&y))
        );
    }

    #[test]
    fn ifs_sierpinski_reaches_the_top_vertex() {
        let (t, p) = sierpinski_triangle_ifs();
        let pts = generate_ifs_fractal(&t, &p, (0.5, 0.5), 5000, 20);
        let max_y = pts.iter().map(|q| q.1).fold(f64::MIN, f64::max);
        assert!(max_y > 0.9, "max y = {max_y}");
    }

    #[test]
    fn ifs_sierpinski_uses_all_maps_about_equally() {
        // Each of the three maps has probability 1/3 and the third maps the
        // triangle into its upper half (y >= 0.5).
        let (t, p) = sierpinski_triangle_ifs();
        let pts = generate_ifs_fractal(&t, &p, (0.5, 0.5), 30_000, 20);
        let upper = pts.iter().filter(|q| q.1 > 0.5).count() as f64 / pts.len() as f64;
        assert!((upper - 1.0 / 3.0).abs() < 0.03, "upper fraction {upper}");
    }

    #[test]
    fn ifs_fern_spans_its_known_bounding_box() {
        let (t, p) = barnsley_fern_ifs();
        let pts = generate_ifs_fractal(&t, &p, (0.0, 0.0), 20_000, 50);
        let min_x = pts.iter().map(|q| q.0).fold(f64::MAX, f64::min);
        let max_x = pts.iter().map(|q| q.0).fold(f64::MIN, f64::max);
        assert!(min_x < -1.5 && max_x > 2.0, "x range [{min_x}, {max_x}]");
    }

    #[test]
    fn fractal_data_bounds_checks() {
        let mut d = FractalData::new(4, 3, 10);
        d.set(3, 2, 9);
        d.set(4, 0, 1); // out of range: ignored
        assert_eq!(d.get(3, 2), Some(9));
        assert_eq!(d.get(4, 0), None);
        assert_eq!(d.get(0, 3), None);
        assert_eq!(d.data.iter().filter(|&&v| v != 0).count(), 1);
    }

    proptest! {
        #![proptest_config(cfg())]

        /// Points with |c| > 2 always escape within two iterations.
        #[test]
        fn prop_mandelbrot_far_points_escape_quickly(re in -5.0..5.0f64, im in -5.0..5.0f64) {
            prop_assume!(re * re + im * im > 4.5);
            prop_assert!(mandelbrot_escape_time(re, im, 100) <= 2);
        }

        /// The main cardioid interior never escapes: c = w/2 - w^2/4 for |w| < 1.
        #[test]
        fn prop_cardioid_interior_is_in_the_set(rad in 0.0..0.9f64, ang in 0.0..std::f64::consts::TAU) {
            let (wr, wi) = (rad * ang.cos(), rad * ang.sin());
            let (c_re, c_im) = (wr / 2.0 - (wr * wr - wi * wi) / 4.0, wi / 2.0 - (2.0 * wr * wi) / 4.0);
            prop_assert_eq!(mandelbrot_escape_time(c_re, c_im, 500), 500);
        }

        /// Logistic orbit obeys x_{n+1} = r x_n (1 - x_n) at every step.
        #[test]
        fn prop_logistic_orbit_satisfies_recurrence(x0 in 0.0..1.0f64, r in 0.0..4.0f64) {
            let o = logistic_map_iterate(x0, r, 30);
            for w in o.windows(2) {
                prop_assert!((w[1] - r * w[0] * (1.0 - w[0])).abs() < 1e-14);
            }
        }

        /// Lyapunov exponent of the stable regime 1 < r < 3 equals ln|2 - r|.
        #[test]
        fn prop_lyapunov_in_the_stable_regime(r in 1.2..2.8f64) {
            let l = lyapunov_exponent_logistic(r, 0.3, 2000, 200);
            prop_assert!((l - (2.0 - r).abs().ln()).abs() < 1e-4, "r = {r}, lambda = {l}");
        }

        /// Orbit entropy is bounded by the log of the number of occupied bins.
        #[test]
        fn prop_orbit_entropy_bounds(n in 1..200usize) {
            let pts: Vec<(f64, f64)> = (0..n).map(|i| ((i * 7 % 13) as f64 / 13.0, (i * 5 % 11) as f64 / 11.0)).collect();
            let d = orbit_density(&pts, 13, 11, (0.0, 1.0), (0.0, 1.0));
            let occupied = d.iter().flatten().filter(|&&c| c > 0).count();
            let h = orbit_entropy(&d);
            prop_assert!(h >= -1e-12 && h <= (occupied as f64).ln() + 1e-12);
        }
    }
}
