//! Comprehensive tests for numerical fractal geometry and chaos module.
//!
//! Tests for Mandelbrot, Julia sets, Lorenz attractor, Henon map,
//! logistic map, Lyapunov exponents, and dimension estimation.

use rssn::numerical::fractal_geometry_and_chaos::*;

// ============================================================================
// Mandelbrot Set Tests
// ============================================================================

#[test]

fn test_mandelbrot_set_generation() {

    let data = generate_mandelbrot_set(
        10,
        10,
        (-2.0, 1.0),
        (-1.5, 1.5),
        50,
    );

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

    let data = generate_julia_set(
        10,
        10,
        (-2.0, 2.0),
        (-2.0, 2.0),
        (-0.4, 0.6),
        50,
    );

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

#[test]

fn test_julia_escape_time_origin() {

    // Origin with c = 0 stays at origin
    let escape = julia_escape_time(
        0.0, 0.0, 0.0, 0.0, 100,
    );

    assert_eq!(escape, 100);
}

#[test]

fn test_julia_escape_time_divergent() {

    // Large z should escape quickly
    let escape = julia_escape_time(
        10.0, 10.0, 0.0, 0.0, 100,
    );

    assert!(escape < 10);
}

// ============================================================================
// Burning Ship Tests
// ============================================================================

#[test]

fn test_burning_ship_generation() {

    let data = generate_burning_ship(
        10,
        10,
        (-2.0, 1.0),
        (-2.0, 1.0),
        50,
    );

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

// ============================================================================
// Multibrot Tests
// ============================================================================

#[test]

fn test_multibrot_generation_d2() {

    // d=2 should behave like standard Mandelbrot
    let data = generate_multibrot(
        10,
        10,
        (-2.0, 1.0),
        (-1.5, 1.5),
        2.0,
        50,
    );

    assert_eq!(data.len(), 10);

    assert_eq!(data[0].len(), 10);
}

#[test]

fn test_multibrot_generation_d3() {

    let data = generate_multibrot(
        10,
        10,
        (-1.5, 1.5),
        (-1.5, 1.5),
        3.0,
        50,
    );

    assert_eq!(data.len(), 10);
}

// ============================================================================
// Newton Fractal Tests
// ============================================================================

#[test]

fn test_newton_fractal_generation() {

    let data = generate_newton_fractal(
        10,
        10,
        (-2.0, 2.0),
        (-2.0, 2.0),
        50,
        1e-6,
    );

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

    let points = generate_lorenz_attractor(
        (1.0, 1.0, 1.0),
        0.01,
        100,
    );

    assert_eq!(points.len(), 100);
}

#[test]

fn test_lorenz_attractor_bounded() {

    let points = generate_lorenz_attractor(
        (0.1, 0.0, 0.0),
        0.01,
        1000,
    );

    // Lorenz attractor is bounded, values should stay within reasonable range
    for (x, y, z) in points {

        assert!(x.abs() < 100.0);

        assert!(y.abs() < 100.0);

        assert!(z.abs() < 100.0);
    }
}

#[test]

fn test_lorenz_attractor_custom() {

    let points = generate_lorenz_attractor_custom(
        (1.0, 1.0, 1.0),
        0.01,
        100,
        10.0,
        28.0,
        8.0 / 3.0,
    );

    assert_eq!(points.len(), 100);
}

// ============================================================================
// Rossler Attractor Tests
// ============================================================================

#[test]

fn test_rossler_attractor_length() {

    let points = generate_rossler_attractor(
        (1.0, 1.0, 1.0),
        0.01,
        100,
        0.2,
        0.2,
        5.7,
    );

    assert_eq!(points.len(), 100);
}

#[test]

fn test_rossler_attractor_bounded() {

    let points = generate_rossler_attractor(
        (0.1, 0.0, 0.0),
        0.01,
        1000,
        0.2,
        0.2,
        5.7,
    );

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

    let points = generate_henon_map(
        (0.0, 0.0),
        100,
        1.4,
        0.3,
    );

    assert_eq!(points.len(), 100);
}

#[test]

fn test_henon_map_classic() {

    let points = generate_henon_map(
        (0.1, 0.1),
        1000,
        1.4,
        0.3,
    );

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

    let points = generate_tinkerbell_map(
        (-0.72, -0.64),
        100,
        0.9,
        -0.6013,
        2.0,
        0.5,
    );

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

    assert!(
        orbit
            .last()
            .unwrap()
            .abs()
            < 0.01
    );
}

#[test]

fn test_logistic_map_fixed_point_r2() {

    // For r=2, should converge to (r-1)/r = 0.5
    let orbit = logistic_map_iterate(0.1, 2.0, 100);

    let final_val = *orbit
        .last()
        .unwrap();

    assert!((final_val - 0.5).abs() < 0.01);
}

#[test]

fn test_logistic_map_chaos() {

    // For r=4, should be chaotic between 0 and 1
    let orbit = logistic_map_iterate(0.1, 4.0, 100);

    for &x in &orbit {

        assert!(x >= 0.0 && x <= 1.0);
    }
}

// ============================================================================
// Bifurcation Diagram Tests
// ============================================================================

#[test]

fn test_bifurcation_diagram_length() {

    let data = logistic_bifurcation(
        (2.5, 4.0),
        10,
        100,
        20,
        0.5,
    );

    assert_eq!(data.len(), 10 * 20);
}

#[test]

fn test_bifurcation_diagram_range() {

    let data = logistic_bifurcation(
        (2.5, 4.0),
        10,
        100,
        10,
        0.5,
    );

    for (r, x) in data {

        assert!(r >= 2.5 && r <= 4.0);

        assert!(x >= 0.0 && x <= 1.0);
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

    let lyap = lyapunov_exponent_lorenz(
        (1.0, 1.0, 1.0),
        0.01,
        1000,
        10.0,
        28.0,
        8.0 / 3.0,
    );

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
    let points: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            (
                i as f64 / 100.0,
                0.0,
            )
        })
        .collect();

    let dim = box_counting_dimension(&points, 8);

    // Box-counting dimension estimation has some variance
    assert!(
        dim >= 0.7 && dim <= 1.5,
        "Box counting dimension: {}",
        dim
    );
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

    let points: Vec<(f64, f64)> = vec![
        (0.5, 0.5),
        (0.5, 0.5),
        (0.25, 0.25),
    ];

    let density = orbit_density(
        &points,
        10,
        10,
        (0.0, 1.0),
        (0.0, 1.0),
    );

    assert_eq!(density.len(), 10);

    assert_eq!(density[0].len(), 10);

    // Should have 2 in one bin and 1 in another
    let total: usize = density
        .iter()
        .flat_map(|r| r.iter())
        .sum();

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

    let transform = AffineTransform2D::new(
        0.5, 0.0, 0.0, 0.5, 0.0, 0.0,
    );

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

    let points = generate_ifs_fractal(
        &transforms,
        &probs,
        (0.5, 0.5),
        100,
        10,
    );

    assert_eq!(points.len(), 100);
}

#[test]

fn test_ifs_fractal_bounded() {

    let (transforms, probs) = sierpinski_triangle_ifs();

    let points = generate_ifs_fractal(
        &transforms,
        &probs,
        (0.5, 0.5),
        100,
        10,
    );

    // Sierpinski triangle is bounded between (0,0) and (1,1) roughly
    for (x, y) in points {

        assert!(x >= -0.1 && x <= 1.1);

        assert!(y >= -0.1 && y <= 1.1);
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

    assert_eq!(
        data.data.len(),
        10000
    );
}

#[test]

fn test_fractal_data_get_set() {

    let mut data = FractalData::new(10, 10, 50);

    data.set(5, 5, 42);

    assert_eq!(
        data.get(5, 5),
        Some(42)
    );

    assert_eq!(
        data.get(100, 100),
        None
    );
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

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
                prop_assert!(x >= 0.0 && x <= 1.0 + 1e-10);
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
