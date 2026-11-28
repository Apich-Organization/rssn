use num_complex::Complex64;
use rssn::kernels::complex_dynamics::koenigs_coordinate;
use rssn::kernels::complex_dynamics::newton_fractal_basin;
use std::cell::Cell;

#[test]
fn koenigs_retains_first_linearized_estimate_for_tiny_multipliers() {
    let origin = Complex64::new(0.0, 0.0);
    let seed = Complex64::new(0.25, 0.0);
    for magnitude in [1e-200, 1e-300] {
        let lambda = Complex64::new(magnitude, 0.0);
        let map = |z| (z + z * z) * magnitude;
        let estimate = koenigs_coordinate(map, origin, lambda, seed, 8).unwrap();
        assert!((estimate - 0.3125).norm() < 1e-14, "{estimate}");
        assert!(
            koenigs_coordinate(|_| Complex64::new(f64::MAX, 0.0), origin, lambda, seed, 1).is_err()
        );
    }
}

#[test]
fn koenigs_scales_displacement_before_subtracting_opposite_large_points() {
    let extent = f64::MAX;
    let fixed = Complex64::new(-0.55 * extent, 0.0);
    let seed = Complex64::new(0.24 * extent, -0.79 * extent);
    let lambda = Complex64::new(0.7, 0.7);
    let map = |z: Complex64| {
        let center = fixed / extent;
        (center + lambda * (z / extent - center)) * extent
    };
    let next = map(seed);
    assert!(next.is_finite());
    assert!(!(next - fixed).is_finite());
    let estimate = koenigs_coordinate(map, fixed, lambda, seed, 1).unwrap();
    assert!((estimate / extent - Complex64::new(0.79, -0.79)).norm() < 1e-14);
}

#[test]
fn newton_uses_finite_quotients_at_large_scales() {
    let origin = Complex64::new(0.0, 0.0);
    let test = newton_fractal_basin(
        |z| (z - 1.0) * 1e200,
        |_| Complex64::new(1e200, 0.0),
        Complex64::new(2.0, 0.0),
        4,
        1e-12,
    );
    assert!(matches!(test, Some((root,1)) if (root - 1.0).norm() < 1e-14));
    let extent = f64::MAX;
    let target = Complex64::new(-0.9 * extent, 0.9 * extent);
    let derivative = Complex64::new(0.5, 0.5);
    let estimate = newton_fractal_basin(
        |z| derivative * (z - target),
        |_| derivative,
        origin,
        4,
        1e-12,
    )
    .unwrap();
    assert!((estimate.0 / extent - target / extent).norm() < 1e-14);
    let calls = Cell::new(0);
    let outside = newton_fractal_basin(
        |_| {
            calls.set(calls.get() + 1);
            Complex64::new(0.9 * extent, 0.9 * extent)
        },
        |_| derivative,
        origin,
        4,
        1e-12,
    );
    assert!(outside.is_none());
    assert_eq!(calls.get(), 1);
}
