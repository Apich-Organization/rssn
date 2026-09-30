//! Interpolation and curve evaluation (ported from `numerical_interpolate_test.rs`
//! and `numerical/interpolate.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::interpolate::{
    b_spline, bezier_curve, cubic_spline_interpolation, lagrange_interpolation,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn lagrange_quadratic_value() {
    let poly = lagrange_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)]).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(poly.eval(1.5), 2.25, 1e-9);
}

#[test]
fn lagrange_quadratic_coefficients() {
    // Coefficients are stored highest degree first: x^2.
    let poly = lagrange_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)]).unwrap_or_else(|e| panic!("{e}"));
    let want = [1.0, 0.0, 0.0];
    assert_eq!(poly.coeffs.len(), 3);
    for (got, want) in poly.coeffs.iter().zip(want) {
        assert_approx_eq!(*got, want, 1e-9);
    }
}

#[test]
fn lagrange_linear_coefficients() {
    // Through (1,2) and (3,4): x + 1.
    let poly = lagrange_interpolation(&[(1.0, 2.0), (3.0, 4.0)]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(poly.coeffs.len(), 2);
    assert_approx_eq!(poly.coeffs[0], 1.0, 1e-9);
    assert_approx_eq!(poly.coeffs[1], 1.0, 1e-9);
}

#[test]
fn lagrange_duplicate_abscissae_is_error() {
    assert!(lagrange_interpolation(&[(1.0, 2.0), (1.0, 3.0)]).is_err());
}

#[test]
fn cubic_spline_symmetric_hat() {
    let spline = cubic_spline_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]).unwrap_or_else(|e| panic!("{e}"));
    let val = spline(0.5);
    assert_approx_eq!(val, 0.6875, 1e-9);
    assert_approx_eq!(spline(1.5), 0.6875, 1e-9);
    assert!(val > 0.0 && val < 1.0);
}

#[test]
fn cubic_spline_passes_through_knots() {
    let points = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)];
    let spline = cubic_spline_interpolation(&points).unwrap_or_else(|e| panic!("{e}"));
    for (x, y) in points {
        assert_approx_eq!(spline(x), y, 1e-9);
    }
}

#[test]
fn cubic_spline_reproduces_a_line() {
    let spline = cubic_spline_interpolation(&[(0.0, 0.0), (1.0, 2.0), (2.0, 4.0), (3.0, 6.0)])
        .unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(spline(1.5), 3.0, 1e-9);
    assert_approx_eq!(spline(0.25), 0.5, 1e-9);
}

#[test]
fn cubic_spline_needs_two_points() {
    assert!(cubic_spline_interpolation(&[(0.0, 1.0)]).is_err());
}

#[test]
fn bezier_quadratic_midpoint() {
    let p = bezier_curve(&[vec![0.0, 0.0], vec![1.0, 2.0], vec![2.0, 0.0]], 0.5);
    assert_approx_eq!(p[0], 1.0, 1e-9);
    assert_approx_eq!(p[1], 1.0, 1e-9);
}

#[test]
fn bezier_interpolates_end_points_and_handles_edge_cases() {
    let cp = [vec![1.0, -1.0], vec![3.0, 5.0], vec![-2.0, 4.0], vec![0.0, 0.0]];
    assert_eq!(bezier_curve(&cp, 0.0), vec![1.0, -1.0]);
    assert_eq!(bezier_curve(&cp, 1.0), vec![0.0, 0.0]);
    assert!(bezier_curve(&[], 0.5).is_empty());
    assert_eq!(bezier_curve(&[vec![7.0]], 0.3), vec![7.0]);
}

#[test]
fn b_spline_quadratic_is_a_bezier_curve() {
    let cp = [vec![0.0], vec![1.0], vec![2.0]];
    let knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let p = b_spline(&cp, 2, &knots, 0.5).unwrap_or_default();
    assert_approx_eq!(p[0], 1.0, 1e-9);
    // Same as the Bezier curve with the same control points.
    let cp2 = [vec![0.0, 0.0], vec![1.0, 2.0], vec![2.0, 0.0]];
    let s = b_spline(&cp2, 2, &knots, 0.3).unwrap_or_default();
    let b = bezier_curve(&cp2, 0.3);
    assert_approx_eq!(s[0], b[0], 1e-12);
    assert_approx_eq!(s[1], b[1], 1e-12);
}

#[test]
fn b_spline_rejects_inconsistent_knot_vector() {
    let cp = [vec![0.0], vec![1.0], vec![2.0]];
    assert!(b_spline(&cp, 2, &[0.0, 0.0, 1.0, 1.0], 0.5).is_none());
    assert!(b_spline(&cp, 3, &[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 0.5).is_none());
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_lagrange_reproduces_lines(a in -10.0..10.0f64, b in -10.0..10.0f64, x in -5.0..5.0f64) {
        let f = |x: f64| a * x + b;
        let poly = lagrange_interpolation(&[(0.0, f(0.0)), (1.0, f(1.0))]).map_err(TestCaseError::fail)?;
        prop_assert!((poly.eval(x) - f(x)).abs() < 1e-8);
    }

    #[test]
    fn prop_lagrange_passes_through_data(ys in prop::collection::vec(-10.0..10.0f64, 2..6)) {
        let pts: Vec<(f64, f64)> = ys.iter().enumerate().map(|(i, &y)| (i as f64, y)).collect();
        let poly = lagrange_interpolation(&pts).map_err(TestCaseError::fail)?;
        for (x, y) in pts {
            prop_assert!((poly.eval(x) - y).abs() < 1e-6 * y.abs().max(1.0));
        }
    }

    #[test]
    fn prop_cubic_spline_passes_through_data(ys in prop::collection::vec(-10.0..10.0f64, 3..8)) {
        let pts: Vec<(f64, f64)> = ys.iter().enumerate().map(|(i, &y)| (i as f64, y)).collect();
        let spline = cubic_spline_interpolation(&pts).map_err(TestCaseError::fail)?;
        for (x, y) in pts {
            prop_assert!((spline(x) - y).abs() < 1e-8);
        }
    }

    #[test]
    fn prop_bezier_stays_in_convex_hull_1d(ps in prop::collection::vec(-10.0..10.0f64, 2..7), t in 0.0..1.0f64) {
        let cp: Vec<Vec<f64>> = ps.iter().map(|&p| vec![p]).collect();
        let v = bezier_curve(&cp, t)[0];
        let lo = ps.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = ps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        prop_assert!(v >= lo - 1e-9 && v <= hi + 1e-9);
    }
}
