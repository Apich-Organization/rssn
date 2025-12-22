use rssn::numerical::interpolate::*;
use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;

#[test]
fn test_lagrange_basic() {
    let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
    let poly = lagrange_interpolation(&points).unwrap();
    assert_approx_eq!(poly.eval(1.5), 2.25, 1e-9);
}

#[test]
fn test_cubic_spline_basic() {
    let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)];
    let spline = cubic_spline_interpolation(&points).unwrap();
    assert_approx_eq!(spline(0.5), 0.6875, 1e-9);
    assert_approx_eq!(spline(1.5), 0.6875, 1e-9);
}

#[test]
fn test_bezier_basic() {
    let control_points = vec![vec![0.0, 0.0], vec![1.0, 2.0], vec![2.0, 0.0]];
    let p = bezier_curve(&control_points, 0.5);
    assert_approx_eq!(p[0], 1.0, 1e-9);
    assert_approx_eq!(p[1], 1.0, 1e-9);
}

#[test]
fn test_b_spline_basic() {
    let control_points = vec![vec![0.0], vec![1.0], vec![2.0]];
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let p = b_spline(&control_points, 2, &knots, 0.5);
    assert_approx_eq!(p.unwrap()[0], 1.0, 1e-9);
}

proptest! {
    #[test]
    fn proptest_lagrange_linear(a in -10.0..10.0f64, b in -10.0..10.0f64) {
        let f = |x: f64| a * x + b;
        let points = vec![(0.0, f(0.0)), (1.0, f(1.0))];
        let poly = lagrange_interpolation(&points).unwrap();
        let x_test = 0.5;
        assert_approx_eq!(poly.eval(x_test), f(x_test), 1e-9);
    }
}
