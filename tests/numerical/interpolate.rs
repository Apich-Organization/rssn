// File: tests/numerical/interpolate.rs

use assert_approx_eq::assert_approx_eq;
use rssn::numerical::interpolate::{cubic_spline_interpolation, lagrange_interpolation};
use rssn::numerical::polynomial::Polynomial;

/// Tests Lagrange interpolation for a simple quadratic function, f(x) = x^2.
/// Given points (0,0), (1,1), and (2,4), the interpolating polynomial should be x^2.
#[test]
fn test_lagrange_interpolation_quadratic() {
    let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
    let poly = lagrange_interpolation(&points).unwrap();

    // The expected polynomial is 1.0 * x^2 + 0.0 * x + 0.0
    // The coefficients are stored from highest degree to lowest.
    let expected_coeffs = vec![1.0, 0.0, 0.0];

    assert_eq!(poly.coeffs.len(), expected_coeffs.len());
    for (c1, c2) in poly.coeffs.iter().zip(expected_coeffs.iter()) {
        assert_approx_eq!(*c1, *c2, 1e-9);
    }
}

/// Tests Lagrange interpolation with a known linear function.
/// Given points (1,2) and (3,4), the interpolating polynomial should be x + 1.
#[test]
fn test_lagrange_interpolation_linear() {
    let points = vec![(1.0, 2.0), (3.0, 4.0)];
    let poly = lagrange_interpolation(&points).unwrap();

    // The expected polynomial is 1.0 * x + 1.0
    let expected_coeffs = vec![1.0, 1.0];

    assert_eq!(poly.coeffs.len(), expected_coeffs.len());
    for (c1, c2) in poly.coeffs.iter().zip(expected_coeffs.iter()) {
        assert_approx_eq!(*c1, *c2, 1e-9);
    }
}

/// Tests cubic spline interpolation with a few points.
/// The test ensures that the spline passes through the given data points.
#[test]
fn test_cubic_spline_interpolation_passes_through_points() {
    let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)];
    let spline = cubic_spline_interpolation(&points).unwrap();

    for (x, y) in &points {
        assert_approx_eq!(spline(*x), *y, 1e-9);
    }
}

/// Tests cubic spline interpolation at an intermediate point.
/// For a simple linear set of points, the spline should behave linearly.
#[test]
fn test_cubic_spline_interpolation_intermediate_point() {
    let points = vec![(0.0, 0.0), (1.0, 2.0), (2.0, 4.0), (3.0, 6.0)];
    let spline = cubic_spline_interpolation(&points).unwrap();

    // Test a point halfway between two data points.
    assert_approx_eq!(spline(1.5), 3.0, 1e-9);
}
