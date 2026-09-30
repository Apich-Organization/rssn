//! Finite-difference calculus kernels (ported from the old
//! `numerical_calculus_test.rs` and `numerical/calculus.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::calculus::{derivative, gradient, hessian, jacobian, partial_derivative};

#[test]
fn derivative_of_x_squared() {
    // d/dx x^2 at 3 = 6
    assert_approx_eq!(derivative(|x| x * x, 3.0), 6.0, 1e-5);
}

#[test]
fn derivative_of_transcendentals() {
    assert_approx_eq!(derivative(f64::sin, 0.0), 1.0, 1e-8);
    assert_approx_eq!(derivative(f64::exp, 1.0), std::f64::consts::E, 1e-7);
    assert_approx_eq!(derivative(f64::ln, 2.0), 0.5, 1e-8);
}

#[test]
fn partial_derivative_single_variable() {
    let f = |p: &[f64]| p[0] * p[0];
    assert_approx_eq!(partial_derivative(f, &[3.0], 0), 6.0, 1e-5);
}

#[test]
fn partial_derivative_out_of_range_index_is_nan() {
    assert!(partial_derivative(|p: &[f64]| p[0], &[1.0], 3).is_nan());
}

#[test]
fn gradient_x_squared() {
    let grad = gradient(|p: &[f64]| p[0] * p[0], &[3.0]);
    assert_eq!(grad.len(), 1);
    assert_approx_eq!(grad[0], 6.0, 1e-6);
}

#[test]
fn gradient_x_squared_plus_y_squared() {
    let grad = gradient(|p: &[f64]| p[0] * p[0] + p[1] * p[1], &[1.0, 2.0]);
    assert_eq!(grad.len(), 2);
    assert_approx_eq!(grad[0], 2.0, 1e-6);
    assert_approx_eq!(grad[1], 4.0, 1e-6);
}

#[test]
fn gradient_sin_x_plus_cos_y() {
    let grad = gradient(
        |p: &[f64]| p[0].sin() + p[1].cos(),
        &[0.0, std::f64::consts::FRAC_PI_2],
    );
    assert_approx_eq!(grad[0], 1.0, 1e-6);
    assert_approx_eq!(grad[1], -1.0, 1e-6);
}

#[test]
fn gradient_x_squared_plus_two_y() {
    let grad = gradient(|p: &[f64]| p[0] * p[0] + 2.0 * p[1], &[2.0, 5.0]);
    assert_approx_eq!(grad[0], 4.0, 1e-5);
    assert_approx_eq!(grad[1], 2.0, 1e-5);
}

#[test]
fn jacobian_of_product_and_sum_of_squares() {
    // f1 = x*y, f2 = x^2 + y^2 ; J = [[y, x], [2x, 2y]] = [[2,1],[2,4]] at (1,2)
    let f = |p: &[f64], out: &mut [f64]| {
        out[0] = p[0] * p[1];
        out[1] = p[0] * p[0] + p[1] * p[1];
    };
    let jac = jacobian(f, &[1.0, 2.0], 2);
    let want = [2.0, 1.0, 2.0, 4.0];
    assert_eq!(jac.len(), 4);
    for (got, want) in jac.iter().zip(want) {
        assert_approx_eq!(*got, want, 1e-5);
    }
}

#[test]
fn jacobian_of_rectangular_map() {
    let f = |p: &[f64], out: &mut [f64]| {
        out[0] = p[0] * p[1];
        out[1] = p[0] + p[1];
        out[2] = p[0].exp();
    };
    let jac = jacobian(f, &[1.0, 2.0], 3);
    let want = [2.0, 1.0, 1.0, 1.0, std::f64::consts::E, 0.0];
    assert_eq!(jac.len(), 6);
    for (got, want) in jac.iter().zip(want) {
        assert_approx_eq!(*got, want, 1e-6);
    }
}

#[test]
fn hessian_of_cubic_form() {
    // f = x^2 y + y^3 ; H = [[2y, 2x], [2x, 6y]] = [[4,2],[2,12]] at (1,2)
    let f = |p: &[f64]| p[0] * p[0] * p[1] + p[1] * p[1] * p[1];
    let h = hessian(f, &[1.0, 2.0]);
    let want = [4.0, 2.0, 2.0, 12.0];
    for (got, want) in h.iter().zip(want) {
        assert_approx_eq!(*got, want, 1e-4);
    }
    assert_eq!(h[1].to_bits(), h[2].to_bits(), "Hessian must be exactly symmetric");
}

/// Fixed seed so the property tests are reproducible.
fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_derivative_of_linear_is_slope(a in -10.0..10.0f64, b in -10.0..10.0f64, x in -10.0..10.0f64) {
        prop_assert!((derivative(|t| a * t + b, x) - a).abs() < 1e-5);
    }

    #[test]
    fn prop_derivative_of_quadratic(a in -5.0..5.0f64, x in -5.0..5.0f64) {
        prop_assert!((derivative(|t| a * t * t, x) - 2.0 * a * x).abs() < 1e-4);
    }

    #[test]
    fn prop_gradient_of_linear_form(
        a in -5.0..5.0f64, b in -5.0..5.0f64, x in -5.0..5.0f64, y in -5.0..5.0f64,
    ) {
        let g = gradient(|p: &[f64]| a * p[0] + b * p[1], &[x, y]);
        prop_assert!((g[0] - a).abs() < 1e-5);
        prop_assert!((g[1] - b).abs() < 1e-5);
    }

    #[test]
    fn prop_gradient_matches_jacobian_of_scalar_map(x in -3.0..3.0f64, y in -3.0..3.0f64) {
        let f = |p: &[f64]| p[0] * p[0] * p[1] + p[1].sin();
        let g = gradient(f, &[x, y]);
        let j = jacobian(|p: &[f64], out: &mut [f64]| out[0] = f(p), &[x, y], 1);
        prop_assert!((g[0] - j[0]).abs() < 1e-6);
        prop_assert!((g[1] - j[1]).abs() < 1e-6);
        // Analytic values.
        prop_assert!((g[0] - 2.0 * x * y).abs() < 1e-5);
        prop_assert!((g[1] - (x * x + y.cos())).abs() < 1e-5);
    }

    #[test]
    fn prop_hessian_of_quadratic_form_is_constant(
        a in -3.0..3.0f64, b in -3.0..3.0f64, c in -3.0..3.0f64,
        x in -2.0..2.0f64, y in -2.0..2.0f64,
    ) {
        // f = a x^2 + b x y + c y^2  =>  H = [[2a, b], [b, 2c]]
        let h = hessian(|p: &[f64]| a * p[0] * p[0] + b * p[0] * p[1] + c * p[1] * p[1], &[x, y]);
        prop_assert!((h[0] - 2.0 * a).abs() < 1e-4);
        prop_assert!((h[1] - b).abs() < 1e-4);
        prop_assert!((h[2] - b).abs() < 1e-4);
        prop_assert!((h[3] - 2.0 * c).abs() < 1e-4);
    }
}
