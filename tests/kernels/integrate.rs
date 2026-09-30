//! Numerical quadrature (ported from the old `numerical_integrate_test.rs`).
//! The `Expr`-based `quadrature` front end no longer exists.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::integrate::{
    QuadratureMethod, adaptive_quadrature, gauss_legendre_quadrature, romberg_integration,
    simpson_rule, trapezoidal_rule,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 100,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn poly_2(c: [f64; 3], x: f64) -> f64 {
    c[0] + c[1] * x + c[2] * x * x
}

fn integral_poly_2(c: [f64; 3], a: f64, b: f64) -> f64 {
    let anti = |x: f64| c[0] * x + 0.5 * c[1] * x * x + c[2] * x * x * x / 3.0;
    anti(b) - anti(a)
}

/// Cubic `c0 + c1 x + c2 x^2 + c3 x^3` and its exact integral over `[a, b]`.
fn cubic(c: [f64; 4], x: f64) -> f64 {
    c[0] + x * (c[1] + x * (c[2] + x * c[3]))
}

fn integral_cubic(c: [f64; 4], a: f64, b: f64) -> f64 {
    let anti = |x: f64| x * (c[0] + x * (c[1] / 2.0 + x * (c[2] / 3.0 + x * c[3] / 4.0)));
    anti(b) - anti(a)
}

#[test]
fn trapezoidal_x_squared() {
    let res = trapezoidal_rule(|x| x * x, (0.0, 1.0), 1000);
    assert_approx_eq!(res, 1.0 / 3.0, 1e-6);
}

#[test]
fn trapezoidal_is_exact_for_linear_functions() {
    // int_0^2 (3x + 1) dx = 6 + 2 = 8, one step suffices.
    assert_approx_eq!(trapezoidal_rule(|x| 3.0 * x + 1.0, (0.0, 2.0), 1), 8.0, 1e-12);
}

#[test]
fn simpson_x_squared() {
    let res = simpson_rule(|x| x * x, (0.0, 1.0), 10).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(res, 1.0 / 3.0, 1e-12);
}

#[test]
fn simpson_degenerate_cases() {
    assert_eq!(simpson_rule(|x| x, (0.0, 1.0), 0), Ok(0.0));
    assert_eq!(simpson_rule(|x| x, (2.0, 2.0), 10), Ok(0.0));
    // Odd step counts are rounded up to even and stay accurate.
    let odd = simpson_rule(|x| x * x * x, (0.0, 2.0), 7).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(odd, 4.0, 1e-12);
}

#[test]
fn adaptive_sin_over_half_period() {
    let res = adaptive_quadrature(f64::sin, (0.0, std::f64::consts::PI), 1e-6);
    assert_approx_eq!(res, 2.0, 1e-6);
}

#[test]
fn adaptive_handles_sharp_peak() {
    // int_{-1}^{1} 1/(1 + 100 x^2) dx = (2/10) atan(10)
    let exact = 0.2 * 10f64.atan();
    let res = adaptive_quadrature(|x| 1.0 / (1.0 + 100.0 * x * x), (-1.0, 1.0), 1e-8);
    assert_approx_eq!(res, exact, 1e-6);
}

#[test]
fn romberg_exp() {
    let exact = std::f64::consts::E - 1.0;
    assert_approx_eq!(romberg_integration(f64::exp, (0.0, 1.0), 5), exact, 1e-8);
    assert_approx_eq!(romberg_integration(f64::exp, (0.0, 1.0), 8), exact, 1e-12);
}

#[test]
fn romberg_degenerate_cases() {
    assert_eq!(romberg_integration(f64::exp, (0.0, 1.0), 0), 0.0);
    assert_eq!(romberg_integration(f64::exp, (1.0, 1.0), 5), 0.0);
}

#[test]
fn gauss_legendre_polynomial_of_degree_three() {
    let exact = 1.0 / 4.0 + 1.0 / 3.0 + 1.0;
    let res = gauss_legendre_quadrature(|x| x.powi(3) + x.powi(2) + 1.0, (0.0, 1.0));
    assert_approx_eq!(res, exact, 1e-10);
}

#[test]
fn gauss_legendre_exact_up_to_degree_nine() {
    // int_{-1}^{2} x^9 dx = (2^10 - 1)/10
    let res = gauss_legendre_quadrature(|x| x.powi(9), (-1.0, 2.0));
    assert_approx_eq!(res, (1024.0 - 1.0) / 10.0, 1e-9);
}

#[test]
fn quadrature_method_enum_is_comparable() {
    assert_eq!(QuadratureMethod::Simpson, QuadratureMethod::Simpson);
    assert_ne!(QuadratureMethod::Simpson, QuadratureMethod::Romberg);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_trapezoid_is_linear(
        c1 in -10.0..10.0f64, c2 in -10.0..10.0f64, a in -5.0..0.0f64, b in 0.0..5.0f64,
    ) {
        let combined = |x: f64| c1 * x + c2 * x * x;
        let i1 = trapezoidal_rule(|x| x, (a, b), 100);
        let i2 = trapezoidal_rule(|x| x * x, (a, b), 100);
        let ic = trapezoidal_rule(combined, (a, b), 100);
        prop_assert!((ic - (c1 * i1 + c2 * i2)).abs() < 1e-9);
    }

    #[test]
    fn prop_reversing_limits_negates_integral(a in -10.0..10.0f64, b in -10.0..10.0f64) {
        let f = |x: f64| x * x - x + 5.0;
        let ab = simpson_rule(f, (a, b), 20).map_err(TestCaseError::fail)?;
        let ba = simpson_rule(f, (b, a), 20).map_err(TestCaseError::fail)?;
        prop_assert!((ab + ba).abs() < 1e-9);
        let aab = adaptive_quadrature(f, (a, b), 1e-6);
        let aba = adaptive_quadrature(f, (b, a), 1e-6);
        prop_assert!((aab + aba).abs() < 1e-6);
    }

    #[test]
    fn prop_gauss_legendre_interval_splitting(a in -5.0..0.0f64, b in 0.0..5.0f64, c in 5.1..10.0f64) {
        let f = |x: f64| x * x + 1.0;
        let whole = gauss_legendre_quadrature(f, (a, c));
        let parts = gauss_legendre_quadrature(f, (a, b)) + gauss_legendre_quadrature(f, (b, c));
        prop_assert!((whole - parts).abs() < 1e-9);
    }

    #[test]
    fn prop_simpson_exact_for_quadratics(
        c0 in -10.0..10.0f64, c1 in -10.0..10.0f64, c2 in -10.0..10.0f64,
        a in -10.0..0.0f64, b in 0.1..10.0f64,
    ) {
        let c = [c0, c1, c2];
        let approx = simpson_rule(|x| poly_2(c, x), (a, b), 10).map_err(TestCaseError::fail)?;
        prop_assert!((approx - integral_poly_2(c, a, b)).abs() < 1e-8);
    }

    #[test]
    fn prop_simpson_exact_for_random_cubics(
        c0 in -10.0..10.0f64, c1 in -10.0..10.0f64, c2 in -10.0..10.0f64, c3 in -10.0..10.0f64,
        a in -5.0..0.0f64, b in 0.1..5.0f64, n in 1..20usize,
    ) {
        let c = [c0, c1, c2, c3];
        let approx = simpson_rule(|x| cubic(c, x), (a, b), 2 * n).map_err(TestCaseError::fail)?;
        let exact = integral_cubic(c, a, b);
        prop_assert!((approx - exact).abs() < 1e-9 * exact.abs().max(1.0), "{approx} vs {exact}");
    }

    #[test]
    fn prop_gauss_legendre_exact_for_random_cubics(
        c0 in -10.0..10.0f64, c1 in -10.0..10.0f64, c2 in -10.0..10.0f64, c3 in -10.0..10.0f64,
        a in -5.0..0.0f64, b in 0.1..5.0f64,
    ) {
        let c = [c0, c1, c2, c3];
        let approx = gauss_legendre_quadrature(|x| cubic(c, x), (a, b));
        let exact = integral_cubic(c, a, b);
        prop_assert!((approx - exact).abs() < 1e-9 * exact.abs().max(1.0));
    }

    #[test]
    fn prop_romberg_exact_for_random_cubics(
        c0 in -10.0..10.0f64, c1 in -10.0..10.0f64, c2 in -10.0..10.0f64, c3 in -10.0..10.0f64,
        a in -5.0..0.0f64, b in 0.1..5.0f64,
    ) {
        let c = [c0, c1, c2, c3];
        let approx = romberg_integration(|x| cubic(c, x), (a, b), 4);
        let exact = integral_cubic(c, a, b);
        prop_assert!((approx - exact).abs() < 1e-8 * exact.abs().max(1.0));
    }

    #[test]
    fn prop_adaptive_matches_exponential_antiderivative(k in 0.1..3.0f64, b in 0.5..3.0f64) {
        let approx = adaptive_quadrature(|x| (k * x).exp(), (0.0, b), 1e-9);
        let exact = ((k * b).exp() - 1.0) / k;
        prop_assert!((approx - exact).abs() < 1e-6 * exact.max(1.0));
    }
}
