// File: tests/numerical/integrate.rs

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::numerical::integrate::adaptive_quadrature;
use rssn::numerical::integrate::gauss_legendre_quadrature;
use rssn::numerical::integrate::romberg_integration;
use rssn::numerical::integrate::simpson_rule;
use rssn::numerical::integrate::trapezoidal_rule;
use rssn::numerical::integrate::QuadratureMethod;
use rssn::numerical::integrate::{
    self,
};

// --- Helper Functions ---
fn poly_2(
    c : [f64; 3],
    x : f64,
) -> f64 {

    c[0] + c[1] * x + c[2] * x * x
}

fn integral_poly_2(
    c : [f64; 3],
    a : f64,
    b : f64,
) -> f64 {

    let antideriv = |x : f64| c[0] * x + 0.5 * c[1] * x * x + (c[2] * x * x * x) / 3.0;

    antideriv(b) - antideriv(a)
}

// --- 1. Standard Unit Tests ---

#[test]

fn test_trapezoidal_basic() {

    let f = |x : f64| x * x;

    let res = trapezoidal_rule(f, (0.0, 1.0), 1000);

    assert_approx_eq!(res, 1.0 / 3.0, 1e-4);
}

#[test]

fn test_simpson_basic() {

    let f = |x : f64| x * x;

    // Simpson's rule is exact for polynomials up to degree 2 (and actually 3).
    let res = simpson_rule(f, (0.0, 1.0), 10).unwrap();

    assert_approx_eq!(
        res,
        1.0 / 3.0,
        1e-10
    );
}

#[test]

fn test_adaptive_basic() {

    let f = |x : f64| x.sin();

    // Integral of sin(x) from 0 to pi is -cos(pi) - (-cos(0)) = 1 - (-1) = 2
    let res = adaptive_quadrature(
        f,
        (
            0.0,
            std::f64::consts::PI,
        ),
        1e-6,
    );

    assert_approx_eq!(res, 2.0, 1e-6);
}

#[test]

fn test_romberg_basic() {

    let f = |x : f64| x.exp();

    // Integral of exp(x) from 0 to 1 is e^1 - e^0 = e - 1
    let exact = std::f64::consts::E - 1.0;

    let res = romberg_integration(f, (0.0, 1.0), 5);

    assert_approx_eq!(res, exact, 1e-8);
}

#[test]

fn test_gauss_legendre_basic() {

    let f = |x : f64| x.powi(3) + x.powi(2) + 1.0;

    // Exact: [x^4/4 + x^3/3 + x]_0^1 = 1/4 + 1/3 + 1 = 0.25 + 0.333... + 1 = 1.58333...
    let exact = 1.0 / 4.0 + 1.0 / 3.0 + 1.0;

    // Gauss-Legendre with n=5 is exact for polynomials up to degree 2n-1 = 9.
    let res = gauss_legendre_quadrature(f, (0.0, 1.0));

    assert_approx_eq!(res, exact, 1e-10);
}

// --- 2. Property-Based Tests (Proptest) ---

proptest! {
    // Limits
    #![proptest_config(ProptestConfig::with_cases(100))]

    #[test]
    fn prop_test_linearity_trapezoid(
        c1 in -10.0..10.0f64,
        c2 in -10.0..10.0f64,
        a in -5.0..0.0f64,
        b in 0.0..5.0f64
    ) {
        let f1 = |x: f64| x;
        let f2 = |x: f64| x * x;
        let combined = |x: f64| c1 * f1(x) + c2 * f2(x);

        let range = (a, b);
        let n = 100;

        let i_f1 = trapezoidal_rule(f1, range, n);
        let i_f2 = trapezoidal_rule(f2, range, n);
        let i_comb = trapezoidal_rule(combined, range, n);

        // Trapezoidal is a linear operator
        assert_approx_eq!(i_comb, c1 * i_f1 + c2 * i_f2, 1e-9);
    }

    #[test]
    fn prop_test_reversal(
        a in -10.0..10.0f64,
        b in -10.0..10.0f64
    ) {
        let f = |x: f64| x * x - x + 5.0;

        // Simpson
        let i_ab = simpson_rule(f, (a, b), 20).unwrap();
        let i_ba = simpson_rule(f, (b, a), 20).unwrap();

        assert_approx_eq!(i_ab, -i_ba, 1e-9);

        // Adaptive
        let i_ab_adap = adaptive_quadrature(f, (a, b), 1e-6);
        let i_ba_adap = adaptive_quadrature(f, (b, a), 1e-6);
        // Adaptive might vary slightly due to recursion path, but should be close.
        assert_approx_eq!(i_ab_adap, -i_ba_adap, 1e-6);
    }

    #[test]
    fn prop_test_interval_splitting_gauss(
        a in -5.0..0.0f64,
        b in 0.0..5.0f64,
        c in 5.1..10.0f64
    ) {
         let f = |x: f64| x * x + 1.0;

         // Integral(a, c) = Integral(a, b) + Integral(b, c)
         let i_ac = gauss_legendre_quadrature(f, (a, c)); // n=5 ok for x^2
         let i_ab = gauss_legendre_quadrature(f, (a, b));
         let i_bc = gauss_legendre_quadrature(f, (b, c));

         assert_approx_eq!(i_ac, i_ab + i_bc, 1e-9);
    }

    #[test]
    fn prop_test_polynomial_exactness_simpson(
        c0 in -10.0..10.0f64,
        c1 in -10.0..10.0f64,
        c2 in -10.0..10.0f64,
        a in -10.0..0.0f64,
        b in 0.1..10.0f64
    ) {
        let coeffs = [c0, c1, c2];
        let f = |x: f64| poly_2(coeffs, x);
        let exact = integral_poly_2(coeffs, a, b);
        let approx = simpson_rule(f, (a, b), 10).unwrap();

        // Simpson is exact for quadratic
        assert_approx_eq!(approx, exact, 1e-9);
    }
}
