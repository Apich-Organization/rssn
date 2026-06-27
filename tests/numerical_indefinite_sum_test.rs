// File: tests/numerical_indefinite_sum_test.rs
//
// Integration tests for the 'rssn' crate's public API in the numerical::indefinite_sum module.

use assert_approx_eq::assert_approx_eq;
use rssn::numerical::indefinite_sum::*;
use rssn::numerical::special::hurwitz_zeta;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]
fn test_exp_closed_form() {
    let f = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let anti = try_closed_form_sum(&f, "x");
    assert!(
        anti.is_some(),
        "e^x should have a closed-form anti-difference"
    );
}

#[test]
fn test_constant_closed_form() {
    let f = Expr::new_constant(3.0);
    let anti = try_closed_form_sum(&f, "x").unwrap();
    // F(x) = 3x; F(2) - F(1) = 3*2 - 3*1 = 3 = f(1). Check!
    let f2 = eval_antidiff(&anti, "x", 2.0).unwrap();
    let f1 = eval_antidiff(&anti, "x", 1.0).unwrap();
    let diff = f2 - f1;
    assert_approx_eq!(diff, 3.0, 1e-10f64);
}

#[test]
fn test_sin_recurrence() {
    let a = 0.7_f64;
    let f = Expr::Sin(Arc::new(Expr::new_mul(
        Expr::new_constant(a),
        Expr::new_variable("x"),
    )));
    let anti = try_closed_form_sum(&f, "x").unwrap();
    // F(x+1) - F(x) should equal sin(ax)
    let x = 3.5_f64;
    let fxp1 = eval_antidiff(&anti, "x", x + 1.0).unwrap();
    let fx = eval_antidiff(&anti, "x", x).unwrap();
    let expected = (a * x).sin();
    assert_approx_eq!(fxp1 - fx, expected, 1e-10f64);
}

#[test]
fn test_series_antidiff_constant() {
    // f(z) = 1 (coeffs = [1.0])  →  Δ⁻¹ 1 = z
    let result = series_antidiff(&[1.0], 0.0, 5.0);
    assert_approx_eq!(result, 5.0, 1e-8f64);
}

#[test]
fn test_hurwitz_x_squared_recurrence() {
    // Δ⁻¹ x² = −ζ(−2, x+1) should satisfy F(x+1)−F(x) = (x+1)²
    let x = 4.0_f64;
    let fx = -hurwitz_zeta(-2.0, x + 1.0);
    let fxp1 = -hurwitz_zeta(-2.0, x + 2.0);
    let expected = (x + 1.0).powi(2);
    assert_approx_eq!(fxp1 - fx, expected, 1e-8f64);
}
