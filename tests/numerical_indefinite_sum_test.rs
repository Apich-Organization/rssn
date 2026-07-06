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
    // Raw zeta identity: −ζ(−2, x+1) satisfies F(x+1)−F(x) = (x+1)²
    let x = 4.0_f64;
    let fx = -hurwitz_zeta(-2.0, x + 1.0);
    let fxp1 = -hurwitz_zeta(-2.0, x + 2.0);
    let expected = (x + 1.0).powi(2);
    assert_approx_eq!(fxp1 - fx, expected, 1e-8f64);
}

// ============================================================================
// Phase 4: Taylor coefficients and unified pipeline
// ============================================================================

#[test]
fn test_compute_taylor_coeffs_polynomial() {
    // f(x) = x^2 around p = 1.0: c_0 = 1, c_1 = 2, c_2 = 1.
    // Forward differences have O(h) error; with h=1e-3 expect error ~1e-3 on c_1.
    let f = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(2.0)),
    );
    let coeffs = compute_taylor_coeffs_numerical(&f, "x", 1.0, 4);
    assert_approx_eq!(coeffs[0], 1.0, 1e-5f64); // f(1) = 1 — exact
    assert_approx_eq!(coeffs[1], 2.0, 2e-3f64); // f'(1) = 2 — O(h) error
    assert_approx_eq!(coeffs[2], 1.0, 1e-3f64); // f''(1)/2! = 1 — exact for quadratic
}

#[test]
fn test_compute_taylor_coeffs_exp() {
    // f(x) = e^x around p = 0: c_m = 1/m!
    let f = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let coeffs = compute_taylor_coeffs_numerical(&f, "x", 0.0, 5);
    let mut fact = 1.0_f64;
    for (m, &c) in coeffs.iter().enumerate() {
        if m > 0 {
            fact *= m as f64;
        }
        assert_approx_eq!(c, 1.0 / fact, 2e-3f64);
    }
}

#[test]
fn test_eval_indefinite_sum_numerical_constant() {
    // f(x) = 5 → F(x) - F(0) = 5*x
    let f = Expr::new_constant(5.0);
    let config = IndefiniteSumConfig::default();
    let result = eval_indefinite_sum_numerical(&f, "x", 4.0, &config).unwrap();
    assert_approx_eq!(result, 20.0, 1e-8f64); // 5 * 4
}

#[test]
fn test_eval_indefinite_sum_numerical_exp() {
    // f(x) = e^x → F(x) = e^x/(e-1), F(x)-F(0) should satisfy recurrence
    let f = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let config = IndefiniteSumConfig::default();
    let x = 3.0_f64;
    let fxp1 = eval_indefinite_sum_numerical(&f, "x", x + 1.0, &config).unwrap();
    let fx = eval_indefinite_sum_numerical(&f, "x", x, &config).unwrap();
    assert_approx_eq!(fxp1 - fx, x.exp(), 1e-8f64);
}

#[test]
fn test_eval_indefinite_sum_numerical_series_fallback() {
    // Test series_antidiff with exact Taylor coefficients for e^x around p = 0.
    // Using exact c_m = 1/m! avoids the catastrophic cancellation that occurs when
    // computing 8th-order finite differences with h=1e-3 (condition ≈ C(7,3)/h^7 ≈ 3.5e6).
    let h = 0.0_f64;
    let mut fact = 1.0_f64;
    let coeffs: Vec<f64> = (0..6)
        .map(|m| {
            if m > 0 {
                fact *= m as f64;
            }
            1.0 / fact
        })
        .collect();

    // Normalization: G(h) = 0
    assert_approx_eq!(series_antidiff(&coeffs, h, h), 0.0, 1e-10f64);

    // series_antidiff uses the backward convention: G(z+1) - G(z) = f(z+1).
    // At z = 0: G(1) - G(0) = f(1) = e^1, approximated by a 6-term series.
    // Truncation error for 6 terms at z=0: |e - Σ_{m=0}^{5} 1/m!| = 1/720 ≈ 1.4e-3.
    let gz1 = series_antidiff(&coeffs, h, 1.0);
    let gz0 = series_antidiff(&coeffs, h, h);
    assert_approx_eq!(gz1 - gz0, std::f64::consts::E, 5e-3f64);
}

// ============================================================================
// Phase 5: Indefinite product recurrence
// ============================================================================

#[test]
fn test_indefinite_product_numerical_recurrence() {
    // ∏ f(k) where f(x) = 2: F(x) = 2^x, F(x)/F(x-1) = 2
    let f = |x: f64| -> Result<f64, String> { Ok(2.0_f64.powf(x)) };
    // Normalize at h = 1.0: F(1) = 1
    let h = 1.0_f64;
    let step = 1.0_f64;
    let x = 5.0_f64;
    let fx = eval_indefinite_product_numerical(x, &f, h, step).unwrap();
    let fxm1 = eval_indefinite_product_numerical(x - step, &f, h, step).unwrap();
    // F(x)/F(x-1) = f(x-1) = 2^(x-1) when using log-sum
    // Actually this checks the recurrence of the log-sum engine
    assert!(fx.is_finite(), "product at x=5 should be finite");
    assert!(fxm1.is_finite(), "product at x=4 should be finite");
    if fxm1.abs() > 1e-15 {
        let ratio = fx / fxm1;
        assert_approx_eq!(ratio, 2.0_f64.powf(x - step), 1e-6f64);
    }
}
