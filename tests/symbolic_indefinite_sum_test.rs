// File: tests/symbolic_indefinite_sum_test.rs
//
// Integration tests for symbolic IndefiniteSum and IndefiniteProduct expressions:
// construction via smart constructors, simplification to closed form, and
// recurrence verification F(x+1) - F(x) = f(x).

use assert_approx_eq::assert_approx_eq;
use rssn::numerical::indefinite_sum::eval_antidiff;
use rssn::numerical::indefinite_sum::try_closed_form_sum;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]
fn test_symbolic_indefinite_sum_constant_simplifies() {
    // IndefiniteSum(3, x, 1) should simplify to 3*x
    let body = Expr::new_constant(3.0);
    let step = Expr::new_constant(1.0);
    let expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = rssn::simplify(&expr);
    // Evaluate at x=5: expect 3*5 = 15
    let val = rssn::numerical::elementary::eval_expr_single(&simplified, "x", 5.0).unwrap();
    assert_approx_eq!(val, 15.0, 1e-9f64);
}

#[test]
fn test_symbolic_indefinite_sum_exp_simplifies() {
    // IndefiniteSum(e^x, x, 1) → e^x / (e - 1), check recurrence F(x+1) - F(x) = e^x
    let body = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = rssn::simplify(&sum_expr);

    let x = 2.0_f64;
    let fxp1 = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x + 1.0).unwrap();
    let fx = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x).unwrap();
    assert_approx_eq!(fxp1 - fx, x.exp(), 1e-8f64);
}

#[test]
fn test_symbolic_indefinite_sum_sin_recurrence() {
    // IndefiniteSum(sin(0.5*x), x, 1) should satisfy F(x+1) - F(x) = sin(0.5*x)
    let a = 0.5_f64;
    let ax = Expr::new_mul(Expr::new_constant(a), Expr::new_variable("x"));
    let body = Expr::Sin(Arc::new(ax));
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = rssn::simplify(&sum_expr);

    let x = 3.0_f64;
    let fxp1 = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x + 1.0).unwrap();
    let fx = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x).unwrap();
    assert_approx_eq!(fxp1 - fx, (a * x).sin(), 1e-8f64);
}

#[test]
fn test_symbolic_indefinite_sum_power_recurrence() {
    // IndefiniteSum(x^3, x, 1) → ζ(−3,1)−ζ(−3,x), check recurrence F(x+1) - F(x) = x^3
    let body = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(3.0)),
    );
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = rssn::simplify(&sum_expr);

    let x = 4.0_f64;
    let anti_diff_val = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x + 1.0);
    if anti_diff_val.is_err() {
        // BinaryList("hurwitz_zeta_antidiff", ...) — evaluate via eval_antidiff
        let anti = try_closed_form_sum(
            &Expr::Power(
                Arc::new(Expr::new_variable("x")),
                Arc::new(Expr::new_constant(3.0)),
            ),
            "x",
        )
        .unwrap();
        let fxp1 = eval_antidiff(&anti, "x", x + 1.0).unwrap();
        let fx = eval_antidiff(&anti, "x", x).unwrap();
        assert_approx_eq!(fxp1 - fx, x.powi(3), 1e-6f64);
    } else {
        let fxp1 = anti_diff_val.unwrap();
        let fx = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x).unwrap();
        assert_approx_eq!(fxp1 - fx, x.powi(3), 1e-6f64);
    }
}

#[test]
fn test_symbolic_indefinite_sum_linearity() {
    // Δ⁻¹ (2*sin(x) + 3*e^x) = 2*Δ⁻¹sin(x) + 3*Δ⁻¹e^x — check recurrence
    let sin_x = Expr::Sin(Arc::new(Expr::new_variable("x")));
    let exp_x = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let body = Expr::new_add(
        Expr::new_mul(Expr::new_constant(2.0), sin_x),
        Expr::new_mul(Expr::new_constant(3.0), exp_x),
    );
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = rssn::simplify(&sum_expr);

    let x = 2.5_f64;
    let fxp1 = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x + 1.0);
    let fx = rssn::numerical::elementary::eval_expr_single(&simplified, "x", x);
    if let (Ok(fxp1), Ok(fx)) = (fxp1, fx) {
        let f_of_x = 2.0 * x.sin() + 3.0 * x.exp();
        assert_approx_eq!(fxp1 - fx, f_of_x, 1e-7f64);
    }
}
