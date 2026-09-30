//! Comprehensive Calculus & Continuous Analysis domain integration tests for RSSN.
//!
//! Validates the unified Identity Transformation Operator architecture under `compute(expr, config)`.

use num_bigint::BigInt;
use rssn::compute::{
    compute, cos, d, definite_integral, exp, gradient, indefinite_integral, integral, limit, sin,
    ComputeConfig,
};
use rssn::symbolic::core::Expr;

#[test]
fn test_symbolic_calculus_derivative_sin_cos() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // d(sin(x))/dx => cos(x)
    let deriv = d(sin(x.clone()), "x");
    let res = compute(&deriv, &config_sym);

    eprintln!("d(sin(x))/dx result: {:?}", res);
    assert_eq!(res, cos(x));
}

#[test]
fn test_symbolic_calculus_derivative_polynomial() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // d(x^3)/dx => 3 * x^2
    let x_cubed = Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(3)));
    let res = compute(&d(x_cubed, "x"), &config_sym);

    eprintln!("d(x^3)/dx result: {:?}", res);
    // Should be a product 3 * x^2
    assert!(matches!(res, Expr::Mul(_, _)));
}

#[test]
fn test_symbolic_calculus_partial_derivative_multivariate() {
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // d(x * y)/dx => y
    let xy = Expr::new_mul(x.clone(), y.clone());
    let res = compute(&d(xy, "x"), &config_sym);

    eprintln!("d(x * y)/dx result: {:?}", res);
    assert_eq!(res, y);
}

#[test]
fn test_numerical_integration_quadrature_x_squared() {
    let x = Expr::new_variable("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(2)));
    let config_num = ComputeConfig::default()
        .with_calculus()
        .target_numerical(1e-7);

    // integral(x^2, "x", Some((0.0, 1.0))) => ~0.333333
    let int_expr = integral(x_sq, "x", Some((0.0, 1.0)));
    let res = compute(&int_expr, &config_num);

    eprintln!("Numerical integral(x^2, 0, 1) result: {:?}", res);
    match res {
        Expr::Constant(val) => {
            assert!(
                (val - 1.0 / 3.0).abs() < 1e-4,
                "Expected ~0.333333, got {val}"
            );
        }
        other => panic!("Expected numerical Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_integration_non_elementary_gaussian() {
    let x = Expr::new_variable("x");
    // exp(-x^2) has no elementary antiderivative; numerical quadrature is required
    let neg_x_sq = Expr::new_neg(Expr::new_pow(x.clone(), Expr::Constant(2.0)));
    let integrand = exp(neg_x_sq);
    let config_num = ComputeConfig::default()
        .with_calculus()
        .target_numerical(1e-7);

    let int_expr = integral(integrand, "x", Some((0.0, 1.0)));
    let res = compute(&int_expr, &config_num);

    eprintln!("Numerical integral(exp(-x^2), 0, 1) result: {:?}", res);
    // Exact value of \int_0^1 e^{-x^2} dx is sqrt(pi)/2 * erf(1) ~ 0.7468241328
    match res {
        Expr::Constant(val) => {
            assert!(
                (val - 0.746824).abs() < 1e-4,
                "Expected ~0.746824, got {val}"
            );
        }
        other => panic!("Expected numerical Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_symbolic_indefinite_integration() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // integral(x, "x", None) => x^2 / 2
    let int_expr = indefinite_integral(x.clone(), "x");
    let res = compute(&int_expr, &config_sym);

    eprintln!("Symbolic integral(x) result: {:?}", res);
    // Should be Div(x^2, 2) or Mul(1/2, x^2)
    assert!(
        matches!(res, Expr::Div(_, _) | Expr::Mul(_, _)),
        "Expected division or multiplication for x^2/2, got: {:?}",
        res
    );
}

#[test]
fn test_symbolic_definite_integration() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // definite_integral(x, "x", 0.0, 1.0) => 0.5 (or Rational 1/2)
    let int_expr = definite_integral(x.clone(), "x", 0.0, 1.0);
    let res = compute(&int_expr, &config_sym);

    eprintln!("Symbolic definite_integral(x, 0, 1) result: {:?}", res);
    if let Some(val) = res.to_f64() {
        assert!((val - 0.5).abs() < 1e-6);
    } else {
        panic!("Expected numerical or rational 0.5, got: {:?}", res);
    }
}

#[test]
fn test_numerical_derivative_with_binding() {
    let x = Expr::new_variable("x");
    // d(x^3)/dx at x=2 => 3 * 2^2 = 12.0
    let x_cubed = Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(3)));
    let config_num = ComputeConfig::default()
        .with_calculus()
        .bind("x", 2.0)
        .target_numerical(1e-7);

    let res = compute(&d(x_cubed, "x"), &config_num);

    eprintln!("Numerical d(x^3)/dx at x=2 result: {:?}", res);
    match res {
        Expr::Constant(val) => assert!((val - 12.0).abs() < 1e-5),
        Expr::BigInt(val) => assert_eq!(val, BigInt::from(12)),
        other => panic!("Expected 12.0, got: {:?}", other),
    }
}

#[test]
fn test_symbolic_limit_evaluation() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // lim_{x -> 2} (x^2 - 1) = 3
    let expr = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(2))),
        Expr::new_bigint(BigInt::from(1)),
    );
    let lim = limit(expr, "x", Expr::new_bigint(BigInt::from(2)));
    let res = compute(&lim, &config_sym);

    eprintln!("Symbolic limit result: {:?}", res);
    match res {
        Expr::BigInt(b) => assert_eq!(b, BigInt::from(3)),
        Expr::Constant(c) => assert!((c - 3.0).abs() < 1e-6),
        other => panic!("Expected 3, got: {:?}", other),
    }
}

#[test]
fn test_vector_gradient_operator() {
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // f(x, y) = x^2 + y^2 => grad = [2x, 2y]
    let f = Expr::new_add(
        Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(2))),
        Expr::new_pow(y.clone(), Expr::new_bigint(BigInt::from(2))),
    );

    let grad_expr = gradient(f, &["x", "y"]);
    let res = compute(&grad_expr, &config_sym);

    eprintln!("Symbolic gradient result: {:?}", res);
    match res {
        Expr::Vector(comps) => {
            assert_eq!(comps.len(), 2);
            eprintln!("Gradient components: {:?}", comps);
        }
        other => panic!("Expected Expr::Vector, got: {:?}", other),
    }
}

#[test]
fn test_dag_structural_sharing_no_explosion() {
    let x = Expr::new_variable("x");
    let config_sym = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    // Build nested expression: f_0 = sin(x), f_{k+1} = d(f_k, "x")
    // Each step should remain compact via EGraph equivalences and Arc sharing.
    let mut current = sin(x.clone());
    for _ in 0..4 {
        current = compute(&d(current, "x"), &config_sym);
    }

    // 4th derivative of sin(x) is sin(x)
    eprintln!("4th derivative of sin(x) result: {:?}", current);
    assert_eq!(current, sin(x));
}
