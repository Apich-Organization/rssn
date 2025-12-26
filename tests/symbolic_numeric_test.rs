use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use rssn::symbolic::core::Expr;
use rssn::symbolic::numeric::evaluate_numerical;

#[test]

fn test_evaluate_constant() {

    let expr = Expr::new_constant(3.14);

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 3.14).abs() < 1e-10);
}

#[test]

fn test_evaluate_bigint() {

    let expr = Expr::BigInt(BigInt::from(42));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 42.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_rational() {

    // 3/4 = 0.75
    let expr = Expr::Rational(BigRational::new(BigInt::from(3), BigInt::from(4)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 0.75).abs() < 1e-10);
}

#[test]

fn test_evaluate_pi() {

    let expr = Expr::Pi;

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::PI).abs() < 1e-10);
}

#[test]

fn test_evaluate_e() {

    let expr = Expr::E;

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::E).abs() < 1e-10);
}

#[test]

fn test_evaluate_arithmetic_add() {

    // 2 + 3 = 5
    let expr = Expr::new_add(Expr::new_constant(2.0), Expr::new_constant(3.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 5.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_arithmetic_sub() {

    // 5 - 3 = 2
    let expr = Expr::new_sub(Expr::new_constant(5.0), Expr::new_constant(3.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 2.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_arithmetic_mul() {

    // 4 * 5 = 20
    let expr = Expr::new_mul(Expr::new_constant(4.0), Expr::new_constant(5.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 20.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_arithmetic_div() {

    // 10 / 2 = 5
    let expr = Expr::new_div(Expr::new_constant(10.0), Expr::new_constant(2.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 5.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_power() {

    // 2^3 = 8
    let expr = Expr::new_pow(Expr::new_constant(2.0), Expr::new_constant(3.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 8.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_sqrt() {

    // sqrt(16) = 4
    let expr = Expr::new_sqrt(Expr::new_constant(16.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 4.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_exp() {

    // e^1 = e
    let expr = Expr::new_exp(Expr::new_constant(1.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::E).abs() < 1e-10);
}

#[test]

fn test_evaluate_log() {

    // ln(e) = 1
    let expr = Expr::new_log(Expr::E);

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 1.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_abs() {

    // abs(-5) = 5
    let expr = Expr::new_abs(Expr::new_constant(-5.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 5.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_neg() {

    // -(-3) = 3
    let expr = Expr::new_neg(Expr::new_constant(-3.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 3.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_sin_pi_over_6() {

    // sin(π/6) = 0.5
    let expr = Expr::new_sin(Expr::new_div(Expr::Pi, Expr::new_constant(6.0)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 0.5).abs() < 1e-10);
}

#[test]

fn test_evaluate_sin_pi_over_2() {

    // sin(π/2) = 1
    let expr = Expr::new_sin(Expr::new_div(Expr::Pi, Expr::new_constant(2.0)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 1.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_cos_pi() {

    // cos(π) = -1
    let expr = Expr::new_cos(Expr::Pi);

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - (-1.0)).abs() < 1e-10);
}

#[test]

fn test_evaluate_tan_pi_over_4() {

    // tan(π/4) = 1
    let expr = Expr::new_tan(Expr::new_div(Expr::Pi, Expr::new_constant(4.0)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 1.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_arcsin() {

    // arcsin(0.5) ≈ π/6
    let expr = Expr::new_arcsin(Expr::new_constant(0.5));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::FRAC_PI_6).abs() < 1e-10);
}

#[test]

fn test_evaluate_arccos() {

    // arccos(0) = π/2
    let expr = Expr::new_arccos(Expr::new_constant(0.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
}

#[test]

fn test_evaluate_arctan() {

    // arctan(1) = π/4
    let expr = Expr::new_arctan(Expr::new_constant(1.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
}

#[test]

fn test_evaluate_sinh() {

    // sinh(0) = 0
    let expr = Expr::new_sinh(Expr::new_constant(0.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!(result.abs() < 1e-10);
}

#[test]

fn test_evaluate_cosh() {

    // cosh(0) = 1
    let expr = Expr::new_cosh(Expr::new_constant(0.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 1.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_tanh() {

    // tanh(0) = 0
    let expr = Expr::new_tanh(Expr::new_constant(0.0));

    let result = evaluate_numerical(&expr).unwrap();

    assert!(result.abs() < 1e-10);
}

#[test]

fn test_evaluate_factorial() {

    // 5! = 120
    use std::sync::Arc;

    let expr = Expr::Factorial(Arc::new(Expr::new_constant(5.0)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 120.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_floor() {

    // floor(3.7) = 3
    use std::sync::Arc;

    let expr = Expr::Floor(Arc::new(Expr::new_constant(3.7)));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 3.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_complex_expression() {

    // (2 + 3) * 4 = 20
    let expr = Expr::new_mul(
        Expr::new_add(Expr::new_constant(2.0), Expr::new_constant(3.0)),
        Expr::new_constant(4.0),
    );

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 20.0).abs() < 1e-10);
}

#[test]

fn test_evaluate_nested_expression() {

    // sqrt(2^2 + 3^2) = sqrt(4 + 9) = sqrt(13)
    let expr = Expr::new_sqrt(Expr::new_add(
        Expr::new_pow(Expr::new_constant(2.0), Expr::new_constant(2.0)),
        Expr::new_pow(Expr::new_constant(3.0), Expr::new_constant(2.0)),
    ));

    let result = evaluate_numerical(&expr).unwrap();

    assert!((result - 13.0f64.sqrt()).abs() < 1e-10);
}

#[test]

fn test_evaluate_variable_returns_none() {

    // Variables cannot be evaluated numerically
    let expr = Expr::new_variable("x");

    let result = evaluate_numerical(&expr);

    assert!(result.is_none());
}

#[test]

fn test_evaluate_symbolic_expression_returns_none() {

    // Symbolic expressions with variables cannot be evaluated
    let expr = Expr::new_add(Expr::new_variable("x"), Expr::new_constant(5.0));

    let result = evaluate_numerical(&expr);

    assert!(result.is_none());
}
