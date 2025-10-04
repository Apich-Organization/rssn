// File: tests/numerical/calculus.rs

use assert_approx_eq::assert_approx_eq;
use rssn::numerical::calculus::gradient;
use rssn::symbolic::core::Expr;

/// Tests the gradient of a simple single-variable function, f(x) = x^2.
/// The gradient of x^2 is 2x. At x=3, the gradient should be 6.
#[test]
fn test_gradient_x_squared() {
    let x = Expr::Variable("x".to_string());
    let x_squared = Expr::Mul(Box::new(x.clone()), Box::new(x.clone()));

    let vars = ["x"];
    let point = [3.0];

    let grad = gradient(&x_squared, &vars, &point).unwrap();

    assert_eq!(grad.len(), 1);
    assert_approx_eq!(grad[0], 6.0, 1e-6);
}

/// Tests the gradient of a multivariate function, f(x, y) = x^2 + y^2.
/// The partial derivative with respect to x is 2x.
/// The partial derivative with respect to y is 2y.
/// At (x, y) = (1, 2), the gradient is (2, 4).
#[test]
fn test_gradient_x_squared_plus_y_squared() {
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());

    let x_squared = Expr::Mul(Box::new(x.clone()), Box::new(x.clone()));
    let y_squared = Expr::Mul(Box::new(y.clone()), Box::new(y.clone()));
    let f = Expr::Add(Box::new(x_squared), Box::new(y_squared));

    let vars = ["x", "y"];
    let point = [1.0, 2.0];

    let grad = gradient(&f, &vars, &point).unwrap();

    assert_eq!(grad.len(), 2);
    assert_approx_eq!(grad[0], 2.0, 1e-6);
    assert_approx_eq!(grad[1], 4.0, 1e-6);
}

/// Tests the gradient of a more complex function, f(x, y) = sin(x) + cos(y).
/// The partial derivative with respect to x is cos(x).
/// The partial derivative with respect to y is -sin(y).
/// At (x, y) = (0, PI/2), the gradient is (cos(0), -sin(PI/2)) = (1, -1).
#[test]
fn test_gradient_sin_x_plus_cos_y() {
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());

    let sin_x = Expr::Sin(Box::new(x.clone()));
    let cos_y = Expr::Cos(Box::new(y.clone()));
    let f = Expr::Add(Box::new(sin_x), Box::new(cos_y));

    let vars = ["x", "y"];
    let point = [0.0, std::f64::consts::PI / 2.0];

    let grad = gradient(&f, &vars, &point).unwrap();

    assert_eq!(grad.len(), 2);
    assert_approx_eq!(grad[0], 1.0, 1e-6);
    assert_approx_eq!(grad[1], -1.0, 1e-6);
}
