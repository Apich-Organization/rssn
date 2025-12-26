use assert_approx_eq::assert_approx_eq;
use rssn::numerical::calculus::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_partial_derivative() {

    let x = Expr::new_variable("x");

    let f = Expr::new_pow(
        x,
        Expr::new_constant(2.0),
    ); // f(x) = x^2
    let val = partial_derivative(
        &f, "x", 3.0,
    )
    .unwrap();

    assert_approx_eq!(
        val, 6.0, 1e-5f64
    );
}

#[test]

fn test_gradient() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // f(x,y) = x^2 + 2y
    let f = Expr::new_add(
        Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        ),
        Expr::new_mul(
            Expr::new_constant(2.0),
            y,
        ),
    );

    let grad = gradient(
        &f,
        &["x", "y"],
        &[2.0, 5.0],
    )
    .unwrap();

    // grad = [2x, 2] at (2,5) = [4, 2]
    assert_approx_eq!(
        grad[0], 4.0, 1e-5f64
    );

    assert_approx_eq!(
        grad[1], 2.0, 1e-5f64
    );
}

#[test]

fn test_jacobian() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // f1 = x*y, f2 = x^2 + y^2
    let f1 = Expr::new_mul(
        x.clone(),
        y.clone(),
    );

    let f2 = Expr::new_add(
        Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        ),
        Expr::new_pow(
            y,
            Expr::new_constant(2.0),
        ),
    );

    let jac = jacobian(
        &[f1, f2],
        &["x", "y"],
        &[1.0, 2.0],
    )
    .unwrap();

    // J = [[y, x], [2x, 2y]] at (1,2) = [[2, 1], [2, 4]]
    assert_approx_eq!(
        jac[0][0], 2.0, 1e-5f64
    );

    assert_approx_eq!(
        jac[0][1], 1.0, 1e-5f64
    );

    assert_approx_eq!(
        jac[1][0], 2.0, 1e-5f64
    );

    assert_approx_eq!(
        jac[1][1], 4.0, 1e-5f64
    );
}

#[test]

fn test_hessian() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // f = x^2 * y + y^3
    let f = Expr::new_add(
        Expr::new_mul(
            Expr::new_pow(
                x,
                Expr::new_constant(2.0),
            ),
            y.clone(),
        ),
        Expr::new_pow(
            y,
            Expr::new_constant(3.0),
        ),
    );

    let hess = hessian(
        &f,
        &["x", "y"],
        &[1.0, 2.0],
    )
    .unwrap();

    // fx = 2xy, fy = x^2 + 3y^2
    // fxx = 2y, fxy = 2x, fyx = 2x, fyy = 6y
    // H at (1,2) = [[4, 2], [2, 12]]
    assert_approx_eq!(
        hess[0][0], 4.0, 1e-4f64
    );

    assert_approx_eq!(
        hess[0][1], 2.0, 1e-4f64
    );

    assert_approx_eq!(
        hess[1][0], 2.0, 1e-4f64
    );

    assert_approx_eq!(
        hess[1][1], 12.0, 1e-4f64
    );
}

use proptest::prelude::*;

proptest! {
    #[test]
    fn proptest_partial_derivative_linear(a in -10.0..10.0f64, b in -10.0..10.0f64, x_val in -10.0..10.0f64) {
        let x = Expr::new_variable("x");
        // f(x) = ax + b
        let f = Expr::new_add(
            Expr::new_mul(Expr::new_constant(a), x),
            Expr::new_constant(b)
        );
        let val = partial_derivative(&f, "x", x_val).unwrap();
        // Derivative should be a
        assert_approx_eq!(val, a, 1e-5f64);
    }

    #[test]
    fn proptest_partial_derivative_quadratic(a in -5.0..5.0f64, x_val in -5.0..5.0f64) {
        let x = Expr::new_variable("x");
        // f(x) = ax^2
        let f = Expr::new_mul(
            Expr::new_constant(a),
            Expr::new_pow(x, Expr::new_constant(2.0))
        );
        let val = partial_derivative(&f, "x", x_val).unwrap();
        // Derivative should be 2ax
        assert_approx_eq!(val, 2.0 * a * x_val, 1e-4f64);
    }

    #[test]
    fn proptest_gradient_linear_combination(
        a in -5.0..5.0f64,
        b in -5.0..5.0f64,
        x_val in -5.0..5.0f64,
        y_val in -5.0..5.0f64
    ) {
        let x = Expr::new_variable("x");
        let y = Expr::new_variable("y");
        // f(x,y) = ax + by
        let f = Expr::new_add(
            Expr::new_mul(Expr::new_constant(a), x),
            Expr::new_mul(Expr::new_constant(b), y)
        );
        let grad = gradient(&f, &["x", "y"], &[x_val, y_val]).unwrap();
        // grad = [a, b]
        assert_approx_eq!(grad[0], a, 1e-5f64);
        assert_approx_eq!(grad[1], b, 1e-5f64);
    }
}
