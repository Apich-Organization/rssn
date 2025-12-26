use rssn::symbolic::core::Expr;
use rssn::symbolic::proof::*;
use std::collections::HashMap;

#[test]

fn test_verify_equation_solution() {

    // x^2 - 4 = 0, solution x = 2
    let eq = Expr::new_sub(
        Expr::new_pow(Expr::new_variable("x"), Expr::Constant(2.0)),
        Expr::Constant(4.0),
    );

    let mut solution = HashMap::new();

    solution.insert("x".to_string(), Expr::Constant(2.0));

    assert!(verify_equation_solution(&[eq], &solution, &[]));
}

#[test]

fn test_verify_indefinite_integral() {

    // int(2x dx) = x^2
    let integrand = Expr::new_mul(Expr::Constant(2.0), Expr::new_variable("x"));

    let result = Expr::new_pow(Expr::new_variable("x"), Expr::Constant(2.0));

    assert!(verify_indefinite_integral(&integrand, &result, "x"));
}

#[test]

fn test_verify_derivative() {

    // d/dx(x^3) = 3x^2
    let f = Expr::new_pow(Expr::new_variable("x"), Expr::Constant(3.0));

    let df = Expr::new_mul(
        Expr::Constant(3.0),
        Expr::new_pow(Expr::new_variable("x"), Expr::Constant(2.0)),
    );

    assert!(verify_derivative(&f, &df, "x"));
}

#[test]

fn test_verify_matrix_inverse() {

    // [[2, 0]; [0, 2]] inverse is [[0.5, 0]; [0, 0.5]]
    let a = Expr::Matrix(vec![
        vec![
            Expr::Constant(2.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(2.0),
        ],
    ]);

    let inv = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.5),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.5),
        ],
    ]);

    assert!(verify_matrix_inverse(&a, &inv));
}

#[test]

fn test_verify_limit() {

    // lim_{x->0} sin(x)/x = 1
    let f = Expr::new_div(
        Expr::new_sin(Expr::new_variable("x")),
        Expr::new_variable("x"),
    );

    let target = Expr::Constant(0.0);

    let l = Expr::Constant(1.0);

    assert!(verify_limit(&f, "x", &target, &l));
}
