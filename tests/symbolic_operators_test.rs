use std::sync::Arc;

use rssn::symbolic::core::Expr;

#[test]

fn test_basic_arithmetic_operators() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let c2 = Expr::new_constant(2.0);

    // Add
    let add = x.clone() + y.clone();

    assert_eq!(
        format!("{:?}", add),
        "(x + y)"
    );

    // Sub
    let sub = x.clone() - y.clone();

    assert_eq!(
        format!("{:?}", sub),
        "(x - y)"
    );

    // Mul
    let mul = x.clone() * y.clone();

    assert_eq!(
        format!("{:?}", mul),
        "(x * y)"
    );

    // Div
    let div = x.clone() / y.clone();

    assert_eq!(
        format!("{:?}", div),
        "(x / y)"
    );

    // Neg
    let neg = -x.clone();

    assert_eq!(
        format!("{:?}", neg),
        "-(x)"
    );

    // Expr + f64
    let add_f64 = x.clone() + 2.0;

    assert_eq!(
        format!("{:?}", add_f64),
        "(2 + x)"
    );

    // f64 + Expr
    let f64_add = 2.0 + x.clone();

    assert_eq!(
        format!("{:?}", f64_add),
        "(2 + x)"
    );

    // Expr * f64
    let mul_f64 = x.clone() * 3.0;

    assert_eq!(
        format!("{:?}", mul_f64),
        "(3 * x)"
    );
}

#[test]

fn test_method_calls() {

    let x = Expr::new_variable("x");

    // Trig
    assert_eq!(
        format!("{:?}", x.sin()),
        "sin(x)"
    );

    assert_eq!(
        format!("{:?}", x.cos()),
        "cos(x)"
    );

    assert_eq!(
        format!("{:?}", x.tan()),
        "tan(x)"
    );

    // Exp/Log
    assert_eq!(
        format!("{:?}", x.exp()),
        "exp(x)"
    );

    assert_eq!(
        format!("{:?}", x.ln()),
        "ln(x)"
    );

    // Pow/Sqrt
    assert_eq!(
        format!(
            "{:?}",
            x.pow(Expr::new_constant(
                2.0
            ))
        ),
        "(x^(2))"
    );

    assert_eq!(
        format!("{:?}", x.sqrt()),
        "sqrt(x)"
    );

    // Abs
    assert_eq!(
        format!("{:?}", x.abs()),
        "|x|"
    );

    // Inverse Trig
    assert_eq!(
        format!("{:?}", x.asin()),
        "asin(x)"
    );

    assert_eq!(
        format!("{:?}", x.acos()),
        "acos(x)"
    );

    assert_eq!(
        format!("{:?}", x.atan()),
        "atan(x)"
    );

    // Hyperbolic
    assert_eq!(
        format!("{:?}", x.sinh()),
        "sinh(x)"
    );

    assert_eq!(
        format!("{:?}", x.cosh()),
        "cosh(x)"
    );

    assert_eq!(
        format!("{:?}", x.tanh()),
        "tanh(x)"
    );
}

#[test]

fn test_complex_expressions() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // (x + y) * 2.0
    let expr1 =
        (x.clone() + y.clone()) * 2.0;

    // Expected structure: ((x + y) * 2)
    assert_eq!(
        format!("{:?}", expr1),
        "(2 * (x + y))"
    );

    // sin(x)^2 + cos(x)^2
    let expr2 = x.sin().pow(
        Expr::new_constant(2.0),
    ) + x.cos().pow(
        Expr::new_constant(2.0),
    );

    // Expected: ((sin(x)^(2)) + (cos(x)^(2)))
    assert_eq!(
        format!("{:?}", expr2),
        "((sin(x)^(2)) + (cos(x)^(2)))"
    );
}

#[test]

fn test_operator_references() {

    let a = Expr::new_constant(1.0);

    let b = Expr::new_constant(2.0);

    // &Expr + &Expr
    let c = &a + &b;

    assert_eq!(
        format!("{:?}", c),
        "(1 + 2)"
    );

    // Expr + &Expr
    let d = a.clone() + &b;

    assert_eq!(
        format!("{:?}", d),
        "(1 + 2)"
    );

    // &Expr + Expr
    let e = &a + b.clone();

    assert_eq!(
        format!("{:?}", e),
        "(1 + 2)"
    );
}
