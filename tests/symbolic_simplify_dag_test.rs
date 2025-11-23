use num_bigint::BigInt;
use rssn::symbolic::core::{Expr, DAG_MANAGER};
use rssn::symbolic::simplify_dag::simplify;

#[test]
fn test_simplify_basic_arithmetic() {
    // 1 + 1 -> 2
    let expr = Expr::new_add(Expr::new_constant(1.0), Expr::new_constant(1.0));
    let simplified = simplify(&expr);
    // simplified should be a DAG node representing Constant(2.0)
    // We can compare it with Expr::new_constant(2.0)
    assert_eq!(simplified, Expr::new_constant(2.0));

    // x + 0 -> x
    let x = Expr::new_variable("x");
    let expr = Expr::new_add(&x, Expr::new_constant(0.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, x);

    // x * 1 -> x
    let expr = Expr::new_mul(&x, Expr::new_constant(1.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, x);

    // x * 0 -> 0
    let expr = Expr::new_mul(&x, Expr::new_constant(0.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, Expr::new_constant(0.0));
}

#[test]
fn test_simplify_trig() {
    // sin(0) -> 0
    let expr = Expr::new_sin(Expr::new_constant(0.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, Expr::new_constant(0.0));

    // cos(0) -> 1
    let expr = Expr::new_cos(Expr::new_constant(0.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, Expr::new_constant(1.0));
}

#[test]
fn test_simplify_power() {
    // x^1 -> x
    let x = Expr::new_variable("x");
    let expr = Expr::new_pow(&x, Expr::new_constant(1.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, x);

    // x^0 -> 1
    let expr = Expr::new_pow(&x, Expr::new_constant(0.0));
    let simplified = simplify(&expr);
    assert_eq!(simplified, Expr::new_constant(1.0));
}

#[test]
fn test_simplify_nested() {
    // (x + x) * 2 -> 4*x
    let x = Expr::new_variable("x");
    let expr = Expr::new_mul(Expr::new_add(&x, &x), Expr::new_constant(2.0));
    let simplified = simplify(&expr);
    // Expect 4*x
    println!("Simplified: {:?}", simplified);

    let expected = Expr::new_mul(Expr::new_constant(4.0), &x);
    assert_eq!(simplified, expected);
}
