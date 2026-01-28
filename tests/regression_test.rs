use std::sync::Arc;

use num_bigint::BigInt;
use rssn::symbolic::calculus::differentiate;
use rssn::symbolic::core::Expr;

#[test]

fn test_differentiate_x_squared_stack_overflow()
 {

    let x =
        Expr::Variable("x".to_string());

    let x2 = Expr::new_mul(
        x.clone(),
        x.clone(),
    );

    let d = differentiate(&x2, "x");

    // The derivative of x^2 is 2*x.
    // The simplification process might result in Constant(2.0) or BigInt(2).
    let two_const = Expr::new_constant(2.0);

    let expected_const = Expr::new_mul(
        two_const,
        x.clone(),
    );

    let two_int =
        Expr::new_bigint(BigInt::from(2));

    let expected_int = Expr::new_mul(
        two_int,
        x.clone(),
    );

    println!(
        "Derivative: {:?}",
        d
    );

    println!(
        "Expected (const): {:?}",
        expected_const
    );

    println!(
        "Expected (int): {:?}",
        expected_int
    );

    assert!(
        d == expected_const
            || d == expected_int
    );
}

#[test]

fn test_normalization() {

    let x = Expr::new_variable("x");

    let two = Expr::new_constant(2.0);

    let a = Expr::new_mul(
        x.clone(),
        two.clone(),
    );

    let b = Expr::new_mul(
        two.clone(),
        x.clone(),
    );

    assert_eq!(
        a, b,
        "DAG normalization (A*B = \
         B*A) failed!"
    );
}
