use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify::simplify;

#[test]
fn test_simplify_one_plus_one_expr() {
    let expr = Expr::new_add(Expr::BigInt(BigInt::from(1)), Expr::BigInt(BigInt::from(1)));
    let simplified = simplify(expr);
    assert_eq!(simplified, Expr::BigInt(BigInt::from(2)));
}

#[test]
fn test_simplify_x_plus_2x_dag() {
    let x = Expr::new_variable("x");
    let two_x = Expr::new_mul(Expr::new_bigint(BigInt::from(2)), x.clone());
    let expr = Expr::new_add(x, two_x);
    let simplified = simplify(expr);
    let expected = Expr::new_mul(Expr::new_bigint(BigInt::from(3)), Expr::new_variable("x"));
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_x_plus_2x_expr() {
    let x = Expr::Variable("x".to_string());
    let two_x = Expr::Mul(Expr::BigInt(BigInt::from(2)).into(), x.clone().into());
    let expr = Expr::Add(x.into(), two_x.into());
    let simplified = simplify(expr);
    let expected = Expr::Mul(
        Expr::BigInt(BigInt::from(3)).into(),
        Expr::Variable("x".to_string()).into(),
    );
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_add() {
    let two = Expr::new_constant(2.0);
    let one = Expr::new_constant(1.0);
    let three = Expr::new_constant(3.0);
    let zero = Expr::new_constant(0.0);
    let expr = Expr::new_add(&one, Expr::new_add(&one, zero));
    let simplified = simplify(expr);
    let expected = two;
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_mul() {
    let two = Expr::new_constant(2.0);
    let one = Expr::new_constant(1.0);
    let three = Expr::new_constant(3.0);
    let five = Expr::new_constant(5.0);
    let zero = Expr::new_constant(0.0);
    let expr = Expr::new_mul(one, Expr::new_add(three, two));
    let simplified = simplify(expr);
    let expected = five;
    assert_eq!(simplified, expected);
}