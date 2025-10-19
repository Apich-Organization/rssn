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