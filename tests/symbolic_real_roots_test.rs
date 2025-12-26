use rssn::symbolic::core::Expr;
use rssn::symbolic::polynomial::expr_to_sparse_poly;
use rssn::symbolic::real_roots::*;

fn create_poly_expr(coeffs: Vec<f64>, var: &str) -> Expr {

    let mut expr = Expr::Constant(0.0);

    let x = Expr::Variable(var.to_string());

    for (i, coeff) in coeffs
        .iter()
        .enumerate()
    {

        if *coeff != 0.0 {

            let term = Expr::new_mul(
                Expr::Constant(*coeff),
                Expr::new_pow(x.clone(), Expr::Constant(i as f64)),
            );

            expr = Expr::new_add(expr, term);
        }
    }

    expr
}

#[test]

fn test_sturm_sequence_simple() {

    // x^2 - 2
    let x = Expr::Variable("x".to_string());

    let expr = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::Constant(2.0),
    );

    let poly = expr_to_sparse_poly(&expr, &["x"]);

    let seq = sturm_sequence(&poly, "x");

    // Sequence should be:
    // P0 = x^2 - 2
    // P1 = 2x
    // P2 = 2 (remainder of (x^2-2)/(2x) is -2, negated is 2)
    // Or similar, depending on normalization.

    assert!(!seq.is_empty());

    assert_eq!(seq.len(), 3);
}

#[test]

fn test_count_real_roots() {

    // x^2 - 2 has roots at sqrt(2) and -sqrt(2)
    // sqrt(2) approx 1.414
    let x = Expr::Variable("x".to_string());

    let expr = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::Constant(2.0),
    );

    let poly = expr_to_sparse_poly(&expr, &["x"]);

    // Interval [0, 2] contains sqrt(2)
    let count = count_real_roots_in_interval(&poly, "x", 0.0, 2.0).unwrap();

    assert_eq!(count, 1);

    // Interval [-2, 0] contains -sqrt(2)
    let count = count_real_roots_in_interval(&poly, "x", -2.0, 0.0).unwrap();

    assert_eq!(count, 1);

    // Interval [-2, 2] contains both
    let count = count_real_roots_in_interval(&poly, "x", -2.0, 2.0).unwrap();

    assert_eq!(count, 2);

    // Interval [2, 3] contains none
    let count = count_real_roots_in_interval(&poly, "x", 2.0, 3.0).unwrap();

    assert_eq!(count, 0);
}

#[test]

fn test_isolate_real_roots() {

    // x^2 - 2
    let x = Expr::Variable("x".to_string());

    let expr = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::Constant(2.0),
    );

    let poly = expr_to_sparse_poly(&expr, &["x"]);

    let roots = isolate_real_roots(&poly, "x", 0.001).unwrap();

    assert_eq!(roots.len(), 2);

    // Check first root (-sqrt(2))
    let (a1, b1) = roots[0];

    assert!(a1 < -1.4 && b1 > -1.5);

    assert!((a1 - (-1.4142)).abs() < 0.1);

    // Check second root (sqrt(2))
    let (a2, b2) = roots[1];

    assert!(a2 < 1.5 && b2 > 1.4);

    assert!((a2 - 1.4142).abs() < 0.1);
}

#[test]

fn test_multiple_roots() {

    // x^3 - x = x(x-1)(x+1) -> roots at -1, 0, 1
    let x = Expr::Variable("x".to_string());

    let expr = Expr::new_sub(Expr::new_pow(x.clone(), Expr::Constant(3.0)), x.clone());

    let poly = expr_to_sparse_poly(&expr, &["x"]);

    let count = count_real_roots_in_interval(&poly, "x", -2.0, 2.0).unwrap();

    assert_eq!(count, 3);

    let roots = isolate_real_roots(&poly, "x", 0.001).unwrap();

    assert_eq!(roots.len(), 3);

    // Roots should be close to -1, 0, 1
    assert!((roots[0].0 - (-1.0)).abs() < 0.1);

    assert!((roots[1].0 - 0.0).abs() < 0.1);

    assert!((roots[2].0 - 1.0).abs() < 0.1);
}
