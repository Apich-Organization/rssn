use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;
use rssn::prelude::*;
use std::collections::HashMap;

#[test]

fn test_basic_eval() {

    let mut vars = HashMap::new();

    vars.insert("x".to_string(), 2.0);

    // x + 3
    let expr = Expr::new_add(
        Expr::new_variable("x"),
        Expr::new_constant(3.0),
    );

    assert_eq!(
        numerical_eval_expr(&expr, &vars).unwrap(),
        5.0
    );

    // x^2
    let expr = Expr::new_pow(
        Expr::new_variable("x"),
        Expr::new_constant(2.0),
    );

    assert_eq!(
        numerical_eval_expr(&expr, &vars).unwrap(),
        4.0
    );
}

#[test]

fn test_trig_eval() {

    let vars = HashMap::<String, f64>::new();

    let expr = Expr::new_sin(Expr::Pi);

    assert_approx_eq!(
        numerical_eval_expr(&expr, &vars).unwrap(),
        0.0
    );

    let expr = Expr::new_cos(Expr::Pi);

    assert_approx_eq!(
        numerical_eval_expr(&expr, &vars).unwrap(),
        -1.0
    );
}

#[test]

fn test_domain_errors() {

    let vars = HashMap::new();

    // 1/0
    let expr = Expr::new_div(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    assert!(numerical_eval_expr(&expr, &vars).is_err());

    // sqrt(-1)
    let expr = Expr::new_sqrt(Expr::new_constant(
        -1.0,
    ));

    assert!(numerical_eval_expr(&expr, &vars).is_err());

    // ln(0)
    let expr = Expr::new_log(Expr::new_constant(
        0.0,
    ));

    assert!(numerical_eval_expr(&expr, &vars).is_err());
}

#[test]

fn test_pure_functions() {

    assert_approx_eq!(
        numerical_sin(std::f64::consts::PI),
        0.0
    );

    assert_approx_eq!(
        numerical_cos(std::f64::consts::PI),
        -1.0
    );

    assert_eq!(
        numerical_abs(-5.0),
        5.0
    );

    assert_eq!(
        numerical_sqrt(16.0),
        4.0
    );

    assert_approx_eq!(
        numerical_exp(1.0),
        std::f64::consts::E
    );
}

proptest! {
    #[test]
    fn prop_sin_range(x in -1000.0..1000.0) {
        let res = numerical_sin(x);
        assert!(res >= -1.0 - 1e-15 && res <= 1.0 + 1e-15);
    }

    #[test]
    fn prop_exp_ln_inverse(x in 1e-5..100.0) {
        let res = numerical_exp(numerical_ln(x));
        assert_approx_eq!(res, x, 1e-7);
    }

    #[test]
    fn prop_square_sqrt(x in 0.0..1000.0) {
        let res = numerical_pow(numerical_sqrt(x), 2.0);
        assert_approx_eq!(res, x, 1e-7);
    }
}
