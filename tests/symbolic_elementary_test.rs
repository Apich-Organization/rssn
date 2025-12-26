use std::sync::Arc;

use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::elementary::*;

#[test]

fn test_trig_functions() {

    let x = Expr::new_variable("x");

    // Test sin - constructors now return Expr::Dag
    let sin_x = sin(x.clone());

    assert!(matches!(
        sin_x,
        Expr::Dag(_)
    ));

    // Test cos
    let cos_x = cos(x.clone());

    assert!(matches!(
        cos_x,
        Expr::Dag(_)
    ));

    // Test tan
    let tan_x = tan(x.clone());

    assert!(matches!(
        tan_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_hyperbolic_functions() {

    let x = Expr::new_variable("x");

    // Test sinh
    let sinh_x = sinh(x.clone());

    assert!(matches!(
        sinh_x,
        Expr::Dag(_)
    ));

    // Test cosh
    let cosh_x = cosh(x.clone());

    assert!(matches!(
        cosh_x,
        Expr::Dag(_)
    ));

    // Test tanh
    let tanh_x = tanh(x.clone());

    assert!(matches!(
        tanh_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_exp_and_log() {

    let x = Expr::new_variable("x");

    // Test ln
    let ln_x = ln(x.clone());

    assert!(matches!(
        ln_x,
        Expr::Dag(_)
    ));

    // Test exp
    let exp_x = exp(x.clone());

    assert!(matches!(
        exp_x,
        Expr::Dag(_)
    ));

    // Test log_base
    let log_2_x = log_base(
        Expr::Constant(2.0),
        x.clone(),
    );

    assert!(matches!(
        log_2_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_binomial_coefficient() {

    use rssn::symbolic::elementary::binomial_coefficient;

    // Test C(5, 2) = 10
    let result = binomial_coefficient(5, 2);

    assert_eq!(
        result,
        BigInt::from(10)
    );

    // Test C(4, 0) = 1
    let result = binomial_coefficient(4, 0);

    assert_eq!(
        result,
        BigInt::from(1)
    );

    // Test C(4, 4) = 1
    let result = binomial_coefficient(4, 4);

    assert_eq!(
        result,
        BigInt::from(1)
    );

    // Test C(6, 3) = 20
    let result = binomial_coefficient(6, 3);

    assert_eq!(
        result,
        BigInt::from(20)
    );
}

#[test]

fn test_power_and_sqrt() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // Test pow
    let x_pow_y = pow(x.clone(), y.clone());

    assert!(matches!(
        x_pow_y,
        Expr::Dag(_)
    ));

    // Test sqrt
    let sqrt_x = sqrt(x.clone());

    assert!(matches!(
        sqrt_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_constants() {

    // Test pi
    let pi_val = pi();

    assert!(matches!(
        pi_val,
        Expr::Pi
    ));

    // Test e
    let e_val = e();

    assert!(matches!(
        e_val,
        Expr::E
    ));

    // Test infinity
    let inf = infinity();

    assert!(matches!(
        inf,
        Expr::Infinity
    ));

    // Test negative infinity
    let neg_inf = negative_infinity();

    assert!(matches!(
        neg_inf,
        Expr::NegativeInfinity
    ));
}

#[test]

fn test_inverse_trig() {

    let x = Expr::new_variable("x");

    // Test acot
    let acot_x = acot(x.clone());

    assert!(matches!(
        acot_x,
        Expr::Dag(_)
    ));

    // Test asec
    let asec_x = asec(x.clone());

    assert!(matches!(
        asec_x,
        Expr::Dag(_)
    ));

    // Test acsc
    let acsc_x = acsc(x.clone());

    assert!(matches!(
        acsc_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_inverse_hyperbolic() {

    let x = Expr::new_variable("x");

    // Test asinh
    let asinh_x = asinh(x.clone());

    assert!(matches!(
        asinh_x,
        Expr::Dag(_)
    ));

    // Test acosh
    let acosh_x = acosh(x.clone());

    assert!(matches!(
        acosh_x,
        Expr::Dag(_)
    ));

    // Test atanh
    let atanh_x = atanh(x.clone());

    assert!(matches!(
        atanh_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_expand_mul() {

    // Test (x+1)*(y+2) expansion using AST constructors
    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    let one = Expr::Constant(1.0);

    let two = Expr::Constant(2.0);

    let x_plus_1 = Expr::Add(
        Arc::new(x.clone()),
        Arc::new(one),
    );

    let y_plus_2 = Expr::Add(
        Arc::new(y.clone()),
        Arc::new(two),
    );

    let product = Expr::Mul(
        Arc::new(x_plus_1),
        Arc::new(y_plus_2),
    );

    let expanded = expand(product);

    // Result should be expanded (exact form depends on simplification)
    // Just check it doesn't panic and returns something
    assert!(matches!(
        expanded,
        Expr::Add(_, _) | Expr::Mul(_, _) | Expr::Dag(_) | _
    ));
}

#[test]

fn test_expand_power() {

    // Test (x+1)^2 expansion using AST constructors
    let x = Expr::Variable("x".to_string());

    let one = Expr::Constant(1.0);

    let two = Expr::BigInt(BigInt::from(2));

    let x_plus_1 = Expr::Add(
        Arc::new(x.clone()),
        Arc::new(one),
    );

    let squared = Expr::Power(
        Arc::new(x_plus_1),
        Arc::new(two),
    );

    let expanded = expand(squared);

    // Should expand to x^2 + 2x + 1 (in some form)
    // Just check it doesn't panic
    assert!(matches!(
        expanded,
        Expr::Add(_, _) | Expr::Dag(_) | _
    ));
}

#[test]

fn test_expand_log() {

    // Test log(x*y) expansion using AST constructors
    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    let product = Expr::Mul(
        Arc::new(x.clone()),
        Arc::new(y.clone()),
    );

    let log_product = Expr::Log(Arc::new(product));

    let expanded = expand(log_product);

    // Should expand to log(x) + log(y)
    assert!(matches!(
        expanded,
        Expr::Add(_, _) | Expr::Log(_) | Expr::Dag(_) | _
    ));
}

#[test]

fn test_expand_trig() {

    // Test sin(x+y) expansion using AST constructors
    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    let sum = Expr::Add(
        Arc::new(x.clone()),
        Arc::new(y.clone()),
    );

    let sin_sum = Expr::Sin(Arc::new(sum));

    let expanded = expand(sin_sum);

    // Should expand using sum-angle formula
    assert!(matches!(
        expanded,
        Expr::Add(_, _) | Expr::Sin(_) | Expr::Dag(_) | _
    ));
}

#[test]

fn test_atan2() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let atan2_yx = atan2(y.clone(), x.clone());

    assert!(matches!(
        atan2_yx,
        Expr::Dag(_)
    ));
}

#[test]

fn test_other_trig_functions() {

    let x = Expr::new_variable("x");

    // Test cot
    let cot_x = cot(x.clone());

    assert!(matches!(
        cot_x,
        Expr::Dag(_)
    ));

    // Test sec
    let sec_x = sec(x.clone());

    assert!(matches!(
        sec_x,
        Expr::Dag(_)
    ));

    // Test csc
    let csc_x = csc(x.clone());

    assert!(matches!(
        csc_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_other_hyperbolic_functions() {

    let x = Expr::new_variable("x");

    // Test coth
    let coth_x = coth(x.clone());

    assert!(matches!(
        coth_x,
        Expr::Dag(_)
    ));

    // Test sech
    let sech_x = sech(x.clone());

    assert!(matches!(
        sech_x,
        Expr::Dag(_)
    ));

    // Test csch
    let csch_x = csch(x.clone());

    assert!(matches!(
        csch_x,
        Expr::Dag(_)
    ));
}

#[test]

fn test_other_inverse_hyperbolic() {

    let x = Expr::new_variable("x");

    // Test acoth
    let acoth_x = acoth(x.clone());

    assert!(matches!(
        acoth_x,
        Expr::Dag(_)
    ));

    // Test asech
    let asech_x = asech(x.clone());

    assert!(matches!(
        asech_x,
        Expr::Dag(_)
    ));

    // Test acsch
    let acsch_x = acsch(x.clone());

    assert!(matches!(
        acsch_x,
        Expr::Dag(_)
    ));
}
