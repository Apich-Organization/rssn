use num_traits::ToPrimitive;
use rssn::symbolic::core::Expr;
use rssn::symbolic::numeric::evaluate_numerical;
use rssn::symbolic::series::{
    analytic_continuation, analyze_convergence, asymptotic_expansion, fourier_series,
    laurent_series, product, summation, taylor_series,
};
use std::sync::Arc;

fn assert_is_value(
    expr: &Expr,
    expected: f64,
) {

    let val = evaluate_numerical(expr).unwrap_or_else(|| {
        panic!(
            "Failed to evaluate expression: {:?}",
            expr
        )
    });

    assert!(
        (val - expected).abs() < 1e-5,
        "Expected {}, got {}",
        expected,
        val
    );
}

#[test]

fn test_taylor_series_exp() {

    // Taylor series of e^x around 0 at order 3: 1 + x + x^2/2 + x^3/6
    let x = Expr::new_variable("x");

    let expr = Expr::new_exp(x.clone());

    let center = Expr::new_constant(0.0);

    let series = taylor_series(
        &expr, "x", &center, 3,
    );

    // Evaluate at x = 0.5
    // e^0.5 = 1.64872
    // Approx: 1 + 0.5 + 0.25/2 + 0.125/6 = 1.5 + 0.125 + 0.020833 = 1.645833

    let val_at_0_5 = rssn::symbolic::calculus::substitute(
        &series,
        "x",
        &Expr::new_constant(0.5),
    );

    assert_is_value(
        &val_at_0_5,
        1.645833,
    );
}

#[test]

fn test_taylor_series_sin() {

    // Taylor series of sin(x) around 0 at order 3: x - x^3/6
    let x = Expr::new_variable("x");

    let expr = Expr::new_sin(x.clone());

    let center = Expr::new_constant(0.0);

    let series = taylor_series(
        &expr, "x", &center, 3,
    );

    // Evaluate at x = 0.1
    // sin(0.1) ≈ 0.0998334
    // Approx: 0.1 - 0.001/6 = 0.1 - 0.0001666 = 0.0998333

    let val_at_0_1 = rssn::symbolic::calculus::substitute(
        &series,
        "x",
        &Expr::new_constant(0.1),
    );

    assert_is_value(
        &val_at_0_1,
        0.0998333,
    );
}

#[test]

fn test_summation_finite() {

    // Sum of i from 1 to 5: 1+2+3+4+5 = 15
    let i = Expr::new_variable("i");

    let sum = summation(
        &i,
        "i",
        &Expr::new_constant(1.0),
        &Expr::new_constant(5.0),
    );

    assert_is_value(&sum, 15.0);
}

#[test]

fn test_product_finite() {

    // Product of i from 1 to 4: 1*2*3*4 = 24
    let i = Expr::new_variable("i");

    let prod = product(
        &i,
        "i",
        &Expr::new_constant(1.0),
        &Expr::new_constant(4.0),
    );

    assert_is_value(&prod, 24.0);
}

#[test]

fn test_asymptotic_expansion() {

    // Expansion of (x + 1) / (x - 1) at infinity
    // = (1 + 1/x) / (1 - 1/x) = (1 + 1/x)(1 + 1/x + 1/x^2 + ...)
    // = 1 + 2/x + 2/x^2 + ...

    let x = Expr::new_variable("x");

    let expr = Expr::new_div(
        Expr::new_add(
            x.clone(),
            Expr::new_constant(1.0),
        ),
        Expr::new_sub(
            x.clone(),
            Expr::new_constant(1.0),
        ),
    );

    let series = asymptotic_expansion(
        &expr,
        "x",
        &Expr::Infinity,
        2,
    );

    // Evaluate at x = 100
    // Exact: 101/99 = 1.020202...
    // Series: 1 + 2/100 + 2/10000 = 1.0202

    let val_at_100 = rssn::symbolic::calculus::substitute(
        &series,
        "x",
        &Expr::new_constant(100.0),
    );

    assert_is_value(&val_at_100, 1.0202);
}
