// File: tests/symbolic\calculus.rs
//
// Integration tests for the 'rssn' crate's public API in the symbolic::calculus module.
//
// Goal: Ensure the public functions and types in this module behave correctly
// when used from an external crate context.
//
// --- IMPORTANT FOR NEW CONTRIBUTORS ---
// 1. Standard Tests (`#[test]`): Use these for known inputs and simple assertions.
// 2. Property Tests (`proptest!`): Use these for invariants and edge cases.
//    Proptest runs the test with thousands of generated inputs.

use assert_approx_eq::assert_approx_eq; /* A useful macro for numerical comparisons */
use num_bigint::BigInt;
use proptest::prelude::*;
use rssn::symbolic::calculus;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify;
use rssn::symbolic::simplify_dag;

// --- 1. Standard Unit/Integration Tests ---
#[test]

fn test_initial_conditions_or_edge_cases()
 {
    // Example: Test a function with input '0' or large, known values.
    // let result = symbolic::calculus::some_function(42.0);
    // assert_approx_eq!(result, 1.0, 1e-6);
}

#[test]

fn test_factorial_overflow() {

    assert_eq!(
        calculus::factorial(170),
        7.257415615307994e306
    );

    assert!(
        calculus::factorial(171)
            .is_infinite()
    );
}

#[test]

fn test_expected_error_behavior() {
    // Example: Test if a function correctly returns an error for invalid input (e.g., division by zero).
    // assert!(symbolic::calculus::divide(1.0, 0.0).is_err());
}

// --- 2. Property-Based Tests (Proptest) ---
proptest! {
    #[test]
    fn prop_test_invariants_hold(
        // Define inputs using strategies (e.g., f64 in a specific range)
        a in any::<f64>(),
        b in -100.0..100.0f64,
    ) {
        // INVARIANT 1: Test an operation and its inverse
        // let val = symbolic::calculus::add(a, b);
        // assert_approx_eq!(symbolic::calculus::subtract(val, b), a, 1e-9);

        // INVARIANT 2: Test basic property (e.g., matrix transpose twice is the original)
        // let matrix = symbolic::calculus::create_random_matrix();
        // assert_eq!(matrix.transpose().transpose(), matrix);
    }
}

#[test]

fn test_differentiate_sin_x() {

    let x = Expr::new_variable("x");

    let sin_x =
        Expr::new_sin(x.clone());

    let derivative =
        calculus::differentiate(
            &sin_x,
            "x",
        );

    let expected = Expr::new_cos(x);

    assert_eq!(derivative, expected);
}

#[test]

fn test_differentiate_cos_x() {

    let x = Expr::new_variable("x");

    let cos_x =
        Expr::new_cos(x.clone());

    let derivative =
        calculus::differentiate(
            &cos_x,
            "x",
        );

    let expected =
        Expr::new_neg(Expr::new_sin(x));

    assert_eq!(derivative, expected);
}

#[test]

fn test_differentiate_x_cubed() {

    let x = Expr::new_variable("x");

    let three = Expr::new_bigint(
        BigInt::from(3),
    );

    let two = Expr::new_bigint(
        BigInt::from(2),
    );

    let x_cubed = Expr::new_pow(
        x.clone(),
        three.clone(),
    );

    let derivative =
        simplify_dag::simplify(
            &calculus::differentiate(
                &x_cubed,
                "x",
            ),
        );

    let expected = Expr::new_mul(
        three,
        Expr::new_pow(x, two),
    );

    assert_eq!(derivative, expected);
}

#[test]

fn test_differentiate_product_rule() {

    let x = Expr::new_variable("x");

    let sin_x =
        Expr::new_sin(x.clone());

    let expr = Expr::new_mul(
        x.clone(),
        sin_x.clone(),
    );

    let derivative =
        calculus::differentiate(
            &expr, "x",
        );

    let expected = Expr::new_add(
        sin_x,
        Expr::new_mul(
            &x,
            Expr::new_cos(x.clone()),
        ),
    );

    assert_eq!(derivative, expected);
}
