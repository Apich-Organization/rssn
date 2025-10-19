// File: tests/symbolic\simplify.rs
//
// Integration tests for the 'rssn' crate's public API in the symbolic::simplify module.
//
// Goal: Ensure the public functions and types in this module behave correctly
// when used from an external crate context.
//
// --- IMPORTANT FOR NEW CONTRIBUTORS ---
// 1. Standard Tests (`#[test]`): Use these for known inputs and simple assertions.
// 2. Property Tests (`proptest!`): Use these for invariants and edge cases.
//    Proptest runs the test with thousands of generated inputs.

use assert_approx_eq::assert_approx_eq; // A useful macro for numerical comparisons
use num_bigint::BigInt;
use proptest::prelude::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify;

// --- 1. Standard Unit/Integration Tests ---
#[test]
fn test_initial_conditions_or_edge_cases() {
    // Example: Test a function with input '0' or large, known values.
    // let result = symbolic::simplify::some_function(42.0);
    // assert_approx_eq!(result, 1.0, 1e-6);
}

#[test]
fn test_expected_error_behavior() {
    // Example: Test if a function correctly returns an error for invalid input (e.g., division by zero).
    // assert!(symbolic::simplify::divide(1.0, 0.0).is_err());
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
        // let val = symbolic::simplify::add(a, b);
        // assert_approx_eq!(symbolic::simplify::subtract(val, b), a, 1e-9);

        // INVARIANT 2: Test basic property (e.g., matrix transpose twice is the original)
        // let matrix = symbolic::simplify::create_random_matrix();
        // assert_eq!(matrix.transpose().transpose(), matrix);
    }
}

#[test]
fn test_simplify_add_x_x() {
    let x = Expr::new_variable("x");
    let expr = Expr::new_add(x.clone(), x.clone());
    let simplified = simplify::simplify(expr);
    let two = Expr::new_bigint(BigInt::from(2));
    let expected = Expr::new_mul(two, x);
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_add_2x_3x() {
    let x = Expr::new_variable("x");
    let two = Expr::new_bigint(BigInt::from(2));
    let three = Expr::new_bigint(BigInt::from(3));
    let five = Expr::new_bigint(BigInt::from(5));
    let expr = Expr::new_add(
        Expr::new_mul(two, x.clone()),
        Expr::new_mul(three, x.clone()),
    );
    let simplified = simplify::simplify(expr);
    let expected = Expr::new_mul(five, x);
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_nested_add() {
    let x = Expr::new_variable("x");
    let one = Expr::new_bigint(BigInt::from(1));
    let two = Expr::new_bigint(BigInt::from(2));
    let three = Expr::new_bigint(BigInt::from(3));
    let expr = Expr::new_add(Expr::new_add(x.clone(), one), Expr::new_add(x.clone(), two));
    let simplified = simplify::simplify(expr);
    let expected = Expr::new_add(Expr::new_mul(Expr::new_bigint(BigInt::from(2)), x), three);
    assert_eq!(simplified, expected);
}

#[test]
fn test_simplify_constants() {
    let one = Expr::new_bigint(BigInt::from(1));
    let two = Expr::new_bigint(BigInt::from(2));
    let three = Expr::new_bigint(BigInt::from(3));
    let expr = Expr::new_add(one, two);
    let simplified = simplify::simplify(expr);
    assert_eq!(simplified, three);
}
