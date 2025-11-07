// File: tests/numerical\polynomial.rs
//
// Integration tests for the 'rssn' crate's public API in the numerical::polynomial module.
//
// Goal: Ensure the public functions and types in this module behave correctly 
// when used from an external crate context.
//
// --- IMPORTANT FOR NEW CONTRIBUTORS ---
// 1. Standard Tests (`#[test]`): Use these for known inputs and simple assertions.
// 2. Property Tests (`proptest!`): Use these for invariants and edge cases.
//    Proptest runs the test with thousands of generated inputs.

use rssn::numerical::polynomial; 
use proptest::prelude::*; 
use assert_approx_eq::assert_approx_eq; // A useful macro for numerical comparisons

// --- 1. Standard Unit/Integration Tests ---
#[test]
fn test_polynomial_eval_and_derivative_high_precision() {
    // P(x) = 3*x^3 - 2*x^2 + 0.5*x - 7
    // coeffs in descending-power order: [3, -2, 0.5, -7]
    let p = polynomial::Polynomial {
        coeffs: vec![3.0, -2.0, 0.5, -7.0],
    };

    // Evaluation point chosen to avoid trivial cancellation
    let x = 1.23456789f64;

    // Ground truth computed with high precision (Python Decimal at 50 digits)
    // Decimal('-3.78600268967063961730879300000000')
    // Rounded to f64: -3.7860026896706396
    let expected = -3.7860026896706396f64;

    let result = p.eval(x);

    // Tolerance appropriate for f64 Horner evaluation
    let tol = 1e-12f64;
    assert!(
        (result - expected).abs() < tol,
        "Polynomial eval high-precision test failed: result={}, expected={}, diff={}",
        result,
        expected,
        (result - expected).abs()
    );

    // Also test derivative coefficients:
    // dP/dx = 9*x^2 - 4*x + 0.5 -> coefficients [9.0, -4.0, 0.5]
    let d = p.derivative();
    let expected_deriv = vec![9.0, -4.0, 0.5];

    assert_eq!(
        d.coeffs.len(),
        expected_deriv.len(),
        "Derivative coefficient length mismatch"
    );
    for (i, (&got, &want)) in d.coeffs.iter().zip(expected_deriv.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-15,
            "Derivative coefficient mismatch at index {}: got={}, want={}",
            i,
            got,
            want
        );
    }
}

#[test]
fn test_expected_error_behavior() {
    // Example: Test if a function correctly returns an error for invalid input (e.g., division by zero).
    // assert!(numerical::polynomial::divide(1.0, 0.0).is_err());

    let p = polynomial::Polynomial{
        coeffs: vec![1.0,2.0,3.0],
    };

    assert!(p.clone().div_scalar(0.0).is_err(), "Expected error on division by zero");

    let p_div = p.clone().div_scalar(2.0).expect("should divide by non-zero");
    let expected = polynomial::Polynomial{
        coeffs: vec![0.5,1.0,1.5],
    };

    for &x in &[0.0_f64, 1.0_f64, -2.5_f64]{
        assert_approx_eq!(p_div.eval(x), expected.eval(x), 1e-12)
    }

}


// --- 2. Property-Based Tests (Proptest) ---
proptest! {
    #[test]
    fn prop_add_sub_inverse(
        coeffs_p in prop::collection::vec(-1e6f64..1e6f64,1..6),
        coeffs_q in prop::collection::vec(-1e6f64..1e6f64,1..6),
    ) {
        
        prop_assume!(coeffs_p.iter().all(|v| v.is_finite()));
        prop_assume!(coeffs_q.iter().all(|v| v.is_finite()));

        let p = polynomial::Polynomial { coeffs: coeffs_p};
        let q = polynomial::Polynomial { coeffs: coeffs_q};

        let sum = p.clone() + q.clone();
        let recovered = sum -q.clone();

        for &x in &[0.0_f64, 0.3141592653_f64, -1.2345_f64] {
            assert_approx_eq!(p.eval(x),recovered.eval(x),1e-9);
        }

    }

    #[test]
    fn prop_scalar_mul_div_inverse(
        coeffs in prop::collection::vec(-1e6f64..1e6f64, 1..6),
        s in -1e3f64..1e3f64,
    ){

        prop_assume!(coeffs.iter().all(|v| v.is_finite()));
        prop_assume!(s.is_finite() &&s.abs() > 1e-12 );

        let p = polynomial::Polynomial{coeffs:coeffs};
        let scaled = p.clone() * s;

        let recovered = scaled.div_scalar(s).except("s non-zero by assume");

        let &x in &[0.0_f64, 0.6180339887_f64, -3.14159_f64] {
            assert_approx_eq!(p.eval(x),recovered.eval(x),1e-9);
        }

    }

    fn prop_derivative_matches_finite_difference(
        coeffs in prop::collection::vec(-1e3f64..1e3f64, 1..6),
        x in -10.0f64..10.0f64,

    ){
        prop_assume!(coeffs.iter().all(|v| v.is_finite()));
        prop_assume!(x.is_finite());

        let p = polynomial::Polynomial{coeffs:coeffs};

        prop_assume!(p.coeff.len()>=2);

        let d = p.derivative();

        let h = 1e-6f64;
        let numeric = (p.eval(x+h) - p.eval(x-h))/(2.0*h);
        let analytic = d.eval(x);

        assert_approx!(analytic,numeric,1e-6);
    }

}
