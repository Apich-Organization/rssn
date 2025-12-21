use rssn::prelude::*;
use rssn::prelude::numerical::*;
use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;

#[test]
fn test_poly_basic() {
    let p = numerical_Polynomial { coeffs: vec![1.0, 2.0, 1.0] }; // x^2 + 2x + 1
    assert_eq!(p.eval(0.0), 1.0);
    assert_eq!(p.eval(1.0), 4.0);
    assert_eq!(p.degree(), 2);
}

#[test]
fn test_poly_arithmetic() {
    let p1 = numerical_Polynomial { coeffs: vec![1.0, 1.0] }; // x + 1
    let p2 = numerical_Polynomial { coeffs: vec![1.0, -1.0] }; // x - 1
    
    let sum = p1.clone() + p2.clone();
    assert_eq!(sum.coeffs, vec![2.0, 0.0]); // 2x
    
    let prod = p1 * p2;
    assert_eq!(prod.coeffs, vec![1.0, 0.0, -1.0]); // x^2 - 1
}

#[test]
fn test_poly_calculus() {
    let p = numerical_Polynomial { coeffs: vec![1.0, 0.0, 0.0] }; // x^2
    let dp = p.derivative();
    assert_eq!(dp.coeffs, vec![2.0, 0.0]); // 2x
    
    let ip = p.integral();
    assert_eq!(ip.coeffs, vec![1.0/3.0, 0.0, 0.0, 0.0]); // (1/3)x^3
}

#[test]
fn test_poly_roots() {
    let p = numerical_Polynomial { coeffs: vec![1.0, 0.0, -1.0] }; // x^2 - 1
    let roots = p.find_roots().unwrap();
    assert_eq!(roots.len(), 2);
    assert_approx_eq!(roots[0].abs(), 1.0, 1e-9);
    assert_approx_eq!(roots[1].abs(), 1.0, 1e-9);
}

proptest! {
    #[test]
    fn prop_poly_eval_sum(coeffs1 in proptest::collection::vec(-100.0..100.0f64, 1..10),
                         coeffs2 in proptest::collection::vec(-100.0..100.0f64, 1..10),
                         x in -10.0..10.0f64) {
        let p1 = numerical_Polynomial { coeffs: coeffs1 };
        let p2 = numerical_Polynomial { coeffs: coeffs2 };
        let sum = p1.clone() + p2.clone();
        assert_approx_eq!(sum.eval(x), p1.eval(x) + p2.eval(x), 1e-4);
    }

    #[test]
    fn prop_poly_eval_mul(coeffs1 in proptest::collection::vec(-10.0..10.0f64, 1..5),
                         coeffs2 in proptest::collection::vec(-10.0..10.0f64, 1..5),
                         x in -2.0..2.0f64) {
        let p1 = numerical_Polynomial { coeffs: coeffs1 };
        let p2 = numerical_Polynomial { coeffs: coeffs2 };
        let prod = p1.clone() * p2.clone();
        assert_approx_eq!(prod.eval(x), p1.eval(x) * p2.eval(x), 1e-4);
    }
}
