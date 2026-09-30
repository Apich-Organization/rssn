use rssn::compute::{compute, d, ComputeConfig};
use rssn::symbolic::core::Expr;
use num_bigint::BigInt;

#[test]
fn test_compute_calculus_symbolic_identity() {
    let x = Expr::new_variable("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(2)));

    // API: let config = default.with_calculus(); compute(d(...)/dx, config)
    let config = ComputeConfig::default()
        .with_calculus()
        .target_symbolic();

    let result = compute(&d(x_sq, "x"), &config);

    // Should reduce d(x^2)/dx to 2*x (or equivalent product)
    eprintln!("Symbolic derivative result: {:?}", result);
    assert!(matches!(result, Expr::Mul(_, _)));
}

#[test]
fn test_compute_calculus_numerical_identity_with_binding() {
    let x = Expr::new_variable("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_bigint(BigInt::from(2)));

    // API: compute with numerical binding
    let config = ComputeConfig::default()
        .with_calculus()
        .bind("x", 3.0)
        .target_numerical(1e-7);

    let result = compute(&d(x_sq, "x"), &config);

    eprintln!("Numerical derivative at x=3: {:?}", result);
    match result {
        Expr::Constant(val) => assert!((val - 6.0).abs() < 1e-6),
        Expr::BigInt(val) => assert_eq!(val, BigInt::from(6)),
        other => panic!("Expected numerical constant, got: {:?}", other),
    }
}
