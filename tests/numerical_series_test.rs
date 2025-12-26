use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::numerical::series::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_taylor_coefficients() {

    let x = Expr::new_variable("x");

    let f = Expr::new_pow(
        x,
        Expr::new_constant(2.0),
    ); // f(x) = x^2
    let coeffs = taylor_coefficients(
        &f, "x", 0.0, 2,
    )
    .unwrap();

    // coeffs = [0, 0, 1]
    assert_approx_eq!(
        coeffs[0], 0.0, 1e-10f64
    );

    assert_approx_eq!(
        coeffs[1], 0.0, 1e-10f64
    );

    assert_approx_eq!(
        coeffs[2], 1.0, 1e-10f64
    );
}

#[test]

fn test_evaluate_power_series() {

    let coeffs = vec![1.0, 1.0, 0.5]; // 1 + x + x^2/2
    let val = evaluate_power_series(
        &coeffs, 0.0, 1.0,
    );

    assert_approx_eq!(
        val, 2.5, 1e-10f64
    );
}

#[test]

fn test_sum_series() {

    let n = Expr::new_variable("n");

    let f = n;

    let sum =
        sum_series(&f, "n", 1, 10)
            .unwrap();

    assert_approx_eq!(
        sum, 55.0, 1e-10f64
    );
}

proptest! {
    #[test]
    fn proptest_sum_linear(n_max in 1..100i64) {
        let n = Expr::new_variable("n");
        let f = n;
        let sum = sum_series(&f, "n", 1, n_max).unwrap();
        let expected = (n_max * (n_max + 1)) as f64 / 2.0;
        assert_approx_eq!(sum, expected, 1e-9f64);
    }
}
