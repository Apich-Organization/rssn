use num_complex::Complex;
use rssn::numerical::multi_valued::*;
use rssn::symbolic::core::Expr;
use std::f64::consts::PI;

#[test]

fn test_newton_method_complex_simple_roots() {

    let z = Expr::Variable("z".to_string());

    let f = Expr::new_sub(
        Expr::new_pow(z.clone(), Expr::Constant(2.0)),
        Expr::Constant(1.0),
    );

    let f_prime = Expr::new_mul(Expr::Constant(2.0), z.clone());

    let root1 = newton_method_complex(&f, &f_prime, Complex::new(2.0, 0.0), 1e-6, 100).unwrap();

    assert!((root1.re - 1.0).abs() < 1e-5);

    let root2 = newton_method_complex(&f, &f_prime, Complex::new(-2.0, 0.0), 1e-6, 100).unwrap();

    assert!((root2.re + 1.0).abs() < 1e-5);
}

#[test]

fn test_complex_log_k() {

    // log(1) = 0 (k=0)
    let z = Complex::new(1.0, 0.0);

    let log0 = complex_log_k(z, 0);

    assert!(log0.re.abs() < 1e-9);

    assert!(log0.im.abs() < 1e-9);

    // log(1) + 2*pi*i (k=1)
    let log1 = complex_log_k(z, 1);

    assert!(log1.re.abs() < 1e-9);

    assert!((log1.im - 2.0 * PI).abs() < 1e-9);
}

#[test]

fn test_complex_sqrt_k() {

    // sqrt(1) = 1 (k=0)
    let z = Complex::new(1.0, 0.0);

    let sqrt0 = complex_sqrt_k(z, 0);

    assert!((sqrt0.re - 1.0).abs() < 1e-9);

    assert!(sqrt0.im.abs() < 1e-9);

    // sqrt(1) = -1 (k=1)
    let sqrt1 = complex_sqrt_k(z, 1);

    assert!((sqrt1.re + 1.0).abs() < 1e-9);

    assert!(sqrt1.im.abs() < 1e-9);

    // sqrt(-1) = i (k=0)
    let z_neg = Complex::new(-1.0, 0.0);

    let sqrt_neg0 = complex_sqrt_k(z_neg, 0);

    assert!(sqrt_neg0.re.abs() < 1e-9);

    assert!((sqrt_neg0.im - 1.0).abs() < 1e-9);
}

#[cfg(test)]

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_log_k_exp(re in -2.0..2.0f64, im in -2.0..2.0f64, k in -5..5i32) {
            let z = Complex::new(re, im);
            if z.norm() > 1e-6 {
                let log_z_k = complex_log_k(z, k);
                let exp_log = log_z_k.exp();
                // exp(log(z) + 2kpi*i) = exp(log(z)) * exp(2kpi*i) = z * 1 = z
                prop_assert!((exp_log - z).norm() < 1e-7);
            }
        }

        #[test]
        fn prop_sqrt_k_sqr(re in -2.0..2.0f64, im in -2.0..2.0f64, k in 0..2i32) {
             let z = Complex::new(re, im);
             let sqrt_z = complex_sqrt_k(z, k);
             let sqr = sqrt_z * sqrt_z;
             prop_assert!((sqr - z).norm() < 1e-7);
        }
    }
}
