use std::collections::HashMap;

use num_complex::Complex;
use rssn::numerical::complex_analysis::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_contour_integral_circle() {

    // ∮ 1/z dz around unit circle = 2πi
    let f = |z : Complex<f64>| 1.0 / z;

    let mut path = Vec::new();

    let n = 1000;

    for i in 0 ..= n {

        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);

        path.push(Complex::new(
            angle.cos(),
            angle.sin(),
        ));
    }

    let res = contour_integral(f, &path);

    println!(
        "Circle integral res: {:?}",
        res
    );

    assert!((res.re - 0.0).abs() < 1e-7);

    assert!((res.im - 2.0 * std::f64::consts::PI).abs() < 1e-7);
}

#[test]

fn test_residue_simple_pole() {

    // f(z) = 1/z has residue 1 at z=0
    let f = |z : Complex<f64>| 1.0 / z;

    let res = residue(
        f,
        Complex::new(0.0, 0.0),
        0.1,
        1000,
    );

    println!(
        "Simple pole residue: {:?}",
        res
    );

    assert!((res.re - 1.0).abs() < 1e-7);

    assert!((res.im - 0.0).abs() < 1e-7);
}

#[test]

fn test_residue_double_pole() {

    // f(z) = 1/z^2 has residue 0 at z=0
    let f = |z : Complex<f64>| 1.0 / (z * z);

    let res = residue(
        f,
        Complex::new(0.0, 0.0),
        0.1,
        1000,
    );

    println!(
        "Double pole residue: {:?}",
        res
    );

    assert!(res.norm() < 1e-7);
}

#[test]

fn test_eval_complex_expr() {

    let mut vars = HashMap::new();

    vars.insert(
        "z".to_string(),
        Complex::new(1.0, 1.0),
    );

    // expr = z^2 + 1
    let z = Expr::Variable("z".to_string());

    let expr = Expr::new_add(
        Expr::new_pow(
            z,
            Expr::Constant(2.0),
        ),
        Expr::Constant(1.0),
    );

    let res = eval_complex_expr(&expr, &vars).unwrap();

    println!(
        "Eval res: {:?}",
        res
    );

    assert!((res.re - 1.0).abs() < 1e-9);

    assert!((res.im - 2.0).abs() < 1e-9);
}

#[test]

fn test_moebius_transform() {

    // f(z) = (z - 1) / (z + 1)
    let m = MobiusTransformation::new(
        Complex::new(1.0, 0.0),
        Complex::new(-1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
    );

    assert_eq!(
        m.apply(Complex::new(
            1.0, 0.0
        )),
        Complex::new(0.0, 0.0)
    );

    assert_eq!(
        m.apply(Complex::new(
            0.0, 0.0
        )),
        Complex::new(-1.0, 0.0)
    );

    let inv = m.inverse();

    let res = inv.apply(
        m.apply(Complex::new(
            2.0, 3.0,
        )),
    );

    assert!((res.re - 2.0).abs() < 1e-9);

    assert!((res.im - 3.0).abs() < 1e-9);
}

#[cfg(test)]

mod proptests {

    use proptest::prelude::*;

    use super::*;

    proptest! {
        #[test]
        fn prop_moebius_inversion(
            a_re in -10.0..10.0, a_im in -10.0..10.0,
            b_re in -10.0..10.0, b_im in -10.0..10.0,
            c_re in -10.0..10.0, c_im in -10.0..10.0,
            d_re in -10.0..10.0, d_im in -10.0..10.0,
            z_re in -10.0..10.0, z_im in -10.0..10.0
        ) {
            let m = MobiusTransformation::new(
                Complex::new(a_re, a_im), Complex::new(b_re, b_im),
                Complex::new(c_re, c_im), Complex::new(d_re, d_im)
            );
            let z = Complex::new(z_re, z_im);

            // Avoid singularities (denominator close to zero)
            if (m.c * z + m.d).norm() > 1e-3 {
                let w = m.apply(z);
                let inv = m.inverse();
                let z_back = inv.apply(w);
                prop_assert!((z - z_back).norm() < 1e-7);
            }
        }

        #[test]
        fn prop_eval_add_symmetric(
            a_re in -100.0..100.0, a_im in -100.0..100.0,
            b_re in -100.0..100.0, b_im in -100.0..100.0
        ) {
            let mut vars = HashMap::new();
            vars.insert("a".to_string(), Complex::new(a_re, a_im));
            vars.insert("b".to_string(), Complex::new(b_re, b_im));

            let expr1 = Expr::new_add(Expr::Variable("a".to_string()), Expr::Variable("b".to_string()));
            let expr2 = Expr::new_add(Expr::Variable("b".to_string()), Expr::Variable("a".to_string()));

            let res1 = eval_complex_expr(&expr1, &vars).unwrap();
            let res2 = eval_complex_expr(&expr2, &vars).unwrap();
            prop_assert!((res1 - res2).norm() < 1e-9);
        }
    }
}
