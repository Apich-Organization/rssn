//! Boundary and limiting arguments of the special functions: no NaN where the function has a
//! value or a limit, no hang at infinity or at huge arguments.

use rssn::kernels::special::airy_ai;
use rssn::kernels::special::airy_ai_prime;
use rssn::kernels::special::airy_bi;
use rssn::kernels::special::airy_bi_prime;
use rssn::kernels::special::bessel_i0;
use rssn::kernels::special::bessel_i1;
use rssn::kernels::special::cosh_integral;
use rssn::kernels::special::cosine_integral;
use rssn::kernels::special::exp_integral_e1;
use rssn::kernels::special::exp_integral_ei;
use rssn::kernels::special::ln_gamma_numerical;
use rssn::kernels::special::log_integral;
use rssn::kernels::special::sinc;
use rssn::kernels::special::sine_integral;
use rssn::kernels::special::sinh_integral;
use rssn::kernels::special_ext::elliptic_e;
use std::f64::consts::PI;

#[test]
fn log_gamma_of_negative_arguments_is_ln_abs_gamma() {
    let close = |a: f64, b: f64| (a - b).abs() < 1e-12 * (1.0 + b.abs());
    assert!(close(ln_gamma_numerical(-0.5), (2.0 * PI.sqrt()).ln()));
    assert!(close(
        ln_gamma_numerical(-1.5),
        (4.0 * PI.sqrt() / 3.0).ln()
    ));
    assert!(close(
        ln_gamma_numerical(-2.5),
        (-(8.0 * PI.sqrt() / 15.0)).abs().ln()
    ));
    // near the pole at 0 from the left: |Gamma(-e)| ~ 1/e
    assert!(close(
        ln_gamma_numerical(-1e-12),
        27.631_021_115_928_55 + 5.772_156_649e-12
    ));
    // reflection consistency: ln|G(x)| + ln|G(1-x)| = ln(pi / |sin(pi x)|)
    for x in [-7.3, -0.2, -12.9] {
        let lhs = ln_gamma_numerical(x) + ln_gamma_numerical(1.0 - x);
        assert!(close(lhs, (PI / (PI * x).sin().abs()).ln()), "{x}");
    }
    for pole in [-1.0, -2.0, -30.0] {
        assert_eq!(ln_gamma_numerical(pole), f64::INFINITY);
    }
    // positive arguments are unchanged
    assert!(close(ln_gamma_numerical(5.0), 24.0_f64.ln()));
}

#[test]
fn infinite_arguments_give_the_limits() {
    let inf = f64::INFINITY;
    assert_eq!(bessel_i0(inf), inf);
    assert_eq!(bessel_i0(-inf), inf);
    assert_eq!(bessel_i1(inf), inf);
    assert_eq!(bessel_i1(-inf), -inf);
    assert_eq!(sinc(inf), 0.0);
    assert_eq!(sinc(-inf), 0.0);
    assert_eq!(exp_integral_e1(0.0), inf);
    assert_eq!(exp_integral_e1(inf), 0.0);
    assert_eq!(exp_integral_ei(inf), inf);
    assert_eq!(exp_integral_ei(-inf), 0.0);
    assert_eq!(log_integral(inf), inf);
    assert!((sine_integral(inf) - PI / 2.0).abs() < 1e-15);
    assert!((sine_integral(-inf) + PI / 2.0).abs() < 1e-15);
    assert_eq!(cosine_integral(0.0), -inf);
    assert_eq!(cosine_integral(inf), 0.0);
    assert_eq!(sinh_integral(inf), inf);
    assert_eq!(sinh_integral(-inf), -inf);
    assert_eq!(cosh_integral(0.0), -inf);
    assert_eq!(cosh_integral(inf), inf);
}

#[test]
fn airy_functions_at_large_arguments_are_finite_and_fast() {
    // the Wronskian Ai Bi' - Ai' Bi = 1/pi on both sides of the switch to the asymptotic series
    for x in [
        -9.9, -10.0, -12.5, -30.0, -100.0, -1e3, -1e5, -1e8, -0.3, 0.7, 5.0, 20.0,
    ] {
        let w = airy_ai(x) * airy_bi_prime(x) - airy_ai_prime(x) * airy_bi(x);
        // the oscillatory functions grow like x^(1/4) in the derivatives: scale the tolerance
        assert!((w * PI - 1.0).abs() < 1e-6, "x = {x}: W pi = {}", w * PI);
    }
    // reference values
    assert!((airy_ai(-10.0) - 0.040_241_238_486_443_19).abs() < 1e-12);
    assert!(
        (airy_bi(-10.0) + 0.314_679_829_643_838).abs() < 1e-12,
        "{}",
        airy_bi(-10.0)
    );
    // the derivative formulas agree with differences of the values
    for x in [-15.0, -50.0, -400.0] {
        let h = 1e-5;
        let (da, db) = (
            (airy_ai(x + h) - airy_ai(x - h)) / (2.0 * h),
            (airy_bi(x + h) - airy_bi(x - h)) / (2.0 * h),
        );
        assert!(
            (da - airy_ai_prime(x)).abs() < 1e-6 * (1.0 + da.abs()),
            "Ai' at {x}"
        );
        assert!(
            (db - airy_bi_prime(x)).abs() < 1e-6 * (1.0 + db.abs()),
            "Bi' at {x}"
        );
    }
    // limits and overflow, without long loops
    assert_eq!(airy_ai(f64::INFINITY), 0.0);
    assert_eq!(airy_ai(f64::NEG_INFINITY), 0.0);
    assert_eq!(airy_bi(f64::NEG_INFINITY), 0.0);
    assert_eq!(airy_bi(f64::INFINITY), f64::INFINITY);
    assert_eq!(airy_bi(1e10), f64::INFINITY);
    assert_eq!(airy_bi_prime(1e10), f64::INFINITY);
    assert!(airy_ai_prime(f64::NEG_INFINITY).is_nan());
    assert!(airy_ai(f64::NAN).is_nan() && airy_bi(f64::NAN).is_nan());
    assert!(airy_bi(100.0).is_finite() && airy_bi(100.0) > 1e100);
    assert!(airy_ai(-1e5).is_finite() && airy_ai(-1e5).abs() < 0.1);
}

#[test]
fn complete_elliptic_e_for_large_negative_parameters() {
    let direct = |m: f64| {
        // Simpson's rule for the defining integral
        let n = 4000;
        let h = PI / 2.0 / f64::from(n);
        let f = |t: f64| (1.0 - m * t.sin().powi(2)).sqrt();
        let mut sum = f(0.0) + f(PI / 2.0);
        for k in 1..n {
            sum += f(f64::from(k) * h) * if k % 2 == 1 { 4.0 } else { 2.0 };
        }
        sum * h / 3.0
    };
    for m in [-1.5, -3.0, -40.0, -1e4] {
        let (v, want) = (elliptic_e(m), direct(m));
        assert!((v - want).abs() < 1e-9 * want, "m = {m}: {v} vs {want}");
    }
    assert_eq!(elliptic_e(f64::NEG_INFINITY), f64::INFINITY);
    assert!((elliptic_e(0.0) - PI / 2.0).abs() < 1e-15);
    assert_eq!(elliptic_e(1.0), 1.0);
}
