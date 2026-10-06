//! Reference-value tests for special functions (new and existing).

use rssn::kernels::special::{
    erf_numerical, erfc_numerical, gamma_numerical, inverse_erf_numerical, ln_gamma_numerical,
    digamma_numerical,
};
use rssn::kernels::special_ext::{
    bessel_in, bessel_jn, carlson_rf, elliptic_e, elliptic_e_inc, elliptic_f, elliptic_k,
    hyp1f1, hyp2f1, jacobi_elliptic,
};

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} vs {b}");
}

#[test]
fn complete_elliptic() {
    close(elliptic_k(0.5), 1.854_074_677_301_372, 1e-14);
    close(elliptic_e(0.5), 1.350_643_881_047_675, 1e-14);
    close(elliptic_k(0.0), std::f64::consts::FRAC_PI_2, 1e-15);
    close(elliptic_e(0.0), std::f64::consts::FRAC_PI_2, 1e-15);
    close(elliptic_k(-1.0), 1.311_028_777_146_059_9, 1e-14);
    close(elliptic_k(0.9), 2.578_092_113_348_173, 1e-14);
    close(elliptic_e(0.9), 1.104_774_732_704_073, 1e-13);
    assert!(elliptic_k(1.0).is_infinite());
    assert!((elliptic_e(1.0) - 1.0).abs() < 1e-15);
    assert!(elliptic_k(2.0).is_nan());
    // Legendre relation: E K' + E' K - K K' = pi/2
    let m = 0.3;
    let r = elliptic_e(m) * elliptic_k(1.0 - m) + elliptic_e(1.0 - m) * elliptic_k(m)
        - elliptic_k(m) * elliptic_k(1.0 - m);
    close(r, std::f64::consts::FRAC_PI_2, 1e-13);
}

#[test]
fn incomplete_elliptic() {
    // F(pi/2, m) = K(m), E(pi/2, m) = E(m)
    close(elliptic_f(std::f64::consts::FRAC_PI_2, 0.5), elliptic_k(0.5), 1e-13);
    close(elliptic_e_inc(std::f64::consts::FRAC_PI_2, 0.5), elliptic_e(0.5), 1e-13);
    // F(phi, 0) = phi
    close(elliptic_f(0.7, 0.0), 0.7, 1e-14);
    // reference: F(pi/4, 0.5) = 0.826017876...
    close(elliptic_f(std::f64::consts::FRAC_PI_4, 0.5), 0.826_017_876_249_245, 1e-12);
    // periodicity
    close(elliptic_f(std::f64::consts::PI + 0.3, 0.4), elliptic_f(0.3, 0.4) + 2.0 * elliptic_k(0.4), 1e-13);
    close(carlson_rf(1.0, 2.0, 3.0), 0.726_945_935_468_9, 1e-11);
}

#[test]
fn jacobi_functions() {
    let (sn, cn, dn) = jacobi_elliptic(0.7, 0.0);
    close(sn, 0.7_f64.sin(), 1e-14);
    close(cn, 0.7_f64.cos(), 1e-14);
    close(dn, 1.0, 1e-14);
    let (sn, cn, dn) = jacobi_elliptic(1.0, 0.5);
    close(sn, 0.803_001_824_895_643_9, 1e-9);
    close(cn, (1.0 - sn * sn).sqrt(), 1e-12);
    close(dn, (1.0 - 0.5 * sn * sn).sqrt(), 1e-12);
    close(sn * sn + cn * cn, 1.0, 1e-14);
}

#[test]
fn hypergeometric() {
    // 2F1(1,1;2;z) = -ln(1-z)/z
    for z in [-3.0, -0.4, 0.3, 0.8, 0.97] {
        close(hyp2f1(1.0, 1.0, 2.0, z), -(1.0 - z as f64).ln() / z, 1e-11);
    }
    // 2F1(a,b;b;z) = (1-z)^-a
    close(hyp2f1(0.5, 2.0, 2.0, 0.6), 0.4_f64.powf(-0.5), 1e-13);
    // 2F1(1/2,1/2;3/2;z^2) = asin(z)/z
    close(hyp2f1(0.5, 0.5, 1.5, 0.25), 0.5_f64.asin() / 0.5, 1e-13);
    assert!(hyp2f1(1.0, 1.0, 2.0, 1.0).is_nan());
    // 1F1(1;2;x) = (e^x - 1)/x ; 1F1(a;a;x) = e^x
    close(hyp1f1(1.0, 2.0, 3.0), (3.0_f64.exp() - 1.0) / 3.0, 1e-13);
    close(hyp1f1(1.0, 2.0, -3.0), ((-3.0_f64).exp() - 1.0) / -3.0, 1e-13);
    close(hyp1f1(2.5, 2.5, 4.0), 4.0_f64.exp(), 1e-13);
    // Kummer's second example: 1F1(1/2;1;x) relation with Bessel: e^{-x/2} I0(x/2)
    close(hyp1f1(0.5, 1.0, 2.0), (1.0_f64).exp() * bessel_in(0, 1.0), 1e-13);
}

#[test]
fn bessel_integer_order() {
    close(bessel_jn(0, 1.0), 0.765_197_686_557_966_6, 1e-15);
    close(bessel_jn(1, 1.0), 0.440_050_585_744_933_5, 1e-15);
    close(bessel_jn(2, 5.0), 0.046_565_116_277_752_2, 1e-14);
    close(bessel_jn(0, 10.0), -0.245_935_764_451_348_3, 1e-14);
    close(bessel_jn(0, 50.0), 0.055_812_327_669_251_9, 1e-13);
    close(bessel_jn(-1, 2.0), -bessel_jn(1, 2.0), 1e-15);
    close(bessel_in(0, 1.0), 1.266_065_877_752_008_4, 1e-14);
    close(bessel_in(1, 2.0), 1.590_636_854_637_329, 1e-14);
    // recurrence J_{n-1} + J_{n+1} = 2n/x J_n
    close(bessel_jn(2, 3.0) + bessel_jn(4, 3.0), 6.0 / 3.0 * bessel_jn(3, 3.0), 1e-14);
}

#[test]
fn existing_gamma_erf_reference_values() {
    close(gamma_numerical(5.0), 24.0, 1e-13);
    close(gamma_numerical(0.5), std::f64::consts::PI.sqrt(), 1e-13);
    close(gamma_numerical(-1.5), 2.363_271_801_207_355, 1e-12);
    close(ln_gamma_numerical(100.0), 359.134_205_369_575_4, 1e-13);
    close(digamma_numerical(1.0), -0.577_215_664_901_532_9, 1e-12);
    close(digamma_numerical(0.5), -1.963_510_026_021_423_5, 1e-12);
    close(erf_numerical(0.5), 0.520_499_877_813_046_5, 1e-14);
    close(erf_numerical(2.0), 0.995_322_265_018_952_7, 1e-14);
    close(erfc_numerical(3.0), 2.209_049_699_858_544e-5, 1e-9);
    close(inverse_erf_numerical(0.5), 0.476_936_276_204_469_9, 1e-12);
    close(erf_numerical(inverse_erf_numerical(-0.9)), -0.9, 1e-13);
}
