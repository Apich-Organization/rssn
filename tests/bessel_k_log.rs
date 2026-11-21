//! Arbitrary real-order K in the tails and at the centered quadrature peak.

use std::f64::consts::LN_2;
use std::f64::consts::PI;

use rssn::kernels::special::bessel_k;
use rssn::kernels::special::bessel_k0;
use rssn::kernels::special::bessel_k1;
use rssn::kernels::special::ln_bessel_k;

#[test]
fn large_order_and_argument_do_not_overflow_the_integrand() {
    // mpmath at 65 digits. The original two factors of the integrand
    // could both leave the representable range despite a finite product.
    for (nu, x, expected_value, expected_log) in [
        (
            800.0,
            400.0,
            5.595_652_515_000_873e111,
            257.308_935_281_967_87,
        ),
        (
            1000.0,
            1000.0,
            1.298_180_251_466_700_9e-233,
            -536.241_363_190_353_7,
        ),
        (
            0.5,
            700.0,
            4.670_609_799_936_133_5e-306,
            -703.049_748_814_877,
        ),
        (
            0.5,
            1e-300,
            1.253_314_137_315_500_2e150,
            345.613_555_301_751_6,
        ),
        (0.0, 1e-130, 299.451_993_604_884_35, 5.701_954_116_243_442),
        (
            0.3,
            f64::from_bits(1),
            1.807_351_518_830_335_4e97,
            223.942_616_544_840_5,
        ),
        (0.0001, 1e-300, 691.441_236_153_746_4, 6.538_778_167_238_901),
    ] {
        let value = bessel_k(nu, x);
        let logarithm = ln_bessel_k(nu, x);
        assert!(
            (value / expected_value - 1.0).abs() < 4e-13,
            "{nu},{x}: {value}"
        );
        assert!(
            (logarithm - expected_log).abs() < 4e-13,
            "{nu},{x}: {logarithm}"
        );
        assert_eq!(value, bessel_k(-nu, x));
        assert_eq!(logarithm, ln_bessel_k(-nu, x));
    }
}

#[test]
fn half_order_and_zero_order_have_independent_limits() {
    let mut probes = vec![
        f64::from_bits(1),
        1e-300,
        0.25,
        0.8,
        1.0,
        16.0,
        16.0_f64.next_up(),
        740.0,
        1000.0,
    ];
    probes.extend([700.0_f64, 745.0_f64, f64::MAX]);
    for x in probes {
        let reference = 0.5 * ((PI / 2.0).ln() - x.ln()) - x;
        let logarithm = ln_bessel_k(0.5, x);
        assert!((logarithm - reference).abs() < 3e-13, "x={x}: {logarithm}");
        if x >= 740.0 {
            let expected = reference.exp().to_bits();
            let actual = bessel_k(0.5, x).to_bits();
            assert!(actual.abs_diff(expected) <= 1);
        }
    }
    // Corrections to K0(x) ~ -ln(x/2) - EulerGamma vanish at these inputs.
    for x in [f64::from_bits(1), 1e-300, 1e-130] {
        let expected = LN_2 - x.ln() - 0.577_215_664_901_532_9;
        let actual = bessel_k0(x);
        assert!((actual / expected - 1.0).abs() < 5e-14);
    }
    let tiny = 1e-300_f64;
    // K3(x) ~ 2² Gamma(3) x^-3 = 8 x^-3.
    assert_eq!(bessel_k(3.0, tiny), f64::INFINITY);
    assert!((ln_bessel_k(3.0, tiny) - (8.0_f64.ln() - 3.0 * tiny.ln())).abs() < 8e-13);
}

#[test]
fn recurrence_and_nonfinite_limits_keep_the_real_domain() {
    let nu = 800.0;
    let x = 600.0;
    let from_recurrence = bessel_k(nu - 1.0, x) + 2.0 * nu / x * bessel_k(nu, x);
    assert!((bessel_k(nu + 1.0, x) / from_recurrence - 1.0).abs() < 4e-13);
    assert_eq!(bessel_k0(0.0), f64::INFINITY);
    assert_eq!(bessel_k1(f64::INFINITY), 0.0);
    assert_eq!(ln_bessel_k(1.0, f64::INFINITY), f64::NEG_INFINITY);
    assert_eq!(ln_bessel_k(f64::INFINITY, 2.0), f64::INFINITY);
    assert!(ln_bessel_k(f64::INFINITY, f64::INFINITY).is_nan());
    assert!(bessel_k(0.0, -2.0).is_nan());
    assert!(ln_bessel_k(f64::NAN, 1.0).is_nan());
}
