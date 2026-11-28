//! Gaussian quantiles at limits of the binary64 location and scale.

use rssn::kernels::distributions::Distribution;

#[test]
fn normal_quantiles_survive_overflow_and_subnormal_scales() {
    let largest = f64::MAX;
    let distribution = Distribution::Normal {
        mu: largest / 4.0,
        sigma: largest,
    };
    assert_eq!(distribution.quantile(0.5).unwrap(), largest / 4.0);

    // References from multiprecision bisection on the Gaussian CDF.
    let result = distribution.quantile(0.51).unwrap();
    assert!(result.to_bits().abs_diff(0x7fd1_9aba_9f47_8691) <= 3);
    let cancelling = Distribution::Normal {
        mu: -largest,
        sigma: largest,
    };
    let result = cancelling.quantile(0.95).unwrap();
    assert!(result.to_bits().abs_diff(0x7fe4_a2a4_12ce_d576) <= 3);

    let smallest = Distribution::Normal {
        mu: 0.0,
        sigma: f64::from_bits(1),
    };
    let quantile = smallest.quantile(6.220_960_574_271_78e-16).unwrap();
    assert_eq!(quantile.to_bits(), (1_u64 << 63) | 8);
}

#[test]
fn gaussian_support_medians_and_probability_validation() {
    let normal = Distribution::Normal {
        mu: -0.0,
        sigma: f64::MAX,
    };
    let lognormal = Distribution::LogNormal {
        mu: std::f64::consts::LN_2,
        sigma: f64::MAX,
    };
    assert_eq!(normal.quantile(0.5).unwrap().to_bits(), (-0.0f64).to_bits());
    assert_eq!(normal.quantile(0.0).unwrap(), f64::NEG_INFINITY);
    assert_eq!(normal.quantile(1.0).unwrap(), f64::INFINITY);
    assert_eq!(lognormal.quantile(0.5).unwrap(), 2.0);
    assert_eq!(lognormal.quantile(-0.0).unwrap(), 0.0);
    assert_eq!(lognormal.quantile(1.0).unwrap(), f64::INFINITY);
    for probability in [f64::NAN, f64::INFINITY, -0.01, 1.01] {
        assert!(normal.quantile(probability).is_err());
        assert!(lognormal.quantile(probability).is_err());
    }
}

#[test]
fn normal_quantiles_resolve_both_representable_probability_tails() {
    let normal = Distribution::Normal { mu: 0.0, sigma: 1.0 };
    let near_one = normal.quantile(1.0 - f64::EPSILON / 2.0).unwrap();
    assert!(near_one.to_bits().abs_diff(0x4020_6b48_528c_ea52) <= 4);

    let near_zero = normal.quantile(f64::from_bits(1)).unwrap();
    let expected = f64::from_bits(0xc043_3bd3_f27f_cd03);
    assert!((near_zero - expected).abs() <= 2.0e-14 * expected.abs());
}
