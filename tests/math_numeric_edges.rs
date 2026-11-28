//! Edge cases of the statistics, quadrature and interpolation kernels: empty
//! and non-finite data, degenerate parameters, singular and infinite
//! integrands, tail accuracy of p-values, and validation of interpolation
//! nodes.

use std::f64::consts::PI;

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::integrate as ig;
use rssn::kernels::interpolate as ip;
use rssn::kernels::special::regularized_beta;
use rssn::kernels::special::regularized_upper_gamma;
use rssn::kernels::stats;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0xED6E5),
        failure_persistence: None,
        cases: 300,
        ..ProptestConfig::default()
    }
}

// ---------------------------------------------------------------- statistics

#[test]
fn empty_and_nan_samples_give_nan_not_panics() {
    let mut empty: Vec<f64> = vec![];
    assert!(stats::median(&mut empty).is_nan());
    assert!(stats::percentile(&mut empty, 50.0).is_nan());
    assert!(stats::min(&mut empty).is_nan());
    assert!(stats::max(&mut empty).is_nan());
    assert!(stats::skewness(&mut empty).is_nan());
    assert!(stats::kurtosis(&mut empty).is_nan());
    assert!(stats::iqr(&mut empty).is_nan());
    assert!(stats::std_dev(&[]).is_nan());
    assert!(stats::range(&[]).is_nan());
    assert!(stats::geometric_mean(&[]).is_nan());
    assert!(stats::harmonic_mean(&[]).is_nan());
    assert!(stats::correlation(&[], &[]).is_nan());
    assert!(stats::covariance(&[1.0], &[1.0]).is_nan());
    assert!(
        stats::covariance(&[1.0, 2.0], &[1.0]).is_nan(),
        "length mismatch"
    );
    assert!(stats::chi_squared_test(&[1.0, 2.0], &[1.0]).0.is_nan());
    assert!(stats::two_sample_t_test(&[], &[]).0.is_nan());
    assert!(stats::welch_t_test(&[1.0], &[2.0]).0.is_nan());
    assert!(stats::simple_linear_regression(&[]).0.is_nan());
    assert!(stats::mode(&[], 2).is_none());
    let mut g: Vec<&mut [f64]> = vec![];
    assert!(stats::one_way_anova(&mut g).0.is_nan());
    let mut with_nan = vec![1.0, f64::NAN, 3.0];
    assert!(stats::median(&mut with_nan).is_nan());
    assert!(stats::percentile(&mut with_nan, 10.0).is_nan());
    assert!(stats::iqr(&mut [1.0, 2.0, f64::NAN, 4.0, 5.0]).is_nan());
    assert!(stats::shannon_entropy(&[0.5, f64::NAN]).is_nan());
}

#[test]
fn percentile_is_fractional_and_bounded() {
    let mut d: Vec<f64> = (1..=10).map(f64::from).collect();
    assert_eq!(stats::percentile(&mut d, 50.0), 5.5);
    let lo = stats::percentile(&mut d, 33.0);
    let mid = stats::percentile(&mut d, 33.3);
    let hi = stats::percentile(&mut d, 34.0);
    assert!(
        lo < mid && mid < hi,
        "fractional percentiles interpolate: {lo} {mid} {hi}"
    );
    assert!(stats::percentile(&mut d, -0.1).is_nan());
    assert!(stats::percentile(&mut d, 100.1).is_nan());
    assert!(stats::percentile(&mut d, f64::NAN).is_nan());
    assert_eq!(stats::percentile(&mut d, 0.0), 1.0);
    assert_eq!(stats::percentile(&mut d, 100.0), 10.0);
}

#[test]
fn mean_variance_are_stable() {
    assert!(
        stats::mean(&[1e308, 1e308]).is_finite(),
        "the sum overflows, the mean does not"
    );
    assert_eq!(stats::mean(&[1e308, 1e308]), 1e308);
    assert!(stats::mean(&[f64::NAN, 1.0]).is_nan());
    assert_eq!(stats::mean(&[f64::INFINITY, 1.0]), f64::INFINITY);
    // A large offset must not destroy the variance (population variance of 4, 7, 13, 16 is 22.5).
    let v = stats::variance(&[1e9 + 4.0, 1e9 + 7.0, 1e9 + 13.0, 1e9 + 16.0]);
    assert!((v - 22.5).abs() < 1e-6, "{v}");
    assert!(stats::correlation(&[1.0, 1.0, 1.0], &[1.0, 2.0, 3.0]).is_nan());
    let mu = stats::mean(&[1e308, -1e308, 5.0]);
    assert!(mu.is_finite());
}

#[test]
fn p_values_keep_their_tails() {
    // A huge t statistic: the two-sided p-value must be tiny but not exactly 0,
    // and must agree with the incomplete beta identity p = I_{df/(df+t^2)}(df/2, 1/2).
    let a: Vec<f64> = (0..30).map(|i| f64::from(i) * 0.01).collect();
    let b: Vec<f64> = (0..30).map(|i| 10.0 + f64::from(i) * 0.01).collect();
    let (t, p) = stats::two_sample_t_test(&a, &b);
    assert!(t < -100.0, "{t}");
    let df = 58.0_f64;
    let want = regularized_beta(df / (df + t * t), df / 2.0, 0.5);
    assert!(p > 0.0 && p < 1e-50, "{p}");
    assert!(
        (p - want).abs() <= 1e-6 * want,
        "p = {p:e}, beta identity {want:e}"
    );
    // Chi-square tail: p = Q(df/2, chi/2)
    let obs = [50.0, 5.0, 5.0, 5.0, 5.0];
    let exp = [10.0; 5];
    let (chi, p) = stats::chi_squared_test(&obs, &exp);
    let want = regularized_upper_gamma(2.0, chi / 2.0);
    assert!(
        p > 0.0 && (p - want).abs() <= 1e-9 * want,
        "p = {p:e}, Q = {want:e}"
    );
    // ANOVA on well separated groups
    let (mut g1, mut g2, mut g3) = (
        vec![1.0, 1.1, 0.9, 1.05],
        vec![50.0, 50.1, 49.9, 50.05],
        vec![100.0, 100.1, 99.9, 100.05],
    );
    let (f, p) = stats::one_way_anova(&mut [&mut g1, &mut g2, &mut g3]);
    assert!(f > 1e6 && p > 0.0 && p < 1e-15, "F = {f}, p = {p:e}");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn stats_never_panic(data in prop::collection::vec(prop_oneof![8 => -1e6..1e6f64, 1 => Just(f64::NAN), 1 => Just(f64::INFINITY), 1 => Just(0.0), 1 => Just(1e308)], 0..12), p in -10.0..110.0f64) {
        let mut a = data.clone();
        let mut b = data.clone();
        let _ = stats::mean(&data);
        let _ = stats::variance(&data);
        let _ = stats::std_dev(&data);
        let _ = stats::median(&mut a);
        let _ = stats::percentile(&mut b, p);
        let _ = stats::iqr(&mut a);
        let _ = stats::min(&mut a);
        let _ = stats::max(&mut a);
        let _ = stats::skewness(&mut a);
        let _ = stats::kurtosis(&mut a);
        let _ = stats::z_scores(&data);
        let _ = stats::mode(&data, 3);
        let _ = stats::geometric_mean(&data);
        let _ = stats::harmonic_mean(&data);
        let _ = stats::shannon_entropy(&data);
        let _ = stats::coefficient_of_variation(&data);
        let _ = stats::standard_error(&data);
        let _ = stats::correlation(&data, &data);
        let _ = stats::two_sample_t_test(&data, &data);
        let _ = stats::welch_t_test(&data, &data);
        let _ = stats::chi_squared_test(&data, &data);
        let pairs: Vec<(f64, f64)> = data.iter().map(|&x| (x, x)).collect();
        let _ = stats::simple_linear_regression(&pairs);
    }
}

// ---------------------------------------------------------------- quadrature

#[test]
fn adaptive_quadrature_is_bounded_on_bad_integrands() {
    // A NaN integrand or a zero tolerance used to recurse 2^100 times.
    let t = std::time::Instant::now();
    assert!(ig::adaptive_quadrature(|_| f64::NAN, (0.0, 1.0), 1e-8).is_nan());
    let v = ig::adaptive_quadrature(|x| x.sin(), (0.0, 1.0), 0.0);
    assert!((v - (1.0 - 1f64.cos())).abs() < 1e-10, "{v}");
    let v = ig::adaptive_quadrature(
        |x| (1000.0 * x).sin().abs() + 1.0 / (x - 0.3).abs().sqrt().max(1e-12),
        (0.0, 1.0),
        1e-9,
    );
    assert!(v.is_finite());
    assert!(ig::adaptive_quadrature(|x| x, (0.0, f64::INFINITY), 1e-8).is_nan());
    assert!(ig::adaptive_quadrature(|x| x, (f64::NAN, 1.0), 1e-8).is_nan());
    assert!(ig::adaptive_quadrature(|x| x, (0.0, 1.0), f64::NAN).is_nan());
    // Endpoint singularities are integrated, not evaluated at the pole.
    let v = ig::adaptive_quadrature(|x| 1.0 / x.sqrt(), (0.0, 1.0), 1e-8);
    assert!((v - 2.0).abs() < 1e-6, "{v}");
    let v = ig::adaptive_quadrature(|x| x.ln(), (0.0, 1.0), 1e-8);
    assert!((v + 1.0).abs() < 1e-6, "{v}");
    assert!(t.elapsed().as_secs_f64() < 10.0, "took {:?}", t.elapsed());
}

#[test]
fn fixed_rules_handle_degenerate_arguments() {
    // A tiny interval is not "zero" because |a - b| < eps.
    let v = ig::gauss_legendre_quadrature(|_| 1e30, (0.0, 1e-20));
    assert!((v - 1e10).abs() < 1.0, "{v}");
    assert!((ig::trapezoidal_rule(|_| 1e30, (0.0, 1e-20), 4) - 1e10).abs() < 1.0);
    assert_eq!(ig::trapezoidal_rule(|x| x, (1.0, 1.0), 10), 0.0);
    assert_eq!(ig::trapezoidal_rule(|x| x, (0.0, 1.0), 0), 0.0);
    assert!(ig::trapezoidal_rule(|x| x, (0.0, f64::NAN), 10).is_nan());
    assert_eq!(ig::simpson_rule(|x| x, (0.0, 1.0), 0).unwrap(), 0.0);
    assert!(ig::simpson_rule(|x| x, (0.0, f64::INFINITY), 10).is_err());
    assert!(
        (ig::simpson_rule(|x| x * x, (0.0, 1.0), 3).unwrap() - 1.0 / 3.0).abs() < 1e-12,
        "odd step count is rounded up"
    );
    // Reversed bounds flip the sign.
    assert!((ig::gauss_kronrod(|x| x, 1.0, 0.0, 1e-12, 100).value + 0.5).abs() < 1e-12);
    // Romberg: huge step counts neither overflow the shifts nor allocate max_steps^2.
    let v = ig::romberg_integration(|x| x.exp(), (0.0, 1.0), 1_000_000);
    assert!((v - (std::f64::consts::E - 1.0)).abs() < 1e-12);
    assert!((ig::romberg_integration(|x| x, (0.0, 1.0), 40) - 0.5).abs() < 1e-12);
    assert!(ig::romberg_integration(|x| x, (f64::NAN, 1.0), 5).is_nan());
}

#[test]
fn infinite_and_singular_ranges() {
    let q = ig::gauss_kronrod_any(
        |x| (-x * x).exp(),
        f64::NEG_INFINITY,
        f64::INFINITY,
        1e-10,
        200,
    );
    assert!((q.value - PI.sqrt()).abs() < 1e-8, "{q:?}");
    let q = ig::gauss_kronrod_any(|x| 1.0 / x.sqrt(), 0.0, 1.0, 1e-10, 200);
    assert!((q.value - 2.0).abs() < 1e-7, "{q:?}");
    let q = ig::gauss_kronrod_any(|x| x.ln(), 0.0, 1.0, 1e-10, 200);
    assert!((q.value + 1.0).abs() < 1e-8, "{q:?}");
    let _ = ig::gauss_kronrod(|x| x, 0.0, 1.0, 1e-10, 0);
    let _ = ig::gauss_kronrod(|x| x, f64::NAN, 1.0, 1e-10, 10);
    let _ = ig::gauss_kronrod(|x| x, 0.0, f64::INFINITY, 1e-10, 10);
    let q = rssn::kernels::quadrature::tanh_sinh(|x| 1.0 / (x * (1.0 - x)).sqrt(), 0.0, 1.0, 1e-10);
    assert!((q.value - PI).abs() < 1e-6, "{q:?}");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn quadrature_never_panics(a in prop_oneof![-5.0..5.0f64, Just(f64::NAN), Just(f64::INFINITY), Just(f64::NEG_INFINITY)], b in prop_oneof![-5.0..5.0f64, Just(f64::NAN), Just(f64::INFINITY), Just(1e308)], n in 0usize..70, tol in prop_oneof![Just(0.0), Just(f64::NAN), 1e-12..1e-2f64]) {
        let f = |x: f64| (3.0 * x).sin() / x;
        let _ = ig::trapezoidal_rule(f, (a, b), n);
        let _ = ig::simpson_rule(f, (a, b), n);
        let _ = ig::romberg_integration(f, (a, b), n);
        let _ = ig::gauss_legendre_quadrature(f, (a, b));
        let _ = ig::adaptive_quadrature(f, (a, b), tol);
        let _ = ig::gauss_kronrod(f, a, b, tol, n);
        let _ = ig::gauss_kronrod_any(f, a, b, tol, n);
    }
}

// ------------------------------------------------------------- interpolation

#[test]
fn interpolation_validates_its_nodes() {
    assert!(
        ip::cubic_spline_interpolation(&[(0.0, 0.0), (0.0, 1.0), (1.0, 2.0)]).is_err(),
        "duplicate x"
    );
    assert!(
        ip::cubic_spline_interpolation(&[(1.0, 0.0), (0.0, 1.0), (2.0, 2.0)]).is_err(),
        "unsorted"
    );
    assert!(
        ip::cubic_spline_interpolation(&[(0.0, 0.0), (f64::NAN, 1.0), (2.0, 2.0)]).is_err(),
        "NaN x"
    );
    assert!(
        ip::cubic_spline_interpolation(&[(0.0, f64::NAN), (1.0, 1.0)]).is_err(),
        "NaN y"
    );
    assert!(ip::cubic_spline_interpolation(&[(0.0, 0.0)]).is_err());
    assert!(ip::cubic_spline_interpolation(&[]).is_err());
    assert!(ip::lagrange_interpolation(&[(0.0, 0.0), (0.0, 1.0)]).is_err());
    assert!(ip::pchip_interpolation(&[(0.0, 0.0), (0.0, 1.0)]).is_err());
    assert!(ip::linear_interpolation(&[(0.0, 0.0), (0.0, 1.0)]).is_err());
    // Close but distinct nodes are not "duplicates" (the old test was |dx| < 1e-9).
    let p = ip::lagrange_interpolation(&[(0.0, 0.0), (1e-10, 1.0)]).unwrap();
    assert!((p.eval(1e-10) - 1.0).abs() < 1e-3);
    // Evaluation at NaN / infinity / far outside is NaN or finite, never a panic.
    let sp = ip::cubic_spline_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]).unwrap();
    let _ = sp(f64::NAN);
    let _ = sp(f64::INFINITY);
    assert!(sp(1e6).is_finite());
    assert!((sp(1.0) - 1.0).abs() < 1e-12);
    assert!((sp(0.5) - 0.6875).abs() < 1e-12);
    // Two points: a straight line.
    let two = ip::cubic_spline_interpolation(&[(0.0, 0.0), (1.0, 1.0)]).unwrap();
    assert!((two(0.25) - 0.25).abs() < 1e-12);
}
