//! Descriptive statistics, hypothesis tests and distributions (ported from
//! `numerical_stats_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::stats::{
    BinomialDist, ExponentialDist, GammaDist, PoissonDist, UniformDist, VarianceType,
    chi_squared_test, coefficient_of_variation, correlation, covariance, geometric_mean,
    harmonic_mean, iqr, kurtosis, max, mean, median, min, mode, one_way_anova, percentile, range,
    shannon_entropy, simple_linear_regression, skewness, standard_error, std_dev,
    two_sample_t_test, variance, variance_with_type, welch_t_test, z_scores,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

const EIGHT: [f64; 8] = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];

#[test]
fn mean_basic_and_empty() {
    assert_approx_eq!(mean(&[1.0, 2.0, 3.0, 4.0, 5.0]), 3.0, 1e-10);
    assert_eq!(mean(&[]), 0.0);
}

#[test]
fn variance_is_population_variance() {
    assert_approx_eq!(variance(&EIGHT), 4.0, 1e-10);
}

#[test]
fn variance_with_type_population_sample_and_edges() {
    assert_approx_eq!(
        variance_with_type(&EIGHT, VarianceType::Population).unwrap_or(f64::NAN),
        4.0,
        1e-10
    );
    assert_approx_eq!(
        variance_with_type(&EIGHT, VarianceType::Sample).unwrap_or(f64::NAN),
        32.0 / 7.0,
        1e-10
    );
    // Welford stability with a large offset.
    let shifted = [1e12 + 1.0, 1e12 + 2.0, 1e12 + 3.0];
    assert_approx_eq!(
        variance_with_type(&shifted, VarianceType::Sample).unwrap_or(f64::NAN),
        1.0,
        1e-10
    );
    assert_eq!(variance_with_type(&[], VarianceType::Sample), None);
    assert_eq!(variance_with_type(&[1.0], VarianceType::Sample), None);
    assert_eq!(
        variance_with_type(&[1.0], VarianceType::Population),
        Some(0.0)
    );
}

#[test]
fn std_dev_is_sample_standard_deviation() {
    assert_approx_eq!(std_dev(&EIGHT), (32.0f64 / 7.0).sqrt(), 1e-10);
}

#[test]
fn geometric_and_harmonic_means() {
    assert_approx_eq!(
        geometric_mean(&[1.0, 2.0, 4.0, 8.0]),
        2.828_427_124_746_190_3,
        1e-10
    );
    assert_approx_eq!(
        harmonic_mean(&[1.0, 2.0, 4.0]),
        1.714_285_714_285_714_2,
        1e-10
    );
    assert!(geometric_mean(&[]).is_nan() && harmonic_mean(&[]).is_nan());
}

#[test]
fn range_min_max_median() {
    assert_approx_eq!(range(&[1.0, 5.0, 3.0, 9.0, 2.0]), 8.0, 1e-10);
    assert!(range(&[]).is_nan());
    assert_eq!(min(&mut [3.0, 1.0, 2.0]), 1.0);
    assert_eq!(max(&mut [3.0, 1.0, 2.0]), 3.0);
    assert_eq!(median(&mut [3.0, 1.0, 2.0]), 2.0);
    assert_eq!(median(&mut [4.0, 1.0, 3.0, 2.0]), 2.5);
}

#[test]
fn percentile_and_iqr() {
    let mut d: Vec<f64> = (1..=10).map(f64::from).collect();
    assert_eq!(percentile(&mut d, 50.0), 5.5);
    assert!(percentile(&mut d, -1.0).is_nan());
    let mut d: Vec<f64> = (1..=9).map(f64::from).collect();
    // Q1 = 3, Q3 = 7 for 1..=9 under the standard-library-free definition used by statrs.
    let spread = iqr(&mut d);
    assert!(spread > 3.0 && spread < 5.0, "iqr = {spread}");
    assert!(iqr(&mut [1.0, 2.0, 3.0]).is_nan());
}

#[test]
fn z_scores_are_centred_and_symmetric() {
    let z = z_scores(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    assert_eq!(z.len(), 5);
    assert!(z[2].abs() < 1e-10);
    assert!((z[0] + z[4]).abs() < 1e-10);
    // With the sample standard deviation sqrt(2.5): z[4] = 2 / sqrt(2.5).
    assert_approx_eq!(z[4], 2.0 / 2.5f64.sqrt(), 1e-10);
    assert!(z_scores(&[]).is_empty());
    assert_eq!(z_scores(&[3.0, 3.0, 3.0]), vec![0.0; 3]);
}

#[test]
fn mode_of_data() {
    assert_eq!(mode(&[1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0], 0), Some(3.0));
    assert_eq!(mode(&[1.0, 2.0, 3.0, 4.0], 0), None);
    assert_eq!(mode(&[], 0), None);
    assert_eq!(mode(&[1.11, 1.12, 5.0], 1), Some(1.1));
}

#[test]
fn covariance_and_correlation() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [2.0, 4.0, 6.0, 8.0, 10.0];
    // Sample covariance: 2 * var_sample(x) = 2 * 2.5.
    assert_approx_eq!(covariance(&x, &y), 5.0, 1e-10);
    assert_approx_eq!(correlation(&x, &y), 1.0, 1e-10);
    assert_approx_eq!(correlation(&x, &[10.0, 8.0, 6.0, 4.0, 2.0]), -1.0, 1e-10);
    assert_approx_eq!(correlation(&x, &x), 1.0, 1e-10);
}

#[test]
fn simple_linear_regression_exact_line() {
    let (slope, intercept) =
        simple_linear_regression(&[(1.0, 3.0), (2.0, 5.0), (3.0, 7.0), (4.0, 9.0)]);
    assert_approx_eq!(slope, 2.0, 1e-10);
    assert_approx_eq!(intercept, 1.0, 1e-10);
    let (s, i) = simple_linear_regression(&[]);
    assert!(s.is_nan() && i.is_nan());
}

#[test]
fn simple_linear_regression_noisy_data_least_squares() {
    // Known least-squares fit: slope 0.6, intercept 2.2 for these five points.
    let (slope, intercept) =
        simple_linear_regression(&[(1.0, 2.0), (2.0, 4.0), (3.0, 5.0), (4.0, 4.0), (5.0, 5.0)]);
    assert_approx_eq!(slope, 0.6, 1e-10);
    assert_approx_eq!(intercept, 2.2, 1e-10);
}

#[test]
fn shannon_entropy_values() {
    assert_approx_eq!(shannon_entropy(&[0.25; 4]), 2.0, 1e-10);
    assert_approx_eq!(shannon_entropy(&[1.0, 0.0, 0.0]), 0.0, 1e-10);
    assert_approx_eq!(shannon_entropy(&[0.5, 0.5]), 1.0, 1e-10);
}

#[test]
fn welch_t_test_identical_samples() {
    let s = [1.0, 2.0, 3.0, 4.0, 5.0];
    let (t, p) = welch_t_test(&s, &s);
    assert!(t.abs() < 1e-10);
    assert_approx_eq!(p, 1.0, 1e-6);
}

#[test]
fn welch_t_test_shifted_samples() {
    // Equal sample variances (2.5) and n = 5 each: t = -1, df = 8, two-sided p = 0.346593507...
    let (t, p) = welch_t_test(&[1.0, 2.0, 3.0, 4.0, 5.0], &[2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_approx_eq!(t, -1.0, 1e-10);
    assert_approx_eq!(p, 0.346_593_507_087_334_16, 1e-6);
    let (t, p) = welch_t_test(&[1.0], &[1.0, 2.0]);
    assert!(t.is_nan() && p.is_nan());
}

#[test]
fn two_sample_t_test_shifted_samples() {
    let (t, p) = two_sample_t_test(&[1.0, 2.0, 3.0, 4.0, 5.0], &[2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_approx_eq!(t, -1.0, 1e-10);
    assert_approx_eq!(p, 0.346_593_507_087_334_16, 1e-6);
}

#[test]
fn two_sample_t_test_identical_samples() {
    let s = [1.0, 2.0, 3.0, 4.0, 5.0];
    let (t, p) = two_sample_t_test(&s, &s);
    assert!(t.abs() < 1e-10);
    assert_approx_eq!(p, 1.0, 1e-6);
}

#[test]
fn chi_squared_perfect_match_and_known_value() {
    let (chi, p) = chi_squared_test(&[10.0, 20.0, 30.0], &[10.0, 20.0, 30.0]);
    assert!(chi.abs() < 1e-10);
    assert_approx_eq!(p, 1.0, 1e-9);
    // chi^2 = 4/12 + 4/18 = 5/9, df = 2 => p = exp(-chi^2 / 2)
    let (chi, p) = chi_squared_test(&[10.0, 20.0, 30.0], &[12.0, 18.0, 30.0]);
    assert_approx_eq!(chi, 5.0 / 9.0, 1e-10);
    assert_approx_eq!(p, (-5.0f64 / 18.0).exp(), 1e-9);
    assert!(chi_squared_test(&[1.0], &[1.0, 2.0]).0.is_nan());
}

#[test]
fn anova_known_f_statistic() {
    let (mut a, mut b, mut c) = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]);
    let (f, p) = one_way_anova(&mut [&mut a, &mut b, &mut c]);
    // SSB = 54 (df 2), SSW = 6 (df 6) => F = 27; p = (6 / (6 + 2F))^3 = 1e-3.
    assert_approx_eq!(f, 27.0, 1e-10);
    assert_approx_eq!(p, 1e-3, 1e-9);
}

#[test]
fn coefficient_of_variation_and_standard_error() {
    let data = [10.0, 20.0, 30.0, 40.0, 50.0];
    // sample sd = sqrt(250), mean = 30
    assert_approx_eq!(coefficient_of_variation(&data), 250f64.sqrt() / 30.0, 1e-10);
    assert!(coefficient_of_variation(&[-1.0, 1.0]).is_nan());
    let d = [1.0, 2.0, 3.0, 4.0, 5.0];
    assert_approx_eq!(standard_error(&d), std_dev(&d) / 5f64.sqrt(), 1e-10);
    assert_approx_eq!(standard_error(&d), (2.5f64 / 5.0).sqrt(), 1e-10);
}

#[test]
fn skewness_sign() {
    assert!(skewness(&mut [1.0, 1.0, 1.0, 2.0, 10.0]) > 0.0);
    assert!(skewness(&mut [-10.0, -2.0, -1.0, -1.0, -1.0]) < 0.0);
    assert!(skewness(&mut [1.0, 2.0, 3.0, 4.0, 5.0]).abs() < 1e-12);
}

#[test]
fn kurtosis_matches_unbiased_excess_kurtosis() {
    let k = kurtosis(&mut EIGHT.clone());
    assert_approx_eq!(k, 0.940_625, 1e-9);
}

#[test]
fn kurtosis_edge_cases() {
    assert!(kurtosis(&mut [1.0, 2.0, 3.0]).is_nan());
    assert_eq!(kurtosis(&mut [2.0, 2.0, 2.0, 2.0]), 0.0);
}

#[test]
fn uniform_distribution() {
    let u = UniformDist::new(2.0, 6.0).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(u.pdf(3.0), 0.25, 1e-12);
    assert_approx_eq!(u.cdf(3.0), 0.25, 1e-12);
    assert_eq!(u.pdf(7.0), 0.0);
    assert_eq!(u.cdf(7.0), 1.0);
    assert!(UniformDist::new(1.0, 0.0).is_err());
}

#[test]
fn binomial_distribution() {
    let b = BinomialDist::new(10, 0.5).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(b.pmf(5), 252.0 / 1024.0, 1e-12);
    assert_approx_eq!(b.cdf(10), 1.0, 1e-12);
    assert_approx_eq!(b.cdf(4) + b.pmf(5) + (1.0 - b.cdf(5)), 1.0, 1e-12);
    assert!(BinomialDist::new(10, 1.5).is_err());
}

#[test]
fn poisson_exponential_gamma_distributions() {
    let p = PoissonDist::new(2.0).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(p.pmf(0), (-2.0f64).exp(), 1e-12);
    assert_approx_eq!(p.pmf(2), 2.0 * (-2.0f64).exp(), 1e-12);
    assert_approx_eq!(p.cdf(1), 3.0 * (-2.0f64).exp(), 1e-12);

    let e = ExponentialDist::new(1.0).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(e.pdf(0.0), 1.0, 1e-12);
    assert_approx_eq!(e.cdf(1.0), 1.0 - (-1.0f64).exp(), 1e-12);

    // Gamma(shape 2, rate 1): cdf(x) = 1 - (1 + x) e^-x
    let g = GammaDist::new(2.0, 1.0).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(g.cdf(3.0), 1.0 - 4.0 * (-3.0f64).exp(), 1e-10);
    assert_approx_eq!(g.pdf(1.0), (-1.0f64).exp(), 1e-12);
    assert!(PoissonDist::new(-1.0).is_err());
    assert!(ExponentialDist::new(0.0).is_err());
    assert!(GammaDist::new(-1.0, 1.0).is_err());
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_mean_and_variance_of_constant(val in -1000.0..1000.0f64, n in 2usize..100) {
        let data = vec![val; n];
        prop_assert!((mean(&data) - val).abs() < 1e-10);
        prop_assert!(variance(&data).abs() < 1e-10);
        prop_assert!(range(&data).abs() < 1e-10);
    }

    #[test]
    fn prop_variance_is_shift_invariant(data in prop::collection::vec(-100.0..100.0f64, 2..30), shift in -1000.0..1000.0f64) {
        let shifted: Vec<f64> = data.iter().map(|x| x + shift).collect();
        prop_assert!((variance(&data) - variance(&shifted)).abs() < 1e-6);
    }

    #[test]
    fn prop_variance_scales_quadratically(data in prop::collection::vec(-100.0..100.0f64, 2..30), s in -5.0..5.0f64) {
        let scaled: Vec<f64> = data.iter().map(|x| x * s).collect();
        prop_assert!((variance(&scaled) - s * s * variance(&data)).abs() < 1e-6 * (1.0 + variance(&data)));
    }

    #[test]
    fn prop_z_scores_have_zero_mean(n in 3usize..50) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let z = z_scores(&data);
        prop_assert!((z.iter().sum::<f64>() / z.len() as f64).abs() < 1e-10);
    }

    #[test]
    fn prop_am_gm_hm_inequality(data in prop::collection::vec(0.1..100.0f64, 2..20)) {
        let (am, gm, hm) = (mean(&data), geometric_mean(&data), harmonic_mean(&data));
        prop_assert!(am >= gm - 1e-9);
        prop_assert!(gm >= hm - 1e-9);
    }

    #[test]
    fn prop_correlation_is_bounded_and_symmetric(
        x in prop::collection::vec(-100.0..100.0f64, 5..20), noise in prop::collection::vec(-10.0..10.0f64, 5..20),
    ) {
        let n = x.len().min(noise.len());
        let (x, y): (Vec<f64>, Vec<f64>) = (x[..n].to_vec(), x[..n].iter().zip(&noise).map(|(a, b)| 0.5 * a + b).collect());
        prop_assume!(std_dev(&x) > 1e-6 && std_dev(&y) > 1e-6);
        let c = correlation(&x, &y);
        prop_assert!(c.abs() <= 1.0 + 1e-9);
        prop_assert!((c - correlation(&y, &x)).abs() < 1e-12);
    }

    #[test]
    fn prop_self_correlation_is_one(n in 3usize..20) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        prop_assert!((correlation(&data, &data) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn prop_regression_recovers_planted_line(a in -10.0..10.0f64, b in -10.0..10.0f64) {
        let pts: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, a * i as f64 + b)).collect();
        let (slope, icpt) = simple_linear_regression(&pts);
        prop_assert!((slope - a).abs() < 1e-8 && (icpt - b).abs() < 1e-8);
    }

    #[test]
    fn prop_entropy_bounds(n in 2usize..10) {
        let h = shannon_entropy(&vec![1.0 / n as f64; n]);
        prop_assert!(h >= 0.0);
        prop_assert!((h - (n as f64).log2()).abs() < 1e-9);
    }

    #[test]
    fn prop_uniform_cdf_is_monotone_and_linear(lo in -10.0..0.0f64, width in 0.5..10.0f64, f in 0.0..1.0f64) {
        let u = UniformDist::new(lo, lo + width).map_err(TestCaseError::fail)?;
        let x = lo + f * width;
        prop_assert!((u.cdf(x) - f).abs() < 1e-9);
        prop_assert!(u.cdf(x + 0.1) >= u.cdf(x));
    }
}

#[test]
fn skewness_is_adjusted_fisher_pearson() {
    // Hand computation: mean 3, m2 = 12.4, m3 = 63.6, g1 = 1.45657, G1 = sqrt(20)/3 g1 = 2.17130.
    let g = skewness(&mut [1.0, 1.0, 1.0, 2.0, 10.0]);
    assert_approx_eq!(g, 20f64.sqrt() / 3.0 * 63.6 / 12.4f64.powf(1.5), 1e-12);
    assert_approx_eq!(g, 2.1713, 1e-4);
    // Sign flips under reflection, invariant under shift and scale.
    let a = skewness(&mut [-10.0, -2.0, -1.0, -1.0, -1.0]);
    assert_approx_eq!(a, -g, 1e-12);
    let b = skewness(&mut [12.0, 12.0, 12.0, 14.0, 30.0]);
    assert_approx_eq!(b, g, 1e-12);
    assert!(skewness(&mut [1.0, 2.0]).is_nan());
    assert_eq!(skewness(&mut [3.0, 3.0, 3.0]), 0.0);
}

#[test]
fn kurtosis_is_shift_and_scale_invariant_and_flat_data_is_negative() {
    let base = kurtosis(&mut EIGHT.clone());
    let mut moved: Vec<f64> = EIGHT.iter().map(|v| 3.0 * v - 7.0).collect();
    assert_approx_eq!(kurtosis(&mut moved), base, 1e-10);
    // Two-point-mass data is the flattest possible: strongly negative excess kurtosis.
    assert!(kurtosis(&mut [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0]) < -1.0);
}

#[test]
fn two_sample_t_test_agrees_with_welch_for_equal_sizes_and_unequal_sizes_hand_value() {
    // For equal n the pooled and Welch t statistics coincide.
    let a = [1.0, 2.0, 4.0, 8.0, 3.0];
    let b = [2.0, 3.0, 9.0, 7.0, 6.0];
    let (t_pool, _) = two_sample_t_test(&a, &b);
    let (t_welch, _) = welch_t_test(&a, &b);
    assert_approx_eq!(t_pool, t_welch, 1e-12);
    // Unequal sizes: x = [1,2,3] (mean 2, s^2 1), y = [4,6,8,10] (mean 7, s^2 20/3);
    // sp^2 = (2*1 + 3*20/3)/5 = 4.4, t = -5 / sqrt(4.4 (1/3 + 1/4)) = -3.1218...
    let (t, p) = two_sample_t_test(&[1.0, 2.0, 3.0], &[4.0, 6.0, 8.0, 10.0]);
    assert_approx_eq!(t, -5.0 / (4.4f64 * (1.0 / 3.0 + 0.25)).sqrt(), 1e-12);
    assert!(p > 0.0 && p < 0.05);
    assert!(two_sample_t_test(&[1.0], &[1.0, 2.0]).0.is_nan());
}
