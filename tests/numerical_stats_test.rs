use rssn::numerical::stats::*;

#[test]
fn test_mean() {
    let data = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    assert!((mean(&data) - 3.0).abs() < 1e-10);

    let empty: Vec<f64> = vec![];
    assert_eq!(mean(&empty), 0.0);
}

#[test]
fn test_variance() {
    let data = vec![
        2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0,
    ];
    let var = variance(&data);
    assert!((var - 4.0).abs() < 1e-10);
}

#[test]
fn test_geometric_mean() {
    let data = vec![1.0, 2.0, 4.0, 8.0];
    let gm = geometric_mean(&data);
    // (1*2*4*8)^(1/4) = 64^(1/4) = 2.828...
    assert!((gm - 2.8284271247461903).abs() < 1e-10);
}

#[test]
fn test_harmonic_mean() {
    let data = vec![1.0, 2.0, 4.0];
    let hm = harmonic_mean(&data);
    // 3 / (1/1 + 1/2 + 1/4) = 3 / 1.75 = 1.714...
    assert!((hm - 1.7142857142857142).abs() < 1e-10);
}

#[test]
fn test_range() {
    let data = vec![
        1.0, 5.0, 3.0, 9.0, 2.0,
    ];
    assert!((range(&data) - 8.0).abs() < 1e-10);
}

#[test]
fn test_z_scores() {
    let data = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let z = z_scores(&data);
    assert_eq!(z.len(), 5);
    // Mean is 3, so z[2] (the value 3) should be 0
    assert!(z[2].abs() < 1e-10);
    // z[0] and z[4] should be symmetric
    assert!((z[0] + z[4]).abs() < 1e-10);
}

#[test]
fn test_mode() {
    let data = vec![
        1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0,
    ];
    let m = mode(&data, 0);
    assert_eq!(m, Some(3.0));

    // No mode when all unique
    let unique = vec![1.0, 2.0, 3.0, 4.0];
    assert_eq!(mode(&unique, 0), None);
}

#[test]
fn test_covariance() {
    let x = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let y = vec![
        2.0, 4.0, 6.0, 8.0, 10.0,
    ];
    let cov = covariance(&x, &y);
    // Perfect positive correlation, cov = 2.5 * 2 = 5 (sample covariance)
    assert!(cov > 0.0);
}

#[test]
fn test_correlation() {
    let x = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let y = vec![
        2.0, 4.0, 6.0, 8.0, 10.0,
    ];
    let corr = correlation(&x, &y);
    assert!((corr - 1.0).abs() < 1e-10); // Perfect positive correlation

    let y_neg = vec![
        10.0, 8.0, 6.0, 4.0, 2.0,
    ];
    let corr_neg = correlation(&x, &y_neg);
    assert!((corr_neg + 1.0).abs() < 1e-10); // Perfect negative correlation
}

#[test]
fn test_simple_linear_regression() {
    // y = 2x + 1
    let data = vec![
        (1.0, 3.0),
        (2.0, 5.0),
        (3.0, 7.0),
        (4.0, 9.0),
    ];
    let (slope, intercept) = simple_linear_regression(&data);
    assert!((slope - 2.0).abs() < 1e-10, "slope was {}", slope);
    assert!(
        (intercept - 1.0).abs() < 1e-10,
        "intercept was {}",
        intercept
    );
}

#[test]
fn test_shannon_entropy() {
    // Uniform distribution: H = log2(n)
    let uniform = vec![
        0.25, 0.25, 0.25, 0.25,
    ];
    let h = shannon_entropy(&uniform);
    assert!((h - 2.0).abs() < 1e-10); // log2(4) = 2

    // Certain event: H = 0
    let certain = vec![1.0, 0.0, 0.0];
    assert!((shannon_entropy(&certain) - 0.0).abs() < 1e-10);
}

#[test]
fn test_welch_t_test() {
    let s1 = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let s2 = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let (t, p) = welch_t_test(&s1, &s2);
    // Same samples should have t ≈ 0 and p ≈ 1
    assert!(t.abs() < 1e-10);
    assert!((p - 1.0).abs() < 0.1);
}

#[test]
fn test_chi_squared() {
    let observed = vec![10.0, 20.0, 30.0];
    let expected = vec![10.0, 20.0, 30.0];
    let (chi, p) = chi_squared_test(&observed, &expected);
    assert!(chi.abs() < 1e-10); // Perfect match
    assert!((p - 1.0).abs() < 0.1);
}

#[test]
fn test_coefficient_of_variation() {
    let data = vec![
        10.0, 20.0, 30.0, 40.0, 50.0,
    ];
    let cv = coefficient_of_variation(&data);
    assert!(cv > 0.0);
}

#[test]
fn test_standard_error() {
    let data = vec![
        1.0, 2.0, 3.0, 4.0, 5.0,
    ];
    let se = standard_error(&data);
    // SE = std_dev / sqrt(n)
    let expected = std_dev(&data) / (5.0_f64).sqrt();
    assert!((se - expected).abs() < 1e-10);
}

// Property-based tests
proptest::proptest! {
    /// Mean of identical values should equal that value
    #[test]
    fn prop_mean_of_constant(val in -1000.0..1000.0, n in 2usize..100) {
        let data: Vec<f64> = vec![val; n];
        proptest::prop_assert!((mean(&data) - val).abs() < 1e-10);
    }

    /// Variance of identical values should be 0
    #[test]
    fn prop_variance_of_constant(val in -1000.0..1000.0, n in 2usize..100) {
        let data: Vec<f64> = vec![val; n];
        proptest::prop_assert!(variance(&data).abs() < 1e-10);
    }

    /// Range of identical values should be 0
    #[test]
    fn prop_range_of_constant(val in -1000.0..1000.0, n in 2usize..100) {
        let data: Vec<f64> = vec![val; n];
        proptest::prop_assert!(range(&data).abs() < 1e-10);
    }

    /// Z-scores should have mean 0
    #[test]
    fn prop_z_scores_mean_zero(n in 3usize..50) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let z = z_scores(&data);
        let z_mean: f64 = z.iter().sum::<f64>() / z.len() as f64;
        proptest::prop_assert!(z_mean.abs() < 1e-10);
    }

    /// Geometric mean of positive values should be positive
    #[test]
    fn prop_geometric_mean_positive(n in 2usize..20) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        proptest::prop_assert!(geometric_mean(&data) > 0.0);
    }

    /// AM >= GM >= HM for positive values
    #[test]
    fn prop_mean_inequality(n in 2usize..20) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let am = mean(&data);
        let gm = geometric_mean(&data);
        let hm = harmonic_mean(&data);
        proptest::prop_assert!(am >= gm - 1e-10);
        proptest::prop_assert!(gm >= hm - 1e-10);
    }

    /// Correlation of a dataset with itself should be 1
    #[test]
    fn prop_self_correlation(n in 3usize..20) {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let corr = correlation(&data, &data);
        proptest::prop_assert!((corr - 1.0).abs() < 1e-10);
    }

    /// Shannon entropy should be non-negative
    #[test]
    fn prop_entropy_nonnegative(n in 2usize..10) {
        let prob = 1.0 / n as f64;
        let data: Vec<f64> = vec![prob; n];
        proptest::prop_assert!(shannon_entropy(&data) >= 0.0);
    }
}
