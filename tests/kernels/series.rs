//! Power series, partial sums and accelerated infinite sums.

use rssn::kernels::series::{evaluate_power_series, sum_range, sum_to_infinity};

#[test]
fn power_series_matches_exponential_taylor_polynomial() {
    let coefficients: Vec<f64> = (0..15)
        .scan(1.0, |fact, n| {
            if n > 0 {
                *fact *= f64::from(n);
            }
            Some(1.0 / *fact)
        })
        .collect();
    assert!((evaluate_power_series(&coefficients, 0.0, 1.0) - std::f64::consts::E).abs() < 1e-10);
    // Centered at 1: e^x = e * sum (x-1)^n / n!.
    let at_one = evaluate_power_series(&coefficients, 1.0, 2.0) * std::f64::consts::E;
    assert!((at_one - std::f64::consts::E.powi(2)).abs() < 1e-9);
}

#[test]
fn power_series_edge_cases() {
    assert_eq!(evaluate_power_series(&[], 0.0, 5.0), 0.0);
    assert_eq!(evaluate_power_series(&[7.0], 3.0, 100.0), 7.0);
    assert_eq!(evaluate_power_series(&[1.0, 2.0, 3.0], 1.0, 1.0), 1.0);
}

#[test]
fn infinite_sums_still_accelerate() {
    let leibniz = sum_to_infinity(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0, 1e-12, 10_000);
    assert!((leibniz.value - std::f64::consts::FRAC_PI_4).abs() < 1e-10, "{leibniz:?}");
    let basel = sum_to_infinity(|k| 1.0 / (k * k), 1, 1e-10, 200_000);
    assert!((basel.value - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-9);
    assert_eq!(sum_range(|k| k, 1, 10), 55.0);
}
