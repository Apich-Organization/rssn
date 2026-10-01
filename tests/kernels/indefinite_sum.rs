//! Indefinite sums and products at non-integer arguments.

use rssn::kernels::indefinite_sum::{
    IndefiniteSumConfig, indefinite_product, indefinite_sum, numeric_taylor_coefficients,
    series_antidiff,
};
use rssn::kernels::special::{digamma_numerical, gamma_numerical};

const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

fn sum(f: impl Fn(f64) -> f64, x: f64) -> f64 {
    indefinite_sum(f, x, &IndefiniteSumConfig::default()).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn sum_of_k_is_x_times_x_minus_one_over_two() {
    for x in [0.5, 2.5, 7.3, -1.4, 12.75] {
        let exact = x * (x - 1.0) / 2.0;
        assert!((sum(|k| k, x) - exact).abs() < 1e-8, "x={x}");
    }
}

#[test]
fn integer_arguments_are_plain_sums() {
    assert_eq!(sum(|k| k * k, 5.0), 30.0);
    assert_eq!(sum(|k| k, 0.0), 0.0);
    // F(x) - F(x-1)... negative direction: F(-2) = -(f(-2) + f(-1)).
    assert_eq!(sum(|k| k, -2.0), 3.0);
}

#[test]
fn polynomials_of_higher_degree_use_the_series_strategy() {
    // sum k^3 = (x(x-1)/2)^2.
    let x = 4.6;
    let exact = (x * (x - 1.0) / 2.0_f64).powi(2);
    assert!((sum(|k| k.powi(3), x) - exact).abs() < 1e-7 * exact);
    // Faulhaber for k^7 at a non-integer point: B_8(x) - B_8(0) over 8.
    let k7 = sum(|k| k.powi(7), 2.5);
    let via_series = series_antidiff(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 0.0, 2.5);
    assert!((k7 - via_series).abs() < 1e-6 * via_series.abs().max(1.0), "{k7} {via_series}");
}

#[test]
fn harmonic_numbers_via_digamma() {
    // sum_{k<x} 1/(k+1) = psi(x+1) + gamma.
    for x in [0.5, 3.25, 10.5] {
        let exact = digamma_numerical(x + 1.0) + EULER_GAMMA;
        assert!((sum(|k| 1.0 / (k + 1.0), x) - exact).abs() < 1e-8, "x={x}");
    }
}

#[test]
fn geometric_sum_is_two_to_the_x_minus_one() {
    for x in [0.5, 2.5, 6.7] {
        let exact = 2.0_f64.powf(x) - 1.0;
        assert!((sum(f64::exp2, x) - exact).abs() < 1e-6 * exact.max(1.0), "x={x}");
    }
}

#[test]
fn step_and_normalisation_point() {
    // F(x+S) - F(x) = f(x) with S = 0.5, F(h) = 0, h = 1: sum of k over 1, 1.5, ...
    let cfg = IndefiniteSumConfig { h: 1.0, step: 0.5, ..IndefiniteSumConfig::default() };
    let value = |x: f64| indefinite_sum(|k| k, x, &cfg).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(value(1.0), 0.0);
    assert!((value(3.0) - (1.0 + 1.5 + 2.0 + 2.5)).abs() < 1e-12);
    assert!((value(3.2) - value(2.7) - 2.7).abs() < 1e-7);
    let bad = IndefiniteSumConfig { step: 0.0, ..IndefiniteSumConfig::default() };
    assert!(indefinite_sum(|k| k, 1.5, &bad).is_err());
}

#[test]
fn recurrence_holds_at_non_integer_points() {
    let f = |k: f64| (k + 1.0).sqrt();
    for x in [0.3, 1.7, 4.4] {
        let d = sum(f, x + 1.0) - sum(f, x);
        assert!((d - f(x)).abs() < 1e-7, "x={x}: {d} vs {}", f(x));
    }
}

#[test]
fn product_of_k_plus_one_is_gamma_of_x_plus_one() {
    let cfg = IndefiniteSumConfig::default();
    for x in [0.5, 3.5, 6.25] {
        let p = indefinite_product(|k| k + 1.0, x, &cfg).unwrap_or_else(|e| panic!("{e}"));
        let exact = gamma_numerical(x + 1.0);
        assert!((p - exact).abs() < 1e-6 * exact, "x={x}: {p} vs {exact}");
    }
    let p = indefinite_product(|k| k + 1.0, 5.0, &cfg).unwrap_or_else(|e| panic!("{e}"));
    assert!((p - 120.0).abs() < 1e-9);
}

#[test]
fn product_rejects_non_positive_factors() {
    let cfg = IndefiniteSumConfig::default();
    assert!(indefinite_product(|k| k - 1.0, 3.0, &cfg).is_err());
}

#[test]
fn taylor_coefficients_of_exponential_and_polynomial() {
    let c = numeric_taylor_coefficients(f64::exp, 0.0, 0.5, 10);
    let mut fact = 1.0;
    for (m, value) in c.iter().enumerate() {
        if m > 0 {
            fact *= m as f64;
        }
        assert!((value - 1.0 / fact).abs() < 1e-6, "m={m}: {value}");
    }
    let c = numeric_taylor_coefficients(|t| 2.0 + 3.0 * (t - 1.0) + (t - 1.0).powi(2), 1.0, 0.5, 5);
    for (value, want) in c.iter().zip([2.0, 3.0, 1.0, 0.0, 0.0]) {
        assert!((value - want).abs() < 1e-10, "{c:?}");
    }
}
