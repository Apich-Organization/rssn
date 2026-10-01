//! Sequence acceleration on classical slowly converging series.

use std::f64::consts::{LN_2, PI};

use rssn::kernels::convergence::{
    aitken_acceleration, find_sequence_limit, richardson_extrapolation,
    richardson_extrapolation_with, sum_series_numerical, wynn_epsilon,
};

fn partial_sums(
    term: impl Fn(f64) -> f64,
    first: f64,
    count: usize,
) -> Vec<f64> {
    let mut sum = 0.0;
    (0..count)
        .map(|i| {
            sum += term(first + i as f64);
            sum
        })
        .collect()
}

fn last(v: &[f64]) -> f64 {
    v.last().copied().unwrap_or(f64::NAN)
}

#[test]
fn aitken_on_leibniz_and_alternating_harmonic() {
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 20);
    let acc = aitken_acceleration(&leibniz);
    assert_eq!(acc.len(), 18);
    assert!((last(&acc) - PI / 4.0).abs() < (leibniz[19] - PI / 4.0).abs() / 100.0);
    let harmonic = partial_sums(|k| (-1.0_f64).powf(k + 1.0) / k, 1.0, 20);
    assert!((last(&aitken_acceleration(&harmonic)) - LN_2).abs() < 1e-3);
    assert!(aitken_acceleration(&[1.0, 2.0]).is_empty());
}

#[test]
fn aitken_is_exact_for_geometric_error() {
    let seq: Vec<f64> = (0..6).map(|n| 3.0 + 0.5_f64.powi(n)).collect();
    for v in aitken_acceleration(&seq) {
        assert!((v - 3.0).abs() < 1e-12, "{v}");
    }
}

#[test]
fn wynn_epsilon_on_alternating_series() {
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 13);
    let est = wynn_epsilon(&leibniz);
    assert!((last(&est) - PI / 4.0).abs() < 1e-8, "{est:?}");
    assert_eq!(est[0], leibniz[12]);
    let harmonic = partial_sums(|k| (-1.0_f64).powf(k + 1.0) / k, 1.0, 13);
    assert!((last(&wynn_epsilon(&harmonic)) - LN_2).abs() < 1e-8);
    assert!(wynn_epsilon(&[]).is_empty());
    // Exactly converged sequences are returned as they are.
    assert_eq!(last(&wynn_epsilon(&[1.0, 1.0, 1.0])), 1.0);
}

#[test]
fn richardson_on_basel_partial_sums_at_doubling_counts() {
    let counts = [8_u32, 16, 32, 64, 128, 256, 512];
    let sums: Vec<f64> = counts
        .iter()
        .map(|&n| (1..=n).map(|k| 1.0 / f64::from(k * k)).sum())
        .collect();
    let rich = richardson_extrapolation_with(&sums, 2.0);
    assert!((last(&rich) - PI * PI / 6.0).abs() < 1e-8, "{}", last(&rich));
    assert!((last(&rich) - PI * PI / 6.0).abs() < (sums[6] - PI * PI / 6.0).abs() * 1e-5);
}

#[test]
fn richardson_removes_h_squared_errors() {
    // Trapezoid-like approximations of 1 with error c2 h^2 + c4 h^4.
    let seq: Vec<f64> = (0..4)
        .map(|i| {
            let h = 0.5_f64.powi(i);
            1.0 + 3.0 * h * h - 2.0 * h.powi(4)
        })
        .collect();
    let rich = richardson_extrapolation(&seq);
    assert_eq!(rich[0], seq[0]);
    assert!((rich[2] - 1.0).abs() < 1e-12, "{rich:?}");
    assert!(richardson_extrapolation(&[]).is_empty());
}

#[test]
fn find_sequence_limit_finds_known_limits() {
    let l = find_sequence_limit(|n| 2.0 + 0.5_f64.powf(n), 20, 1e-8)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((l - 2.0).abs() < 1e-8);
    // Partial sums of the Leibniz series.
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 40);
    let l = find_sequence_limit(|n| leibniz[n as usize], 40, 1e-6)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((l - PI / 4.0).abs() < 1e-5, "{l}");
    assert!(find_sequence_limit(|_| f64::NAN, 10, 1e-8).is_err());
    assert!(find_sequence_limit(|n| n, 10, 1e-8).is_err());
}

#[test]
fn sum_series_numerical_stops_at_small_terms() {
    let s = sum_series_numerical(|n| 0.5_f64.powf(n), 0, 200, 1e-14);
    assert!((s - 2.0).abs() < 1e-12);
    let s = sum_series_numerical(|n| n, 1, 10, 0.0);
    assert_eq!(s, 55.0);
}
