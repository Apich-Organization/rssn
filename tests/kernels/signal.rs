//! Signal processing (ported from `numerical_signal_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::signal::{convolve, cross_correlation, fft, hamming_window, hann_window};
use rustfft::num_complex::Complex;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn fft_of_constant_signal_is_a_dc_spike() {
    let mut input = vec![Complex::new(1.0, 0.0); 4];
    let out = fft(&mut input);
    assert!((out[0].re - 4.0).abs() < 1e-9);
    for k in 1..4 {
        assert!(out[k].norm() < 1e-9, "bin {k} = {}", out[k]);
    }
}

#[test]
fn fft_of_cosine_has_two_symmetric_bins() {
    let n = 8;
    let mut x: Vec<Complex<f64>> = (0..n)
        .map(|i| Complex::new((2.0 * std::f64::consts::PI * i as f64 / n as f64).cos(), 0.0))
        .collect();
    let out = fft(&mut x);
    for (k, v) in out.iter().enumerate() {
        let want = if k == 1 || k == n - 1 { n as f64 / 2.0 } else { 0.0 };
        assert!((v.re - want).abs() < 1e-9 && v.im.abs() < 1e-9, "bin {k}: {v}");
    }
}

#[test]
fn fft_does_not_modify_its_input() {
    let mut x = vec![Complex::new(1.0, 2.0), Complex::new(3.0, -1.0)];
    let before = x.clone();
    let out = fft(&mut x);
    assert_eq!(x, before);
    assert!((out[0] - Complex::new(4.0, 1.0)).norm() < 1e-12);
    assert!((out[1] - Complex::new(-2.0, 3.0)).norm() < 1e-12);
}

#[test]
fn convolve_known_result() {
    let result = convolve(&[1.0, 2.0, 3.0], &[0.0, 1.0, 0.5]);
    assert_eq!(result, vec![0.0, 1.0, 2.5, 4.0, 1.5]);
}

#[test]
fn convolve_is_polynomial_multiplication() {
    // (1 + x)(1 - x) = 1 - x^2
    assert_eq!(convolve(&[1.0, 1.0], &[1.0, -1.0]), vec![1.0, 0.0, -1.0]);
}

#[test]
fn cross_correlation_known_result() {
    // v reversed is [0.5, 1, 0]; convolving with a = [1, 2, 3].
    let result = cross_correlation(&[1.0, 2.0, 3.0], &[0.0, 1.0, 0.5]);
    assert_eq!(result, vec![0.5, 2.0, 3.5, 3.0, 0.0]);
}

#[test]
fn hann_and_hamming_windows() {
    let hann = hann_window(5);
    assert_eq!(hann.len(), 5);
    assert!(hann[0].abs() < 1e-9);
    assert!((hann[1] - 0.5).abs() < 1e-9);
    assert!((hann[2] - 1.0).abs() < 1e-9);
    assert!((hann[3] - 0.5).abs() < 1e-9);
    assert!(hann[4].abs() < 1e-9);

    let hamming = hamming_window(5);
    assert!((hamming[0] - 0.08).abs() < 1e-9);
    assert!((hamming[1] - 0.54).abs() < 1e-9);
    assert!((hamming[2] - 1.0).abs() < 1e-9);
    assert!((hamming[4] - 0.08).abs() < 1e-9);
}

#[test]
fn window_edge_lengths() {
    assert!(hann_window(0).is_empty());
    assert_eq!(hann_window(1), vec![1.0]);
    assert!(hamming_window(0).is_empty());
    assert_eq!(hamming_window(1), vec![1.0]);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_convolve_with_unit_impulse_is_identity(a in prop::collection::vec(-10.0..10.0f64, 1..10)) {
        prop_assert_eq!(convolve(&a, &[1.0]), a);
    }

    #[test]
    fn prop_convolve_commutes(a in prop::collection::vec(-10.0..10.0f64, 1..8),
                              b in prop::collection::vec(-10.0..10.0f64, 1..8)) {
        let (ab, ba) = (convolve(&a, &b), convolve(&b, &a));
        prop_assert_eq!(ab.len(), a.len() + b.len() - 1);
        for (x, y) in ab.iter().zip(&ba) {
            prop_assert!((x - y).abs() < 1e-9);
        }
    }

    #[test]
    fn prop_fft_parseval(re in prop::collection::vec(-10.0..10.0f64, 16)) {
        let mut x: Vec<Complex<f64>> = re.iter().map(|&r| Complex::new(r, 0.0)).collect();
        let energy_time: f64 = x.iter().map(|c| c.norm_sqr()).sum();
        let out = fft(&mut x);
        let energy_freq: f64 = out.iter().map(|c| c.norm_sqr()).sum::<f64>() / 16.0;
        prop_assert!((energy_time - energy_freq).abs() < 1e-8 * energy_time.max(1.0));
    }
}
