//! Signal processing (ported from `numerical_signal_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::signal::blackman_window;
use rssn::kernels::signal::convolve;
use rssn::kernels::signal::cross_correlation;
use rssn::kernels::signal::fft;
use rssn::kernels::signal::hamming_window;
use rssn::kernels::signal::hann_window;
use rssn::kernels::signal::ifft;
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
        .map(|i| {
            Complex::new(
                (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos(),
                0.0,
            )
        })
        .collect();
    let out = fft(&mut x);
    for (k, v) in out.iter().enumerate() {
        let want = if k == 1 || k == n - 1 {
            n as f64 / 2.0
        } else {
            0.0
        };
        assert!(
            (v.re - want).abs() < 1e-9 && v.im.abs() < 1e-9,
            "bin {k}: {v}"
        );
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
fn short_convolution_matches_exact_dyadic_coefficients() {
    // Every coefficient fits in i64 and divides by a power of two, so the
    // integer oracle tests exact values rather than a tolerance.
    for (length, width) in [
        (64, 1),
        (259, 8),
        (2_048, 16),
        (2_049, 2),
        (32_768, 8),
        (33_011, 16),
    ] {
        let left: Vec<i64> = (0..length).map(|i| (i * 47) % 67 - 33).collect();
        let taps: Vec<i64> = (0..width).map(|i| (i * 29) % 37 - 18).collect();
        let samples: Vec<_> = left.iter().map(|&x| x as f64 / 16.0).collect();
        let filter: Vec<_> = taps.iter().map(|&x| x as f64 / 32.0).collect();
        for reverse in [false, true] {
            let values = if reverse {
                cross_correlation(&samples, &filter)
            } else {
                convolve(&samples, &filter)
            };
            assert_eq!(values.len(), left.len() + taps.len() - 1);
            for (k, value) in values.into_iter().enumerate() {
                let exact: i64 = taps
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &tap)| {
                        let offset = if reverse {
                            taps.len() - 1 - j
                        } else {
                            j
                        };
                        left.get(k.checked_sub(offset)?).map(|&x| x * tap)
                    })
                    .sum();
                assert_eq!(value, exact as f64 / 512.0, "coefficient {k}");
            }
        }
    }
}

#[test]
fn short_convolution_keeps_legacy_float_order_and_non_finite_locality() {
    let mut samples = vec![-0.0; 32_771];
    let half_sqrt_epsilon = 2.0_f64.powi(-27);
    samples[2_047] = 1.0;
    samples[2_048] = 1.0 + half_sqrt_epsilon;
    samples[4_095] = f64::INFINITY;
    samples[8_191] = f64::NEG_INFINITY;
    samples[16_383] = f64::NAN;
    samples[24_575] = f64::MAX;
    samples[24_576] = f64::MAX;
    samples[32_768] = f64::from_bits(1);
    let taps = [1.0 - half_sqrt_epsilon, -1.0, -0.0, 2.0, -2.0, 0.5];
    let mut expected = vec![0.0; samples.len() + taps.len() - 1];
    for (i, &x) in samples.iter().enumerate() {
        for (j, &y) in taps.iter().enumerate() {
            expected[i + j] += x * y;
        }
    }

    // This coefficient would change if a fused multiply-add were introduced.
    assert_eq!(expected[2_048].to_bits(), 0.0f64.to_bits());
    assert!((1.0 + half_sqrt_epsilon).mul_add(taps[0], -1.0) < 0.0);
    for threads in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let values = pool.install(|| convolve(&samples, &taps));
        assert_eq!(values.len(), expected.len());
        for (k, (&value, &want)) in values.iter().zip(&expected).enumerate() {
            if want.is_nan() {
                assert!(value.is_nan(), "coefficient {k} with {threads} workers");
            } else {
                assert_eq!(
                    value.to_bits(),
                    want.to_bits(),
                    "coefficient {k} with {threads} workers"
                );
            }
        }
        assert_eq!(values[4_094].to_bits(), 0.0f64.to_bits());
        assert_eq!(values[4_095], f64::INFINITY);
        assert_eq!(values[4_096], f64::NEG_INFINITY);
        assert!(values[4_097].is_nan());
        assert_eq!(values[4_101].to_bits(), 0.0f64.to_bits());
    }
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

#[test]
fn ifft_inverts_fft_and_normalises_by_n() {
    let n = 32;
    let mut input: Vec<Complex<f64>> = (0..n)
        .map(|i| {
            let t = i as f64;
            Complex::new((0.3 * t).sin() + 0.5, (0.7 * t).cos() - 0.2)
        })
        .collect();
    let original = input.clone();
    let spectrum = fft(&mut input);
    let mut spectrum_mut = spectrum.clone();
    let back = ifft(&mut spectrum_mut);
    for (got, want) in back.iter().zip(&original) {
        assert!((got - want).norm() < 1e-12, "{got} vs {want}");
    }
    // The transform does not modify its input, like `fft`.
    assert_eq!(spectrum_mut, spectrum);
    // A DC spectrum inverses to a constant of value sum/n.
    let mut dc = vec![Complex::new(0.0, 0.0); 4];
    dc[0] = Complex::new(8.0, 0.0);
    for value in ifft(&mut dc) {
        assert!((value.re - 2.0).abs() < 1e-12 && value.im.abs() < 1e-12);
    }
    // Edge lengths.
    let mut empty: Vec<Complex<f64>> = vec![];
    assert!(ifft(&mut empty).is_empty());
    let one = vec![Complex::new(3.0, -2.0)];
    assert_eq!(ifft(&mut one.clone()), one);
}

#[test]
fn blackman_window_shape_and_edges() {
    let w = blackman_window(5);
    assert_eq!(w.len(), 5);
    assert!(w[0].abs() < 1e-12);
    assert!(w[4].abs() < 1e-12);
    // 0.42 - 0.5 cos(pi/2) + 0.08 cos(pi) = 0.34 at i = 1, symmetric at i = 3.
    assert!((w[1] - 0.34).abs() < 1e-12, "{w:?}");
    assert!((w[3] - 0.34).abs() < 1e-12, "{w:?}");
    assert!((w[2] - 1.0).abs() < 1e-12, "{w:?}");
    for (a, b) in w.iter().zip(w.iter().rev()) {
        assert!((a - b).abs() < 1e-12, "not symmetric: {w:?}");
    }
    assert!(blackman_window(0).is_empty());
    assert_eq!(blackman_window(1), vec![1.0]);
    // A longer window stays in [0, 1] up to the round-off of the cosine sum
    // (the exact Blackman window is non-negative on [0, n-1]).
    let w = blackman_window(64);
    assert!(w.iter().all(|v| (-1e-12..=1.0 + 1e-12).contains(v)));
    assert!(w[0].abs() < 1e-12 && w[63].abs() < 1e-12);
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
