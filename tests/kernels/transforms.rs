//! Radix-2 FFT / IFFT (ported from `numerical_transforms_test.rs`).

use num_complex::Complex;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::transforms::{fft, fft_slice, ifft, ifft_slice};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn close(a: Complex<f64>, b: Complex<f64>) -> bool {
    (a - b).norm() < 1e-9
}

#[test]
fn fft_known_spectrum() {
    let mut data = vec![
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];
    fft(&mut data);
    // [2, 1 - i, 0, 1 + i]
    let want = [Complex::new(2.0, 0.0), Complex::new(1.0, -1.0), Complex::new(0.0, 0.0), Complex::new(1.0, 1.0)];
    for (g, w) in data.iter().zip(want) {
        assert!(close(*g, w), "{g} vs {w}");
    }
}

#[test]
fn ifft_known_signal() {
    let mut data = vec![
        Complex::new(2.0, 0.0),
        Complex::new(1.0, -1.0),
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 1.0),
    ];
    ifft(&mut data);
    let want = [1.0, 1.0, 0.0, 0.0];
    for (g, w) in data.iter().zip(want) {
        assert!(close(*g, Complex::new(w, 0.0)), "{g} vs {w}");
    }
}

#[test]
fn round_trip() {
    let mut data: Vec<Complex<f64>> = (1..=4).map(|i| Complex::new(2.0 * i as f64 - 1.0, 2.0 * i as f64)).collect();
    let original = data.clone();
    fft(&mut data);
    ifft(&mut data);
    for (g, w) in data.iter().zip(&original) {
        assert!(close(*g, *w));
    }
}

#[test]
fn non_power_of_two_input_is_zero_padded() {
    let mut data = vec![Complex::new(1.0, 0.0); 3];
    fft(&mut data);
    assert_eq!(data.len(), 4);
    // DC bin is the sum of the samples.
    assert!(close(data[0], Complex::new(3.0, 0.0)));
}

#[test]
fn slice_variants_match_vec_variants() {
    let src: Vec<Complex<f64>> = (0..8).map(|i| Complex::new(i as f64, (i * i) as f64 * 0.1)).collect();
    let mut a = src.clone();
    let mut b = src.clone();
    fft(&mut a);
    fft_slice(&mut b);
    for (x, y) in a.iter().zip(&b) {
        assert!(close(*x, *y));
    }
    ifft_slice(&mut b);
    for (x, y) in b.iter().zip(&src) {
        assert!(close(*x, *y));
    }
}

#[test]
fn empty_and_single_element_are_untouched() {
    let mut empty: Vec<Complex<f64>> = vec![];
    fft(&mut empty);
    ifft(&mut empty);
    assert!(empty.is_empty());
    let mut one = vec![Complex::new(3.0, -2.0)];
    fft(&mut one);
    assert_eq!(one, vec![Complex::new(3.0, -2.0)]);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_round_trip(re in prop::collection::vec(-100.0..100.0f64, 16), im in prop::collection::vec(-100.0..100.0f64, 16)) {
        let mut data: Vec<Complex<f64>> = re.iter().zip(&im).map(|(&r, &i)| Complex::new(r, i)).collect();
        let original = data.clone();
        fft(&mut data);
        ifft(&mut data);
        for (g, w) in data.iter().zip(&original) {
            prop_assert!((g - w).norm() < 1e-7);
        }
    }

    #[test]
    fn prop_fft_is_linear(a in prop::collection::vec(-10.0..10.0f64, 8), b in prop::collection::vec(-10.0..10.0f64, 8), s in -3.0..3.0f64) {
        let to_c = |v: &[f64]| -> Vec<Complex<f64>> { v.iter().map(|&x| Complex::new(x, 0.0)).collect() };
        let (mut fa, mut fb) = (to_c(&a), to_c(&b));
        let mut fsum: Vec<Complex<f64>> = a.iter().zip(&b).map(|(x, y)| Complex::new(s * x + y, 0.0)).collect();
        fft(&mut fa);
        fft(&mut fb);
        fft(&mut fsum);
        for k in 0..8 {
            prop_assert!((fsum[k] - (fa[k] * s + fb[k])).norm() < 1e-8);
        }
    }

    #[test]
    fn prop_parseval(re in prop::collection::vec(-10.0..10.0f64, 16)) {
        let mut data: Vec<Complex<f64>> = re.iter().map(|&r| Complex::new(r, 0.0)).collect();
        let t: f64 = data.iter().map(|c| c.norm_sqr()).sum();
        fft(&mut data);
        let f: f64 = data.iter().map(|c| c.norm_sqr()).sum::<f64>() / 16.0;
        prop_assert!((t - f).abs() < 1e-8 * t.max(1.0));
    }
}
