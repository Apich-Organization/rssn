//! FFT / IFFT at power-of-two and arbitrary lengths.

use std::f64::consts::PI;
use std::time::Instant;

use num_complex::Complex;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::transforms::fft;
use rssn::kernels::transforms::fft_slice;
use rssn::kernels::transforms::ifft;
use rssn::kernels::transforms::ifft_slice;
use rssn::kernels::transforms::irfft;
use rssn::kernels::transforms::rfft;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn close(
    a: Complex<f64>,
    b: Complex<f64>,
) -> bool {
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
    let want = [
        Complex::new(2.0, 0.0),
        Complex::new(1.0, -1.0),
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 1.0),
    ];
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
    let mut data: Vec<Complex<f64>> = (1..=4)
        .map(|i| Complex::new(2.0 * i as f64 - 1.0, 2.0 * i as f64))
        .collect();
    let original = data.clone();
    fft(&mut data);
    ifft(&mut data);
    for (g, w) in data.iter().zip(&original) {
        assert!(close(*g, *w));
    }
}

#[test]
fn non_power_of_two_input_keeps_its_length() {
    let mut data = vec![Complex::new(1.0, 0.0); 3];
    fft(&mut data);
    assert_eq!(data.len(), 3);
    // A constant signal of length 3: [3, 0, 0], not the 4-point padded spectrum.
    assert!(close(data[0], Complex::new(3.0, 0.0)));
    assert!(data[1].norm() < 1e-12);
    assert!(data[2].norm() < 1e-12);
}

#[test]
fn slice_variants_match_vec_variants() {
    let src: Vec<Complex<f64>> = (0..8)
        .map(|i| Complex::new(i as f64, (i * i) as f64 * 0.1))
        .collect();
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

// -------------------------------------------------------------------------
// Arbitrary lengths (Bluestein chirp-z) against a direct-DFT oracle
// -------------------------------------------------------------------------

/// Naive `O(n^2)` DFT, the reference every fast path is checked against.
fn direct_dft(x: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = x.len();
    (0..n)
        .map(|k| {
            x.iter()
                .enumerate()
                .fold(Complex::new(0.0, 0.0), |acc, (j, xj)| {
                    let ang = -2.0 * PI * (j as f64) * (k as f64) / (n as f64);
                    acc + xj * Complex::from_polar(1.0, ang)
                })
        })
        .collect()
}

/// A deterministic complex signal with both parts of order one.
fn sample(n: usize) -> Vec<Complex<f64>> {
    (0..n)
        .map(|i| {
            let t = i as f64;
            Complex::new((0.3 * t).sin() + 0.5, (0.7 * t).cos() - 0.2)
        })
        .collect()
}

/// `1..=64`, the three primes, and the two highly composite lengths.
fn lengths() -> Vec<usize> {
    let mut v: Vec<usize> = (1..=64).collect();
    v.extend([97, 101, 127, 120, 210]);
    v
}

#[test]
fn arbitrary_lengths_round_trip_and_satisfy_parseval() {
    for n in lengths() {
        let original = sample(n);
        let mut data = original.clone();
        fft(&mut data);
        assert_eq!(data.len(), n, "length changed at n = {n}");
        let time_energy: f64 = original.iter().map(|c| c.norm_sqr()).sum();
        let freq_energy: f64 = data.iter().map(|c| c.norm_sqr()).sum::<f64>() / n as f64;
        assert!(
            (time_energy - freq_energy).abs() <= 1e-12 * time_energy.max(1.0),
            "Parseval at n = {n}: {time_energy} vs {freq_energy}"
        );
        ifft(&mut data);
        for (got, want) in data.iter().zip(&original) {
            assert!(
                (got - want).norm() <= 1e-12,
                "round trip at n = {n}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn arbitrary_lengths_match_the_direct_dft() {
    for n in lengths() {
        let x = sample(n);
        let want = direct_dft(&x);
        let scale = want.iter().map(|c| c.norm()).fold(1.0_f64, f64::max);
        let mut got = x.clone();
        fft(&mut got);
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!(
                (g - w).norm() <= 1e-12 * scale,
                "bin {k} at n = {n}: {g} vs {w}"
            );
        }
    }
}

#[test]
fn prime_lengths_match_the_direct_dft() {
    for n in [97usize, 101, 127] {
        let x = sample(n);
        let want = direct_dft(&x);
        let scale = want.iter().map(|c| c.norm()).fold(1.0_f64, f64::max);
        let mut got = x.clone();
        fft(&mut got);
        assert_eq!(got.len(), n);
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!(
                (g - w).norm() <= 1e-12 * scale,
                "prime bin {k} at n = {n}: {g} vs {w}"
            );
        }
        // The slice entry point must agree with the Vec one exactly.
        let mut via_slice = x.clone();
        fft_slice(&mut via_slice);
        for (a, b) in got.iter().zip(&via_slice) {
            assert_eq!(a.re.to_bits(), b.re.to_bits(), "slice re at n = {n}");
            assert_eq!(a.im.to_bits(), b.im.to_bits(), "slice im at n = {n}");
        }
        // And the inverse of that spectrum is the signal.
        ifft(&mut got);
        for (g, w) in got.iter().zip(&x) {
            assert!((g - w).norm() <= 1e-12, "prime inverse at n = {n}");
        }
    }
}

/// The radix-2 kernel exactly as it stood before the arbitrary-length path
/// was added: the power-of-two route must keep reproducing it bit for bit.
fn radix2_reference(data: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let mut data = data.to_vec();
    let n = data.len();
    if n <= 1 {
        return data;
    }
    let mut j = 0;
    for i in 0..n {
        if i < j {
            data.swap(i, j);
        }
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
    }
    let mut len = 2;
    while len <= n {
        let half_len = len / 2;
        let w_m = Complex::from_polar(1.0, -2.0 * PI / (len as f64));
        for i in (0..n).step_by(len) {
            let mut w = Complex::new(1.0, 0.0);
            for k in 0..half_len {
                let even_idx = i + k;
                let odd_idx = i + k + half_len;
                let t = w * data[odd_idx];
                let u = data[even_idx];
                data[even_idx] = u + t;
                data[odd_idx] = u - t;
                w *= w_m;
            }
        }
        len *= 2;
    }
    data
}

#[test]
fn power_of_two_results_are_bit_identical_to_the_radix2_kernel() {
    for n in [1usize, 2, 4, 8, 16, 64, 128, 1024] {
        let x = sample(n);
        let want = radix2_reference(&x);
        let mut got = x.clone();
        fft(&mut got);
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g.re.to_bits(), w.re.to_bits(), "re bin {k} at n = {n}");
            assert_eq!(g.im.to_bits(), w.im.to_bits(), "im bin {k} at n = {n}");
        }
    }
}

#[test]
fn cost_of_the_general_path() {
    let m = 2048_usize;
    let log_m = 11_usize;
    println!(
        "n = 1000 (not a power of two): Bluestein, m = {m}, \
         3 transforms x {log_m} stages of {}/2 = {} butterflies",
        m,
        3 * m * log_m / 2
    );
    println!(
        "n = 1024 (power of two): 10 stages of 512 = {} butterflies",
        1024 * 10 / 2
    );

    let reps = 32_u32;
    for (label, n) in [("n = 1000", 1000usize), ("n = 1024", 1024)] {
        let pattern = sample(n);
        let mut data = pattern.clone();
        fft(&mut data);
        let start = Instant::now();
        for _ in 0..reps {
            data.copy_from_slice(&pattern);
            fft(&mut data);
        }
        let per = start.elapsed() / reps;
        println!("{label}: {per:?} per forward transform ({reps} reps)");
        assert!(data.iter().all(|c| c.re.is_finite() && c.im.is_finite()));
    }
}

// -------------------------------------------------------------------------
// Real-input transform pair (rfft / irfft)
// -------------------------------------------------------------------------

/// A deterministic real signal that excites every bin.
fn real_sample(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let t = i as f64;
            0.5 + (0.3 * t).sin() - 0.25 * (0.11 * t).cos() + 0.1 * ((i * i) as f64 * 0.01).sin()
        })
        .collect()
}

/// Even, odd, prime, highly composite, and the two smallest lengths.
fn real_lengths() -> Vec<usize> {
    vec![1, 2, 3, 4, 5, 8, 15, 16, 100, 1000]
}

#[test]
fn rfft_matches_the_complex_fft() {
    for n in real_lengths() {
        let x = real_sample(n);
        let half = rfft(&x);
        assert_eq!(half.len(), n / 2 + 1, "half-spectrum length at n = {n}");

        let mut full: Vec<Complex<f64>> = x.iter().map(|&value| Complex::new(value, 0.0)).collect();
        fft(&mut full);
        let scale = full.iter().map(|c| c.norm()).fold(1.0_f64, f64::max);

        for (k, (got, want)) in half.iter().zip(&full).enumerate() {
            assert!(
                (got - want).norm() <= 1e-12 * scale,
                "bin {k} at n = {n}: {got} vs {want}"
            );
        }
        // Bins above n/2 are the conjugate mirror of the half spectrum.
        for k in 1..=(n - 1) / 2 {
            let mirrored = half[k].conj();
            assert!(
                (full[n - k] - mirrored).norm() <= 1e-12 * scale,
                "mirror bin {k} at n = {n}: {} vs {mirrored}",
                full[n - k]
            );
        }
    }
}

#[test]
fn real_round_trip() {
    for n in real_lengths() {
        let x = real_sample(n);
        let back = irfft(&rfft(&x), n);
        assert_eq!(back.len(), n, "round-trip length at n = {n}");
        for (got, want) in back.iter().zip(&x) {
            assert!(
                (got - want).abs() <= 1e-12 * want.abs().max(1.0),
                "round trip at n = {n}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn real_parseval() {
    for n in real_lengths() {
        let x = real_sample(n);
        let spectrum = rfft(&x);
        let time_energy: f64 = x.iter().map(|value| value * value).sum();
        // The interior bins stand for their conjugate partner too; the DC and
        // (even n) Nyquist bins are not repeated.
        let last = spectrum.len() - 1;
        let mut freq_energy = spectrum[0].norm_sqr();
        for (k, bin) in spectrum.iter().enumerate().skip(1) {
            let weight = if n % 2 == 0 && k == last {
                1.0
            } else {
                2.0
            };
            freq_energy += weight * bin.norm_sqr();
        }
        let freq_energy = freq_energy / n as f64;
        assert!(
            (time_energy - freq_energy).abs() <= 1e-12 * time_energy.max(1.0),
            "Parseval at n = {n}: {time_energy} vs {freq_energy}"
        );
    }
}

#[test]
fn real_spectra_of_impulse_constant_and_sine() {
    let impulse = rfft(&[3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    for bin in &impulse {
        assert!((bin.re - 3.0).abs() < 1e-12, "{bin}");
        assert!(bin.im.abs() < 1e-12, "{bin}");
    }

    let constant = rfft(&[2.0; 8]);
    assert!((constant[0].re - 16.0).abs() < 1e-12);
    assert!(constant[1..].iter().all(|bin| bin.norm() < 1e-12));

    let n = 8_usize;
    let sine: Vec<f64> = (0..n)
        .map(|j| (2.0 * PI * j as f64 / n as f64).sin())
        .collect();
    let spectrum = rfft(&sine);
    // sin(2 pi j / 8) has X_1 = n/(2i) = -4i and X_7 = +4i.
    assert!((spectrum[1] - Complex::new(0.0, -4.0)).norm() < 1e-12);
    assert!(spectrum[0].norm() < 1e-12);
    assert!(spectrum[2].norm() < 1e-12);
    assert!(spectrum[3].norm() < 1e-12);
    assert!(spectrum[4].norm() < 1e-12);

    let cosine: Vec<f64> = (0..n)
        .map(|j| (2.0 * PI * 2.0 * j as f64 / n as f64).cos())
        .collect();
    let spectrum = rfft(&cosine);
    // cos(2 pi * 2 j / 8) has X_2 = n/2 = 4.
    assert!((spectrum[2] - Complex::new(4.0, 0.0)).norm() < 1e-12);
    assert!(spectrum[1].norm() < 1e-12);
    assert!(spectrum[3].norm() < 1e-12);
}

#[test]
fn real_prime_lengths_match_the_direct_dft() {
    for n in [3usize, 5, 101, 127] {
        let x = real_sample(n);
        let full: Vec<Complex<f64>> = x.iter().map(|&value| Complex::new(value, 0.0)).collect();
        let want = direct_dft(&full);
        let got = rfft(&x);
        let scale = want.iter().map(|c| c.norm()).fold(1.0_f64, f64::max);
        assert_eq!(got.len(), n / 2 + 1);
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!(
                (g - w).norm() <= 1e-12 * scale,
                "prime bin {k} at n = {n}: {g} vs {w}"
            );
        }
    }
}

#[test]
fn real_transform_edge_and_error_cases() {
    assert!(rfft(&[]).is_empty());
    assert_eq!(rfft(&[2.5]), vec![Complex::new(2.5, 0.0)]);
    assert_eq!(irfft(&[Complex::new(2.5, 0.0)], 1), vec![2.5]);

    // A half spectrum must have exactly n/2 + 1 bins.
    assert!(irfft(&[Complex::new(1.0, 0.0); 2], 8).is_empty());
    assert!(irfft(&[], 8).is_empty());
    assert!(irfft(&[Complex::new(1.0, 0.0)], 0).is_empty());

    // DC and Nyquist imaginary parts are ignored: a real signal has none.
    let x = [1.0, 2.0, 3.0, 4.0];
    let clean = irfft(&rfft(&x), 4);
    let mut noisy = rfft(&x);
    noisy[0].im = 0.5;
    noisy[2].im = -0.25;
    let from_noisy = irfft(&noisy, 4);
    for (a, b) in clean.iter().zip(&from_noisy) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }
    for (got, want) in from_noisy.iter().zip(&x) {
        assert!((got - want).abs() < 1e-12);
    }
}
