//! `kernels::fftconv`: the arbitrary-length FFT,
//! convolution, correlation and exact integer polynomial multiplication,
//! each against the textbook O(n²) definition.

use rssn::kernels::fftconv::circular_convolve;
use rssn::kernels::fftconv::convolve;
use rssn::kernels::fftconv::correlate;
use rssn::kernels::fftconv::fft;
use rssn::kernels::fftconv::poly_mul_i64;
use rssn::kernels::fftconv::poly_mul_i64_checked;
use rssn::kernels::random::Rng;
use std::f64::consts::PI;

type C = (f64, f64);

/// `X_k = Σ_j x_j e^{∓2πi jk/n}` by definition (`-` forward, `+` inverse,
/// the inverse carrying the `1/n`).
fn dft(
    x: &[C],
    inverse: bool,
) -> Vec<C> {
    let n = x.len();
    let sign = if inverse { 1.0 } else { -1.0 };
    (0..n)
        .map(|k| {
            let mut s = (0.0, 0.0);
            for (j, v) in x.iter().enumerate() {
                let ang = sign * 2.0 * PI * ((j * k) % n) as f64 / n as f64;
                s.0 += v.0 * ang.cos() - v.1 * ang.sin();
                s.1 += v.0 * ang.sin() + v.1 * ang.cos();
            }
            if inverse {
                (s.0 / n as f64, s.1 / n as f64)
            } else {
                s
            }
        })
        .collect()
}

fn random_signal(
    n: usize,
    seed: u64,
) -> Vec<C> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| (rng.normal(), rng.normal())).collect()
}

fn max_error(
    a: &[C],
    b: &[C],
) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| ((x.0 - y.0).powi(2) + (x.1 - y.1).powi(2)).sqrt())
        .fold(0.0, f64::max)
}

#[test]
fn fft_matches_the_definition_for_power_of_two_and_arbitrary_lengths() {
    // Cover power-of-two, prime and awkward composite lengths.
    for n in [1_usize, 2, 3, 5, 8, 12, 13, 64, 97, 100] {
        let x = random_signal(n, 100 + n as u64);
        let want = dft(&x, false);
        let got = fft(&x, false);
        assert!(
            max_error(&got, &want) < 1e-9 * (n as f64).max(1.0),
            "forward n = {n}: {}",
            max_error(&got, &want)
        );
        let want = dft(&x, true);
        let got = fft(&x, true);
        assert!(
            max_error(&got, &want) < 1e-9,
            "inverse n = {n}: {}",
            max_error(&got, &want)
        );
    }
}

#[test]
fn fft_round_trip_parseval_and_known_spectra() {
    for n in [7_usize, 16, 30, 127, 256, 1000] {
        let x = random_signal(n, 7 + n as u64);
        let back = fft(&fft(&x, false), true);
        assert!(
            max_error(&back, &x) < 1e-10,
            "round trip n = {n}: {}",
            max_error(&back, &x)
        );
        // Parseval: Σ|X_k|² = n Σ|x_j|².
        let spectrum = fft(&x, false);
        let energy: f64 = x.iter().map(|v| v.0 * v.0 + v.1 * v.1).sum();
        let spectral: f64 = spectrum.iter().map(|v| v.0 * v.0 + v.1 * v.1).sum();
        assert!(
            (spectral - n as f64 * energy).abs() < 1e-9 * spectral,
            "Parseval n = {n}"
        );
    }
    // The impulse has a flat spectrum, the constant a single line at zero.
    for n in [6_usize, 8, 11] {
        let mut impulse = vec![(0.0, 0.0); n];
        impulse[0] = (1.0, 0.0);
        for v in fft(&impulse, false) {
            assert!((v.0 - 1.0).abs() < 1e-12 && v.1.abs() < 1e-12);
        }
        let flat = fft(&vec![(1.0, 0.0); n], false);
        assert!((flat[0].0 - n as f64).abs() < 1e-12);
        for v in &flat[1..] {
            assert!(v.0.abs() < 1e-12 && v.1.abs() < 1e-12, "n = {n}: {v:?}");
        }
    }
    // A pure harmonic e^{2πi 3 j/10} (length 10) has one line at k = 3.
    let harmonic: Vec<C> = (0..10)
        .map(|j| {
            (
                (2.0 * PI * 3.0 * j as f64 / 10.0).cos(),
                (2.0 * PI * 3.0 * j as f64 / 10.0).sin(),
            )
        })
        .collect();
    let spectrum = fft(&harmonic, false);
    for (k, v) in spectrum.iter().enumerate() {
        let want = if k == 3 { 10.0 } else { 0.0 };
        assert!(
            (v.0 - want).abs() < 1e-10 && v.1.abs() < 1e-10,
            "k = {k}: {v:?}"
        );
    }
    // Lengths 0 and 1 are returned unchanged.
    assert!(fft(&[], false).is_empty());
    assert_eq!(fft(&[(2.5, -1.0)], true), vec![(2.5, -1.0)]);
}

fn naive_convolution(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

#[test]
fn convolution_matches_the_definition_on_both_code_paths() {
    let mut rng = Rng::new(5);
    // (3, 5) and (16, 40) use the direct loop (shorter factor ≤ 16), the
    // others the FFT; 17 is the first length that switches.
    for (la, lb) in [
        (1_usize, 1_usize),
        (3, 5),
        (16, 40),
        (17, 17),
        (17, 100),
        (64, 63),
        (200, 301),
    ] {
        let a: Vec<f64> = (0..la).map(|_| rng.normal()).collect();
        let b: Vec<f64> = (0..lb).map(|_| rng.normal()).collect();
        let got = convolve(&a, &b);
        let want = naive_convolution(&a, &b);
        assert_eq!(got.len(), la + lb - 1);
        let err = got
            .iter()
            .zip(&want)
            .map(|(g, w)| (g - w).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-10 * (la.max(lb) as f64), "({la}, {lb}): {err}");
    }
    assert!(convolve(&[], &[1.0]).is_empty() && convolve(&[1.0], &[]).is_empty());
    // Commutative, and a delta is the identity.
    let a: Vec<f64> = (0..40).map(|_| rng.normal()).collect();
    let b: Vec<f64> = (0..33).map(|_| rng.normal()).collect();
    let (ab, ba) = (convolve(&a, &b), convolve(&b, &a));
    assert!(ab.iter().zip(&ba).all(|(x, y)| (x - y).abs() < 1e-11));
    let same = convolve(&a, &[1.0]);
    assert_eq!(same, a);
}

#[test]
fn circular_convolution_wraps_around() {
    let mut rng = Rng::new(9);
    for n in [1_usize, 4, 7, 12, 31, 64] {
        let a: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        let b: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        let want: Vec<f64> = (0..n)
            .map(|k| (0..n).map(|j| a[j] * b[(k + n - j) % n]).sum())
            .collect();
        let got = circular_convolve(&a, &b);
        let err = got
            .iter()
            .zip(&want)
            .map(|(g, w)| (g - w).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-10, "n = {n}: {err}");
    }
}

#[test]
fn correlation_is_convolution_with_the_reversed_sequence() {
    let mut rng = Rng::new(11);
    for (la, lb) in [(5_usize, 3_usize), (30, 20), (64, 40)] {
        let a: Vec<f64> = (0..la).map(|_| rng.normal()).collect();
        let b: Vec<f64> = (0..lb).map(|_| rng.normal()).collect();
        let got = correlate(&a, &b);
        assert_eq!(got.len(), la + lb - 1);
        // c[k] = Σ_i a[i + k - (lb - 1)] b[i], zero outside the range of a.
        for k in 0..got.len() {
            let mut want = 0.0;
            for i in 0..lb {
                let idx = (i + k) as isize - (lb as isize - 1);
                if (0..la as isize).contains(&idx) {
                    want += a[idx as usize] * b[i];
                }
            }
            assert!(
                (got[k] - want).abs() < 1e-10,
                "({la}, {lb}) k = {k}: {} vs {want}",
                got[k]
            );
        }
    }
    // The autocorrelation peaks at zero lag with the energy.
    let a: Vec<f64> = (0..50).map(|_| rng.normal()).collect();
    let auto = correlate(&a, &a);
    let peak = (0..auto.len())
        .max_by(|&i, &j| auto[i].total_cmp(&auto[j]))
        .unwrap();
    assert_eq!(peak, a.len() - 1);
    assert!((auto[peak] - a.iter().map(|v| v * v).sum::<f64>()).abs() < 1e-10);
}

#[test]
fn integer_polynomial_products_are_exact() {
    // (1 + x)^2 and a Pascal row: (1 + x)^10 = (1 + x)^4 (1 + x)^6.
    assert_eq!(poly_mul_i64(&[1, 1], &[1, 1]), vec![1, 2, 1]);
    let binom = |n: usize| -> Vec<i64> {
        let mut row = vec![1_i64];
        for _ in 0..n {
            let mut next = vec![1_i64; row.len() + 1];
            for i in 1..row.len() {
                next[i] = row[i - 1] + row[i];
            }
            row = next;
        }
        row
    };
    assert_eq!(poly_mul_i64(&binom(4), &binom(6)), binom(10));
    // Signed coefficients of 24 bits (two 15-bit limbs each) against exact
    // 128-bit arithmetic, long enough to take the FFT path.
    let mut rng = Rng::new(21);
    for (la, lb) in [(1_usize, 1_usize), (7, 9), (60, 80), (130, 200)] {
        let draw = |rng: &mut Rng, n: usize| -> Vec<i64> {
            (0..n)
                .map(|_| (rng.below(1 << 24) as i64) - (1 << 23))
                .collect()
        };
        let a = draw(&mut rng, la);
        let b = draw(&mut rng, lb);
        let got = poly_mul_i64(&a, &b);
        let mut want = vec![0_i128; la + lb - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                want[i + j] += i128::from(*x) * i128::from(*y);
            }
        }
        assert_eq!(got.len(), la + lb - 1);
        for (k, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(i128::from(*g), *w, "({la}, {lb}) coefficient {k}");
        }
    }
    // Small coefficients: everything is exact and fits.
    let a: Vec<i64> = (0..300).map(|i| (i % 7) - 3).collect();
    let b: Vec<i64> = (0..257).map(|i| 5 - (i % 11)).collect();
    let got = poly_mul_i64(&a, &b);
    for k in 0..got.len() {
        let mut want = 0_i64;
        for i in 0..a.len() {
            if k >= i && k - i < b.len() {
                want += a[i] * b[k - i];
            }
        }
        assert_eq!(got[k], want, "coefficient {k}");
    }
    assert!(poly_mul_i64(&[], &[1]).is_empty());
}

#[test]
fn full_width_integer_convolution_is_exact_and_deterministic() {
    assert_eq!(
        poly_mul_i64_checked(&[i64::MIN], &[1]),
        Some(vec![i64::MIN])
    );
    assert_eq!(poly_mul_i64_checked(&[i64::MIN], &[-1]), None);
    assert_eq!(poly_mul_i64(&[i64::MIN], &[-1]), vec![i64::MIN]);
    let mut wide = vec![0; 35];
    wide[20] = i64::MIN;
    let narrow: Vec<_> = (0..27).map(|i| i64::from(i % 2)).collect();
    let mut product = vec![0; 61];
    for (dest, &factor) in product[20..47].iter_mut().zip(&narrow) {
        *dest = i64::MIN * factor;
    }
    assert_eq!(poly_mul_i64_checked(&wide, &narrow), Some(product.clone()));
    assert_eq!(poly_mul_i64(&wide, &narrow), product);

    let mut rng = Rng::new(3418);
    for (m, n) in [(18, 18), (50, 37), (300, 250)] {
        let a: Vec<_> = (0..m).map(|_| rng.next_u64() as i64).collect();
        let b: Vec<_> = (0..n).map(|_| rng.next_u64() as i64).collect();
        let mut expected = vec![0_i64; m + n - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                expected[i + j] = expected[i + j].wrapping_add(x.wrapping_mul(y));
            }
        }
        assert_eq!(poly_mul_i64(&a, &b), expected);
        assert!(poly_mul_i64_checked(&a, &b).is_none());
    }
}

#[test]
fn fft_axes_batch_strided_and_empty_shapes() {
    use rssn::kernels::fft::Complex64;
    use rssn::kernels::fft::Norm;
    use rssn::kernels::fft::fft_axes;
    use rssn::kernels::fft::fftfreq;
    use rssn::kernels::fft::irfft_axes;
    use rssn::kernels::fft::rfft_axes;
    use rssn::kernels::fft::rfftfreq;
    let samples: Vec<_> = (0..12).map(|i| (f64::from(i) / 2.0).sin()).collect();
    let complex: Vec<_> = samples.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    for axes in [[0, 1], [1, 0]] {
        let f = fft_axes(&complex, &[3, 4], &axes, false, Norm::Unitary).unwrap();
        let recovered = fft_axes(&f, &[3, 4], &axes, true, Norm::Unitary).unwrap();
        for (actual, expected) in recovered.iter().zip(&complex) {
            assert!((actual - expected).norm() < 1e-13);
        }

        let (f, shape) = rfft_axes(&samples, &[3, 4], &axes).unwrap();
        let n = [3, 4][axes[1]];
        let (recovered, full_shape) = irfft_axes(&f, &shape, &axes, n).unwrap();
        assert_eq!(full_shape, [3, 4]);
        for (actual, expected) in recovered.iter().zip(&samples) {
            assert!((actual - expected).abs() < 1e-13);
        }
    }
    let empty_shape = [usize::MAX, 2, 0];
    let (f, shape) = rfft_axes(&[], &empty_shape, &[1]).unwrap();
    assert!(f.is_empty());
    assert_eq!(
        irfft_axes(&f, &shape, &[1], 2).unwrap(),
        (vec![], empty_shape.to_vec())
    );
    assert!(rfftfreq(4, f64::MAX)[1] > 0.0);
    assert!(fftfreq(4, f64::MAX)[2] < 0.0);
    assert!(fftfreq(usize::MAX, 1.0).is_empty());
    assert!(rfftfreq(usize::MAX, 1.0).is_empty());
}
