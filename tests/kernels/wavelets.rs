//! Discrete wavelet transforms and the discrete fractional Fourier transform.

use num_complex::Complex64;
use rssn::kernels::wavelets::Wavelet;
use rssn::kernels::wavelets::WaveletError;
use rssn::kernels::wavelets::dwt;
use rssn::kernels::wavelets::frft;
use rssn::kernels::wavelets::frft_matrix;
use rssn::kernels::wavelets::idwt;
use rssn::kernels::wavelets::wavedec;
use rssn::kernels::wavelets::waverec;

const FAMILIES: [Wavelet; 4] = [Wavelet::Haar, Wavelet::Db2, Wavelet::Db3, Wavelet::Db4];

fn signal(n: usize) -> Vec<f64> {
    let mut s = 12345_u64;
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        })
        .collect()
}

#[test]
fn filters_are_orthonormal_with_the_advertised_vanishing_moments() {
    for (w, moments) in FAMILIES.iter().zip([1, 2, 3, 4]) {
        let (h, g) = (w.low_pass(), w.high_pass());
        let sum: f64 = h.iter().sum();
        assert!((sum - 2.0_f64.sqrt()).abs() < 1e-12, "{w:?}: sum h = {sum}");
        for shift in 0..h.len() / 2 {
            let dot: f64 = (0..h.len())
                .filter(|m| m + 2 * shift < h.len())
                .map(|m| h[m] * h[m + 2 * shift])
                .sum();
            let want = if shift == 0 { 1.0 } else { 0.0 };
            assert!(
                (dot - want).abs() < 1e-10,
                "{w:?}: double-shift orthogonality {shift}: {dot}"
            );
            let cross: f64 = (0..h.len())
                .filter(|m| m + 2 * shift < h.len())
                .map(|m| g[m] * h[m + 2 * shift])
                .sum();
            assert!(cross.abs() < 1e-10, "{w:?}: g not orthogonal to h: {cross}");
        }
        for p in 0..moments {
            let mom: f64 = g
                .iter()
                .enumerate()
                .map(|(m, v)| v * (m as f64).powi(p))
                .sum();
            assert!(mom.abs() < 1e-8, "{w:?}: moment {p} of g is {mom}");
        }
    }
}

#[test]
fn one_level_is_an_orthogonal_transform_with_exact_reconstruction() {
    for w in FAMILIES {
        let x = signal(64);
        let (a, d) = dwt(&x, w).unwrap();
        assert_eq!((a.len(), d.len()), (32, 32));
        let energy: f64 = x.iter().map(|v| v * v).sum();
        let coeff: f64 = a.iter().chain(&d).map(|v| v * v).sum();
        assert!((energy - coeff).abs() < 1e-12 * energy, "{w:?}: Parseval");
        let back = idwt(&a, &d, w).unwrap();
        let err = x
            .iter()
            .zip(&back)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-13, "{w:?}: reconstruction {err}");
    }
    // a constant has no detail; an interior linear ramp has none for db2 and up
    let (_, d) = dwt(&[3.0; 16], Wavelet::Db2).unwrap();
    assert!(d.iter().all(|v| v.abs() < 1e-12));
    let ramp: Vec<f64> = (0..32).map(f64::from).collect();
    let (_, d) = dwt(&ramp, Wavelet::Db2).unwrap();
    assert!(d[..d.len() - 2].iter().all(|v| v.abs() < 1e-10), "{d:?}");
}

#[test]
fn multilevel_decomposition_round_trips() {
    for w in FAMILIES {
        let x = signal(128);
        let c = wavedec(&x, w, 4).unwrap();
        assert_eq!(c.len(), 5);
        assert_eq!(c[0].len(), 8);
        assert_eq!(c[1].len(), 8);
        assert_eq!(c[4].len(), 64);
        let back = waverec(&c, w).unwrap();
        assert!(
            x.iter().zip(&back).all(|(p, q)| (p - q).abs() < 1e-12),
            "{w:?}"
        );
    }
}

#[test]
fn invalid_inputs_are_errors() {
    assert_eq!(
        dwt(&[1.0, 2.0, 3.0], Wavelet::Haar),
        Err(WaveletError::OddLength)
    );
    assert_eq!(dwt(&[], Wavelet::Haar), Err(WaveletError::OddLength));
    assert_eq!(
        idwt(&[1.0], &[1.0, 2.0], Wavelet::Haar),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        wavedec(&signal(12), Wavelet::Haar, 3),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavedec(&signal(16), Wavelet::Db4, 3),
        Err(WaveletError::TooManyLevels)
    );
    assert!(waverec(&[], Wavelet::Haar).is_err());
}

fn cplx(n: usize) -> Vec<Complex64> {
    let (a, b) = (signal(n), signal(n + 7));
    a.iter()
        .zip(&b[7..])
        .map(|(p, q)| Complex64::new(*p, *q))
        .collect()
}

fn max_diff(
    a: &[Complex64],
    b: &[Complex64],
) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| (p - q).norm())
        .fold(0.0, f64::max)
}

#[test]
fn the_fractional_fourier_transform_is_a_unitary_group_through_dft_and_reversal() {
    for n in [8, 15, 16, 24] {
        let x = cplx(n);
        // a = 0: identity; a = 1: unitary DFT; a = 2: reversal; a = 4: identity again
        assert!(
            max_diff(&frft(&x, 0.0).unwrap(), &x) < 1e-9,
            "n = {n}: a = 0"
        );
        let dft: Vec<Complex64> = (0..n)
            .map(|j| {
                (0..n)
                    .map(|k| {
                        Complex64::from_polar(
                            1.0,
                            -2.0 * std::f64::consts::PI * (j * k) as f64 / n as f64,
                        ) * x[k]
                    })
                    .sum::<Complex64>()
                    / (n as f64).sqrt()
            })
            .collect();
        assert!(
            max_diff(&frft(&x, 1.0).unwrap(), &dft) < 1e-9,
            "n = {n}: a = 1"
        );
        let rev: Vec<Complex64> = (0..n).map(|j| x[(n - j) % n]).collect();
        assert!(
            max_diff(&frft(&x, 2.0).unwrap(), &rev) < 1e-9,
            "n = {n}: a = 2"
        );
        assert!(
            max_diff(&frft(&x, 4.0).unwrap(), &x) < 1e-9,
            "n = {n}: a = 4"
        );
        // index additivity and unitarity at fractional orders
        let step = frft(&frft(&x, 0.3).unwrap(), 0.45).unwrap();
        assert!(
            max_diff(&step, &frft(&x, 0.75).unwrap()) < 1e-9,
            "n = {n}: additivity"
        );
        let y = frft(&x, 0.37).unwrap();
        let (e0, e1): (f64, f64) = (
            x.iter().map(|v| v.norm_sqr()).sum(),
            y.iter().map(|v| v.norm_sqr()).sum(),
        );
        assert!((e0 - e1).abs() < 1e-9 * e0, "n = {n}: norm");
    }
}

#[test]
fn the_fractional_matrix_is_unitary_and_validates_the_size() {
    let n = 12;
    let m = frft_matrix(n, 0.6).unwrap();
    for i in 0..n {
        for j in 0..n {
            let dot: Complex64 = (0..n).map(|k| m[i * n + k] * m[j * n + k].conj()).sum();
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((dot - want).norm() < 1e-9, "({i}, {j}): {dot}");
        }
    }
    assert!(frft_matrix(1, 0.5).is_err());
}
