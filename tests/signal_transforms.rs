//! Independent references for the FFT-backed signal transforms and filter banks.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;
use rssn::kernels::fft;
use rssn::kernels::signal;
use rssn::kernels::signal::SignalError;
use rssn::kernels::wavelets;
use rssn::kernels::wavelets::FrftPlan;
use rssn::kernels::wavelets::Wavelet;
use rssn::kernels::wavelets::WaveletError;

fn close(
    got: f64,
    want: f64,
) {
    assert!(
        (got - want).abs() <= 5e-11 * want.abs().max(1.0),
        "{got} vs {want}"
    );
}

fn close_complex(
    got: Complex64,
    want: Complex64,
) {
    close(got.re, want.re);
    close(got.im, want.im);
}

#[test]
fn cosine_and_sine_transforms_match_direct_orthonormal_sums() {
    for n in [0, 1, 2, 3, 8, 15, 16, 17, 64, 97] {
        let values: Vec<_> = (0..n)
            .map(|i| (f64::from(i) * 0.73).cos() + f64::from(i % 3))
            .collect();
        let cosine = signal::dct(&values);
        let sine = signal::dst(&values);
        let scale = (2.0 / f64::from(n)).sqrt();
        for (k, (&cosine, &sine)) in cosine.iter().zip(&sine).enumerate() {
            let mut expected_cosine = 0.0;
            let mut expected_sine = 0.0;
            for (i, &value) in values.iter().enumerate() {
                let phase = PI * (i as f64 + 0.5) / f64::from(n);
                expected_cosine += value * (phase * k as f64).cos();
                expected_sine += value * (phase * (k + 1) as f64).sin();
            }
            let dc = if k == 0 {
                2.0_f64.sqrt()
            } else {
                1.0
            };
            let nyquist = if k + 1 == values.len() {
                2.0_f64.sqrt()
            } else {
                1.0
            };
            close(cosine, expected_cosine * scale / dc);
            close(sine, expected_sine * scale / nyquist);
        }
        for result in [signal::idct(&cosine), signal::idst(&sine)] {
            assert_eq!(result.len(), values.len());
            for (&got, &want) in result.iter().zip(&values) {
                close(got, want);
            }
        }
    }
}

#[test]
fn analytic_signal_preserves_dc_and_nyquist_and_rotates_positive_frequencies() {
    for n in [1_u32, 2, 5, 12, 17] {
        let (signal, expected): (Vec<_>, Vec<_>) = (0..n)
            .map(|i| {
                let phase = TAU * f64::from(i) / f64::from(n);
                let nyquist = if n.is_multiple_of(2) && i.is_multiple_of(2) {
                    1.5
                } else if n.is_multiple_of(2) {
                    -1.5
                } else {
                    0.0
                };
                (
                    2.0 + phase.cos() + nyquist,
                    if n <= 2 { 0.0 } else { phase.sin() },
                )
            })
            .unzip();
        let analytic = signal::analytic_signal(&signal).unwrap();
        let hilbert = signal::hilbert(&signal).unwrap();
        let envelope = signal::analytic_envelope(&signal).unwrap();
        for (index, value) in analytic.iter().enumerate() {
            assert_eq!(value.re, signal[index]);
            close(value.im, expected[index]);
            close(hilbert[index], expected[index]);
            close(envelope[index], signal[index].hypot(expected[index]));
        }
    }
    assert!(signal::analytic_signal(&[]).unwrap().is_empty());
    assert_eq!(signal::hilbert(&[f64::MAX; 4]).unwrap(), [0.0; 4]);
    let tiny = f64::from_bits(1);
    assert_eq!(
        signal::hilbert(&[tiny, 0.0, -tiny, 0.0]).unwrap(),
        [0.0, tiny, 0.0, -tiny]
    );
    assert_eq!(
        signal::hilbert(&[f64::NAN]),
        Err(SignalError::NonFiniteInput)
    );
}

#[test]
fn zoom_fft_matches_direct_complex_sums_at_arbitrary_frequencies() {
    let samples: Vec<_> = (0..37)
        .map(|n| {
            let value = f64::from(n);
            Complex64::new((0.63 * value).cos(), (0.77 * value).sin())
        })
        .collect();
    for (start, step, bins) in [
        (0.21, 0.0035, 71),
        (-0.34, -0.071, 13),
        (1.43, 0.0, 23),
        (-1.0, 0.5, 1),
    ] {
        let actual = signal::zoom_fft(&samples, start, step, bins).unwrap();
        assert_eq!(actual.len(), bins);
        for (k, got) in actual.into_iter().enumerate() {
            let frequency = start + k as f64 * step;
            let want = samples
                .iter()
                .enumerate()
                .map(|(i, &sample)| {
                    sample * Complex64::from_polar(1.0, -TAU * frequency * i as f64)
                })
                .sum();
            close_complex(got, want);
        }
    }
    let regular = fft::fft(&samples);
    let zoomed =
        signal::zoom_fft(&samples, 0.0, 1.0 / samples.len() as f64, samples.len()).unwrap();
    for (actual, expected) in zoomed.into_iter().zip(regular) {
        close_complex(actual, expected);
    }
    assert_eq!(
        signal::zoom_fft(&[], 0.0, 0.0, 3).unwrap(),
        vec![Complex64::ZERO; 3]
    );
    assert!(signal::zoom_fft(&samples, 0.0, 1.0, 0).unwrap().is_empty());
    assert_eq!(
        signal::zoom_fft(&samples, f64::INFINITY, 0.4, 1),
        Err(SignalError::InvalidFrequency)
    );
    assert_eq!(
        signal::zoom_fft(&samples, 0.0, 0.4, usize::MAX),
        Err(SignalError::OutputSize)
    );
}

#[test]
fn paired_welch_spectra_keep_phase_and_cancel_incoherent_frames() {
    let x = [1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
    let sin = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0];
    let switched = [0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 1.0];
    let window = [1.0; 4];
    let result = signal::welch_cross_spectrum(&x, &sin, &window, 4, 4.0).unwrap();
    assert_eq!(result.frequencies, [0.0, 1.0, 2.0]);
    close_complex(result.density[1], Complex64::new(0.0, -0.5));
    close(result.coherence[1], 1.0);
    for index in [0, 2] {
        close_complex(result.density[index], Complex64::ZERO);
        close(result.coherence[index], 0.0);
    }
    let result = signal::welch_cross_spectrum(&x, &switched, &window, 4, 4.0).unwrap();
    close_complex(result.density[1], Complex64::ZERO);
    close(result.coherence[1], 0.0);
    let auto = signal::welch_cross_spectrum(&x, &x, &window, 1, 9.0).unwrap();
    let power = signal::welch_psd(&x, &window, 1, 9.0).unwrap();
    for (&cross, power) in auto.density.iter().zip(power.density) {
        close_complex(cross, Complex64::new(power, 0.0));
    }
    assert_eq!(
        signal::welch_cross_spectrum(&x, &[], &window, 4, 4.0),
        Err(SignalError::SampleShape)
    );
    assert_eq!(
        signal::welch_cross_spectrum(&[f64::NAN, 0.0], &[0.0; 2], &[1.0], 1, 1.0),
        Err(SignalError::NonFiniteInput)
    );
}

#[test]
fn real_tensor_fft_matches_the_full_complex_transform_for_mixed_axes() {
    let shape = [3, 4, 5];
    let samples: Vec<_> = (0..60)
        .map(|i| (f64::from(i) * 0.7).cos() + f64::from(i % 3))
        .collect();
    let complex: Vec<_> = samples.iter().map(|&x| Complex64::new(x, 0.0)).collect();
    for axes in [&[0, 1, 2][..], &[2, 0], &[1], &[0, 2, 1], &[2, 1, 0]] {
        let last = axes[axes.len() - 1];
        let full = fft::fft_axes(&complex, &shape, axes, false, fft::Norm::Backward).unwrap();
        let (half, spectrum_shape) = fft::rfft_axes(&samples, &shape, axes).unwrap();
        let stride: usize = shape[last + 1..].iter().product();
        let expected: Vec<_> = full
            .into_iter()
            .enumerate()
            .filter_map(|(i, bin)| {
                ((i / stride) % shape[last] < spectrum_shape[last]).then_some(bin)
            })
            .collect();
        assert_eq!(half.len(), expected.len());
        for (&got, want) in half.iter().zip(expected) {
            close_complex(got, want);
        }
        let (back, recovered_shape) =
            fft::irfft_axes(&half, &spectrum_shape, axes, shape[last]).unwrap();
        assert_eq!(recovered_shape, shape);
        assert_eq!(back.len(), samples.len());
        for (got, &want) in back.into_iter().zip(&samples) {
            close(got, want);
        }
    }
    let (empty, dims) = fft::rfft_axes(&[], &[0, usize::MAX, 3], &[2]).unwrap();
    assert!(empty.is_empty());
    assert_eq!(dims, [0, usize::MAX, 2]);
    assert_eq!(
        fft::irfft_axes(&empty, &dims, &[2], 3).unwrap(),
        (Vec::new(), vec![0, usize::MAX, 3])
    );
    assert!(fft::rfft_axes(&samples, &shape, &[0, 0]).is_err());
    assert!(fft::rfft_axes(&samples, &shape, &[]).is_err());
}

#[test]
fn fft_batches_are_bit_identical_across_thread_counts() {
    let shape = [16, 64, 257];
    let signal: Vec<_> = (0..(shape.iter().product::<usize>() + 1))
        .map(|i| (i as f64 * 0.031).sin())
        .collect();
    let window = signal::hamming_window(257);
    let compute = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let fft = fft::rfft_axes(&signal[1..], &shape, &[2, 0, 1]).unwrap();
                let stft = signal::stft(&signal, &window, 127);
                (fft, stft)
            })
    };
    assert_eq!(compute(1), compute(4));
}

#[test]
fn packet_bands_follow_tree_order_and_reconstruct_all_families() {
    let coefficients = wavelets::wavelet_packet(&[1.0, 2.0, 3.0, 4.0], Wavelet::Haar, 2).unwrap();
    for (band, value) in coefficients.iter().zip([5.0, -2.0, -1.0, 0.0]) {
        assert_eq!(band.len(), 1);
        close(band[0], value);
    }
    for wavelet in Wavelet::ALL {
        let samples: Vec<_> = (0..96).map(|i| (f64::from(i) * 0.55).cos()).collect();
        let leaf = wavelets::wavelet_packet(&samples, wavelet, 5).unwrap();
        assert_eq!(leaf.len(), 32);
        assert!(leaf.iter().all(|band| band.len() == 3));
        let back = wavelets::inverse_wavelet_packet(&leaf, wavelet).unwrap();
        for (got, &want) in back.into_iter().zip(&samples) {
            close(got, want);
        }
    }
    assert_eq!(
        wavelets::wavelet_packet(&[7.0], Wavelet::Haar, 0).unwrap(),
        [vec![7.0]]
    );
    assert_eq!(
        wavelets::wavelet_packet(&[1.0; 3], Wavelet::Haar, 1),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavelets::wavedec(&[1.0; 8], Wavelet::Haar, usize::BITS as usize),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavelets::inverse_wavelet_packet(&[vec![1.0], vec![0.0; 2]], Wavelet::Haar),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        wavelets::inverse_wavelet_packet(&[vec![0.0]; 0], Wavelet::Haar),
        Err(WaveletError::Shape)
    );
}

#[test]
fn reusable_fractional_fourier_plan_matches_matrix_fft_and_group_law() {
    for n in [2, 3, 4, 6, 11, 12, 17] {
        let signal: Vec<_> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.34).cos(), (i as f64 * 0.62).sin()))
            .collect();
        let plan = FrftPlan::new(n).unwrap();
        let a = -0.625;
        let b = 1.625;
        let direct = fft::fft_norm(&signal, fft::Norm::Unitary);
        let transformed = plan.transform(&signal, a).unwrap();
        let matrix = plan.matrix(a).unwrap();
        for (got, row) in transformed.iter().zip(matrix.chunks_exact(n)) {
            let expected = signal.iter().zip(row).map(|(x, m)| x * m).sum();
            close_complex(*got, expected);
        }
        let composed = plan.transform(&transformed, b).unwrap();
        for (got, want) in composed.into_iter().zip(direct) {
            close_complex(got, want);
        }
        let periodic = plan.transform(&signal, a + 400.0).unwrap();
        for (got, want) in periodic.into_iter().zip(transformed) {
            close_complex(got, want);
        }
        assert!(plan.transform(&signal[..1], 0.1).is_err());
        assert!(plan.matrix(f64::NAN).is_err());
        assert_eq!(wavelets::frft(&signal, -4.0).unwrap(), signal);
    }
    assert!(FrftPlan::new(1).is_err());
    assert!(FrftPlan::new(usize::MAX).is_err());
}
