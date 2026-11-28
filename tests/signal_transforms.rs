//! Independent references for the FFT-backed signal transforms and filter banks.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;
use rssn::kernels::fft;
use rssn::kernels::fftconv;
use rssn::kernels::fftconv::OverlapSavePlan;
use rssn::kernels::signal;
use rssn::kernels::signal::FilterError;
use rssn::kernels::signal::SignalError;
use rssn::kernels::wavelets;
use rssn::kernels::wavelets::FrftPlan;
use rssn::kernels::wavelets::Wavelet;
use rssn::kernels::wavelets::WaveletError;
use rssn::kernels::wavelets::WaveletProcessingError;

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
fn zoom_fft_preserves_phase_and_normalization_across_cache_sizes() {
    for len in [257, 32_769, 65_536] {
        let index = len / 3;
        let pulse = Complex64::new(0.75, -0.25);
        let mut samples = vec![Complex64::ZERO; len];
        samples[index] = pulse;
        let step = -1.0 / 256.0;
        for start in [0.0, 1.0 / 32.0] {
            let spectrum = signal::zoom_fft(&samples, start, step, 7).unwrap();
            assert_eq!(spectrum.len(), 7);
            for (bin, actual) in spectrum.into_iter().enumerate() {
                let cycles = ((start + bin as f64 * step) * index as f64).rem_euclid(1.0);
                let expected = pulse * Complex64::from_polar(1.0, -TAU * cycles);
                close_complex(actual, expected);
            }
        }
    }
}

#[test]
fn concurrent_zoom_grids_with_the_same_fft_extent_keep_separate_kernels() {
    let samples: Vec<_> = (0..11)
        .map(|n| Complex64::new(f64::from(n).sin(), f64::from(n).cos()))
        .collect();
    std::thread::scope(|scope| {
        for (len, start, step, bins) in [
            (11, 0.07, 0.003, 5),
            (11, 0.11, 0.003, 5),
            (10, 0.07, 0.003, 5),
            (11, 0.07, -0.003, 5),
            (11, 0.07, 0.003, 6),
        ] {
            let input = &samples[..len];
            scope.spawn(move || {
                let result = signal::zoom_fft(input, start, step, bins).unwrap();
                assert_eq!(result.len(), bins);
                for (bin, actual) in result.into_iter().enumerate() {
                    let cycles = start + bin as f64 * step;
                    let expected = input
                        .iter()
                        .enumerate()
                        .map(|(i, sample)| {
                            sample * Complex64::from_polar(1.0, -TAU * cycles * i as f64)
                        })
                        .sum();
                    close_complex(actual, expected);
                }
            });
        }
    });
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
fn welch_retains_combined_subnormal_power() {
    let tiny = f64::from_bits(1);
    let sample = 1.28 * tiny.sqrt();
    let spectrum = signal::welch_psd(&[sample, sample, 0.0, 0.0], &[1.0; 4], 4, 2.0).unwrap();
    assert_eq!(spectrum.density[1], tiny);
}

#[test]
fn welch_parallel_blocks_match_weighted_frames_and_each_other() {
    for frame_len in [128, 129] {
        let count = 2051;
        let bin = 7;
        let window = vec![1.0; frame_len];
        let rate = frame_len as f64;
        let basis: Vec<_> = (0..frame_len)
            .map(|index| (TAU * bin as f64 * index as f64 / rate).sin_cos())
            .collect();
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;
        for frame in 0..count {
            let ax = f64::from(frame % 7 + 1) / 7.0;
            let ay = if frame % 5 < 3 { 0.75 } else { -0.25 };
            x.extend(basis.iter().map(|&(_, cosine)| ax * cosine));
            y.extend(basis.iter().map(|&(sine, _)| ay * sine));
            sum_xx += ax * ax;
            sum_yy += ay * ay;
            sum_xy += ax * ay;
        }
        // Trailing samples are outside every complete frame.
        x.push(f64::NAN);
        y.push(f64::NAN);
        let compute = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    (
                        signal::welch_psd(&x, &window, frame_len, rate).unwrap(),
                        signal::welch_cross_spectrum(&x, &y, &window, frame_len, rate).unwrap(),
                    )
                })
        };
        let (power, cross) = compute(1);
        assert_eq!((power.clone(), cross.clone()), compute(4));
        close(power.density[bin], sum_xx / (2.0 * f64::from(count)));
        close_complex(
            cross.density[bin],
            Complex64::new(0.0, -sum_xy / (2.0 * f64::from(count))),
        );
        close(cross.coherence[bin], sum_xy * sum_xy / (sum_xx * sum_yy));
        close(power.density.iter().sum(), power.density[bin]);
        x[(count as usize - 1) * frame_len] = f64::NAN;
        assert_eq!(
            signal::welch_psd(&x, &window, frame_len, rate),
            Err(SignalError::NonFiniteInput)
        );
    }
}

#[test]
fn time_series_lags_match_direct_moments_and_yule_walker_solutions() {
    let samples = [2.0, -1.0, 4.0, 0.0, 3.0, -2.0];
    let expected_covariance = [
        14.0 / 3.0,
        -19.0 / 6.0,
        7.0 / 3.0,
        -7.0 / 3.0,
        4.0 / 3.0,
        -0.5,
    ];
    let covariances = signal::autocovariance(&samples, samples.len() - 1).unwrap();
    let autocorrelations = signal::autocorrelation(&samples, samples.len() - 1).unwrap();
    for ((covariance, correlation), expected) in covariances
        .into_iter()
        .zip(autocorrelations)
        .zip(expected_covariance)
    {
        close(covariance, expected);
        close(correlation, expected / expected_covariance[0]);
    }

    // Rational last coefficients of the independent Toeplitz linear systems.
    let expected_partial = [
        1.0,
        -19.0 / 28.0,
        31.0 / 423.0,
        -1609.0 / 6356.0,
        -1759.0 / 5959.0,
        2056.0 / 19125.0,
    ];
    let partial = signal::partial_autocorrelation(&samples, samples.len() - 1).unwrap();
    for (actual, expected) in partial.into_iter().zip(expected_partial) {
        close(actual, expected);
    }
    assert_eq!(signal::partial_autocorrelation(&samples, 0).unwrap(), [1.0]);

    // Exercise the FFT crossover and its prefix zero-padding bound.
    let samples: Vec<_> = (0..1800)
        .map(|index| (f64::from(index) * 0.31).sin() + f64::from(index % 3))
        .collect();
    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    for (lag, covariance) in signal::autocovariance(&samples, 200)
        .unwrap()
        .into_iter()
        .enumerate()
    {
        let direct = samples
            .iter()
            .zip(&samples[lag..])
            .map(|(a, b)| (a - mean) * (b - mean))
            .sum::<f64>()
            / samples.len() as f64;
        close(covariance, direct);
    }
}

#[test]
fn time_series_lags_validate_domain_and_scale_before_squaring() {
    assert_eq!(signal::autocovariance(&[f64::MAX; 3], 2).unwrap(), [0.0; 3]);
    assert_eq!(signal::autocovariance(&[5.0], 0).unwrap(), [0.0]);
    assert_eq!(
        signal::autocorrelation(&[5.0; 3], 2),
        Err(SignalError::DegenerateVariance)
    );
    assert_eq!(
        signal::partial_autocorrelation(&[0.0; 2], 1),
        Err(SignalError::DegenerateVariance)
    );
    assert_eq!(
        signal::autocovariance(&[], 0),
        Err(SignalError::InsufficientSamples)
    );
    assert_eq!(
        signal::autocorrelation(&[1.0], usize::MAX),
        Err(SignalError::InvalidLag)
    );
    assert_eq!(
        signal::partial_autocorrelation(&[1.0, f64::INFINITY], 1),
        Err(SignalError::NonFiniteInput)
    );
    assert_eq!(
        signal::autocovariance(&[f64::MAX, -f64::MAX], 1),
        Err(SignalError::NonFiniteResult)
    );
    for amplitude in [f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1)] {
        let samples = [amplitude, -amplitude];
        assert_eq!(signal::autocorrelation(&samples, 1).unwrap(), [1.0, -0.5]);
        assert_eq!(
            signal::partial_autocorrelation(&samples, 1).unwrap(),
            [1.0, -0.5]
        );
    }

    let offset = 2.0_f64.powi(450);
    let unit = 2.0_f64.powi(398);
    let centered = [0.0, 1.0, 3.0, -1.0];
    let shifted = centered.map(|sample| offset + unit * sample);
    let want = signal::autocovariance(&centered, 2).unwrap();
    for (covariance, expected) in signal::autocovariance(&shifted, 2)
        .unwrap()
        .into_iter()
        .zip(want)
    {
        close(covariance / unit / unit, expected);
    }
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
        wavelets::inverse_wavelet_packet(&[], Wavelet::Haar),
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

#[test]
fn packed_fft_handles_degenerate_inputs_and_real_endpoints() {
    assert!(fft::fft(&[]).is_empty());
    assert!(fft::rfft(&[]).is_empty());
    assert!(fft::irfft(&[], 0).is_empty());
    assert_eq!(
        fft::fft(&[Complex64::new(3.0, -0.0)])[0].im.to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(fft::rfft(&[3.0]), [Complex64::new(3.0, 0.0)]);
    assert_eq!(fft::irfft(&[Complex64::new(3.0, f64::NAN)], 1), [3.0]);

    let original = [3.0, -2.0, 1.5, 2.0, 4.0, 0.5];
    let mut spectrum = fft::rfft(&original);
    assert_eq!(spectrum[0].im.to_bits(), 0);
    assert_eq!(spectrum[3].im.to_bits(), 0);
    spectrum[0].im = f64::NAN;
    spectrum[3].im = f64::INFINITY;
    for (got, want) in fft::irfft(&spectrum, original.len())
        .into_iter()
        .zip(original)
    {
        close(got, want);
    }
    let filled = fft::irfft(&[Complex64::ONE], 4);
    assert_eq!(filled, [0.25; 4]);
}

#[test]
fn wide_strided_fft_matches_independent_lanes() {
    let shape = [257, 1_024];
    let [rows, cols] = shape;
    let data = (0..rows * cols)
        .map(|index| {
            let index = index as f64;
            Complex64::new((index * 0.037).cos(), (index * 0.107).sin())
        })
        .collect::<Vec<_>>();
    let compute = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| fft::fft_axes(&data, &shape, &[0], false, fft::Norm::Backward).unwrap())
    };
    let parallel = compute(4);
    assert_eq!(compute(1), parallel);
    for column in [0, 1, 512, 1_023] {
        let lane = data
            .chunks_exact(cols)
            .map(|row| row[column])
            .collect::<Vec<_>>();
        for (row, expected) in parallel.chunks_exact(cols).zip(fft::fft(&lane)) {
            close_complex(row[column], expected);
        }
    }
}

#[test]
fn overlap_save_matches_direct_convolution_across_blocks_and_tail() {
    let taps = (0..37)
        .map(|index| (f64::from(index) * 0.13).sin() * 0.125)
        .collect::<Vec<_>>();
    let plan = OverlapSavePlan::new(&taps).unwrap();
    for len in [1, 8, 214, 215, 216, 537, 70_001] {
        let x = (0..len)
            .map(|index| (f64::from(index) * 0.11).cos())
            .collect::<Vec<_>>();
        let expected = signal::convolve(&x, &taps);
        let got = plan.apply(&x).unwrap();
        let automatic = fftconv::convolve(&taps, &x);
        assert_eq!(got.len(), expected.len());
        assert_eq!(automatic.len(), expected.len());
        for (got, expected) in got.into_iter().zip(&expected) {
            close(got, *expected);
        }
        for (got, expected) in automatic.into_iter().zip(expected) {
            close(got, expected);
        }
    }
    assert!(plan.apply(&[]).unwrap().is_empty());
    assert_eq!(
        plan.apply(&[f64::NAN]),
        Err(FilterError::Signal(SignalError::NonFiniteInput))
    );
    assert!(matches!(
        OverlapSavePlan::new(&[]),
        Err(FilterError::InvalidParameters)
    ));
}

#[test]
fn fft_convolution_retains_inputs_of_very_different_scales() {
    let x = [1e80; 30];
    let y = [1e-80; 31];
    for (got, expected) in fftconv::convolve(&x, &y)
        .into_iter()
        .zip(signal::convolve(&x, &y))
    {
        close(got, expected);
    }
    for period in [1, 2, 3, 4, 31, 32] {
        let a = (0..period)
            .map(|index| (index as f64 * 0.271).sin())
            .collect::<Vec<_>>();
        let b = (0..period)
            .map(|index| (index as f64 * 0.089).cos())
            .collect::<Vec<_>>();
        let output = fftconv::circular_convolve(&a, &b);
        for (lag, value) in output.into_iter().enumerate() {
            let expected = a
                .iter()
                .enumerate()
                .map(|(index, sample)| sample * b[(lag + period - index) % period])
                .sum();
            close(value, expected);
        }
    }
}

#[test]
fn fft_autocorrelation_is_linear_and_uses_nonnegative_lags() {
    for n in [0, 1, 2, 19, 48] {
        let samples = (0..n)
            .map(|index| (f64::from(index) * 0.4).cos())
            .collect::<Vec<_>>();
        let output = fftconv::autocorrelate(&samples);
        assert_eq!(output.len(), samples.len());
        for (lag, got) in output.into_iter().enumerate() {
            let expected = samples
                .iter()
                .zip(&samples[lag..])
                .map(|(a, b)| a * b)
                .sum();
            close(got, expected);
        }
    }
}

#[test]
fn windowed_fir_filters_have_symmetric_taps_and_expected_passbands() {
    let lowpass = signal::firwin_lowpass(65, 0.25, 6.0).unwrap();
    let bandpass = signal::firwin_bandpass(64, 0.3, 0.7, 6.0).unwrap();
    for taps in [&lowpass, &bandpass] {
        for (&left, &right) in taps.iter().zip(taps.iter().rev()) {
            close(left, right);
        }
    }
    let low = signal::frequency_response(&lowpass, &[1.0], &[0.0, 0.25, 1.0]).unwrap();
    close(low[0].norm(), 1.0);
    assert!(low[1].norm() < 0.001);
    close_complex(low[0], low[2]);
    let band = signal::frequency_response(&bandpass, &[1.0], &[0.0, 0.25, 0.5]).unwrap();
    assert!(band[0].norm() < 0.001);
    close(band[1].norm(), 1.0);
    assert!(band[2].norm() < 0.001);
    assert_eq!(
        signal::firwin_lowpass(0, 0.5, 2.0),
        Err(FilterError::InvalidParameters)
    );
    assert_eq!(
        signal::firwin_bandpass(10, 0.6, 0.5, 2.0),
        Err(FilterError::InvalidParameters)
    );
    assert_eq!(
        signal::firwin_lowpass(usize::MAX, 0.5, 2.0),
        Err(FilterError::Signal(SignalError::OutputSize))
    );
}

#[test]
fn second_order_sections_agree_with_direct_filtering_and_response() {
    let sections = [
        [2.0, 1.0, -0.5, 2.0, -0.25, 0.125],
        [1.0, 0.0, 0.5, 1.0, 0.2, 0.0],
    ];
    let x = (0..150)
        .map(|i| (f64::from(i) * 0.31).sin())
        .collect::<Vec<_>>();
    let expected = sections.iter().fold(x.clone(), |state, row| {
        signal::lfilter(&row[..3], &row[3..], &state)
    });
    let got = signal::sos_filter(&x, &sections).unwrap();
    for (got, expected) in got.into_iter().zip(expected) {
        close(got, expected);
    }

    assert_eq!(signal::sos_filter(&x, &[]).unwrap(), x);
    assert_eq!(
        signal::sos_filter(&x, &[[0.0; 6]]),
        Err(FilterError::InvalidParameters)
    );
    let response = signal::frequency_response(&[1.0], &[1.0, -0.5], &[0.0, 0.25, 0.5]).unwrap();
    close_complex(response[0], Complex64::new(2.0, 0.0));
    close_complex(response[1], Complex64::new(0.8, -0.4));
    close_complex(response[2], Complex64::new(2.0 / 3.0, 0.0));
    assert_eq!(
        signal::frequency_response(&[1.0], &[1.0, -1.0], &[0.0]),
        Err(FilterError::Signal(SignalError::NonFiniteResult))
    );
}

#[test]
fn polyphase_resampling_is_centered_and_rejects_aliases() {
    let impulse = (0..51)
        .map(|index| f64::from(index == 25))
        .collect::<Vec<_>>();
    let output = signal::resample_poly(&impulse, 4, 2).unwrap();
    assert_eq!(output, signal::resample_poly(&impulse, 2, 1).unwrap());
    assert_eq!(output.len(), 102);
    assert!((output[50] - 1.0).abs() < 0.002);
    for offset in 1..=30 {
        close(output[50 + offset], output[50 - offset]);
    }
    let signal = (0..600)
        .map(|i| {
            let i = f64::from(i);
            (TAU * 0.05 * i).cos() + (TAU * 0.35 * i).sin()
        })
        .collect::<Vec<_>>();
    let down = signal::resample_poly(&signal, 1, 3).unwrap();
    assert_eq!(down.len(), 200);
    for (index, value) in down.into_iter().enumerate().skip(15).take(170) {
        let expected = (TAU * 0.15 * index as f64).cos();
        assert!((value - expected).abs() < 0.005, "{value} vs {expected}");
    }
    assert_eq!(
        signal::resample_poly(&[1.0, 2.0], 3, 3).unwrap(),
        [1.0, 2.0]
    );
    assert_eq!(
        signal::resample_poly(&[1.0], usize::MAX, 1),
        Err(FilterError::Signal(SignalError::OutputSize))
    );
    assert_eq!(
        signal::resample_poly(&[1.0], 1, 0),
        Err(FilterError::Signal(SignalError::InvalidRate))
    );
}

#[test]
fn stationary_wavelets_reconstruct_and_commute_with_periodic_shifts() {
    for wavelet in Wavelet::ALL {
        let original = (0..29)
            .map(|index| (f64::from(index) * 0.26).cos())
            .collect::<Vec<_>>();
        let bands = wavelets::swt(&original, wavelet, 4).unwrap();
        let back = wavelets::iswt(&bands, wavelet).unwrap();
        let mut rotated = original.clone();
        rotated.rotate_right(3);
        let shifted = wavelets::swt(&rotated, wavelet, 4).unwrap();
        assert_eq!(bands.len(), 5);
        for (band, shifted) in bands.into_iter().zip(shifted) {
            assert_eq!(band.len(), original.len());
            for (value, shifted) in band.iter().cycle().skip(26).take(29).zip(shifted) {
                close(*value, shifted);
            }
        }
        for (got, expected) in back.into_iter().zip(original) {
            close(got, expected);
        }
    }
    assert_eq!(
        wavelets::swt(&[1.0], Wavelet::Haar, 1),
        Err(WaveletProcessingError::Transform(
            WaveletError::TooManyLevels
        ))
    );
    assert_eq!(
        wavelets::iswt(&[vec![1.0; 4], vec![0.0; 3]], Wavelet::Haar),
        Err(WaveletProcessingError::Transform(WaveletError::Shape))
    );
}

#[test]
fn stationary_wavelet_parallel_buffers_preserve_order_and_overflow_checks() {
    let signal: Vec<_> = (0..4099)
        .map(|i| (f64::from(i) * 0.079).sin() - (f64::from(i) * 0.173).cos())
        .collect();
    let run = |workers| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(|| {
                let bands = wavelets::swt(&signal, Wavelet::Db6, 5).unwrap();
                let restored = wavelets::iswt(&bands, Wavelet::Db6).unwrap();
                let overflow = vec![f64::MAX; signal.len()];
                assert_eq!(
                    wavelets::swt(&overflow, Wavelet::Haar, 1),
                    Err(WaveletProcessingError::NonFiniteResult)
                );
                let bits: Vec<_> = bands.into_iter().flatten().map(f64::to_bits).collect();
                (bits, restored)
            })
    };
    let (serial, signal_serial) = run(1);
    let (parallel, signal_parallel) = run(4);
    assert_eq!(serial, parallel);
    assert_eq!(serial.len(), 6 * signal.len());
    for ((a, b), expected) in signal_serial.into_iter().zip(signal_parallel).zip(signal) {
        assert_eq!(a.to_bits(), b.to_bits());
        close(a, expected);
    }
}

#[test]
fn wavelet_denoising_treats_approximation_and_details_separately() {
    use wavelets::Boundary;

    let signal = [3.0, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0];
    let hard = |threshold| {
        wavelets::wavelet_hard_threshold(&signal, threshold, Wavelet::Haar, 1, Boundary::Periodic)
    };
    let soft = |threshold| {
        wavelets::wavelet_soft_threshold(&signal, threshold, Wavelet::Haar, 1, Boundary::Periodic)
    };
    for (got, expected) in hard(0.5).unwrap().into_iter().zip(signal) {
        close(got, expected);
    }
    for got in hard(2.0).unwrap() {
        close(got, 2.0);
    }
    for (got, expected) in soft(0.5).unwrap().into_iter().zip(signal) {
        let shrink = std::f64::consts::FRAC_1_SQRT_2 * 0.5;
        close(got, expected - (expected - 2.0).signum() * shrink);
    }
    assert_eq!(soft(-0.1), Err(WaveletProcessingError::InvalidCoefficient));
    assert_eq!(Wavelet::Db10.coeff_len(0, Boundary::Zero), 0);
    assert_eq!(
        wavelets::idwt2(&[], usize::MAX, 1, Wavelet::Haar, Boundary::Periodic),
        Err(WaveletError::Shape)
    );
}

#[test]
fn frequency_response_keeps_exact_sweep_order_across_worker_counts() {
    let numerator: Vec<_> = (0..4096)
        .map(|index| f64::from(index % 11 - 5) / 4096.0)
        .collect();
    let frequencies: Vec<_> = (0..64).map(|index| f64::from(index) / 64.0).collect();
    let compute = |workers| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(|| {
                signal::frequency_response(&numerator, &[1.0, -0.25], &frequencies).unwrap()
            })
    };
    assert_eq!(compute(1), compute(4));
}

#[test]
fn savgol_reproduces_quadratics_in_the_interior_and_constants_at_edges() {
    let expected = [-3.0, 12.0, 17.0, 12.0, -3.0];
    let coefficients = signal::savgol_coeffs(5, 2);
    assert_eq!(coefficients.len(), expected.len());
    for (actual, expected) in coefficients.into_iter().zip(expected) {
        assert!((actual - expected / 35.0).abs() < 1e-12);
    }
    let pulse: Vec<_> = (0_i32..13).map(|x| f64::from(x * x - 3 * x + 2)).collect();
    let smooth = signal::savgol_filter(&pulse, 5, 2);
    assert_eq!(smooth.len(), pulse.len());
    for (&actual, &expected) in smooth[2..11].iter().zip(&pulse[2..11]) {
        assert!((actual - expected).abs() < 1e-10);
    }
    let flat = signal::savgol_filter(&[7.0; 8], 5, 2);
    assert_eq!(flat.len(), 8);
    assert!(flat.iter().all(|&value| (value - 7.0).abs() < 1e-12));
    assert!(signal::savgol_coeffs(4, 2).is_empty());
    assert!(signal::savgol_coeffs(5, 5).is_empty());
    assert!(signal::savgol_filter(&[], 5, 2).is_empty());
}

#[test]
fn inverse_stft_retains_representable_subnormal_inverse() {
    let tiny = f64::from_bits(1);
    let frame = [vec![Complex64::new(tiny, 0.0); 2]];
    let output = signal::inverse_stft(&frame, &[1.0, 1.0], 2).unwrap();
    assert_eq!(output, [tiny, 0.0]);
}

#[test]
fn cross_welch_retains_subnormal_self_density_and_coherence() {
    let tiny = f64::from_bits(1);
    let sample = 1.28 * tiny.sqrt();
    let inputs = [sample, sample, 0.0, 0.0];
    let spectrum = signal::welch_cross_spectrum(&inputs, &inputs, &[1.0; 4], 4, 2.0).unwrap();
    assert_eq!(spectrum.density[1], Complex64::new(tiny, 0.0));
    assert!((spectrum.coherence[1] - 1.0).abs() < 1e-15);
}

#[test]
fn cross_welch_normalizes_independent_signal_amplitudes() {
    let x = [1e-200, 0.0, -1e-200, 0.0];
    let y = [0.0, 1e200, 0.0, -1e200];
    let spectrum = signal::welch_cross_spectrum(&x, &y, &[1.0; 4], 4, 4.0).unwrap();
    assert!((spectrum.density[1] - Complex64::new(0.0, -0.5)).norm() < 1e-15);
    assert!((spectrum.coherence[1] - 1.0).abs() < 1e-15);
}

#[test]
fn cross_welch_uses_only_covered_samples_when_finding_scales() {
    let x = [1.0, f64::NAN, 2.0, f64::INFINITY];
    let result = signal::welch_cross_spectrum(&x, &x, &[1.0], 2, 1.0).unwrap();
    assert_eq!(result.density, [Complex64::new(2.5, 0.0)]);
    assert!((result.coherence[0] - 1.0).abs() < 1e-15);
}

#[test]
fn cross_welch_keeps_coherence_when_density_is_too_small_to_represent() {
    let x = [1e-200, 0.0, -1e-200, 0.0];
    let result = signal::welch_cross_spectrum(&x, &x, &[1.0; 4], 4, 4.0).unwrap();
    assert_eq!(result.density[1], Complex64::ZERO);
    assert!((result.coherence[1] - 1.0).abs() < 1e-15);
}

#[test]
fn cross_welch_refuses_true_density_overflow() {
    let x = [1e200, 0.0, -1e200, 0.0];
    assert_eq!(
        signal::welch_cross_spectrum(&x, &x, &[1.0; 4], 4, 4.0),
        Err(signal::SignalError::NonFiniteResult)
    );
}

#[test]
fn cross_welch_recovers_when_window_products_round_to_zero() {
    let tiny = f64::from_bits(1);
    let x = [tiny, tiny, 0.0, 0.0];
    let y = [1e150, 1e150, 0.0, 0.0];
    let spectrum = signal::welch_cross_spectrum(&x, &y, &[1.0; 4], 4, 2.0).unwrap();
    let want = tiny * 5e149;
    assert!((spectrum.density[0].re / want - 1.0).abs() < 1e-14);
    assert!((spectrum.coherence[0] - 1.0).abs() < 1e-15);
}

#[test]
fn cross_welch_combines_subnormal_samples_and_a_subnormal_rate() {
    let tiny = f64::from_bits(1);
    let x = [tiny; 4];
    let spectrum = signal::welch_cross_spectrum(&x, &x, &[1.0; 4], 4, tiny).unwrap();
    assert_eq!(spectrum.density[0], Complex64::new(f64::from_bits(4), 0.0));
    assert!((spectrum.coherence[0] - 1.0).abs() < 1e-15);
}
