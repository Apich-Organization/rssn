//! Independent references and boundary checks for framed signals and digraphs.

use std::f64::consts::TAU;

use rssn::kernels::graph::Graph;
use rssn::kernels::graph::GraphStructureError;
use rssn::kernels::graph::strongly_connected_components;
use rssn::kernels::graph::topological_sort;
use rssn::kernels::signal::SignalError;
use rssn::kernels::signal::butterworth_lowpass;
use rssn::kernels::signal::dct;
use rssn::kernels::signal::dst;
use rssn::kernels::signal::fft;
use rssn::kernels::signal::hamming_window;
use rssn::kernels::signal::hann_window;
use rssn::kernels::signal::idct;
use rssn::kernels::signal::idst;
use rssn::kernels::signal::inverse_stft;
use rssn::kernels::signal::kaiser_window;
use rssn::kernels::signal::lfilter;
use rssn::kernels::signal::stft;
use rssn::kernels::signal::upfirdn;
use rssn::kernels::signal::welch_psd;
use rssn::kernels::special::bessel_i0;
use rssn::kernels::special::bessel_i0e;
use rustfft::num_complex::Complex;

fn close(
    got: f64,
    expected: f64,
) {
    let scale = 1.0_f64.max(got.abs()).max(expected.abs());
    assert!(
        (got - expected).abs() <= 2e-12 * scale,
        "got {got}, expected {expected}"
    );
}

fn check_values(
    got: &[f64],
    expected: &[f64],
) {
    assert_eq!(got.len(), expected.len());
    for (&got, &expected) in got.iter().zip(expected) {
        close(got, expected);
    }
}

fn direct_spectrum(
    samples: &[f64],
    window: &[f64],
    bin: usize,
) -> Complex<f64> {
    samples
        .iter()
        .zip(window)
        .enumerate()
        .map(|(time, (sample, weight))| {
            let phase = -TAU * (bin * time) as f64 / window.len() as f64;
            Complex::from_polar(sample * weight, phase)
        })
        .sum()
}

#[test]
fn welch_matches_averaged_direct_periodograms_for_odd_and_even_lengths() {
    let data = [1.5, -2.0, 3.0, 9.0, 7.0, -4.0, 2.0, -1.0, 0.5, 7.0, 99.0];
    let rate = 24.0;
    let hop = 3;
    for width in [5, 6] {
        let window = hamming_window(width);
        let frames: Vec<_> = data.windows(width).step_by(hop).collect();
        let normalization = rate * window.iter().map(|weight| weight * weight).sum::<f64>();
        let result = welch_psd(&data, &window, hop, rate).unwrap();
        assert_eq!(result.frequencies.len(), width / 2 + 1);
        for (bin, (&freq, &density)) in result.frequencies.iter().zip(&result.density).enumerate() {
            close(freq, bin as f64 * rate / width as f64);
            let phase_partners = if bin == 0 || bin * 2 == width {
                1.0
            } else {
                2.0
            };
            let sum: f64 = frames
                .iter()
                .map(|frame| direct_spectrum(frame, &window, bin).norm_sqr())
                .sum();
            close(
                density,
                phase_partners * sum / (normalization * frames.len() as f64),
            );
        }
    }
}

#[test]
fn welch_integrates_to_mean_weighted_power() {
    let data = [-1.0, 5.0, 3.0, 9.0, -1.0, -2.0, 2.0, 8.0];
    let window = [0.0, 1.0, 3.0, 0.25, -2.0];
    let hop = 2;
    let rate = 9.0;
    let result = welch_psd(&data, &window, hop, rate).unwrap();
    let norm: f64 = window.iter().map(|w| w * w).sum();
    let frames: Vec<_> = data.windows(window.len()).step_by(hop).collect();
    let total: f64 = frames
        .iter()
        .flat_map(|frame| frame.iter().zip(window).map(|(x, w)| (x * w).powi(2)))
        .sum();
    close(
        result.density.iter().sum::<f64>() * rate / window.len() as f64,
        total / (norm * frames.len() as f64),
    );
}

#[test]
fn welch_distinguishes_dc_nyquist_and_odd_last_bins() {
    let dc = welch_psd(&[3.0; 4], &[1.0; 4], 4, 4.0).unwrap();
    check_values(&dc.density, &[9.0, 0.0, 0.0]);
    let nyquist = welch_psd(&[3.0, -3.0, 3.0, -3.0], &[1.0; 4], 4, 4.0).unwrap();
    check_values(&nyquist.density, &[0.0, 0.0, 9.0]);
    let odd: Vec<_> = (0..5)
        .map(|i| (2.0 * TAU * f64::from(i) / 5.0).cos())
        .collect();
    let odd = welch_psd(&odd, &[1.0; 5], 5, 5.0).unwrap();
    check_values(&odd.density, &[0.0, 0.0, 0.5]);
    let singleton = welch_psd(&[1.0, 2.0], &[9.0], 1, 0.5).unwrap();
    check_values(&singleton.frequencies, &[0.0]);
    check_values(&singleton.density, &[5.0]);
}

#[test]
fn welch_is_invariant_to_scaling_and_ignores_an_incomplete_tail() {
    let data = [5.0, -1.0, 4.0, 1.0, -4.0, f64::NAN];
    let regular = welch_psd(&data, &[1.0, -0.5, 0.25], 2, 3.0).unwrap();
    for scale in [1e300, 1e-300] {
        let window = [scale, -0.5 * scale, 0.25 * scale];
        let rescaled = welch_psd(&data, &window, 2, 3.0).unwrap();
        check_values(&regular.density, &rescaled.density);
    }
}

#[test]
fn welch_rejects_bad_domains_and_unrepresentable_power() {
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            welch_psd(&[1.0], &[1.0], 1, invalid),
            Err(SignalError::InvalidSampleRate)
        );
    }
    for window in [&[][..], &[0.0, -0.0], &[1.0, f64::NAN]] {
        assert_eq!(
            welch_psd(&[1.0, 1.0], window, 1, 1.0),
            Err(SignalError::InvalidWindow)
        );
    }
    assert_eq!(
        welch_psd(&[1.0], &[1.0], 0, 1.0),
        Err(SignalError::InvalidHop)
    );
    assert_eq!(
        welch_psd(&[], &[1.0], 1, 1.0),
        Err(SignalError::InsufficientSamples)
    );
    assert_eq!(
        welch_psd(&[f64::INFINITY], &[1.0], 1, 1.0),
        Err(SignalError::NonFiniteInput)
    );
    assert_eq!(
        welch_psd(&[f64::MAX], &[1.0], 1, 1.0),
        Err(SignalError::NonFiniteResult)
    );
}

#[test]
fn stft_matches_direct_dft_and_preserves_full_frame_order() {
    let data = [6.0, 0.0, -1.0, 10.0, -4.0, 2.5, 7.0, -91.0];
    for window in [&[1.0, 0.2, -0.75][..], &[0.5, -1.0, 3.0, 2.0]] {
        let frames = stft(&data, window, 3);
        let inputs: Vec<_> = data.windows(window.len()).step_by(3).collect();
        assert_eq!(frames.len(), inputs.len());
        for (input, spectrum) in inputs.into_iter().zip(frames) {
            assert_eq!(spectrum.len(), window.len());
            for (bin, value) in spectrum.into_iter().enumerate() {
                let reference = direct_spectrum(input, window, bin);
                close(value.re, reference.re);
                close(value.im, reference.im);
            }
        }
    }
}

#[test]
fn empty_and_incompatible_stft_inputs_do_not_enter_the_fft_planner() {
    assert!(fft(&mut []).is_empty());
    assert!(stft(&[], &[1.0], 1).is_empty());
    assert!(stft(&[1.0], &[], 1).is_empty());
    assert!(stft(&[1.0], &[1.0], 0).is_empty());
    assert!(stft(&[1.0], &[1.0, 2.0], 1).is_empty());
}

#[test]
fn inverse_stft_recovers_signed_and_overlapping_windows() {
    let input: Vec<_> = (0..27)
        .map(|i| (0.2 * f64::from(i)).cos() + 0.3 * f64::from(i))
        .collect();
    for window in [
        vec![1.5, -1.0, 0.0, 2.5, 1.0],
        hamming_window(8),
        hamming_window(7),
    ] {
        let hop = 2;
        let frames = stft(&input, &window, hop);
        let recovered = inverse_stft(&frames, &window, hop).unwrap();
        let covered = (frames.len() - 1) * hop + window.len();
        check_values(&recovered, &input[..covered]);
    }
}

#[test]
fn inverse_stft_retains_nonzero_weights_with_underflowing_squares() {
    let window = [1.0, 1e-170];
    let input = [0.0, 1.0];
    let frames = stft(&input, &window, 2);
    check_values(&inverse_stft(&frames, &window, 2).unwrap(), &input);
}

#[test]
fn inverse_stft_reports_uncovered_endpoints_and_gaps() {
    let window = hann_window(5);
    let frames = stft(&[1.0; 7], &window, 2);
    assert_eq!(
        inverse_stft(&frames, &window, 2),
        Err(SignalError::UncoveredSample(0))
    );
    let window = [1.0; 2];
    let frames = stft(&[1.0; 8], &window, 3);
    assert_eq!(
        inverse_stft(&frames, &window, 3),
        Err(SignalError::UncoveredSample(2))
    );
}

#[test]
fn inverse_stft_validates_shape_rate_and_finiteness() {
    let bins = [vec![Complex::ZERO]];
    assert_eq!(inverse_stft(&bins, &[1.0], 0), Err(SignalError::InvalidHop));
    assert_eq!(inverse_stft(&bins, &[], 1), Err(SignalError::InvalidWindow));
    assert_eq!(
        inverse_stft(&[], &[1.0], 1),
        Err(SignalError::InsufficientSamples)
    );
    assert_eq!(
        inverse_stft(&bins, &[1.0; 2], 1),
        Err(SignalError::FrameShape)
    );
    let bad = [vec![Complex::new(0.0, f64::NAN)]];
    assert_eq!(
        inverse_stft(&bad, &[1.0], 1),
        Err(SignalError::NonFiniteInput)
    );
    let oversized = [bins[0].clone(), bins[0].clone()];
    assert_eq!(
        inverse_stft(&oversized, &[1.0], usize::MAX),
        Err(SignalError::OutputSize)
    );
}

fn direct_resample(
    input: &[f64],
    taps: &[f64],
    up: usize,
    down: usize,
) -> Vec<f64> {
    let mut expanded = vec![0.0; (input.len() - 1) * up + 1];
    for (index, &x) in input.iter().enumerate() {
        expanded[index * up] = x;
    }
    let mut linear = vec![0.0; expanded.len() + taps.len() - 1];
    for (i, x) in expanded.into_iter().enumerate() {
        for (j, &h) in taps.iter().enumerate() {
            linear[i + j] += x * h;
        }
    }
    linear.into_iter().step_by(down).collect()
}

#[test]
fn polyphase_matches_explicit_zero_insertion_and_convolution() {
    let data = [-1.5, 3.0, 1.25, -4.5];
    let taps = [1.0, 3.0, 1.0, -3.0, 0.0, -5.0, 4.0];
    for len in 1..=data.len() {
        for width in 1..=taps.len() {
            for up in 1..=8 {
                for down in 1..=8 {
                    let expected = direct_resample(&data[..len], &taps[..width], up, down);
                    let output = upfirdn(&data[..len], &taps[..width], up, down).unwrap();
                    check_values(&output, &expected);
                }
            }
        }
    }
}

#[test]
fn polyphase_can_skip_an_unallocatable_intermediate() {
    let output = upfirdn(&[1.0, 2.0, 3.0], &[5.0, 6.0], usize::MAX, usize::MAX).unwrap();
    check_values(&output, &[5.0, 10.0, 15.0]);
    assert_eq!(
        upfirdn(&[1.0; 2], &[1.0], usize::MAX, 1),
        Err(SignalError::OutputSize)
    );
    assert_eq!(
        upfirdn(&[1.0; 2], &[1.0], isize::MAX as usize, 1),
        Err(SignalError::OutputSize)
    );
}

#[test]
fn polyphase_checks_domain_and_handles_empty_inputs() {
    assert!(upfirdn(&[], &[1.0], 1, 1).unwrap().is_empty());
    assert!(upfirdn(&[1.0], &[], 1, 1).unwrap().is_empty());
    assert_eq!(upfirdn(&[1.0], &[1.0], 0, 1), Err(SignalError::InvalidRate));
    assert_eq!(upfirdn(&[1.0], &[1.0], 1, 0), Err(SignalError::InvalidRate));
    assert_eq!(
        upfirdn(&[1.0], &[f64::NAN], 1, 1),
        Err(SignalError::NonFiniteInput)
    );
    assert_eq!(
        upfirdn(&[f64::MAX], &[2.0], 1, 1),
        Err(SignalError::NonFiniteResult)
    );
}

#[test]
fn topology_and_scc_agree_with_exhaustive_three_vertex_reachability() {
    const N: usize = 3;
    for mask in 0_u16..(1 << (N * N)) {
        let mut graph = Graph::new(N);
        let mut reachable = [[false; N]; N];
        for (from, row) in reachable.iter_mut().enumerate() {
            for (to, path) in row.iter_mut().enumerate() {
                if mask & (1 << (from * N + to)) != 0 {
                    *path = true;
                    graph.add_edge(from, to, f64::NAN);
                }
            }
        }
        for via in 0..N {
            for from in 0..N {
                for to in 0..N {
                    reachable[from][to] |= reachable[from][via] && reachable[via][to];
                }
            }
        }

        let labels = strongly_connected_components(&graph).unwrap();
        assert_eq!(labels.len(), N);
        for from in 0..N {
            for to in 0..N {
                let mutually_reachable = from == to || (reachable[from][to] && reachable[to][from]);
                assert_eq!(
                    labels[from] == labels[to],
                    mutually_reachable,
                    "mask {mask}"
                );
                if reachable[from][to] && !mutually_reachable {
                    assert!(labels[from] < labels[to], "mask {mask}");
                }
            }
        }
        if (0..N).any(|vertex| reachable[vertex][vertex]) {
            assert_eq!(
                topological_sort(&graph),
                Err(GraphStructureError::DirectedCycle)
            );
        } else {
            let order = topological_sort(&graph).unwrap();
            assert_eq!(order.len(), N);
            for (index, &vertex) in order.iter().enumerate() {
                assert!(!order[..index].contains(&vertex), "mask {mask}");
                assert!(
                    order[index + 1..]
                        .iter()
                        .all(|&after| !reachable[after][vertex]),
                    "mask {mask}"
                );
            }
        }
    }
}

#[test]
fn topology_counts_parallel_arcs_and_keeps_isolated_vertices() {
    let mut graph = Graph::new(5);
    graph.add_edge(0, 1, 1.0);
    graph.add_edge(0, 1, 2.0);
    graph.add_edge(1, 3, 1.0);
    graph.add_edge(2, 3, 1.0);
    assert_eq!(topological_sort(&graph).unwrap(), [0, 2, 4, 1, 3]);
    assert_eq!(topological_sort(&Graph::new(0)).unwrap(), []);
    assert_eq!(strongly_connected_components(&Graph::new(0)).unwrap(), []);
}

#[test]
fn graph_operations_reject_out_of_range_edge_destinations() {
    let mut graph = Graph::new(3);
    graph.add_edge(0, 1, 1.0);
    graph.add_edge(2, usize::MAX, 1.0);
    assert_eq!(
        topological_sort(&graph),
        Err(GraphStructureError::InvalidVertex)
    );
    assert_eq!(
        strongly_connected_components(&graph),
        Err(GraphStructureError::InvalidVertex)
    );
}

#[test]
fn deep_directed_path_and_cycle_use_no_call_stack_recursion() {
    let size = 65_536;
    let mut graph = Graph::new(size);
    for u in 1..size {
        graph.add_edge(u - 1, u, 0.0);
    }
    assert!(topological_sort(&graph).unwrap().into_iter().eq(0..size));
    assert!(
        strongly_connected_components(&graph)
            .unwrap()
            .into_iter()
            .eq(0..size)
    );

    graph.add_edge(size - 1, 0, 0.0);
    assert_eq!(
        topological_sort(&graph),
        Err(GraphStructureError::DirectedCycle)
    );
    assert!(
        strongly_connected_components(&graph)
            .unwrap()
            .into_iter()
            .all(|label| label == 0)
    );
}


#[test]
fn real_sine_and_cosine_transforms_are_orthonormal() {
    for width in [1, 2, 5, 12, 19] {
        let samples: Vec<_> = (0..width)
            .map(|index| (f64::from(index) * 0.7).cos() + f64::from(index * index))
            .collect();
        let sine = dst(&samples);
        let cosine = dct(&samples);
        check_values(&idst(&sine), &samples);
        check_values(&idct(&cosine), &samples);
        let power = samples.iter().map(|value| value * value).sum::<f64>();
        close(sine.iter().map(|value| value * value).sum(), power);
        close(cosine.iter().map(|value| value * value).sum(), power);
    }
    for transform in [dct, idct, dst, idst] {
        assert!(transform(&[]).is_empty());
    }
}

#[test]
fn filter_matches_the_analytic_first_order_step_response() {
    let step = lfilter(&[0.25], &[0.5, -0.25], &[1.0; 16]);
    for (n, response) in (1..).zip(step) {
        close(response, 1.0 - 0.5_f64.powi(n));
    }
    assert!(
        lfilter(&[1.0], &[0.0], &[1.0; 3])
            .iter()
            .all(|v| v.is_nan())
    );
    assert!(lfilter(&[], &[1.0], &[1.0]).is_empty());
}

#[test]
fn butterworth_preserves_dc_and_has_the_requested_half_power_cutoff() {
    let cutoff = 0.3;
    for order in 1..=6 {
        let (b, a) = butterworth_lowpass(order, cutoff);
        assert_eq!(b.len(), order + 1);
        assert_eq!(a.len(), order + 1);
        let response = |frequency: f64| {
            let evaluate = |coefficients: &[f64]| {
                coefficients
                    .iter()
                    .enumerate()
                    .map(|(index, value)| Complex::from_polar(*value, -frequency * index as f64))
                    .sum::<Complex<f64>>()
            };
            evaluate(&b) / evaluate(&a)
        };
        close(response(0.0).norm(), 1.0);
        close(response(cutoff * TAU / 2.0).norm_sqr(), 0.5);
        close(response(TAU / 2.0).norm(), 0.0);
    }
}

#[test]
fn scaled_bessel_i0_matches_high_precision_series_across_switch_and_overflow() {
    for (x, expected) in [
        (0.0, 1.0),
        (1.0, 0.465_759_607_593_640_43),
        (10.0, 0.127_833_337_163_428_6),
        (30.0, 0.073_145_946_482_237_3),
        (31.0, 0.071_946_496_696_983_83),
        (710.0, 0.014_974_675_005_024_156),
        (1_000.0, 0.012_617_240_455_891_257),
    ] {
        close(bessel_i0e(x) / expected, 1.0);
        close(bessel_i0e(-x) / expected, 1.0);
    }
    close(bessel_i0e(f64::MAX) * f64::MAX.sqrt(), 1.0 / TAU.sqrt());
    assert_eq!(bessel_i0e(f64::INFINITY), 0.0);
    assert!(bessel_i0e(f64::NAN).is_nan());
    assert!(bessel_i0(710.0).is_finite());
    assert_eq!(bessel_i0(f64::MAX), f64::INFINITY);
}

#[test]
fn kaiser_window_preserves_finite_weights_when_bessel_i0_overflows() {
    let shape = kaiser_window(7, 1_000.0);
    let expected = [
        0.0,
        2.974_057_683_214_161_4e-111,
        1.496_489_835_603_844_6e-25,
        1.0,
        1.496_489_835_603_844_6e-25,
        2.974_057_683_214_161_4e-111,
        0.0,
    ];
    for (value, known) in shape.iter().zip(expected) {
        if known == 0.0 {
            assert_eq!(*value, known);
        } else {
            close(value / known, 1.0);
        }
    }
    assert_eq!(shape, kaiser_window(7, -1_000.0));
    assert_eq!(
        kaiser_window(3, 748.0),
        vec![f64::from_bits(2), 1.0, f64::from_bits(2)]
    );
    assert_eq!(kaiser_window(5, f64::MAX), vec![0.0, 0.0, 1.0, 0.0, 0.0]);
    assert!(kaiser_window(2, f64::NAN).iter().all(|v| v.is_nan()));
}
