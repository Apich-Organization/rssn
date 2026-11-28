//! Seeded and statistical checks for repeated jump mark sampling.

use rssn::sim::stochastic::JumpDistribution;
use rssn::sim::stochastic::JumpSpec;
use rssn::sim::stochastic::Method;
use rssn::sim::stochastic::Options;
use rssn::sim::stochastic::SdeError;
use rssn::sim::stochastic::Stream;
use rssn::sim::stochastic::Structure;
use rssn::sim::stochastic::ensemble_seed;
use rssn::sim::stochastic::integrate_jump_diffusion;
use std::f64::consts::PI;
use std::f64::consts::TAU;

struct ReferenceStream(u64);

impl ReferenceStream {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        (self.0 >> 11) as f64 * 2.0_f64.powi(-53)
    }

    fn normal(&mut self) -> f64 {
        let radius = (-2.0 * self.uniform().max(f64::MIN_POSITIVE).ln()).sqrt();
        radius * (TAU * self.uniform()).cos()
    }

    fn poisson_small(
        &mut self,
        expectation: f64,
    ) -> usize {
        let limit = (-expectation).exp();
        let mut count = 0;
        let mut product = self.uniform();
        while product > limit {
            product *= self.uniform();
            count += 1;
        }
        count
    }
}

fn assert_reference(
    distribution: JumpDistribution,
    mut next_mark: impl FnMut(&mut ReferenceStream) -> [f64; 2],
) {
    let seed = 123_456;
    let mut stream = Stream::new(seed);
    let mut scalar_stream = stream.clone();
    let mut reference = ReferenceStream(seed);
    let mut output = [0.0; 74];
    distribution.fill_marks(&mut stream, &mut output).unwrap();
    for pair in output.as_chunks::<2>().0 {
        let expected = next_mark(&mut reference);
        assert_eq!(
            pair.as_slice(),
            distribution.sample(&mut scalar_stream).unwrap()
        );
        for (&a, e) in pair.iter().zip(expected) {
            assert!((a - e).abs() <= 8.0e-15 * e.abs().max(1.0), "{a} vs {e}");
        }
    }
    assert_eq!(stream, scalar_stream);
    assert_eq!(stream.next_unit(), reference.uniform());
}

#[test]
fn batch_sampling_preserves_each_distribution_stream_order() {
    assert_reference(JumpDistribution::Constant(vec![1.3, -2.1]), |_| [1.3, -2.1]);
    assert_reference(
        JumpDistribution::Normal {
            mean: vec![0.7, -1.2],
            std: vec![0.0, 2.5],
        },
        |stream| [0.7 + 0.0 * stream.normal(), -1.2 + 2.5 * stream.normal()],
    );
    assert_reference(
        JumpDistribution::Uniform {
            low: vec![-2.0, 1.2],
            high: vec![2.0, 1.2],
        },
        |stream| [-2.0 + 4.0 * stream.uniform(), 1.2 + 0.0 * stream.uniform()],
    );
    for p_up in [0.0, 0.3, 1.0] {
        assert_reference(
            JumpDistribution::DoubleExponential {
                p_up,
                eta_up: 2.0,
                eta_down: 4.0,
                dim: 2,
            },
            |stream| {
                std::array::from_fn(|_| {
                    let upward = stream.uniform() < p_up;
                    let magnitude = -stream.uniform().max(f64::MIN_POSITIVE).ln();
                    if upward {
                        magnitude / 2.0
                    } else {
                        -magnitude / 4.0
                    }
                })
            },
        );
    }
    assert_reference(
        JumpDistribution::AlphaStable {
            alpha: 2.0,
            beta: -0.8,
            scale: 0.4,
            location: 0.2,
            dim: 2,
        },
        |stream| std::array::from_fn(|_| 0.2 + 0.4 * std::f64::consts::SQRT_2 * stream.normal()),
    );
    assert_reference(
        JumpDistribution::AlphaStable {
            alpha: 1.0,
            beta: 0.0,
            scale: 0.6,
            location: -0.2,
            dim: 2,
        },
        |stream| {
            std::array::from_fn(|_| {
                let angle = PI * (stream.uniform() - 0.5);
                let _exponential = stream.uniform();
                0.6 * angle.tan() - 0.2
            })
        },
    );
}

#[test]
fn uniform_crossing_the_float_range_has_finite_marks_and_the_correct_mean() {
    let limit = f64::MAX;
    let wide = JumpDistribution::Uniform {
        low: vec![-limit, -limit],
        high: vec![limit, limit / 2.0],
    };
    assert_eq!(wide.mean_mark(), Some(vec![0.0, -limit / 4.0]));
    let mut sampled = [0.0; 512];
    let mut stream = Stream::new(0);
    wide.fill_marks(&mut stream, &mut sampled).unwrap();
    assert_eq!(sampled[0], -limit); // The first LCG draw for seed zero is exactly zero.
    let mut reference = ReferenceStream(0);
    for mark in sampled.as_chunks::<2>().0 {
        assert!(mark[0].is_finite() && mark[1].is_finite());
        assert!((mark[0] / limit - (2.0 * reference.uniform() - 1.0)).abs() < 4.0e-16);
        assert!((mark[1] / limit - (1.5 * reference.uniform() - 1.0)).abs() < 4.0e-16);
    }
    assert_eq!(stream.next_unit(), reference.uniform());
}

#[test]
fn empirical_moments_match_normal_uniform_and_asymmetric_exponential_laws() {
    let mixture_mean = 0.3_f64 / 2.0 - 0.7_f64 / 4.0;
    let mixture_var = 0.3 * 2.0 / 4.0 + 0.7 * 2.0 / 16.0 - mixture_mean.powi(2);
    let cases = [
        (
            JumpDistribution::Normal {
                mean: vec![1.3],
                std: vec![0.7],
            },
            1.3,
            0.49,
        ),
        (
            JumpDistribution::Uniform {
                low: vec![-2.0],
                high: vec![1.0],
            },
            -0.5,
            0.75,
        ),
        (
            JumpDistribution::DoubleExponential {
                p_up: 0.3,
                eta_up: 2.0,
                eta_down: 4.0,
                dim: 1,
            },
            mixture_mean,
            mixture_var,
        ),
    ];
    let mut sampled = vec![0.0; 50_000];
    for (distribution, expected_mean, expected_var) in cases {
        distribution
            .fill_marks(&mut Stream::new(19), &mut sampled)
            .unwrap();
        let mean = sampled.iter().sum::<f64>() / sampled.len() as f64;
        let variance =
            sampled.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / sampled.len() as f64;
        assert!(
            (mean - expected_mean).abs() < 0.012,
            "mean {mean} vs {expected_mean}"
        );
        assert!(
            (variance - expected_var).abs() < 0.018,
            "variance {variance} vs {expected_var}"
        );
    }
}

#[test]
fn malformed_mark_requests_are_rejected_before_draws() {
    let distribution = JumpDistribution::Constant(vec![1.0; 2]);
    let mut stream = Stream::new(42);
    let mut output = [5.0; 3];
    assert!(matches!(
        distribution.fill_marks(&mut stream, &mut output),
        Err(SdeError::Shape(_))
    ));
    assert_eq!(output, [5.0; 3]);
    let invalid = [
        JumpDistribution::Constant(vec![]),
        JumpDistribution::Normal {
            mean: vec![0.0],
            std: vec![-1.0],
        },
        JumpDistribution::Uniform {
            low: vec![1.0],
            high: vec![0.0],
        },
        JumpDistribution::DoubleExponential {
            p_up: 0.5,
            eta_up: f64::INFINITY,
            eta_down: 1.0,
            dim: 1,
        },
        JumpDistribution::AlphaStable {
            alpha: 0.0,
            beta: 0.0,
            scale: 1.0,
            location: 0.0,
            dim: 1,
        },
    ];
    for distribution in invalid {
        assert!(matches!(
            distribution.fill_marks(&mut stream, &mut output),
            Err(SdeError::Invalid(_))
        ));
        assert_eq!(output, [5.0; 3]);
    }
    assert_eq!(stream, Stream::new(42));
    distribution.fill_marks(&mut stream, &mut []).unwrap();
    assert_eq!(stream, Stream::new(42));
    let impossible = JumpDistribution::DoubleExponential {
        p_up: 0.5,
        eta_up: 1.0,
        eta_down: 1.0,
        dim: usize::MAX,
    };
    assert!(matches!(
        impossible.sample(&mut stream),
        Err(SdeError::Invalid(_))
    ));
    assert_eq!(stream, Stream::new(42));
}

#[test]
fn fixed_jump_diffusion_reuses_marks_without_reordering_draws_or_callbacks() {
    let (seed, dt, steps) = (29, 0.25, 80);
    let mut spec = JumpSpec {
        intensity: 14.0,
        distribution: JumpDistribution::Normal {
            mean: vec![0.8, 0.2],
            std: vec![0.45, 0.9],
        },
        compensated: false,
    };
    let options = Options::fixed(
        Method::EulerMaruyama,
        Structure::Diagonal,
        seed,
        0.0,
        dt,
        steps,
    );
    let mut callbacks = Vec::new();
    let zero = |_: f64, _: &[f64]| Ok(vec![0.0; 2]);
    let initial = [1.0, 2.0];
    let path = integrate_jump_diffusion(
        zero,
        zero,
        |t, y, mark| {
            callbacks.push((t, y.to_vec(), mark.to_vec()));
            Ok(mark.to_vec())
        },
        &initial,
        &spec,
        &options,
    )
    .unwrap();
    let mut reference = ReferenceStream(ensemble_seed(seed, 0x4A55_4D50));
    let mut state = initial;
    let mut calls = callbacks.iter();
    for (step, &actual_jumps) in path.jump_counts.iter().enumerate() {
        let expected_jumps = reference.poisson_small(spec.intensity * dt);
        assert_eq!(actual_jumps, expected_jumps);
        for _ in 0..expected_jumps {
            let expected_mark = [
                0.8 + 0.45 * reference.normal(),
                0.2 + 0.9 * reference.normal(),
            ];
            let (time, pre_jump_state, actual_mark) = calls.next().unwrap();
            assert_eq!(*time, step as f64 * dt);
            assert_eq!(pre_jump_state, &state);
            assert_eq!(actual_mark, &expected_mark);
            state[0] += expected_mark[0];
            state[1] += expected_mark[1];
        }
        assert_eq!(path.path.y[step + 1], state);
    }
    assert!(calls.next().is_none());
    assert_eq!(path.total_jumps, callbacks.len());

    spec.intensity = 0.0;
    spec.distribution = JumpDistribution::Constant(vec![]);
    let no_jumps = integrate_jump_diffusion(
        zero,
        zero,
        |_, _, _| unreachable!(),
        &initial,
        &spec,
        &options,
    )
    .unwrap();
    assert_eq!(no_jumps.total_jumps, 0);
    assert_eq!(no_jumps.path.last().unwrap(), &initial);
}
