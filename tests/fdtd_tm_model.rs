//! Analytic and deterministic checks of the fixed two-dimensional TM model.

use ndarray::Array2;
use rssn::sim::models_fdtd_electrodynamics::FdtdParameters;
use rssn::sim::models_fdtd_electrodynamics::run_fdtd_final_field;
use rssn::sim::models_fdtd_electrodynamics::run_fdtd_simulation;
use std::f64::consts::TAU;

fn parameters() -> FdtdParameters {
    FdtdParameters {
        width: 45,
        height: 49,
        time_steps: 12,
        source_pos: (22, 24),
        source_freq: 0.125,
    }
}

fn pulse(
    step: usize,
    freq: f64,
) -> f64 {
    let t = step as f64;
    (-((t - 30.0).powi(2)) / 100.0).exp() * (TAU * freq * t).sin()
}

fn equal_bits(
    a: &Array2<f64>,
    b: &Array2<f64>,
) {
    assert_eq!(a.dim(), b.dim());
    for (&left, &right) in a.iter().zip(b) {
        assert_eq!(left.to_bits(), right.to_bits(), "{left} vs {right}");
    }
}

#[test]
fn thin_grids_have_only_the_damped_soft_source() {
    let steps = 18;
    let freq = 0.125;
    let mut scalar = 0.0;
    let timeline: Vec<_> = (0..steps)
        .map(|t| {
            scalar = 0.95 * (scalar + pulse(t, freq));
            scalar
        })
        .collect();

    for (width, height, source_pos) in [
        (1, 1, (0, 0)),
        (1, 43, (0, 22)),
        (43, 1, (22, 0)),
        (2, 32, (1, 31)),
        (32, 2, (31, 1)),
    ] {
        let p = FdtdParameters {
            width,
            height,
            time_steps: steps,
            source_pos,
            source_freq: freq,
        };
        let final_field = run_fdtd_final_field(&p).unwrap();
        for (site, &field) in final_field.indexed_iter() {
            let expected = if site == source_pos {
                timeline[steps - 1]
            } else {
                0.0
            };
            assert_eq!(field.to_bits(), expected.to_bits());
        }
        let sampled = run_fdtd_simulation(&p);
        assert_eq!(sampled.len(), steps.div_ceil(5));
        for (k, snapshot) in sampled.iter().enumerate() {
            assert_eq!(snapshot[source_pos].to_bits(), timeline[k * 5].to_bits());
        }
    }
}

#[test]
fn the_first_wavefront_undergoes_one_damping_factor_at_the_belt() {
    let p = FdtdParameters {
        width: 25,
        height: 26,
        time_steps: 3,
        source_pos: (10, 10),
        ..parameters()
    };
    let e = run_fdtd_final_field(&p).unwrap();
    let neighbor = 0.25 * pulse(1, p.source_freq);
    assert_eq!(e[[9, 10]].to_bits(), (neighbor * 0.95).to_bits());
    assert_eq!(e[[10, 9]].to_bits(), (neighbor * 0.95).to_bits());
    assert_eq!(e[[11, 10]].to_bits(), neighbor.to_bits());
    assert_eq!(e[[10, 11]].to_bits(), neighbor.to_bits());
    assert_eq!(e[[9, 9]].to_bits(), 0.0_f64.to_bits());
    assert!((e[[10, 10]] - pulse(2, p.source_freq)).abs() < 1e-18);
}

#[test]
fn the_interior_matches_the_discrete_wave_equation_after_eliminating_h() {
    let p = parameters();
    let shape = (p.width, p.height);
    let mut current = Array2::<f64>::zeros(shape);
    let mut previous = current.clone();
    let mut last_pulse = 0.0;
    for t in 0..p.time_steps {
        let mut next = Array2::<f64>::zeros(shape);
        for i in 1..p.width - 1 {
            for j in 1..p.height - 1 {
                let center = current[[i, j]];
                let laplacian = current[[i - 1, j]]
                    + current[[i + 1, j]]
                    + current[[i, j - 1]]
                    + current[[i, j + 1]]
                    - 4.0 * center;
                next[[i, j]] = 2.0 * center - previous[[i, j]] + 0.25 * laplacian;
            }
        }
        let excitation = pulse(t, p.source_freq);
        next[p.source_pos] += excitation - last_pulse;
        last_pulse = excitation;
        previous = current;
        current = next;
    }
    let observed = run_fdtd_final_field(&p).unwrap();
    for (&got, &want) in observed.iter().zip(&current) {
        assert!((got - want).abs() < 3e-16, "{got} vs {want}");
    }
}

#[test]
fn sampled_times_and_final_update_keep_their_independent_contracts() {
    let mut p = parameters();
    p.time_steps = 17;
    let snapshots = run_fdtd_simulation(&p);
    assert_eq!(snapshots.len(), 4);
    for (sample, array) in snapshots.iter().enumerate() {
        let at = FdtdParameters {
            time_steps: 5 * sample + 1,
            ..p
        };
        equal_bits(array, &run_fdtd_final_field(&at).unwrap());
    }
    let final_field = run_fdtd_final_field(&p).unwrap();
    assert!(
        final_field
            .iter()
            .zip(&snapshots[3])
            .any(|(&a, &b)| a.to_bits() != b.to_bits())
    );

    p.time_steps = 0;
    assert!(run_fdtd_simulation(&p).is_empty());
    let no_steps = run_fdtd_final_field(&p).unwrap();
    assert_eq!(no_steps.dim(), (p.width, p.height));
    assert!(no_steps.iter().all(|v| v.to_bits() == 0));
}

#[test]
fn discrete_frequency_aliases_and_final_input_errors_are_explicit() {
    let p = parameters();
    let ordinary = run_fdtd_final_field(&p).unwrap();
    let alias = FdtdParameters {
        source_freq: 1_099_511_627_776.125,
        ..p
    };
    equal_bits(&ordinary, &run_fdtd_final_field(&alias).unwrap());
    let integer = FdtdParameters {
        source_freq: 1e308,
        ..p
    };
    assert!(
        run_fdtd_final_field(&integer)
            .unwrap()
            .iter()
            .all(|v| v.to_bits() == 0)
    );
    for e in run_fdtd_simulation(&integer) {
        assert!(e.iter().all(|v| v.to_bits() == 0));
    }
    for bad in [
        FdtdParameters { width: 0, ..p },
        FdtdParameters { height: 0, ..p },
        FdtdParameters {
            source_pos: (p.width, 0),
            ..p
        },
        FdtdParameters {
            width: usize::MAX,
            height: 3,
            ..p
        },
        FdtdParameters {
            width: isize::MAX as usize,
            height: 1,
            source_pos: (0, 0),
            ..p
        },
        FdtdParameters {
            source_freq: f64::NAN,
            ..p
        },
        FdtdParameters {
            source_freq: f64::INFINITY,
            ..p
        },
    ] {
        assert!(run_fdtd_final_field(&bad).is_err(), "{bad:?}");
    }
}

#[test]
fn active_workers_produce_bitwise_equal_fields_and_samples() {
    let p = FdtdParameters {
        width: 125,
        height: 135,
        source_pos: (12, 11),
        time_steps: 23,
        ..parameters()
    };
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| (run_fdtd_simulation(&p), run_fdtd_final_field(&p).unwrap()))
    };
    let (one, final_one) = run(1);
    let (two, final_two) = run(2);
    equal_bits(&final_one, &final_two);
    assert_eq!(one.len(), two.len());
    for (a, b) in one.iter().zip(&two) {
        equal_bits(a, b);
    }
}
