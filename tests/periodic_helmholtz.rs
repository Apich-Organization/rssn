//! Manufactured solutions and scale checks for periodic Helmholtz plans.

use std::f64::consts::TAU;

use num_complex::Complex;
use rssn::kernels::fft::Norm;
use rssn::kernels::fft::fft_axes;
use rssn::sim::spectral::SpectralGrid;
use rssn::sim::spectral::fft_nd;

#[test]
fn positive_small_shift_and_long_period_modes_remain_invertible() {
    let grid = SpectralGrid::new(&[32], &[TAU]).unwrap();
    for sigma in [1e-300, f64::from_bits(1)] {
        let constant = grid.helmholtz(sigma, &vec![sigma; grid.size()]).unwrap();
        assert!(constant.iter().all(|&value| (value - 1.0).abs() < 2e-14));
    }

    let lambda = grid.k2(1);
    let sigma = (-lambda).next_up();
    let mode: Vec<_> = (0..grid.size()).map(|p| grid.coord(p)[0].cos()).collect();
    let forcing: Vec<_> = mode.iter().map(|value| (lambda + sigma) * value).collect();
    let near_resonant = grid.helmholtz(sigma, &forcing).unwrap();
    assert!(
        near_resonant
            .iter()
            .zip(&mode)
            .all(|(a, b)| (a - b).abs() < 2e-14)
    );

    let long = SpectralGrid::new(&[32], &[TAU * 1e8]).unwrap();
    let source: Vec<_> = mode.iter().map(|value| long.k2(1) * value).collect();
    let poisson = long.helmholtz(0.0, &source).unwrap();
    assert!(
        poisson
            .iter()
            .zip(&mode)
            .all(|(a, b)| (a - b).abs() < 2e-14)
    );
}

#[test]
fn singular_compatibility_scales_with_the_forcing() {
    let grid = SpectralGrid::new(&[64], &[TAU]).unwrap();
    let plan = grid.helmholtz_plan(0.0).unwrap();
    for scale in [1e-250, 1.0, 1e250] {
        let tiny_mean = vec![scale; grid.size()];
        assert!(plan.solve(&tiny_mean).unwrap_err().contains("singular"));
        let source: Vec<_> = (0..grid.size())
            .map(|p| scale * grid.coord(p)[0].sin())
            .collect();
        let compatible = plan.solve(&source).unwrap();
        assert!(
            compatible.iter().zip(source).all(|(obtained, expected)| {
                (obtained / scale - expected / scale).abs() < 3e-14
            })
        );
    }
    let resonance = grid.helmholtz_plan(-grid.k2(1)).unwrap();
    let source: Vec<_> = (0..grid.size())
        .map(|p| grid.coord(p)[0].cos() * 1e-250)
        .collect();
    assert!(resonance.solve(&source).unwrap_err().contains("singular"));
    let compatible: Vec<_> = (0..grid.size())
        .map(|p| (2.0 * grid.coord(p)[0]).sin())
        .collect();
    let solution = resonance.solve(&compatible).unwrap();
    assert!(
        solution
            .iter()
            .zip(compatible)
            .all(|(actual, source)| (3.0 * actual - source).abs() < 2e-14)
    );
    assert_eq!(
        resonance.solve(&vec![0.0; grid.size()]).unwrap(),
        vec![0.0; grid.size()]
    );
}

#[test]
fn plan_keys_track_box_length_and_shift_and_keep_existing_plans_valid() {
    let mut grid = SpectralGrid::new(&[32], &[TAU]).unwrap();
    let source: Vec<_> = (0..grid.size()).map(|p| grid.coord(p)[0].cos()).collect();
    let unit = grid.helmholtz_plan(1.0).unwrap();
    let shifted = grid.helmholtz_plan(2.0).unwrap();
    grid.len[0] *= 2.0;
    let extended = grid.helmholtz_plan(1.0).unwrap();
    for (plan, divisor) in [(unit, 2.0), (shifted, 3.0), (extended, 1.25)] {
        assert!(
            plan.solve(&source)
                .unwrap()
                .iter()
                .zip(&source)
                .all(|(actual, source)| (divisor * actual - source).abs() < 2e-14)
        );
    }
}

#[test]
fn large_uniform_forcing_uses_a_representable_transform() {
    let grid = SpectralGrid::new(&[16, 8], &[TAU, TAU]).unwrap();
    let forcing = 0.75 * f64::MAX;
    let solution = grid.helmholtz(4.0, &vec![forcing; grid.size()]).unwrap();
    assert!(
        solution
            .iter()
            .all(|value| value.to_bits() == (forcing / 4.0).to_bits())
    );
    assert!(grid.helmholtz(0.25, &vec![forcing; grid.size()]).is_err());
}

#[test]
fn batched_three_dimensional_modes_keep_order_across_thread_counts() {
    let grid = SpectralGrid::new(&[16, 8, 4], &[TAU, 2.0 * TAU, TAU]).unwrap();
    let sigma = 1.75;
    let mode: Vec<_> = (0..grid.size())
        .map(|p| {
            let [x, y, z] = grid.coord(p);
            (2.0 * x).sin() * (0.5 * y).cos() * z.sin()
        })
        .collect();
    let samples: Vec<Vec<f64>> = (1..=16)
        .map(|sample| {
            mode.iter()
                .map(|&value| 7.0 * f64::from(sample) * value + 0.25 * sigma)
                .collect()
        })
        .collect();
    let plan = grid.helmholtz_plan(sigma).unwrap();
    let run_at = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| plan.solve_batch(&samples).unwrap())
    };
    let serial = run_at(1);
    let parallel = run_at(2);
    let mut scale = 1.0;
    for (source, (single, multi)) in samples.iter().zip(serial.iter().zip(&parallel)) {
        let one = plan.solve(source).unwrap();
        for (((a, b), c), expected) in single.iter().zip(multi).zip(one).zip(&mode) {
            assert_eq!(a.to_bits(), b.to_bits());
            assert_eq!(a.to_bits(), c.to_bits());
            assert!((a - (scale * expected + 0.25)).abs() < 4e-13);
        }
        scale += 1.0;
    }
    assert!(plan.solve_batch(&[]).unwrap().is_empty());
    assert!(
        plan.solve_batch(&[vec![f64::NAN; grid.size()], vec![]])
            .unwrap_err()
            .contains("finite")
    );
    assert!(
        plan.solve_batch(&[vec![], vec![f64::NAN; grid.size()]])
            .unwrap_err()
            .contains("length")
    );
}

#[test]
fn malformed_and_unrepresentable_plans_return_errors() {
    let grid = SpectralGrid::new(&[4], &[TAU]).unwrap();
    assert!(grid.helmholtz_plan(f64::INFINITY).is_err());
    assert!(grid.helmholtz(1.0, &[1.0, 2.0]).is_err());
    for length in [0.0, -1.0, f64::NAN, f64::MAX, f64::MIN_POSITIVE] {
        let malformed = SpectralGrid::new(&[4], &[length]).unwrap();
        assert!(malformed.helmholtz_plan(0.0).is_err());
    }
    let mut malformed = grid;
    malformed.dim = 4;
    assert!(malformed.helmholtz_plan(0.0).is_err());
    malformed.dim = 3;
    malformed.n = [1 << (usize::BITS - 2); 3];
    assert!(malformed.helmholtz_plan(0.0).is_err());
}

#[test]
fn shared_axis_fft_matches_the_legacy_x_fastest_transform() {
    let input: Vec<_> = (0..8 * 4 * 16)
        .map(|index: i32| Complex::new(f64::from(index % 17) * 0.03, f64::from(index % 13) * 0.02))
        .collect();
    let mut legacy = input.clone();
    fft_nd(&mut legacy, &[8, 4, 16], false).unwrap();
    let planned = fft_axes(&input, &[16, 4, 8], &[2, 1, 0], false, Norm::Backward).unwrap();
    for (actual, reference) in planned.iter().zip(&legacy) {
        assert!((*actual - *reference).norm() < 3e-13);
    }
    let recovered = fft_axes(&planned, &[16, 4, 8], &[2, 1, 0], true, Norm::Backward).unwrap();
    for (actual, reference) in recovered.iter().zip(&input) {
        assert!((*actual - *reference).norm() < 2e-14);
    }
    assert_eq!(
        fft_axes(&input, &[1, input.len()], &[0], false, Norm::Backward).unwrap(),
        input
    );
    assert_eq!(
        fft_axes(&input, &[input.len()], &[], true, Norm::Backward).unwrap(),
        input
    );
}
