//! Split-step Schrodinger solver (ported from `physics_sim_schrodinger_test.rs`).

use num_complex::Complex;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::schrodinger_quantum::*;
use std::f64::consts::PI;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 12,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn params(
    n: usize,
    l: f64,
    dt: f64,
    steps: usize,
    potential: Vec<f64>,
) -> SchrodingerParameters {
    SchrodingerParameters {
        nx: n,
        ny: n,
        lx: l,
        ly: l,
        dt,
        time_steps: steps,
        hbar: 1.0,
        mass: 1.0,
        potential,
    }
}

fn gaussian(
    n: usize,
    l: f64,
    x0: f64,
    y0: f64,
    sigma: f64,
    kx: f64,
) -> Vec<Complex<f64>> {
    let d = l / n as f64;
    let mut psi = vec![Complex::default(); n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (i as f64 * d, j as f64 * d);
            let env = (-((x - x0).powi(2) + (y - y0).powi(2)) / (2.0 * sigma * sigma)).exp();
            psi[j * n + i] = Complex::from_polar(env, kx * x);
        }
    }
    psi
}

#[test]
fn snapshots_are_taken_every_tenth_step_with_the_grid_shape() {
    let n = 32;
    let p = params(n, 10.0, 0.1, 25, vec![0.0; n * n]);
    let mut psi = vec![Complex::new(1.0, 0.0); n * n];
    let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(snaps.len(), 3); // steps 0, 10, 20
    for s in &snaps {
        assert_eq!((s.nrows(), s.ncols()), (n, n));
    }
}

#[test]
fn uniform_wavefunction_stays_uniform_with_unit_density() {
    let n = 16;
    let p = params(n, 5.0, 0.05, 11, vec![0.0; n * n]);
    let mut psi = vec![Complex::new(1.0, 0.0); n * n];
    let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
    for s in &snaps {
        assert!(s.iter().all(|&v| (v - 1.0).abs() < 1e-9));
    }
}

#[test]
fn total_probability_is_conserved_in_a_potential() {
    let n = 32;
    let l = 32.0;
    let d = l / n as f64;
    let mut pot = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (i as f64 * d - 16.0, j as f64 * d - 16.0);
            pot[j * n + i] = 0.02 * (x * x + y * y);
        }
    }
    let p = params(n, l, 0.05, 41, pot);
    let mut psi = gaussian(n, l, 16.0, 16.0, 2.0, 0.5);
    let norm0: f64 = psi.iter().map(|c| c.norm_sqr()).sum();
    let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
    for s in &snaps {
        let norm: f64 = s.iter().sum();
        assert!((norm - norm0).abs() < 1e-8 * norm0, "{norm} vs {norm0}");
    }
}

#[test]
fn free_packet_moves_with_group_velocity_k_over_m() {
    let (n, l) = (64usize, 64.0);
    let kx = 2.0 * PI * 4.0 / l; // resolved on the periodic grid
    let p = params(n, l, 1.0, 21, vec![0.0; n * n]);
    let mut psi = gaussian(n, l, 20.0, 32.0, 3.0, kx);
    let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(snaps.len(), 3);
    let centroid = |s: &ndarray::Array2<f64>| {
        let (mut sx, mut sy, mut m) = (0.0, 0.0, 0.0);
        for ((j, i), &v) in s.indexed_iter() {
            sx += v * i as f64;
            sy += v * j as f64;
            m += v;
        }
        (sx / m, sy / m)
    };
    let (x1, y1) = centroid(&snaps[1]);
    let (x2, y2) = centroid(&snaps[2]);
    // Ten time units of unit-mass motion at k = 2 pi 4 / 64 = 0.3927.
    assert!((x2 - x1 - 10.0 * kx).abs() < 0.05, "dx = {}", x2 - x1);
    assert!((y2 - y1).abs() < 1e-6);
}

#[test]
fn a_hard_barrier_reflects_most_of_the_packet() {
    let (n, l) = (64usize, 64.0);
    let mut pot = vec![0.0; n * n];
    for j in 0..n {
        for i in 40..44 {
            pot[j * n + i] = 1e3;
        }
    }
    let p = params(n, l, 0.5, 41, pot);
    let mut psi = gaussian(n, l, 20.0, 32.0, 3.0, 0.4);
    let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
    let last = snaps.last().unwrap_or_else(|| panic!("no snapshots"));
    let behind: f64 = last
        .indexed_iter()
        .filter(|((_, i), _)| *i >= 44 && *i < 60)
        .map(|(_, &v)| v)
        .sum();
    let total: f64 = last.iter().sum();
    assert!(
        behind / total < 0.05,
        "transmitted fraction {}",
        behind / total
    );
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_density_is_finite_and_nonnegative(steps in 1usize..10, dt in 0.01f64..0.1) {
        let n = 16;
        let p = params(n, 5.0, dt, steps, vec![0.0; n * n]);
        let mut psi = vec![Complex::new(0.5, 0.5); n * n];
        let snaps = run_schrodinger_simulation(&p, &mut psi).unwrap_or_else(|e| panic!("{e}"));
        for s in &snaps {
            prop_assert!(s.iter().all(|v| v.is_finite() && *v >= 0.0));
        }
    }
}
