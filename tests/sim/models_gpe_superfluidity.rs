//! Gross-Pitaevskii ground state (ported from `physics_sim_gpe_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::gpe_superfluidity::*;

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
    d_tau: f64,
    steps: usize,
    g: f64,
) -> GpeParameters {
    GpeParameters {
        nx: n,
        ny: n,
        lx: l,
        ly: l,
        d_tau,
        time_steps: steps,
        g,
        trap_strength: 1.0,
    }
}

#[test]
fn density_has_the_grid_shape_and_is_normalised() {
    let p = params(32, 10.0, 0.1, 20, 100.0);
    let rho = run_gpe_ground_state_finder(&p).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!((rho.nrows(), rho.ncols()), (32, 32));
    let dx = p.lx / p.nx as f64;
    // The code normalises sum(|psi|^2) dx dy to nx * ny.
    let total = rho.iter().sum::<f64>() * dx * dx;
    assert!((total - 1024.0).abs() < 1e-6, "total = {total}");
}

#[test]
fn noninteracting_ground_state_is_the_harmonic_oscillator_gaussian() {
    // g = 0, omega = 1: rho(x, y) = A exp(-(x^2 + y^2)), A = N / pi with the code's normalisation.
    let n = 64;
    let p = params(n, 10.0, 0.05, 400, 0.0);
    let rho = run_gpe_ground_state_finder(&p).unwrap_or_else(|e| panic!("{e}"));
    let dx = 10.0 / n as f64;
    let a = (n * n) as f64 / std::f64::consts::PI;
    for (j, i) in [(32usize, 32usize), (32, 40), (40, 32), (36, 36), (32, 45)] {
        let (x, y) = ((i as f64 - 32.0) * dx, (j as f64 - 32.0) * dx);
        let exact = a * (-(x * x + y * y)).exp();
        let rel = (rho[[j, i]] - exact).abs() / exact;
        assert!(rel < 0.02, "({j},{i}): {} vs {exact}", rho[[j, i]]);
    }
}

#[test]
fn interaction_broadens_the_condensate() {
    let n = 32;
    let free = run_gpe_ground_state_finder(&params(n, 10.0, 0.05, 100, 0.0))
        .unwrap_or_else(|e| panic!("{e}"));
    let strong = run_gpe_ground_state_finder(&params(n, 10.0, 0.05, 100, 100.0))
        .unwrap_or_else(|e| panic!("{e}"));
    let second_moment = |a: &ndarray::Array2<f64>| {
        let mut s = 0.0;
        for ((j, i), &v) in a.indexed_iter() {
            let (x, y) = (i as f64 - 16.0, j as f64 - 16.0);
            s += v * (x * x + y * y);
        }
        s / a.iter().sum::<f64>()
    };
    assert!(second_moment(&strong) > second_moment(&free));
    // Peak density is lower for the repulsive condensate.
    let max = |a: &ndarray::Array2<f64>| a.iter().cloned().fold(0.0, f64::max);
    assert!(max(&strong) < max(&free));
}

#[test]
fn density_is_symmetric_about_the_trap_centre() {
    let n = 32;
    let rho = run_gpe_ground_state_finder(&params(n, 10.0, 0.05, 100, 0.5))
        .unwrap_or_else(|e| panic!("{e}"));
    for j in 0..n {
        for i in 0..n {
            let (mi, mj) = ((n - i) % n, (n - j) % n);
            assert!((rho[[j, i]] - rho[[j, mi]]).abs() < 1e-8 * (1.0 + rho[[j, i]]));
            assert!((rho[[j, i]] - rho[[i, j]]).abs() < 1e-8 * (1.0 + rho[[j, i]]));
            assert!((rho[[j, i]] - rho[[mj, i]]).abs() < 1e-8 * (1.0 + rho[[j, i]]));
        }
    }
    let max = rho.iter().cloned().fold(0.0, f64::max);
    assert!((rho[[16, 16]] - max).abs() < 1e-9 * max);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_gpe_density_is_finite_and_nonnegative(steps in 1usize..5, d_tau in 0.01f64..0.05) {
        let rho = run_gpe_ground_state_finder(&params(16, 5.0, d_tau, steps, 50.0)).unwrap_or_else(|e| panic!("{e}"));
        for &v in rho.iter() {
            prop_assert!(v.is_finite() && v >= 0.0);
        }
    }
}
