//! 2D Ising model Monte Carlo (ported from `physics_sim_ising_test.rs`).
//!
//! `run_ising_simulation` uses a fixed default seed; the `_seeded` variant takes an explicit one.

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::ising_statistical::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 16,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn run(
    w: usize,
    h: usize,
    t: f64,
    steps: usize,
) -> (Vec<i8>, f64) {
    run_ising_simulation(&IsingParameters {
        width: w,
        height: h,
        temperature: t,
        mc_steps: steps,
    })
}

/// Energy per spin with periodic boundaries (J = 1).
fn energy_per_spin(
    g: &[i8],
    w: usize,
    h: usize,
) -> f64 {
    let mut e = 0.0;
    for i in 0..h {
        for j in 0..w {
            let s = f64::from(g[i * w + j]);
            e -= s * f64::from(g[i * w + (j + 1) % w]);
            e -= s * f64::from(g[((i + 1) % h) * w + j]);
        }
    }
    e / (w * h) as f64
}

#[test]
fn spins_are_plus_or_minus_one_and_magnetisation_matches_the_grid() {
    let (grid, mag) = run(10, 12, 2.0, 20);
    assert_eq!(grid.len(), 120);
    assert!(grid.iter().all(|&s| s == 1 || s == -1));
    let m = grid.iter().map(|&s| f64::from(s)).sum::<f64>() / 120.0;
    assert!((mag - m.abs()).abs() < 1e-12);
}

#[test]
fn zero_steps_returns_a_random_initial_configuration() {
    let (grid, mag) = run(30, 30, 1.0, 0);
    assert_eq!(grid.len(), 900);
    // 900 fair coin flips: |m| < 0.25 with probability > 1 - 1e-8.
    assert!(mag < 0.25, "mag = {mag}");
}

#[test]
fn quench_to_low_temperature_lowers_the_energy() {
    // A random state has E / N ~ 0; the ground state has -2.  Even metastable stripe states have
    // E / N = -1.6 on a 10 x 10 lattice.
    let (grid, _) = run(10, 10, 0.1, 100);
    let e = energy_per_spin(&grid, 10, 10);
    assert!(e < -1.0, "E/N = {e}");
}

#[test]
fn runs_are_reproducible_and_seed_dependent() {
    let params = IsingParameters {
        width: 16,
        height: 16,
        temperature: 2.0,
        mc_steps: 30,
    };
    assert_eq!(run_ising_simulation(&params), run_ising_simulation(&params));
    assert_eq!(
        run_ising_simulation(&params),
        run_ising_simulation_seeded(&params, DEFAULT_ISING_SEED)
    );
    let a = run_ising_simulation_seeded(&params, 1);
    assert_eq!(a, run_ising_simulation_seeded(&params, 1));
    assert_ne!(a.0, run_ising_simulation_seeded(&params, 2).0);
}

#[test]
fn scenario_writes_its_csv_into_the_given_directory() {
    let dir = std::env::temp_dir().join(format!("rssn_ising_{}", std::process::id()));
    // The .npy outputs need the `npy` feature; without it the call reports an error.
    let res = simulate_ising_phase_transition_scenario(&dir);
    if cfg!(feature = "npy") {
        res.unwrap_or_else(|e| panic!("{e}"));
        assert!(dir.join("ising_magnetization_vs_temp.csv").is_file());
    } else {
        assert!(res.is_err());
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn high_temperature_is_disordered() {
    let (_, mag) = run(20, 20, 100.0, 50);
    assert!(mag < 0.3, "mag = {mag}");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_ising_state_is_valid(t in 0.1f64..5.0, steps in 1usize..10) {
        let (grid, mag) = run(8, 8, t, steps);
        prop_assert_eq!(grid.len(), 64);
        prop_assert!((0.0..=1.0).contains(&mag));
        prop_assert!(grid.iter().all(|&s| s == 1 || s == -1));
    }
}
