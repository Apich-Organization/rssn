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

#[test]
fn lattice_simulation_covers_any_dimension_and_every_model() {
    // 1-D Ising ring: e = -tanh(beta) per site
    let r = run_lattice_simulation(
        [200],
        true,
        SpinModel::Ising,
        SpinAlgorithm::HeatBath,
        1.5,
        (500, 6000),
        11,
    )
    .unwrap();
    assert!((r.energy + (1.0f64 / 1.5).tanh()).abs() < 0.02, "{r:?}");
    // 2-D Ising above Tc against Onsager
    let r = run_lattice_simulation(
        [16, 16],
        true,
        SpinModel::Ising,
        SpinAlgorithm::Wolff,
        3.0,
        (500, 8000),
        12,
    )
    .unwrap();
    let exact = rssn::sim::spin_models::onsager_energy_per_site(3.0, 1.0);
    assert!((r.energy - exact).abs() < 0.02, "{} vs {exact}", r.energy);
    // 4-D Ising is ordered at T = 3 (Tc(4D) = 6.68)
    let r = run_lattice_simulation(
        [4, 4, 4, 4],
        true,
        SpinModel::Ising,
        SpinAlgorithm::Wolff,
        3.0,
        (300, 600),
        13,
    )
    .unwrap();
    assert!(r.magnetization > 0.8 && r.energy < -3.0, "{r:?}");
    // XY and Heisenberg heat bath in 1-D match their exact chains
    let beta = 1.5;
    let r = run_lattice_simulation(
        [200],
        true,
        SpinModel::Heisenberg,
        SpinAlgorithm::HeatBath,
        1.0 / beta,
        (500, 6000),
        14,
    )
    .unwrap();
    assert!(
        (r.energy + (1.0 / beta.tanh() - 1.0 / beta)).abs() < 0.02,
        "{r:?}"
    );
    // Potts q = 3 low temperature is ordered; cluster algorithm with the cold side of Tc
    let r = run_lattice_simulation(
        [12, 12],
        true,
        SpinModel::Potts(3),
        SpinAlgorithm::SwendsenWang,
        0.5,
        (200, 800),
        15,
    )
    .unwrap();
    assert!(r.magnetization > 0.95, "{r:?}");
    // open boundaries and invalid input
    assert!(
        run_lattice_simulation(
            [6, 6],
            false,
            SpinModel::Xy,
            SpinAlgorithm::Wolff,
            0.5,
            (10, 10),
            1
        )
        .is_ok()
    );
    assert!(
        run_lattice_simulation(
            [6, 6],
            true,
            SpinModel::Ising,
            SpinAlgorithm::HeatBath,
            -1.0,
            (1, 1),
            1
        )
        .is_err()
    );
}

#[test]
fn wolff_observables_come_from_the_generic_engine_and_reject_tiny_lattices() {
    let o = run_wolff_simulation(2, 2, 2.0, 10, 10, 1);
    assert!(o.energy.is_nan());
    let o = run_wolff_simulation(16, 16, 1.2, 100, 500, 2);
    assert!(o.magnetization > 0.97 && o.mean_cluster > 0.9);
}

#[test]
fn checkerboard_is_correct_on_odd_periodic_lattices() {
    use rssn::sim::spin_models::enumerate_ising_dos;
    // exact canonical energy per spin by enumeration vs the long-run average of the checkerboard
    for (w, h) in [(3usize, 3usize), (5, 3), (3, 4)] {
        let t = 2.5;
        let dos = enumerate_ising_dos::<2>([w, h], true).unwrap();
        let z: f64 = dos.iter().map(|&(e, c)| c * (-(e as f64) / t).exp()).sum();
        let exact = dos
            .iter()
            .map(|&(e, c)| c * e as f64 * (-(e as f64) / t).exp())
            .sum::<f64>()
            / z
            / (w * h) as f64;
        let mut acc = 0.0;
        let runs = 1500;
        for s in 0..runs {
            let (g, _) = run_ising_simulation_seeded(
                &IsingParameters {
                    width: w,
                    height: h,
                    temperature: t,
                    mc_steps: 15,
                },
                1000 + s,
            );
            acc += energy_per_spin(&g, w, h);
        }
        let mean = acc / runs as f64;
        assert!(
            (mean - exact).abs() < 0.04,
            "{w}x{h}: {mean} vs exact {exact}"
        );
    }
    // invalid input does not panic
    let (g, m) = run_ising_simulation_seeded(
        &IsingParameters {
            width: 0,
            height: 4,
            temperature: 1.0,
            mc_steps: 3,
        },
        1,
    );
    assert!(g.is_empty() && m.is_nan());
    let (_, m) = run_ising_simulation_seeded(
        &IsingParameters {
            width: 4,
            height: 4,
            temperature: -1.0,
            mc_steps: 3,
        },
        1,
    );
    assert!(m.is_nan());
}
