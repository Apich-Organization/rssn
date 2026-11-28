//! 2D FDTD electrodynamics model (ported from `physics_sim_fdtd_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::fdtd_electrodynamics::*;

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
    steps: usize,
    src: (usize, usize),
    f: f64,
) -> FdtdParameters {
    FdtdParameters {
        width: n,
        height: n,
        time_steps: steps,
        source_pos: src,
        source_freq: f,
    }
}

#[test]
fn snapshots_are_taken_every_fifth_step() {
    let snaps = run_fdtd_simulation(&params(50, 100, (25, 25), 0.1));
    assert_eq!(snaps.len(), 20);
    for s in &snaps {
        assert_eq!(s.shape(), &[50, 50]);
    }
    let sum_abs: f64 = snaps
        .last()
        .map_or(0.0, |s| s.iter().map(|v| v.abs()).sum());
    assert!(sum_abs > 0.0);
}

#[test]
fn no_steps_gives_no_snapshots() {
    assert!(run_fdtd_simulation(&params(30, 0, (15, 15), 0.1)).is_empty());
}

#[test]
fn field_respects_causality() {
    // One cell per half-step: after 20 steps nothing can be more than 40 cells from the source.
    let p = FdtdParameters {
        width: 120,
        height: 40,
        time_steps: 21,
        source_pos: (10, 20),
        source_freq: 0.1,
    };
    let snaps = run_fdtd_simulation(&p);
    let last = snaps.last().unwrap_or_else(|| panic!("no snapshots"));
    for i in 60..120 {
        for j in 0..40 {
            assert_eq!(last[[i, j]], 0.0, "field at ({i},{j})");
        }
    }
    // ...but the source region itself is excited.
    assert!(last.iter().any(|&v| v != 0.0));
}

#[test]
fn field_is_symmetric_under_transposition_for_a_diagonal_source() {
    let snaps = run_fdtd_simulation(&params(40, 60, (20, 20), 0.05));
    let last = snaps.last().unwrap_or_else(|| panic!("no snapshots"));
    for i in 0..40 {
        for j in 0..40 {
            assert!((last[[i, j]] - last[[j, i]]).abs() < 1e-9, "({i},{j})");
        }
    }
}

#[test]
fn boundary_damping_lets_the_field_decay() {
    let snaps = run_fdtd_simulation(&params(40, 400, (20, 20), 0.1));
    let energy = |a: &ndarray::Array2<f64>| a.iter().map(|v| v * v).sum::<f64>();
    let peak = snaps.iter().map(energy).fold(0.0, f64::max);
    let end = energy(snaps.last().unwrap_or_else(|| panic!("no snapshots")));
    assert!(peak > 0.0);
    assert!(end < peak, "end energy {end} vs peak {peak}");
}

#[cfg(feature = "npy")]
#[test]
fn simulate_and_save_final_state_writes_a_zero_column() {
    // The 1D helper has no source term, so the saved field is identically zero.
    let path = std::env::temp_dir().join(format!("rssn_fdtd_{}.npy", std::process::id()));
    simulate_and_save_final_state(16, 10, path.to_str().unwrap_or_else(|| panic!("path")))
        .unwrap_or_else(|e| panic!("{e}"));
    let arr = rssn::io::read_npy_file(&path).unwrap_or_else(|e| panic!("{e}"));
    let _ = std::fs::remove_file(&path);
    assert_eq!(arr.shape(), &[16, 1]);
    assert!(arr.iter().all(|&v| v == 0.0));
}

#[cfg(not(feature = "npy"))]
#[test]
fn simulate_and_save_final_state_requires_the_npy_feature() {
    let err = simulate_and_save_final_state(16, 10, "unused.npy").err();
    assert!(err.is_some_and(|e| e.contains("npy")));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_fdtd_fields_stay_finite(width in 20usize..40, freq in 0.01f64..0.5) {
        let snaps = run_fdtd_simulation(&params(width, 20, (width / 2, width / 2), freq));
        for s in &snaps {
            prop_assert!(s.iter().all(|v| v.is_finite()));
        }
        prop_assert_eq!(snaps.len(), 4);
    }
}
