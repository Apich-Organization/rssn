//! Unit and property-based tests for the physics sim FDTD electrodynamics module.

use ndarray::Array2;
use proptest::prelude::*;
use rssn::physics::physics_sim::fdtd_electrodynamics::*;

#[test]
fn test_run_fdtd_simulation_basic() {
    let params = FdtdParameters {
        width: 50,
        height: 50,
        time_steps: 100,
        source_pos: (25, 25),
        source_freq: 0.1,
    };

    let snapshots = run_fdtd_simulation(&params);
    // 100 / 5 = 20 snapshots (actually 0, 5, ..., 95, so 20 snapshots)
    assert!(snapshots.len() >= 19);

    let final_ez = snapshots.last().unwrap();
    assert_eq!(final_ez.shape(), &[50, 50]);

    // Check if there is some activity in the field
    let sum_abs: f64 = final_ez.iter().map(|&x| x.abs()).sum();
    assert!(sum_abs > 0.0);
}

#[test]
fn test_simulate_and_save_final_state_smoke() {
    let res = simulate_and_save_final_state(50, 50, "test_fdtd.npy");
    assert!(res.is_ok());
    // Cleanup if needed, but we ignored .npy in gitignore
}

proptest! {
    #[test]
    fn prop_fdtd_stability(width in 20..40usize, freq in 0.01..0.5f64) {
        let params = FdtdParameters {
            width,
            height: width,
            time_steps: 20,
            source_pos: (width/2, width/2),
            source_freq: freq,
        };

        let snapshots = run_fdtd_simulation(&params);
        if let Some(ez) = snapshots.last() {
            for &val in ez.iter() {
                prop_assert!(val.is_finite());
            }
        }
    }
}
