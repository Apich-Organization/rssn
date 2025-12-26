//! Unit and property-based tests for the physics sim Navier-Stokes module.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::physics::physics_sim::navier_stokes_fluid::*;

#[test]

fn test_simulate_lid_driven_cavity_scenario(
) {

    simulate_lid_driven_cavity_scenario(
    );
}

#[test]

fn test_lid_driven_cavity_basic() {

    const K : usize = 4;

    const N : usize =
        2_usize.pow(K as u32) + 1; // 17
    let params =
        NavierStokesParameters {
            nx : N,
            ny : N,
            re : 10.0,
            dt : 0.001,
            n_iter : 10,
            lid_velocity : 1.0,
        };

    let res =
        run_lid_driven_cavity(&params)
            .unwrap();

    let (u, v, p) = res;

    assert_eq!(u.shape(), &[N, N]);

    assert_eq!(v.shape(), &[N, N]);

    assert_eq!(p.shape(), &[N, N]);

    // Check if lid velocity is reflected in top row of U
    // Top row indices in centered U are [N-1, :]
    // Since lid velocity is 1.0, we expect some positive value in the top row
    let top_row_avg =
        u.row(N - 1).sum() / N as f64;

    assert!(top_row_avg > 0.0);
}

proptest! {
    #[test]
    fn prop_navier_stokes_stability(lid_v in 0.1..2.0f64, dt in 0.001..0.01f64) {
        const N: usize = 9; // small grid for speed
        let params = NavierStokesParameters {
            nx: N,
            ny: N,
            re: 100.0,
            dt,
            n_iter: 5,
            lid_velocity: lid_v,
        };

        let res = run_lid_driven_cavity(&params).unwrap();
        let (u, v, p) = res;

        for val in u.iter() { prop_assert!(val.is_finite()); }
        for val in v.iter() { prop_assert!(val.is_finite()); }
        for val in p.iter() { prop_assert!(val.is_finite()); }
    }
}
