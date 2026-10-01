//! Projection-method Navier-Stokes solvers (ported from `physics_sim_navier_stokes_test.rs`).

use ndarray::Array2;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::navier_stokes_fluid::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 12,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn cavity(
    n: usize,
    dt: f64,
    iters: usize,
    lid: f64,
) -> NavierStokesParameters {
    NavierStokesParameters {
        nx: n,
        ny: n,
        re: 10.0,
        dt,
        n_iter: iters,
        lid_velocity: lid,
    }
}

#[test]
fn cavity_output_shapes_and_lid_row() {
    const N: usize = 17;
    let (u, v, p) =
        run_lid_driven_cavity(&cavity(N, 0.001, 10, 1.0)).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.shape(), &[N, N]);
    assert_eq!(v.shape(), &[N, N]);
    assert_eq!(p.shape(), &[N, N]);
    // The lid row is never touched by the pressure projection.
    for i in 0..N {
        assert!((u[[N - 1, i]] - 1.0).abs() < 1e-12);
    }
    // The bottom wall row stays at rest.
    for i in 0..N {
        assert_eq!(u[[0, i]], 0.0);
    }
}

#[test]
fn zero_lid_velocity_gives_a_quiescent_fluid() {
    let (u, v, p) =
        run_lid_driven_cavity(&cavity(9, 0.01, 5, 0.0)).unwrap_or_else(|e| panic!("{e}"));
    for a in [&u, &v, &p] {
        assert!(a.iter().all(|&x| x == 0.0));
    }
}

#[test]
fn moving_lid_drags_the_fluid_beneath_it() {
    let (u, _, _) =
        run_lid_driven_cavity(&cavity(17, 0.01, 30, 1.0)).unwrap_or_else(|e| panic!("{e}"));
    let below: f64 = (1..16).map(|i| u[[15, i]].abs()).sum();
    assert!(below > 0.0, "no flow induced below the lid");
}

#[test]
fn channel_flow_respects_inflow_wall_and_obstacle_conditions() {
    let n = 17;
    let mut mask = Array2::<bool>::from_elem((n, n), false);
    for j in 6..11 {
        for i in 6..9 {
            mask[[j, i]] = true;
        }
    }
    let (u, v, p) =
        run_channel_flow(n, n, 100.0, 0.001, 20, &mask).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(u.shape(), &[n, n]);
    assert_eq!(p.shape(), &[n, n]);
    for j in 0..n {
        assert_eq!(u[[j, 0]], 1.0, "inflow");
        assert_eq!(v[[j, 0]], 0.0);
    }
    for i in 0..n {
        assert_eq!(v[[0, i]], 0.0, "wall");
        assert_eq!(v[[n - 1, i]], 0.0, "wall");
    }
    for j in 6..11 {
        for i in 6..9 {
            assert_eq!((u[[j, i]], v[[j, i]]), (0.0, 0.0), "obstacle at ({j},{i})");
        }
    }
    for a in [&u, &v, &p] {
        assert!(a.iter().all(|x| x.is_finite()));
    }
}

#[test]
fn unobstructed_channel_velocity_stays_between_zero_and_the_inflow_speed() {
    let n = 17;
    let mask = Array2::<bool>::from_elem((n, n), false);
    let (u, _, _) =
        run_channel_flow(n, n, 100.0, 0.001, 300, &mask).unwrap_or_else(|e| panic!("{e}"));
    for i in 0..n {
        assert!(
            (-1e-9..=1.0 + 1e-9).contains(&u[[8, i]]),
            "u[8,{i}] = {}",
            u[[8, i]]
        );
    }
    // The flow has developed along the channel axis.
    assert!(u[[8, n - 1]] > 0.2);
}

#[test]
fn lid_cavity_scenario_writes_into_the_given_directory() {
    let dir = std::env::temp_dir().join(format!("rssn_cavity_{}", std::process::id()));
    simulate_lid_driven_cavity_scenario(&dir);
    if cfg!(feature = "npy") {
        assert!(dir.join("cavity_u_velocity.npy").is_file());
        assert!(dir.join("cavity_pressure.npy").is_file());
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn channel_flow_reports_a_wrongly_sized_mask_instead_of_panicking() {
    let n = 17;
    let small = Array2::<bool>::from_elem((9, 9), false);
    assert!(run_channel_flow(n, n, 100.0, 0.001, 1, &small).is_err());
    let wide = Array2::<bool>::from_elem((n, n + 4), false);
    assert!(run_channel_flow(n, n, 100.0, 0.001, 1, &wide).is_err());
    let ok = Array2::<bool>::from_elem((n, n), false);
    assert!(
        run_channel_flow(n, n + 1, 100.0, 0.001, 1, &ok).is_err(),
        "non-square grid"
    );
    assert!(
        run_channel_flow(
            2,
            2,
            100.0,
            0.001,
            1,
            &Array2::<bool>::from_elem((2, 2), false)
        )
        .is_err()
    );
}

#[test]
fn viscosity_controls_how_fast_the_lid_shear_penetrates() {
    // After one step u just below the lid is ~ dt * nu * U / h^2 (pure diffusion, no advection yet),
    // so it scales like 1 / re and is independent of the (unused before) pressure.
    let n = 17;
    let mut p10 = cavity(n, 0.001, 1, 1.0);
    p10.re = 10.0;
    let mut p100 = p10.clone();
    p100.re = 100.0;
    let (u10, _, _) = run_lid_driven_cavity(&p10).unwrap_or_else(|e| panic!("{e}"));
    let (u100, _, _) = run_lid_driven_cavity(&p100).unwrap_or_else(|e| panic!("{e}"));
    let (a, b) = (u10[[n - 2, 8]], u100[[n - 2, 8]]);
    assert!(a > 0.0 && b > 0.0);
    assert!((a / b - 10.0).abs() < 1.0, "ratio {}", a / b);
    let h = 1.0 / (n as f64 - 1.0);
    assert!((a - 0.001 * 0.1 / (h * h)).abs() < 0.1 * a, "u = {a}");
}

#[test]
fn projection_leaves_the_interior_divergence_free() {
    // div(u) over interior cells must be much smaller than the divergence the explicit step creates.
    let n = 17;
    let (u, v, _) =
        run_lid_driven_cavity(&cavity(n, 0.002, 40, 1.0)).unwrap_or_else(|e| panic!("{e}"));
    let h = 1.0 / (n as f64 - 1.0);
    // Cell-centred fields: central differences approximate the staggered divergence.
    let mut max_div = 0.0f64;
    let mut max_grad = 0.0f64;
    for j in 3..n - 3 {
        for i in 3..n - 3 {
            let div = (u[[j, i + 1]] - u[[j, i - 1]]) / (2.0 * h)
                + (v[[j + 1, i]] - v[[j - 1, i]]) / (2.0 * h);
            max_div = max_div.max(div.abs());
            max_grad = max_grad.max(((u[[j, i + 1]] - u[[j, i - 1]]) / (2.0 * h)).abs());
        }
    }
    assert!(max_grad > 1e-3, "flow too weak to test: {max_grad}");
    assert!(
        max_div < 0.5 * max_grad,
        "div {max_div} vs |du/dx| {max_grad}"
    );
}

#[test]
fn fluid_beneath_the_lid_moves_with_the_lid_and_recirculates() {
    let n = 17;
    let (u, _, _) =
        run_lid_driven_cavity(&cavity(n, 0.002, 400, 1.0)).unwrap_or_else(|e| panic!("{e}"));
    assert!(
        u[[n - 2, n / 2]] > 0.1,
        "near-lid u = {}",
        u[[n - 2, n / 2]]
    );
    // Return flow in the lower half of the cavity.
    let below: f64 = (2..n / 2).map(|j| u[[j, n / 2]]).fold(f64::MAX, f64::min);
    assert!(below < 0.0, "no recirculation: min u = {below}");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_cavity_fields_are_finite(lid in 0.1f64..2.0, dt in 0.001f64..0.01) {
        let n = 9;
        let params = NavierStokesParameters { nx: n, ny: n, re: 100.0, dt, n_iter: 5, lid_velocity: lid };
        let (u, v, p) = run_lid_driven_cavity(&params).unwrap_or_else(|e| panic!("{e}"));
        for a in [&u, &v, &p] {
            prop_assert!(a.iter().all(|x| x.is_finite()));
        }
    }
}
