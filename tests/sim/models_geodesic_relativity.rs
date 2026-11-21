//! Schwarzschild geodesics (ported from `physics_sim_geodesic_test.rs`).

use std::num::NonZeroUsize;

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::geodesic_relativity::*;
use rssn::sim::relativity::KerrNewman;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 12,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn params(
    mass: f64,
    state: [f64; 4],
    end: f64,
) -> GeodesicParameters {
    GeodesicParameters {
        black_hole_mass: mass,
        initial_state: state,
        proper_time_end: end,
        initial_dt: 0.1,
    }
}

#[test]
fn effective_potential_matches_the_hand_computation() {
    // L = r^2 phi_dot = 3.5:  V = -1/10 + L^2 / 200 - L^2 / 1000 = -0.051.
    let p = params(1.0, [10.0, 0.0, 0.0, 0.035], 10.0);
    assert!((p.effective_potential(10.0, 3.5) + 0.051).abs() < 1e-12);
}

#[test]
fn effective_potential_without_angular_momentum_is_newtonian() {
    for (m, r) in [(1.0, 5.0), (2.0, 7.0), (0.5, 100.0)] {
        let p = params(m, [r, 0.0, 0.0, 0.0], 1.0);
        assert!((p.effective_potential(r, 0.0) + m / r).abs() < 1e-15);
    }
}

#[test]
fn schwarzschild_isco_is_an_inflection_of_the_effective_potential() {
    // At r = 6M the circular-orbit angular momentum is L^2 = 12 M^2, and V'(r) = V''(r) = 0.
    let p = params(1.0, [6.0, 0.0, 0.0, 0.0], 1.0);
    let l = 12f64.sqrt();
    let h = 1e-4;
    let d1 = (p.effective_potential(6.0 + h, l) - p.effective_potential(6.0 - h, l)) / (2.0 * h);
    let d2 = (p.effective_potential(6.0 + h, l) - 2.0 * p.effective_potential(6.0, l)
        + p.effective_potential(6.0 - h, l))
        / (h * h);
    assert!(d1.abs() < 1e-8, "V' = {d1}");
    assert!(d2.abs() < 1e-5, "V'' = {d2}");
}

#[test]
fn path_starts_at_the_initial_cartesian_position() {
    let path = run_geodesic_simulation(&params(1.0, [10.0, 0.0, 0.0, 0.035], 100.0));
    assert!(path.len() > 10);
    let (x0, y0) = path[0];
    assert!((x0 - 10.0).abs() < 1e-12);
    assert!(y0.abs() < 1e-12);
}

#[test]
fn circular_orbit_keeps_a_constant_radius() {
    // Circular geodesic at r = 10 M: phi_dot = sqrt(M / (r^2 (r - 3M))).
    let w = (1.0f64 / (100.0 * 7.0)).sqrt();
    let path = run_geodesic_simulation(&params(1.0, [10.0, 0.0, 0.0, w], 300.0));
    for (x, y) in &path {
        let r = x.hypot(*y);
        assert!((r - 10.0).abs() < 1e-3, "r = {r}");
    }
    // It really goes round: 300 * 0.0378 = 11.3 rad, so it reaches y < -9.
    assert!(path.iter().any(|&(_, y)| y < -9.0));
}

#[test]
fn radial_infall_stays_on_the_x_axis() {
    let path = run_geodesic_simulation(&params(1.0, [20.0, 0.0, 0.0, 0.0], 50.0));
    for (_, y) in &path {
        assert!(y.abs() < 1e-12);
    }
    let (x_end, _) = path[path.len() - 1];
    assert!(x_end < 20.0);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_geodesic_path_is_finite(mass in 0.5f64..2.0, r0 in 6.0f64..20.0) {
        let path = run_geodesic_simulation(&params(mass, [r0, 0.0, 0.0, 0.01], 50.0));
        prop_assert!(!path.is_empty());
        for (x, y) in path {
            prop_assert!(x.is_finite() && y.is_finite());
        }
    }

    #[test]
    fn prop_potential_scales_linearly_with_mass_at_fixed_l_zero(m in 0.1f64..10.0, r in 3.0f64..50.0) {
        let p = params(m, [r, 0.0, 0.0, 0.0], 1.0);
        prop_assert!((p.effective_potential(r, 0.0) + m / r).abs() < 1e-12);
    }
}

#[test]
fn black_hole_scenario_writes_one_csv_per_orbit() {
    let dir = std::env::temp_dir().join(format!("rssn_geodesic_{}", std::process::id()));
    simulate_black_hole_orbits_scenario(&dir).unwrap_or_else(|e| panic!("{e}"));
    for name in ["stable_orbit", "plunging_orbit", "photon_orbit"] {
        let file = dir.join(format!("orbit_{name}.csv"));
        let text = std::fs::read_to_string(&file).unwrap_or_else(|e| panic!("{file:?}: {e}"));
        assert!(text.starts_with("x,y\n"), "{name}");
        assert!(text.lines().count() > 2, "{name} has no path points");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn schwarzschild_engine_path_agrees_with_the_adaptive_rk45_reference() {
    let p = params(1.0, [12.0, 0.0, 0.0, 0.03], 150.0);
    let a = run_geodesic_simulation(&p);
    let b = run_geodesic_simulation_rk45(&p);
    let (xa, ya) = a[a.len() - 1];
    let (xb, yb) = b[b.len() - 1];
    assert!(
        ((xa - xb).powi(2) + (ya - yb).powi(2)).sqrt() < 0.2,
        "({xa},{ya}) vs ({xb},{yb})"
    );
}

#[test]
fn kerr_frame_dragging_twists_a_zero_angular_velocity_infall() {
    let p = params(1.0, [15.0, 0.0, 0.0, 0.0], 60.0);
    let kerr = run_kerr_newman_geodesic(&p, 0.9, 0.0).unwrap();
    let schw = run_kerr_newman_geodesic(&p, 0.0, 0.0).unwrap();
    let rn = run_kerr_newman_geodesic(&p, 0.0, 0.6).unwrap();
    assert!(schw.iter().all(|&(_, y)| y.abs() < 1e-12));
    assert!(rn.iter().all(|&(_, y)| y.abs() < 1e-12));
    // frame dragging: the azimuth grows in the sense of the spin although the initial L is not zero
    let (xk, yk) = kerr[kerr.len() - 1];
    assert!(yk.abs() > 1e-4 && xk.hypot(yk) < 15.0, "({xk},{yk})");
    // charge changes the radial motion
    let (xs, _) = schw[schw.len() - 1];
    let (xr, _) = rn[rn.len() - 1];
    assert!((xs - xr).abs() > 1e-3, "{xs} vs {xr}");
}

#[test]
fn kerr_circular_orbit_keeps_its_radius() {
    let (m, a, r) = (1.0, 0.9, 8.0);
    let k = KerrNewman { m, a, q: 0.0 };
    let (_, _, omega) = k.circular_orbit(r, true).unwrap();
    let w = 2.0 * m / r;
    let (gtt, gtp, gpp) = (-(1.0 - w), -a * w, r * r + a * a + a * a * w);
    let ut = 1.0 / (-(gtt + 2.0 * gtp * omega + gpp * omega * omega)).sqrt();
    let path =
        run_kerr_newman_geodesic(&params(m, [r, 0.0, 0.0, omega * ut], 200.0), a, 0.0).unwrap();
    for (x, y) in &path {
        assert!((x.hypot(*y) - r).abs() < 1e-3, "r = {}", x.hypot(*y));
    }
    assert!(path.iter().any(|&(_, y)| y < -r * 0.9));
    // a naked-singularity-style initial radius inside the horizon is rejected by the engine
    assert!(run_kerr_newman_geodesic(&params(1.0, [1.5, 0.0, 0.0, 0.01], 1.0), 0.0, 0.0).is_err());
}

#[test]
fn analytic_prograde_circular_orbits_are_valid_across_the_stationary_limit() {
    // Specialize Bardeen–Press–Teukolsky (1972), equations 2.12–2.13, to
    // unit mass and zero charge. The proper angular rate follows from
    // Ω = 1 / (r^(3/2) + a) and the circular-orbit normalization.
    let a = 0.99;
    let metric = KerrNewman { m: 1.0, a, q: 0.0 };
    for r in [1.9, 2.0_f64.next_down(), 2.0, 2.0_f64.next_up(), 2.1] {
        let root_r = r.sqrt();
        let r32 = r * root_r;
        let den = root_r * (root_r * (r32 - 3.0 * root_r + 2.0 * a)).sqrt();
        let phid = den.recip();
        let energy = (r32 - 2.0 * root_r + a) / den;
        let angular = (r * r - 2.0 * a * root_r + a * a) / den;
        for sign in [1.0, -1.0] {
            let g = equatorial_timelike_geodesic(
                KerrNewman {
                    a: sign * a,
                    ..metric
                },
                [r, 0.0, 0.125, sign * phid],
                EquatorialEnergyBranch::Nonnegative,
            )
            .unwrap();
            assert!(
                (g.energy - energy).abs() < 2e-14,
                "{r}: {} vs {energy}",
                g.energy
            );
            assert!((g.lz - sign * angular).abs() < 4e-14, "{r}: {}", g.lz);
            assert!(g.constraint().abs() < 4e-14, "{r}: {}", g.constraint());
        }
        if r < 1.95 {
            let mut p = params(1.0, [r, 0.0, 0.125, phid], 1.0);
            p.initial_dt = 0.005;
            let path = run_kerr_newman_geodesic(&p, a, 0.0).unwrap();
            assert!(path.iter().all(|&(x, y)| (x.hypot(y) - r).abs() < 2e-7));
            let &(x, y) = path.last().unwrap();
            let phase = 0.125 + phid;
            let error = (x - r * phase.cos()).hypot(y - r * phase.sin());
            assert!(error < 2e-7, "{error}");
        }
    }
}

#[test]
fn both_ergoregion_energy_branches_are_future_directed() {
    let metric = KerrNewman {
        m: 1.0,
        a: 1.0,
        q: 0.0,
    };
    // At an equatorial radius 3/2 the metric coefficients are rational:
    // g_tt=1/3, g_tphi=-4/3, g_phiphi=55/12 and g_rr=9.
    let r_dot = 0.1;
    let ep = (1.0_f64 - (1.0 + 9.0 * r_dot * r_dot) / 3.0).sqrt();
    for (branch, energy) in [
        (EquatorialEnergyBranch::Nonnegative, ep),
        (EquatorialEnergyBranch::Nonpositive, -ep),
    ] {
        let state = [1.5, r_dot, 0.0, 2.0];
        let geodesic = equatorial_timelike_geodesic(metric, state, branch).unwrap();
        let ut = 8.0 - 3.0 * energy;
        let want_lz = -4.0 / 3.0 * ut + 55.0 / 6.0;
        assert!((geodesic.energy - energy).abs() < 3e-15);
        assert!((geodesic.lz - want_lz).abs() < 5e-15);
        assert!((geodesic.pr - (9.0 * r_dot)).abs() < 1e-15);
        assert!(geodesic.constraint().abs() < 1e-12);
        // Hamiltonian evolution must advance coordinate time on either root.
        let mut advanced = geodesic;
        advanced.step(0.0001, true).unwrap();
        assert!((advanced.t / 0.0001 - ut).abs() < 0.01);

        let dt = 1.0 / 1024.0;
        let mut p = params(1.0, state, 2.0 * dt);
        p.initial_dt = dt;
        let sampled = run_kerr_newman_geodesic_sampled(
            &p,
            1.0,
            0.0,
            branch,
            NonZeroUsize::new(usize::MAX).unwrap(),
        )
        .unwrap();
        let mut terminal = geodesic;
        terminal.integrate(dt, 2, true, 1.001, 1e4).unwrap();
        assert_eq!(
            sampled,
            [
                (1.5, 0.0),
                (
                    terminal.r * terminal.phi.cos(),
                    terminal.r * terminal.phi.sin()
                )
            ]
        );
    }

    for state in [
        [0.5, 0.0, 0.0, 2.0],
        [1.0, 0.0, 0.0, 2.0],
        [1.5, 0.0, 0.0, -2.0],
        [1.5, 0.0, 0.0, 0.125],
        [1.5, 0.0, f64::NAN, 2.0],
    ] {
        assert!(
            equatorial_timelike_geodesic(metric, state, EquatorialEnergyBranch::Nonnegative)
                .is_err()
        );
    }
    for r in [2.0, 2.1] {
        assert!(
            equatorial_timelike_geodesic(
                metric,
                [r, 0.0, 0.0, 2.0],
                EquatorialEnergyBranch::Nonpositive
            )
            .is_err()
        );
    }
}

#[test]
fn sampled_trajectory_retains_only_stride_points_and_the_reached_end() {
    let branch = EquatorialEnergyBranch::Nonnegative;
    let mut plunge = params(1.0, [2.01, -0.1, 0.0, 0.0], 0.3);
    plunge.initial_dt = 0.001;
    let cases = [
        params(1.0, [10.0, 0.0, 0.25, 0.05], 0.0),
        params(1.0, [10.0, 0.0, 0.25, 0.05], 3.17),
        params(1.0, [9_999.5, 8.0, 0.0, 0.0], 2.0),
        plunge,
    ];
    for p in cases {
        let full = run_kerr_newman_geodesic(&p, 0.0, 0.0).unwrap();
        for stride in [1, 3, 7, usize::MAX] {
            let sampled = run_kerr_newman_geodesic_sampled(
                &p,
                0.0,
                0.0,
                branch,
                NonZeroUsize::new(stride).unwrap(),
            )
            .unwrap();
            let mut expected: Vec<_> = full.iter().step_by(stride).copied().collect();
            if (full.len() - 1) % stride != 0 {
                expected.push(*full.last().unwrap());
            }
            assert_eq!(sampled, expected);
        }
        if p.initial_state[0] > 9_000.0 {
            let &(x, y) = full.last().unwrap();
            assert!(x.hypot(y) >= 1e4);
            assert!(full.len() < 5);
        } else if p.initial_state[0] < 3.0 {
            let &(x, y) = full.last().unwrap();
            assert!((x.hypot(y) - 2.002).abs() < 1e-12);
            assert!(full.len() < 300);
        }
    }
}

#[test]
fn invalid_time_bounds_terminate_without_invoking_the_legacy_integrator() {
    for (end, dt) in [
        (f64::INFINITY, 0.1),
        (f64::NAN, 0.1),
        (-1.0, 0.1),
        (1.0, f64::INFINITY),
        (1.0, f64::NAN),
        (1.0, 0.0),
        (1.0, -0.1),
    ] {
        let mut p = params(1.0, [10.0, 0.0, 0.0, 0.045], end);
        p.initial_dt = dt;
        assert!(run_kerr_newman_geodesic(&p, 0.0, 0.0).is_err());
        assert_eq!(run_geodesic_simulation(&p), [(10.0, 0.0)]);
        assert_eq!(run_geodesic_simulation_rk45(&p), [(10.0, 0.0)]);
    }
}
