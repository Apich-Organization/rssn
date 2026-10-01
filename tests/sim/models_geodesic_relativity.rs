//! Schwarzschild geodesics (ported from `physics_sim_geodesic_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::geodesic_relativity::*;

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
