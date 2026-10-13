//! Kerr-Newman geodesics: characteristic radii, circular orbits, capture, convergence.

use rssn::sim::relativity::*;

fn bardeen_isco(a: f64, prograde: bool) -> f64 {
    let z1 = 1.0 + (1.0 - a * a).cbrt() * ((1.0 + a).cbrt() + (1.0 - a).cbrt());
    let z2 = (3.0 * a * a + z1 * z1).sqrt();
    let s = (3.0 - z1) * (3.0 + z1 + 2.0 * z2);
    if prograde { 3.0 + z2 - s.sqrt() } else { 3.0 + z2 + s.sqrt() }
}

#[test]
fn schwarzschild_isco_and_photon_sphere() {
    let bh = KerrNewman::schwarzschild(1.0);
    assert!((bh.isco(true).unwrap() - 6.0).abs() < 1e-6);
    assert!((bh.photon_orbit(true).unwrap() - 3.0).abs() < 1e-6);
    assert!((bh.horizon().unwrap() - 2.0).abs() < 1e-14);
    let (e, l, om) = bh.circular_orbit(6.0, true).unwrap();
    assert!((e - (8.0_f64 / 9.0).sqrt()).abs() < 1e-12);
    assert!((l - 12.0_f64.sqrt()).abs() < 1e-12);
    assert!((om - 216.0_f64.powf(-0.5)).abs() < 1e-12);
}

#[test]
fn kerr_isco_and_photon_orbits_match_closed_forms() {
    for a in [0.3, 0.6, 0.9] {
        let bh = KerrNewman { m: 1.0, a, q: 0.0 };
        for prograde in [true, false] {
            let isco = bh.isco(prograde).unwrap();
            assert!((isco - bardeen_isco(a, prograde)).abs() < 1e-5, "a={a} prograde={prograde}: {isco}");
        }
        let ph = bh.photon_orbit(true).unwrap();
        let exact = 2.0 * (1.0 + (2.0 / 3.0 * (-a).acos()).cos());
        assert!((ph - exact).abs() < 1e-6, "photon orbit {ph} vs {exact}");
        let phr = bh.photon_orbit(false).unwrap();
        let exact_r = 2.0 * (1.0 + (2.0 / 3.0 * a.acos()).cos());
        assert!((phr - exact_r).abs() < 1e-6, "retro photon orbit {phr} vs {exact_r}");
    }
}

#[test]
fn reissner_nordstrom_photon_sphere_and_isco_shift() {
    for q in [0.3, 0.6, 0.9] {
        let bh = KerrNewman { m: 1.0, a: 0.0, q };
        let exact = (3.0 + (9.0 - 8.0 * q * q).sqrt()) / 2.0;
        assert!((bh.photon_orbit(true).unwrap() - exact).abs() < 1e-6);
        let isco = bh.isco(true).unwrap();
        assert!(isco < 6.0 && isco > 4.0, "q={q}: ISCO {isco}");
    }
    // charge-to-mass ratio above one: naked singularity, no horizon
    assert!(KerrNewman { m: 1.0, a: 0.0, q: 1.2 }.horizon().is_none());
}

#[test]
fn circular_orbit_stays_circular_and_has_keplerian_frequency() {
    let bh = KerrNewman::schwarzschild(1.0);
    let mut g = Geodesic::circular(bh, 10.0, true).unwrap();
    let omega = 10.0_f64.powf(-1.5);
    for _ in 0..2000 {
        g.step(0.5, true).unwrap();
    }
    assert!((g.r - 10.0).abs() < 1e-5, "radius drift {}", g.r - 10.0);
    assert!((g.phi / g.t - omega).abs() < 1e-6 * omega, "angular frequency {}", g.phi / g.t);
    assert!(g.constraint().abs() < 1e-12);
}

#[test]
fn kerr_circular_orbit_and_constraint_preservation() {
    let bh = KerrNewman { m: 1.0, a: 0.9, q: 0.0 };
    let (_, _, om) = bh.circular_orbit(8.0, true).unwrap();
    assert!((om - 1.0 / (8.0_f64.powf(1.5) + 0.9)).abs() < 1e-12);
    let mut g = Geodesic::circular(bh, 8.0, true).unwrap();
    for _ in 0..1000 {
        g.step(0.5, true).unwrap();
    }
    assert!((g.r - 8.0).abs() < 1e-5);
    assert!((g.phi / g.t - om).abs() < 1e-4 * om);
}

fn eccentric_r(
    h: f64,
    project: bool,
) -> f64 {
    let bh = KerrNewman::schwarzschild(1.0);
    // bound orbit with E < 1, started at apocentre-like point with some radial velocity
    let mut g = Geodesic::new(bh, 20.0, std::f64::consts::FRAC_PI_2, 0.0, 0.97, 4.2, 0.0, 1.0, 1.0).unwrap();
    let n = (400.0 / h).round() as usize;
    g.integrate(h, n, project, 2.5, 1e6).unwrap();
    g.r
}

#[test]
fn implicit_midpoint_is_second_order_on_a_precessing_orbit() {
    let reference = eccentric_r(0.25, true);
    let e1 = (eccentric_r(4.0, true) - reference).abs();
    let e2 = (eccentric_r(2.0, true) - reference).abs();
    let order = (e1 / e2).log2();
    assert!((order - 2.0).abs() < 0.4, "order {order} ({e1:e}, {e2:e})");
}

#[test]
fn constraint_stays_on_the_mass_shell_with_projection() {
    let bh = KerrNewman { m: 1.0, a: 0.5, q: 0.2 };
    let mut g = Geodesic::new(bh, 15.0, 1.2, 0.0, 0.98, 3.5, 0.4, -1.0, 1.0).unwrap();
    assert!(g.constraint().abs() < 1e-12);
    for _ in 0..500 {
        g.step(0.1, true).unwrap();
        assert!(g.constraint().abs() < 1e-11, "constraint {}", g.constraint());
    }
    // E and L are exactly conserved (cyclic coordinates)
    assert_eq!(g.energy, 0.98);
    assert_eq!(g.lz, 3.5);
}

#[test]
fn radial_null_infall_reports_horizon_crossing_with_exact_affine_length() {
    let bh = KerrNewman::schwarzschild(1.0);
    let mut g = Geodesic::new(bh, 10.0, std::f64::consts::FRAC_PI_2, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0).unwrap();
    let r_stop = 2.01;
    let end = g.integrate(0.002, 20000, true, r_stop, 1e6).unwrap();
    assert_eq!(end, Termination::Horizon);
    assert!((g.r - r_stop).abs() < 1e-9);
    // dr/dlambda = -E exactly for radial null rays
    assert!((g.lambda - (10.0 - r_stop)).abs() < 1e-3, "lambda {}", g.lambda);
    assert!(g.t > 15.0, "coordinate time diverges toward the horizon: {}", g.t);
}

#[test]
fn photon_capture_threshold_is_the_critical_impact_parameter() {
    let bh = KerrNewman::schwarzschild(1.0);
    let bc = 27.0_f64.sqrt();
    let run = |b: f64| {
        let mut g = Geodesic::new(bh, 30.0, std::f64::consts::FRAC_PI_2, 0.0, 1.0, b, 0.0, -1.0, 0.0).unwrap();
        g.integrate(0.05, 40_000, true, 2.05, 31.0).unwrap()
    };
    assert_eq!(run(bc * 0.97), Termination::Horizon);
    assert_eq!(run(bc * 1.03), Termination::Escaped);
}

#[test]
fn forbidden_initial_data_is_rejected() {
    let bh = KerrNewman::schwarzschild(1.0);
    // L far too large for a bound start at r = 3
    assert!(Geodesic::new(bh, 3.0, std::f64::consts::FRAC_PI_2, 0.0, 0.9, 6.0, 0.0, 1.0, 1.0).is_err());
    assert!(Geodesic::circular(bh, 2.5, true).is_err());
}
