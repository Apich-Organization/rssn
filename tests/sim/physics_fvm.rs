//! Finite volume method solvers (ported from `physics_fvm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_fvm::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn top_hat(x: f64) -> f64 {
    if x > 0.4 && x < 0.6 {
        1.0
    } else {
        0.0
    }
}

#[test]
fn mesh_places_cell_centres_and_initial_values() {
    let mesh = Mesh::new(4, 2.0, |x| x);
    assert_eq!(mesh.num_cells(), 4);
    assert_approx_eq!(mesh.dx, 0.5);
    let expect = [0.25, 0.75, 1.25, 1.75];
    for (c, e) in mesh.cells.iter().zip(expect) {
        assert_approx_eq!(c.value, e);
    }
}

#[test]
fn limiters_have_known_values() {
    assert_eq!(minmod(1.0, 2.0), 1.0);
    assert_eq!(minmod(-3.0, -2.0), -2.0);
    assert_eq!(minmod(1.0, -2.0), 0.0);
    assert_eq!(minmod(0.0, 5.0), 0.0);
    assert_approx_eq!(van_leer(1.0, 1.0), 1.0);
    assert_approx_eq!(van_leer(1.0, 3.0), 1.5);
    assert_eq!(van_leer(-1.0, 2.0), 0.0);
}

#[test]
fn lax_friedrichs_flux_of_a_uniform_state_is_the_physical_flux() {
    let f = |u: f64| 0.5 * u * u;
    assert_approx_eq!(lax_friedrichs_flux(2.0, 2.0, 0.1, 0.1, f), 2.0);
    // Jump: 0.5 (f(1) + f(3)) - 0.5 (dx/dt)(3 - 1) = 2.5 - 1.0 with dx/dt = 1.
    assert_approx_eq!(lax_friedrichs_flux(1.0, 3.0, 0.1, 0.1, f), 1.5);
}

#[test]
fn advection_1d_with_unit_cfl_shifts_exactly() {
    let n = 100;
    let mut mesh = Mesh::new(n, 1.0, top_hat);
    let before: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();
    let dx = 1.0 / n as f64;
    let res = solve_advection_1d(&mut mesh, 1.0, dx, 10, || (0.0, 0.0));
    assert_eq!(res.len(), n);
    for i in 10..n {
        assert_approx_eq!(res[i], before[i - 10], 1e-12);
    }
}

#[test]
fn advection_1d_stays_bounded_at_half_cfl() {
    let mut mesh = Mesh::new(100, 1.0, top_hat);
    let dt = 0.5 * 0.01;
    let res = solve_advection_1d(&mut mesh, 1.0, dt, 10, || (0.0, 0.0));
    for &v in &res {
        assert!((-1e-12..=1.0 + 1e-12).contains(&v));
    }
    // Mass moves right but is conserved while nothing has left the domain.
    let mass: f64 = res.iter().sum();
    assert_approx_eq!(mass, 20.0, 1e-9);
}

#[test]
fn advection_1d_scenario_is_bounded_and_moved_right() {
    let res = simulate_1d_advection_scenario();
    assert_eq!(res.len(), 200);
    assert!(res.iter().all(|&v| (-1e-9..=1.0 + 1e-9).contains(&v)));
    // Pulse started on (0.2, 0.4); after t = 0.5 it is centred near 0.8.
    let (imax, _) = res
        .iter()
        .enumerate()
        .fold((0, f64::MIN), |a, (i, &v)| if v > a.1 { (i, v) } else { a });
    let x = (imax as f64 + 0.5) / 200.0;
    assert!((x - 0.8).abs() < 0.08, "peak at {x}");
}

#[test]
fn burgers_uniform_state_is_stationary() {
    let mut mesh = Mesh::new(50, 1.0, |_| 0.7);
    let res = solve_burgers_1d(&mut mesh, 0.001, 30);
    for &v in &res {
        assert_approx_eq!(v, 0.7, 1e-12);
    }
}

#[test]
fn burgers_shock_moves_at_the_rankine_hugoniot_speed() {
    // u_L = 1, u_R = 0  =>  shock speed 1/2.
    let mut mesh = Mesh::new(200, 1.0, |x| if x < 0.25 { 1.0 } else { 0.0 });
    let dt = 0.001;
    let res = solve_burgers_1d(&mut mesh, dt, 400); // t = 0.4, shock near 0.45
    let total: f64 = res.iter().sum::<f64>() / 200.0;
    // Integral grows by the inflow flux u_L^2 / 2 = 0.5 over t = 0.4.
    assert_approx_eq!(total, 0.45, 1e-3);
    assert!(res[60] > 0.9, "behind the shock: {}", res[60]);
    assert!(res[120] < 0.1, "ahead of the shock: {}", res[120]);
    let mid = res.iter().position(|&v| v < 0.5).unwrap_or(0);
    let x = mid as f64 / 200.0;
    assert!((x - 0.45).abs() < 0.03, "shock at {x}");
}

#[test]
fn shallow_water_still_water_stays_still() {
    let n = 40;
    let res = solve_shallow_water_1d(vec![1.0; n], vec![0.0; n], 0.025, 0.001, 50, 9.81);
    for s in &res {
        assert_approx_eq!(s.h, 1.0, 1e-12);
        assert_approx_eq!(s.hu, 0.0, 1e-12);
    }
}

#[test]
fn shallow_water_dam_break_conserves_mass_and_forms_a_bore() {
    let n = 100;
    let mut h = vec![1.0; n];
    for v in h.iter_mut().skip(50) {
        *v = 0.5;
    }
    let dx = 1.0 / n as f64;
    let mass0: f64 = h.iter().sum::<f64>() * dx;
    let res = solve_shallow_water_1d(h, vec![0.0; n], dx, 0.001, 50, 9.81);
    let mass: f64 = res.iter().map(|s| s.h).sum::<f64>() * dx;
    assert_approx_eq!(mass, mass0, 1e-9);
    assert!(res[45].h < 1.0); // rarefaction
    assert!(res[55].h > 0.5); // bore
    assert!(res[50].hu > 0.0); // flow toward the shallow side
    assert!(res.iter().all(|s| s.h > 0.0 && s.h.is_finite()));
}

fn border2(
    i: usize,
    j: usize,
    w: usize,
    h: usize,
) -> bool {
    i == 0 || j == 0 || i == w - 1 || j == h - 1
}

#[test]
fn advection_2d_unit_cfl_shifts_exactly_in_x() {
    let (w, h) = (30, 20);
    let mut mesh = Mesh2D::new(w, h, (1.0, 1.0), |x, y| (3.0 * x + y).sin() + 2.0);
    let before: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();
    let dt = mesh.dx;
    let res = solve_advection_2d(&mut mesh, (1.0, 0.0), dt, 3, border2);
    assert_eq!(res.len(), w * h);
    // Interior cells more than 3 steps from the left wall hold the shifted profile.
    for j in 1..h - 1 {
        for i in 5..w - 1 {
            assert_approx_eq!(res[j * w + i], before[j * w + i - 3], 1e-12);
        }
    }
}

#[test]
fn advection_2d_scenario_stays_finite_and_bounded() {
    let res = simulate_2d_advection_scenario();
    assert_eq!(res.len(), 100 * 100);
    assert!(
        res.iter()
            .all(|v| v.is_finite() && *v > -1e-9 && *v < 1.0 + 1e-9)
    );
    assert!(res.iter().any(|&v| v > 0.3));
}

#[test]
fn advection_3d_unit_cfl_shifts_exactly_in_z() {
    let (w, h, d) = (8, 8, 12);
    let mut mesh = Mesh3D::new(w, h, d, (1.0, 1.0, 1.0), |_, _, z| z * z);
    let before: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();
    let dt = mesh.dz;
    let border = |i: usize, j: usize, k: usize, w: usize, h: usize, d: usize| {
        i == 0 || j == 0 || k == 0 || i == w - 1 || j == h - 1 || k == d - 1
    };
    let res = solve_advection_3d(&mut mesh, (0.0, 0.0, 1.0), dt, 2, border);
    let plane = w * h;
    for k in 4..d - 1 {
        for j in 1..h - 1 {
            for i in 1..w - 1 {
                let idx = k * plane + j * w + i;
                assert_approx_eq!(res[idx], before[idx - 2 * plane], 1e-12);
            }
        }
    }
}

#[test]
fn advection_3d_scenario_is_bounded() {
    let res = simulate_3d_advection_scenario();
    assert_eq!(res.len(), 30 * 30 * 30);
    assert!(
        res.iter()
            .all(|v| v.is_finite() && *v > -1e-9 && *v < 1.0 + 1e-9)
    );
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_advection_1d_with_matching_inflow_conserves_a_uniform_state(
        v in 0.1f64..2.0,
        steps in 1usize..20,
    ) {
        let n = 50;
        let mut mesh = Mesh::new(n, 1.0, |_| 1.0);
        let dt = 0.1 / n as f64 / v;
        let res = solve_advection_1d(&mut mesh, v, dt, steps, || (1.0, 1.0));
        for x in res {
            prop_assert!((x - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn prop_upwind_advection_obeys_the_maximum_principle(
        v in 0.1f64..1.0,
        steps in 1usize..30,
        seed in 0u64..1000,
    ) {
        let n = 40;
        let mut s = seed;
        let init: Vec<f64> = (0..n).map(|_| { s = s.wrapping_mul(6364136223846793005).wrapping_add(1); ((s >> 33) % 100) as f64 / 100.0 }).collect();
        let mut mesh = Mesh::new(n, 1.0, |x| init[((x * n as f64) as usize).min(n - 1)]);
        let dt = 0.9 / n as f64 / v;
        let res = solve_advection_1d(&mut mesh, v, dt, steps, || (0.0, 0.0));
        for x in res {
            prop_assert!((-1e-12..=1.0 + 1e-12).contains(&x));
        }
    }
}
