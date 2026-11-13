//! Coverage round 2: finite-volume Euler solver (grid geometry, conservation across every
//! scheme, mirror symmetry, contact preservation, agreement between Riemann solvers) and the
//! AMR driver (validation, equivalence with the uniform solver, conservation).

use std::f64::consts::PI;

use rssn::sim::euler_amr::AmrConfig;
use rssn::sim::euler_amr::solve_amr;
use rssn::sim::euler_nd::*;

const GAMMA: f64 = 1.4;

fn sod(x: [f64; 3]) -> State {
    if x[0] < 0.5 {
        [1.0, 0.0, 0.0, 0.0, 1.0]
    } else {
        [0.125, 0.0, 0.0, 0.0, 0.1]
    }
}

fn sod_mirrored(x: [f64; 3]) -> State {
    // x -> 1 - x with u -> -u
    let [r, u, v, w, p] = sod([1.0 - x[0], 0.0, 0.0]);
    [r, -u, v, w, p]
}

fn all_schemes() -> Vec<(Reconstruction, RiemannSolver)> {
    let mut v = Vec::new();
    let recons = [
        Reconstruction::FirstOrder,
        Reconstruction::Muscl(SlopeLimiter::Minmod),
        Reconstruction::Muscl(SlopeLimiter::VanLeer),
        Reconstruction::Muscl(SlopeLimiter::Mc),
        Reconstruction::Muscl(SlopeLimiter::Superbee),
        Reconstruction::Weno5,
    ];
    for r in recons {
        for s in [
            RiemannSolver::Hllc,
            RiemannSolver::Roe,
            RiemannSolver::Rusanov,
        ] {
            v.push((r, s));
        }
    }
    v
}

// ------------------------------------------------------------- geometry

#[test]
fn grid_geometry_is_consistent() {
    assert!(FvGrid::new(&[], &[], &[]).is_err());
    assert!(FvGrid::new(&[4, 4, 4, 4], &[0.0; 4], &[1.0; 4]).is_err());
    assert!(FvGrid::new(&[1], &[0.0], &[1.0]).is_err());
    assert!(FvGrid::new(&[4, 4], &[0.0], &[1.0, 1.0]).is_err());
    assert!(FvGrid::new(&[4, 4], &[0.0, 0.0], &[1.0]).is_err());
    let g = FvGrid::new(&[4, 6, 5], &[0.0, -1.0, 2.0], &[2.0, 2.0, 3.0]).unwrap();
    assert_eq!(g.len(), 120);
    assert!(!g.is_empty());
    assert!(
        (g.h(0) - 0.5).abs() < 1e-15
            && (g.h(1) - 0.5).abs() < 1e-15
            && (g.h(2) - 0.2).abs() < 1e-15
    );
    assert!((g.cell_volume() - 0.5 * 0.5 * 0.2).abs() < 1e-15);
    // total volume of all cells = domain volume
    assert!((g.cell_volume() * g.len() as f64 - 2.0 * 3.0 * 1.0).abs() < 1e-12);
    for p in 0..g.len() {
        let ijk = g.ijk(p);
        assert_eq!(g.index(ijk), p);
        let c = g.centre(p);
        for a in 0..3 {
            assert!(c[a] > g.lo[a] && c[a] < g.hi[a]);
        }
    }
    assert_eq!(g.centre(0), [0.25, -0.75, 2.1]);
    // an unused axis reports a zero coordinate in lower dimensions
    let g1 = FvGrid::new(&[8], &[1.0], &[2.0]).unwrap();
    assert_eq!(g1.centre(3)[1], 0.0);
    assert_eq!(g1.centre(3)[2], 0.0);
}

// ------------------------------------------------------- Euler solutions

#[test]
fn every_scheme_conserves_mass_momentum_and_energy_on_a_periodic_wave() {
    let grid = FvGrid::new(&[32], &[0.0], &[1.0]).unwrap();
    let init = |x: [f64; 3]| {
        [
            1.0 + 0.3 * (2.0 * PI * x[0]).sin(),
            0.7,
            0.0,
            0.0,
            1.0 + 0.1 * (4.0 * PI * x[0]).cos(),
        ]
    };
    for (recon, riemann) in all_schemes() {
        let mut cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Periodic);
        cfg.recon = recon;
        cfg.riemann = riemann;
        let s0 = solve(&grid, &cfg, &init, 0.0).unwrap();
        let s1 = solve(&grid, &cfg, &init, 0.25).unwrap();
        assert!(
            s1.steps > 0 && (s1.t - 0.25).abs() < 1e-12,
            "{recon:?}/{riemann:?}"
        );
        let (a, b) = (s0.totals(), s1.totals());
        for c in 0..5 {
            assert!(
                (a[c] - b[c]).abs() < 1e-12 * (1.0 + a[c].abs()),
                "{recon:?}/{riemann:?}: component {c}: {} -> {}",
                a[c],
                b[c]
            );
        }
        // positivity
        for p in 0..grid.len() {
            let w = s1.primitive(p);
            assert!(w[0] > 0.0 && w[4] > 0.0, "{recon:?}/{riemann:?}: {w:?}");
        }
    }
}

#[test]
fn t_end_zero_returns_the_initial_condition_untouched() {
    let grid = FvGrid::new(&[16], &[0.0], &[1.0]).unwrap();
    let cfg = FvConfig::euler(GAMMA);
    let s = solve(&grid, &cfg, &sod, 0.0).unwrap();
    assert_eq!(s.steps, 0);
    for p in 0..grid.len() {
        let w = s.primitive(p);
        let e = sod(grid.centre(p));
        for c in 0..5 {
            assert!((w[c] - e[c]).abs() < 1e-14);
        }
    }
}

#[test]
fn sod_tube_is_mirror_symmetric_for_every_scheme() {
    let grid = FvGrid::new(&[48], &[0.0], &[1.0]).unwrap();
    for (recon, riemann) in all_schemes() {
        let mut cfg = FvConfig::euler(GAMMA);
        cfg.recon = recon;
        cfg.riemann = riemann;
        let a = solve(&grid, &cfg, &sod, 0.12).unwrap();
        let b = solve(&grid, &cfg, &sod_mirrored, 0.12).unwrap();
        for i in 0..48 {
            let wa = a.primitive(i);
            let wb = b.primitive(47 - i);
            assert!(
                (wa[0] - wb[0]).abs() < 1e-9,
                "{recon:?}/{riemann:?} rho at {i}"
            );
            assert!(
                (wa[1] + wb[1]).abs() < 1e-9,
                "{recon:?}/{riemann:?} u at {i}"
            );
            assert!(
                (wa[4] - wb[4]).abs() < 1e-9,
                "{recon:?}/{riemann:?} p at {i}"
            );
        }
    }
}

#[test]
fn a_moving_contact_keeps_pressure_and_velocity_exactly_with_hllc_and_roe() {
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let init = |x: [f64; 3]| {
        let rho = if (0.25..0.5).contains(&x[0]) {
            2.0
        } else {
            1.0
        };
        [rho, 0.6, 0.0, 0.0, 1.0]
    };
    for riemann in [RiemannSolver::Hllc, RiemannSolver::Roe] {
        let mut cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Periodic);
        cfg.recon = Reconstruction::FirstOrder;
        cfg.riemann = riemann;
        let s = solve(&grid, &cfg, &init, 0.3).unwrap();
        for p in 0..grid.len() {
            let w = s.primitive(p);
            assert!((w[1] - 0.6).abs() < 1e-9, "{riemann:?}: u = {}", w[1]);
            assert!((w[4] - 1.0).abs() < 1e-9, "{riemann:?}: p = {}", w[4]);
        }
        // mass is moved, never created
        assert!((s.totals()[0] - (1.0 + 0.25)).abs() < 1e-12);
    }
}

#[test]
fn riemann_solvers_agree_on_a_smooth_wave_and_all_match_the_exact_advection() {
    let n = 64;
    let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
    let rho = |x: f64, t: f64| 1.0 + 0.2 * (2.0 * PI * (x - t)).sin();
    let init = |x: [f64; 3]| [rho(x[0], 0.0), 1.0, 0.0, 0.0, 1.0];
    let mut sols = Vec::new();
    for riemann in [
        RiemannSolver::Hllc,
        RiemannSolver::Roe,
        RiemannSolver::Rusanov,
    ] {
        let mut cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Periodic);
        cfg.recon = Reconstruction::Weno5;
        cfg.riemann = riemann;
        let s = solve(&grid, &cfg, &init, 0.5).unwrap();
        let err = s.l1_error(0, |x| rho(x[0], 0.5));
        assert!(err < 2e-3, "{riemann:?}: L1 error {err}");
        sols.push(s);
    }
    for k in 1..3 {
        for p in 0..n {
            assert!(
                (sols[0].primitive(p)[0] - sols[k].primitive(p)[0]).abs() < 5e-3,
                "solver {k} differs from HLLC at {p}"
            );
        }
    }
}

#[test]
fn invalid_euler_setups_are_errors() {
    let grid = FvGrid::new(&[16], &[0.0], &[1.0]).unwrap();
    for cfl in [0.0, 1.0, -0.2, 2.0, f64::NAN] {
        let mut cfg = FvConfig::euler(GAMMA);
        cfg.cfl = cfl;
        assert!(solve(&grid, &cfg, &sod, 0.1).is_err(), "cfl {cfl}");
    }
    // a vacuum/negative-density initial state must not run (or must at least not panic)
    let cfg = FvConfig::euler(GAMMA);
    let bad = |_x: [f64; 3]| [-1.0, 0.0, 0.0, 0.0, 1.0];
    let r = std::panic::catch_unwind(|| {
        solve(
            &FvGrid::new(&[16], &[0.0], &[1.0]).unwrap(),
            &cfg,
            &bad,
            0.01,
        )
    });
    assert!(r.is_ok(), "negative density must not panic");
    if let Ok(res) = r {
        assert!(res.is_err(), "negative density should be rejected");
    }
    // gravity on shallow water is rejected
    let sw = FvConfig::shallow_water(9.81);
    let src = Sources {
        gravity: [0.0, -1.0, 0.0],
        ..Sources::default()
    };
    let flow = |_x: [f64; 3]| [1.0, 0.0, 0.0, 0.0, 0.0];
    assert!(solve_with_sources(&grid, &sw, &flow, &src, 0.1).is_err());
}

// ------------------------------------------------------------------- AMR

fn amr_scheme() -> FvConfig {
    let mut c = FvConfig::euler(GAMMA);
    c.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
    c.with_axis(0, GasBc::Reflective)
}

#[test]
fn amr_validation_errors() {
    let grid = FvGrid::new(&[32], &[0.0], &[1.0]).unwrap();
    let good = AmrConfig::new(amr_scheme());
    assert!(solve_amr(&grid, &good, &sod, 0.01).is_ok());
    // cfl out of range
    let mut cfg = good;
    cfg.fv.cfl = 1.5;
    assert!(solve_amr(&grid, &cfg, &sod, 0.01).is_err());
    // odd base grid
    let odd = FvGrid::new(&[31], &[0.0], &[1.0]).unwrap();
    assert!(solve_amr(&odd, &good, &sod, 0.01).is_err());
    // one periodic face only
    let mut cfg = good;
    cfg.fv = cfg.fv.with_face(0, GasBc::Periodic);
    assert!(solve_amr(&grid, &cfg, &sod, 0.01).is_err());
    // non-reflecting boundaries unsupported
    let mut cfg = good;
    cfg.fv = cfg
        .fv
        .with_axis(0, GasBc::NonReflecting([1.0, 0.0, 0.0, 0.0, 1.0]));
    assert!(solve_amr(&grid, &cfg, &sod, 0.01).is_err());
    // inadmissible initial state
    let bad = |_x: [f64; 3]| [1.0, 0.0, 0.0, 0.0, -1.0];
    assert!(solve_amr(&grid, &good, &bad, 0.01).is_err());
    // shallow water in 3D
    let g3 = FvGrid::new(&[4, 4, 4], &[0.0; 3], &[1.0; 3]).unwrap();
    let sw = AmrConfig::new(FvConfig::shallow_water(9.81));
    assert!(solve_amr(&g3, &sw, &|_x: [f64; 3]| [1.0, 0.0, 0.0, 0.0, 0.0], 0.01).is_err());
}

#[test]
fn amr_with_no_refinement_levels_equals_the_uniform_solver() {
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let mut cfg = AmrConfig::new(amr_scheme());
    cfg.max_level = 0;
    let amr = solve_amr(&grid, &cfg, &sod, 0.1).unwrap();
    assert_eq!(amr.cells_per_level(), vec![40]);
    let uni = solve(&grid, &amr_scheme(), &sod, 0.1).unwrap();
    let prim = amr.primitive_on_level(0);
    assert_eq!(prim.len(), 40);
    for p in 0..40 {
        let a = prim[p];
        let b = uni.primitive(p);
        for c in 0..5 {
            assert!(
                (a[c] - b[c]).abs() < 1e-9,
                "cell {p} comp {c}: {} vs {}",
                a[c],
                b[c]
            );
        }
    }
}

#[test]
fn amr_does_not_refine_a_uniform_state_and_is_deterministic() {
    let grid = FvGrid::new(&[16, 16], &[0.0; 2], &[1.0; 2]).unwrap();
    let uniform = |_x: [f64; 3]| [1.0, 0.2, -0.1, 0.0, 1.0];
    let mut cfg = AmrConfig::new(
        FvConfig::euler(GAMMA)
            .with_axis(0, GasBc::Periodic)
            .with_axis(1, GasBc::Periodic),
    );
    cfg.max_level = 2;
    let sol = solve_amr(&grid, &cfg, &uniform, 0.1).unwrap();
    assert_eq!(sol.cells_per_level()[0], 256);
    assert!(
        sol.cells_per_level().iter().skip(1).all(|&c| c == 0),
        "{:?}",
        sol.cells_per_level()
    );
    // a uniform flow stays uniform
    for w in sol.primitive_on_level(0) {
        assert!(
            (w[0] - 1.0).abs() < 1e-12 && (w[1] - 0.2).abs() < 1e-12 && (w[4] - 1.0).abs() < 1e-12
        );
    }
    // determinism with refinement
    let blob = |x: [f64; 3]| {
        let r2 = (x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2);
        [1.0 + (-60.0 * r2).exp(), 0.0, 0.0, 0.0, 1.0]
    };
    let a = solve_amr(&grid, &cfg, &blob, 0.05).unwrap();
    let b = solve_amr(&grid, &cfg, &blob, 0.05).unwrap();
    assert_eq!(a.cells_per_level(), b.cells_per_level());
    assert_eq!(a.steps, b.steps);
    for (x, y) in a.totals().iter().zip(b.totals().iter()) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    // periodic refluxed run conserves all totals
    let t0 = solve_amr(&grid, &cfg, &blob, 0.0).unwrap().totals();
    let t1 = a.totals();
    for c in 0..5 {
        assert!(
            (t0[c] - t1[c]).abs() < 1e-10 * (1.0 + t0[c].abs()),
            "component {c}: {} -> {}",
            t0[c],
            t1[c]
        );
    }
    assert!(
        a.cells_per_level()[1] > 0,
        "the density blob must trigger refinement"
    );
}
