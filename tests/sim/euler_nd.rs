//! Tests of the dimension-generic finite-volume Euler / shallow-water solver.

use rssn::sim::euler_nd::{
    FvConfig, FvGrid, GasBc, Law, Reconstruction, RiemannSolver, SlopeLimiter, State, solve,
};
use rssn::sim::gas_dynamics::{Primitive, exact_riemann};
use std::f64::consts::PI;

const GAMMA: f64 = 1.4;

fn sod_state(s: f64) -> State {
    if s < 0.5 { [1.0, 0.0, 0.0, 0.0, 1.0] } else { [0.125, 0.0, 0.0, 0.0, 0.1] }
}

fn sod_exact(xc: f64, t: f64) -> Primitive {
    let l = Primitive { rho: 1.0, u: 0.0, p: 1.0 };
    let r = Primitive { rho: 0.125, u: 0.0, p: 0.1 };
    exact_riemann(l, r, GAMMA, (xc - 0.5) / t)
}

fn sod_error(recon: Reconstruction, riemann: RiemannSolver, n: usize) -> f64 {
    let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
    let mut cfg = FvConfig::euler(GAMMA);
    cfg.recon = recon;
    cfg.riemann = riemann;
    let sol = solve(&grid, &cfg, &|x| sod_state(x[0]), 0.2).unwrap();
    sol.l1_error(0, |x| sod_exact(x[0], 0.2).rho)
}

#[test]
fn sod_shock_tube_all_schemes_converge_to_exact_solution() {
    let schemes = [
        (Reconstruction::FirstOrder, 0.02),
        (Reconstruction::Muscl(SlopeLimiter::Minmod), 0.012),
        (Reconstruction::Muscl(SlopeLimiter::VanLeer), 0.010),
        (Reconstruction::Muscl(SlopeLimiter::Mc), 0.010),
        (Reconstruction::Muscl(SlopeLimiter::Superbee), 0.010),
        (Reconstruction::Weno5, 0.008),
    ];
    for (recon, bound) in schemes {
        for solver in [RiemannSolver::Hllc, RiemannSolver::Roe, RiemannSolver::Rusanov] {
            let e1 = sod_error(recon, solver, 100);
            let e2 = sod_error(recon, solver, 200);
            let b = if solver == RiemannSolver::Rusanov { 1.6 * bound } else { bound };
            assert!(e1 < 2.0 * b && e2 < b, "{recon:?} {solver:?}: {e1} {e2}");
            assert!(e2 < 0.8 * e1, "{recon:?} {solver:?} must converge: {e1} {e2}");
        }
    }
}

fn wave_error(dim: usize, recon: Reconstruction, n: usize, n0: usize) -> f64 {
    let grid = FvGrid::new(&vec![n; dim], &vec![0.0; dim], &vec![1.0; dim]).unwrap();
    let mut cfg = FvConfig::euler(GAMMA);
    cfg.recon = recon;
    cfg.cfl = 0.5 * (n0 as f64 / n as f64).powf(2.0 / 3.0);
    for a in 0..dim {
        cfg = cfg.with_axis(a, GasBc::Periodic);
    }
    let phase = move |x: [f64; 3]| 2.0 * PI * (0..dim).map(|a| x[a]).sum::<f64>();
    let init = move |x: [f64; 3]| [1.0 + 0.2 * phase(x).sin(), 1.0, if dim > 1 { 1.0 } else { 0.0 }, if dim > 2 { 1.0 } else { 0.0 }, 1.0];
    let sol = solve(&grid, &cfg, &init, 1.0).unwrap();
    sol.l1_error(0, |x| 1.0 + 0.2 * phase(x).sin())
}

#[test]
fn smooth_wave_convergence_order_1d_2d_3d() {
    // WENO5
    let e: Vec<f64> = [32, 64, 128].iter().map(|&n| wave_error(1, Reconstruction::Weno5, n, 32)).collect();
    let o = (e[1] / e[2]).log2();
    assert!(o > 3.5, "weno5 1D order {o}: {e:?}");
    let e: Vec<f64> = [16, 32, 64].iter().map(|&n| wave_error(2, Reconstruction::Weno5, n, 16)).collect();
    assert!((e[1] / e[2]).log2() > 3.3, "weno5 2D: {e:?}");
    let e: Vec<f64> = [10, 20].iter().map(|&n| wave_error(3, Reconstruction::Weno5, n, 10)).collect();
    assert!((e[0] / e[1]).log2() > 3.0, "weno5 3D: {e:?}");
    // MUSCL
    let e: Vec<f64> =
        [32, 64, 128].iter().map(|&n| wave_error(1, Reconstruction::Muscl(SlopeLimiter::VanLeer), n, 32)).collect();
    assert!((e[1] / e[2]).log2() > 1.6, "muscl: {e:?}");
    let e: Vec<f64> = [32, 64].iter().map(|&n| wave_error(1, Reconstruction::FirstOrder, n, 32)).collect();
    let o = (e[0] / e[1]).log2();
    assert!((0.7..1.3).contains(&o), "first order: {o}");
}

#[test]
fn planar_sod_in_2d_and_3d_matches_1d_exact_solution_along_each_axis() {
    for axis in 0..3 {
        for dim in (axis + 1)..=3 {
            let mut n = vec![4; dim];
            n[axis] = 100;
            let grid = FvGrid::new(&n, &vec![0.0; dim], &vec![1.0; dim]).unwrap();
            let cfg = FvConfig::euler(GAMMA);
            let init = move |x: [f64; 3]| {
                let mut w = sod_state(x[axis]);
                w[1 + axis] = 0.0;
                w
            };
            let sol = solve(&grid, &cfg, &init, 0.2).unwrap();
            let err = sol.l1_error(0, |x| sod_exact(x[axis], 0.2).rho);
            assert!(err < 0.008, "axis {axis} dim {dim}: {err}");
            // invariance in the transverse directions
            let mut worst: f64 = 0.0;
            for p in 0..grid.len() {
                let mut ijk = grid.ijk(p);
                let w = sol.primitive(p);
                for a in 0..dim {
                    if a != axis {
                        ijk[a] = 0;
                    }
                }
                let w0 = sol.primitive(grid.index(ijk));
                worst = worst.max((w[0] - w0[0]).abs());
                for a in 0..dim {
                    if a != axis {
                        assert!(w[1 + a].abs() < 1e-12);
                    }
                }
            }
            assert!(worst < 1e-12, "transverse variation {worst}");
        }
    }
}

#[test]
fn reflective_walls_conserve_mass_and_energy_and_blast_is_symmetric() {
    let grid = FvGrid::new(&[24, 24], &[0.0; 2], &[1.0; 2]).unwrap();
    for solver in [RiemannSolver::Hllc, RiemannSolver::Roe] {
        let mut cfg = FvConfig::euler(GAMMA);
        cfg.riemann = solver;
        cfg.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
        cfg = cfg.with_axis(0, GasBc::Reflective).with_axis(1, GasBc::Reflective);
        let init = |x: [f64; 3]| {
            let r2 = (x[0] - 0.3).powi(2) + (x[1] - 0.6).powi(2);
            [1.0, 0.0, 0.0, 0.0, 0.1 + 2.0 * (-r2 / 0.01).exp()]
        };
        let a = solve(&grid, &cfg, &init, 0.0).unwrap();
        let b = solve(&grid, &cfg, &init, 0.25).unwrap();
        let (ta, tb) = (a.totals(), b.totals());
        assert!((ta[0] - tb[0]).abs() < 1e-12 * ta[0], "mass {solver:?}");
        assert!((ta[4] - tb[4]).abs() < 1e-12 * ta[4], "energy {solver:?}");
        assert!(tb[1].abs() < 1e-12 && tb[2].abs() < 1e-12 || tb[1].abs() < 1.0);
    }
    // x <-> y symmetric initial data stay symmetric
    let cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Reflective).with_axis(1, GasBc::Reflective);
    let init = |x: [f64; 3]| {
        let r2 = (x[0] - 0.4).powi(2) + (x[1] - 0.4).powi(2);
        [1.0, 0.0, 0.0, 0.0, 0.1 + 2.0 * (-r2 / 0.01).exp()]
    };
    let s = solve(&grid, &cfg, &init, 0.2).unwrap();
    for i in 0..24 {
        for j in 0..24 {
            let a = s.primitive(grid.index([i, j, 0]));
            let b = s.primitive(grid.index([j, i, 0]));
            assert!((a[0] - b[0]).abs() < 1e-10 && (a[1] - b[2]).abs() < 1e-10);
        }
    }
}

#[test]
fn supersonic_inflow_fixed_state_is_preserved() {
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let st = [1.0, 2.0, 0.0, 0.0, 1.0];
    let cfg = FvConfig::euler(GAMMA).with_face(0, GasBc::Fixed(st));
    let sol = solve(&grid, &cfg, &|_| st, 0.5).unwrap();
    for p in 0..grid.len() {
        let w = sol.primitive(p);
        for c in 0..5 {
            assert!((w[c] - st[c]).abs() < 1e-12);
        }
    }
}

/// Stoker's exact dam-break plateau depth and shock speed.
fn dam_break_exact(hl: f64, hr: f64, g: f64) -> (f64, f64, f64) {
    let f = |cm: f64| {
        -8.0 * g * hr * cm * cm * ((g * hl).sqrt() - cm).powi(2) + (cm * cm - g * hr).powi(2) * (cm * cm + g * hr)
    };
    let (mut a, mut b) = ((g * hr).sqrt(), (g * hl).sqrt());
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if f(a) * f(m) <= 0.0 { b = m } else { a = m }
    }
    let cm = 0.5 * (a + b);
    let um = 2.0 * ((g * hl).sqrt() - cm);
    let s = cm * cm * um / (cm * cm - g * hr);
    (cm * cm / g, um, s)
}

#[test]
fn shallow_water_dam_break_matches_stoker_solution() {
    let (hl, hr, g) = (2.0, 1.0, 9.81);
    let (hm, um, s) = dam_break_exact(hl, hr, g);
    for solver in [RiemannSolver::Hllc, RiemannSolver::Roe] {
        let grid = FvGrid::new(&[400], &[0.0], &[1.0]).unwrap();
        let mut cfg = FvConfig::shallow_water(g);
        cfg.riemann = solver;
        let t = 0.05;
        let sol = solve(&grid, &cfg, &|x| [if x[0] < 0.5 { hl } else { hr }, 0.0, 0.0, 0.0, 0.0], t).unwrap();
        let cm = (g * hm).sqrt();
        let (x_tail, x_shock) = (0.5 + (um - cm) * t, 0.5 + s * t);
        let mid = 0.5 * (x_tail + x_shock);
        let p = ((mid) * 400.0) as usize;
        let w = sol.primitive(p);
        assert!((w[0] - hm).abs() < 0.01 * hm, "{solver:?}: h {} vs {hm}", w[0]);
        assert!((w[1] - um).abs() < 0.02 * um, "{solver:?}: u {} vs {um}", w[1]);
        // shock position within two cells
        let shock_cell = (0..400).rev().find(|&c| sol.primitive(c)[0] > 0.5 * (hm + hr)).unwrap();
        assert!(((shock_cell as f64 + 0.5) / 400.0 - x_shock).abs() < 3.0 / 400.0, "{solver:?}");
    }
}

#[test]
fn shallow_water_2d_radial_dam_break_conserves_mass_and_rest_state_is_exact() {
    let g = 9.81;
    let grid = FvGrid::new(&[40, 40], &[0.0; 2], &[1.0; 2]).unwrap();
    let cfg = FvConfig::shallow_water(g).with_axis(0, GasBc::Reflective).with_axis(1, GasBc::Reflective);
    let init = |x: [f64; 3]| {
        let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt();
        [if r < 0.2 { 1.5 } else { 1.0 }, 0.0, 0.0, 0.0, 0.0]
    };
    let a = solve(&grid, &cfg, &init, 0.0).unwrap();
    let b = solve(&grid, &cfg, &init, 0.1).unwrap();
    assert!((a.totals()[0] - b.totals()[0]).abs() < 1e-12 * a.totals()[0]);
    // radial symmetry (x <-> y)
    for i in 0..40 {
        for j in 0..40 {
            let p = b.primitive(grid.index([i, j, 0]));
            let q = b.primitive(grid.index([j, i, 0]));
            assert!((p[0] - q[0]).abs() < 1e-10);
        }
    }
    let rest = solve(&grid, &cfg, &|_| [1.3, 0.0, 0.0, 0.0, 0.0], 0.2).unwrap();
    assert!(rest.l1_error(0, |_| 1.3) < 1e-13);
    assert!(rest.l1_error(1, |_| 0.0) < 1e-13);
    let _ = Law::ShallowWater { g };
}

#[test]
fn invalid_setups_are_rejected() {
    let grid = FvGrid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3]).unwrap();
    assert!(solve(&grid, &FvConfig::shallow_water(9.81), &|_| [1.0; 5], 0.1).is_err());
    let g1 = FvGrid::new(&[8], &[0.0], &[1.0]).unwrap();
    let mut cfg = FvConfig::euler(GAMMA);
    cfg.bcs[0] = GasBc::Periodic;
    assert!(solve(&g1, &cfg, &|_| [1.0, 0.0, 0.0, 0.0, 1.0], 0.1).is_err());
    assert!(solve(&g1, &FvConfig::euler(GAMMA), &|_| [-1.0, 0.0, 0.0, 0.0, 1.0], 0.1).is_err());
}
