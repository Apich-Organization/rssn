//! Tests of the dimension-generic finite-volume Euler / shallow-water solver.

use rssn::sim::euler_nd::FvConfig;
use rssn::sim::euler_nd::FvGrid;
use rssn::sim::euler_nd::GasBc;
use rssn::sim::euler_nd::Law;
use rssn::sim::euler_nd::Reconstruction;
use rssn::sim::euler_nd::RiemannSolver;
use rssn::sim::euler_nd::SlopeLimiter;
use rssn::sim::euler_nd::Sources;
use rssn::sim::euler_nd::State;
use rssn::sim::euler_nd::solve;
use rssn::sim::euler_nd::solve_with_sources;
use rssn::sim::gas_dynamics::Primitive;
use rssn::sim::gas_dynamics::exact_riemann;
use std::f64::consts::PI;

const GAMMA: f64 = 1.4;

fn sod_state(s: f64) -> State {
    if s < 0.5 {
        [1.0, 0.0, 0.0, 0.0, 1.0]
    } else {
        [0.125, 0.0, 0.0, 0.0, 0.1]
    }
}

fn sod_exact(
    xc: f64,
    t: f64,
) -> Primitive {
    let l = Primitive {
        rho: 1.0,
        u: 0.0,
        p: 1.0,
    };
    let r = Primitive {
        rho: 0.125,
        u: 0.0,
        p: 0.1,
    };
    exact_riemann(l, r, GAMMA, (xc - 0.5) / t)
}

fn sod_error(
    recon: Reconstruction,
    riemann: RiemannSolver,
    n: usize,
) -> f64 {
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
        for solver in [
            RiemannSolver::Hllc,
            RiemannSolver::Roe,
            RiemannSolver::Rusanov,
        ] {
            let e1 = sod_error(recon, solver, 100);
            let e2 = sod_error(recon, solver, 200);
            let b = if solver == RiemannSolver::Rusanov {
                1.6 * bound
            } else {
                bound
            };
            assert!(e1 < 2.0 * b && e2 < b, "{recon:?} {solver:?}: {e1} {e2}");
            assert!(
                e2 < 0.8 * e1,
                "{recon:?} {solver:?} must converge: {e1} {e2}"
            );
        }
    }
}

fn wave_error(
    dim: usize,
    recon: Reconstruction,
    n: usize,
    n0: usize,
) -> f64 {
    let grid = FvGrid::new(&vec![n; dim], &vec![0.0; dim], &vec![1.0; dim]).unwrap();
    let mut cfg = FvConfig::euler(GAMMA);
    cfg.recon = recon;
    cfg.cfl = 0.5 * (n0 as f64 / n as f64).powf(2.0 / 3.0);
    for a in 0..dim {
        cfg = cfg.with_axis(a, GasBc::Periodic);
    }
    let phase = move |x: [f64; 3]| 2.0 * PI * (0..dim).map(|a| x[a]).sum::<f64>();
    let init = move |x: [f64; 3]| {
        [
            1.0 + 0.2 * phase(x).sin(),
            1.0,
            if dim > 1 { 1.0 } else { 0.0 },
            if dim > 2 { 1.0 } else { 0.0 },
            1.0,
        ]
    };
    let sol = solve(&grid, &cfg, &init, 1.0).unwrap();
    sol.l1_error(0, |x| 1.0 + 0.2 * phase(x).sin())
}

#[test]
fn smooth_wave_convergence_order_1d_2d_3d() {
    // WENO5
    let e: Vec<f64> = [32, 64, 128]
        .iter()
        .map(|&n| wave_error(1, Reconstruction::Weno5, n, 32))
        .collect();
    let o = (e[1] / e[2]).log2();
    assert!(o > 3.5, "weno5 1D order {o}: {e:?}");
    let e: Vec<f64> = [16, 32, 64]
        .iter()
        .map(|&n| wave_error(2, Reconstruction::Weno5, n, 16))
        .collect();
    assert!((e[1] / e[2]).log2() > 3.3, "weno5 2D: {e:?}");
    let e: Vec<f64> = [10, 20]
        .iter()
        .map(|&n| wave_error(3, Reconstruction::Weno5, n, 10))
        .collect();
    assert!((e[0] / e[1]).log2() > 3.0, "weno5 3D: {e:?}");
    // MUSCL
    let e: Vec<f64> = [32, 64, 128]
        .iter()
        .map(|&n| wave_error(1, Reconstruction::Muscl(SlopeLimiter::VanLeer), n, 32))
        .collect();
    assert!((e[1] / e[2]).log2() > 1.6, "muscl: {e:?}");
    let e: Vec<f64> = [32, 64]
        .iter()
        .map(|&n| wave_error(1, Reconstruction::FirstOrder, n, 32))
        .collect();
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
        cfg = cfg
            .with_axis(0, GasBc::Reflective)
            .with_axis(1, GasBc::Reflective);
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
    let cfg = FvConfig::euler(GAMMA)
        .with_axis(0, GasBc::Reflective)
        .with_axis(1, GasBc::Reflective);
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
fn dam_break_exact(
    hl: f64,
    hr: f64,
    g: f64,
) -> (f64, f64, f64) {
    let f = |cm: f64| {
        -8.0 * g * hr * cm * cm * ((g * hl).sqrt() - cm).powi(2)
            + (cm * cm - g * hr).powi(2) * (cm * cm + g * hr)
    };
    let (mut a, mut b) = ((g * hr).sqrt(), (g * hl).sqrt());
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if f(a) * f(m) <= 0.0 {
            b = m
        } else {
            a = m
        }
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
        let sol = solve(
            &grid,
            &cfg,
            &|x| [if x[0] < 0.5 { hl } else { hr }, 0.0, 0.0, 0.0, 0.0],
            t,
        )
        .unwrap();
        let cm = (g * hm).sqrt();
        let (x_tail, x_shock) = (0.5 + (um - cm) * t, 0.5 + s * t);
        let mid = 0.5 * (x_tail + x_shock);
        let p = ((mid) * 400.0) as usize;
        let w = sol.primitive(p);
        assert!(
            (w[0] - hm).abs() < 0.01 * hm,
            "{solver:?}: h {} vs {hm}",
            w[0]
        );
        assert!(
            (w[1] - um).abs() < 0.02 * um,
            "{solver:?}: u {} vs {um}",
            w[1]
        );
        // shock position within two cells
        let shock_cell = (0..400)
            .rev()
            .find(|&c| sol.primitive(c)[0] > 0.5 * (hm + hr))
            .unwrap();
        assert!(
            ((shock_cell as f64 + 0.5) / 400.0 - x_shock).abs() < 3.0 / 400.0,
            "{solver:?}"
        );
    }
}

#[test]
fn shallow_water_2d_radial_dam_break_conserves_mass_and_rest_state_is_exact() {
    let g = 9.81;
    let grid = FvGrid::new(&[40, 40], &[0.0; 2], &[1.0; 2]).unwrap();
    let cfg = FvConfig::shallow_water(g)
        .with_axis(0, GasBc::Reflective)
        .with_axis(1, GasBc::Reflective);
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
    assert!(
        solve(
            &g1,
            &FvConfig::euler(GAMMA),
            &|_| [-1.0, 0.0, 0.0, 0.0, 1.0],
            0.1
        )
        .is_err()
    );
}

// ------------------------------------------------------------ source terms

const G: f64 = 9.81;

/// A smooth bump and a discontinuous step, both of amplitude of order one.
fn smooth_bed(x: [f64; 3]) -> f64 {
    0.4 * (-(x[0] - 0.5).powi(2) / 0.01).exp() + 0.1 * (2.0 * PI * x[0]).sin()
}

fn step_bed(x: [f64; 3]) -> f64 {
    if x[0] > 0.45 && x[0] < 0.7 {
        0.5
    } else {
        0.0
    }
}

fn bed_2d(x: [f64; 3]) -> f64 {
    0.3 * (-((x[0] - 0.5).powi(2) + (x[1] - 0.4).powi(2)) / 0.03).exp()
}

/// Largest speed and largest free-surface deviation from `eta0` of a run.
fn lake_at_rest(
    dim: usize,
    n: usize,
    bed: &dyn Fn([f64; 3]) -> f64,
    eta0: f64,
    cfg: &FvConfig,
    t_end: f64,
) -> (f64, f64) {
    let grid = FvGrid::new(&vec![n; dim], &vec![0.0; dim], &vec![1.0; dim]).unwrap();
    let sources = Sources {
        bed: Some(bed),
        ..Sources::default()
    };
    let sol = solve_with_sources(
        &grid,
        cfg,
        &|x| [eta0 - bed(x), 0.0, 0.0, 0.0, 0.0],
        &sources,
        t_end,
    )
    .unwrap();
    let (mut speed, mut surface) = (0.0_f64, 0.0_f64);
    for p in 0..grid.len() {
        let w = sol.primitive(p);
        speed = speed.max(w[1].abs()).max(w[2].abs());
        surface = surface.max((w[0] + bed(grid.centre(p)) - eta0).abs());
    }
    (speed, surface)
}

#[test]
fn the_lake_at_rest_is_preserved_to_round_off_by_every_scheme() {
    let reconstructions = [
        Reconstruction::FirstOrder,
        Reconstruction::Muscl(SlopeLimiter::Minmod),
        Reconstruction::Muscl(SlopeLimiter::Mc),
        Reconstruction::Muscl(SlopeLimiter::Superbee),
        Reconstruction::Weno5,
    ];
    for recon in reconstructions {
        for riemann in [
            RiemannSolver::Hllc,
            RiemannSolver::Roe,
            RiemannSolver::Rusanov,
        ] {
            let mut cfg = FvConfig::shallow_water(G).with_axis(0, GasBc::Reflective);
            cfg.recon = recon;
            cfg.riemann = riemann;
            for (name, bed) in [
                ("smooth", &smooth_bed as &dyn Fn([f64; 3]) -> f64),
                ("step", &step_bed),
            ] {
                let (speed, surface) = lake_at_rest(1, 80, bed, 1.0, &cfg, 0.5);
                assert!(
                    speed < 1e-12 && surface < 1e-12,
                    "{recon:?} {riemann:?} {name}: speed {speed:e}, surface {surface:e}"
                );
            }
        }
    }
}

#[test]
fn the_two_dimensional_lake_at_rest_and_the_naive_scheme_it_replaces() {
    let mut cfg = FvConfig::shallow_water(G)
        .with_axis(0, GasBc::Reflective)
        .with_axis(1, GasBc::Reflective);
    cfg.recon = Reconstruction::Muscl(SlopeLimiter::Mc);
    let (speed, surface) = lake_at_rest(2, 40, &bed_2d, 1.0, &cfg, 0.3);
    assert!(
        speed < 1e-12 && surface < 1e-12,
        "speed {speed:e}, surface {surface:e}"
    );
    // Periodic in one direction, wall in the other.
    let mut periodic = FvConfig::shallow_water(G)
        .with_axis(0, GasBc::Periodic)
        .with_axis(1, GasBc::Reflective);
    periodic.recon = Reconstruction::Weno5;
    let (speed, surface) = lake_at_rest(
        2,
        40,
        &|x| 0.2 * (2.0 * PI * x[0]).sin() * x[1],
        1.0,
        &periodic,
        0.2,
    );
    assert!(
        speed < 1e-12 && surface < 1e-12,
        "periodic: speed {speed:e}, surface {surface:e}"
    );
    // Ignoring the bed (same depth profile, flat-bottom equations) is not at
    // rest: the pressure imbalance drives velocities of order one.
    let grid = FvGrid::new(&[80], &[0.0], &[1.0]).unwrap();
    let mut flat = FvConfig::shallow_water(G).with_axis(0, GasBc::Reflective);
    flat.recon = Reconstruction::Muscl(SlopeLimiter::Mc);
    let sol = solve(
        &grid,
        &flat,
        &|x| [1.0 - smooth_bed(x), 0.0, 0.0, 0.0, 0.0],
        0.1,
    )
    .unwrap();
    let speed = (0..grid.len())
        .map(|p| sol.primitive(p)[1].abs())
        .fold(0.0, f64::max);
    assert!(
        speed > 0.3,
        "the flat-bottom scheme moves the lake: {speed}"
    );
}

#[test]
fn a_perturbation_over_a_bed_conserves_mass_and_converges() {
    // A small hump on the free surface over a smooth bed, reflecting walls.
    let hump = |x: [f64; 3]| 1.0 + 0.01 * (-(x[0] - 0.25).powi(2) / 0.004).exp() - smooth_bed(x);
    let run = |n: usize| {
        let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let mut cfg = FvConfig::shallow_water(G).with_axis(0, GasBc::Reflective);
        cfg.recon = Reconstruction::Muscl(SlopeLimiter::Mc);
        let sources = Sources {
            bed: Some(&smooth_bed),
            ..Sources::default()
        };
        solve_with_sources(
            &grid,
            &cfg,
            &|x| [hump(x), 0.0, 0.0, 0.0, 0.0],
            &sources,
            0.04,
        )
        .unwrap()
    };
    let coarse = run(100);
    // Mass is conserved to round-off.
    let initial_mass: f64 = (0..100)
        .map(|p| hump(coarse.grid.centre(p)) * coarse.grid.cell_volume())
        .sum();
    assert!(
        (coarse.totals()[0] - initial_mass).abs() < 1e-12,
        "{} vs {initial_mass}",
        coarse.totals()[0]
    );
    // Self-convergence of the free surface.
    let surface = |sol: &rssn::sim::euler_nd::FvSolution| -> Vec<f64> {
        (0..sol.grid.len())
            .map(|p| sol.primitive(p)[0] + smooth_bed(sol.grid.centre(p)))
            .collect()
    };
    let (s100, s200, s400, s800) = (
        surface(&coarse),
        surface(&run(200)),
        surface(&run(400)),
        surface(&run(800)),
    );
    // L1 differences of successive refinements (averaging the finer grid onto the coarser).
    let l1 = |coarse: &[f64], fine: &[f64]| -> f64 {
        let ratio = fine.len() / coarse.len();
        coarse
            .iter()
            .enumerate()
            .map(|(i, c)| {
                (c - fine[ratio * i..ratio * (i + 1)].iter().sum::<f64>() / ratio as f64).abs()
            })
            .sum::<f64>()
            / coarse.len() as f64
    };
    let (e1, e2, e3) = (l1(&s100, &s200), l1(&s200, &s400), l1(&s400, &s800));
    assert!(
        e1 / e2 > 1.6 && e2 / e3 > 1.6,
        "successive L1 differences {e1:e} {e2:e} {e3:e}"
    );
    // The hump splits into two waves of about half the amplitude.
    let top = s800.iter().copied().fold(f64::NEG_INFINITY, f64::max) - 1.0;
    assert!(top > 0.003 && top < 0.008, "peak {top}");
}

#[test]
fn invalid_source_combinations_are_rejected() {
    let grid = FvGrid::new(&[20], &[0.0], &[1.0]).unwrap();
    let init = |_: [f64; 3]| [1.0, 0.0, 0.0, 0.0, 1.0];
    let flat = |_: [f64; 3]| 0.0;
    let bed = Sources {
        bed: Some(&flat),
        ..Sources::default()
    };
    assert!(
        solve_with_sources(&grid, &FvConfig::euler(GAMMA), &init, &bed, 0.01).is_err(),
        "a bed with Euler"
    );
    let gravity = Sources {
        gravity: [0.0, -1.0, 0.0],
        ..Sources::default()
    };
    assert!(
        solve_with_sources(&grid, &FvConfig::shallow_water(G), &init, &gravity, 0.01).is_err(),
        "gravity with shallow water"
    );
    // The defaults are the plain solver.
    let a = solve(&grid, &FvConfig::euler(GAMMA), &|x| sod_state(x[0]), 0.05).unwrap();
    let b = solve_with_sources(
        &grid,
        &FvConfig::euler(GAMMA),
        &|x| sod_state(x[0]),
        &Sources::default(),
        0.05,
    )
    .unwrap();
    assert_eq!(a.q, b.q);
}

#[test]
fn a_uniform_gas_in_free_fall_gains_exactly_g_t() {
    // Periodic along the field: the state stays uniform, v = g t, and the
    // energy gain is the kinetic energy, so the pressure does not change.
    let grid = FvGrid::new(&[8, 8], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
    let cfg = FvConfig::euler(GAMMA)
        .with_axis(0, GasBc::Periodic)
        .with_axis(1, GasBc::Periodic);
    let sources = Sources {
        gravity: [0.0, -2.0, 0.0],
        ..Sources::default()
    };
    let sol =
        solve_with_sources(&grid, &cfg, &|_| [1.3, 0.2, 0.0, 0.0, 0.9], &sources, 0.3).unwrap();
    for p in 0..grid.len() {
        let w = sol.primitive(p);
        assert!(
            (w[0] - 1.3).abs() < 1e-12 && (w[1] - 0.2).abs() < 1e-12,
            "{w:?}"
        );
        assert!((w[2] + 0.6).abs() < 1e-12, "v = g t: {}", w[2]);
        assert!((w[4] - 0.9).abs() < 1e-11, "pressure {}", w[4]);
    }
    // The total momentum grows by rho g V t.
    let momentum = sol.totals()[2];
    assert!((momentum - 1.3 * -2.0 * 0.3).abs() < 1e-11, "{momentum}");
}

#[test]
fn an_isothermal_atmosphere_is_held_to_the_order_of_the_scheme() {
    // rho = exp(-g x / T), p = rho T (T = g = 1) between two walls: the
    // residual motion is a truncation effect, first order because of the
    // zero-gradient ghost cells at the walls.
    let residual = |n: usize| {
        let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let mut cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Reflective);
        cfg.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
        let sources = Sources {
            gravity: [-1.0, 0.0, 0.0],
            ..Sources::default()
        };
        let init = |x: [f64; 3]| {
            let rho = (-x[0]).exp();
            [rho, 0.0, 0.0, 0.0, rho]
        };
        let sol = solve_with_sources(&grid, &cfg, &init, &sources, 0.5).unwrap();
        (0..n)
            .map(|p| sol.primitive(p)[1].abs())
            .fold(0.0, f64::max)
    };
    let (a, b, c) = (residual(40), residual(80), residual(160));
    assert!(a < 0.05, "the equilibrium is nearly held: {a}");
    assert!(
        a / b > 1.8 && b / c > 1.8,
        "residual speeds {a:e} {b:e} {c:e}"
    );
    // A uniform column is not in equilibrium: gravity accelerates its interior.
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let sources = Sources {
        gravity: [-1.0, 0.0, 0.0],
        ..Sources::default()
    };
    let cfg = FvConfig::euler(GAMMA).with_axis(0, GasBc::Reflective);
    let sol =
        solve_with_sources(&grid, &cfg, &|_| [1.0, 0.0, 0.0, 0.0, 1.0], &sources, 0.2).unwrap();
    assert!(
        sol.primitive(20)[1] < -0.05,
        "free fall of the interior: {}",
        sol.primitive(20)[1]
    );
}
