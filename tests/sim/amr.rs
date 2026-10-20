//! Block-structured AMR for the Euler equations.

use rssn::sim::euler_amr::{AmrConfig, solve_amr};
use rssn::sim::euler_nd::{FvConfig, FvGrid, GasBc, Reconstruction, SlopeLimiter, solve};
use rssn::sim::gas_dynamics::{Primitive, exact_riemann};
use std::f64::consts::PI;

fn sod(x: [f64; 3]) -> [f64; 5] {
    if x[0] < 0.5 { [1.0, 0.0, 0.0, 0.0, 1.0] } else { [0.125, 0.0, 0.0, 0.0, 0.1] }
}

fn scheme() -> FvConfig {
    let mut c = FvConfig::euler(1.4);
    c.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
    c.cfl = 0.4;
    c.with_axis(0, GasBc::Reflective)
}

fn density_l1_error(rho: &[[f64; 5]], n: usize) -> f64 {
    let l = Primitive { rho: 1.0, u: 0.0, p: 1.0 };
    let r = Primitive { rho: 0.125, u: 0.0, p: 0.1 };
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            (rho[i][0] - exact_riemann(l, r, 1.4, (x - 0.5) / 0.2).rho).abs() / n as f64
        })
        .sum()
}

#[test]
fn amr_sod_tube_matches_the_uniform_fine_grid_with_fewer_cells_and_conserves() {
    let grid = FvGrid::new(&[64], &[0.0], &[1.0]).unwrap();
    let mut cfg = AmrConfig::new(scheme());
    cfg.max_level = 2;
    cfg.threshold = 0.03;
    let sol = solve_amr(&grid, &cfg, &sod, 0.2).unwrap();
    let cells = sol.cells_per_level();
    println!("cells per level {cells:?}");
    assert!(cells[1] > 0 && cells[2] > 0, "refinement was created: {cells:?}");
    let fine_n = 256;
    let amr_err = density_l1_error(&sol.primitive_on_level(2), fine_n);
    let uni = |n: usize| {
        let g = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let s = solve(&g, &scheme(), &sod, 0.2).unwrap();
        (0..n).map(|i| s.primitive(i)).collect::<Vec<_>>()
    };
    let coarse_err = density_l1_error(&uni(64), 64);
    let fine_err = density_l1_error(&uni(fine_n), fine_n);
    println!("L1: coarse {coarse_err:.5}, amr {amr_err:.5}, uniform fine {fine_err:.5}");
    assert!(amr_err < 0.5 * coarse_err, "AMR {amr_err} vs coarse {coarse_err}");
    assert!(amr_err < 1.5 * fine_err, "AMR {amr_err} vs fine {fine_err}");
    let total_cells: usize = cells.iter().sum();
    assert!(total_cells < fine_n, "AMR uses {total_cells} cells against {fine_n}");
    // mass and energy are conserved to round-off (the waves do not reach the walls)
    let t = sol.totals();
    let (m0, e0) = (0.5 * 1.0 + 0.5 * 0.125, 0.5 * 1.0 / 0.4 + 0.5 * 0.1 / 0.4);
    assert!((t[0] - m0).abs() < 1e-12 && (t[4] - e0).abs() < 1e-12, "totals {t:?}");
}

#[test]
fn amr_conserves_with_refluxing_in_a_closed_2d_box() {
    // blast wave in a reflective box: refinement follows the front; totals stay exact
    let grid = FvGrid::new(&[32, 32], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
    let mut c = FvConfig::euler(1.4);
    c.recon = Reconstruction::Muscl(SlopeLimiter::Minmod);
    c.cfl = 0.4;
    let c = c.with_axis(0, GasBc::Reflective).with_axis(1, GasBc::Reflective);
    let mut cfg = AmrConfig::new(c);
    cfg.max_level = 2;
    cfg.threshold = 0.05;
    cfg.regrid_interval = 2;
    let blast = |x: [f64; 3]| {
        let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt();
        if r < 0.15 { [1.0, 0.0, 0.0, 0.0, 10.0] } else { [1.0, 0.0, 0.0, 0.0, 0.1] }
    };
    let sol = solve_amr(&grid, &cfg, &blast, 0.12).unwrap();
    let cells = sol.cells_per_level();
    assert!(cells[1] > 0, "{cells:?}");
    assert!(cells[0] + cells[1] + cells[2] < 32 * 32 * 16, "fewer cells than the uniform fine grid: {cells:?}");
    // initial totals from the fine-sampled initial data on the finest patches are close; compare against a restart at t = 0
    let s0 = solve_amr(&grid, &cfg, &blast, 1e-12).unwrap().totals();
    let s1 = sol.totals();
    for c in [0, 4] {
        assert!((s1[c] - s0[c]).abs() < 1e-11 * s0[c].abs(), "component {c}: {} vs {}", s1[c], s0[c]);
    }
    // x and y symmetry of the radially symmetric blast is kept by the hierarchy
    let rho = sol.primitive_on_level(2);
    let n = 128;
    let mut asym: f64 = 0.0;
    for j in 0..n {
        for i in 0..n {
            asym = asym.max((rho[i + n * j][0] - rho[j + n * i][0]).abs());
        }
    }
    assert!(asym < 0.05, "x/y asymmetry {asym}");
    let uniform = {
        let g = FvGrid::new(&[128, 128], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
        solve(&g, &cfg.fv, &blast, 0.12).unwrap()
    };
    let mut l1 = 0.0;
    for p in 0..128 * 128 {
        l1 += (uniform.primitive(p)[0] - rho[p][0]).abs() / (128.0 * 128.0);
    }
    assert!(l1 < 0.02, "AMR vs uniform fine L1 density difference {l1}");
}

#[test]
fn amr_periodic_density_wave_keeps_second_order_accuracy_across_refinement_boundaries() {
    let exact = |x: f64, t: f64| 1.0 + 0.2 * (2.0 * PI * (x - t)).sin();
    let init = |x: [f64; 3]| [exact(x[0], 0.0), 1.0, 0.0, 0.0, 1.0];
    let err = |n: usize| {
        let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let mut c = FvConfig::euler(1.4);
        c.recon = Reconstruction::Muscl(SlopeLimiter::Mc);
        c.cfl = 0.4;
        let c = c.with_axis(0, GasBc::Periodic);
        let mut cfg = AmrConfig::new(c);
        cfg.max_level = 1;
        cfg.threshold = 0.12; // refine only the steep parts of the wave
        let sol = solve_amr(&grid, &cfg, &init, 0.5).unwrap();
        let cells = sol.cells_per_level();
        assert!(cells[1] > 0 && cells[1] < 2 * n, "partial refinement: {cells:?}");
        let rho = sol.primitive_on_level(1);
        let m = 2 * n;
        (0..m).map(|i| (rho[i][0] - exact((i as f64 + 0.5) / m as f64, 0.5)).abs() / m as f64).sum::<f64>()
    };
    let (e1, e2) = (err(40), err(80));
    assert!((e1 / e2).log2() > 1.6, "AMR order {} ({e1}, {e2})", (e1 / e2).log2());
    // periodic mass conservation through the refinement interfaces
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let mut c = FvConfig::euler(1.4);
    c.cfl = 0.4;
    let mut cfg = AmrConfig::new(c.with_axis(0, GasBc::Periodic));
    cfg.max_level = 1;
    cfg.threshold = 0.12;
    let t0 = solve_amr(&grid, &cfg, &init, 1e-12).unwrap().totals();
    let t1 = solve_amr(&grid, &cfg, &init, 0.5).unwrap().totals();
    assert!((t1[0] - t0[0]).abs() < 1e-12 && (t1[4] - t0[4]).abs() < 1e-11 * t0[4], "{t0:?} {t1:?}");
    assert!(solve_amr(&FvGrid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3]).unwrap(), &cfg, &init, 0.1).is_err());
}
