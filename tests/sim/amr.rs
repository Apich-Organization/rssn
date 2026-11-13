//! Block-structured AMR for the Euler equations.

use rssn::sim::euler_amr::AmrConfig;
use rssn::sim::euler_amr::AmrSolution;
use rssn::sim::euler_amr::solve_amr;
use rssn::sim::euler_nd::FvConfig;
use rssn::sim::euler_nd::FvGrid;
use rssn::sim::euler_nd::GasBc;
use rssn::sim::euler_nd::Reconstruction;
use rssn::sim::euler_nd::SlopeLimiter;
use rssn::sim::euler_nd::solve;
use rssn::sim::gas_dynamics::Primitive;
use rssn::sim::gas_dynamics::exact_riemann;
use std::f64::consts::PI;

const GAMMA: f64 = 1.4;

fn sod(x: [f64; 3]) -> [f64; 5] {
    if x[0] < 0.5 {
        [1.0, 0.0, 0.0, 0.0, 1.0]
    } else {
        [0.125, 0.0, 0.0, 0.0, 0.1]
    }
}

fn sod_exact_density(
    x: f64,
    t: f64,
) -> f64 {
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
    exact_riemann(l, r, GAMMA, (x - 0.5) / t).rho
}

fn scheme() -> FvConfig {
    let mut c = FvConfig::euler(GAMMA);
    c.recon = Reconstruction::Muscl(SlopeLimiter::VanLeer);
    c.cfl = 0.4;
    c.with_axis(0, GasBc::Reflective)
}

/// 2D blast of radius 0.15 in a reflective box.
fn blast2(x: [f64; 3]) -> [f64; 5] {
    let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt();
    if r < 0.15 {
        [1.0, 0.0, 0.0, 0.0, 10.0]
    } else {
        [1.0, 0.0, 0.0, 0.0, 0.1]
    }
}

/// Density blob advected at `u = 1` in a periodic box (an exact solution of Euler).
fn blob(
    x: [f64; 3],
    t: f64,
) -> f64 {
    let r2 = (x[0] - 0.3 - t).powi(2) + (x[1] - 0.5).powi(2) + (x[2] - 0.5).powi(2);
    1.0 + 0.5 * (-r2 / 0.01).exp()
}

fn periodic3() -> FvConfig {
    let mut c = FvConfig::euler(GAMMA);
    c.recon = Reconstruction::Muscl(SlopeLimiter::Minmod);
    c.cfl = 0.4;
    c.with_axis(0, GasBc::Periodic)
        .with_axis(1, GasBc::Periodic)
        .with_axis(2, GasBc::Periodic)
}

/// FNV-1a over every conserved value of the hierarchy, level by level.
fn checksum(sol: &AmrSolution) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for lv in &sol.levels {
        for (_, _, q) in &lv.patches {
            for s in q {
                for &v in s {
                    h = (h ^ v.to_bits()).wrapping_mul(0x100_0000_01b3);
                }
            }
        }
    }
    h
}

fn density_l1_error(
    rho: &[[f64; 5]],
    n: usize,
) -> f64 {
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            (rho[i][0] - sod_exact_density(x, 0.2)).abs() / n as f64
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
    assert!(
        cells[1] > 0 && cells[2] > 0,
        "refinement was created: {cells:?}"
    );
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
    assert!(
        amr_err < 0.5 * coarse_err,
        "AMR {amr_err} vs coarse {coarse_err}"
    );
    assert!(amr_err < 1.5 * fine_err, "AMR {amr_err} vs fine {fine_err}");
    let total_cells: usize = cells.iter().sum();
    assert!(
        total_cells < fine_n,
        "AMR uses {total_cells} cells against {fine_n}"
    );
    // mass and energy are conserved to round-off (the waves do not reach the walls)
    let t = sol.totals();
    let (m0, e0) = (0.5 * 1.0 + 0.5 * 0.125, 0.5 * 1.0 / 0.4 + 0.5 * 0.1 / 0.4);
    assert!(
        (t[0] - m0).abs() < 1e-12 && (t[4] - e0).abs() < 1e-12,
        "totals {t:?}"
    );
}

#[test]
fn amr_sod_l1_beats_a_uniform_grid_with_the_same_number_of_cells() {
    // the adapted run is compared with a uniform grid holding exactly the
    // number of cells the hierarchy uses (levels summed)
    let grid = FvGrid::new(&[64], &[0.0], &[1.0]).unwrap();
    let mut cfg = AmrConfig::new(scheme());
    cfg.max_level = 2;
    cfg.threshold = 0.03;
    cfg.buffer = 1;
    let sol = solve_amr(&grid, &cfg, &sod, 0.2).unwrap();
    let cells: usize = sol.cells_per_level().iter().sum();
    let amr_err = sol.composite_l1_error(0, |x| sod_exact_density(x[0], 0.2));
    let uniform = solve(
        &FvGrid::new(&[cells], &[0.0], &[1.0]).unwrap(),
        &scheme(),
        &sod,
        0.2,
    )
    .unwrap();
    let uni_err = uniform.l1_error(0, |x| sod_exact_density(x[0], 0.2));
    println!(
        "Sod L1 at {cells} cells: adapted {amr_err:.6} ({:?}), uniform {uni_err:.6}",
        sol.cells_per_level()
    );
    assert!(
        amr_err < 0.95 * uni_err,
        "adapted {amr_err} vs uniform {uni_err} at {cells} cells"
    );
}

#[test]
fn amr_3d_conserves_mass_momentum_and_energy_only_with_refluxing() {
    // periodic box, hierarchy frozen at t = 0: the only transfer between the
    // levels is the flux across the coarse-fine interfaces, which the registers
    // correct. Without them the mismatch of the coarse and the fine fluxes
    // shows up directly in the totals.
    let grid = FvGrid::new(&[12, 12, 12], &[0.0; 3], &[1.0; 3]).unwrap();
    let init = |x: [f64; 3]| [blob(x, 0.0), 1.0, 0.0, 0.0, 1.0];
    let mut cfg = AmrConfig::new(periodic3());
    cfg.max_level = 2;
    cfg.threshold = 0.05;
    cfg.buffer = 1;
    cfg.regrid_interval = 1_000_000;
    let start = solve_amr(&grid, &cfg, &init, 0.0).unwrap();
    let cells = start.cells_per_level();
    assert!(cells[1] > 0 && cells[2] > 0, "3D refinement: {cells:?}");
    let with = solve_amr(&grid, &cfg, &init, 0.05).unwrap();
    cfg.reflux = false;
    let without = solve_amr(&grid, &cfg, &init, 0.05).unwrap();
    let (t0, t1, t2) = (start.totals(), with.totals(), without.totals());
    let scale = t0[0].abs();
    println!("cells {cells:?}, initial {t0:?}");
    println!("reflux on {t1:?}");
    println!("reflux off {t2:?}");
    // the flow is along x: mass, x-momentum and energy are carried, the two
    // transverse momenta stay zero by symmetry
    for c in [0, 1, 4] {
        let on = (t1[c] - t0[c]).abs();
        let off = (t2[c] - t0[c]).abs();
        println!("component {c}: |d| reflux on {on:.3e}, off {off:.3e}");
        assert!(
            on <= 1e-12 * scale,
            "component {c} with refluxing: {on:.3e} (scale {scale:.3e})"
        );
        assert!(
            off > 1e3 * on,
            "component {c} without refluxing: {off:.3e} against {on:.3e}"
        );
        assert!(
            off > 1e-8 * scale,
            "component {c} without refluxing must lose conservation: {off:.3e}"
        );
    }
    for c in [2, 3] {
        assert!(
            (t1[c] - t0[c]).abs() <= 1e-12 * scale && (t2[c] - t0[c]).abs() <= 1e-12 * scale,
            "transverse momentum {c} must stay zero: {} {}",
            t1[c],
            t2[c]
        );
    }
}

#[test]
fn amr_3d_advected_blob_resolves_across_refinement_boundaries() {
    let grid = FvGrid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3]).unwrap();
    let init = |x: [f64; 3]| [blob(x, 0.0), 1.0, 0.0, 0.0, 1.0];
    let mut cfg = AmrConfig::new(periodic3());
    cfg.max_level = 2;
    cfg.threshold = 0.02;
    cfg.buffer = 1;
    cfg.regrid_interval = 2;
    let sol = solve_amr(&grid, &cfg, &init, 0.1).unwrap();
    let cells = sol.cells_per_level();
    assert!(cells[1] > 0 && cells[2] > 0, "3D refinement: {cells:?}");
    assert!(
        cells[2] < 32 * 32 * 32,
        "level 2 covers only part of the box: {cells:?}"
    );
    let amr_err = sol.composite_l1_error(0, |x| blob(x, 0.1));
    let coarse = solve(&grid, &cfg.fv, &init, 0.1).unwrap();
    let coarse_err = coarse.l1_error(0, |x| blob(x, 0.1));
    println!("3D blob L1: adapted {amr_err:.6} ({cells:?}), coarse {coarse_err:.6}");
    assert!(
        amr_err < 0.5 * coarse_err,
        "adapted {amr_err} vs coarse {coarse_err}"
    );
    // the hierarchy keeps the totals through repeated regridding
    let t0 = solve_amr(&grid, &cfg, &init, 0.0).unwrap().totals();
    let t1 = sol.totals();
    for c in 0..5 {
        assert!(
            (t1[c] - t0[c]).abs() < 1e-11 * t0[0],
            "component {c}: {} vs {}",
            t1[c],
            t0[c]
        );
    }
}

/// The 1D and 2D results of the pre-3D module, bit for bit.
///
/// The checksums were taken from the 1D/2D-only implementation; the third axis
/// contributes no flux, no ghost cell and a factor-one volume, so every value
/// of a 1D or 2D run has to stay identical.
#[test]
fn amr_1d_and_2d_results_are_bit_identical_to_the_1d_2d_only_implementation() {
    let g1 = FvGrid::new(&[64], &[0.0], &[1.0]).unwrap();
    let g2 = FvGrid::new(&[32, 32], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
    let mut c2 = FvConfig::euler(GAMMA);
    c2.recon = Reconstruction::Muscl(SlopeLimiter::Minmod);
    c2.cfl = 0.4;
    let c2 = c2
        .with_axis(0, GasBc::Reflective)
        .with_axis(1, GasBc::Reflective);
    let check = |tag: &str, sol: &AmrSolution, sum: u64, cells: &[usize]| {
        assert_eq!(sol.cells_per_level(), cells, "{tag}: cells");
        assert_eq!(checksum(sol), sum, "{tag}: states");
    };
    let mut c1 = AmrConfig::new(scheme());
    c1.threshold = 0.03;
    c1.max_level = 0;
    check(
        "1d uniform",
        &solve_amr(&g1, &c1, &sod, 0.2).unwrap(),
        0xdd48_ecc4_1323_ca88,
        &[64],
    );
    c1.max_level = 1;
    check(
        "1d one level",
        &solve_amr(&g1, &c1, &sod, 0.2).unwrap(),
        0x516f_75e3_6dec_29c4,
        &[64, 86],
    );
    c1.max_level = 2;
    check(
        "1d two levels",
        &solve_amr(&g1, &c1, &sod, 0.2).unwrap(),
        0x53fc_336a_32be_d4eb,
        &[64, 86, 66],
    );
    let mut c = AmrConfig::new(c2);
    c.threshold = 0.25;
    c.buffer = 1;
    c.regrid_interval = 2;
    c.max_level = 0;
    check(
        "2d uniform",
        &solve_amr(&g2, &c, &blast2, 0.05).unwrap(),
        0x8e84_2c5e_3683_9ead,
        &[1024],
    );
    c.max_level = 1;
    check(
        "2d one level",
        &solve_amr(&g2, &c, &blast2, 0.05).unwrap(),
        0x3aa6_3087_f18f_2de1,
        &[1024, 1936],
    );
    c.max_level = 2;
    check(
        "2d two levels",
        &solve_amr(&g2, &c, &blast2, 0.05).unwrap(),
        0xe32b_5171_38f3_bd74,
        &[1024, 1936, 4340],
    );
}

#[test]
fn amr_conserves_with_refluxing_in_a_closed_2d_box() {
    // blast wave in a reflective box: refinement follows the front; totals stay exact
    let grid = FvGrid::new(&[32, 32], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
    let mut c = FvConfig::euler(GAMMA);
    c.recon = Reconstruction::Muscl(SlopeLimiter::Minmod);
    c.cfl = 0.4;
    let c = c
        .with_axis(0, GasBc::Reflective)
        .with_axis(1, GasBc::Reflective);
    let mut cfg = AmrConfig::new(c);
    cfg.max_level = 2;
    cfg.threshold = 0.25;
    cfg.regrid_interval = 2;
    cfg.buffer = 1;
    let sol = solve_amr(&grid, &cfg, &blast2, 0.12).unwrap();
    let cells = sol.cells_per_level();
    assert!(cells[1] > 0, "{cells:?}");
    assert!(
        cells[2] < 32 * 32 * 16 * 3 / 4,
        "level 2 covers only part of the box: {cells:?}"
    );
    // initial totals from the fine-sampled initial data on the finest patches are close; compare against a restart at t = 0
    let s0 = solve_amr(&grid, &cfg, &blast2, 1e-12).unwrap().totals();
    let s1 = sol.totals();
    for c in [0, 4] {
        assert!(
            (s1[c] - s0[c]).abs() < 1e-11 * s0[c].abs(),
            "component {c}: {} vs {}",
            s1[c],
            s0[c]
        );
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
    assert!(asym < 0.25, "x/y asymmetry {asym}");
    let uniform = {
        let g = FvGrid::new(&[128, 128], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
        solve(&g, &cfg.fv, &blast2, 0.12).unwrap()
    };
    let mut l1 = 0.0;
    for p in 0..128 * 128 {
        l1 += (uniform.primitive(p)[0] - rho[p][0]).abs() / (128.0 * 128.0);
    }
    // the same comparison for the plain 32 x 32 grid prolonged to the fine grid
    let coarse = solve(&grid, &cfg.fv, &blast2, 0.12).unwrap();
    let mut l1_coarse = 0.0;
    for j in 0..128 {
        for i in 0..128 {
            l1_coarse += (uniform.primitive(i + 128 * j)[0]
                - coarse.primitive(i / 4 + 32 * (j / 4))[0])
                .abs()
                / (128.0 * 128.0);
        }
    }
    println!("blast L1 vs fine uniform: amr {l1:.4}, coarse {l1_coarse:.4}");
    assert!(
        l1 < 0.6 * l1_coarse,
        "AMR {l1} should be much closer to the fine solution than the coarse grid {l1_coarse}"
    );
}

#[test]
fn amr_periodic_density_wave_keeps_second_order_accuracy_across_refinement_boundaries() {
    let exact = |x: f64, t: f64| 1.0 + 0.2 * (2.0 * PI * (x - t)).sin();
    let init = |x: [f64; 3]| [exact(x[0], 0.0), 1.0, 0.0, 0.0, 1.0];
    let err = |n: usize| {
        let grid = FvGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let mut c = FvConfig::euler(GAMMA);
        c.recon = Reconstruction::Muscl(SlopeLimiter::Mc);
        c.cfl = 0.4;
        let c = c.with_axis(0, GasBc::Periodic);
        let mut cfg = AmrConfig::new(c);
        cfg.max_level = 1;
        cfg.threshold = 1.1 / n as f64; // refine only the steep parts of the wave
        cfg.buffer = 1;
        let sol = solve_amr(&grid, &cfg, &init, 0.5).unwrap();
        let cells = sol.cells_per_level();
        assert!(
            cells[1] > 0 && cells[1] < 2 * n,
            "partial refinement: {cells:?}"
        );
        sol.composite_l1_error(0, |x| exact(x[0], 0.5))
    };
    let (e1, e2) = (err(40), err(80));
    assert!(
        (e1 / e2).log2() > 1.6,
        "AMR order {} ({e1}, {e2})",
        (e1 / e2).log2()
    );
    // periodic mass conservation through the refinement interfaces
    let grid = FvGrid::new(&[40], &[0.0], &[1.0]).unwrap();
    let mut c = FvConfig::euler(GAMMA);
    c.cfl = 0.4;
    let mut cfg = AmrConfig::new(c.with_axis(0, GasBc::Periodic));
    cfg.max_level = 1;
    cfg.threshold = 1.1 / 40.0;
    cfg.buffer = 1;
    let t0 = solve_amr(&grid, &cfg, &init, 1e-12).unwrap().totals();
    let t1 = solve_amr(&grid, &cfg, &init, 0.5).unwrap().totals();
    assert!(
        (t1[0] - t0[0]).abs() < 1e-12 && (t1[4] - t0[4]).abs() < 1e-11 * t0[4],
        "{t0:?} {t1:?}"
    );
    // 3D shallow water (barotropic) is accepted: a lake at rest stays at rest
    let mut sw = AmrConfig::new(FvConfig::shallow_water(9.81));
    sw.max_level = 1;
    let rest = solve_amr(
        &FvGrid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3]).unwrap(),
        &sw,
        &|_| [1.0, 0.0, 0.0, 0.0, 0.0],
        0.05,
    )
    .unwrap();
    assert!((rest.totals()[0] - 1.0).abs() < 1e-12);
}

#[test]
fn amr_shallow_water_3d_refines_and_conserves_mass_and_momentum() {
    let grid = FvGrid::new(&[12, 12, 12], &[0.0; 3], &[1.0; 3]).unwrap();
    let init = |x: [f64; 3]| {
        let r2 = (x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2) + (x[2] - 0.5).powi(2);
        [1.0 + 0.5 * (-r2 / 0.02).exp(), 0.3, 0.0, -0.2, 0.0]
    };
    let fv = FvConfig::shallow_water(9.81)
        .with_axis(0, GasBc::Periodic)
        .with_axis(1, GasBc::Periodic)
        .with_axis(2, GasBc::Periodic);
    let mut cfg = AmrConfig::new(fv);
    cfg.max_level = 1;
    cfg.buffer = 1;
    let a = solve_amr(&grid, &cfg, &init, 1e-12).unwrap().totals();
    let sol = solve_amr(&grid, &cfg, &init, 0.05).unwrap();
    let b = sol.totals();
    assert!(sol.cells_per_level()[1] > 0, "the blob must be refined");
    for c in 0..4 {
        assert!(
            (a[c] - b[c]).abs() < 1e-11 * (1.0 + a[c].abs()),
            "slot {c}: {} -> {}",
            a[c],
            b[c]
        );
    }
    assert!(b[3].abs() > 1e-3, "the vertical momentum must be carried");
}
