//! Determinism and timing checks of the parallel simulation sites (Part 6 of
//! `docs/ALGORITHM_COMPLETION_CHECKLIST.md`). Each site runs in a one-thread
//! rayon pool and in the global pool; the results must be equal (Debug text,
//! so floats are compared bit-exactly). Timings are printed, never asserted:
//! run with `--nocapture` and grep `PAR`.

use std::time::Instant;

fn compare<T: std::fmt::Debug + Send>(
    site: &str,
    size: &str,
    f: impl Fn() -> T + Send + Sync,
) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let t0 = Instant::now();
    let a = pool.install(&f);
    let ts = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let b = f();
    let tp = t1.elapsed().as_secs_f64();
    assert_eq!(
        format!("{a:?}"),
        format!("{b:?}"),
        "{site}: parallel result differs from serial"
    );
    println!(
        "PAR | {site} | {size} | serial {:.1} ms | parallel {:.1} ms | x{:.2} | threads {}",
        ts * 1e3,
        tp * 1e3,
        ts / tp.max(1e-12),
        rayon::current_num_threads()
    );
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

#[test]
fn par_dg3d() {
    use rssn::sim::dg3d::Dg3Bc;
    use rssn::sim::dg3d::Dg3Limiter;
    use rssn::sim::dg3d::Dg3Mesh;
    use rssn::sim::dg3d::Dg3Options;
    use rssn::sim::dg3d::LinearAdvection3d;
    use rssn::sim::dg3d::solve_dg3d;
    use std::f64::consts::PI;
    let law = LinearAdvection3d { a: [1.0, 0.5, -0.5] };
    let per: [Dg3Bc; 6] = std::array::from_fn(|_| Dg3Bc::Periodic);
    let mesh = Dg3Mesh {
        ne: [8; 3],
        lo: [0.0; 3],
        hi: [1.0; 3],
        p: 2,
    };
    let init = |x: [f64; 3]| {
        vec![(2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos() * (2.0 * PI * x[2]).cos()]
    };
    for limiter in [Dg3Limiter::None, Dg3Limiter::Hierarchical] {
        let opt = Dg3Options {
            limiter,
            ..Dg3Options::default()
        };
        compare(
            "P11 solve_dg3d",
            &format!("8^3 elements p=2, {limiter:?}"),
            || {
                solve_dg3d(&mesh, &law, &per, &init, 0.05, 0.3, &opt)
                    .unwrap()
                    .q
            },
        );
    }
}

#[test]
fn par_fem_simplex() {
    use rssn::sim::boundary::BcValue;
    use rssn::sim::boundary::BoundaryCondition as Bc;
    use rssn::sim::boundary::MarkerBcs;
    use rssn::sim::fem_simplex::FemProblem;
    use rssn::sim::fem_simplex::FemSpace;
    use rssn::sim::fem_simplex::SimplexMesh;
    use rssn::sim::fem_simplex::error_indicators;
    use rssn::sim::fem_simplex::marker;
    use rssn::sim::fem_simplex::solve_static;
    use rssn::sim::structured::LinearSolver;
    let space = FemSpace::new(SimplexMesh::rectangle(100, 100, [0.0; 2], [1.0; 2]), 2).unwrap();
    let mut bcs = MarkerBcs::new();
    for m in [marker::XMIN, marker::XMAX, marker::YMIN, marker::YMAX] {
        bcs = bcs.set(m, Bc::Dirichlet(BcValue::constant(0.0)));
    }
    let f = |_t: f64, x: [f64; 3]| {
        2.0 * std::f64::consts::PI.powi(2)
            * (std::f64::consts::PI * x[0]).sin()
            * (std::f64::consts::PI * x[1]).sin()
    };
    let prob = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &f,
        bcs: &bcs,
    };
    compare(
        "P12 solve_static (assembly + solve)",
        "100x100 triangles, P2",
        || solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap(),
    );
    let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
    compare("P12 error_indicators", "100x100 triangles, P2", || {
        error_indicators(&space, &u, &prob, 0.0)
    });
}

#[test]
fn par_euler_nd_and_ns_mac() {
    use rssn::sim::boundary::Face;
    use rssn::sim::boundary::FluidBoundary;
    use rssn::sim::boundary::FluidFaceBcs;
    use rssn::sim::euler_nd::FvConfig;
    use rssn::sim::euler_nd::FvGrid;
    use rssn::sim::euler_nd::solve;
    use rssn::sim::ns_mac::NsConfig;
    use rssn::sim::ns_mac::NsSolver;
    let grid = FvGrid::new(&[192, 192], &[0.0, 0.0], &[1.0, 1.0]).unwrap();
    let cfg = FvConfig::euler(1.4);
    let init = |x: [f64; 3]| {
        let r2 = (x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2);
        if r2 < 0.04 {
            [1.0, 0.0, 0.0, 0.0, 1.0]
        } else {
            [0.125, 0.0, 0.0, 0.0, 0.1]
        }
    };
    compare("P13 euler_nd::solve", "192x192 2-D blast, t=0.02", || {
        solve(&grid, &cfg, &init, 0.02).unwrap().totals()
    });
    compare(
        "P13 ns_mac cavity (only Csr::matvec parallel; stencils sequential)",
        "96x96, 15 steps",
        || {
            let mut cfg = NsConfig::new(&[96, 96], &[1.0, 1.0], 0.01).unwrap();
            cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
                .with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
            let mut s = NsSolver::new(cfg, 0.005, &|_| [0.0; 3]).unwrap();
            s.run(15).unwrap();
            (s.kinetic_energy(), s.max_speed())
        },
    );
}

#[test]
fn par_fmm_nbody_ewald() {
    use rssn::sim::ewald::ewald_energy_forces;
    use rssn::sim::ewald::ewald_reciprocal;
    use rssn::sim::fmm::FmmOptions;
    use rssn::sim::fmm::fmm3d;
    use rssn::sim::nbody::direct_field;
    use rssn::sim::nbody::gravity_accelerations;
    let mut st = 12345_u64;
    let n = 30_000;
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|_| [lcg(&mut st), lcg(&mut st), lcg(&mut st)])
        .collect();
    let src: Vec<f64> = (0..n).map(|_| lcg(&mut st) - 0.5).collect();
    let opts = FmmOptions {
        order: 4,
        leaf_size: 40,
        theta: 0.5,
        eps: 0.0,
    };
    compare("P14 fmm3d", "30000 points, order 4", || {
        let r = fmm3d(&pos, &src, &opts).unwrap();
        (r.potential, r.field)
    });
    let m = 3000;
    let p2: Vec<[f64; 3]> = pos[..m].to_vec();
    let s2: Vec<f64> = src[..m].to_vec();
    compare("P18 direct_field", "3000 bodies O(N^2)", || {
        direct_field(&p2, &s2, 1e-3)
    });
    let mass: Vec<f64> = src.iter().map(|v| v.abs() + 0.1).collect();
    compare(
        "P18 gravity_accelerations (Barnes-Hut)",
        "30000 bodies",
        || gravity_accelerations(&pos, &mass, 1.0, 1e-3, 0.5).unwrap(),
    );
    // charges: neutral set in a box of length 10
    let nq = 2500;
    let len = [10.0; 3];
    let pq: Vec<[f64; 3]> = (0..nq)
        .map(|_| {
            [
                lcg(&mut st) * 10.0,
                lcg(&mut st) * 10.0,
                lcg(&mut st) * 10.0,
            ]
        })
        .collect();
    let q: Vec<f64> = (0..nq)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    compare("P17 ewald_reciprocal", "2500 charges, kmax 8", || {
        let r = ewald_reciprocal(&pq, &q, len, 0.8, 8, 1.0);
        (r.energy, r.forces, r.virial)
    });
    compare(
        "P17 ewald_energy_forces (real space)",
        "2500 charges",
        || ewald_energy_forces(&pq[..2500], &q, len, 0.8, 2, 4.9, 1.0),
    );
}
