//! Geometric multigrid tests (V, W, FMG; 1D/2D/3D; mixed boundary conditions).

use rssn::sim::boundary::{BcValue, BoundaryCondition as Bc, Face, FaceBcs};
use rssn::sim::multigrid::{MgCycle, MgOptions, MultigridSolver, solve_helmholtz_multigrid};
use rssn::sim::structured::{Grid, LinearSolver, solve_helmholtz};
use std::f64::consts::PI;

fn exact(dim: usize, x: [f64; 3]) -> f64 {
    let e = (x[0] + if dim > 1 { x[1] } else { 0.0 }).exp();
    if dim == 3 { e * (2.0 * PI * x[2]).cos() } else { e }
}

fn minus_lap(dim: usize, x: [f64; 3]) -> f64 {
    let c = match dim {
        1 => 1.0,
        2 => 2.0,
        _ => 2.0 - 4.0 * PI * PI,
    };
    -c * exact(dim, x)
}

fn mixed_bcs(dim: usize) -> FaceBcs {
    let u = move |x: [f64; 3]| exact(dim, x);
    if dim == 1 {
        return FaceBcs::dirichlet_zero()
            .with(Face::XMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
            .with(Face::XMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u(x)) });
    }
    let mut b = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::XMax, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::YMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
        .with(Face::YMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u(x)) });
    if dim == 3 {
        b = b.with_axis(2, &Bc::Periodic);
    }
    b
}

fn grid(dim: usize, n: usize, bcs: &FaceBcs) -> Grid {
    let mut ns = vec![n + 1; dim];
    if dim == 3 {
        ns[2] = n;
    }
    Grid::new(&ns, &vec![0.0; dim], &vec![1.0; dim], bcs).unwrap()
}

#[test]
fn mg_matches_direct_solution_with_mixed_bcs_in_all_dimensions() {
    for (dim, n) in [(1, 64), (2, 32), (3, 16)] {
        let bcs = mixed_bcs(dim);
        let g = grid(dim, n, &bcs);
        let f = |x: [f64; 3]| minus_lap(dim, x);
        let direct = solve_helmholtz(&g, &bcs, 0.0, f, LinearSolver::Direct).unwrap();
        for cycle in [MgCycle::V, MgCycle::W, MgCycle::Fmg] {
            let opts = MgOptions { cycle, tol: 1e-12, ..MgOptions::default() };
            let r = solve_helmholtz_multigrid(&g, &bcs, 0.0, f, opts).unwrap();
            assert!(r.converged, "dim {dim} {cycle:?}: {:?}", r.residuals);
            let diff = r.u.iter().zip(&direct).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(diff < 1e-8, "dim {dim} {cycle:?}: diff {diff}");
            assert!(g.max_error(&r.u, |x| exact(dim, x)) < 0.05);
        }
    }
}

#[test]
fn mg_convergence_factor_independent_of_mesh_size() {
    let bcs = FaceBcs::dirichlet_zero();
    let f = |x: [f64; 3]| 2.0 * PI * PI * (PI * x[0]).sin() * (PI * x[1]).sin();
    let mut factors = Vec::new();
    for n in [17, 33, 65, 129] {
        let g = Grid::unit(2, n, &bcs).unwrap();
        let opts = MgOptions { max_cycles: 8, tol: 1e-14, ..MgOptions::default() };
        let r = solve_helmholtz_multigrid(&g, &bcs, 0.0, f, opts).unwrap();
        factors.push(r.convergence_factor());
    }
    for &c in &factors {
        assert!(c < 0.2, "factors {factors:?}");
    }
    let spread = factors.iter().fold(0.0_f64, |m, v| m.max(*v)) - factors.iter().fold(1.0_f64, |m, v| m.min(*v));
    assert!(spread < 0.1, "factors {factors:?}");
}

#[test]
fn mg_w_cycle_not_worse_than_v_and_fmg_reaches_truncation_error_in_one_pass() {
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::unit(2, 65, &bcs).unwrap();
    let u_ex = |x: [f64; 3]| (PI * x[0]).sin() * (2.0 * PI * x[1]).sin();
    let f = |x: [f64; 3]| 5.0 * PI * PI * u_ex(x);
    let run = |cycle, max_cycles| {
        let opts = MgOptions { cycle, max_cycles, tol: 1e-14, ..MgOptions::default() };
        solve_helmholtz_multigrid(&g, &bcs, 0.0, f, opts).unwrap()
    };
    let v = run(MgCycle::V, 6);
    let w = run(MgCycle::W, 6);
    assert!(w.convergence_factor() <= v.convergence_factor() * 1.05);
    let fmg = run(MgCycle::Fmg, 0);
    let converged = run(MgCycle::V, 40);
    let e_fmg = g.max_error(&fmg.u, u_ex);
    let e_conv = g.max_error(&converged.u, u_ex);
    assert!(e_fmg < 1.5 * e_conv, "fmg {e_fmg} vs converged {e_conv}");
    assert!(fmg.residuals[1] < 0.05 * fmg.residuals[0]);
}

#[test]
fn mg_periodic_helmholtz_and_second_order() {
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let sigma = 3.0;
    let u_ex = |x: [f64; 3]| (2.0 * PI * x[0]).sin() * (4.0 * PI * x[1]).cos();
    let f = |x: [f64; 3]| (20.0 * PI * PI + sigma) * u_ex(x);
    let mut errs = Vec::new();
    for n in [16, 32, 64] {
        let g = Grid::new(&[n, n], &[0.0; 2], &[1.0; 2], &bcs).unwrap();
        let r = solve_helmholtz_multigrid(&g, &bcs, sigma, f, MgOptions::default()).unwrap();
        assert!(r.converged);
        errs.push(g.max_error(&r.u, u_ex));
    }
    assert!((errs[0] / errs[1]).log2() > 1.8 && (errs[1] / errs[2]).log2() > 1.8, "{errs:?}");
    let g = Grid::new(&[32, 32], &[0.0; 2], &[1.0; 2], &bcs).unwrap();
    assert!(MultigridSolver::new(&g, &bcs, sigma, MgOptions::default()).unwrap().num_levels() >= 3);
}

#[test]
fn mg_time_dependent_dirichlet_data() {
    let bcs = FaceBcs::uniform(&Bc::Dirichlet(BcValue::spacetime(|t, x| t * (x[0] + 2.0 * x[1]))));
    let g = Grid::unit(2, 33, &bcs).unwrap();
    let mg = MultigridSolver::new(&g, &bcs, 0.0, MgOptions::default()).unwrap();
    let f = vec![0.0; g.len()];
    let r = mg.solve(&f, 2.0).unwrap();
    assert!(r.converged);
    // harmonic linear function reproduced exactly
    assert!(g.max_error(&r.u, |x| 2.0 * (x[0] + 2.0 * x[1])) < 1e-8);
}
