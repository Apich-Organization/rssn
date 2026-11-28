//! Continuum multigrid relative residuals and scratch isolation.

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::multigrid::MgCycle;
use rssn::sim::multigrid::MgOptions;
use rssn::sim::multigrid::MgResult;
use rssn::sim::multigrid::MultigridSolver;
use rssn::sim::structured::Grid;
use std::sync::Barrier;

#[test]
fn periodic_constant_solution_scales_across_norm_underflow_and_overflow() {
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let grid = Grid::new(&[16], &[0.0], &[1.0], &bcs).unwrap();
    for cycle in [MgCycle::V, MgCycle::W, MgCycle::Fmg] {
        let solver = MultigridSolver::new(
            &grid,
            &bcs,
            2.0,
            MgOptions {
                cycle,
                coarsest_nodes: 4,
                ..MgOptions::default()
            },
        )
        .unwrap();

        // The periodic Laplacian annihilates constants and the shift is two.
        for scale in [1e-200, 1.0, 1e200] {
            let result = solver.solve(&[2.0 * scale; 16], 0.0).unwrap();
            assert!((result.residuals[0] - 1.0).abs() <= 2e-15);
            assert!(result.residuals.iter().all(|r| r.is_finite()));
            assert!(result.converged, "{cycle:?} scale {scale:e}");
            assert!(result.u.iter().all(|u| (u / scale - 1.0).abs() <= 1e-8));
        }
    }
}

#[test]
fn zero_full_multigrid_residual_reports_a_finite_reduction_factor() {
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let grid = Grid::new(&[16], &[0.0], &[1.0], &bcs).unwrap();
    let solver = MultigridSolver::new(
        &grid,
        &bcs,
        2.0,
        MgOptions {
            cycle: MgCycle::Fmg,
            ..MgOptions::default()
        },
    )
    .unwrap();
    let result = solver.solve(&[0.0; 16], 0.0).unwrap();
    assert!(result.converged);
    assert_eq!(
        result
            .residuals
            .iter()
            .map(|r| r.to_bits())
            .collect::<Vec<_>>(),
        [0, 0]
    );
    assert_eq!(result.convergence_factor().to_bits(), 0);
}

#[test]
fn geometric_reduction_can_be_finite_when_the_endpoint_ratio_is_not() {
    for (start, end, expected) in [(1e200, 1e-200, 1e-200), (1e-200, 1e200, 1e200)] {
        let result = MgResult {
            u: Vec::new(),
            residuals: vec![start, 1.0, end],
            converged: false,
        };
        assert!((result.convergence_factor() / expected - 1.0).abs() < 5e-15);
    }
}

fn same_bits(
    a: &MgResult,
    b: &MgResult,
) {
    assert_eq!(a.converged, b.converged);
    assert_eq!(a.u.len(), b.u.len());
    assert_eq!(a.residuals.len(), b.residuals.len());
    for (x, y) in
        a.u.iter()
            .chain(&a.residuals)
            .zip(b.u.iter().chain(&b.residuals))
    {
        assert_eq!(x.to_bits(), y.to_bits());
    }
}

#[test]
fn a_shared_hierarchy_isolates_simultaneous_solves_and_time_dependent_boundaries() {
    let linear = |t: f64, x: [f64; 3]| t * (x[0] + 2.0 * x[1]);
    let bcs = FaceBcs::uniform(&Bc::Dirichlet(BcValue::spacetime(linear)));
    // Only y can coarsen. This exercises a strided line throughout the hierarchy.
    let grid = Grid::new(&[12, 17], &[0.0, 0.0], &[1.0, 1.0], &bcs).unwrap();
    let sigma = 0.5;
    let times = [0.25, -2.0, 1.5, 5.0];

    for cycle in [MgCycle::V, MgCycle::W, MgCycle::Fmg] {
        let solver = MultigridSolver::new(
            &grid,
            &bcs,
            sigma,
            MgOptions {
                cycle,
                coarsest_nodes: 30,
                ..MgOptions::default()
            },
        )
        .unwrap();
        assert!(solver.num_levels() > 1);
        let run = |t| {
            solver
                .solve(&grid.sample(|x| sigma * linear(t, x)), t)
                .unwrap()
        };
        let baseline = times.map(run);
        let barrier = Barrier::new(times.len());
        std::thread::scope(|scope| {
            let barrier = &barrier;
            let run = &run;
            let jobs = times.map(|t| {
                scope.spawn(move || {
                    barrier.wait();
                    run(t)
                })
            });
            for ((t, expected), job) in times.into_iter().zip(&baseline).zip(jobs) {
                let result = job.join().unwrap();
                same_bits(&result, expected);
                assert!(result.converged);
                assert!(grid.max_error(&result.u, |x| linear(t, x)) < 2e-8);
            }
        });
    }
}

#[test]
fn tensor_workspace_covers_levels_that_coarsen_different_axis_sets() {
    let affine = |x: [f64; 3]| 1.0 + x[0] + 2.0 * x[2];
    let bcs = FaceBcs::uniform(&Bc::Dirichlet(BcValue::space(affine))).with_axis(1, &Bc::Periodic);
    // Coarsens x/y/z, then x/y, then only x.
    let grid = Grid::new(&[17, 12, 7], &[0.0; 3], &[1.0; 3], &bcs).unwrap();
    let f = grid.sample(|x| 2.0 * affine(x));
    for cycle in [MgCycle::V, MgCycle::W, MgCycle::Fmg] {
        let solver = MultigridSolver::new(
            &grid,
            &bcs,
            2.0,
            MgOptions {
                cycle,
                coarsest_nodes: 16,
                ..MgOptions::default()
            },
        )
        .unwrap();
        assert_eq!(solver.num_levels(), 4);
        let result = solver.solve(&f, 0.0).unwrap();
        assert!(result.converged, "{cycle:?}");
        assert!(grid.max_error(&result.u, affine) < 1e-8);
    }
}

#[test]
fn rejects_invalid_active_data_and_leaves_dirichlet_source_slots_unused() {
    let bcs = FaceBcs::dirichlet_zero();
    let grid = Grid::new(&[9], &[0.0], &[1.0], &bcs).unwrap();
    assert!(MultigridSolver::new(&grid, &bcs, f64::NAN, MgOptions::default()).is_err());
    assert!(
        MultigridSolver::new(
            &grid,
            &bcs,
            0.0,
            MgOptions {
                tol: -1.0,
                ..MgOptions::default()
            }
        )
        .is_err()
    );
    let solver = MultigridSolver::new(&grid, &bcs, 0.0, MgOptions::default()).unwrap();
    assert!(solver.solve(&[0.0; 8], 0.0).is_err());
    assert!(solver.solve(&[0.0; 9], f64::INFINITY).is_err());
    let mut rhs = [2.0; 9];
    rhs[0] = f64::NAN;
    rhs[8] = f64::INFINITY;
    assert!(solver.solve(&rhs, 0.0).unwrap().converged);
    rhs[4] = f64::NAN;
    assert!(solver.solve(&rhs, 0.0).is_err());
}
