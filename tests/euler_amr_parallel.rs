//! Refluxed AMR conservation, deterministic flux sweeps, and checked time spans.

use rssn::sim::euler_amr::AmrConfig;
use rssn::sim::euler_amr::solve_amr;
use rssn::sim::euler_nd::FvConfig;
use rssn::sim::euler_nd::FvGrid;
use rssn::sim::euler_nd::GasBc;
use rssn::sim::euler_nd::Reconstruction;
use rssn::sim::euler_nd::RiemannSolver;

#[test]
fn refluxed_weno_sweeps_match_between_pools_on_a_rectangular_grid() {
    let grid = FvGrid::new(&[32, 24, 16], &[0.0; 3], &[1.0; 3]).unwrap();
    let mut fv = FvConfig::euler(1.4)
        .with_axis(0, GasBc::Periodic)
        .with_axis(1, GasBc::Periodic)
        .with_axis(2, GasBc::Periodic);
    fv.recon = Reconstruction::Weno5;
    fv.riemann = RiemannSolver::Hlle;
    let mut config = AmrConfig::new(fv);
    config.max_level = 1;
    config.threshold = 0.055;
    let initial = |[x, y, z]: [f64; 3]| {
        let r = (x - 0.3).powi(2) + (y - 0.5).powi(2) + (z - 0.7).powi(2);
        [1.0 + 0.7 * (-r / 0.012).exp(), 0.45, -0.3, 0.2, 1.0]
    };
    let run = |workers, time| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(|| solve_amr(&grid, &config, &initial, time).unwrap())
    };
    let initial = run(1, 0.0);
    let end = 0.017;
    let single = run(1, end);
    let multiple = run(2, end);
    assert!(single.cells_per_level()[1] > 0);
    assert_eq!(single.t.to_bits(), end.to_bits());
    assert_eq!(multiple.t.to_bits(), end.to_bits());
    assert_eq!(single.steps, multiple.steps);
    assert_eq!(single.levels.len(), multiple.levels.len());
    for (serial, parallel) in single.levels.iter().zip(&multiple.levels) {
        assert_eq!(serial.patches.len(), parallel.patches.len());
        for ((lo, n, q), (other_lo, other_n, other_q)) in
            serial.patches.iter().zip(&parallel.patches)
        {
            assert_eq!((lo, n), (other_lo, other_n));
            assert_eq!(q.len(), other_q.len());
            for (a, b) in q.iter().zip(other_q) {
                assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
            }
        }
    }
    for (start, finish) in initial.totals().into_iter().zip(single.totals()) {
        assert!((start - finish).abs() <= 5e-12 * start.abs().max(1.0));
    }
}

#[test]
fn invalid_spans_and_refinement_are_rejected_before_sampling() {
    let grid = FvGrid::new(&[4], &[0.0], &[1.0]).unwrap();
    let mut config = AmrConfig::new(FvConfig::euler(1.4));
    let unused = |_| panic!("invalid request must fail before sampling initial data");
    for time in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(solve_amr(&grid, &config, &unused, time).is_err());
    }
    config.max_level = usize::BITS as usize;
    assert!(solve_amr(&grid, &config, &unused, 0.0).is_err());
    config.max_level = 0;
    let invalid_grid = FvGrid { dim: 4, ..grid };
    assert!(solve_amr(&invalid_grid, &config, &unused, 0.1).is_err());

    let end = f64::from_bits(1);
    let uniform = |_| [1.0, 0.0, 0.0, 0.0, 1.0];
    let result = solve_amr(&grid, &config, &uniform, end).unwrap();
    assert_eq!(result.steps, 1);
    assert_eq!(result.t.to_bits(), end.to_bits());
}
