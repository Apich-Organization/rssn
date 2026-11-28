//! Discrete divergence and volume-flux identities for staggered MAC grids.

use std::f64::consts::TAU;

use rssn::sim::boundary::FluidBoundary;
use rssn::sim::boundary::FluidFaceBcs;
use rssn::sim::ns_mac::Density;
use rssn::sim::ns_mac::ImmersedBody;
use rssn::sim::ns_mac::NsConfig;
use rssn::sim::ns_mac::NsSolver;

fn reference(
    solver: &NsSolver,
    dim: usize,
) -> Vec<f64> {
    let n = solver.n();
    let mut result = Vec::new();
    for z in 0..n[2] {
        for y in 0..n[1] {
            for x in 0..n[0] {
                let i = [x, y, z];
                let mut sum = 0.0;
                for c in 0..dim {
                    let mut ext = n;
                    ext[c] += 1;
                    let mut j = i;
                    j[c] += 1;
                    let flat = |[x, y, z]: [usize; 3]| x + ext[0] * (y + ext[1] * z);
                    let lower = solver.face_aperture(c, i) * solver.u[c][flat(i)];
                    let upper = solver.face_aperture(c, j) * solver.u[c][flat(j)];
                    sum += (upper - lower) / solver.cell_size(c);
                }
                result.push(sum);
            }
        }
    }
    result
}

fn assert_bits(
    a: &[f64],
    b: &[f64],
) {
    assert_eq!(a.len(), b.len());
    for (p, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "cell {p}");
    }
}

#[test]
fn affine_flow_matches_flux_divergence_and_norms_in_2d_and_3d() {
    let n = [5, 4, 3];
    let lengths = [2.0, 3.0, 4.0];
    for dim in 2..=3 {
        let mut cfg = NsConfig::new(&n[..dim], &lengths[..dim], 0.05).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::outflow());
        let solver = NsSolver::new(cfg, 0.02, &|[x, y, z]| {
            [0.25 + 2.0 * x, -0.75 - y, 0.5 + 0.5 * z]
        })
        .unwrap();
        let expected = if dim == 2 { 1.0 } else { 1.5 };
        let volume: f64 = lengths[..dim].iter().product();
        let mut into = vec![f64::NAN; n[..dim].iter().product()];
        solver.divergence_into(&mut into).unwrap();
        assert_bits(&into, &reference(&solver, dim));
        assert_bits(&into, &solver.divergence());
        assert!(into.iter().all(|x| (x - expected).abs() < 3.0e-15));
        let summary = solver.divergence_summary().unwrap();
        assert!((summary.max_abs - expected).abs() < 3.0e-15);
        assert!((summary.rms - expected).abs() < 3.0e-15);
        assert!((summary.fluid_volume - volume).abs() < 1.0e-14);
        let boundary: f64 = (0..dim)
            .map(|c| solver.flux_through(c, n[c]) - solver.flux_through(c, 0))
            .sum();
        assert!((summary.net_outflow - boundary).abs() < 3.0e-14);
        assert!((summary.net_outflow - expected * volume).abs() < 3.0e-14);
    }
}

#[test]
fn immersed_planar_cut_uses_apertures_and_actual_fluid_volume() {
    let n = [5, 4, 3];
    let lengths = [5.0, 4.0, 3.0];
    for dim in 2..=3 {
        let mut cfg = NsConfig::new(&n[..dim], &lengths[..dim], 0.05).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::outflow());
        cfg.body = Some(ImmersedBody {
            level_set: Box::new(|x| x[0] - 1.2),
            velocity: Box::new(|_| [0.0; 3]),
        });
        let solver = NsSolver::new(cfg, 0.02, &|x| [0.0, x[1], 0.0]).unwrap();
        let volume = 3.8 * lengths[1..dim].iter().product::<f64>();
        let values = solver.divergence();
        assert_bits(&values, &reference(&solver, dim));
        assert_eq!(values[0], 0.0);
        assert!((values[1] - 0.8).abs() < 1.0e-12);
        assert!((values[2] - 1.0).abs() < 1.0e-12);
        let summary = solver.divergence_summary().unwrap();
        assert!((summary.fluid_volume - volume).abs() < 1.0e-11);
        assert!((summary.net_outflow - volume).abs() < 1.0e-11);
        assert!((summary.max_abs - 1.0).abs() < 1.0e-11);
        assert!((summary.rms - 1.0).abs() < 1.0e-11);
    }
}

#[test]
fn staircase_solid_cell_contributes_no_fluid_volume() {
    let mut cfg = NsConfig::new(&[5, 5], &[5.0, 5.0], 0.1).unwrap();
    let mut mask = vec![false; 25];
    mask[12] = true;
    cfg.solid = Some(mask);
    let solver = NsSolver::new(cfg, 0.01, &|_| [0.0; 3]).unwrap();
    let summary = solver.divergence_summary().unwrap();
    assert_eq!(summary.fluid_volume, 24.0);
    assert_eq!(
        (summary.net_outflow, summary.rms, summary.max_abs),
        (0.0, 0.0, 0.0)
    );
}

#[test]
fn summary_keeps_large_representable_norms_and_checks_buffer_shapes() {
    let mut cfg = NsConfig::new(&[3, 3], &[3.0, 3.0], 0.1).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::outflow());
    let mut solver = NsSolver::new(cfg, 0.01, &|x| [1.0e200 * x[0], 0.0, 0.0]).unwrap();
    let mut wrong = [4.0; 8];
    assert!(solver.divergence_into(&mut wrong).is_err());
    assert_eq!(wrong, [4.0; 8]);
    let summary = solver.divergence_summary().unwrap();
    assert!((summary.rms / 1.0e200 - 1.0).abs() < 1.0e-15);
    assert!((summary.max_abs / 1.0e200 - 1.0).abs() < 1.0e-15);
    assert!((summary.net_outflow / 9.0e200 - 1.0).abs() < 1.0e-15);
    solver.u[0][1] = f64::NAN;
    assert!(solver.divergence_summary().is_err());
    solver.u[0].pop();
    assert!(solver.divergence_into(&mut [0.0; 9]).is_err());
}

#[test]
fn pressure_projection_removes_divergence_at_constant_and_variable_density() {
    for variable in [false, true] {
        let mut cfg = NsConfig::new(&[11, 9], &[1.0, 1.0], 0.005).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
        if variable {
            cfg.density = Some(Density {
                init: Box::new(|x| 2.0 + 0.5 * (TAU * x[1]).cos()),
                transport: true,
                gravity: [0.0; 3],
            });
        }
        let mut solver =
            NsSolver::new(cfg, 0.001, &|x| [0.25 * (TAU * x[0]).sin(), 0.0, 0.0]).unwrap();
        assert!(solver.divergence_summary().unwrap().rms > 1.0);
        for _ in 0..3 {
            solver.step().unwrap();
            let summary = solver.divergence_summary().unwrap();
            assert!(summary.max_abs < 2.0e-12, "{summary:?}");
            assert!(summary.net_outflow.abs() < 2.0e-15, "{summary:?}");
        }
    }
}

#[test]
fn row_blocks_keep_the_same_bits_across_worker_counts() {
    let mut cfg = NsConfig::new(&[137, 251], &[7.0, 9.0], 0.01).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
    let solver = NsSolver::new(cfg, 0.01, &|[x, y, _]| [x.sin() * y, y.sin() * x, 0.0]).unwrap();
    let run = |workers| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| (solver.divergence(), solver.divergence_summary().unwrap()))
    };
    let (a, sa) = run(1);
    let (b, sb) = run(2);
    assert_bits(&a, &b);
    assert_bits(&a, &reference(&solver, 2));
    assert_eq!(sa, sb);
}
