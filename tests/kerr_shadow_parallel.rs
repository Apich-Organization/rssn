//! Numerical Kerr shadow ordering, cancellation, and analytic boundary checks.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use rssn::sim::kerr::shadow_boundary_analytic;
use rssn::sim::kerr::trace_shadow;

#[test]
fn ray_tracing_preserves_angles_and_bits_across_worker_counts() {
    let trace = |workers| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| trace_shadow(1.0, -0.6, 1.2, 1.0e4, 12).unwrap())
    };
    let serial = trace(1);
    let parallel = trace(4);
    for (index, ((x, y), (parallel_x, parallel_y))) in serial.into_iter().zip(parallel).enumerate()
    {
        let angle = std::f64::consts::TAU * (index as f64 + 0.5) / 12.0;
        assert!((x * angle.sin() - y * angle.cos()).abs() < 1e-14);
        assert!(x * angle.cos() + y * angle.sin() > 0.0);
        assert_eq!(x.to_bits(), parallel_x.to_bits());
        assert_eq!(y.to_bits(), parallel_y.to_bits());
    }
}

#[test]
fn rays_observe_the_cancellation_scope_in_both_pool_modes() {
    for workers in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| {
            let cancellation = Arc::new(AtomicBool::new(true));
            let _guard = rssn::cancellation::install(cancellation);
            let error = trace_shadow(1.0, 0.0, 1.2, 1.0e4, 12).unwrap_err();
            assert_eq!(error, "computation cancelled");
        });
    }
    assert!(trace_shadow(1.0, 0.0, 1.2, 1.0e4, 0).unwrap().is_empty());
}

#[test]
fn analytic_shadow_checks_inputs_and_preserves_mirrored_order() {
    assert!(shadow_boundary_analytic(1.0, 0.0, f64::NAN, 4).is_err());
    assert!(shadow_boundary_analytic(f64::INFINITY, 0.0, 1.2, 4).is_err());
    assert!(shadow_boundary_analytic(1.0, 0.6, 1.2, usize::MAX).is_err());
    assert!(shadow_boundary_analytic(1.0, 0.6, 1.2, 1).is_err());
    assert!(shadow_boundary_analytic(1.0, 0.6, 0.0, 4).is_err());
    let mirrored = shadow_boundary_analytic(1.0, 0.6, 1.2, 65).unwrap();
    assert_eq!(mirrored.len(), 130);
    for (&(upper_x, upper_y), &(lower_x, lower_y)) in
        mirrored[..65].iter().zip(mirrored[65..].iter().rev())
    {
        assert_eq!(upper_x.to_bits(), lower_x.to_bits());
        assert_eq!((-upper_y).to_bits(), lower_y.to_bits());
    }
}
