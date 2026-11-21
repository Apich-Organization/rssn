//! Exterior acoustic fields after restarted complex GMRES.

use rssn::sim::bem_helmholtz::HelmholtzAcaOptions;
use rssn::sim::bem_helmholtz::PlaneWave;
use rssn::sim::bem_helmholtz::Scatterer;
use rssn::sim::bem_helmholtz::solve_helmholtz_exterior;
use rssn::sim::bem_helmholtz::solve_helmholtz_exterior_aca;
use rssn::sim::physics_bem::SurfaceMesh3D;

#[test]
fn compressed_and_dense_exterior_fields_agree_after_restarts() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let wave = PlaneWave {
        direction: [2.0, -1.0, 4.0],
    };
    let opts = HelmholtzAcaOptions {
        restart: 5,
        max_iter: 300,
        gmres_tol: 1e-10,
        leaf_size: 8,
        ..HelmholtzAcaOptions::default()
    };
    let field_points = [[1.0, -0.4, 3.0], [-2.4, 1.3, -0.25], [0.2, 2.2, -1.8]];
    for kind in [Scatterer::SoundSoft, Scatterer::SoundHard] {
        let dense = solve_helmholtz_exterior(&mesh, 1.7, kind, wave, None).unwrap();
        let (iterative, stats) =
            solve_helmholtz_exterior_aca(&mesh, 1.7, kind, wave, None, &opts).unwrap();
        assert!(stats.iterations > opts.restart, "{kind:?}: {stats:?}");
        for p in field_points {
            let reference = dense.scattered_field(&mesh, p);
            let field = iterative.scattered_field(&mesh, p);
            assert!(
                (reference - field).norm() < 2e-8,
                "{kind:?} {p:?}: {field} vs {reference}"
            );
        }
    }
}

#[test]
fn exterior_solver_rejects_invalid_restart_and_tolerance() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 0);
    let wave = PlaneWave {
        direction: [0.0, 0.0, 1.0],
    };
    for (restart, gmres_tol) in [(0, 1e-12), (8, f64::INFINITY), (8, -1e-12)] {
        let opts = HelmholtzAcaOptions {
            restart,
            gmres_tol,
            ..HelmholtzAcaOptions::default()
        };
        let error =
            solve_helmholtz_exterior_aca(&mesh, 1.5, Scatterer::SoundSoft, wave, None, &opts)
                .unwrap_err();
        assert!(error.contains("GMRES"));
    }
}
