//! Regression tests for the iterative ACA boundary operator.

use rssn::sim::bem_aca::AcaOptions;
use rssn::sim::bem_aca::solve_laplace_bem_3d_aca;
use rssn::sim::bem_aca::solve_laplace_bem_3d_p1;
use rssn::sim::bem_aca::solve_laplace_bem_3d_p1_aca;
use rssn::sim::physics_bem::BoundaryCondition;
use rssn::sim::physics_bem::SurfaceMesh3D;
use rssn::sim::physics_bem::evaluate_potential_3d;
use rssn::sim::physics_bem::solve_laplace_bem_3d_mesh;

#[test]
fn mixed_p0_matches_dense_and_retains_relative_accuracy_at_extreme_scales() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let bcs = |scale| {
        (0..mesh.triangles.len())
            .map(|i| {
                let point = mesh.centroid(i).unwrap();
                if point[2] < 0.0 {
                    BoundaryCondition::Potential(scale * point[0])
                } else {
                    BoundaryCondition::Flux(scale * mesh.normal(i).unwrap()[0])
                }
            })
            .collect::<Vec<_>>()
    };
    let dense = solve_laplace_bem_3d_mesh(&mesh, &bcs(1.0)).unwrap();
    // Keep the same discrete system as the dense reference across all three scales.
    let opts = AcaOptions {
        leaf_size: 4,
        eta: 0.0,
        ..AcaOptions::default()
    };
    for scale in [2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
        let (compressed, stats) = solve_laplace_bem_3d_aca(&mesh, &bcs(scale), &opts).unwrap();
        assert!(
            stats.iterations > 0 && stats.residual <= opts.gmres_tol,
            "{stats:?}"
        );
        for ((u, q), (expected_u, expected_q)) in compressed
            .u
            .iter()
            .zip(&compressed.q)
            .zip(dense.u.iter().zip(&dense.q))
        {
            assert!(
                (u / scale - expected_u).abs() < 2e-8,
                "potential at scale {scale:e}"
            );
            assert!(
                (q / scale - expected_q).abs() < 2e-8,
                "flux at scale {scale:e}"
            );
        }
    }
}

#[test]
fn nonzero_neumann_data_keep_the_zero_mean_gauge() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let bcs = (0..mesh.triangles.len())
        .map(|i| BoundaryCondition::Flux(mesh.normal(i).unwrap()[0]))
        .collect::<Vec<_>>();
    let opts = AcaOptions {
        leaf_size: 4,
        aca_tol: 1e-12,
        ..AcaOptions::default()
    };
    let (solution, stats) = solve_laplace_bem_3d_aca(&mesh, &bcs, &opts).unwrap();
    assert!(
        stats.iterations > 0 && stats.residual <= opts.gmres_tol,
        "{stats:?}"
    );
    let weighted_sum: f64 = solution
        .u
        .iter()
        .enumerate()
        .map(|(i, u)| mesh.area(i).unwrap() * u)
        .sum();
    assert!(weighted_sum.abs() < 2e-10, "weighted mean {weighted_sum}");
    let point = [0.2, 0.1, -0.3];
    let value = evaluate_potential_3d(point, &mesh, &solution).unwrap();
    assert!(
        (value - point[0]).abs() < 0.035,
        "interior harmonic potential {value}"
    );
}

#[test]
fn mixed_p1_retains_the_rigid_diagonal_and_exterior_flux_sign() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let opts = AcaOptions {
        leaf_size: 4,
        eta: 0.0,
        ..AcaOptions::default()
    };
    for exterior in [false, true] {
        let bcs = mesh
            .vertices
            .iter()
            .map(|p| {
                let u = if exterior {
                    1.0 + 0.125 * p[0]
                } else {
                    p[0]
                };
                let q = if exterior {
                    -1.0 - 0.25 * p[0]
                } else {
                    p[0]
                };
                if p[2] < 0.0 {
                    BoundaryCondition::Potential(u)
                } else {
                    BoundaryCondition::Flux(q)
                }
            })
            .collect::<Vec<_>>();
        let dense = solve_laplace_bem_3d_p1(&mesh, &bcs, exterior).unwrap();
        let (compressed, stats) =
            solve_laplace_bem_3d_p1_aca(&mesh, &bcs, exterior, &opts).unwrap();
        assert!(
            stats.iterations > 0 && stats.residual <= opts.gmres_tol,
            "{stats:?}"
        );
        for ((u, q), (expected_u, expected_q)) in compressed
            .u
            .iter()
            .zip(&compressed.q)
            .zip(dense.u.iter().zip(&dense.q))
        {
            assert!(
                (u - expected_u).abs() < 2e-8,
                "potential, exterior = {exterior}"
            );
            assert!((q - expected_q).abs() < 2e-8, "flux, exterior = {exterior}");
        }
    }
}
