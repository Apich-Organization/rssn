//! Repeated boundary data for dense three-dimensional Laplace BEM.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use rssn::sim::physics_bem::BoundaryCondition;
use rssn::sim::physics_bem::BoundaryKind;
use rssn::sim::physics_bem::PreparedLaplaceBem3d;
use rssn::sim::physics_bem::SurfaceMesh3D;
use rssn::sim::physics_bem::evaluate_potential_3d;
use rssn::sim::physics_bem::solve_laplace_bem_3d_mesh;

#[test]
fn reuses_mixed_boundary_kinds_for_two_affine_fields() {
    let mesh = SurfaceMesh3D::cube(2.0, [0.0; 3], 4);
    let n = mesh.triangles.len();
    let kinds: Vec<_> = (0..n)
        .map(|t| {
            if mesh.normal(t).unwrap()[0].abs() > 0.5 {
                BoundaryKind::Potential
            } else {
                BoundaryKind::Flux
            }
        })
        .collect();
    let prepared = PreparedLaplaceBem3d::new(&mesh, &kinds).unwrap();
    assert_eq!(prepared.panel_count(), n);
    assert!(prepared.retained_bytes() >= 16 * n * n);
    for (offset, slope) in [(1.25, 0.75), (-2.5, -0.5)] {
        let values: Vec<_> = kinds
            .iter()
            .enumerate()
            .map(|(t, kind)| {
                match kind {
                    | BoundaryKind::Potential => offset + slope * mesh.centroid(t).unwrap()[0],
                    | BoundaryKind::Flux => 0.0,
                }
            })
            .collect();
        let bcs: Vec<_> = kinds
            .iter()
            .zip(&values)
            .map(|(kind, &value)| {
                match kind {
                    | BoundaryKind::Potential => BoundaryCondition::Potential(value),
                    | BoundaryKind::Flux => BoundaryCondition::Flux(value),
                }
            })
            .collect();
        let warm = prepared.solve(&values).unwrap();
        let cold = solve_laplace_bem_3d_mesh(&mesh, &bcs).unwrap();
        for (a, b) in warm
            .u
            .iter()
            .chain(&warm.q)
            .zip(cold.u.iter().chain(&cold.q))
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        for (t, kind) in kinds.iter().enumerate() {
            let normal = mesh.normal(t).unwrap()[0];
            match kind {
                | BoundaryKind::Potential => {
                    assert_eq!(warm.u[t].to_bits(), values[t].to_bits());
                    assert!(warm.q[t] * (slope * normal) > 0.0);
                },
                | BoundaryKind::Flux => assert_eq!(warm.q[t].to_bits(), values[t].to_bits()),
            }
        }
        for x in [-0.4, 0.2] {
            let u = evaluate_potential_3d([x, 0.1, -0.05], &mesh, &warm).unwrap();
            assert!((u - (offset + slope * x)).abs() < 0.025, "{x}: {u}");
        }
    }
}

#[test]
fn pure_neumann_uses_the_zero_area_weighted_mean() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let n = mesh.triangles.len();
    let prepared = PreparedLaplaceBem3d::new(&mesh, &vec![BoundaryKind::Flux; n]).unwrap();
    let outward: Vec<_> = (0..n).map(|t| mesh.normal(t).unwrap()[0]).collect();
    let result = prepared.solve(&outward).unwrap();
    let mean = result
        .u
        .iter()
        .enumerate()
        .map(|(t, u)| mesh.area(t).unwrap() * u)
        .sum::<f64>()
        / mesh.total_area();
    assert!(mean.abs() < 2e-15);
    let err = result
        .u
        .iter()
        .enumerate()
        .map(|(t, u)| (u - mesh.centroid(t).unwrap()[0]).powi(2))
        .sum::<f64>();
    assert!((err / n as f64).sqrt() < 0.08);
    assert!(
        prepared
            .solve(&vec![1.0; n])
            .unwrap_err()
            .contains("net flux")
    );
}

#[test]
fn flux_compatibility_survives_product_overflow_and_underflow() {
    for (size, flux) in [(3.0, 1e308), (0.25, f64::from_bits(1))] {
        let mesh = SurfaceMesh3D::cube(size, [0.0; 3], 1);
        let n = mesh.triangles.len();
        let prepared = PreparedLaplaceBem3d::new(&mesh, &vec![BoundaryKind::Flux; n]).unwrap();
        let mut values = vec![0.0; n];
        values[0] = flux;
        values[1] = -0.5 * flux;
        assert!(prepared.solve(&values).unwrap_err().contains("net flux"));
        let bcs: Vec<_> = values
            .iter()
            .copied()
            .map(BoundaryCondition::Flux)
            .collect();
        assert!(
            solve_laplace_bem_3d_mesh(&mesh, &bcs)
                .unwrap_err()
                .contains("net flux")
        );

        values[1] = -flux;
        let solution = prepared.solve(&values).unwrap();
        assert!(solution.u.iter().all(|u| u.is_finite()));
        if flux > 1.0 {
            let unit_values: Vec<_> = values.iter().map(|q| q / flux).collect();
            let unit_solution = prepared.solve(&unit_values).unwrap();
            for (big, unit) in solution.u.iter().zip(&unit_solution.u) {
                assert!((big / flux - unit).abs() < 1e-14);
            }
        }
    }
}

#[test]
fn prepared_shape_finiteness_and_cancellation_are_checked() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 0);
    let n = mesh.triangles.len();
    let kinds = vec![BoundaryKind::Potential; n];
    assert!(PreparedLaplaceBem3d::new(&mesh, &kinds[..n - 1]).is_err());
    let prepared = PreparedLaplaceBem3d::new(&mesh, &kinds).unwrap();
    assert!(prepared.solve(&[]).is_err());
    assert!(prepared.solve(&vec![f64::NAN; n]).is_err());
    let values = vec![2.0; n];
    {
        let _cancel = rssn::cancellation::install(Arc::new(AtomicBool::new(true)));
        assert!(PreparedLaplaceBem3d::new(&mesh, &kinds).is_err());
        assert!(prepared.solve(&values).unwrap_err().contains("cancelled"));
    }
    let result = prepared.solve(&values).unwrap();
    assert!(result.q.iter().all(|q| q.abs() < 1e-12));
}
