//! Translated-domain patch tests shared by the simplicial elasticity assembly.

use std::cell::RefCell;

use rssn::sim::boundary::ElasticBoundary;
use rssn::sim::boundary::MarkerBcs;
use rssn::sim::elasticity_fem::ElasticityProblemVar;
use rssn::sim::elasticity_fem::Model;
use rssn::sim::elasticity_fem::solve_heterogeneous;
use rssn::sim::fem_simplex::FemSpace;
use rssn::sim::fem_simplex::SimplexMesh;
use rssn::sim::fem_simplex::marker;
use rssn::sim::structured::LinearSolver;

#[test]
fn translated_heterogeneous_patch_uses_only_quadrature_and_centroid_samples() {
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, ElasticBoundary::clamped())
        .set(marker::XMAX, ElasticBoundary::traction([2.5, 0.0, 0.0]));
    for degree in [1, 2] {
        for (mesh, model) in [
            (
                SimplexMesh::rectangle(3, 2, [2.0, -3.0], [3.0, -2.0]),
                Model::PlaneStress,
            ),
            (
                SimplexMesh::boxed([2, 1, 1], [2.0, -3.0, 4.0], [3.0, -2.0, 5.0]),
                Model::Solid3d,
            ),
        ] {
            let space = FemSpace::new(mesh, degree).unwrap();
            let young_samples = RefCell::new(Vec::new());
            let poisson_samples = RefCell::new(Vec::new());
            let body_samples = RefCell::new(Vec::new());
            let young = |x: [f64; 3]| {
                young_samples.borrow_mut().push(x);
                if (2.0..=3.0).contains(&x[0]) && (-3.0..=-2.0).contains(&x[1]) {
                    64.0 + 16.0 * (x[0] - 2.0)
                } else {
                    f64::NAN
                }
            };
            let poisson = |x: [f64; 3]| {
                poisson_samples.borrow_mut().push(x);
                0.0
            };
            let force = |x: [f64; 3]| {
                body_samples.borrow_mut().push(x);
                [-0.5, 0.0, 0.0]
            };
            let prob = ElasticityProblemVar {
                model,
                young: &young,
                poisson: &poisson,
                force: &force,
                thickness: 0.625,
                bcs: &bcs,
            };
            let sol = solve_heterogeneous(&space, &prob, LinearSolver::Direct).unwrap();
            for (d, x) in space.dofs.iter().enumerate() {
                let u = sol.displacement(d);
                assert!(
                    (u[0] - 0.03125 * (x[0] - 2.0)).abs() < 2e-12,
                    "P{degree} {model:?} {d}: {u:?}"
                );
                assert!(u[1..].iter().all(|v| v.abs() < 2e-12));
            }
            for (cell, nodes) in space.mesh.cells.iter().enumerate() {
                let x_mean =
                    nodes.iter().map(|&j| space.mesh.nodes[j][0]).sum::<f64>() / nodes.len() as f64;
                let exact_xx = 2.0 + 0.5 * (x_mean - 2.0);
                let stress = sol.stress_at_cell(cell);
                assert!(
                    (stress[0] - exact_xx).abs() < 8e-11,
                    "P{degree} {model:?} {cell}: {stress:?}"
                );
                assert!(stress[1..].iter().all(|v| v.abs() < 8e-11));
            }
            let young = young_samples.borrow();
            let body = body_samples.borrow();
            assert_eq!(young.len(), body.len() + space.mesh.cells.len());
            assert_eq!(*young, *poisson_samples.borrow());
            let body_per_cell = body.len() / space.mesh.cells.len();
            for (material, body) in young
                .chunks_exact(body_per_cell + 1)
                .zip(body.chunks_exact(body_per_cell))
            {
                assert_eq!(&material[..body_per_cell], body);
            }
        }
    }
}
