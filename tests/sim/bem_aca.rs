//! Iterative ACA-compressed 3D BEM tests.

use rssn::sim::bem_aca::{AcaOptions, solve_laplace_bem_3d_aca};
use rssn::sim::physics_bem::{BoundaryCondition, SurfaceMesh3D, evaluate_potential_3d, solve_laplace_bem_3d_mesh};

fn dirichlet_bcs(mesh: &SurfaceMesh3D, u: impl Fn([f64; 3]) -> f64) -> Vec<BoundaryCondition<f64>> {
    (0..mesh.triangles.len()).map(|t| BoundaryCondition::Potential(u(mesh.centroid(t).unwrap()))).collect()
}

#[test]
fn aca_gmres_matches_dense_lu_on_the_sphere_for_dirichlet_and_mixed_data() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);
    let opts = AcaOptions { leaf_size: 16, aca_tol: 1e-8, ..AcaOptions::default() };
    // Dirichlet u = x, exact flux q = n_x
    let bcs = dirichlet_bcs(&mesh, |c| c[0]);
    let dense = solve_laplace_bem_3d_mesh(&mesh, &bcs).unwrap();
    let (fast, stats) = solve_laplace_bem_3d_aca(&mesh, &bcs, &opts).unwrap();
    let diff = (0..bcs.len()).map(|t| (dense.q[t] - fast.q[t]).abs()).fold(0.0, f64::max);
    assert!(diff < 1e-5, "ACA vs dense {diff}");
    assert!(stats.lowrank_blocks > 0 && stats.iterations < 40, "{stats:?}");
    let exact = (0..bcs.len()).map(|t| (fast.q[t] - mesh.normal(t).unwrap()[0]).abs()).fold(0.0, f64::max);
    assert!(exact < 0.05, "flux error {exact}");

    // mixed: u = x y + z, Dirichlet on the lower hemisphere, Neumann (q = grad u . n) on the upper
    let u_ex = |c: [f64; 3]| c[0] * c[1] + c[2];
    let mixed: Vec<BoundaryCondition<f64>> = (0..mesh.triangles.len())
        .map(|t| {
            let c = mesh.centroid(t).unwrap();
            let n = mesh.normal(t).unwrap();
            if c[2] < 0.0 {
                BoundaryCondition::Potential(u_ex(c))
            } else {
                BoundaryCondition::Flux(c[1] * n[0] + c[0] * n[1] + n[2])
            }
        })
        .collect();
    let dense = solve_laplace_bem_3d_mesh(&mesh, &mixed).unwrap();
    let (fast, _) = solve_laplace_bem_3d_aca(&mesh, &mixed, &opts).unwrap();
    let diff = (0..mixed.len()).map(|t| (dense.u[t] - fast.u[t]).abs().max((dense.q[t] - fast.q[t]).abs())).fold(0.0, f64::max);
    assert!(diff < 1e-5, "mixed ACA vs dense {diff}");
    let worst = (0..mixed.len()).map(|t| (fast.u[t] - u_ex(mesh.centroid(t).unwrap())).abs()).fold(0.0, f64::max);
    assert!(worst < 0.06, "potential error {worst}");
}

#[test]
fn aca_compresses_the_far_field_and_converges_on_a_finer_sphere() {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 3);
    let bcs = dirichlet_bcs(&mesh, |c| c[0] * c[1]);
    let (sol, stats) = solve_laplace_bem_3d_aca(&mesh, &bcs, &AcaOptions::default()).unwrap();
    assert!(stats.compression() < 0.85, "compression {}", stats.compression());
    assert!(stats.lowrank_blocks > stats.dense_blocks / 4);
    // u = x y has q = 2 x y on the unit sphere
    let worst = (0..bcs.len())
        .map(|t| {
            let c = mesh.centroid(t).unwrap();
            let n = mesh.normal(t).unwrap();
            (sol.q[t] - (c[1] * n[0] + c[0] * n[1])).abs()
        })
        .fold(0.0, f64::max);
    assert!(worst < 0.04, "flux error {worst}");
    // interior evaluation agrees with the exact harmonic function
    let p = evaluate_potential_3d([0.2, 0.3, 0.1], &mesh, &sol).unwrap();
    assert!((p - 0.06).abs() < 5e-3, "interior value {p}");
}

#[test]
fn cube_with_point_source_data_and_error_cases() {
    let mesh = SurfaceMesh3D::cube(1.0, [0.0; 3], 6);
    let s = [2.0, 0.5, 0.3];
    let phi = |c: [f64; 3]| 1.0 / ((c[0] - s[0]).powi(2) + (c[1] - s[1]).powi(2) + (c[2] - s[2]).powi(2)).sqrt();
    let bcs = dirichlet_bcs(&mesh, phi);
    let (sol, _) = solve_laplace_bem_3d_aca(&mesh, &bcs, &AcaOptions { leaf_size: 8, ..AcaOptions::default() }).unwrap();
    let mut err: f64 = 0.0;
    for t in 0..mesh.triangles.len() {
        let c = mesh.centroid(t).unwrap();
        let n = mesh.normal(t).unwrap();
        let r2 = (c[0] - s[0]).powi(2) + (c[1] - s[1]).powi(2) + (c[2] - s[2]).powi(2);
        let g = [-(c[0] - s[0]) / r2.powf(1.5), -(c[1] - s[1]) / r2.powf(1.5), -(c[2] - s[2]) / r2.powf(1.5)];
        err = err.max((sol.q[t] - (g[0] * n[0] + g[1] * n[1] + g[2] * n[2])).abs());
    }
    assert!(err < 0.08, "cube flux error {err}");
    // pure Neumann is rejected, wrong data length is rejected
    let neu: Vec<BoundaryCondition<f64>> = (0..mesh.triangles.len()).map(|_| BoundaryCondition::Flux(0.0)).collect();
    assert!(solve_laplace_bem_3d_aca(&mesh, &neu, &AcaOptions::default()).is_err());
    assert!(solve_laplace_bem_3d_aca(&mesh, &neu[..3], &AcaOptions::default()).is_err());
}
