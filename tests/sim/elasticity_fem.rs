//! Linear elasticity on triangle / tetrahedron meshes.

use rssn::sim::boundary::{BcValue, ElasticBoundary, MarkerBcs};
use rssn::sim::elasticity_fem::{ElasticSolution, ElasticityProblem, Material, Model, constitutive_matrix, solve};
use rssn::sim::fem_simplex::{FemSpace, SimplexMesh, marker};
use rssn::sim::structured::LinearSolver;
use std::f64::consts::PI;

fn no_force(_x: [f64; 3]) -> [f64; 3] {
    [0.0; 3]
}

fn dir_exact(f: impl Fn([f64; 3]) -> [f64; 3] + Clone + Send + Sync + 'static) -> ElasticBoundary {
    ElasticBoundary::Displacement([0, 1, 2].map(|a| {
        let f = f.clone();
        Some(BcValue::space(move |x| f(x)[a]))
    }))
}

#[test]
fn patch_test_reproduces_linear_displacement_exactly_p1_p2_with_mixed_bcs() {
    let mat = Material { young: 100.0, poisson: 0.25, thickness: 1.0 };
    let u_ex = |x: [f64; 3]| [0.001 * x[0] + 0.0005 * x[1], -0.0003 * x[0] + 0.002 * x[1], 0.0];
    let c = constitutive_matrix(Model::PlaneStress, &mat);
    let eps = [0.001, 0.002, 0.0002];
    let s: Vec<f64> = (0..3).map(|r| (0..3).map(|k| c[r][k] * eps[k]).sum()).collect();
    for degree in [1, 2] {
        let bcs = MarkerBcs::new()
            .set(marker::XMIN, dir_exact(u_ex))
            .set(marker::YMIN, dir_exact(u_ex))
            .set(marker::XMAX, ElasticBoundary::traction([s[0], s[2], 0.0]))
            .set(marker::YMAX, ElasticBoundary::traction([s[2], s[1], 0.0]));
        let space = FemSpace::new(SimplexMesh::rectangle(5, 3, [0.0; 2], [2.0, 1.0]), degree).unwrap();
        let prob = ElasticityProblem { model: Model::PlaneStress, material: mat, force: &no_force, bcs: &bcs };
        let sol = solve(&space, &prob, LinearSolver::Auto).unwrap();
        for d in 0..space.size() {
            let ex = u_ex(space.dofs[d]);
            let v = sol.displacement(d);
            assert!((v[0] - ex[0]).abs() < 1e-12 && (v[1] - ex[1]).abs() < 1e-12, "P{degree} dof {d}");
        }
        let st = sol.stress_at_cell(3);
        for r in 0..3 {
            assert!((st[r] - s[r]).abs() < 1e-9);
        }
    }
}

#[test]
fn uniaxial_tension_with_rollers_in_3d_is_exact_for_tetrahedra() {
    let mat = Material { young: 200.0, poisson: 0.3, thickness: 1.0 };
    let sig = 1.0;
    let u_ex = move |x: [f64; 3]| [sig * x[0] / 200.0, -0.3 * sig * x[1] / 200.0, -0.3 * sig * x[2] / 200.0];
    for degree in [1, 2] {
        let bcs = MarkerBcs::new()
            .set(marker::XMIN, ElasticBoundary::roller(0))
            .set(marker::YMIN, ElasticBoundary::roller(1))
            .set(marker::ZMIN, ElasticBoundary::roller(2))
            .set(marker::XMAX, ElasticBoundary::traction([sig, 0.0, 0.0]));
        let space = FemSpace::new(SimplexMesh::boxed([2, 1, 1], [0.0; 3], [2.0, 1.0, 1.0]), degree).unwrap();
        let prob = ElasticityProblem { model: Model::Solid3d, material: mat, force: &no_force, bcs: &bcs };
        let sol = solve(&space, &prob, LinearSolver::Auto).unwrap();
        for d in 0..space.size() {
            let ex = u_ex(space.dofs[d]);
            let v = sol.displacement(d);
            for a in 0..3 {
                assert!((v[a] - ex[a]).abs() < 1e-10, "P{degree} dof {d} comp {a}: {} vs {}", v[a], ex[a]);
            }
        }
        let st = sol.stress_at_cell(0);
        assert!((st[0] - sig).abs() < 1e-8 && st[1].abs() < 1e-8);
        assert!((ElasticSolution::von_mises(&st) - sig).abs() < 1e-8);
    }
}

#[test]
fn cantilever_beam_2d_tip_deflection_matches_timoshenko_beam_theory() {
    let (l, h, p) = (10.0, 1.0, 1.0);
    let mat = Material { young: 1000.0, poisson: 0.3, thickness: 1.0 };
    let shear = move |x: [f64; 3]| {
        let s = 2.0 * (x[1] - 0.5 * h) / h;
        [0.0, -1.5 * p / h * (1.0 - s * s), 0.0]
    };
    let bcs = MarkerBcs::new().set(marker::XMIN, ElasticBoundary::clamped()).set(
        marker::XMAX,
        ElasticBoundary::Traction([BcValue::constant(0.0), BcValue::space(move |x| shear(x)[1]), BcValue::constant(0.0)]),
    );
    let space = FemSpace::new(SimplexMesh::rectangle(40, 4, [0.0; 2], [l, h]), 2).unwrap();
    let prob = ElasticityProblem { model: Model::PlaneStress, material: mat, force: &no_force, bcs: &bcs };
    let sol = solve(&space, &prob, LinearSolver::Auto).unwrap();
    let tip = sol.eval([l, 0.5 * h, 0.0]).unwrap()[1];
    let i = h.powi(3) / 12.0;
    let g = mat.young / (2.0 * (1.0 + mat.poisson));
    let exact = p * l.powi(3) / (3.0 * mat.young * i) + p * l / (5.0 / 6.0 * g * h);
    assert!(tip < 0.0 && ((-tip) / exact - 1.0).abs() < 0.04, "tip {tip} vs {exact}");
    // bending stress at the clamped end: sigma_xx = M y / I, M = P L
    let c = (0..sol.space.mesh.cells.len())
        .find(|&c| {
            let v = sol.space.mesh.cells[c].iter().map(|&q| sol.space.mesh.nodes[q]).collect::<Vec<_>>();
            v.iter().all(|n| n[0] < 0.3) && v.iter().all(|n| n[1] > 0.7)
        })
        .unwrap();
    let st = sol.stress_at_cell(c);
    assert!(st[0] > 0.0 && st[0] > 0.5 * p * l * 0.5 / i * 0.5, "top-fibre stress {}", st[0]);
}

#[test]
fn cantilever_beam_3d_tetrahedra() {
    let (l, p) = (8.0, 1.0);
    let mat = Material { young: 1000.0, poisson: 0.3, thickness: 1.0 };
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, ElasticBoundary::clamped())
        .set(marker::XMAX, ElasticBoundary::traction([0.0, 0.0, -p]));
    let space = FemSpace::new(SimplexMesh::boxed([16, 2, 2], [0.0; 3], [l, 1.0, 1.0]), 2).unwrap();
    let prob = ElasticityProblem { model: Model::Solid3d, material: mat, force: &no_force, bcs: &bcs };
    let sol = solve(&space, &prob, LinearSolver::Auto).unwrap();
    let tip = sol.eval([l, 0.5, 0.5]).unwrap()[2];
    let g = mat.young / (2.0 * (1.0 + mat.poisson));
    let exact = p * l.powi(3) / (3.0 * mat.young / 12.0) + p * l / (0.85 * g);
    assert!(tip < 0.0 && ((-tip) / exact - 1.0).abs() < 0.06, "tip {tip} vs {exact}");
}

fn manufactured_error(degree: usize, n: usize) -> f64 {
    let mat = Material { young: 1.0, poisson: 0.3, thickness: 1.0 };
    let (lam, mu) = (0.3 / (1.3 * 0.4), 1.0 / 2.6);
    let u1 = |x: [f64; 3]| (PI * x[0]).sin() * (PI * x[1]).sin();
    let a = |x: f64| x * (1.0 - x);
    let u2 = move |x: [f64; 3]| a(x[0]) * a(x[1]);
    let f = move |x: [f64; 3]| {
        let (sx, sy, cx, cy) = ((PI * x[0]).sin(), (PI * x[1]).sin(), (PI * x[0]).cos(), (PI * x[1]).cos());
        let (da, db) = (1.0 - 2.0 * x[0], 1.0 - 2.0 * x[1]);
        let u1v = sx * sy;
        let f1 = -((lam + mu) * (-PI * PI * u1v + da * db) + mu * (-2.0 * PI * PI * u1v));
        let f2 = -((lam + mu) * (PI * PI * cx * cy - 2.0 * a(x[0])) + mu * (-2.0 * a(x[1])));
        [f1, f2, 0.0]
    };
    let bcs = MarkerBcs::with_default(ElasticBoundary::clamped());
    let space = FemSpace::new(SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]), degree).unwrap();
    let prob = ElasticityProblem { model: Model::PlaneStrain, material: mat, force: &f, bcs: &bcs };
    let sol = solve(&space, &prob, LinearSolver::Auto).unwrap();
    let c1: Vec<f64> = (0..space.size()).map(|d| sol.u[2 * d]).collect();
    let c2: Vec<f64> = (0..space.size()).map(|d| sol.u[2 * d + 1]).collect();
    space.l2_error(&c1, u1).hypot(space.l2_error(&c2, u2))
}

#[test]
fn manufactured_plane_strain_solution_converges_with_optimal_order() {
    let (a, b) = (manufactured_error(1, 8), manufactured_error(1, 16));
    assert!((a / b).log2() > 1.8, "P1 {a} {b}");
    let (a, b) = (manufactured_error(2, 4), manufactured_error(2, 8));
    assert!((a / b).log2() > 2.7, "P2 {a} {b}");
}

#[test]
fn model_dimension_mismatch_is_rejected() {
    let space = FemSpace::new(SimplexMesh::rectangle(2, 2, [0.0; 2], [1.0; 2]), 1).unwrap();
    let bcs = MarkerBcs::new();
    let mat = Material { young: 1.0, poisson: 0.3, thickness: 1.0 };
    let prob = ElasticityProblem { model: Model::Solid3d, material: mat, force: &no_force, bcs: &bcs };
    assert!(solve(&space, &prob, LinearSolver::Auto).is_err());
}
