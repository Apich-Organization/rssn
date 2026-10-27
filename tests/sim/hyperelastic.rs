//! Large-deformation elasticity, elastodynamics and penalty contact.

use rssn::sim::boundary::{BcValue, ElasticBoundary, MarkerBcs};
use rssn::sim::elasticity_fem::{ElasticityProblem, Material, Model, solve as solve_linear};
use rssn::sim::fem_simplex::{FemSpace, SimplexMesh, marker};
use rssn::sim::hyperelastic::{
    DynamicState, HyperMaterial, HyperModel, HyperProblem, PenaltyContact, TimeScheme, first_piola, kinetic_energy, rest_state, solve_dynamic,
    solve_static, strain_energy, tangent_moduli,
};
use rssn::sim::structured::LinearSolver;

fn material(model: HyperModel) -> HyperMaterial {
    HyperMaterial { model, young: 1000.0, poisson: 0.3, density: 2.0 }
}

const F_SAMPLE: [[f64; 3]; 3] = [[1.2, 0.3, -0.1], [0.1, 0.9, 0.2], [-0.2, 0.15, 1.1]];

#[test]
fn piola_stress_tangent_matches_finite_differences_and_is_objective() {
    for model in [HyperModel::StVenantKirchhoff, HyperModel::NeoHookean] {
        let mat = material(model);
        let a = tangent_moduli(&mat, &F_SAMPLE).unwrap();
        let eps = 1e-6;
        for k in 0..3 {
            for l in 0..3 {
                let (mut fp, mut fm) = (F_SAMPLE, F_SAMPLE);
                fp[k][l] += eps;
                fm[k][l] -= eps;
                let (pp, pm) = (first_piola(&mat, &fp).unwrap(), first_piola(&mat, &fm).unwrap());
                for i in 0..3 {
                    for j in 0..3 {
                        let fd = (pp[i][j] - pm[i][j]) / (2.0 * eps);
                        assert!((fd - a[i][j][k][l]).abs() < 1e-5 * (1.0 + fd.abs()), "{model:?} A[{i}{j}{k}{l}] = {} vs {fd}", a[i][j][k][l]);
                    }
                }
            }
        }
        // a rigid rotation carries no stress (frame indifference), unlike small-strain theory
        let (s, c) = 1.0_f64.sin_cos();
        let r = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let p = first_piola(&mat, &r).unwrap();
        assert!(p.iter().flatten().all(|v| v.abs() < 1e-12), "{model:?}: stress of a rotation {p:?}");
        assert!(first_piola(&mat, &[[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]).is_err());
    }
}

fn uniaxial_case(model: HyperModel, stretch: f64, steps: usize) -> (f64, f64, Vec<f64>, usize) {
    let mesh = SimplexMesh::boxed([2, 2, 2], [0.0; 3], [1.0; 3]);
    let mat = material(model);
    let mut bcs = MarkerBcs::new();
    bcs = bcs.set(marker::XMIN, ElasticBoundary::roller(0)).set(marker::YMIN, ElasticBoundary::roller(1)).set(marker::ZMIN, ElasticBoundary::roller(2));
    bcs = bcs.set(marker::XMAX, ElasticBoundary::Displacement([Some(BcValue::constant(stretch - 1.0)), None, None]));
    let zero = |_: [f64; 3]| [0.0; 3];
    let prob = HyperProblem { mesh: &mesh, material: mat, bcs: &bcs, force: &zero, contacts: Vec::new() };
    let sol = solve_static(&prob, steps, 1e-12, 30, LinearSolver::Direct).unwrap();
    // lateral stretch at the far corner (node with coordinates (1,1,1))
    let node = mesh.nodes.iter().position(|x| (x[0] - 1.0).abs() < 1e-12 && (x[1] - 1.0).abs() < 1e-12 && (x[2] - 1.0).abs() < 1e-12).unwrap();
    let lam2 = 1.0 + sol.u[node * 3 + 1];
    let lam3 = 1.0 + sol.u[node * 3 + 2];
    assert!((lam2 - lam3).abs() < 1e-9);
    (lam2, mat.lame().0, sol.residuals, sol.iterations.iter().sum())
}

#[test]
fn uniaxial_stretch_reproduces_the_analytic_lateral_contraction() {
    let l1 = 1.4;
    // St. Venant-Kirchhoff: S_22 = 0  =>  E_22 = -lambda E_11 / (2 (lambda + mu))
    let mat = material(HyperModel::StVenantKirchhoff);
    let (lam, mu) = mat.lame();
    let e11 = 0.5 * (l1 * l1 - 1.0);
    let e22 = -lam * e11 / (2.0 * (lam + mu));
    let exact_svk = (1.0 + 2.0 * e22).sqrt();
    let (got, _, _, _) = uniaxial_case(HyperModel::StVenantKirchhoff, l1, 4);
    assert!((got - exact_svk).abs() < 1e-9, "SVK lateral stretch {got} vs {exact_svk}");
    // neo-Hookean: P_22 = mu (l2 - 1/l2) + lambda ln(J) / l2 = 0 with J = l1 l2^2 (bisection)
    let mat = material(HyperModel::NeoHookean);
    let (lam, mu) = mat.lame();
    let f = |l2: f64| mu * (l2 - 1.0 / l2) + lam * (l1 * l2 * l2).ln() / l2;
    let (mut lo, mut hi) = (0.3, 1.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid } else { lo = mid }
    }
    let exact_nh = 0.5 * (lo + hi);
    let (got, _, res, iters) = uniaxial_case(HyperModel::NeoHookean, l1, 4);
    assert!((got - exact_nh).abs() < 1e-9, "neo-Hookean lateral stretch {got} vs {exact_nh}");
    assert!(iters <= 4 * 6, "Newton iterations {iters}");
    assert!(*res.last().unwrap() < 1e-6, "final residual {res:?}");
}

fn cantilever(model: HyperModel, load: f64, steps: usize) -> (f64, Vec<f64>) {
    let mesh = SimplexMesh::rectangle(24, 6, [0.0, 0.0], [4.0, 1.0]);
    let mat = material(model);
    let bcs = MarkerBcs::new().set(marker::XMIN, ElasticBoundary::clamped()).set(marker::XMAX, ElasticBoundary::traction([0.0, -load, 0.0]));
    let zero = |_: [f64; 3]| [0.0; 3];
    let prob = HyperProblem { mesh: &mesh, material: mat, bcs: &bcs, force: &zero, contacts: Vec::new() };
    let sol = solve_static(&prob, steps, 1e-12, 40, LinearSolver::Direct).unwrap();
    // mean tip deflection of the right edge
    let tip: Vec<usize> = (0..mesh.nodes.len()).filter(|&i| (mesh.nodes[i][0] - 4.0).abs() < 1e-12).collect();
    let d = tip.iter().map(|&i| sol.u[2 * i + 1]).sum::<f64>() / tip.len() as f64;
    (d, sol.residuals)
}

#[test]
fn newton_converges_quadratically_and_reduces_to_linear_elasticity() {
    // small-strain limit against the linear solver (same P1 mesh, plane strain)
    let load = 1e-3;
    let (d_nl, _) = cantilever(HyperModel::StVenantKirchhoff, load, 1);
    let mesh = SimplexMesh::rectangle(24, 6, [0.0, 0.0], [4.0, 1.0]);
    let space = FemSpace::new(mesh.clone(), 1).unwrap();
    let bcs = MarkerBcs::new().set(marker::XMIN, ElasticBoundary::clamped()).set(marker::XMAX, ElasticBoundary::traction([0.0, -load, 0.0]));
    let zero = |_: [f64; 3]| [0.0; 3];
    let prob = ElasticityProblem { model: Model::PlaneStrain, material: Material { young: 1000.0, poisson: 0.3, thickness: 1.0 }, force: &zero, bcs: &bcs };
    let lin = solve_linear(&space, &prob, LinearSolver::Direct).unwrap();
    let tip: Vec<usize> = (0..mesh.nodes.len()).filter(|&i| (mesh.nodes[i][0] - 4.0).abs() < 1e-12).collect();
    let d_lin = tip.iter().map(|&i| lin.displacement(i)[1]).sum::<f64>() / tip.len() as f64;
    assert!((d_nl - d_lin).abs() < 2e-3 * d_lin.abs(), "small-strain limit {d_nl} vs linear {d_lin}");
    // large deflection: the cantilever stiffens geometrically relative to the linear prediction
    let big = 10.0;
    let (d_big, res) = cantilever(HyperModel::StVenantKirchhoff, big, 20);
    let d_big_lin = d_lin * big / load;
    assert!(d_big.abs() < 0.9 * d_big_lin.abs() && d_big.abs() > 0.3 * d_big_lin.abs(), "large deflection {d_big} vs linear {d_big_lin}");
    // quadratic convergence of the last load step: r_{k+1} ~ r_k^2
    assert!(res.len() >= 3 && *res.last().unwrap() < 1e-8 * res[0].max(1.0), "{res:?}");
    for w in res.windows(2).skip(1) {
        if w[0] > 1e-8 * res[0] {
            assert!(w[1] < 5.0 * w[0] * w[0] / res[0] + 1e-14, "not quadratic: {res:?}");
        }
    }
    // neo-Hookean gives a similar response
    let (d_nh, _) = cantilever(HyperModel::NeoHookean, big, 20);
    assert!((d_nh - d_big).abs() < 0.15 * d_big.abs(), "NH {d_nh} vs SVK {d_big}");
}

fn bar_problem<'a>(mesh: &'a SimplexMesh, bcs: &'a MarkerBcs<ElasticBoundary>, force: &'a dyn Fn([f64; 3]) -> [f64; 3]) -> HyperProblem<'a> {
    HyperProblem { mesh, material: material(HyperModel::StVenantKirchhoff), bcs, force, contacts: Vec::new() }
}

fn bar_bcs(traction: f64) -> MarkerBcs<ElasticBoundary> {
    MarkerBcs::new()
        .set(marker::XMIN, ElasticBoundary::clamped())
        .set(marker::YMIN, ElasticBoundary::roller(1))
        .set(marker::YMAX, ElasticBoundary::roller(1))
        .set(marker::ZMIN, ElasticBoundary::roller(2))
        .set(marker::ZMAX, ElasticBoundary::roller(2))
        .set(marker::XMAX, ElasticBoundary::traction([traction, 0.0, 0.0]))
}

fn run_bar(scheme: TimeScheme, dt_factor: f64, periods: f64) -> (f64, Vec<f64>, Vec<f64>) {
    // returns (period ratio, total energy / (F u_peak) history, [scale])
    // tension bar clamped at x = 0, lateral rollers: 1D wave speed c = sqrt((lambda + 2 mu) / rho)
    let length = 1.0;
    let mesh = SimplexMesh::boxed([20, 1, 1], [0.0; 3], [length, 0.2, 0.2]);
    let bcs = bar_bcs(0.01);
    let zero = |_: [f64; 3]| [0.0; 3];
    let prob = bar_problem(&mesh, &bcs, &zero);
    let mat = prob.material;
    let (lam, mu) = mat.lame();
    let c = ((lam + 2.0 * mu) / mat.density).sqrt();
    let period = 4.0 * length / c;
    let dt = period * dt_factor;
    let n = mesh.nodes.len() * 3;
    let init = rest_state(&prob, &vec![0.0; n], LinearSolver::Direct).unwrap();
    let tip = mesh.nodes.iter().position(|x| (x[0] - length).abs() < 1e-12 && x[1] < 1e-12 && x[2] < 1e-12).unwrap();
    let mut history = Vec::new();
    let mut energy = Vec::new();
    let mut obs = |s: &DynamicState| {
        history.push((s.t, s.u[tip * 3]));
        let work: f64 = {
            // external work of the constant end load: traction * area * u averaged over the end face nodes
            let end: Vec<usize> = (0..mesh.nodes.len()).filter(|&i| (mesh.nodes[i][0] - length).abs() < 1e-12).collect();
            end.iter().map(|&i| s.u[i * 3]).sum::<f64>() / end.len() as f64 * 0.01 * 0.04
        };
        energy.push(kinetic_energy(&prob, &s.v) + strain_energy(&prob, &s.u).unwrap() - work);
    };
    let steps = (periods / dt_factor) as usize;
    solve_dynamic(&prob, scheme, &init, dt, steps, 1e-12, LinearSolver::Direct, &mut obs).unwrap();
    // period from two successive upward crossings of the mean level (linear interpolation)
    let (hi, lo) = history.iter().fold((f64::MIN, f64::MAX), |(h, l), v| (h.max(v.1), l.min(v.1)));
    let level = 0.5 * (hi + lo);
    let mut crossings = Vec::new();
    for w in history.windows(2) {
        if w[0].1 < level && w[1].1 >= level {
            crossings.push(w[0].0 + (w[1].0 - w[0].0) * (level - w[0].1) / (w[1].1 - w[0].1));
        }
    }
    assert!(crossings.len() >= 2, "need two crossings, got {crossings:?}");
    // energy scale: external force times the peak displacement
    let scale = 0.01 * 0.04 * hi;
    ((crossings[1] - crossings[0]) / period, energy.iter().map(|e| e / scale).collect(), vec![scale])
}

#[test]
fn newmark_bar_vibration_has_the_analytic_period_and_conserves_energy() {
    let (ratio, energy, _) = run_bar(TimeScheme::Newmark { beta: 0.25, gamma: 0.5 }, 1.0 / 80.0, 2.6);
    assert!((ratio - 1.0).abs() < 0.015, "period / (4 L / c) = {ratio}");
    // the total energy (kinetic + strain - work of the load) starts at zero and stays there
    let drift = energy.iter().map(|e| e.abs()).fold(0.0_f64, f64::max);
    assert!(drift < 6e-3, "relative energy drift {drift}");
}

#[test]
fn generalized_alpha_dissipates_high_frequencies_but_keeps_the_fundamental_period() {
    let (ratio, energy, _) = run_bar(TimeScheme::GeneralizedAlpha { rho_inf: 0.4 }, 1.0 / 40.0, 2.6);
    assert!((ratio - 1.0).abs() < 0.04, "period ratio {ratio}");
    let (_, conserved, _) = run_bar(TimeScheme::Newmark { beta: 0.25, gamma: 0.5 }, 1.0 / 40.0, 2.6);
    let worst = |e: &[f64]| e.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    assert!(worst(&energy) > 5.0 * worst(&conserved) + 1e-4, "dissipation {} vs {}", worst(&energy), worst(&conserved));
    // energy is only ever removed (unconditional stability)
    assert!(energy.windows(2).all(|w| w[1] <= w[0] + 1e-3), "energy increased");
}

#[test]
fn penalty_contact_with_a_rigid_plane_matches_the_analytic_penetration() {
    let delta = 2e-4;
    let mesh = SimplexMesh::rectangle(4, 4, [0.0, 0.0], [1.0, 1.0]);
    let mat = material(HyperModel::NeoHookean);
    let bcs = MarkerBcs::new().set(marker::XMIN, ElasticBoundary::roller(0)).set(marker::XMAX, ElasticBoundary::roller(0)).set(marker::YMIN, ElasticBoundary::roller(1));
    let zero = |_: [f64; 3]| [0.0; 3];
    let eps = 1e4;
    let contact = PenaltyContact {
        marker: marker::YMAX,
        gap: Box::new(move |x| (1.0 - delta) - x[1]),
        normal: Box::new(|_| [0.0, -1.0, 0.0]),
        penalty: eps,
    };
    let prob = HyperProblem { mesh: &mesh, material: mat, bcs: &bcs, force: &zero, contacts: vec![contact] };
    let sol = solve_static(&prob, 1, 1e-12, 40, LinearSolver::Direct).unwrap();
    let (e, nu) = (mat.young, mat.poisson);
    let m = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let p = eps * delta / (1.0 + eps / m);
    let u_top = -p / m;
    for (i, x) in mesh.nodes.iter().enumerate() {
        if (x[1] - 1.0).abs() < 1e-12 {
            let got = sol.u[2 * i + 1];
            assert!((got - u_top).abs() < 3e-3 * u_top.abs(), "top displacement {got} vs {u_top}");
        }
    }
    // without the obstacle the block does not move
    let free = HyperProblem { mesh: &mesh, material: mat, bcs: &bcs, force: &zero, contacts: Vec::new() };
    let s0 = solve_static(&free, 1, 1e-12, 10, LinearSolver::Direct).unwrap();
    assert!(s0.u.iter().all(|v| v.abs() < 1e-12));
}
