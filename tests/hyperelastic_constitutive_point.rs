//! Prepared response of the 3D material and its plane-strain use in a simplex.

use rssn::sim::boundary::ElasticBoundary;
use rssn::sim::boundary::MarkerBcs;
use rssn::sim::fem_simplex::SimplexMesh;
use rssn::sim::fem_simplex::marker;
use rssn::sim::hyperelastic::HyperMaterial;
use rssn::sim::hyperelastic::HyperModel;
use rssn::sim::hyperelastic::HyperProblem;
use rssn::sim::hyperelastic::HyperelasticPoint;
use rssn::sim::hyperelastic::first_piola;
use rssn::sim::hyperelastic::solve_static;
use rssn::sim::hyperelastic::strain_energy;
use rssn::sim::hyperelastic::tangent_moduli;

use rssn::sim::structured::LinearSolver;

type M3 = [[f64; 3]; 3];
const I: M3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const F: M3 = [[1.17, 0.13, -0.08], [0.15, 0.93, 0.04], [0.2, -0.1, 1.08]];
const H: M3 = [[0.3, -0.2, 0.1], [0.15, -0.4, 0.45], [-0.1, 0.05, 0.25]];
const ROTATION: M3 = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
const MODELS: [HyperModel; 2] = [HyperModel::StVenantKirchhoff, HyperModel::NeoHookean];

fn material(model: HyperModel) -> HyperMaterial {
    HyperMaterial {
        model,
        young: 1.0,
        poisson: 0.25,
        density: 2.0,
    }
}

fn next(
    f: M3,
    direction: M3,
    scale: f64,
) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| f[i][j] + scale * direction[i][j]))
}

fn mul(
    a: M3,
    b: M3,
) -> M3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j])
    })
}

fn near(
    x: f64,
    y: f64,
    tol: f64,
) {
    assert!(
        (x - y).abs() <= tol * y.abs().max(1.0),
        "{x:.17e} != {y:.17e}"
    );
}

fn near_matrix(
    a: M3,
    b: M3,
    tol: f64,
) {
    for (x, y) in a.into_iter().flatten().zip(b.into_iter().flatten()) {
        near(x, y, tol);
    }
}

#[test]
fn infinitesimal_limit_is_isotropic_linear_elasticity() {
    let epsilon = 2e-5;
    let tr = H[0][0] + H[1][1] + H[2][2];
    let symmetric: M3 = std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (H[i][j] + H[j][i])));
    for model in MODELS {
        let mat = material(model);
        let (lambda, mu) = mat.lame();
        let linear: M3 = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (if i == j { lambda * tr } else { 0.0 }) + 2.0 * mu * symmetric[i][j]
            })
        });
        let linear_energy =
            0.5 * lambda * tr * tr + mu * symmetric.iter().flatten().map(|x| x * x).sum::<f64>();
        let rest = HyperelasticPoint::new(&mat, &I).unwrap();
        near_matrix(rest.tangent_action(&H).unwrap(), linear, 2e-15);
        near(rest.energy_density().unwrap(), 0.0, 2e-15);
        let p = HyperelasticPoint::new(&mat, &next(I, H, epsilon)).unwrap();
        near_matrix(
            p.first_piola().unwrap().map(|row| row.map(|x| x / epsilon)),
            linear,
            3e-5,
        );
        near(
            p.energy_density().unwrap() / epsilon.powi(2),
            linear_energy,
            2e-5,
        );
    }
}

#[test]
fn directional_tangent_agrees_with_stress_difference_and_energy_derivative() {
    let step = 1e-5;
    for model in MODELS {
        let mat = material(model);
        let point = HyperelasticPoint::new(&mat, &F).unwrap();
        let tangent = point.tangent_action(&H).unwrap();
        let p1 = HyperelasticPoint::new(&mat, &next(F, H, step)).unwrap();
        let p2 = HyperelasticPoint::new(&mat, &next(F, H, -step)).unwrap();
        let (s1, s2) = (p1.first_piola().unwrap(), p2.first_piola().unwrap());
        let difference =
            std::array::from_fn(|i| std::array::from_fn(|j| (s1[i][j] - s2[i][j]) / (2.0 * step)));
        near_matrix(tangent, difference, 3e-10);

        let tensor = point.tangent_moduli().unwrap();
        let contraction: M3 = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                tensor[i][j]
                    .iter()
                    .flatten()
                    .zip(H.iter().flatten())
                    .map(|(x, h)| x * h)
                    .sum()
            })
        });
        near_matrix(tangent, contraction, 3e-15);

        let stress = point.first_piola().unwrap();
        let power = stress
            .iter()
            .flatten()
            .zip(H.iter().flatten())
            .map(|(p, h)| p * h)
            .sum::<f64>();
        let w1 = p1.energy_density().unwrap();
        let w2 = p2.energy_density().unwrap();
        near(power, (w1 - w2) / (2.0 * step), 3e-10);
    }
}

#[test]
fn proper_rotation_transforms_piola_and_directional_tangent() {
    for model in MODELS {
        let mat = material(model);
        let a = HyperelasticPoint::new(&mat, &F).unwrap();
        let b = HyperelasticPoint::new(&mat, &mul(ROTATION, F)).unwrap();
        near_matrix(
            b.first_piola().unwrap(),
            mul(ROTATION, a.first_piola().unwrap()),
            3e-15,
        );
        near_matrix(
            b.tangent_action(&mul(ROTATION, H)).unwrap(),
            mul(ROTATION, a.tangent_action(&H).unwrap()),
            3e-15,
        );
        near(
            b.energy_density().unwrap(),
            a.energy_density().unwrap(),
            3e-15,
        );
        let rigid = HyperelasticPoint::new(&mat, &ROTATION).unwrap();
        near_matrix(rigid.first_piola().unwrap(), [[0.0; 3]; 3], 3e-15);
        near(rigid.energy_density().unwrap(), 0.0, 3e-15);
    }
}

#[test]
fn neo_hooke_handles_large_positive_jacobians_and_unequal_rotated_stretches() {
    let mat = material(HyperModel::NeoHookean);
    let (lambda, mu) = mat.lame();
    let f = [[1e-150, 0.0, 0.0], [0.0, 1e200, 0.0], [0.0, 0.0, 1e150]];
    let ln_j = 200.0 * 10.0_f64.ln();
    let expected: M3 = std::array::from_fn(|i| {
        let v: f64 = f[i][i];
        let mut p = [0.0; 3];
        p[i] = mu * (v - 1.0 / v) + lambda * ln_j / v;
        p
    });

    for rotation in [I, ROTATION] {
        let gradient = mul(rotation, f);
        let point = HyperelasticPoint::new(&mat, &gradient).unwrap();
        near_matrix(point.first_piola().unwrap(), mul(rotation, expected), 2e-14);
        assert!(
            point
                .tangent_moduli()
                .unwrap()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .all(|x| x.is_finite())
        );

        // For a radial increment H=F/4, dlnJ=3/4 and Q H^T Q=Q/4.
        let h = gradient.map(|row| row.map(|x| x / 4.0));
        let radial: M3 = std::array::from_fn(|i| {
            let mut d = [0.0; 3];
            d[i] = 0.25 * (mu * f[i][i] + (mu - lambda * ln_j + 3.0 * lambda) / f[i][i]);
            d
        });
        near_matrix(
            point.tangent_action(&h).unwrap(),
            mul(rotation, radial),
            3e-14,
        );
    }

    // The true J is 1e450 while the strain energy density is still finite.
    let isotropic = I.map(|row| row.map(|x| 1e150 * x));
    let point = HyperelasticPoint::new(&mat, &isotropic).unwrap();
    near(point.energy_density().unwrap(), 1.5 * mu * 1e300, 3e-15);

    let soft = HyperMaterial { young: 1e-200, ..mat };
    let scaled = HyperelasticPoint::new(&soft, &f).unwrap();
    near(scaled.energy_density().unwrap(), 2e199, 5e-15);
    assert!(
        HyperelasticPoint::new(&mat, &f)
            .unwrap()
            .energy_density()
            .is_err()
    );
}

#[test]
fn invalid_jacobians_and_directions_are_errors() {
    for model in MODELS {
        let mat = material(model);
        for value in [0.0, -1.0, 1e-12, f64::INFINITY, f64::NAN] {
            let mut f = I;
            f[0][0] = value;
            assert!(first_piola(&mat, &f).is_err());
            assert!(tangent_moduli(&mat, &f).is_err());
        }
        let mut f = I;
        f[0][0] = 1e-12f64.next_up();
        assert!(first_piola(&mat, &f).is_ok());
        f[0][1] = f64::NAN;
        assert!(
            HyperelasticPoint::new(&mat, &I)
                .unwrap()
                .tangent_action(&f)
                .is_err()
        );
        let invalid_material = HyperMaterial { poisson: 0.5, ..mat };
        assert!(HyperelasticPoint::new(&invalid_material, &I).is_err());
    }
}

#[test]
fn affine_plane_strain_energy_matches_reference_volume_times_point_energy() {
    let mesh = SimplexMesh::rectangle(3, 2, [0.0; 2], [2.0, 3.0]);
    let bcs = MarkerBcs::new();
    let force = |_| [0.0; 3];
    let f = [[1.125, 0.0625, 0.0], [0.125, 0.875, 0.0], [0.0, 0.0, 1.0]];
    let displacement: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|x| [0.125 * x[0] + 0.0625 * x[1], 0.125 * x[0] - 0.125 * x[1]])
        .collect();
    for model in MODELS {
        let mat = material(model);
        let problem = HyperProblem {
            mesh: &mesh,
            material: mat,
            bcs: &bcs,
            force: &force,
            contacts: Vec::new(),
        };
        let local = HyperelasticPoint::new(&mat, &f)
            .unwrap()
            .energy_density()
            .unwrap();
        near(
            strain_energy(&problem, &displacement).unwrap(),
            6.0 * local,
            3e-15,
        );
        assert!(strain_energy(&problem, &displacement[..1]).is_err());
    }
}

#[test]
fn nonlinear_assembly_reports_inversion_and_iteration_exhaustion() {
    let mesh = SimplexMesh::rectangle(3, 2, [0.0; 2], [1.0; 2]);
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, ElasticBoundary::clamped())
        .set(marker::XMAX, ElasticBoundary::traction([0.1, 0.01, 0.0]));
    let force = |_| [0.0; 3];
    for model in MODELS {
        let problem = HyperProblem {
            mesh: &mesh,
            material: material(model),
            bcs: &bcs,
            force: &force,
            contacts: Vec::new(),
        };
        let err = solve_static(&problem, 1, 1e-12, 1, LinearSolver::Direct).unwrap_err();
        assert!(err.contains("Newton did not converge"), "{err}");
        let inverted: Vec<_> = mesh.nodes.iter().flat_map(|x| [-2.0 * x[0], 0.0]).collect();
        let err = strain_energy(&problem, &inverted).unwrap_err();
        assert!(err.contains("det F"), "{err}");
    }
}
