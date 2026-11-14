//! Plane-stress linear elasticity with Q4 elements (ported from `physics_sim_linear_elasticity_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::models::linear_elasticity::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

const SQ: [(f64, f64); 4] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];

fn k_of(
    e: f64,
    nu: f64,
) -> ndarray::Array2<f64> {
    element_stiffness_matrix(SQ[0], SQ[1], SQ[2], SQ[3], e, nu)
}

#[test]
fn stiffness_matrix_is_8x8_and_symmetric() {
    let k = k_of(1e7, 0.3);
    assert_eq!(k.shape(), &[8, 8]);
    for i in 0..8 {
        for j in 0..8 {
            assert!((k[[i, j]] - k[[j, i]]).abs() < 1e-9);
        }
    }
}

#[test]
fn rigid_translations_produce_no_force() {
    let k = k_of(1e7, 0.3);
    let tx = ndarray::arr1(&[1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]);
    let ty = ndarray::arr1(&[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    for t in [tx, ty] {
        let f = k.dot(&t);
        assert!(f.iter().all(|v| v.abs() < 1e-6), "{f:?}");
    }
}

#[test]
fn stiffness_matrix_is_positive_semidefinite() {
    let k = k_of(1e7, 0.3);
    // Deterministic probes of the quadratic form u^T K u.
    for s in 0..50 {
        let u: Vec<f64> = (0..8).map(|i| ((s * 8 + i) as f64 * 1.7).sin()).collect();
        let ku = k.dot(&ndarray::arr1(&u));
        let q: f64 = u.iter().zip(ku.iter()).map(|(a, b)| a * b).sum();
        assert!(q > -1e-6, "u^T K u = {q}");
    }
}

#[test]
fn stiffness_diagonal_has_the_analytic_q4_value() {
    // Unit square, N1 = (1-x)(1-y): K00 = c (int (1-y)^2 + (1-nu)/2 int (1-x)^2) = c/3 (1 + (1-nu)/2),
    // with c = E / (1 - nu^2). For nu = 0, E = 1 this is exactly 1/2.
    assert!((k_of(1.0, 0.0)[[0, 0]] - 0.5).abs() < 1e-12);
    let (e, nu) = (2.0, 0.25);
    let c = e / (1.0 - nu * nu);
    let expected = c / 3.0 * (1.0 + (1.0 - nu) / 2.0);
    assert!((k_of(e, nu)[[0, 0]] - expected).abs() < 1e-12);
}

#[test]
fn stiffness_uses_the_node_coordinates() {
    // Unit thickness plane stress: scaling the element in-plane leaves K unchanged...
    let (e, nu) = (1e5, 0.3);
    let k1 = k_of(e, nu);
    let k3 = element_stiffness_matrix((0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), e, nu);
    for i in 0..8 {
        for j in 0..8 {
            assert!((k1[[i, j]] - k3[[i, j]]).abs() < 1e-8 * (1.0 + k1[[i, j]].abs()));
        }
    }
    // ...but a stretched element is not the unit one.
    let kr = element_stiffness_matrix((0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (0.0, 1.0), e, nu);
    assert!((kr[[0, 0]] - k1[[0, 0]]).abs() > 1.0);
}

#[test]
fn uniform_strain_energy_equals_half_c_times_area_on_a_trapezoid() {
    // u = (x, 0) is in the bilinear space: strain (1, 0, 0), energy = 1/2 c11 A.
    let (e, nu) = (7.0, 0.2);
    let c11 = e / (1.0 - nu * nu);
    let quad = [(0.0, 0.0), (4.0, 0.0), (3.0, 2.0), (1.0, 2.0)];
    let area = 0.5 * (4.0 + 2.0) * 2.0;
    let k = element_stiffness_matrix(quad[0], quad[1], quad[2], quad[3], e, nu);
    let mut u = ndarray::Array1::<f64>::zeros(8);
    for n in 0..4 {
        u[2 * n] = quad[n].0;
    }
    let energy = 0.5 * u.dot(&k.dot(&u));
    assert!(
        (energy - 0.5 * c11 * area).abs() < 1e-9 * energy.abs(),
        "{energy}"
    );
    // Rigid rotation u = (-y, x) carries no energy.
    let mut r = ndarray::Array1::<f64>::zeros(8);
    for n in 0..4 {
        r[2 * n] = -quad[n].1;
        r[2 * n + 1] = quad[n].0;
    }
    assert!(k.dot(&r).iter().all(|v| v.abs() < 1e-8));
}

#[test]
fn cantilever_scenario_writes_into_the_given_directory() {
    let dir = std::env::temp_dir().join(format!("rssn_beam_{}", std::process::id()));
    simulate_cantilever_beam_scenario(&dir).unwrap_or_else(|e| panic!("{e}"));
    assert!(dir.join("beam_original.csv").is_file());
    assert!(dir.join("beam_deformed.csv").is_file());
    let _ = std::fs::remove_dir_all(&dir);
}

fn one_element(
    fixed: Vec<usize>,
    loads: Vec<(usize, f64, f64)>,
) -> ElasticityParameters {
    ElasticityParameters {
        nodes: SQ.to_vec(),
        elements: vec![[0, 1, 2, 3]],
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes: fixed,
        loads,
    }
}

#[test]
fn single_element_pulled_to_the_right() {
    let p = one_element(vec![0, 3], vec![(1, 1000.0, 0.0), (2, 1000.0, 0.0)]);
    let d = run_elasticity_simulation(&p).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(d.len(), 8);
    for i in [0, 1, 6, 7] {
        assert_eq!(d[i], 0.0, "fixed dof {i}");
    }
    assert!(d[2] > 0.0 && d[4] > 0.0);
    // Symmetric loading about y = 1/2 gives equal x-displacements on both loaded nodes.
    assert!((d[2] - d[4]).abs() < 1e-9 * d[2].abs());
}

#[test]
fn displacement_is_linear_in_the_load() {
    let d1 = run_elasticity_simulation(&one_element(
        vec![0, 3],
        vec![(1, 100.0, 0.0), (2, 100.0, 0.0)],
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    let d3 = run_elasticity_simulation(&one_element(
        vec![0, 3],
        vec![(1, 300.0, 0.0), (2, 300.0, 0.0)],
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    for (a, b) in d1.iter().zip(&d3) {
        assert!((3.0 * a - b).abs() < 1e-9 * (1.0 + b.abs()));
    }
}

#[test]
fn a_softer_material_deflects_more() {
    let mut soft = one_element(vec![0, 3], vec![(1, 1000.0, 0.0), (2, 1000.0, 0.0)]);
    let stiff = soft.clone();
    soft.youngs_modulus = 1e6;
    let ds = run_elasticity_simulation(&soft).unwrap_or_else(|e| panic!("{e}"));
    let dk = run_elasticity_simulation(&stiff).unwrap_or_else(|e| panic!("{e}"));
    assert!((ds[2] / dk[2] - 10.0).abs() < 1e-6);
}

#[test]
fn cantilever_tip_deflects_in_the_load_direction() {
    // Reduced version of `simulate_cantilever_beam_scenario` that writes no files.
    let (nx, ny) = (8usize, 2usize);
    let mut nodes: Nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((i as f64 * 4.0 / nx as f64, j as f64 * 1.0 / ny as f64));
        }
    }
    let mut elements: Elements = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let n1 = j * (nx + 1) + i;
            elements.push([n1, n1 + 1, n1 + nx + 2, n1 + nx + 1]);
        }
    }
    let fixed_nodes: Vec<usize> = (0..=ny).map(|j| j * (nx + 1)).collect();
    let tip = (ny / 2) * (nx + 1) + nx;
    let params = ElasticityParameters {
        nodes,
        elements,
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes: fixed_nodes.clone(),
        loads: vec![(tip, 0.0, -1e3)],
    };
    let d = run_elasticity_simulation(&params).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(d.len(), 2 * (nx + 1) * (ny + 1));
    for n in fixed_nodes {
        assert_eq!((d[2 * n], d[2 * n + 1]), (0.0, 0.0));
    }
    assert!(d.iter().all(|v| v.is_finite()));
    assert!(d[2 * tip + 1] < 0.0);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_stiffness_scales_with_youngs_modulus(scale in 0.1f64..10.0) {
        let k1 = k_of(1e7, 0.3);
        let k2 = k_of(1e7 * scale, 0.3);
        for i in 0..8 {
            for j in 0..8 {
                prop_assert!((k1[[i, j]] * scale - k2[[i, j]]).abs() < 1e-5);
            }
        }
    }
}
