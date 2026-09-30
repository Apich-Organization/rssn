//! Finite-element building blocks (ported from `numerical_physics_fea_test.rs`).
//!
//! Tests for materials, elements, stress analysis, and mesh generation.

use std::f64::consts::PI;

use rssn::sim::physics_fea::*;

// ============================================================================
// Material Tests
// ============================================================================

#[test]

fn test_material_new() {
    let mat = Material::new(200e9, 0.3, 7850.0, 50.0, 12e-6, 250e6);

    assert_eq!(mat.youngs_modulus, 200e9);

    assert_eq!(mat.poissons_ratio, 0.3);
}

#[test]

fn test_material_steel() {
    let steel = Material::steel();

    assert_eq!(steel.youngs_modulus, 200e9);

    assert_eq!(steel.poissons_ratio, 0.3);

    assert_eq!(steel.density, 7850.0);
}

#[test]

fn test_material_aluminum() {
    let al = Material::aluminum();

    assert_eq!(al.youngs_modulus, 70e9);

    assert!((al.poissons_ratio - 0.33).abs() < 1e-10);
}

#[test]

fn test_material_copper() {
    let cu = Material::copper();

    assert_eq!(cu.youngs_modulus, 117e9);

    assert!((cu.poissons_ratio - 0.34).abs() < 1e-10);
}

#[test]

fn test_shear_modulus() {
    let steel = Material::steel();

    // G = E / (2(1+ν)) = 200e9 / (2 * 1.3) ≈ 76.9 GPa
    let g = steel.shear_modulus();

    assert!(g > 76e9 && g < 78e9);
}

#[test]

fn test_bulk_modulus() {
    let steel = Material::steel();

    // K = E / (3(1-2ν)) = 200e9 / (3 * 0.4) ≈ 166.7 GPa
    let k = steel.bulk_modulus();

    assert!(k > 165e9 && k < 168e9);
}

// ============================================================================
// Node Tests
// ============================================================================

#[test]

fn test_node2d_new() {
    let node = Node2D::new(0, 1.0, 2.0);

    assert_eq!(node.id, 0);

    assert_eq!(node.x, 1.0);

    assert_eq!(node.y, 2.0);
}

#[test]

fn test_node2d_distance() {
    let n1 = Node2D::new(0, 0.0, 0.0);

    let n2 = Node2D::new(1, 3.0, 4.0);

    assert!((n1.distance_to(&n2) - 5.0).abs() < 1e-10);
}

#[test]

fn test_node3d_new() {
    let node = Node3D::new(0, 1.0, 2.0, 3.0);

    assert_eq!(node.id, 0);

    assert_eq!(node.x, 1.0);

    assert_eq!(node.y, 2.0);

    assert_eq!(node.z, 3.0);
}

// ============================================================================
// 1D Element Tests
// ============================================================================

#[test]

fn test_linear_element_1d() {
    let elem = LinearElement1D {
        length: 1.0,
        youngs_modulus: 200e9,
        area: 0.001,
    };

    let k = elem.local_stiffness_matrix();

    // k = EA/L = 200e9 * 0.001 / 1.0 = 200e6
    assert_eq!(k.rows(), 2);

    assert_eq!(k.cols(), 2);

    assert!((*k.get(0, 0) - 200e6).abs() < 1e-6);

    assert!((*k.get(0, 1) + 200e6).abs() < 1e-6);
}

// ============================================================================
// 2D Triangle Element Tests
// ============================================================================

#[test]

fn test_triangle_element_area() {
    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        0.01,
        Material::steel(),
        true,
    );

    assert!((elem.area() - 0.5).abs() < 1e-10);
}

#[test]

fn test_triangle_element_constitutive_matrix_plane_stress() {
    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        0.01,
        Material::steel(),
        true,
    );

    let d = elem.constitutive_matrix();

    assert_eq!(d.rows(), 3);

    assert_eq!(d.cols(), 3);

    // D11 = E/(1-ν²) = 200e9 / (1-0.09) = 219.78e9
    assert!(*d.get(0, 0) > 200e9);
}

#[test]

fn test_triangle_element_b_matrix() {
    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        0.01,
        Material::steel(),
        true,
    );

    let b = elem.b_matrix();

    assert_eq!(b.rows(), 3);

    assert_eq!(b.cols(), 6);
}

#[test]

fn test_triangle_element_stiffness_matrix() {
    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        0.01,
        Material::steel(),
        true,
    );

    let k = elem.local_stiffness_matrix();

    assert_eq!(k.rows(), 6);

    assert_eq!(k.cols(), 6);

    // Stiffness matrix should be symmetric
    for i in 0..6 {
        for j in 0..6 {
            assert!((*k.get(i, j) - *k.get(j, i)).abs() < 1e-6);
        }
    }
}

#[test]

fn test_von_mises_stress() {
    // Pure tension in x
    let stress = [100e6, 0.0, 0.0];

    let vm = TriangleElement2D::von_mises_stress(&stress);

    assert!((vm - 100e6).abs() < 1e-6);

    // Pure shear
    let stress = [0.0, 0.0, 100e6];

    let vm = TriangleElement2D::von_mises_stress(&stress);

    // von Mises for pure shear = sqrt(3) * τ
    assert!((vm - 3.0_f64.sqrt() * 100e6).abs() < 1e-6);
}

// ============================================================================
// Beam Element Tests
// ============================================================================

#[test]

fn test_beam_element_2d_new() {
    let beam = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, 0.0);

    assert_eq!(beam.length, 1.0);

    assert_eq!(beam.youngs_modulus, 200e9);
}

#[test]

fn test_beam_element_stiffness_matrix() {
    let beam = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, 0.0);

    let k = beam.local_stiffness_matrix();

    assert_eq!(k.rows(), 6);

    assert_eq!(k.cols(), 6);

    // Stiffness matrix should be symmetric
    for i in 0..6 {
        for j in 0..6 {
            assert!((*k.get(i, j) - *k.get(j, i)).abs() < 1.0);
        }
    }
}

#[test]

fn test_beam_element_transformation_matrix() {
    // Horizontal beam (angle = 0)
    let beam = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, 0.0);

    let t = beam.transformation_matrix();

    // Should be identity for cos(0)=1, sin(0)=0
    assert!((*t.get(0, 0) - 1.0).abs() < 1e-10);

    // 90 degree beam
    let beam90 = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, PI / 2.0);

    let t90 = beam90.transformation_matrix();

    // cos(90°)=0, sin(90°)=1
    assert!(t90.get(0, 0).abs() < 1e-10);

    assert!((*t90.get(0, 1) - 1.0).abs() < 1e-10);
}

#[test]

fn test_beam_element_mass_matrix() {
    let beam = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, 0.0);

    let m = beam.mass_matrix(7850.0);

    assert_eq!(m.rows(), 6);

    assert_eq!(m.cols(), 6);

    // Mass matrix should be symmetric
    for i in 0..6 {
        for j in 0..6 {
            assert!((*m.get(i, j) - *m.get(j, i)).abs() < 1e-10);
        }
    }
}

// ============================================================================
// Thermal Element Tests
// ============================================================================

#[test]

fn test_thermal_element_1d() {
    let elem = ThermalElement1D::new(1.0, 50.0, 0.001);

    let k = elem.conductivity_matrix();

    // k = κA/L = 50 * 0.001 / 1.0 = 0.05
    assert!((*k.get(0, 0) - 0.05).abs() < 1e-10);

    assert!((*k.get(0, 1) + 0.05).abs() < 1e-10);
}

#[test]

fn test_thermal_triangle_2d() {
    let elem = ThermalTriangle2D::new([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 0.01, 50.0);

    assert!((elem.area() - 0.5).abs() < 1e-10);

    let k = elem.conductivity_matrix();

    assert_eq!(k.rows(), 3);

    assert_eq!(k.cols(), 3);
}

// ============================================================================
// Stress Analysis Tests
// ============================================================================

#[test]

fn test_principal_stresses_uniaxial() {
    // Pure tension in x
    let (s1, s2, angle) = principal_stresses(&[100e6, 0.0, 0.0]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!(s2.abs() < 1e-6);

    assert!(angle.abs() < 1e-10);
}

#[test]

fn test_principal_stresses_biaxial() {
    // Equal biaxial tension
    let (s1, s2, _angle) = principal_stresses(&[100e6, 100e6, 0.0]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!((s2 - 100e6).abs() < 1e-6);
}

#[test]

fn test_principal_stresses_pure_shear() {
    // Pure shear
    let (s1, s2, _angle) = principal_stresses(&[0.0, 0.0, 100e6]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!((s2 + 100e6).abs() < 1e-6);
}

#[test]

fn test_max_shear_stress() {
    let tau = max_shear_stress(100e6, -100e6);

    assert!((tau - 100e6).abs() < 1e-6);
}

#[test]

fn test_safety_factor_von_mises() {
    let stress = [100e6, 0.0, 0.0];

    let sf = safety_factor_von_mises(&stress, 250e6);

    assert!((sf - 2.5).abs() < 1e-6);
}

// ============================================================================
// Mesh Generation Tests
// ============================================================================

#[test]

fn test_create_rectangular_mesh() {
    let (nodes, elements) = create_rectangular_mesh(1.0, 1.0, 2, 2);

    // 3x3 = 9 nodes
    assert_eq!(nodes.len(), 9);

    // 2x2 rectangles × 2 triangles each = 8 triangles
    assert_eq!(elements.len(), 8);
}

#[test]

fn test_create_rectangular_mesh_nodes() {
    let (nodes, _elements) = create_rectangular_mesh(2.0, 1.0, 2, 1);

    // 3×2 = 6 nodes
    assert_eq!(nodes.len(), 6);

    // Check corner nodes
    assert!(nodes[0].x.abs() < 1e-10);

    assert!(nodes[0].y.abs() < 1e-10);

    assert!((nodes[2].x - 2.0).abs() < 1e-10);
}

#[test]

fn test_refine_mesh() {
    let (nodes, elements) = create_rectangular_mesh(1.0, 1.0, 1, 1);

    // Initial: 4 nodes, 2 triangles
    assert_eq!(nodes.len(), 4);

    assert_eq!(elements.len(), 2);

    let (new_nodes, new_elements) = refine_mesh(&nodes, &elements);

    // Each triangle becomes 4 triangles
    assert_eq!(new_elements.len(), 8);

    // New nodes from edge midpoints
    assert!(new_nodes.len() > 4);
}

// ============================================================================
// Assembly Tests
// ============================================================================

#[test]

fn test_assemble_global_stiffness_matrix() {
    let elem1 = LinearElement1D {
        length: 1.0,
        youngs_modulus: 200e9,
        area: 0.001,
    };

    let k1 = elem1.local_stiffness_matrix();

    let elements = vec![(k1.clone(), 0, 1), (k1.clone(), 1, 2)];

    let global_k = assemble_global_stiffness_matrix(3, &elements);

    assert_eq!(global_k.rows(), 3);

    assert_eq!(global_k.cols(), 3);

    // Interior node (node 1) should have contribution from both elements
    assert!((*global_k.get(1, 1) - 2.0 * 200e6).abs() < 1e-6);
}

#[test]

fn test_solve_static_structural() {
    // Simple 2-element bar
    let elem = LinearElement1D {
        length: 1.0,
        youngs_modulus: 200e9,
        area: 0.001,
    };

    let k_local = elem.local_stiffness_matrix();

    let elements = vec![(k_local.clone(), 0, 1), (k_local.clone(), 1, 2)];

    let global_k = assemble_global_stiffness_matrix(3, &elements);

    // Fixed at node 0, force at node 2
    let forces = vec![0.0, 0.0, 1000.0];

    let fixed_dofs = vec![(0, 0.0)];

    let result = solve_static_structural(global_k, forces, &fixed_dofs);

    assert!(result.is_ok());

    let u = result.unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(u.len(), 3);

    assert!(u[0].abs() < 1e-20); // Fixed
    assert!(u[2] > u[1]); // Displacement increases toward load
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use proptest::prelude::*;

    use super::*;

    proptest! {
        #[test]
        fn prop_von_mises_non_negative(sx in -1e9..1e9f64, sy in -1e9..1e9f64, txy in -1e9..1e9f64) {
            let vm = TriangleElement2D::von_mises_stress(&[sx, sy, txy]);
            prop_assert!(vm >= 0.0);
        }

        #[test]
        fn prop_shear_modulus_positive(e in 1e6..1e12f64, nu in 0.0..0.49f64) {
            let mat = Material::new(e, nu, 1000.0, 1.0, 1e-6, 1e6);
            prop_assert!(mat.shear_modulus() > 0.0);
        }

        #[test]
        fn prop_triangle_area_positive(
            x1 in 0.0..10.0f64, y1 in 0.0..10.0f64,
            x2 in 0.0..10.0f64, y2 in 0.0..10.0f64,
            x3 in 0.0..10.0f64, y3 in 0.0..10.0f64
        ) {
            let elem = TriangleElement2D::new(
                [0, 1, 2],
                [(x1, y1), (x2, y2), (x3, y3)],
                0.01,
                Material::steel(),
                true,
            );
            prop_assert!(elem.area() >= 0.0);
        }

        #[test]
        fn prop_principal_stress_ordering(sx in -1e9..1e9f64, sy in -1e9..1e9f64, txy in -1e9..1e9f64) {
            let (s1, s2, _) = principal_stresses(&[sx, sy, txy]);
            prop_assert!(s1 >= s2);
        }
    }
}

// ============================================================================
// Added: exact structural-mechanics results
// ============================================================================

mod strengthened {
    use std::f64::consts::PI;

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::matrix::Matrix;
    use rssn::sim::physics_fea::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            cases: 48,
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-30)
    }

    fn mat_vec(k: &Matrix<f64>, x: &[f64]) -> Vec<f64> {
        (0..k.rows()).map(|i| (0..k.cols()).map(|j| k.get(i, j) * x[j]).sum()).collect()
    }

    fn tri(coords: [(f64, f64); 3], plane_stress: bool) -> TriangleElement2D {
        TriangleElement2D::new([0, 1, 2], coords, 0.01, Material::steel(), plane_stress)
    }

    #[test]
    fn material_relations() {
        for m in [Material::steel(), Material::aluminum(), Material::copper()] {
            let (e, nu) = (m.youngs_modulus, m.poissons_ratio);
            let (g, k) = (m.shear_modulus(), m.bulk_modulus());
            assert!(close(g, e / (2.0 * (1.0 + nu)), 1e-12));
            // E = 9 K G / (3 K + G) and nu = (3K - 2G) / (2 (3K + G))
            assert!(close(e, 9.0 * k * g / (3.0 * k + g), 1e-12));
            assert!(close(nu, (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g)), 1e-12));
        }
        let steel = Material::steel();
        assert!(close(steel.shear_modulus(), 76.923_076_923e9, 1e-9));
        assert!(close(steel.bulk_modulus(), 166.666_666_667e9, 1e-9));
        assert_eq!((steel.yield_strength, steel.thermal_conductivity, steel.thermal_expansion), (250e6, 50.0, 12e-6));
    }

    #[test]
    fn node_helpers() {
        assert!((Node2D::new(1, 1.0, 1.0).distance_to(&Node2D::new(2, 4.0, 5.0)) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn linear_bar_stiffness_is_ea_over_l_with_unit_row_sums_zero() {
        let k = LinearElement1D { length: 2.0, youngs_modulus: 100e9, area: 0.002 }.local_stiffness_matrix();
        let ea_l = 100e9 * 0.002 / 2.0;
        assert_eq!(k.data(), &vec![ea_l, -ea_l, -ea_l, ea_l]);
    }

    #[test]
    fn constitutive_matrices_match_textbook_formulas() {
        let (e, nu) = (200e9, 0.3);
        let ps = tri([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], true).constitutive_matrix();
        let f = e / (1.0 - nu * nu);
        for (idx, want) in [f, f * nu, 0.0, f * nu, f, 0.0, 0.0, 0.0, f * (1.0 - nu) / 2.0].iter().enumerate() {
            assert!(close(ps.data()[idx], *want, 1e-12) || (ps.data()[idx] == 0.0 && *want == 0.0), "plane stress entry {idx}");
        }
        let pe = tri([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], false).constitutive_matrix();
        let f = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        assert!(close(*pe.get(0, 0), f * (1.0 - nu), 1e-12));
        assert!(close(*pe.get(0, 1), f * nu, 1e-12));
        assert!(close(*pe.get(2, 2), Material::steel().shear_modulus(), 1e-12), "D33 must equal G in both formulations");
        assert!(close(*ps.get(2, 2), Material::steel().shear_modulus(), 1e-12));
    }

    #[test]
    fn b_matrix_of_the_unit_right_triangle() {
        let b = tri([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], true).b_matrix();
        // shape function gradients: N1 = 1 - x - y, N2 = x, N3 = y
        let want = [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.0];
        assert_eq!(b.data().len(), 18);
        for (g, w) in b.data().iter().zip(want) {
            assert!((g - w).abs() < 1e-12, "{:?}", b.data());
        }
    }

    #[test]
    fn cst_reproduces_constant_strain_fields_exactly() {
        let coords = [(0.1, 0.2), (1.3, 0.1), (0.4, 0.9)];
        let el = tri(coords, true);
        // u = a x + b y + c, v = d x + e y + f
        let (a, b, c, d, e, f) = (1e-4, 2e-5, 3e-6, -4e-5, 5e-5, 7e-6);
        let u: Vec<f64> = coords.iter().flat_map(|&(x, y)| [a * x + b * y + c, d * x + e * y + f]).collect();
        let strain = compute_element_strain(&el.b_matrix(), &u);
        assert!(close(strain[0], a, 1e-9) && close(strain[1], e, 1e-9) && close(strain[2], b + d, 1e-9), "{strain:?}");
        let stress = el.compute_stress(&u);
        let dm = el.constitutive_matrix();
        for i in 0..3 {
            let want: f64 = (0..3).map(|j| dm.get(i, j) * strain[j]).sum();
            assert!(close(stress[i], want, 1e-9));
        }
        // sigma_x = E/(1-nu^2) (eps_x + nu eps_y)
        assert!(close(stress[0], 200e9 / (1.0 - 0.09) * (a + 0.3 * e), 1e-9));
    }

    #[test]
    #[should_panic(expected = "Need 6 displacement")]
    fn compute_stress_requires_six_displacements() {
        let _ = tri([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], true).compute_stress(&[0.0; 4]);
    }

    #[test]
    fn cst_stiffness_annihilates_rigid_body_motions_and_stores_the_strain_energy() {
        let coords = [(0.0, 0.0), (2.0, 0.0), (0.5, 1.5)];
        let el = tri(coords, true);
        let k = el.local_stiffness_matrix();
        let scale = k.data().iter().fold(0.0f64, |m, v| m.max(v.abs()));
        // Translations in x and y.
        for tv in [[1.0, 0.0], [0.0, 1.0]] {
            let u: Vec<f64> = (0..3).flat_map(|_| tv).collect();
            assert!(mat_vec(&k, &u).iter().all(|f| f.abs() < 1e-9 * scale));
        }
        // Small rigid rotation u = -theta y, v = theta x.
        let u: Vec<f64> = coords.iter().flat_map(|&(x, y)| [-y, x]).collect();
        assert!(mat_vec(&k, &u).iter().all(|f| f.abs() < 1e-9 * scale));
        // Uniform axial strain eps: energy = 1/2 u^T K u = 1/2 t A eps^T D eps
        let eps = 1e-3;
        let u: Vec<f64> = coords.iter().flat_map(|&(x, _)| [eps * x, 0.0]).collect();
        let ku = mat_vec(&k, &u);
        let energy = 0.5 * u.iter().zip(&ku).map(|(a, b)| a * b).sum::<f64>();
        let d11 = el.constitutive_matrix().get(0, 0).to_owned();
        assert!(close(energy, 0.5 * el.thickness * el.area() * d11 * eps * eps, 1e-9));
        // Symmetric positive semi-definite.
        assert!(energy > 0.0);
    }

    #[test]
    fn two_triangle_plate_under_uniform_strain_is_in_equilibrium() {
        let (nodes, elems) = create_rectangular_mesh(2.0, 1.0, 1, 1);
        let mat = Material::steel();
        let t = 0.01;
        let mut list = Vec::new();
        for e in &elems {
            let coords = [
                (nodes[e[0]].x, nodes[e[0]].y),
                (nodes[e[1]].x, nodes[e[1]].y),
                (nodes[e[2]].x, nodes[e[2]].y),
            ];
            let el = TriangleElement2D::new(*e, coords, t, mat, true);
            let dofs = [2 * e[0], 2 * e[0] + 1, 2 * e[1], 2 * e[1] + 1, 2 * e[2], 2 * e[2] + 1];
            list.push((el.local_stiffness_matrix(), dofs));
        }
        let k = assemble_2d_stiffness_matrix(2 * nodes.len(), &list);
        assert_eq!((k.rows(), k.cols()), (8, 8));
        let eps = 1e-4;
        let u: Vec<f64> = nodes.iter().flat_map(|n| [eps * n.x, 0.0]).collect();
        let f = mat_vec(&k, &u);
        // Net force vanishes, and the total axial force is sigma_x * t * height on each end.
        assert!(f.iter().step_by(2).sum::<f64>().abs() < 1e-3);
        assert!(f.iter().skip(1).step_by(2).sum::<f64>().abs() < 1e-3);
        let sigma_x = mat.youngs_modulus / (1.0 - 0.09) * eps;
        let right_end: f64 = nodes.iter().filter(|n| (n.x - 2.0).abs() < 1e-12).map(|n| f[2 * n.id]).sum();
        assert!(close(right_end, sigma_x * t * 1.0, 1e-9), "{right_end} vs {}", sigma_x * t);
        // Energy 1/2 u^T K u = 1/2 sigma_x eps V.
        let energy = 0.5 * u.iter().zip(&f).map(|(a, b)| a * b).sum::<f64>();
        assert!(close(energy, 0.5 * sigma_x * eps * 2.0 * t, 1e-9));
    }

    #[test]
    fn stress_invariants() {
        let vm = |s: [f64; 3]| TriangleElement2D::von_mises_stress(&s);
        assert!(close(vm([100e6, 100e6, 0.0]), 100e6, 1e-12), "equibiaxial: sigma_vm = sigma");
        assert!(close(vm([100e6, -100e6, 0.0]), 3f64.sqrt() * 100e6, 1e-12));
        let (s1, s2, ang) = principal_stresses(&[80e6, 20e6, 40e6]);
        // centre 50, radius sqrt(30^2 + 40^2) = 50
        assert!((s1 - 100e6).abs() < 1e-3 && s2.abs() < 1e-3, "s1 = {s1}, s2 = {s2}");
        assert!((ang - 0.5 * (40.0f64).atan2(30.0)).abs() < 1e-12);
        assert!(close(max_shear_stress(s1, s2), 50e6, 1e-9));
        // Invariants: s1 + s2 = sx + sy ; s1 s2 = sx sy - txy^2
        assert!((s1 + s2 - 100e6).abs() < 1.0);
        assert!((s1 * s2 - (80e6 * 20e6 - 40e6 * 40e6)).abs() < 1e3);
        assert_eq!(safety_factor_von_mises(&[0.0, 0.0, 0.0], 250e6), f64::INFINITY);
        assert!((safety_factor_von_mises(&[0.0, 0.0, 100e6], 250e6) - 250.0 / (3f64.sqrt() * 100.0)).abs() < 1e-9);
    }

    #[test]
    fn beam_local_stiffness_reference_entries_and_symmetry() {
        let (e, a, i, l) = (200e9, 0.001, 1e-6, 2.0);
        let k = BeamElement2D::new(l, e, a, i, 0.0).local_stiffness_matrix();
        assert!(close(*k.get(0, 0), e * a / l, 1e-12));
        assert!(close(*k.get(1, 1), 12.0 * e * i / l.powi(3), 1e-12));
        assert!(close(*k.get(1, 2), 6.0 * e * i / (l * l), 1e-12));
        assert!(close(*k.get(2, 2), 4.0 * e * i / l, 1e-12));
        assert!(close(*k.get(2, 5), 2.0 * e * i / l, 1e-12));
        for r in 0..6 {
            for c in 0..6 {
                assert!(close(*k.get(r, c), *k.get(c, r), 1e-12) || (k.get(r, c) - k.get(c, r)).abs() < 1e-9);
            }
        }
        // Translation and rigid rotation produce no force: u = [1,0,0,1,0,0], [0,1,0,0,1,0], [0,0,1,0,l... ] etc.
        let scale = k.data().iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for u in [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0, l, 1.0]] {
            assert!(mat_vec(&k, &u).iter().all(|f| f.abs() < 1e-9 * scale), "{u:?}");
        }
    }

    #[test]
    fn beam_transformation_is_orthogonal_and_rotates_the_stiffness() {
        let beam = BeamElement2D::new(1.5, 70e9, 0.002, 3e-6, PI / 3.0);
        let t = beam.transformation_matrix();
        assert!((t.clone() * t.transpose()).is_identity(1e-12), "T T^T = I");
        let kg = beam.global_stiffness_matrix();
        let kl = beam.local_stiffness_matrix();
        // Rotation preserves the Frobenius norm and the trace.
        assert!(close(kg.frobenius_norm(), kl.frobenius_norm(), 1e-12));
        assert!(close(kg.trace().unwrap_or(f64::NAN), kl.trace().unwrap_or(f64::NAN), 1e-12));
        assert!(kg.is_symmetric());
        // A vertical bar (90 degrees) has its axial stiffness on the y dof.
        let vk = BeamElement2D::new(1.0, 200e9, 0.001, 1e-6, PI / 2.0).global_stiffness_matrix();
        assert!(close(*vk.get(1, 1), 200e9 * 0.001, 1e-9));
        // ...and the bending stiffness 12 E I / L^3 acts along global x.
        assert!(close(*vk.get(0, 0), 12.0 * 200e9 * 1e-6, 1e-9));
    }

    #[test]
    fn beam_mass_matrix_conserves_total_mass() {
        let (rho, a, l) = (7850.0, 0.001, 2.0);
        let m = BeamElement2D::new(l, 200e9, a, 1e-6, 0.0).mass_matrix(rho);
        let total = rho * a * l;
        // Rigid translation in y: u^T M u = total mass (and likewise for x).
        let uy = [0.0, 1.0, 0.0, 0.0, 1.0, 0.0];
        let ux = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        for u in [uy, ux] {
            let mu = mat_vec(&m, &u);
            assert!(close(u.iter().zip(&mu).map(|(a, b)| a * b).sum::<f64>(), total, 1e-12));
        }
    }

    #[test]
    fn cantilever_tip_load_matches_euler_bernoulli() {
        // One beam element is exact for a tip load: delta = P L^3 / (3 E I), theta = P L^2 / (2 E I).
        let (e, a, i, l, p) = (200e9, 0.001, 8.0e-6, 2.5, 5000.0);
        let k = BeamElement2D::new(l, e, a, i, 0.0).global_stiffness_matrix();
        let mut f = vec![0.0; 6];
        f[4] = p;
        let u = solve_static_structural(k, f, &[(0, 0.0), (1, 0.0), (2, 0.0)]).unwrap_or_else(|e| panic!("{e}"));
        assert!(close(u[4], p * l.powi(3) / (3.0 * e * i), 1e-9), "tip deflection {}", u[4]);
        assert!(close(u[5], p * l * l / (2.0 * e * i), 1e-9), "tip rotation {}", u[5]);
        assert!(u[3].abs() < 1e-15, "no axial extension for a transverse load");
    }

    #[test]
    fn bar_chain_displacements_and_prescribed_motion() {
        let k = LinearElement1D { length: 1.0, youngs_modulus: 200e9, area: 0.001 }.local_stiffness_matrix();
        let kk = 200e6;
        let global = assemble_global_stiffness_matrix(3, &[(k.clone(), 0, 1), (k.clone(), 1, 2)]);
        assert_eq!(global.data(), &vec![kk, -kk, 0.0, -kk, 2.0 * kk, -kk, 0.0, -kk, kk]);
        let u = solve_static_structural(global.clone(), vec![0.0, 0.0, 1000.0], &[(0, 0.0)]).unwrap_or_else(|e| panic!("{e}"));
        assert!(close(u[1], 1000.0 / kk, 1e-9) && close(u[2], 2000.0 / kk, 1e-9), "{u:?}");
        // Prescribed end displacement moves the whole unloaded bar rigidly.
        let u = solve_static_structural(global.clone(), vec![0.0; 3], &[(0, 0.001)]).unwrap_or_else(|e| panic!("{e}"));
        assert!(u.iter().all(|&v| (v - 0.001).abs() < 1e-15), "{u:?}");
        // Prescribing both ends fixes the middle node by symmetry.
        let u = solve_static_structural(global.clone(), vec![0.0; 3], &[(0, 0.0), (2, 0.002)]).unwrap_or_else(|e| panic!("{e}"));
        assert!((u[1] - 0.001).abs() < 1e-12, "{u:?}");
        // Errors: wrong force length, and an unconstrained (singular) system.
        assert!(solve_static_structural(global.clone(), vec![0.0; 2], &[(0, 0.0)]).is_err());
        assert!(solve_static_structural(global, vec![1.0, 0.0, 0.0], &[]).is_err());
    }

    #[test]
    fn penalty_boundary_conditions_approximate_the_exact_solution() {
        let k = LinearElement1D { length: 1.0, youngs_modulus: 200e9, area: 0.001 }.local_stiffness_matrix();
        let mut global = assemble_global_stiffness_matrix(3, &[(k.clone(), 0, 1), (k, 1, 2)]);
        let mut f = vec![0.0, 0.0, 1000.0];
        apply_boundary_conditions_penalty(&mut global, &mut f, &[(0, 0.0)], 1e20);
        let sol = solve_static_structural(global, f, &[]).unwrap_or_else(|e| panic!("{e}"));
        assert!((sol[1] - 5e-6).abs() < 1e-9 && (sol[2] - 1e-5).abs() < 1e-9, "{sol:?}");
    }

    #[test]
    fn thermal_elements() {
        let k1 = ThermalElement1D::new(2.0, 10.0, 0.5).conductivity_matrix();
        assert_eq!(k1.data(), &vec![2.5, -2.5, -2.5, 2.5]);
        let tri = ThermalTriangle2D::new([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 0.1, 4.0);
        let k = tri.conductivity_matrix();
        // k t / (4 A) = 4 * 0.1 / 2 = 0.2 ; grad N1 = (-1, -1), N2 = (1, 0), N3 = (0, 1)
        let want = [0.4, -0.2, -0.2, -0.2, 0.2, 0.0, -0.2, 0.0, 0.2];
        for (g, w) in k.data().iter().zip(want) {
            assert!((g - w).abs() < 1e-12, "{:?}", k.data());
        }
        // A uniform temperature produces no heat flow: rows sum to zero.
        for r in 0..3 {
            assert!((0..3).map(|c| k.get(r, c)).sum::<f64>().abs() < 1e-12);
        }
    }

    #[test]
    fn rectangular_mesh_geometry() {
        let (nodes, elems) = create_rectangular_mesh(3.0, 2.0, 3, 2);
        assert_eq!((nodes.len(), elems.len()), (12, 12));
        assert!(nodes.iter().enumerate().all(|(i, n)| n.id == i));
        assert_eq!((nodes[11].x, nodes[11].y), (3.0, 2.0));
        let area = |e: &[usize; 3]| {
            let (a, b, c) = (&nodes[e[0]], &nodes[e[1]], &nodes[e[2]]);
            0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y))
        };
        // Counter-clockwise triangles that tile the rectangle.
        assert!(elems.iter().all(|e| area(e) > 0.0));
        assert!((elems.iter().map(area).sum::<f64>() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn mesh_refinement_preserves_area_and_shares_edge_midpoints() {
        let (nodes, elems) = create_rectangular_mesh(1.0, 1.0, 1, 1);
        let (n2, e2) = refine_mesh(&nodes, &elems);
        // 4 corners + 5 distinct edge midpoints (the diagonal is shared).
        assert_eq!((n2.len(), e2.len()), (9, 8));
        assert!(n2.iter().enumerate().all(|(i, n)| n.id == i));
        let area = |e: &[usize; 3]| {
            let (a, b, c) = (&n2[e[0]], &n2[e[1]], &n2[e[2]]);
            0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)).abs()
        };
        assert!(e2.iter().all(|e| (area(e) - 0.125).abs() < 1e-12));
        assert!((e2.iter().map(area).sum::<f64>() - 1.0).abs() < 1e-12);
        // The diagonal midpoint is the centre of the square.
        assert!(n2.iter().any(|n| (n.x - 0.5).abs() < 1e-12 && (n.y - 0.5).abs() < 1e-12));
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_principal_stresses_are_invariant_and_ordered(sx in -1e9..1e9f64, sy in -1e9..1e9f64, txy in -1e9..1e9f64) {
            let (s1, s2, _) = principal_stresses(&[sx, sy, txy]);
            let scale = sx.abs().max(sy.abs()).max(txy.abs()).max(1.0);
            prop_assert!((s1 + s2 - (sx + sy)).abs() < 1e-9 * scale);
            prop_assert!((s1 * s2 - (sx * sy - txy * txy)).abs() < 1e-9 * scale * scale);
            // The maximum shear is the Mohr circle radius.
            prop_assert!((max_shear_stress(s1, s2) - (((sx - sy) / 2.0).powi(2) + txy * txy).sqrt()).abs() < 1e-9 * scale);
        }

        #[test]
        fn prop_von_mises_from_principal_stresses(sx in -1e9..1e9f64, sy in -1e9..1e9f64, txy in -1e9..1e9f64) {
            let vm = TriangleElement2D::von_mises_stress(&[sx, sy, txy]);
            let (s1, s2, _) = principal_stresses(&[sx, sy, txy]);
            let from_principal = (s1 * s1 - s1 * s2 + s2 * s2).sqrt();
            prop_assert!((vm - from_principal).abs() <= 1e-9 * vm.max(1.0));
        }

        #[test]
        fn prop_cst_stiffness_is_symmetric_with_zero_row_sums_for_translations(
            x2 in 0.5..3.0f64, x3 in 0.0..3.0f64, y3 in 0.5..3.0f64,
        ) {
            let el = tri([(0.0, 0.0), (x2, 0.0), (x3, y3)], true);
            let k = el.local_stiffness_matrix();
            let scale = k.data().iter().fold(0.0f64, |m, v| m.max(v.abs()));
            for r in 0..6 {
                for c in 0..6 {
                    prop_assert!((k.get(r, c) - k.get(c, r)).abs() < 1e-9 * scale);
                }
                // x-translation and y-translation
                let fx: f64 = (0..3).map(|n| k.get(r, 2 * n)).sum();
                let fy: f64 = (0..3).map(|n| k.get(r, 2 * n + 1)).sum();
                prop_assert!(fx.abs() < 1e-8 * scale && fy.abs() < 1e-8 * scale);
            }
        }

        #[test]
        fn prop_cantilever_deflection_scales_linearly_with_load(p in 100.0..1e5f64, l in 0.5..5.0f64) {
            let (e, i) = (200e9, 1e-5);
            let k = BeamElement2D::new(l, e, 0.01, i, 0.0).global_stiffness_matrix();
            let mut f = vec![0.0; 6];
            f[4] = p;
            let u = solve_static_structural(k, f, &[(0, 0.0), (1, 0.0), (2, 0.0)]).map_err(TestCaseError::fail)?;
            prop_assert!((u[4] - p * l.powi(3) / (3.0 * e * i)).abs() < 1e-8 * u[4].abs());
        }
    }
}
