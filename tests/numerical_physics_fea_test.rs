//! Comprehensive tests for numerical FEA module.
//!
//! Tests for materials, elements, stress analysis, and mesh generation.

use std::f64::consts::PI;

use rssn::numerical::physics_fea::*;

// ============================================================================
// Material Tests
// ============================================================================

#[test]

fn test_material_new() {

    let mat = Material::new(
        200e9,
        0.3,
        7850.0,
        50.0,
        12e-6,
        250e6,
    );

    assert_eq!(
        mat.youngs_modulus,
        200e9
    );

    assert_eq!(
        mat.poissons_ratio,
        0.3
    );
}

#[test]

fn test_material_steel() {

    let steel = Material::steel();

    assert_eq!(
        steel.youngs_modulus,
        200e9
    );

    assert_eq!(
        steel.poissons_ratio,
        0.3
    );

    assert_eq!(
        steel.density,
        7850.0
    );
}

#[test]

fn test_material_aluminum() {

    let al = Material::aluminum();

    assert_eq!(
        al.youngs_modulus,
        70e9
    );

    assert!(
        (al.poissons_ratio - 0.33)
            .abs()
            < 1e-10
    );
}

#[test]

fn test_material_copper() {

    let cu = Material::copper();

    assert_eq!(
        cu.youngs_modulus,
        117e9
    );

    assert!(
        (cu.poissons_ratio - 0.34)
            .abs()
            < 1e-10
    );
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

    assert!(
        (n1.distance_to(&n2) - 5.0)
            .abs()
            < 1e-10
    );
}

#[test]

fn test_node3d_new() {

    let node =
        Node3D::new(0, 1.0, 2.0, 3.0);

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

    let k =
        elem.local_stiffness_matrix();

    // k = EA/L = 200e9 * 0.001 / 1.0 = 200e6
    assert_eq!(k.rows(), 2);

    assert_eq!(k.cols(), 2);

    assert!(
        (*k.get(0, 0) - 200e6).abs()
            < 1e-6
    );

    assert!(
        (*k.get(0, 1) + 200e6).abs()
            < 1e-6
    );
}

// ============================================================================
// 2D Triangle Element Tests
// ============================================================================

#[test]

fn test_triangle_element_area() {

    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ],
        0.01,
        Material::steel(),
        true,
    );

    assert!(
        (elem.area() - 0.5).abs()
            < 1e-10
    );
}

#[test]

fn test_triangle_element_constitutive_matrix_plane_stress(
) {

    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ],
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
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ],
        0.01,
        Material::steel(),
        true,
    );

    let b = elem.b_matrix();

    assert_eq!(b.rows(), 3);

    assert_eq!(b.cols(), 6);
}

#[test]

fn test_triangle_element_stiffness_matrix(
) {

    let elem = TriangleElement2D::new(
        [0, 1, 2],
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ],
        0.01,
        Material::steel(),
        true,
    );

    let k =
        elem.local_stiffness_matrix();

    assert_eq!(k.rows(), 6);

    assert_eq!(k.cols(), 6);

    // Stiffness matrix should be symmetric
    for i in 0 .. 6 {

        for j in 0 .. 6 {

            assert!(
                (*k.get(i, j)
                    - *k.get(j, i))
                .abs()
                    < 1e-6
            );
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
    assert!(
        (vm - 3.0_f64.sqrt() * 100e6)
            .abs()
            < 1e-6
    );
}

// ============================================================================
// Beam Element Tests
// ============================================================================

#[test]

fn test_beam_element_2d_new() {

    let beam = BeamElement2D::new(
        1.0, 200e9, 0.001, 1e-6, 0.0,
    );

    assert_eq!(beam.length, 1.0);

    assert_eq!(
        beam.youngs_modulus,
        200e9
    );
}

#[test]

fn test_beam_element_stiffness_matrix()
{

    let beam = BeamElement2D::new(
        1.0, 200e9, 0.001, 1e-6, 0.0,
    );

    let k =
        beam.local_stiffness_matrix();

    assert_eq!(k.rows(), 6);

    assert_eq!(k.cols(), 6);

    // Stiffness matrix should be symmetric
    for i in 0 .. 6 {

        for j in 0 .. 6 {

            assert!(
                (*k.get(i, j)
                    - *k.get(j, i))
                .abs()
                    < 1.0
            );
        }
    }
}

#[test]

fn test_beam_element_transformation_matrix(
) {

    // Horizontal beam (angle = 0)
    let beam = BeamElement2D::new(
        1.0, 200e9, 0.001, 1e-6, 0.0,
    );

    let t =
        beam.transformation_matrix();

    // Should be identity for cos(0)=1, sin(0)=0
    assert!(
        (*t.get(0, 0) - 1.0).abs()
            < 1e-10
    );

    // 90 degree beam
    let beam90 = BeamElement2D::new(
        1.0,
        200e9,
        0.001,
        1e-6,
        PI / 2.0,
    );

    let t90 =
        beam90.transformation_matrix();

    // cos(90°)=0, sin(90°)=1
    assert!(
        t90.get(0, 0).abs() < 1e-10
    );

    assert!(
        (*t90.get(0, 1) - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_beam_element_mass_matrix() {

    let beam = BeamElement2D::new(
        1.0, 200e9, 0.001, 1e-6, 0.0,
    );

    let m = beam.mass_matrix(7850.0);

    assert_eq!(m.rows(), 6);

    assert_eq!(m.cols(), 6);

    // Mass matrix should be symmetric
    for i in 0 .. 6 {

        for j in 0 .. 6 {

            assert!(
                (*m.get(i, j)
                    - *m.get(j, i))
                .abs()
                    < 1e-10
            );
        }
    }
}

// ============================================================================
// Thermal Element Tests
// ============================================================================

#[test]

fn test_thermal_element_1d() {

    let elem = ThermalElement1D::new(
        1.0, 50.0, 0.001,
    );

    let k = elem.conductivity_matrix();

    // k = κA/L = 50 * 0.001 / 1.0 = 0.05
    assert!(
        (*k.get(0, 0) - 0.05).abs()
            < 1e-10
    );

    assert!(
        (*k.get(0, 1) + 0.05).abs()
            < 1e-10
    );
}

#[test]

fn test_thermal_triangle_2d() {

    let elem = ThermalTriangle2D::new(
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ],
        0.01,
        50.0,
    );

    assert!(
        (elem.area() - 0.5).abs()
            < 1e-10
    );

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
    let (s1, s2, angle) =
        principal_stresses(&[
            100e6, 0.0, 0.0,
        ]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!(s2.abs() < 1e-6);

    assert!(angle.abs() < 1e-10);
}

#[test]

fn test_principal_stresses_biaxial() {

    // Equal biaxial tension
    let (s1, s2, _angle) =
        principal_stresses(&[
            100e6, 100e6, 0.0,
        ]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!((s2 - 100e6).abs() < 1e-6);
}

#[test]

fn test_principal_stresses_pure_shear()
{

    // Pure shear
    let (s1, s2, _angle) =
        principal_stresses(&[
            0.0, 0.0, 100e6,
        ]);

    assert!((s1 - 100e6).abs() < 1e-6);

    assert!((s2 + 100e6).abs() < 1e-6);
}

#[test]

fn test_max_shear_stress() {

    let tau =
        max_shear_stress(100e6, -100e6);

    assert!((tau - 100e6).abs() < 1e-6);
}

#[test]

fn test_safety_factor_von_mises() {

    let stress = [100e6, 0.0, 0.0];

    let sf = safety_factor_von_mises(
        &stress,
        250e6,
    );

    assert!((sf - 2.5).abs() < 1e-6);
}

// ============================================================================
// Mesh Generation Tests
// ============================================================================

#[test]

fn test_create_rectangular_mesh() {

    let (nodes, elements) =
        create_rectangular_mesh(
            1.0, 1.0, 2, 2,
        );

    // 3x3 = 9 nodes
    assert_eq!(nodes.len(), 9);

    // 2x2 rectangles × 2 triangles each = 8 triangles
    assert_eq!(elements.len(), 8);
}

#[test]

fn test_create_rectangular_mesh_nodes()
{

    let (nodes, _elements) =
        create_rectangular_mesh(
            2.0, 1.0, 2, 1,
        );

    // 3×2 = 6 nodes
    assert_eq!(nodes.len(), 6);

    // Check corner nodes
    assert!(nodes[0].x.abs() < 1e-10);

    assert!(nodes[0].y.abs() < 1e-10);

    assert!(
        (nodes[2].x - 2.0).abs()
            < 1e-10
    );
}

#[test]

fn test_refine_mesh() {

    let (nodes, elements) =
        create_rectangular_mesh(
            1.0, 1.0, 1, 1,
        );

    // Initial: 4 nodes, 2 triangles
    assert_eq!(nodes.len(), 4);

    assert_eq!(elements.len(), 2);

    let (new_nodes, new_elements) =
        refine_mesh(&nodes, &elements);

    // Each triangle becomes 4 triangles
    assert_eq!(
        new_elements.len(),
        8
    );

    // New nodes from edge midpoints
    assert!(new_nodes.len() > 4);
}

// ============================================================================
// Assembly Tests
// ============================================================================

#[test]

fn test_assemble_global_stiffness_matrix(
) {

    let elem1 = LinearElement1D {
        length: 1.0,
        youngs_modulus: 200e9,
        area: 0.001,
    };

    let k1 =
        elem1.local_stiffness_matrix();

    let elements = vec![
        (k1.clone(), 0, 1),
        (k1.clone(), 1, 2),
    ];

    let global_k = assemble_global_stiffness_matrix(3, &elements);

    assert_eq!(global_k.rows(), 3);

    assert_eq!(global_k.cols(), 3);

    // Interior node (node 1) should have contribution from both elements
    assert!(
        (*global_k.get(1, 1)
            - 2.0 * 200e6)
            .abs()
            < 1e-6
    );
}

#[test]

fn test_solve_static_structural() {

    // Simple 2-element bar
    let elem = LinearElement1D {
        length: 1.0,
        youngs_modulus: 200e9,
        area: 0.001,
    };

    let k_local =
        elem.local_stiffness_matrix();

    let elements = vec![
        (
            k_local.clone(),
            0,
            1,
        ),
        (
            k_local.clone(),
            1,
            2,
        ),
    ];

    let global_k = assemble_global_stiffness_matrix(3, &elements);

    // Fixed at node 0, force at node 2
    let forces = vec![0.0, 0.0, 1000.0];

    let fixed_dofs = vec![(0, 0.0)];

    let result =
        solve_static_structural(
            global_k,
            forces,
            &fixed_dofs,
        );

    assert!(result.is_ok());

    let u = result.unwrap();

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
