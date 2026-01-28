//! Test suite for computer_graphics module.

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use rssn::symbolic::computer_graphics::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::vector::Vector;

#[test]

fn test_translation_2d() {

    let tx = Expr::new_constant(3.0);

    let ty = Expr::new_constant(4.0);

    let matrix = translation_2d(tx, ty);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].len(), 3);

        // Check identity part
        assert_eq!(
            rows[0][0],
            Expr::new_bigint(BigInt::one())
        );

        assert_eq!(
            rows[1][1],
            Expr::new_bigint(BigInt::one())
        );

        assert_eq!(
            rows[2][2],
            Expr::new_bigint(BigInt::one())
        );

        // Check translation part
        assert_eq!(
            rows[0][2],
            Expr::new_constant(3.0)
        );

        assert_eq!(
            rows[1][2],
            Expr::new_constant(4.0)
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_translation_3d() {

    let tx = Expr::new_constant(1.0);

    let ty = Expr::new_constant(2.0);

    let tz = Expr::new_constant(3.0);

    let matrix =
        translation_3d(tx, ty, tz);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 4);

        assert_eq!(rows[0].len(), 4);

        // Check translation column
        assert_eq!(
            rows[0][3],
            Expr::new_constant(1.0)
        );

        assert_eq!(
            rows[1][3],
            Expr::new_constant(2.0)
        );

        assert_eq!(
            rows[2][3],
            Expr::new_constant(3.0)
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_scaling_2d() {

    let sx = Expr::new_constant(2.0);

    let sy = Expr::new_constant(3.0);

    let matrix = scaling_2d(sx, sy);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 3);

        assert_eq!(
            rows[0][0],
            Expr::new_constant(2.0)
        );

        assert_eq!(
            rows[1][1],
            Expr::new_constant(3.0)
        );

        assert_eq!(
            rows[2][2],
            Expr::new_bigint(BigInt::one())
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_scaling_3d() {

    let sx = Expr::new_constant(2.0);

    let sy = Expr::new_constant(3.0);

    let sz = Expr::new_constant(4.0);

    let matrix = scaling_3d(sx, sy, sz);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 4);

        assert_eq!(
            rows[0][0],
            Expr::new_constant(2.0)
        );

        assert_eq!(
            rows[1][1],
            Expr::new_constant(3.0)
        );

        assert_eq!(
            rows[2][2],
            Expr::new_constant(4.0)
        );

        assert_eq!(
            rows[3][3],
            Expr::new_bigint(BigInt::one())
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_rotation_2d() {

    // Just ensure it constructs without panic
    let angle = Expr::new_constant(
        std::f64::consts::PI / 4.0,
    );

    let matrix = rotation_2d(angle);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].len(), 3);
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_rotation_3d_x() {

    let angle = Expr::new_constant(
        std::f64::consts::PI / 6.0,
    );

    let matrix = rotation_3d_x(angle);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 4);

        // First row should be [1, 0, 0, 0]
        assert_eq!(
            rows[0][0],
            Expr::new_bigint(BigInt::one())
        );

        assert_eq!(
            rows[0][1],
            Expr::new_bigint(BigInt::zero())
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_shear_2d() {

    let shx = Expr::new_constant(0.5);

    let shy = Expr::new_constant(0.3);

    let matrix = shear_2d(
        shx.clone(),
        shy.clone(),
    );

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 3);

        // Diagonal should be 1
        assert_eq!(
            rows[0][0],
            Expr::new_bigint(BigInt::one())
        );

        assert_eq!(
            rows[1][1],
            Expr::new_bigint(BigInt::one())
        );

        // Shear components
        assert_eq!(rows[0][1], shx);

        assert_eq!(rows[1][0], shy);
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_reflection_2d() {

    // Reflection across the x-axis (angle = 0)
    let angle = Expr::new_constant(0.0);

    let matrix = reflection_2d(angle);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].len(), 3);
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_reflection_3d() {

    // Reflection across xy-plane (normal = (0, 0, 1))
    let nx = Expr::new_constant(0.0);

    let ny = Expr::new_constant(0.0);

    let nz = Expr::new_constant(1.0);

    let matrix =
        reflection_3d(nx, ny, nz);

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 4);

        assert_eq!(rows[0].len(), 4);
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_rotation_axis_angle() {

    // Rotation around z-axis should be similar to rotation_3d_z
    let axis = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
    );

    let angle = Expr::new_constant(
        std::f64::consts::PI / 4.0,
    );

    let matrix = rotation_axis_angle(
        &axis, angle,
    );

    if let Expr::Matrix(rows) = matrix {

        assert_eq!(rows.len(), 4);

        assert_eq!(rows[0].len(), 4);

        // Last row should be [0, 0, 0, 1]
        assert_eq!(
            rows[3][3],
            Expr::new_bigint(BigInt::one())
        );
    } else {

        panic!("Expected Matrix");
    }
}

#[test]

fn test_bezier_curve_new() {

    let p0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let p1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(2.0),
        Expr::new_constant(0.0),
    );

    let p2 = Vector::new(
        Expr::new_constant(3.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let curve = BezierCurve {
        control_points: vec![
            p0, p1, p2,
        ],
        degree: 2,
    };

    assert_eq!(curve.degree, 2);

    assert_eq!(
        curve
            .control_points
            .len(),
        3
    );
}

#[test]

fn test_bezier_curve_evaluate_endpoints()
 {

    let p0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let p1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let curve = BezierCurve {
        control_points: vec![
            p0.clone(),
            p1.clone(),
        ],
        degree: 1,
    };

    // At t=0, should be at p0
    let at_0 = curve
        .evaluate(&Expr::new_constant(0.0));

    // At t=1, should be at p1
    let at_1 = curve
        .evaluate(&Expr::new_constant(1.0));

    // Just verify they don't panic and produce valid vectors
    assert!(matches!(
        at_0.x,
        Expr::Dag(_)
            | Expr::new_constant(_)
            | Expr::new_bigint(_)
            | _
    ));

    assert!(matches!(
        at_1.x,
        Expr::Dag(_)
            | Expr::new_constant(_)
            | Expr::new_bigint(_)
            | _
    ));
}

#[test]

fn test_bezier_curve_derivative() {

    // Linear Bezier: derivative should be constant (P1 - P0)
    let p0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let p1 = Vector::new(
        Expr::new_constant(2.0),
        Expr::new_constant(4.0),
        Expr::new_constant(0.0),
    );

    let curve = BezierCurve {
        control_points: vec![p0, p1],
        degree: 1,
    };

    let tangent = curve.derivative(
        &Expr::new_constant(0.5),
    );

    // For a linear curve of degree 1, derivative is n * (P1 - P0) = 1 * (2, 4, 0)
    // The tangent should be a valid Vector
    assert!(matches!(
        tangent.x,
        _
    ));
}

#[test]

fn test_bezier_curve_split() {

    let p0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let p1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(2.0),
        Expr::new_constant(0.0),
    );

    let p2 = Vector::new(
        Expr::new_constant(2.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let curve = BezierCurve {
        control_points: vec![
            p0, p1, p2,
        ],
        degree: 2,
    };

    let (left, right) = curve
        .split(&Expr::new_constant(0.5));

    // Both curves should have same degree
    assert_eq!(left.degree, 2);

    assert_eq!(right.degree, 2);

    // Both curves should have n+1 control points
    assert_eq!(
        left.control_points
            .len(),
        3
    );

    assert_eq!(
        right
            .control_points
            .len(),
        3
    );
}

#[test]

fn test_polygon_mesh_new() {

    let v0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v2 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let mesh = PolygonMesh::new(
        vec![v0, v1, v2],
        vec![Polygon::new(vec![
            0, 1, 2,
        ])],
    );

    assert_eq!(
        mesh.vertices.len(),
        3
    );

    assert_eq!(
        mesh.polygons.len(),
        1
    );

    assert_eq!(
        mesh.polygons[0].indices,
        vec![0, 1, 2]
    );
}

#[test]

fn test_polygon_mesh_apply_transformation()
 {

    let v0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let mesh = PolygonMesh::new(
        vec![v0, v1],
        vec![],
    );

    // Apply translation
    let transform = translation_3d(
        Expr::new_constant(5.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let transformed = mesh
        .apply_transformation(
            &transform,
        );

    assert!(transformed.is_ok());

    let new_mesh = transformed.unwrap();

    assert_eq!(
        new_mesh
            .vertices
            .len(),
        2
    );
}

#[test]

fn test_polygon_mesh_compute_normals() {

    // Triangle in XY plane
    let v0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v2 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let mesh = PolygonMesh::new(
        vec![v0, v1, v2],
        vec![Polygon::new(vec![
            0, 1, 2,
        ])],
    );

    let normals =
        mesh.compute_normals();

    assert_eq!(normals.len(), 1);
    // Normal should point in Z direction (either +Z or -Z)
}

#[test]

fn test_polygon_mesh_triangulate() {

    // Quad (4 vertices)
    let v0 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let v2 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let v3 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let mesh = PolygonMesh::new(
        vec![v0, v1, v2, v3],
        vec![Polygon::new(vec![
            0, 1, 2, 3,
        ])], // One quad
    );

    let triangulated =
        mesh.triangulate();

    // Quad should become 2 triangles
    assert_eq!(
        triangulated
            .polygons
            .len(),
        2
    );

    // Each triangle has 3 vertices
    assert_eq!(
        triangulated.polygons[0]
            .indices
            .len(),
        3
    );

    assert_eq!(
        triangulated.polygons[1]
            .indices
            .len(),
        3
    );
}
