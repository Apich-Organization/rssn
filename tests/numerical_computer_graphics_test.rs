//! Comprehensive tests for numerical computer graphics module.
//!
//! Tests for vectors, matrices, transformations, quaternions, ray tracing, and curves.

use rssn::numerical::computer_graphics::*;
use std::f64::consts::PI;

// ============================================================================
// Vector2D Tests
// ============================================================================

#[test]

fn test_vector2d_new() {

    let v = Vector2D::new(1.0, 2.0);

    assert_eq!(v.x, 1.0);

    assert_eq!(v.y, 2.0);
}

#[test]

fn test_vector2d_magnitude() {

    let v = Vector2D::new(3.0, 4.0);

    assert!(
        (v.magnitude() - 5.0).abs()
            < 1e-10
    );
}

#[test]

fn test_vector2d_normalize() {

    let v = Vector2D::new(3.0, 4.0);

    let n = v.normalize();

    assert!(
        (n.magnitude() - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_vector2d_rotate() {

    let v = Vector2D::new(1.0, 0.0);

    let rotated = v.rotate(PI / 2.0);

    assert!(rotated.x.abs() < 1e-10);

    assert!(
        (rotated.y - 1.0).abs() < 1e-10
    );
}

#[test]

fn test_vector2d_perpendicular() {

    let v = Vector2D::new(1.0, 0.0);

    let perp = v.perpendicular();

    assert_eq!(perp.x, 0.0);

    assert_eq!(perp.y, 1.0);
}

// ============================================================================
// Vector3D Tests
// ============================================================================

#[test]

fn test_vector3d_new() {

    let v =
        Vector3D::new(1.0, 2.0, 3.0);

    assert_eq!(v.x, 1.0);

    assert_eq!(v.y, 2.0);

    assert_eq!(v.z, 3.0);
}

#[test]

fn test_vector3d_magnitude() {

    let v =
        Vector3D::new(1.0, 2.0, 2.0);

    assert!(
        (v.magnitude() - 3.0).abs()
            < 1e-10
    );
}

#[test]

fn test_vector3d_normalize() {

    let v =
        Vector3D::new(1.0, 2.0, 2.0);

    let n = v.normalize();

    assert!(
        (n.magnitude() - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_vector3d_add() {

    let v1 =
        Vector3D::new(1.0, 2.0, 3.0);

    let v2 =
        Vector3D::new(4.0, 5.0, 6.0);

    let sum = v1 + v2;

    assert_eq!(sum.x, 5.0);

    assert_eq!(sum.y, 7.0);

    assert_eq!(sum.z, 9.0);
}

#[test]

fn test_vector3d_sub() {

    let v1 =
        Vector3D::new(4.0, 5.0, 6.0);

    let v2 =
        Vector3D::new(1.0, 2.0, 3.0);

    let diff = v1 - v2;

    assert_eq!(diff.x, 3.0);

    assert_eq!(diff.y, 3.0);

    assert_eq!(diff.z, 3.0);
}

#[test]

fn test_vector3d_scalar_mul() {

    let v =
        Vector3D::new(1.0, 2.0, 3.0);

    let result = v * 2.0;

    assert_eq!(result.x, 2.0);

    assert_eq!(result.y, 4.0);

    assert_eq!(result.z, 6.0);
}

// ============================================================================
// Dot and Cross Product Tests
// ============================================================================

#[test]

fn test_dot_product() {

    let v1 =
        Vector3D::new(1.0, 0.0, 0.0);

    let v2 =
        Vector3D::new(0.0, 1.0, 0.0);

    assert_eq!(
        dot_product(&v1, &v2),
        0.0
    );
}

#[test]

fn test_dot_product_parallel() {

    let v1 =
        Vector3D::new(1.0, 0.0, 0.0);

    let v2 =
        Vector3D::new(2.0, 0.0, 0.0);

    assert_eq!(
        dot_product(&v1, &v2),
        2.0
    );
}

#[test]

fn test_cross_product() {

    let v1 =
        Vector3D::new(1.0, 0.0, 0.0);

    let v2 =
        Vector3D::new(0.0, 1.0, 0.0);

    let cross = cross_product(&v1, &v2);

    assert_eq!(cross.x, 0.0);

    assert_eq!(cross.y, 0.0);

    assert_eq!(cross.z, 1.0);
}

#[test]

fn test_cross_product_anticommutative()
{

    let v1 =
        Vector3D::new(1.0, 2.0, 3.0);

    let v2 =
        Vector3D::new(4.0, 5.0, 6.0);

    let c1 = cross_product(&v1, &v2);

    let c2 = cross_product(&v2, &v1);

    assert!(
        (c1.x + c2.x).abs() < 1e-10
    );

    assert!(
        (c1.y + c2.y).abs() < 1e-10
    );

    assert!(
        (c1.z + c2.z).abs() < 1e-10
    );
}

// ============================================================================
// Reflection and Refraction Tests
// ============================================================================

#[test]

fn test_reflect() {

    let incident =
        Vector3D::new(1.0, -1.0, 0.0);

    let normal =
        Vector3D::new(0.0, 1.0, 0.0);

    let reflected =
        reflect(&incident, &normal);

    assert!(
        (reflected.x - 1.0).abs()
            < 1e-10
    );

    assert!(
        (reflected.y - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_refract_straight() {

    let incident =
        Vector3D::new(0.0, -1.0, 0.0);

    let normal =
        Vector3D::new(0.0, 1.0, 0.0);

    let refracted = refract(
        &incident, &normal, 1.0,
    );

    assert!(refracted.is_some());

    let r = refracted.unwrap();

    assert!((r.x).abs() < 1e-10);

    assert!((r.y + 1.0).abs() < 1e-10);
}

// ============================================================================
// Interpolation Tests
// ============================================================================

#[test]

fn test_lerp() {

    let v1 =
        Vector3D::new(0.0, 0.0, 0.0);

    let v2 =
        Vector3D::new(2.0, 4.0, 6.0);

    let mid = lerp(&v1, &v2, 0.5);

    assert!(
        (mid.x - 1.0).abs() < 1e-10
    );

    assert!(
        (mid.y - 2.0).abs() < 1e-10
    );

    assert!(
        (mid.z - 3.0).abs() < 1e-10
    );
}

#[test]

fn test_lerp_endpoints() {

    let v1 =
        Vector3D::new(1.0, 2.0, 3.0);

    let v2 =
        Vector3D::new(4.0, 5.0, 6.0);

    let start = lerp(&v1, &v2, 0.0);

    assert!(
        (start.x - v1.x).abs() < 1e-10
    );

    let end = lerp(&v1, &v2, 1.0);

    assert!(
        (end.x - v2.x).abs() < 1e-10
    );
}

#[test]

fn test_angle_between() {

    let v1 =
        Vector3D::new(1.0, 0.0, 0.0);

    let v2 =
        Vector3D::new(0.0, 1.0, 0.0);

    let angle = angle_between(&v1, &v2);

    assert!(
        (angle - PI / 2.0).abs()
            < 1e-10
    );
}

// ============================================================================
// Color Tests
// ============================================================================

#[test]

fn test_color_constants() {

    assert_eq!(Color::BLACK.r, 0.0);

    assert_eq!(Color::WHITE.r, 1.0);

    assert_eq!(Color::RED.r, 1.0);

    assert_eq!(Color::RED.g, 0.0);
}

#[test]

fn test_color_lerp() {

    let c1 = Color::BLACK;

    let c2 = Color::WHITE;

    let mid = c1.lerp(&c2, 0.5);

    assert!(
        (mid.r - 0.5).abs() < 1e-10
    );

    assert!(
        (mid.g - 0.5).abs() < 1e-10
    );
}

#[test]

fn test_color_clamp() {

    let c =
        Color::new(1.5, -0.5, 0.5, 1.0);

    let clamped = c.clamp();

    assert_eq!(clamped.r, 1.0);

    assert_eq!(clamped.g, 0.0);

    assert_eq!(clamped.b, 0.5);
}

// ============================================================================
// Matrix Tests
// ============================================================================

#[test]

fn test_translation_matrix() {

    let m = translation_matrix(
        1.0, 2.0, 3.0,
    );

    assert_eq!(*m.get(0, 3), 1.0);

    assert_eq!(*m.get(1, 3), 2.0);

    assert_eq!(*m.get(2, 3), 3.0);
}

#[test]

fn test_scaling_matrix() {

    let m =
        scaling_matrix(2.0, 3.0, 4.0);

    assert_eq!(*m.get(0, 0), 2.0);

    assert_eq!(*m.get(1, 1), 3.0);

    assert_eq!(*m.get(2, 2), 4.0);
}

#[test]

fn test_identity_matrix() {

    let m = identity_matrix();

    for i in 0..4 {

        for j in 0..4 {

            if i == j {

                assert_eq!(
                    *m.get(i, j),
                    1.0
                );
            } else {

                assert_eq!(
                    *m.get(i, j),
                    0.0
                );
            }
        }
    }
}

#[test]

fn test_rotation_matrix_x() {

    let m = rotation_matrix_x(0.0);

    // Should be identity for 0 rotation
    assert!(
        (*m.get(1, 1) - 1.0).abs()
            < 1e-10
    );

    assert!(
        (*m.get(2, 2) - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_rotation_matrix_y() {

    let m = rotation_matrix_y(0.0);

    assert!(
        (*m.get(0, 0) - 1.0).abs()
            < 1e-10
    );

    assert!(
        (*m.get(2, 2) - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_rotation_matrix_z() {

    let m = rotation_matrix_z(0.0);

    assert!(
        (*m.get(0, 0) - 1.0).abs()
            < 1e-10
    );

    assert!(
        (*m.get(1, 1) - 1.0).abs()
            < 1e-10
    );
}

// ============================================================================
// Quaternion Tests
// ============================================================================

#[test]

fn test_quaternion_identity() {

    let q = Quaternion::identity();

    assert_eq!(q.w, 1.0);

    assert_eq!(q.x, 0.0);

    assert_eq!(q.y, 0.0);

    assert_eq!(q.z, 0.0);
}

#[test]

fn test_quaternion_magnitude() {

    let q = Quaternion::new(
        1.0, 0.0, 0.0, 0.0,
    );

    assert!(
        (q.magnitude() - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_quaternion_normalize() {

    let q = Quaternion::new(
        2.0, 0.0, 0.0, 0.0,
    );

    let n = q.normalize();

    assert!(
        (n.magnitude() - 1.0).abs()
            < 1e-10
    );
}

#[test]

fn test_quaternion_conjugate() {

    let q = Quaternion::new(
        1.0, 2.0, 3.0, 4.0,
    );

    let c = q.conjugate();

    assert_eq!(c.w, 1.0);

    assert_eq!(c.x, -2.0);

    assert_eq!(c.y, -3.0);

    assert_eq!(c.z, -4.0);
}

#[test]

fn test_quaternion_multiply_identity() {

    let q = Quaternion::new(
        1.0, 2.0, 3.0, 4.0,
    )
    .normalize();

    let identity =
        Quaternion::identity();

    let result = q.multiply(&identity);

    assert!(
        (result.w - q.w).abs() < 1e-10
    );

    assert!(
        (result.x - q.x).abs() < 1e-10
    );
}

#[test]

fn test_quaternion_from_axis_angle() {

    let q = Quaternion::from_axis_angle(
        &Vector3D::new(0.0, 0.0, 1.0),
        PI / 2.0,
    );

    // Should be approximately (cos(pi/4), 0, 0, sin(pi/4))
    assert!(
        (q.w - (PI / 4.0).cos()).abs()
            < 1e-10
    );

    assert!(
        (q.z - (PI / 4.0).sin()).abs()
            < 1e-10
    );
}

#[test]

fn test_quaternion_rotate_vector() {

    // 90 degree rotation around Z axis
    let q = Quaternion::from_axis_angle(
        &Vector3D::new(0.0, 0.0, 1.0),
        PI / 2.0,
    );

    let v =
        Vector3D::new(1.0, 0.0, 0.0);

    let rotated = q.rotate_vector(&v);

    assert!(rotated.x.abs() < 1e-10);

    assert!(
        (rotated.y - 1.0).abs() < 1e-10
    );
}

// ============================================================================
// Ray Tracing Tests
// ============================================================================

#[test]

fn test_ray_at() {

    let ray = Ray::new(
        Point3D::new(0.0, 0.0, 0.0),
        Vector3D::new(1.0, 0.0, 0.0),
    );

    let p = ray.at(2.0);

    assert_eq!(p.x, 2.0);

    assert_eq!(p.y, 0.0);

    assert_eq!(p.z, 0.0);
}

#[test]

fn test_ray_sphere_intersection_hit() {

    let ray = Ray::new(
        Point3D::new(0.0, 0.0, -5.0),
        Vector3D::new(0.0, 0.0, 1.0),
    );

    let sphere = Sphere::new(
        Point3D::new(0.0, 0.0, 0.0),
        1.0,
    );

    let hit = ray_sphere_intersection(
        &ray, &sphere,
    );

    assert!(hit.is_some());

    let h = hit.unwrap();

    assert!((h.t - 4.0).abs() < 1e-10);
}

#[test]

fn test_ray_sphere_intersection_miss() {

    let ray = Ray::new(
        Point3D::new(0.0, 5.0, -5.0),
        Vector3D::new(0.0, 0.0, 1.0),
    );

    let sphere = Sphere::new(
        Point3D::new(0.0, 0.0, 0.0),
        1.0,
    );

    let hit = ray_sphere_intersection(
        &ray, &sphere,
    );

    assert!(hit.is_none());
}

#[test]

fn test_ray_plane_intersection_hit() {

    let ray = Ray::new(
        Point3D::new(0.0, 1.0, 0.0),
        Vector3D::new(0.0, -1.0, 0.0),
    );

    let plane = Plane::new(
        Point3D::new(0.0, 0.0, 0.0),
        Vector3D::new(0.0, 1.0, 0.0),
    );

    let hit = ray_plane_intersection(
        &ray, &plane,
    );

    assert!(hit.is_some());

    let h = hit.unwrap();

    assert!((h.t - 1.0).abs() < 1e-10);
}

#[test]

fn test_ray_triangle_intersection_hit()
{

    let ray = Ray::new(
        Point3D::new(0.25, 0.25, -1.0),
        Vector3D::new(0.0, 0.0, 1.0),
    );

    let v0 =
        Point3D::new(0.0, 0.0, 0.0);

    let v1 =
        Point3D::new(1.0, 0.0, 0.0);

    let v2 =
        Point3D::new(0.0, 1.0, 0.0);

    let hit = ray_triangle_intersection(
        &ray, &v0, &v1, &v2,
    );

    assert!(hit.is_some());

    let h = hit.unwrap();

    assert!((h.t - 1.0).abs() < 1e-10);
}

// ============================================================================
// Curve Tests
// ============================================================================

#[test]

fn test_bezier_quadratic_endpoints() {

    let p0 =
        Point3D::new(0.0, 0.0, 0.0);

    let p1 =
        Point3D::new(0.5, 1.0, 0.0);

    let p2 =
        Point3D::new(1.0, 0.0, 0.0);

    let start = bezier_quadratic(
        &p0, &p1, &p2, 0.0,
    );

    assert!(
        (start.x - p0.x).abs() < 1e-10
    );

    let end = bezier_quadratic(
        &p0, &p1, &p2, 1.0,
    );

    assert!(
        (end.x - p2.x).abs() < 1e-10
    );
}

#[test]

fn test_bezier_cubic_endpoints() {

    let p0 =
        Point3D::new(0.0, 0.0, 0.0);

    let p1 =
        Point3D::new(0.25, 1.0, 0.0);

    let p2 =
        Point3D::new(0.75, 1.0, 0.0);

    let p3 =
        Point3D::new(1.0, 0.0, 0.0);

    let start = bezier_cubic(
        &p0, &p1, &p2, &p3, 0.0,
    );

    assert!(
        (start.x - p0.x).abs() < 1e-10
    );

    let end = bezier_cubic(
        &p0, &p1, &p2, &p3, 1.0,
    );

    assert!(
        (end.x - p3.x).abs() < 1e-10
    );
}

#[test]

fn test_catmull_rom_through_points() {

    let p0 =
        Point3D::new(-1.0, 0.0, 0.0);

    let p1 =
        Point3D::new(0.0, 0.0, 0.0);

    let p2 =
        Point3D::new(1.0, 0.0, 0.0);

    let p3 =
        Point3D::new(2.0, 0.0, 0.0);

    // At t=0, should be at p1
    let start = catmull_rom(
        &p0, &p1, &p2, &p3, 0.0,
    );

    assert!(
        (start.x - p1.x).abs() < 1e-10
    );

    // At t=1, should be at p2
    let end = catmull_rom(
        &p0, &p1, &p2, &p3, 1.0,
    );

    assert!(
        (end.x - p2.x).abs() < 1e-10
    );
}

// ============================================================================
// Utility Function Tests
// ============================================================================

#[test]

fn test_degrees_to_radians() {

    let rad = degrees_to_radians(180.0);

    assert!((rad - PI).abs() < 1e-10);
}

#[test]

fn test_radians_to_degrees() {

    let deg = radians_to_degrees(PI);

    assert!(
        (deg - 180.0).abs() < 1e-10
    );
}

#[test]

fn test_transform_point_identity() {

    let m = identity_matrix();

    let p = Point3D::new(1.0, 2.0, 3.0);

    let result =
        transform_point(&m, &p);

    assert!(
        (result.x - 1.0).abs() < 1e-10
    );

    assert!(
        (result.y - 2.0).abs() < 1e-10
    );

    assert!(
        (result.z - 3.0).abs() < 1e-10
    );
}

#[test]

fn test_transform_point_translation() {

    let m = translation_matrix(
        1.0, 2.0, 3.0,
    );

    let p = Point3D::new(0.0, 0.0, 0.0);

    let result =
        transform_point(&m, &p);

    assert!(
        (result.x - 1.0).abs() < 1e-10
    );

    assert!(
        (result.y - 2.0).abs() < 1e-10
    );

    assert!(
        (result.z - 3.0).abs() < 1e-10
    );
}

#[test]

fn test_barycentric_coordinates() {

    let v0 =
        Point3D::new(0.0, 0.0, 0.0);

    let v1 =
        Point3D::new(1.0, 0.0, 0.0);

    let v2 =
        Point3D::new(0.0, 1.0, 0.0);

    // Centroid
    let center = Point3D::new(
        1.0 / 3.0,
        1.0 / 3.0,
        0.0,
    );

    let (u, v, w) =
        barycentric_coordinates(
            &center, &v0, &v1, &v2,
        );

    assert!(
        (u - 1.0 / 3.0).abs() < 1e-10
    );

    assert!(
        (v - 1.0 / 3.0).abs() < 1e-10
    );

    assert!(
        (w - 1.0 / 3.0).abs() < 1e-10
    );
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_normalize_magnitude_one(x in -100.0..100.0f64, y in -100.0..100.0f64, z in -100.0..100.0f64) {
            if x != 0.0 || y != 0.0 || z != 0.0 {
                let v = Vector3D::new(x, y, z);
                let n = v.normalize();
                prop_assert!((n.magnitude() - 1.0).abs() < 1e-10);
            }
        }

        #[test]
        fn prop_dot_product_parallel(scale in 0.1..10.0f64) {
            let v = Vector3D::new(1.0, 0.0, 0.0);
            let scaled = v * scale;
            prop_assert!((dot_product(&v, &scaled) - scale).abs() < 1e-10);
        }

        #[test]
        fn prop_cross_product_perpendicular(
            x1 in -10.0..10.0f64, y1 in -10.0..10.0f64, z1 in -10.0..10.0f64,
            x2 in -10.0..10.0f64, y2 in -10.0..10.0f64, z2 in -10.0..10.0f64,
        ) {
            let v1 = Vector3D::new(x1, y1, z1);
            let v2 = Vector3D::new(x2, y2, z2);
            let cross = cross_product(&v1, &v2);
            // Cross product is perpendicular to both inputs
            prop_assert!(dot_product(&cross, &v1).abs() < 1e-6);
            prop_assert!(dot_product(&cross, &v2).abs() < 1e-6);
        }

        #[test]
        fn prop_quaternion_normalize(w in -10.0..10.0f64, x in -10.0..10.0f64, y in -10.0..10.0f64, z in -10.0..10.0f64) {
            if w != 0.0 || x != 0.0 || y != 0.0 || z != 0.0 {
                let q = Quaternion::new(w, x, y, z);
                let n = q.normalize();
                prop_assert!((n.magnitude() - 1.0).abs() < 1e-10);
            }
        }

        #[test]
        fn prop_lerp_endpoints(t in 0.0..1.0f64) {
            let v1 = Vector3D::new(0.0, 0.0, 0.0);
            let v2 = Vector3D::new(1.0, 2.0, 3.0);
            let result = lerp(&v1, &v2, t);
            prop_assert!(result.x >= 0.0 && result.x <= 1.0);
            prop_assert!(result.y >= 0.0 && result.y <= 2.0);
            prop_assert!(result.z >= 0.0 && result.z <= 3.0);
        }

        #[test]
        fn prop_conversion_roundtrip(degrees in 0.0..360.0f64) {
            let result = radians_to_degrees(degrees_to_radians(degrees));
            prop_assert!((result - degrees).abs() < 1e-10);
        }
    }
}
