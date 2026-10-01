//! Computer-graphics kernels (ported from `numerical_computer_graphics_test.rs`).
//!
//! Tests for vectors, matrices, transformations, quaternions, ray tracing, and curves.

use std::f64::consts::PI;

use rssn::kernels::computer_graphics::*;

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

    assert!((v.magnitude() - 5.0).abs() < 1e-10);
}

#[test]

fn test_vector2d_normalize() {
    let v = Vector2D::new(3.0, 4.0);

    let n = v.normalize();

    assert!((n.magnitude() - 1.0).abs() < 1e-10);
}

#[test]

fn test_vector2d_rotate() {
    let v = Vector2D::new(1.0, 0.0);

    let rotated = v.rotate(PI / 2.0);

    assert!(rotated.x.abs() < 1e-10);

    assert!((rotated.y - 1.0).abs() < 1e-10);
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
    let v = Vector3D::new(1.0, 2.0, 3.0);

    assert_eq!(v.x, 1.0);

    assert_eq!(v.y, 2.0);

    assert_eq!(v.z, 3.0);
}

#[test]

fn test_vector3d_magnitude() {
    let v = Vector3D::new(1.0, 2.0, 2.0);

    assert!((v.magnitude() - 3.0).abs() < 1e-10);
}

#[test]

fn test_vector3d_normalize() {
    let v = Vector3D::new(1.0, 2.0, 2.0);

    let n = v.normalize();

    assert!((n.magnitude() - 1.0).abs() < 1e-10);
}

#[test]

fn test_vector3d_add() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);

    let v2 = Vector3D::new(4.0, 5.0, 6.0);

    let sum = v1 + v2;

    assert_eq!(sum.x, 5.0);

    assert_eq!(sum.y, 7.0);

    assert_eq!(sum.z, 9.0);
}

#[test]

fn test_vector3d_sub() {
    let v1 = Vector3D::new(4.0, 5.0, 6.0);

    let v2 = Vector3D::new(1.0, 2.0, 3.0);

    let diff = v1 - v2;

    assert_eq!(diff.x, 3.0);

    assert_eq!(diff.y, 3.0);

    assert_eq!(diff.z, 3.0);
}

#[test]

fn test_vector3d_scalar_mul() {
    let v = Vector3D::new(1.0, 2.0, 3.0);

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
    let v1 = Vector3D::new(1.0, 0.0, 0.0);

    let v2 = Vector3D::new(0.0, 1.0, 0.0);

    assert_eq!(dot_product(&v1, &v2), 0.0);
}

#[test]

fn test_dot_product_parallel() {
    let v1 = Vector3D::new(1.0, 0.0, 0.0);

    let v2 = Vector3D::new(2.0, 0.0, 0.0);

    assert_eq!(dot_product(&v1, &v2), 2.0);
}

#[test]

fn test_cross_product() {
    let v1 = Vector3D::new(1.0, 0.0, 0.0);

    let v2 = Vector3D::new(0.0, 1.0, 0.0);

    let cross = cross_product(&v1, &v2);

    assert_eq!(cross.x, 0.0);

    assert_eq!(cross.y, 0.0);

    assert_eq!(cross.z, 1.0);
}

#[test]

fn test_cross_product_anticommutative() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);

    let v2 = Vector3D::new(4.0, 5.0, 6.0);

    let c1 = cross_product(&v1, &v2);

    let c2 = cross_product(&v2, &v1);

    assert!((c1.x + c2.x).abs() < 1e-10);

    assert!((c1.y + c2.y).abs() < 1e-10);

    assert!((c1.z + c2.z).abs() < 1e-10);
}

// ============================================================================
// Reflection and Refraction Tests
// ============================================================================

#[test]

fn test_reflect() {
    let incident = Vector3D::new(1.0, -1.0, 0.0);

    let normal = Vector3D::new(0.0, 1.0, 0.0);

    let reflected = reflect(&incident, &normal);

    assert!((reflected.x - 1.0).abs() < 1e-10);

    assert!((reflected.y - 1.0).abs() < 1e-10);
}

#[test]

fn test_refract_straight() {
    let incident = Vector3D::new(0.0, -1.0, 0.0);

    let normal = Vector3D::new(0.0, 1.0, 0.0);

    let refracted = refract(&incident, &normal, 1.0);

    assert!(refracted.is_some());

    let r = refracted.unwrap_or_else(|| panic!("expected Some"));

    assert!((r.x).abs() < 1e-10);

    assert!((r.y + 1.0).abs() < 1e-10);
}

// ============================================================================
// Interpolation Tests
// ============================================================================

#[test]

fn test_lerp() {
    let v1 = Vector3D::new(0.0, 0.0, 0.0);

    let v2 = Vector3D::new(2.0, 4.0, 6.0);

    let mid = lerp(&v1, &v2, 0.5);

    assert!((mid.x - 1.0).abs() < 1e-10);

    assert!((mid.y - 2.0).abs() < 1e-10);

    assert!((mid.z - 3.0).abs() < 1e-10);
}

#[test]

fn test_lerp_endpoints() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);

    let v2 = Vector3D::new(4.0, 5.0, 6.0);

    let start = lerp(&v1, &v2, 0.0);

    assert!((start.x - v1.x).abs() < 1e-10);

    let end = lerp(&v1, &v2, 1.0);

    assert!((end.x - v2.x).abs() < 1e-10);
}

#[test]

fn test_angle_between() {
    let v1 = Vector3D::new(1.0, 0.0, 0.0);

    let v2 = Vector3D::new(0.0, 1.0, 0.0);

    let angle = angle_between(&v1, &v2);

    assert!((angle - PI / 2.0).abs() < 1e-10);
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

    assert!((mid.r - 0.5).abs() < 1e-10);

    assert!((mid.g - 0.5).abs() < 1e-10);
}

#[test]

fn test_color_clamp() {
    let c = Color::new(1.5, -0.5, 0.5, 1.0);

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
    let m = translation_matrix(1.0, 2.0, 3.0);

    assert_eq!(*m.get(0, 3), 1.0);

    assert_eq!(*m.get(1, 3), 2.0);

    assert_eq!(*m.get(2, 3), 3.0);
}

#[test]

fn test_scaling_matrix() {
    let m = scaling_matrix(2.0, 3.0, 4.0);

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
                assert_eq!(*m.get(i, j), 1.0);
            } else {
                assert_eq!(*m.get(i, j), 0.0);
            }
        }
    }
}

#[test]

fn test_rotation_matrix_x() {
    let m = rotation_matrix_x(0.0);

    // Should be identity for 0 rotation
    assert!((*m.get(1, 1) - 1.0).abs() < 1e-10);

    assert!((*m.get(2, 2) - 1.0).abs() < 1e-10);
}

#[test]

fn test_rotation_matrix_y() {
    let m = rotation_matrix_y(0.0);

    assert!((*m.get(0, 0) - 1.0).abs() < 1e-10);

    assert!((*m.get(2, 2) - 1.0).abs() < 1e-10);
}

#[test]

fn test_rotation_matrix_z() {
    let m = rotation_matrix_z(0.0);

    assert!((*m.get(0, 0) - 1.0).abs() < 1e-10);

    assert!((*m.get(1, 1) - 1.0).abs() < 1e-10);
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
    let q = Quaternion::new(1.0, 0.0, 0.0, 0.0);

    assert!((q.magnitude() - 1.0).abs() < 1e-10);
}

#[test]

fn test_quaternion_normalize() {
    let q = Quaternion::new(2.0, 0.0, 0.0, 0.0);

    let n = q.normalize();

    assert!((n.magnitude() - 1.0).abs() < 1e-10);
}

#[test]

fn test_quaternion_conjugate() {
    let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);

    let c = q.conjugate();

    assert_eq!(c.w, 1.0);

    assert_eq!(c.x, -2.0);

    assert_eq!(c.y, -3.0);

    assert_eq!(c.z, -4.0);
}

#[test]

fn test_quaternion_multiply_identity() {
    let q = Quaternion::new(1.0, 2.0, 3.0, 4.0).normalize();

    let identity = Quaternion::identity();

    let result = q.multiply(&identity);

    assert!((result.w - q.w).abs() < 1e-10);

    assert!((result.x - q.x).abs() < 1e-10);
}

#[test]

fn test_quaternion_from_axis_angle() {
    let q = Quaternion::from_axis_angle(&Vector3D::new(0.0, 0.0, 1.0), PI / 2.0);

    // Should be approximately (cos(pi/4), 0, 0, sin(pi/4))
    assert!((q.w - (PI / 4.0).cos()).abs() < 1e-10);

    assert!((q.z - (PI / 4.0).sin()).abs() < 1e-10);
}

#[test]

fn test_quaternion_rotate_vector() {
    // 90 degree rotation around Z axis
    let q = Quaternion::from_axis_angle(&Vector3D::new(0.0, 0.0, 1.0), PI / 2.0);

    let v = Vector3D::new(1.0, 0.0, 0.0);

    let rotated = q.rotate_vector(&v);

    assert!(rotated.x.abs() < 1e-10);

    assert!((rotated.y - 1.0).abs() < 1e-10);
}

// ============================================================================
// Ray Tracing Tests
// ============================================================================

#[test]

fn test_ray_at() {
    let ray = Ray::new(Point3D::new(0.0, 0.0, 0.0), Vector3D::new(1.0, 0.0, 0.0));

    let p = ray.at(2.0);

    assert_eq!(p.x, 2.0);

    assert_eq!(p.y, 0.0);

    assert_eq!(p.z, 0.0);
}

#[test]

fn test_ray_sphere_intersection_hit() {
    let ray = Ray::new(Point3D::new(0.0, 0.0, -5.0), Vector3D::new(0.0, 0.0, 1.0));

    let sphere = Sphere::new(Point3D::new(0.0, 0.0, 0.0), 1.0);

    let hit = ray_sphere_intersection(&ray, &sphere);

    assert!(hit.is_some());

    let h = hit.unwrap_or_else(|| panic!("expected Some"));

    assert!((h.t - 4.0).abs() < 1e-10);
}

#[test]

fn test_ray_sphere_intersection_miss() {
    let ray = Ray::new(Point3D::new(0.0, 5.0, -5.0), Vector3D::new(0.0, 0.0, 1.0));

    let sphere = Sphere::new(Point3D::new(0.0, 0.0, 0.0), 1.0);

    let hit = ray_sphere_intersection(&ray, &sphere);

    assert!(hit.is_none());
}

#[test]

fn test_ray_plane_intersection_hit() {
    let ray = Ray::new(Point3D::new(0.0, 1.0, 0.0), Vector3D::new(0.0, -1.0, 0.0));

    let plane = Plane::new(Point3D::new(0.0, 0.0, 0.0), Vector3D::new(0.0, 1.0, 0.0));

    let hit = ray_plane_intersection(&ray, &plane);

    assert!(hit.is_some());

    let h = hit.unwrap_or_else(|| panic!("expected Some"));

    assert!((h.t - 1.0).abs() < 1e-10);
}

#[test]

fn test_ray_triangle_intersection_hit() {
    let ray = Ray::new(Point3D::new(0.25, 0.25, -1.0), Vector3D::new(0.0, 0.0, 1.0));

    let v0 = Point3D::new(0.0, 0.0, 0.0);

    let v1 = Point3D::new(1.0, 0.0, 0.0);

    let v2 = Point3D::new(0.0, 1.0, 0.0);

    let hit = ray_triangle_intersection(&ray, &v0, &v1, &v2);

    assert!(hit.is_some());

    let h = hit.unwrap_or_else(|| panic!("expected Some"));

    assert!((h.t - 1.0).abs() < 1e-10);
}

// ============================================================================
// Curve Tests
// ============================================================================

#[test]

fn test_bezier_quadratic_endpoints() {
    let p0 = Point3D::new(0.0, 0.0, 0.0);

    let p1 = Point3D::new(0.5, 1.0, 0.0);

    let p2 = Point3D::new(1.0, 0.0, 0.0);

    let start = bezier_quadratic(&p0, &p1, &p2, 0.0);

    assert!((start.x - p0.x).abs() < 1e-10);

    let end = bezier_quadratic(&p0, &p1, &p2, 1.0);

    assert!((end.x - p2.x).abs() < 1e-10);
}

#[test]

fn test_bezier_cubic_endpoints() {
    let p0 = Point3D::new(0.0, 0.0, 0.0);

    let p1 = Point3D::new(0.25, 1.0, 0.0);

    let p2 = Point3D::new(0.75, 1.0, 0.0);

    let p3 = Point3D::new(1.0, 0.0, 0.0);

    let start = bezier_cubic(&p0, &p1, &p2, &p3, 0.0);

    assert!((start.x - p0.x).abs() < 1e-10);

    let end = bezier_cubic(&p0, &p1, &p2, &p3, 1.0);

    assert!((end.x - p3.x).abs() < 1e-10);
}

#[test]

fn test_catmull_rom_through_points() {
    let p0 = Point3D::new(-1.0, 0.0, 0.0);

    let p1 = Point3D::new(0.0, 0.0, 0.0);

    let p2 = Point3D::new(1.0, 0.0, 0.0);

    let p3 = Point3D::new(2.0, 0.0, 0.0);

    // At t=0, should be at p1
    let start = catmull_rom(&p0, &p1, &p2, &p3, 0.0);

    assert!((start.x - p1.x).abs() < 1e-10);

    // At t=1, should be at p2
    let end = catmull_rom(&p0, &p1, &p2, &p3, 1.0);

    assert!((end.x - p2.x).abs() < 1e-10);
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

    assert!((deg - 180.0).abs() < 1e-10);
}

#[test]

fn test_transform_point_identity() {
    let m = identity_matrix();

    let p = Point3D::new(1.0, 2.0, 3.0);

    let result = transform_point(&m, &p);

    assert!((result.x - 1.0).abs() < 1e-10);

    assert!((result.y - 2.0).abs() < 1e-10);

    assert!((result.z - 3.0).abs() < 1e-10);
}

#[test]

fn test_transform_point_translation() {
    let m = translation_matrix(1.0, 2.0, 3.0);

    let p = Point3D::new(0.0, 0.0, 0.0);

    let result = transform_point(&m, &p);

    assert!((result.x - 1.0).abs() < 1e-10);

    assert!((result.y - 2.0).abs() < 1e-10);

    assert!((result.z - 3.0).abs() < 1e-10);
}

#[test]

fn test_barycentric_coordinates() {
    let v0 = Point3D::new(0.0, 0.0, 0.0);

    let v1 = Point3D::new(1.0, 0.0, 0.0);

    let v2 = Point3D::new(0.0, 1.0, 0.0);

    // Centroid
    let center = Point3D::new(1.0 / 3.0, 1.0 / 3.0, 0.0);

    let (u, v, w) = barycentric_coordinates(&center, &v0, &v1, &v2);

    assert!((u - 1.0 / 3.0).abs() < 1e-10);

    assert!((v - 1.0 / 3.0).abs() < 1e-10);

    assert!((w - 1.0 / 3.0).abs() < 1e-10);
}

// ============================================================================
// Property Tests
// ============================================================================

mod proptests {

    use proptest::prelude::*;

    use super::*;

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

// ============================================================================
// Added: geometric identities with independently known results
// ============================================================================

mod strengthened {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::computer_graphics::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    fn near(
        a: Vector3D,
        b: Vector3D,
        tol: f64,
    ) -> bool {
        (a - b).magnitude() < tol
    }

    fn near_p(
        a: Point3D,
        x: f64,
        y: f64,
        z: f64,
    ) -> bool {
        (a.x - x).abs() < 1e-9 && (a.y - y).abs() < 1e-9 && (a.z - z).abs() < 1e-9
    }

    #[test]
    fn vector_algebra_reference_values() {
        let a = Vector3D::new(1.0, 2.0, 3.0);
        let b = Vector3D::new(4.0, 5.0, 6.0);
        assert_eq!(dot_product(&a, &b), 32.0);
        assert!(near(
            cross_product(&a, &b),
            Vector3D::new(-3.0, 6.0, -3.0),
            1e-12
        ));
        assert!(near(-a, Vector3D::new(-1.0, -2.0, -3.0), 0.0 + 1e-15));
        assert!(near(a / 2.0, Vector3D::new(0.5, 1.0, 1.5), 1e-15));
        assert_eq!(a.magnitude_squared(), 14.0);
        assert!(near(
            project(&a, &Vector3D::new(2.0, 0.0, 0.0)),
            Vector3D::new(1.0, 0.0, 0.0),
            1e-12
        ));
        assert!(near(
            project(&a, &Vector3D::new(0.0, 0.0, 0.0)),
            Vector3D::new(0.0, 0.0, 0.0),
            0.0 + 1e-15
        ));
        assert_eq!(
            dot_product_2d(&Vector2D::new(1.0, 2.0), &Vector2D::new(3.0, 4.0)),
            11.0
        );
        assert_eq!(angle_between(&a, &Vector3D::new(0.0, 0.0, 0.0)), 0.0);
        assert!((angle_between(&a, &-a) - PI).abs() < 1e-7);
    }

    #[test]
    fn points_and_2d_helpers() {
        assert!((Point2D::new(0.0, 0.0).distance_to(&Point2D::new(3.0, 4.0)) - 5.0).abs() < 1e-12);
        assert!(
            (Point3D::new(1.0, 2.0, 2.0).distance_to(&Point3D::new(0.0, 0.0, 0.0)) - 3.0).abs()
                < 1e-12
        );
        let v = Vector2D::new(3.0, 4.0);
        assert_eq!((v + v).x, 6.0);
        assert_eq!((v - v).y, 0.0);
        assert_eq!((v * 2.0).y, 8.0);
        assert_eq!((v / 2.0).x, 1.5);
        assert_eq!((-v).x, -3.0);
        // Rotating by 2 pi is the identity.
        let r = v.rotate(2.0 * PI);
        assert!((r.x - 3.0).abs() < 1e-12 && (r.y - 4.0).abs() < 1e-12);
    }

    #[test]
    fn rotation_matrices_rotate_basis_vectors() {
        let x = Point3D::new(1.0, 0.0, 0.0);
        let y = Point3D::new(0.0, 1.0, 0.0);
        let z = Point3D::new(0.0, 0.0, 1.0);
        assert!(near_p(
            transform_point(&rotation_matrix_z(FRAC_PI_2), &x),
            0.0,
            1.0,
            0.0
        ));
        assert!(near_p(
            transform_point(&rotation_matrix_x(FRAC_PI_2), &y),
            0.0,
            0.0,
            1.0
        ));
        assert!(near_p(
            transform_point(&rotation_matrix_y(FRAC_PI_2), &z),
            1.0,
            0.0,
            0.0
        ));
        let axis = rotation_matrix_axis(&Vector3D::new(0.0, 0.0, 1.0), FRAC_PI_2);
        assert!(near_p(transform_point(&axis, &x), 0.0, 1.0, 0.0));
    }

    #[test]
    fn composite_transforms_apply_right_to_left() {
        // scale by 2, then translate by (1, 0, 0): (1,1,1) -> (2,2,2) -> (3,2,2)
        let m = translation_matrix(1.0, 0.0, 0.0) * scaling_matrix(2.0, 2.0, 2.0);
        assert!(near_p(
            transform_point(&m, &Point3D::new(1.0, 1.0, 1.0)),
            3.0,
            2.0,
            2.0
        ));
        // Vectors ignore translation.
        let v = transform_vector(&m, &Vector3D::new(1.0, 1.0, 1.0));
        assert!(near(v, Vector3D::new(2.0, 2.0, 2.0), 1e-12));
        assert!(near_p(
            transform_point(&uniform_scaling_matrix(3.0), &Point3D::new(1.0, 2.0, 3.0)),
            3.0,
            6.0,
            9.0
        ));
    }

    #[test]
    fn quaternion_rotation_agrees_with_rotation_matrix() {
        let axis = Vector3D::new(0.0, 0.0, 1.0);
        let q = Quaternion::from_axis_angle(&axis, 0.9);
        let m = q.to_matrix();
        let p = Point3D::new(1.0, 2.0, 3.0);
        let via_q = q.rotate_vector(&Vector3D::new(1.0, 2.0, 3.0));
        let via_m = transform_point(&m, &p);
        assert!(near_p(via_m, via_q.x, via_q.y, via_q.z));
        assert!(near_p(
            via_m,
            transform_point(&rotation_matrix_z(0.9), &p).x,
            transform_point(&rotation_matrix_z(0.9), &p).y,
            3.0
        ));
    }

    #[test]
    fn quaternion_inverse_and_product() {
        let q = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        let id = q.multiply(&q.inverse());
        assert!(
            (id.w - 1.0).abs() < 1e-12
                && id.x.abs() < 1e-12
                && id.y.abs() < 1e-12
                && id.z.abs() < 1e-12
        );
        // i * j = k
        let i = Quaternion::new(0.0, 1.0, 0.0, 0.0);
        let j = Quaternion::new(0.0, 0.0, 1.0, 0.0);
        let k = i.multiply(&j);
        assert_eq!((k.w, k.x, k.y, k.z), (0.0, 0.0, 0.0, 1.0));
        // i * i = -1
        assert_eq!(i.multiply(&i).w, -1.0);
    }

    #[test]
    fn quaternion_slerp_halfway_is_half_the_rotation() {
        let axis = Vector3D::new(0.0, 0.0, 1.0);
        let a = Quaternion::identity();
        let b = Quaternion::from_axis_angle(&axis, FRAC_PI_2);
        let mid = a.slerp(&b, 0.5);
        let want = Quaternion::from_axis_angle(&axis, FRAC_PI_4);
        assert!((mid.w - want.w).abs() < 1e-9, "w = {} vs {}", mid.w, want.w);
        assert!((mid.z - want.z).abs() < 1e-9, "z = {} vs {}", mid.z, want.z);
        assert!(
            (mid.magnitude() - 1.0).abs() < 1e-9,
            "|slerp| = {}",
            mid.magnitude()
        );
    }

    #[test]
    fn quaternion_slerp_quarter_way() {
        let axis = Vector3D::new(0.0, 0.0, 1.0);
        let a = Quaternion::identity();
        let b = Quaternion::from_axis_angle(&axis, FRAC_PI_2);
        let q = a.slerp(&b, 0.25);
        let want = Quaternion::from_axis_angle(&axis, FRAC_PI_2 / 4.0);
        assert!(
            (q.w - want.w).abs() < 1e-9,
            "t=0.25: w = {} vs {}",
            q.w,
            want.w
        );
        assert!(
            (q.z - want.z).abs() < 1e-9,
            "t=0.25: z = {} vs {}",
            q.z,
            want.z
        );
    }

    #[test]
    fn vector_slerp_halfway_between_orthogonal_unit_vectors() {
        let m = slerp(
            &Vector3D::new(1.0, 0.0, 0.0),
            &Vector3D::new(0.0, 1.0, 0.0),
            0.5,
        );
        let s = std::f64::consts::FRAC_1_SQRT_2;
        assert!(near(m, Vector3D::new(s, s, 0.0), 1e-12));
    }

    #[test]
    fn snell_refraction_angle() {
        // 45 degrees into a medium with eta = 1/1.5 (air -> glass).
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let incident = Vector3D::new(s, -s, 0.0);
        let normal = Vector3D::new(0.0, 1.0, 0.0);
        let eta = 1.0 / 1.5;
        let r = refract(&incident, &normal, eta)
            .unwrap_or_else(|| panic!("no total internal reflection here"));
        let sin_out = r.x / r.magnitude();
        assert!((sin_out - eta * s).abs() < 1e-9, "sin(theta_t) = {sin_out}");
        assert!(r.y < 0.0);
        // Total internal reflection from glass to air at 60 degrees.
        let steep = Vector3D::new(0.866_025_403_784_438_6, -0.5, 0.0);
        assert!(refract(&steep, &normal, 1.5).is_none());
    }

    #[test]
    fn ray_sphere_hit_reports_point_and_normal() {
        let ray = Ray::new(Point3D::new(0.0, 0.0, -5.0), Vector3D::new(0.0, 0.0, 1.0));
        let sphere = Sphere::new(Point3D::new(0.0, 0.0, 0.0), 1.0);
        let h = ray_sphere_intersection(&ray, &sphere).unwrap_or_else(|| panic!("expected a hit"));
        assert!(near_p(h.point, 0.0, 0.0, -1.0));
        assert!(near(h.normal, Vector3D::new(0.0, 0.0, -1.0), 1e-12));
        // Starting inside the sphere the far side is returned.
        let inside = Ray::new(Point3D::new(0.0, 0.0, 0.0), Vector3D::new(0.0, 0.0, 1.0));
        let h =
            ray_sphere_intersection(&inside, &sphere).unwrap_or_else(|| panic!("expected a hit"));
        assert!((h.t - 1.0).abs() < 1e-9);
        // Sphere behind the ray.
        let behind = Ray::new(Point3D::new(0.0, 0.0, 5.0), Vector3D::new(0.0, 0.0, 1.0));
        assert!(ray_sphere_intersection(&behind, &sphere).is_none());
    }

    #[test]
    fn ray_plane_and_triangle_misses() {
        let parallel = Ray::new(Point3D::new(0.0, 1.0, 0.0), Vector3D::new(1.0, 0.0, 0.0));
        let plane = Plane::new(Point3D::new(0.0, 0.0, 0.0), Vector3D::new(0.0, 1.0, 0.0));
        assert!(ray_plane_intersection(&parallel, &plane).is_none());
        let (v0, v1, v2) = (
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
        );
        let outside = Ray::new(Point3D::new(2.0, 2.0, -1.0), Vector3D::new(0.0, 0.0, 1.0));
        assert!(ray_triangle_intersection(&outside, &v0, &v1, &v2).is_none());
    }

    #[test]
    fn curves_at_midpoint() {
        let (p0, p1, p2, p3) = (
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(0.0, 2.0, 0.0),
            Point3D::new(2.0, 2.0, 0.0),
            Point3D::new(2.0, 0.0, 0.0),
        );
        // Cubic Bezier at t = 1/2: (p0 + 3 p1 + 3 p2 + p3) / 8 = (1, 1.5, 0)
        assert!(near_p(bezier_cubic(&p0, &p1, &p2, &p3, 0.5), 1.0, 1.5, 0.0));
        // Quadratic Bezier at t = 1/2: (p0 + 2 p1 + p2)/4 = (0.5, 1, 0)
        assert!(near_p(bezier_quadratic(&p0, &p1, &p3, 0.5), 0.5, 1.0, 0.0));
        // Catmull-Rom of collinear equispaced points is linear.
        let (a, b, c, d) = (
            Point3D::new(-1.0, 0.0, 0.0),
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(2.0, 0.0, 0.0),
        );
        assert!(near_p(catmull_rom(&a, &b, &c, &d, 0.25), 0.25, 0.0, 0.0));
    }

    #[test]
    fn projection_matrices_have_expected_structure() {
        let persp = perspective_matrix(FRAC_PI_2, 1.0, 0.1, 100.0);
        // A point on the optical axis stays on the axis.
        let p = transform_point(&persp, &Point3D::new(0.0, 0.0, -10.0));
        assert!(p.x.abs() < 1e-12 && p.y.abs() < 1e-12);
        // With a 90 degree field of view and aspect 1, the point (z, z, -z) maps to the corner x = y = 1.
        let corner = transform_point(&persp, &Point3D::new(5.0, 5.0, -5.0));
        assert!(
            (corner.x - 1.0).abs() < 1e-9 && (corner.y - 1.0).abs() < 1e-9,
            "{corner:?}"
        );
        let ortho = orthographic_matrix(-2.0, 2.0, -1.0, 1.0, 0.1, 10.0);
        let q = transform_point(&ortho, &Point3D::new(2.0, 1.0, -0.1));
        assert!((q.x - 1.0).abs() < 1e-9 && (q.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn look_at_maps_eye_to_origin_and_target_onto_negative_z() {
        let eye = Point3D::new(1.0, 2.0, 3.0);
        let target = Point3D::new(1.0, 2.0, -7.0);
        let view = look_at_matrix(
            &eye.to_vector(),
            &target.to_vector(),
            &Vector3D::new(0.0, 1.0, 0.0),
        );
        assert!(near_p(transform_point(&view, &eye), 0.0, 0.0, 0.0));
        assert!(near_p(transform_point(&view, &target), 0.0, 0.0, -10.0));
    }

    #[test]
    fn barycentric_of_vertices_and_color_helpers() {
        let (v0, v1, v2) = (
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
        );
        let (u, v, w) = barycentric_coordinates(&v1, &v0, &v1, &v2);
        assert!((u - 0.0).abs() < 1e-12 && (v - 1.0).abs() < 1e-12 && w.abs() < 1e-12);
        let c = Color::rgb(0.2, 0.4, 0.6);
        assert_eq!((c.r, c.g, c.b, c.a), (0.2, 0.4, 0.6, 1.0));
        let l = Color::RED.lerp(&Color::BLUE, 0.25);
        assert!((l.r - 0.75).abs() < 1e-12 && (l.b - 0.25).abs() < 1e-12);
    }

    proptest! {
        #![proptest_config(cfg())]

        #[test]
        fn prop_rotation_matrices_preserve_length(
            x in -10.0..10.0f64, y in -10.0..10.0f64, z in -10.0..10.0f64, ang in -6.3..6.3f64,
        ) {
            let p = Point3D::new(x, y, z);
            let len = p.to_vector().magnitude();
            for m in [rotation_matrix_x(ang), rotation_matrix_y(ang), rotation_matrix_z(ang)] {
                let r = transform_point(&m, &p);
                prop_assert!((r.to_vector().magnitude() - len).abs() < 1e-9 * len.max(1.0));
            }
        }

        #[test]
        fn prop_quaternion_rotation_preserves_length_and_matches_matrix(
            ax in -1.0..1.0f64, ay in -1.0..1.0f64, az in -1.0..1.0f64, ang in -3.0..3.0f64,
            x in -10.0..10.0f64, y in -10.0..10.0f64, z in -10.0..10.0f64,
        ) {
            let axis = Vector3D::new(ax, ay, az);
            prop_assume!(axis.magnitude() > 1e-3);
            let q = Quaternion::from_axis_angle(&axis.normalize(), ang);
            let v = Vector3D::new(x, y, z);
            let r = q.rotate_vector(&v);
            prop_assert!((r.magnitude() - v.magnitude()).abs() < 1e-9 * v.magnitude().max(1.0));
            let m = transform_point(&q.to_matrix(), &Point3D::new(x, y, z));
            prop_assert!((m.x - r.x).abs() < 1e-9 && (m.y - r.y).abs() < 1e-9 && (m.z - r.z).abs() < 1e-9);
        }

        #[test]
        fn prop_angle_between_and_dot_agree(
            a in proptest::collection::vec(-5.0..5.0f64, 3), b in proptest::collection::vec(-5.0..5.0f64, 3),
        ) {
            let (u, v) = (Vector3D::new(a[0], a[1], a[2]), Vector3D::new(b[0], b[1], b[2]));
            prop_assume!(u.magnitude() > 1e-3 && v.magnitude() > 1e-3);
            let theta = angle_between(&u, &v);
            prop_assert!(((u.magnitude() * v.magnitude() * theta.cos()) - dot_product(&u, &v)).abs() < 1e-9);
        }

        #[test]
        fn prop_ray_sphere_hit_lies_on_the_sphere(
            ox in -0.5..0.5f64, oy in -0.5..0.5f64, radius in 1.0..3.0f64,
        ) {
            let ray = Ray::new(Point3D::new(ox, oy, -10.0), Vector3D::new(0.0, 0.0, 1.0));
            let sphere = Sphere::new(Point3D::new(0.0, 0.0, 0.0), radius);
            let h = ray_sphere_intersection(&ray, &sphere).ok_or_else(|| TestCaseError::fail("expected a hit"))?;
            prop_assert!((h.point.to_vector().magnitude() - radius).abs() < 1e-9);
            prop_assert!((h.normal.magnitude() - 1.0).abs() < 1e-9);
        }
    }
}

mod ledger_fill {
    use super::*;

    #[test]
    fn shearing_matrix_places_the_six_shear_factors() {
        let m = shearing_matrix(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        let expect = [
            [1.0, 1.0, 2.0, 0.0],
            [3.0, 1.0, 4.0, 0.0],
            [5.0, 6.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for (i, row) in expect.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert_eq!(*m.get(i, j), v, "({i},{j})");
            }
        }
    }

    #[test]
    fn from_euler_gives_identity_and_axis_rotations() {
        let q = Quaternion::from_euler(0.0, 0.0, 0.0);
        assert_eq!((q.w, q.x, q.y, q.z), (1.0, 0.0, 0.0, 0.0));
        let yaw = Quaternion::from_euler(0.0, 0.0, PI);
        assert!(yaw.w.abs() < 1e-12 && (yaw.z - 1.0).abs() < 1e-12);
        let roll = Quaternion::from_euler(PI / 2.0, 0.0, 0.0);
        assert!(
            (roll.w - (PI / 4.0).cos()).abs() < 1e-12 && (roll.x - (PI / 4.0).sin()).abs() < 1e-12
        );
    }
}

#[test]
fn quaternion_slerp_has_constant_angular_velocity() {
    // Rotation angle of slerp(t) relative to the start must be t * total angle.
    let axis = Vector3D::new(1.0, 2.0, 2.0).normalize();
    let a = Quaternion::identity();
    let total = 1.9;
    let b = Quaternion::from_axis_angle(&axis, total);
    for &t in &[0.1, 0.3, 0.7, 0.9] {
        let q = a.slerp(&b, t);
        let want = Quaternion::from_axis_angle(&axis, total * t);
        assert!(
            (q.w - want.w).abs() < 1e-9,
            "t={t}: w {} vs {}",
            q.w,
            want.w
        );
        assert!(
            (q.x - want.x).abs() < 1e-9,
            "t={t}: x {} vs {}",
            q.x,
            want.x
        );
        assert!((q.magnitude() - 1.0).abs() < 1e-9);
    }
}
