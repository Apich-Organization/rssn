//! Geometric and range regressions for prepared ray-triangle intersections.

#![allow(clippy::unwrap_used)]

use rssn::kernels::computer_graphics::Point3D;
use rssn::kernels::computer_graphics::PreparedTriangle;
use rssn::kernels::computer_graphics::Ray;
use rssn::kernels::computer_graphics::Vector3D;
use rssn::kernels::computer_graphics::ray_triangle_intersection;

const UNBOUNDED: f64 = f64::INFINITY;

fn xy_triangle(
    x: f64,
    y: f64,
) -> [Point3D; 3] {
    [
        Point3D::new(0.0, 0.0, 0.0),
        Point3D::new(x, 0.0, 0.0),
        Point3D::new(0.0, y, 0.0),
    ]
}

fn relative(
    got: f64,
    expected: f64,
) {
    assert!(
        (got / expected - 1.0).abs() < 5e-14,
        "got={got:e}, expected={expected:e}"
    );
}

#[test]
fn exact_plane_hit_survives_independently_scaled_triangle_edges() {
    let subnormal = f64::from_bits(4);
    for (x, y) in [
        (1e-4, 1e-4),
        (1e-200, 1e-200),
        (1e200, 1e200),
        (1e200, 1e-200),
        (subnormal, subnormal),
    ] {
        let [v0, v1, v2] = xy_triangle(x, y);
        let point = Point3D::new(x / 4.0, y / 4.0, 0.0);
        let ray = Ray::new(
            Point3D::new(point.x, point.y, 1.0),
            Vector3D::new(0.0, 0.0, -1.0),
        );
        let hit = ray_triangle_intersection(&ray, &v0, &v1, &v2).unwrap();
        relative(hit.t, 1.0);
        assert_eq!(hit.point, point);
        assert_eq!(hit.normal, Vector3D::new(0.0, 0.0, 1.0));

        let flipped = PreparedTriangle::new(&v0, &v2, &v1).unwrap();
        let hit = flipped.intersect(&ray, 1.0, 1.0).unwrap();
        assert_eq!(hit.point, point);
        assert_eq!(hit.normal, Vector3D::new(0.0, 0.0, -1.0));
    }
}

#[test]
fn prepared_ranges_are_inclusive_and_independent_of_direction_scale() {
    let [v0, v1, v2] = xy_triangle(1.0, 1.0);
    let prepared = PreparedTriangle::new(&v0, &v1, &v2).unwrap();
    let at_origin = Ray::new(Point3D::new(0.25, 0.25, 0.0), Vector3D::new(0.0, 0.0, -1.0));
    assert!(prepared.intersect(&at_origin, 0.0, 0.0).is_some());
    assert!(ray_triangle_intersection(&at_origin, &v0, &v1, &v2).is_none());

    for scale in [2.0_f64.powi(-700), 2.0_f64.powi(700)] {
        let ray = Ray::new(
            Point3D::new(0.25, 0.25, 1.0),
            Vector3D::new(0.0, 0.0, -scale),
        );
        let want = 1.0 / scale;
        let hit = prepared.intersect(&ray, 0.0, UNBOUNDED).unwrap();
        relative(hit.t, want);
        assert!((hit.point.z).abs() < 1e-15);
        assert!(prepared.intersect(&ray, want, want).is_some());
        assert!(prepared.intersect(&ray, 0.0, want.next_down()).is_none());
        assert!(
            prepared
                .intersect(&ray, want.next_up(), UNBOUNDED)
                .is_none()
        );
    }
}

#[test]
fn rejects_bad_ranges_nonfinite_rays_and_degenerate_triangles() {
    let [v0, v1, v2] = xy_triangle(1.0, 1.0);
    let prepared = PreparedTriangle::new(&v0, &v1, &v2).unwrap();
    let origin = Point3D::new(0.25, 0.25, 1.0);
    let ray = Ray::new(origin, Vector3D::new(0.0, 0.0, -1.0));
    for (min_t, max_t) in [
        (-1.0, 2.0),
        (2.0, 1.0),
        (f64::NAN, 2.0),
        (0.0, f64::NAN),
        (UNBOUNDED, UNBOUNDED),
    ] {
        assert!(prepared.intersect(&ray, min_t, max_t).is_none());
    }
    for direction in [
        Vector3D::new(1.0, 0.0, 0.0),
        Vector3D::new(0.0, 0.0, 0.0),
        Vector3D::new(f64::NAN, 0.0, -1.0),
    ] {
        assert!(
            prepared
                .intersect(&Ray::new(origin, direction), 0.0, UNBOUNDED)
                .is_none()
        );
    }
    assert!(PreparedTriangle::new(&v0, &v1, &v0).is_none());
    let p = Point3D::new(4.0, 16.0, 5.0);
    let q = Point3D::new(8.0, 32.0, 10.0);
    assert!(PreparedTriangle::new(&v0, &p, &q).is_none());
    assert!(PreparedTriangle::new(&v0, &v1, &Point3D::new(0.0, 1.0, f64::INFINITY)).is_none());
}

#[test]
fn tilted_plane_and_inclusive_edges_match_analytic_barycentric_points() {
    let normal_component = 3.0_f64.sqrt().recip();
    let [v0, v1, v2] = [
        Point3D::new(5.0, 7.0, 11.0),
        Point3D::new(7.0, 7.0, 13.0),
        Point3D::new(5.0, 10.0, 14.0),
    ];
    let triangle = PreparedTriangle::new(&v0, &v1, &v2).unwrap();
    let dir = Vector3D::new(0.25, -0.125, -1.0);
    for (u, v, expected) in [
        (0.25, 0.25, true),
        (0.0, 0.0, true),
        (1.0, 0.0, true),
        (0.5, 0.5, true),
        (0.75, 0.5, false),
        (-0.25, 0.25, false),
    ] {
        let point = Point3D::new(5.0 + 2.0 * u, 7.0 + 3.0 * v, 11.0 + 2.0 * u + 3.0 * v);
        let ray = Ray::new(
            Point3D::new(
                point.x - 2.0 * dir.x,
                point.y - 2.0 * dir.y,
                point.z - 2.0 * dir.z,
            ),
            dir,
        );
        let hit = triangle.intersect(&ray, 0.0, UNBOUNDED);
        assert_eq!(hit.is_some(), expected, "barycentric ({u}, {v})");
        if let Some(hit) = hit {
            relative(hit.t, 2.0);
            assert!((hit.point.x - point.x).abs() < 1e-13);
            assert!((hit.point.y - point.y).abs() < 1e-13);
            assert!((hit.point.z - point.z).abs() < 1e-13);
            relative(hit.normal.x, -normal_component);
            relative(hit.normal.y, -normal_component);
            relative(hit.normal.z, normal_component);
        }
    }
}

#[test]
fn exact_boundary_fallback_accepts_vertices_and_excludes_adjacent_outside_rays() {
    let v0 = Point3D::new(0.0, 0.017578125, 0.0078125);
    let v1 = Point3D::new(-0.0078125, 0.0078125, -0.001953125);
    let v2 = Point3D::new(-0.009765625, 0.025390625, -0.00390625);
    let ray = Ray::new(
        Point3D::new(-0.0048828125, 0.005859375, 0.0),
        Vector3D::new(-0.00146484375, 0.0009765625, -0.0009765625),
    );
    let hit = PreparedTriangle::new(&v0, &v1, &v2)
        .unwrap()
        .intersect(&ray, 2.0, 2.0)
        .unwrap();
    assert_eq!(hit.point, v1);
    assert_eq!(
        ray_triangle_intersection(&ray, &v0, &v1, &v2)
            .unwrap()
            .point,
        v1
    );

    let [v0, v1, v2] = xy_triangle(1.0, 1.0);
    let triangle = PreparedTriangle::new(&v0, &v1, &v2).unwrap();
    let direction = Vector3D::new(0.0, 0.0, -1.0);
    let on_edge = Ray::new(Point3D::new(0.5, 0.5, 1.0), direction);
    let above_edge = Ray::new(Point3D::new(0.5, 0.5_f64.next_up(), 1.0), direction);
    let below_edge = Ray::new(Point3D::new(0.5, 0.5_f64.next_down(), 1.0), direction);
    assert!(triangle.intersect(&on_edge, 0.0, UNBOUNDED).is_some());
    assert!(triangle.intersect(&above_edge, 0.0, UNBOUNDED).is_none());
    assert!(triangle.intersect(&below_edge, 0.0, UNBOUNDED).is_some());
}

#[test]
fn parallel_rays_keep_input_order_and_individual_range_results() {
    let [v0, v1, v2] = xy_triangle(1.0, 1.0);
    let prepared = PreparedTriangle::new(&v0, &v1, &v2).unwrap();
    let direction = Vector3D::new(0.0, 0.0, -1.0);
    let origins = [
        Point3D::new(0.25, 0.25, 1.0),
        Point3D::new(0.5, 0.5_f64.next_up(), 1.0),
        Point3D::new(0.5, 0.5, 2.0),
        Point3D::new(0.25, 0.25, 3.0),
    ];
    let expected: Vec<_> = origins
        .iter()
        .map(|&p| prepared.intersect(&Ray::new(p, direction), 0.5, 2.0))
        .collect();
    assert_eq!(
        prepared.intersect_parallel(&origins, direction, 0.5, 2.0),
        expected
    );
    assert_eq!(
        expected.iter().map(Option::is_some).collect::<Vec<_>>(),
        [true, false, true, false]
    );
    assert_eq!(
        prepared.intersect_parallel(&origins, direction, 3.0, 1.0),
        [None; 4]
    );
}
