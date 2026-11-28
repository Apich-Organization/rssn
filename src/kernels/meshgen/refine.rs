//! Quality-driven refinement of a planar triangulation.
//!
//! [`refine_quality`] follows Ruppert's scheme: while a triangle has a small
//! angle or a large area, its circumcentre is inserted as a new vertex. A
//! circumcentre that falls outside the mesh triggers a split of the encroached
//! convex-hull edge instead, which is what keeps the boundary triangles from
//! degenerating into slivers.
//!
//! # Termination
//!
//! The number of inserted points is capped at `20 * points.len() + 1000`; when
//! the cap is reached the mesh built so far is returned. Each pass inserts at
//! least one point or stops, so the routine always terminates.

use super::delaunay::delaunay_2d;
use super::util::TriMesh2;
use super::util::centroid;
use super::util::circumcircle;
use super::util::distance2;
use super::util::is_hull_edge;
use super::util::midpoint;
use super::util::point_in_triangle;
use super::util::triangle_area;
use super::util::triangle_min_angle;

/// Refines `triangles` until every triangle satisfies an angle and area bound.
///
/// `min_angle_deg` is the smallest angle a triangle may have, and `max_area`
/// the largest area it may have. New vertices are inserted at circumcentres
/// (or at convex-hull edge midpoints when the circumcentre leaves the mesh),
/// and the triangulation is rebuilt as a Delaunay mesh after every pass.
///
/// # Errors
///
/// Returns `Err` when fewer than three points are given, when a triangle index
/// is out of range, when `min_angle_deg` is not in `(0, 60)`, or when
/// `max_area` is not positive and finite.
pub fn refine_quality(
    points: &[[f64; 2]],
    triangles: &[[usize; 3]],
    min_angle_deg: f64,
    max_area: f64,
) -> Result<TriMesh2, String> {
    if points.len() < 3 {
        return Err("at least three points are required".to_string());
    }
    if min_angle_deg.is_nan() || min_angle_deg <= 0.0 || min_angle_deg >= 60.0 {
        return Err("minimum angle must lie strictly between 0 and 60 degrees".to_string());
    }
    if max_area.is_nan() || max_area <= 0.0 || !max_area.is_finite() {
        return Err("maximum area must be positive and finite".to_string());
    }
    for t in triangles {
        if t[0] >= points.len() || t[1] >= points.len() || t[2] >= points.len() {
            return Err("triangle index out of range".to_string());
        }
    }

    let mut pts = points.to_vec();
    let mut tris = triangles.to_vec();

    if tris.is_empty() {
        return Ok(TriMesh2 {
            points: pts,
            triangles: tris,
        });
    }

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in &pts {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let scale = (max_x - min_x).max(max_y - min_y).max(1.0);
    let tol = 1e-9 * scale;

    let min_rad = min_angle_deg.to_radians();
    let cap = 20 * points.len() + 1000;
    let mut passes = 0_usize;

    loop {
        if pts.len() >= cap || passes > 2 * cap {
            break;
        }
        passes += 1;

        let mut candidates: Vec<[f64; 2]> = Vec::new();
        for t in &tris {
            let a = pts[t[0]];
            let b = pts[t[1]];
            let c = pts[t[2]];
            let angle = triangle_min_angle(a, b, c);
            let area = triangle_area(a, b, c);
            if angle >= min_rad && area <= max_area {
                continue;
            }
            if let Some(p) = candidate(&pts, &tris, *t, min_rad, tol) {
                candidates.push(p);
            }
        }

        if candidates.is_empty() {
            break;
        }

        let mut inserted = false;
        for p in candidates {
            if pts.len() >= cap {
                break;
            }
            if far_enough(&pts, p, tol) {
                pts.push(p);
                inserted = true;
            }
        }

        if !inserted {
            break;
        }

        tris = delaunay_2d(&pts)?;
    }

    Ok(TriMesh2 {
        points: pts,
        triangles: tris,
    })
}

/// Chooses an insertion point for the triangle `t`.
///
/// The circumcentre is preferred; when it lies outside the mesh and `t` has a
/// convex-hull edge, the midpoint of that edge is used instead; the centroid is
/// the last resort.
fn candidate(
    pts: &[[f64; 2]],
    tris: &[[usize; 3]],
    t: [usize; 3],
    min_rad: f64,
    tol: f64,
) -> Option<[f64; 2]> {
    let a = pts[t[0]];
    let b = pts[t[1]];
    let c = pts[t[2]];

    if let Some((centre, _)) = circumcircle(a, b, c) {
        if inside_mesh(pts, tris, centre) {
            if far_enough(pts, centre, tol) {
                return Some(centre);
            }
        } else {
            for &(p, q) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                if is_hull_edge(pts[p], pts[q], pts, tol) {
                    let m = midpoint(pts[p], pts[q]);
                    if far_enough(pts, m, tol) {
                        return Some(m);
                    }
                }
            }
        }
    }

    // A triangle that is already well shaped needs no further work even when
    // the circumcentre is unavailable.
    if triangle_min_angle(a, b, c) >= min_rad {
        return None;
    }

    let g = centroid(a, b, c);
    if far_enough(pts, g, tol) {
        Some(g)
    } else {
        None
    }
}

/// Whether `p` lies strictly inside at least one triangle of `tris`.
fn inside_mesh(
    pts: &[[f64; 2]],
    tris: &[[usize; 3]],
    p: [f64; 2],
) -> bool {
    tris.iter()
        .any(|t| point_in_triangle(p, pts[t[0]], pts[t[1]], pts[t[2]], f64::EPSILON))
}

/// Whether `p` is at least `tol` away from every point already in `pts`.
fn far_enough(
    pts: &[[f64; 2]],
    p: [f64; 2],
    tol: f64,
) -> bool {
    let limit = tol * tol;
    pts.iter().all(|&q| distance2(p, q) > limit)
}

#[cfg(test)]
mod tests {
    use super::refine_quality;
    use crate::kernels::meshgen::util::min_angle_deg;
    use crate::kernels::meshgen::util::triangle_area;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    #[test]
    fn coarse_square_is_refined_below_the_area_bound() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = [[0, 1, 2], [0, 2, 3]];

        let mesh = expect_ok(refine_quality(&points, &triangles, 20.0, 0.01));

        assert!(mesh.triangles.len() > triangles.len());
        assert!(mesh.points.len() > points.len());

        let mut worst_area = 0.0_f64;
        for t in &mesh.triangles {
            let area = triangle_area(mesh.points[t[0]], mesh.points[t[1]], mesh.points[t[2]]);
            worst_area = worst_area.max(area);
        }
        assert!(
            worst_area <= 0.01 + 1e-9,
            "largest triangle area {worst_area}"
        );

        let angle = min_angle_deg(&mesh.points, &mesh.triangles);
        assert!(angle >= 15.0 - 1e-9, "smallest angle {angle} degrees");
    }

    #[test]
    fn already_good_mesh_is_returned_unchanged() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = [[0, 1, 2], [0, 2, 3]];
        let mesh = expect_ok(refine_quality(&points, &triangles, 40.0, 10.0));
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.triangles.len(), 2);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = [[0, 1, 2]];
        assert!(refine_quality(&points, &triangles, 0.0, 1.0).is_err());
        assert!(refine_quality(&points, &triangles, 70.0, 1.0).is_err());
        assert!(refine_quality(&points, &triangles, 20.0, -1.0).is_err());
        assert!(refine_quality(&points, &triangles, 20.0, f64::INFINITY).is_err());
        assert!(refine_quality(&points, &[[0, 1, 7]], 20.0, 1.0).is_err());
        assert!(refine_quality(&[[0.0, 0.0], [1.0, 0.0]], &[], 20.0, 1.0).is_err());
    }
}
