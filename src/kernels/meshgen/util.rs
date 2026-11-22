//! Shared geometry for the mesh-generation kernels.
//!
//! Predicates (orientation, in-circle, in-sphere), circumcircle geometry,
//! triangle quality measures, polygon helpers and the [`TriMesh2`] container
//! used by [`super::delaunay`], [`super::refine`], [`super::advancing_front`]
//! and [`super::tetra`].

/// An unstructured triangle mesh in the plane.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TriMesh2 {
    /// Vertex coordinates, one `[x, y]` per vertex.
    pub points: Vec<[f64; 2]>,
    /// Counter-clockwise triangles, given as triples of vertex indices.
    pub triangles: Vec<[usize; 3]>,
}

/// Twice the signed area of the triangle `a`, `b`, `c`.
///
/// The value is positive when `a`, `b`, `c` wind counter-clockwise and
/// negative when they wind clockwise.
pub(super) fn orient2d(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// The `z` component of `(b - a) x (c - a)` for three points in space.
pub(super) fn orient3d(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> f64 {
    let ux = b[0] - a[0];
    let uy = b[1] - a[1];
    let uz = b[2] - a[2];
    let vx = c[0] - a[0];
    let vy = c[1] - a[1];
    let vz = c[2] - a[2];
    let wx = d[0] - a[0];
    let wy = d[1] - a[1];
    let wz = d[2] - a[2];
    ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)
}

/// The in-circle determinant of `d` against the counter-clockwise triangle
/// `a`, `b`, `c`.
///
/// The value is positive when `d` lies strictly inside the circumcircle of the
/// triangle, zero when it lies on it and negative when it lies outside.
pub(super) fn incircle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> f64 {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];

    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady)
}

/// The circumcentre and squared circumradius of `a`, `b`, `c`.
///
/// Returns [`None`] when the three points are collinear to working precision,
/// in which case the circumcircle is not defined.
pub(super) fn circumcircle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> Option<([f64; 2], f64)> {
    let bx = b[0] - a[0];
    let by = b[1] - a[1];
    let cx = c[0] - a[0];
    let cy = c[1] - a[1];

    let det = 2.0 * (bx * cy - by * cx);
    let scale = (bx * bx + by * by)
        .max(cx * cx + cy * cy)
        .max(f64::MIN_POSITIVE);
    if det.abs() <= 1e-14 * scale {
        return None;
    }

    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    let ux = (cy * b2 - by * c2) / det;
    let uy = (bx * c2 - cx * b2) / det;

    Some(([a[0] + ux, a[1] + uy], ux * ux + uy * uy))
}

/// The normalised key of the undirected edge between `a` and `b`.
pub(super) fn edge_key(
    a: usize,
    b: usize,
) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

/// The area of the triangle `a`, `b`, `c`, always non-negative.
pub(super) fn triangle_area(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    0.5 * orient2d(a, b, c).abs()
}

/// The angle in radians at `p` in the corner `q`-`p`-`r`.
pub(super) fn angle_at(
    p: [f64; 2],
    q: [f64; 2],
    r: [f64; 2],
) -> f64 {
    let ux = q[0] - p[0];
    let uy = q[1] - p[1];
    let vx = r[0] - p[0];
    let vy = r[1] - p[1];

    let du = ux.hypot(uy);
    let dv = vx.hypot(vy);
    if du <= 0.0 || dv <= 0.0 {
        return 0.0;
    }

    let cosine = ((ux * vx + uy * vy) / (du * dv)).clamp(-1.0, 1.0);
    cosine.acos()
}

/// The smallest angle in radians of the triangle `a`, `b`, `c`.
pub(super) fn triangle_min_angle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    angle_at(a, b, c)
        .min(angle_at(b, c, a))
        .min(angle_at(c, a, b))
}

/// The centroid of the triangle `a`, `b`, `c`.
pub(super) fn centroid(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> [f64; 2] {
    [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0]
}

/// The midpoint of the segment `a`-`b`.
pub(super) fn midpoint(
    a: [f64; 2],
    b: [f64; 2],
) -> [f64; 2] {
    [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])]
}

/// The squared distance between `a` and `b`.
pub(super) fn distance2(
    a: [f64; 2],
    b: [f64; 2],
) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

/// Returns `t` with its vertices reordered so that it winds counter-clockwise.
pub(super) fn as_ccw(
    points: &[[f64; 2]],
    mut t: [usize; 3],
) -> [usize; 3] {
    if orient2d(points[t[0]], points[t[1]], points[t[2]]) < 0.0 {
        t.swap(1, 2);
    }
    t
}

/// Whether `p` lies strictly inside the counter-clockwise triangle `a`, `b`, `c`.
pub(super) fn point_in_triangle(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    eps: f64,
) -> bool {
    orient2d(a, b, p) > eps && orient2d(b, c, p) > eps && orient2d(c, a, p) > eps
}

/// Whether the two open segments `a`-`b` and `c`-`d` cross properly.
///
/// Endpoint contact is not reported; only a transversal crossing of the two
/// interiors counts.
pub(super) fn segments_cross(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> bool {
    let o1 = orient2d(a, b, c);
    let o2 = orient2d(a, b, d);
    let o3 = orient2d(c, d, a);
    let o4 = orient2d(c, d, b);

    let opposite = |u: f64, v: f64| (u > 0.0 && v < 0.0) || (u < 0.0 && v > 0.0);
    opposite(o1, o2) && opposite(o3, o4)
}

/// The distance from `p` to the segment `a`-`b`.
pub(super) fn point_segment_distance(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
) -> f64 {
    let abx = b[0] - a[0];
    let aby = b[1] - a[1];
    let len2 = abx * abx + aby * aby;
    if len2 <= 0.0 {
        return distance2(p, a).sqrt();
    }

    let t = (((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / len2).clamp(0.0, 1.0);
    let q = [a[0] + t * abx, a[1] + t * aby];
    distance2(p, q).sqrt()
}

/// Whether `p` lies inside the closed polygon `poly` by the even-odd rule.
pub(super) fn point_in_polygon(
    p: [f64; 2],
    poly: &[[f64; 2]],
) -> bool {
    if poly.len() < 3 {
        return false;
    }

    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let pi = poly[i];
        let pj = poly[j];
        if (pi[1] > p[1]) != (pj[1] > p[1]) {
            let t = (p[1] - pi[1]) / (pj[1] - pi[1]);
            if p[0] < pi[0] + t * (pj[0] - pi[0]) {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// The signed area of the polygon `poly`, positive when it winds counter-clockwise.
pub(super) fn polygon_signed_area(poly: &[[f64; 2]]) -> f64 {
    if poly.len() < 3 {
        return 0.0;
    }

    let mut acc = 0.0;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        acc += poly[j][0] * poly[i][1] - poly[i][0] * poly[j][1];
        j = i;
    }
    acc / 2.0
}

/// Whether `p` lies on the convex-hull boundary of `points` up to `tol`.
///
/// A point is on the hull when every other point sits on one side of the line
/// through it and at least two other points are not collinear with it.
pub(super) fn is_hull_edge(
    p: [f64; 2],
    q: [f64; 2],
    points: &[[f64; 2]],
    tol: f64,
) -> bool {
    let mut has_left = false;
    let mut has_right = false;
    for &r in points {
        let o = orient2d(p, q, r);
        if o > tol {
            has_left = true;
        } else if o < -tol {
            has_right = true;
        }
        if has_left && has_right {
            return false;
        }
    }
    has_left || has_right
}

/// The area of every triangle in `triangles`.
///
/// Triangles whose indices do not address `points` are reported as `0.0`.
#[must_use]
pub fn triangle_areas(
    points: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> Vec<f64> {
    let n = points.len();
    triangles
        .iter()
        .map(|t| {
            if t[0] < n && t[1] < n && t[2] < n {
                triangle_area(points[t[0]], points[t[1]], points[t[2]])
            } else {
                0.0
            }
        })
        .collect()
}

/// The smallest angle in degrees over all `triangles`.
///
/// Returns `0.0` for an empty mesh or a degenerate triangle.
#[must_use]
pub fn min_angle_deg(
    points: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> f64 {
    let n = points.len();
    let mut best = f64::INFINITY;
    for t in triangles {
        if t[0] >= n || t[1] >= n || t[2] >= n {
            continue;
        }
        let angle = triangle_min_angle(points[t[0]], points[t[1]], points[t[2]]);
        best = best.min(angle);
    }
    if best.is_finite() {
        best.to_degrees()
    } else {
        0.0
    }
}

/// The number of unique undirected edges spanned by `triangles`.
#[cfg(test)]
pub(super) fn count_edges(triangles: &[[usize; 3]]) -> usize {
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for t in triangles {
        for k in 0..3 {
            seen.insert(edge_key(t[k], t[(k + 1) % 3]));
        }
    }
    seen.len()
}

/// Triangulates the simple counter-clockwise polygon `poly` by ear clipping.
///
/// The result is a list of counter-clockwise triangles that cover the polygon.
/// A polygon that cannot be reduced to ears (because it is not simple) yields
/// whatever triangles could be clipped before the guard ran out.
pub(super) fn ear_clip(
    points: &[[f64; 2]],
    poly: &[usize],
) -> Vec<[usize; 3]> {
    let mut ring: Vec<usize> = poly.to_vec();
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        ring.pop();
    }
    if ring.len() < 3 {
        return Vec::new();
    }

    let mut out: Vec<[usize; 3]> = Vec::new();
    let mut guard = 0;
    let limit = 4 * ring.len() * ring.len() + 16;

    while ring.len() > 3 {
        guard += 1;
        if guard > limit {
            return out;
        }

        let n = ring.len();
        let mut clipped = false;
        for i in 0..n {
            let prev = ring[(i + n - 1) % n];
            let cur = ring[i];
            let next = ring[(i + 1) % n];

            if orient2d(points[prev], points[cur], points[next]) <= 1e-14 {
                continue;
            }

            let mut blocked = false;
            for &q in &ring {
                if q == prev || q == cur || q == next {
                    continue;
                }
                if point_in_triangle(points[q], points[prev], points[cur], points[next], 1e-14) {
                    blocked = true;
                    break;
                }
            }
            if blocked {
                continue;
            }

            out.push([prev, cur, next]);
            ring.remove(i);
            clipped = true;
            break;
        }

        if !clipped {
            // Drop a collinear vertex so the loop can still make progress.
            let n = ring.len();
            let mut removed = false;
            for i in 0..n {
                let prev = ring[(i + n - 1) % n];
                let cur = ring[i];
                let next = ring[(i + 1) % n];
                if orient2d(points[prev], points[cur], points[next]).abs() <= 1e-14 {
                    ring.remove(i);
                    removed = true;
                    break;
                }
            }
            if !removed {
                return out;
            }
        }
    }

    if ring.len() == 3 && orient2d(points[ring[0]], points[ring[1]], points[ring[2]]) > 0.0 {
        out.push([ring[0], ring[1], ring[2]]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::TriMesh2;
    use super::count_edges;
    use super::incircle;
    use super::min_angle_deg;
    use super::orient2d;
    use super::point_in_polygon;
    use super::triangle_areas;
    use super::triangle_min_angle;

    const UNIT_SQUARE: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

    #[test]
    fn orient2d_signs_follow_the_winding() {
        assert!(orient2d(UNIT_SQUARE[0], UNIT_SQUARE[1], UNIT_SQUARE[2]) > 0.0);
        assert!(orient2d(UNIT_SQUARE[0], UNIT_SQUARE[2], UNIT_SQUARE[1]) < 0.0);
        assert!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]).abs() < 1e-15);
    }

    #[test]
    fn incircle_detects_the_interior_of_the_circumcircle() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        // The circumcentre of the unit right triangle is (0.5, 0.5) with
        // squared radius 0.5, so (0.25, 0.25) is inside and (2, 2) outside.
        assert!(incircle(a, b, c, [0.25, 0.25]) > 0.0);
        assert!(incircle(a, b, c, [1.0, 1.0]).abs() < 1e-15);
        assert!(incircle(a, b, c, [2.0, 2.0]) < 0.0);
    }

    #[test]
    fn triangle_areas_and_min_angle_match_the_square() {
        let points = UNIT_SQUARE.to_vec();
        let triangles = [[0, 1, 2], [0, 2, 3]];
        let areas = triangle_areas(&points, &triangles);
        assert!((areas[0] - 0.5).abs() < 1e-12);
        assert!((areas[1] - 0.5).abs() < 1e-12);
        assert!((min_angle_deg(&points, &triangles) - 45.0).abs() < 1e-9);
        assert_eq!(count_edges(&triangles), 5);
    }

    #[test]
    fn degenerate_triangle_has_zero_minimum_angle() {
        let angle = triangle_min_angle([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]);
        assert!(angle < 1e-12);
        assert_eq!(min_angle_deg(&UNIT_SQUARE, &[]), 0.0);
    }

    #[test]
    fn point_in_polygon_accepts_the_interior_only() {
        let poly = UNIT_SQUARE.to_vec();
        assert!(point_in_polygon([0.5, 0.5], &poly));
        assert!(!point_in_polygon([1.5, 0.5], &poly));
        assert!(!point_in_polygon([-0.1, -0.1], &poly));
    }

    #[test]
    fn empty_mesh_is_the_default() {
        let mesh = TriMesh2::default();
        assert_eq!(mesh.points.len(), 0);
        assert_eq!(mesh.triangles.len(), 0);
    }
}
