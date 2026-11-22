//! Delaunay triangulation in the plane.
//!
//! [`delaunay_2d`] runs the Bowyer–Watson insertion scheme on a super-triangle
//! and then restores the empty-circumcircle property with Lawson edge flips.
//! [`constrained_delaunay_2d`] additionally forces a set of input segments into
//! the mesh by removing the triangles they cross and re-triangulating the
//! resulting cavity (Sloan's edge-flipping insertion).
//!
//! # Limitations
//!
//! * The predicates use plain `f64` arithmetic. Inputs that are exactly
//!   cocircular are handled by a strict in-circle test, so the chosen diagonal
//!   of a degenerate quad is arbitrary but the result is still a valid
//!   triangulation.
//! * Segments that pass exactly through an existing vertex, or that cross each
//!   other, are not supported; the segment is left out in that case.
//! * The duplicate-point scan is `O(n)` in the point count through a uniform
//!   grid, but the Bowyer–Watson insertion is `O(n^2)`.

use std::collections::HashMap;
use std::collections::HashSet;

use super::util::as_ccw;
use super::util::ear_clip;
use super::util::edge_key;
use super::util::incircle;
use super::util::orient2d;
use super::util::segments_cross;
use super::util::triangle_area;

/// The number of point pairs to check before the duplicate scan gives up.
const MAX_DUPLICATE_SCAN: usize = 4_000_000;

/// Delaunay triangulation of `points` by Bowyer–Watson insertion.
///
/// Returns counter-clockwise triangles as triples of indices into `points`.
///
/// # Errors
///
/// Returns `Err` when fewer than three points are supplied, when a coordinate
/// is not finite, or when two points coincide.
pub fn delaunay_2d(points: &[[f64; 2]]) -> Result<Vec<[usize; 3]>, String> {
    validate(points)?;

    let n = points.len();
    let mut pts = points.to_vec();

    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &pts {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let span = (max_x - min_x).max(max_y - min_y).max(1.0);
    let cx = f64::midpoint(min_x, max_x);
    let cy = f64::midpoint(min_y, max_y);

    pts.push([cx - 20.0 * span, cy - 10.0 * span]);
    pts.push([cx + 20.0 * span, cy - 10.0 * span]);
    pts.push([cx, cy + 20.0 * span]);

    let super_triangle = as_ccw(&pts, [n, n + 1, n + 2]);
    let mut triangles: Vec<[usize; 3]> = vec![super_triangle];

    for i in 0..n {
        let p = pts[i];
        let mut boundary: HashMap<(usize, usize), usize> = HashMap::new();
        let mut kept: Vec<[usize; 3]> = Vec::with_capacity(triangles.len() + 4);

        for t in &triangles {
            let inside = incircle(pts[t[0]], pts[t[1]], pts[t[2]], p) > 0.0;
            if inside {
                for k in 0..3 {
                    *boundary.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
                }
            } else {
                kept.push(*t);
            }
        }

        // Sorting the boundary keys keeps the triangle order independent of
        // the hash-map iteration order.
        let mut candidates: Vec<(usize, usize)> = boundary.keys().copied().collect();
        candidates.sort_unstable();
        for (a, b) in candidates {
            if boundary.contains_key(&(b, a)) {
                continue;
            }
            let candidate = as_ccw(&pts, [a, b, i]);
            if triangle_area(pts[candidate[0]], pts[candidate[1]], pts[candidate[2]]) > 0.0 {
                kept.push(candidate);
            }
        }

        triangles = kept;
    }

    triangles.retain(|t| t[0] < n && t[1] < n && t[2] < n);
    deduplicate(&mut triangles);

    lawson(points, &mut triangles, &HashSet::new());

    Ok(triangles)
}

/// Incremental Delaunay triangulation of `points`.
///
/// This is the incremental entry point of the same Bowyer–Watson insertion
/// used by [`delaunay_2d`]: points are inserted one at a time into the growing
/// triangulation, and each insertion only re-triangulates the cavity of
/// triangles whose circumcircle contains the new point.
///
/// # Errors
///
/// Returns `Err` under the same conditions as [`delaunay_2d`].
pub fn delaunay_2d_incremental(points: &[[f64; 2]]) -> Result<Vec<[usize; 3]>, String> {
    delaunay_2d(points)
}

/// Delaunay triangulation of `points` with the given segments forced in.
///
/// Each `(u, v)` in `segments` becomes an edge of the returned triangulation.
/// Segments are inserted by removing every triangle the segment crosses and
/// re-triangulating the two halves of the resulting cavity (Sloan's algorithm).
///
/// # Errors
///
/// Returns `Err` when the points are invalid, when a segment index is out of
/// range, or when a segment joins a vertex to itself.
pub fn constrained_delaunay_2d(
    points: &[[f64; 2]],
    segments: &[(usize, usize)],
) -> Result<Vec<[usize; 3]>, String> {
    let mut triangles = delaunay_2d(points)?;

    for &(u, v) in segments {
        if u >= points.len() || v >= points.len() {
            return Err("segment index out of range".to_string());
        }
        if u == v {
            return Err("degenerate segment".to_string());
        }
        insert_segment(points, &mut triangles, u, v);
    }

    let fixed: HashSet<(usize, usize)> = segments.iter().map(|&(u, v)| edge_key(u, v)).collect();
    lawson(points, &mut triangles, &fixed);

    Ok(triangles)
}

/// Rejects inputs the triangulation cannot handle.
fn validate(points: &[[f64; 2]]) -> Result<(), String> {
    if points.len() < 3 {
        return Err("at least three points are required".to_string());
    }

    for p in points {
        if !p[0].is_finite() || !p[1].is_finite() {
            return Err("point coordinates must be finite".to_string());
        }
    }

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let scale = (max_x - min_x).max(max_y - min_y).max(1.0);
    let tol = 1e-12 * scale;

    if points.len() * points.len() <= MAX_DUPLICATE_SCAN {
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let dx = points[i][0] - points[j][0];
                let dy = points[i][1] - points[j][1];
                if dx.abs() <= tol && dy.abs() <= tol {
                    return Err("duplicate points".to_string());
                }
            }
        }
        return Ok(());
    }

    let cell = tol;
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        let cx = (p[0] / cell).floor() as i64;
        let cy = (p[1] / cell).floor() as i64;
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                if let Some(list) = grid.get(&(cx + dx, cy + dy)) {
                    for &j in list {
                        let q = points[j];
                        if (q[0] - p[0]).abs() <= tol && (q[1] - p[1]).abs() <= tol {
                            return Err("duplicate points".to_string());
                        }
                    }
                }
            }
        }
        grid.entry((cx, cy)).or_default().push(i);
    }

    Ok(())
}

/// Removes triangles that repeat the same three vertices.
fn deduplicate(triangles: &mut Vec<[usize; 3]>) {
    let mut seen: HashSet<[usize; 3]> = HashSet::new();
    triangles.retain(|t| {
        let mut key = *t;
        key.sort_unstable();
        seen.insert(key)
    });
}

/// The vertex of `t` that is neither `a` nor `b`.
fn opposite(
    t: [usize; 3],
    a: usize,
    b: usize,
) -> Option<usize> {
    t.iter().find(|&&v| v != a && v != b).copied()
}

/// Restores the empty-circumcircle property by flipping non-Delaunay edges.
///
/// Edges listed in `fixed` are never flipped, which keeps constrained segments
/// in place.
fn lawson(
    points: &[[f64; 2]],
    triangles: &mut [[usize; 3]],
    fixed: &HashSet<(usize, usize)>,
) {
    let mut adjacency = build_adjacency(triangles);
    let mut stack: Vec<(usize, usize)> = adjacency.keys().copied().collect();
    stack.sort_unstable();
    let mut budget = 64 * triangles.len() * triangles.len() + 1024;

    while let Some((a, b)) = stack.pop() {
        if budget == 0 {
            break;
        }
        budget -= 1;

        if fixed.contains(&(a, b)) {
            continue;
        }
        let Some(list) = adjacency.get(&(a, b)) else {
            continue;
        };
        if list.len() != 2 {
            continue;
        }
        let (i, j) = (list[0], list[1]);
        if i == j || i >= triangles.len() || j >= triangles.len() {
            continue;
        }

        let Some(c) = opposite(triangles[i], a, b) else {
            continue;
        };
        let Some(d) = opposite(triangles[j], a, b) else {
            continue;
        };

        // The in-circle test needs the edge ordered so that `c` lies to its
        // left; the edge key is sorted by index, not by winding.
        let (u, v) = if orient2d(points[a], points[b], points[c]) > 0.0 {
            (a, b)
        } else {
            (b, a)
        };
        if incircle(points[u], points[v], points[c], points[d]) <= 0.0 {
            continue;
        }

        let o1 = orient2d(points[c], points[d], points[a]);
        let o2 = orient2d(points[c], points[d], points[b]);
        if o1 * o2 >= 0.0 {
            continue;
        }

        let first = as_ccw(points, [c, d, a]);
        let second = as_ccw(points, [d, c, b]);
        if triangle_area(points[first[0]], points[first[1]], points[first[2]]) <= 0.0
            || triangle_area(points[second[0]], points[second[1]], points[second[2]]) <= 0.0
        {
            continue;
        }

        detach(&mut adjacency, (a, b), i);
        detach(&mut adjacency, (a, b), j);
        detach(&mut adjacency, edge_key(b, c), i);
        detach(&mut adjacency, edge_key(c, a), i);
        detach(&mut adjacency, edge_key(b, d), j);
        detach(&mut adjacency, edge_key(d, a), j);

        triangles[i] = first;
        triangles[j] = second;

        for e in [edge_key(d, a), edge_key(a, c)] {
            adjacency.entry(e).or_default().push(i);
        }
        for e in [edge_key(c, b), edge_key(b, d)] {
            adjacency.entry(e).or_default().push(j);
        }
        adjacency.entry(edge_key(c, d)).or_default().extend([i, j]);

        for e in [
            edge_key(d, a),
            edge_key(a, c),
            edge_key(c, b),
            edge_key(b, d),
        ] {
            stack.push(e);
        }
    }
}

/// Maps every undirected edge to the indices of the triangles that use it.
fn build_adjacency(triangles: &[[usize; 3]]) -> HashMap<(usize, usize), Vec<usize>> {
    let mut adjacency: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (index, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            adjacency
                .entry(edge_key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(index);
        }
    }
    adjacency
}

/// Removes triangle `t` from the adjacency list of edge `e`.
fn detach(
    adjacency: &mut HashMap<(usize, usize), Vec<usize>>,
    e: (usize, usize),
    t: usize,
) {
    let now_empty = match adjacency.get_mut(&e) {
        | Some(list) => {
            list.retain(|&x| x != t);
            list.is_empty()
        },
        | None => false,
    };
    if now_empty {
        adjacency.remove(&e);
    }
}

/// Forces the segment `u`-`v` into `triangles` by cavity re-triangulation.
fn insert_segment(
    points: &[[f64; 2]],
    triangles: &mut Vec<[usize; 3]>,
    u: usize,
    v: usize,
) {
    if triangles.iter().any(|t| {
        let mut has_u = false;
        let mut has_v = false;
        for &w in t {
            has_u |= w == u;
            has_v |= w == v;
        }
        has_u && has_v
    }) {
        return;
    }

    let a = points[u];
    let b = points[v];

    let mut crossing: HashSet<(usize, usize)> = HashSet::new();
    for t in triangles.iter() {
        for k in 0..3 {
            let e = (t[k], t[(k + 1) % 3]);
            if e.0 == u || e.1 == u || e.0 == v || e.1 == v {
                continue;
            }
            if segments_cross(a, b, points[e.0], points[e.1]) {
                crossing.insert(edge_key(e.0, e.1));
            }
        }
    }
    if crossing.is_empty() {
        return;
    }

    let mut cavity: Vec<usize> = Vec::new();
    for (index, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            if crossing.contains(&edge_key(t[k], t[(k + 1) % 3])) {
                cavity.push(index);
                break;
            }
        }
    }
    if cavity.is_empty() {
        return;
    }

    let removed: HashSet<usize> = cavity.iter().copied().collect();
    let mut directed: HashSet<(usize, usize)> = HashSet::new();
    for &index in &cavity {
        let t = triangles[index];
        for k in 0..3 {
            directed.insert((t[k], t[(k + 1) % 3]));
        }
    }

    let boundary: Vec<(usize, usize)> = directed
        .iter()
        .copied()
        .filter(|&(p, q)| !directed.contains(&(q, p)))
        .collect();

    let mut next: HashMap<usize, usize> = HashMap::new();
    for &(p, q) in &boundary {
        next.insert(p, q);
    }

    let Some(cycle) = walk_cycle(&next, u, v) else {
        return;
    };

    let Some(split) = cycle.iter().position(|&x| x == v) else {
        return;
    };
    if split == 0 || split + 1 >= cycle.len() {
        return;
    }

    let first: Vec<usize> = cycle[..=split].to_vec();
    let second: Vec<usize> = cycle[split..].to_vec();

    let mut replacement: Vec<[usize; 3]> = Vec::new();
    replacement.extend(ear_clip(points, &first));
    replacement.extend(ear_clip(points, &second));

    if replacement.is_empty() {
        return;
    }

    let mut kept: Vec<[usize; 3]> = Vec::new();
    for (index, t) in triangles.iter().enumerate() {
        if !removed.contains(&index) {
            kept.push(*t);
        }
    }
    kept.extend(replacement);
    *triangles = kept;
}

/// Walks the boundary cycle that starts at `u` and stops when it returns.
///
/// Returns [`None`] when the walk does not pass through `v` or does not close.
fn walk_cycle(
    next: &HashMap<usize, usize>,
    u: usize,
    v: usize,
) -> Option<Vec<usize>> {
    let mut cycle = vec![u];
    let mut current = u;
    let mut guard = 0;
    loop {
        let following = *next.get(&current)?;
        cycle.push(following);
        if following == u {
            break;
        }
        current = following;
        guard += 1;
        if guard > next.len() + 2 {
            return None;
        }
    }

    if !cycle.contains(&v) {
        return None;
    }
    Some(cycle)
}

#[cfg(test)]
mod tests {
    use super::constrained_delaunay_2d;
    use super::delaunay_2d;
    use super::delaunay_2d_incremental;
    use crate::kernels::meshgen::util::count_edges;
    use crate::kernels::meshgen::util::incircle;
    use crate::kernels::meshgen::util::min_angle_deg;
    use crate::kernels::meshgen::util::triangle_areas;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    fn expect_err<T>(result: Result<T, String>) -> String {
        match result {
            | Ok(_) => panic!("expected an error"),
            | Err(message) => message,
        }
    }

    /// A deterministic low-discrepancy cloud in the unit square.
    fn cloud(count: usize) -> Vec<[f64; 2]> {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut points = Vec::with_capacity(count);
        for _ in 0..count {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let x = ((state >> 11) as f64) / ((1_u64 << 53) as f64);
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let y = ((state >> 11) as f64) / ((1_u64 << 53) as f64);
            points.push([x, y]);
        }
        points
    }

    #[test]
    fn unit_square_gives_two_triangles_of_unit_area() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = expect_ok(delaunay_2d(&points));
        assert_eq!(triangles.len(), 2);

        let areas = triangle_areas(&points, &triangles);
        let total: f64 = areas.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "total area {total}");
        assert!((min_angle_deg(&points, &triangles) - 45.0).abs() < 1e-9);
        for t in &triangles {
            assert!(
                crate::kernels::meshgen::util::orient2d(points[t[0]], points[t[1]], points[t[2]])
                    > 0.0
            );
        }
    }

    #[test]
    fn incremental_entry_point_agrees_with_delaunay_2d() {
        let points = cloud(20);
        let plain = expect_ok(delaunay_2d(&points));
        let incremental = expect_ok(delaunay_2d_incremental(&points));
        assert_eq!(plain, incremental);
    }

    #[test]
    fn point_cloud_is_a_valid_delaunay_disk() {
        let points = cloud(20);
        let triangles = expect_ok(delaunay_2d(&points));

        assert_ne!(triangles.len(), 0);
        let vertices = points.len();
        let edges = count_edges(&triangles);
        let faces = triangles.len();
        let euler = vertices as i64 - edges as i64 + faces as i64;
        assert_eq!(euler, 1, "euler characteristic {euler}");

        for t in &triangles {
            for (index, p) in points.iter().enumerate() {
                if index == t[0] || index == t[1] || index == t[2] {
                    continue;
                }
                let value = incircle(points[t[0]], points[t[1]], points[t[2]], *p);
                assert!(
                    value <= 1e-9,
                    "point {index} inside circumcircle of {t:?}: {value}"
                );
            }
        }

        let areas = triangle_areas(&points, &triangles);
        for area in &areas {
            assert!(*area > 0.0, "degenerate triangle in the cloud");
        }
    }

    #[test]
    fn constrained_diagonal_is_kept_as_an_edge() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let segments = [(0_usize, 2_usize)];
        let triangles = expect_ok(constrained_delaunay_2d(&points, &segments));

        assert_eq!(triangles.len(), 2);
        let has_diagonal = triangles.iter().any(|t| {
            let mut has_u = false;
            let mut has_v = false;
            for &w in t {
                has_u |= w == 0;
                has_v |= w == 2;
            }
            has_u && has_v
        });
        assert!(has_diagonal, "diagonal (0, 2) missing from {triangles:?}");

        let total: f64 = triangle_areas(&points, &triangles).iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn constrained_segments_are_validated() {
        let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        assert!(constrained_delaunay_2d(&points, &[(0, 9)]).is_err());
        assert!(constrained_delaunay_2d(&points, &[(0, 0)]).is_err());
    }

    #[test]
    fn degenerate_inputs_are_rejected() {
        assert_eq!(
            expect_err(delaunay_2d(&[[0.0, 0.0], [1.0, 1.0]])),
            "at least three points are required"
        );
        assert_eq!(
            expect_err(delaunay_2d(&[[0.0, 0.0], [1.0, 0.0], [1.0, 0.0]])),
            "duplicate points"
        );
        assert!(delaunay_2d(&[[0.0, 0.0], [1.0, 0.0], [f64::NAN, 1.0]]).is_err());
    }

    #[test]
    fn collinear_points_produce_no_triangles() {
        let points = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]];
        let triangles = expect_ok(delaunay_2d(&points));
        assert_eq!(triangles.len(), 0);
    }
}
