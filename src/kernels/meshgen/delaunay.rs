//! Delaunay triangulation in the plane.
//!
//! [`delaunay_2d`] runs the Bowyer–Watson insertion scheme on a super-triangle
//! using exact adaptive-precision geometric predicates with Simulation of
//! Simplicity (SoS) tie-breaking for cocircular points, and then restores the
//! empty-circumcircle property with Lawson edge flips.
//!
//! [`constrained_delaunay_2d`] and [`constrained_delaunay_2d_mesh`] force a
//! set of input segments into the mesh by first preprocessing the segments
//! (splitting segments at any existing collinear vertices, and in
//! [`constrained_delaunay_2d_mesh`] also splitting mutually crossing segments
//! at newly inserted Steiner intersection vertices) before enforcing every
//! sub-segment by cavity re-triangulation (Sloan's algorithm).
//!
//! # Limitations
//!
//! * The duplicate-point scan is `O(n)` in the point count through a uniform
//!   grid, while the 2D Bowyer–Watson insertion is `O(n^2)`.
//! * [`constrained_delaunay_2d`] returns triangle indices into the caller's
//!   fixed point slice, so mutually crossing segments require
//!   [`constrained_delaunay_2d_mesh`], which returns the augmented [`TriMesh2`]
//!   including any inserted Steiner vertices.

use std::collections::HashMap;
use std::collections::HashSet;

use super::util::TriMesh2;
use super::util::as_ccw;
use super::util::distance2;
use super::util::ear_clip;
use super::util::edge_key;
use super::util::incircle_sos;
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
    // Circumcircle of every triangle, in the same order as `triangles`. It is
    // a floating-point pre-filter: only points near or on a circumcircle are
    // decided by the exact predicate, which is what keeps the result identical
    // to a pure `incircle_sos` scan while cutting its cost by an order of
    // magnitude on a large cloud.
    let mut circles: Vec<Option<(f64, f64, f64)>> = vec![circumcircle(&pts, super_triangle)];

    let mut boundary: HashMap<(usize, usize), usize> = HashMap::new();
    let mut candidates: Vec<(usize, usize)> = Vec::new();
    let mut kept: Vec<[usize; 3]> = Vec::with_capacity(triangles.len() + 4);
    let mut kept_circles: Vec<Option<(f64, f64, f64)>> = Vec::with_capacity(triangles.len() + 4);

    for i in 0..n {
        boundary.clear();
        kept.clear();
        kept_circles.clear();

        for (index, t) in triangles.iter().enumerate() {
            let inside = match circles[index] {
                | Some(circle) => point_in_circumcircle(circle, pts[i], &pts, *t, i),
                | None => incircle_sos(&pts, t[0], t[1], t[2], i),
            };
            if inside {
                for k in 0..3 {
                    *boundary.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
                }
            } else {
                kept.push(*t);
                kept_circles.push(circles[index]);
            }
        }

        // Sorting the boundary keys keeps the triangle order independent of
        // the hash-map iteration order.
        candidates.clear();
        candidates.extend(boundary.keys().copied());
        candidates.sort_unstable();
        for &(a, b) in &candidates {
            if boundary.contains_key(&(b, a)) {
                continue;
            }
            let candidate = as_ccw(&pts, [a, b, i]);
            if triangle_area(pts[candidate[0]], pts[candidate[1]], pts[candidate[2]]) > 0.0 {
                kept.push(candidate);
                kept_circles.push(circumcircle(&pts, candidate));
            }
        }

        std::mem::swap(&mut triangles, &mut kept);
        std::mem::swap(&mut circles, &mut kept_circles);
    }

    triangles.retain(|t| t[0] < n && t[1] < n && t[2] < n);
    deduplicate(&mut triangles);

    lawson(points, &mut triangles, &HashSet::new());

    Ok(triangles)
}

/// Floating-point circumcircle `(centre_x, centre_y, radius^2)` of a triangle,
/// or `None` when the triangle is too close to degenerate for the circle to be
/// trusted.
fn circumcircle(
    points: &[[f64; 2]],
    t: [usize; 3],
) -> Option<(f64, f64, f64)> {
    let (ax, ay) = (points[t[0]][0], points[t[0]][1]);
    let (bx, by) = (points[t[1]][0], points[t[1]][1]);
    let (cx, cy) = (points[t[2]][0], points[t[2]][1]);
    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    let edge = ((bx - ax) * (bx - ax) + (by - ay) * (by - ay))
        .max((cx - bx) * (cx - bx) + (cy - by) * (cy - by))
        .max((ax - cx) * (ax - cx) + (ay - cy) * (ay - cy));
    // `d` is twice the signed area; below this ratio the centre is
    // ill-conditioned and the exact predicate is cheaper than trusting it.
    if !d.is_finite() || d.abs() <= 1e-4 * edge {
        return None;
    }
    let a2 = ax * ax + ay * ay;
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    let ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    let uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    let r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy);
    if ux.is_finite() && uy.is_finite() && r2.is_finite() && r2 >= 0.0 {
        Some((ux, uy, r2))
    } else {
        None
    }
}

/// Whether `p` lies inside the cached circumcircle, falling back to the exact
/// predicate only when the floating-point decision is within its error bound.
fn point_in_circumcircle(
    circle: (f64, f64, f64),
    p: [f64; 2],
    points: &[[f64; 2]],
    t: [usize; 3],
    index: usize,
) -> bool {
    let (ux, uy, r2) = circle;
    let dx = p[0] - ux;
    let dy = p[1] - uy;
    let d2 = dx * dx + dy * dy;
    let ux_squared = ux * ux;
    let uy_squared = uy * uy;
    let px_squared = p[0] * p[0];
    let py_squared = p[1] * p[1];
    let magnitude = (d2 + r2 + ux_squared + uy_squared + px_squared + py_squared).max(1.0);
    let tolerance = 1e-10 * magnitude;
    if d2 < r2 - tolerance {
        return true;
    }
    if d2 > r2 + tolerance {
        return false;
    }
    incircle_sos(points, t[0], t[1], t[2], index)
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
/// Each `(u, v)` in `segments` becomes a chain of constrained edges in the
/// returned triangulation: any segment passing through existing collinear
/// vertices is automatically split at those vertices before enforcement.
///
/// # Errors
///
/// Returns `Err` when the points are invalid, when a segment index is out of
/// range, when a segment joins a vertex to itself, or when two segments cross
/// transversely (in which case [`constrained_delaunay_2d_mesh`] should be used
/// so the intersection Steiner points can be returned).
pub fn constrained_delaunay_2d(
    points: &[[f64; 2]],
    segments: &[(usize, usize)],
) -> Result<Vec<[usize; 3]>, String> {
    validate(points)?;
    let (_, sub_segments) = preprocess_segments(points, segments, false, usize::MAX)?;
    enforce_constraints(points, &sub_segments)
}

/// Constrained Delaunay triangulation of `points` and `segments`, returning
/// the full [`TriMesh2`] including any Steiner vertices inserted where input
/// segments cross.
///
/// Segments passing through existing collinear vertices are split at those
/// vertices, and mutually crossing segments are split at their intersection
/// points (deduplicated against existing vertices within tolerance). Every
/// input segment is guaranteed to be represented as a chain of constrained
/// edges in the output mesh.
///
/// # Errors
///
/// Returns `Err` when the points are invalid, when a segment index is out of
/// range, or when a segment joins a vertex to itself.
pub fn constrained_delaunay_2d_mesh(
    points: &[[f64; 2]],
    segments: &[(usize, usize)],
) -> Result<TriMesh2, String> {
    constrained_delaunay_bounded(points, segments, usize::MAX)
}

pub(crate) fn constrained_delaunay_bounded(
    points: &[[f64; 2]],
    segments: &[(usize, usize)],
    vertex_budget: usize,
) -> Result<TriMesh2, String> {
    if points.len() > vertex_budget {
        return Err("constrained Delaunay vertex budget exceeded".into());
    }
    validate(points)?;
    let (points, sub_segments) = preprocess_segments(points, segments, true, vertex_budget)?;
    let triangles = enforce_constraints(&points, &sub_segments)?;
    Ok(TriMesh2 { points, triangles })
}

fn enforce_constraints(
    points: &[[f64; 2]],
    constraints: &[(usize, usize)],
) -> Result<Vec<[usize; 3]>, String> {
    let mut triangles = delaunay_2d(points)?;
    for &(u, v) in constraints {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        insert_segment(points, &mut triangles, u, v);
    }
    let fixed: HashSet<_> = constraints.iter().map(|&(u, v)| edge_key(u, v)).collect();
    lawson(points, &mut triangles, &fixed);
    let edges = build_adjacency(&triangles);
    if fixed.iter().any(|edge| !edges.contains_key(edge)) {
        return Err("constrained Delaunay could not retain every input constraint".into());
    }
    Ok(triangles)
}

/// Preprocesses constraint segments by splitting segments at existing collinear
/// vertices and (when `allow_steiner` is `true`) at mutual segment crossings.
fn preprocess_segments(
    points: &[[f64; 2]],
    segments: &[(usize, usize)],
    allow_steiner: bool,
    vertex_budget: usize,
) -> Result<(Vec<[f64; 2]>, Vec<(usize, usize)>), String> {
    for &(u, v) in segments {
        if u >= points.len() || v >= points.len() {
            return Err("segment index out of range".to_string());
        }
        if u == v {
            return Err("degenerate segment".to_string());
        }
    }

    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in points {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let scale = (max_x - min_x).max(max_y - min_y).max(1.0);
    let pt_tol = 1e-11 * scale;
    let pt_tol2 = pt_tol * pt_tol;
    let collinear_tol = 1e-12 * scale;

    let mut pts = points.to_vec();
    let mut segs = dedup_segments(segments);

    for _ in 0..64 {
        // 1. Split every segment at any existing vertex lying strictly on it.
        let mut after_collinear: Vec<(usize, usize)> = Vec::with_capacity(segs.len());
        for &(u, v) in &segs {
            let a = pts[u];
            let b = pts[v];
            let len2 = distance2(a, b);
            if len2 <= 0.0 {
                continue;
            }
            let len = len2.sqrt();
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];

            let mut interior: Vec<usize> = Vec::new();
            for (w, &p) in pts.iter().enumerate() {
                if w == u || w == v {
                    continue;
                }
                let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2;
                if t <= 0.0 || t >= 1.0 {
                    continue;
                }
                if distance2(p, a) <= pt_tol2 || distance2(p, b) <= pt_tol2 {
                    continue;
                }
                if orient2d(a, b, p).abs() <= collinear_tol * len {
                    interior.push(w);
                }
            }

            if interior.is_empty() {
                after_collinear.push((u, v));
            } else {
                after_collinear.extend(split_segment_by_vertices(&pts, u, v, &interior));
            }
        }
        segs = dedup_segments(&after_collinear);

        // 2. Detect and split mutually crossing segments.
        let mut splits: Vec<Vec<usize>> = vec![Vec::new(); segs.len()];
        let mut any_crossing = false;

        for i in 0..segs.len() {
            if i.is_multiple_of(16) {
                crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            }
            let (a, b) = segs[i];
            for j in (i + 1)..segs.len() {
                let (c, d) = segs[j];
                if a == c || a == d || b == c || b == d {
                    continue;
                }
                if segments_cross(pts[a], pts[b], pts[c], pts[d]) {
                    if !allow_steiner {
                        return Err(
                            "crossing segments require Steiner points; use constrained_delaunay_2d_mesh"
                                .to_string(),
                        );
                    }
                    any_crossing = true;
                    let inter = segment_intersection(pts[a], pts[b], pts[c], pts[d]);
                    let idx = find_or_insert_vertex(&mut pts, inter, pt_tol2, vertex_budget)
                        .ok_or("constrained Delaunay vertex budget exceeded")?;
                    if idx != a && idx != b {
                        splits[i].push(idx);
                    }
                    if idx != c && idx != d {
                        splits[j].push(idx);
                    }
                }
            }
        }

        if !any_crossing {
            break;
        }

        let mut next_segs: Vec<(usize, usize)> = Vec::with_capacity(segs.len() * 2);
        for (i, &(u, v)) in segs.iter().enumerate() {
            if splits[i].is_empty() {
                next_segs.push((u, v));
            } else {
                next_segs.extend(split_segment_by_vertices(&pts, u, v, &splits[i]));
            }
        }
        segs = dedup_segments(&next_segs);
    }

    Ok((pts, segs))
}

/// Deduplicates undirected segments and removes zero-length `(u, u)` entries.
fn dedup_segments(raw: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut seen: HashSet<(usize, usize)> = HashSet::with_capacity(raw.len());
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(raw.len());
    for &(u, v) in raw {
        if u == v {
            continue;
        }
        let key = edge_key(u, v);
        if seen.insert(key) {
            out.push(key);
        }
    }
    out
}

/// Splits the segment `(u, v)` along `u -> v` at all given interior vertices.
fn split_segment_by_vertices(
    pts: &[[f64; 2]],
    u: usize,
    v: usize,
    verts: &[usize],
) -> Vec<(usize, usize)> {
    let a = pts[u];
    let b = pts[v];
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return Vec::new();
    }

    let mut ordered: Vec<(f64, usize)> = verts
        .iter()
        .copied()
        .filter(|&w| w != u && w != v)
        .map(|w| {
            let p = pts[w];
            let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2;
            (t, w)
        })
        .collect();
    ordered.sort_by(|x, y| x.0.total_cmp(&y.0).then_with(|| x.1.cmp(&y.1)));

    let mut chain: Vec<usize> = Vec::with_capacity(ordered.len() + 2);
    chain.push(u);
    for (_, w) in ordered {
        if chain.last().copied() != Some(w) {
            chain.push(w);
        }
    }
    if chain.last().copied() != Some(v) {
        chain.push(v);
    }

    chain.windows(2).map(|pair| (pair[0], pair[1])).collect()
}

/// Computes the intersection point of two transversely crossing segments.
fn segment_intersection(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> [f64; 2] {
    let abx = b[0] - a[0];
    let aby = b[1] - a[1];
    let cdx = d[0] - c[0];
    let cdy = d[1] - c[1];
    let denom = abx * cdy - aby * cdx;
    let t = (((c[0] - a[0]) * cdy - (c[1] - a[1]) * cdx) / denom).clamp(0.0, 1.0);
    [a[0] + t * abx, a[1] + t * aby]
}

/// Returns the index of an existing vertex within `tol2` of `p`, or appends `p`.
fn find_or_insert_vertex(
    pts: &mut Vec<[f64; 2]>,
    p: [f64; 2],
    tol2: f64,
    vertex_budget: usize,
) -> Option<usize> {
    for (idx, &q) in pts.iter().enumerate() {
        if distance2(q, p) <= tol2 {
            return Some(idx);
        }
    }
    let idx = pts.len();
    if idx >= vertex_budget {
        return None;
    }
    pts.push(p);
    Some(idx)
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
        if !incircle_sos(points, u, v, c, d) {
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
    use super::constrained_delaunay_2d_mesh;
    use super::delaunay_2d;
    use super::delaunay_2d_incremental;
    use crate::kernels::meshgen::util::count_edges;
    use crate::kernels::meshgen::util::edge_key;
    use crate::kernels::meshgen::util::incircle;
    use crate::kernels::meshgen::util::min_angle_deg;
    use crate::kernels::meshgen::util::orient2d;
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

    fn has_edge(
        triangles: &[[usize; 3]],
        u: usize,
        v: usize,
    ) -> bool {
        let target = edge_key(u, v);
        triangles.iter().any(|t| {
            edge_key(t[0], t[1]) == target
                || edge_key(t[1], t[2]) == target
                || edge_key(t[2], t[0]) == target
        })
    }

    #[test]
    fn constrained_segment_splits_at_collinear_vertices() {
        // Segment (0, 3) passes through two existing interior vertices 1 and 2.
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [1.5, 2.0],
            [1.5, -2.0],
        ];
        let triangles = expect_ok(constrained_delaunay_2d(&points, &[(0, 3)]));

        assert!(has_edge(&triangles, 0, 1), "sub-edge (0, 1) missing");
        assert!(has_edge(&triangles, 1, 2), "sub-edge (1, 2) missing");
        assert!(has_edge(&triangles, 2, 3), "sub-edge (2, 3) missing");

        for t in &triangles {
            assert!(orient2d(points[t[0]], points[t[1]], points[t[2]]) > 0.0);
        }
        let total: f64 = triangle_areas(&points, &triangles).iter().sum();
        assert!((total - 6.0).abs() < 1e-12, "total area {total}");
    }

    #[test]
    fn constrained_crossing_x_inserts_steiner_point() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let segments = [(0, 2), (1, 3)];
        let mesh = expect_ok(constrained_delaunay_2d_mesh(&points, &segments));

        assert_eq!(mesh.points.len(), 5, "expected one Steiner intersection");
        let center = 4_usize;
        assert!((mesh.points[center][0] - 0.5).abs() < 1e-12);
        assert!((mesh.points[center][1] - 0.5).abs() < 1e-12);

        for corner in 0..4 {
            assert!(
                has_edge(&mesh.triangles, corner, center),
                "sub-edge ({corner}, {center}) missing from {:?}",
                mesh.triangles
            );
        }

        for t in &mesh.triangles {
            assert!(orient2d(mesh.points[t[0]], mesh.points[t[1]], mesh.points[t[2]]) > 0.0);
        }
        let total: f64 = triangle_areas(&mesh.points, &mesh.triangles).iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "total area {total}");
    }

    #[test]
    fn constrained_multiple_crossing_grid_segments() {
        // Unit square corners plus midpoints forming a 2x2 crossing grid of
        // two horizontal and two vertical segments (4 interior intersections).
        let points = [
            [0.0, 0.0], // 0
            [3.0, 0.0], // 1
            [3.0, 3.0], // 2
            [0.0, 3.0], // 3
            [0.0, 1.0], // 4
            [3.0, 1.0], // 5
            [0.0, 2.0], // 6
            [3.0, 2.0], // 7
            [1.0, 0.0], // 8
            [1.0, 3.0], // 9
            [2.0, 0.0], // 10
            [2.0, 3.0], // 11
        ];
        let segments = [(4, 5), (6, 7), (8, 9), (10, 11)];
        let mesh = expect_ok(constrained_delaunay_2d_mesh(&points, &segments));

        // 12 original vertices + 4 crossing Steiner vertices at (1,1),(2,1),(1,2),(2,2).
        assert_eq!(mesh.points.len(), 16, "expected 4 Steiner intersections");

        for t in &mesh.triangles {
            assert!(orient2d(mesh.points[t[0]], mesh.points[t[1]], mesh.points[t[2]]) > 0.0);
        }
        let total: f64 = triangle_areas(&mesh.points, &mesh.triangles).iter().sum();
        assert!((total - 9.0).abs() < 1e-10, "total area {total}");

        // Verify every crossing point exists and every horizontal/vertical chain
        // of 3 sub-segments is present in the output triangulation.
        let find_pt = |x: f64, y: f64| -> usize {
            mesh.points
                .iter()
                .position(|p| (p[0] - x).abs() < 1e-9 && (p[1] - y).abs() < 1e-9)
                .unwrap_or_else(|| panic!("missing vertex ({x}, {y})"))
        };
        let p11 = find_pt(1.0, 1.0);
        let p21 = find_pt(2.0, 1.0);
        let p12 = find_pt(1.0, 2.0);
        let p22 = find_pt(2.0, 2.0);

        for (u, v) in [
            (4, p11),
            (p11, p21),
            (p21, 5),
            (6, p12),
            (p12, p22),
            (p22, 7),
            (8, p11),
            (p11, p12),
            (p12, 9),
            (10, p21),
            (p21, p22),
            (p22, 11),
        ] {
            assert!(
                has_edge(&mesh.triangles, u, v),
                "constrained sub-edge ({u}, {v}) missing"
            );
        }
    }
}
