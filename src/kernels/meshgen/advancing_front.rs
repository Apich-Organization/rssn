//! Advancing-front triangle generation from a closed boundary loop.
//!
//! The boundary is first resampled at the requested spacing, then the front is
//! advanced inwards one edge at a time: for the shortest usable front edge a
//! new vertex is placed on the interior side at the height of an equilateral
//! triangle and accepted only when the resulting triangle stays inside the
//! domain, does not cross the front and does not swallow an existing vertex.
//! Edges that cannot advance are retired and the remaining region is closed by
//! ear clipping, which guarantees that the mesh covers the whole domain.

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use super::util::TriMesh2;
use super::util::distance2;
use super::util::ear_clip;
use super::util::orient2d;
use super::util::point_in_polygon;
use super::util::point_in_triangle;
use super::util::point_segment_distance;
use super::util::polygon_signed_area;
use super::util::segments_cross;
use super::util::triangle_area;

/// The height of an equilateral triangle of unit side.
const EQUILATERAL_HEIGHT: f64 = 0.866_025_403_784_438_6;

/// Meshes the interior of the closed loop `boundary` by front advancement.
///
/// `spacing` is the target edge length; boundary edges are resampled to it and
/// new interior vertices are placed at the height of an equilateral triangle
/// over the front. The returned mesh is counter-clockwise and covers the
/// polygon.
///
/// # Errors
///
/// Returns `Err` when the boundary has fewer than three points, contains
/// non-finite coordinates or repeated consecutive points, when the enclosed
/// area is zero, or when `spacing` is not positive and finite.
pub fn advancing_front_2d(
    boundary: &[[f64; 2]],
    spacing: f64,
) -> Result<TriMesh2, String> {
    advancing_front_bounded(boundary, spacing, usize::MAX)
}

pub(crate) fn advancing_front_bounded(
    boundary: &[[f64; 2]],
    spacing: f64,
    vertex_budget: usize,
) -> Result<TriMesh2, String> {
    if boundary.len() > vertex_budget {
        return Err("advancing-front vertex budget exceeded".into());
    }
    if boundary.len() < 3 {
        return Err("at least three boundary points are required".to_string());
    }
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err("spacing must be positive and finite".to_string());
    }

    let mut loop_points: Vec<[f64; 2]> = Vec::with_capacity(boundary.len());
    for p in boundary {
        if !p[0].is_finite() || !p[1].is_finite() {
            return Err("boundary coordinates must be finite".to_string());
        }
        if let Some(&last) = loop_points.last() {
            if distance2(last, *p) <= 0.0 {
                return Err("consecutive boundary points must differ".to_string());
            }
        }
        loop_points.push(*p);
    }
    if loop_points.len() >= 2
        && distance2(loop_points[0], loop_points[loop_points.len() - 1]) <= 0.0
    {
        loop_points.pop();
    }

    let area = polygon_signed_area(&loop_points);
    if area.abs() <= 1e-14 {
        return Err("boundary encloses no area".to_string());
    }
    if area < 0.0 {
        loop_points.reverse();
    }

    let mut points = resample(&loop_points, spacing, vertex_budget)?;
    if points.len() < 3 {
        return Err("boundary is too short for the requested spacing".to_string());
    }

    let count = points.len();
    let mut front: Vec<(usize, usize)> = (0..count).map(|i| (i, (i + 1) % count)).collect();
    let mut retired: Vec<bool> = vec![false; front.len()];
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let mut queue: VecDeque<usize> = (0..front.len()).collect();
    let mut stuck: Vec<u32> = vec![0; front.len()];

    let max_steps = front.len().saturating_mul(100).saturating_add(10_000);
    let mut steps = 0_usize;

    while let Some(edge) = queue.pop_front() {
        steps += 1;
        if steps.is_multiple_of(64) {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        }
        if steps > max_steps {
            break;
        }
        if retired[edge] || edge >= front.len() {
            continue;
        }

        let (a, b) = front[edge];
        let apex = match try_advance(&points, &front, &triangles, &loop_points, a, b, spacing) {
            | Some(c) => c,
            | None => {
                stuck[edge] += 1;
                if stuck[edge] < 3 {
                    queue.push_back(edge);
                } else {
                    retired[edge] = true;
                }
                continue;
            },
        };

        let c = points.len();
        if c >= vertex_budget {
            return Err("advancing-front vertex budget exceeded".into());
        }
        points.push(apex);
        triangles.push([a, b, c]);

        front[edge] = (a, c);
        retired[edge] = false;
        stuck[edge] = 0;
        queue.push_back(edge);

        front.push((c, b));
        retired.push(false);
        stuck.push(0);
        queue.push_back(front.len() - 1);
    }

    close_remaining(&points, &front, &mut triangles);

    Ok(TriMesh2 { points, triangles })
}

/// Resamples the closed loop so that every edge is close to `spacing` long.
fn resample(
    loop_points: &[[f64; 2]],
    spacing: f64,
    vertex_budget: usize,
) -> Result<Vec<[f64; 2]>, String> {
    let n = loop_points.len();
    let mut count = Vec::with_capacity(n);
    let mut total = 0_usize;
    let edges = loop_points
        .iter()
        .copied()
        .zip(loop_points.iter().copied().cycle().skip(1));
    for (p, q) in edges.clone() {
        let length = (q[0] - p[0]).hypot(q[1] - p[1]);
        let intervals = (length / spacing).round().max(1.0);
        if !intervals.is_finite() || intervals >= usize::MAX as f64 {
            return Err("advancing-front vertex budget exceeded".into());
        }
        let intervals = intervals as usize;
        total = total
            .checked_add(intervals)
            .filter(|&count| count <= vertex_budget)
            .ok_or("advancing-front vertex budget exceeded")?;
        count.push(intervals);
    }
    let mut out = Vec::with_capacity(total);
    for ((p, q), intervals) in edges.zip(count) {
        for j in 0..intervals {
            let t = j as f64 / intervals as f64;
            out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
        }
    }
    Ok(out)
}

/// Tries to place a new vertex over the front edge `a`-`b`.
///
/// Candidates are tried from the full equilateral height down to a flatter
/// triangle; the first candidate that passes every validity check wins.
fn try_advance(
    points: &[[f64; 2]],
    front: &[(usize, usize)],
    triangles: &[[usize; 3]],
    boundary: &[[f64; 2]],
    a: usize,
    b: usize,
    spacing: f64,
) -> Option<[f64; 2]> {
    let pa = points[a];
    let pb = points[b];
    let dx = pb[0] - pa[0];
    let dy = pb[1] - pa[1];
    let length = dx.hypot(dy);
    if length <= 0.0 {
        return None;
    }

    let nx = -dy / length;
    let ny = dx / length;
    let mid = [f64::midpoint(pa[0], pb[0]), f64::midpoint(pa[1], pb[1])];
    let base = EQUILATERAL_HEIGHT * spacing.max(length * 0.5);

    for factor in [1.0, 0.75, 0.55, 1.25, 0.35] {
        let h = base * factor;
        let candidate = [mid[0] + nx * h, mid[1] + ny * h];
        if valid_apex(points, front, triangles, boundary, a, b, candidate, spacing) {
            return Some(candidate);
        }
    }
    None
}

/// Whether `candidate` can be used as the apex of the triangle over `a`-`b`.
fn valid_apex(
    points: &[[f64; 2]],
    front: &[(usize, usize)],
    triangles: &[[usize; 3]],
    boundary: &[[f64; 2]],
    a: usize,
    b: usize,
    candidate: [f64; 2],
    spacing: f64,
) -> bool {
    if !point_in_polygon(candidate, boundary) {
        return false;
    }

    let min_distance = 0.35 * spacing;
    let min_distance2 = min_distance * min_distance;
    for &p in points {
        if distance2(p, candidate) <= min_distance2 {
            return false;
        }
    }

    let pa = points[a];
    let pb = points[b];
    for t in triangles {
        if point_in_triangle(candidate, points[t[0]], points[t[1]], points[t[2]], 1e-14) {
            return false;
        }
    }

    // No existing vertex may end up strictly inside the new triangle.
    for (index, &p) in points.iter().enumerate() {
        if index == a || index == b {
            continue;
        }
        if point_in_triangle(p, pa, pb, candidate, 1e-14) {
            return false;
        }
    }

    // The two new sides must not cross any front edge that is not incident.
    for &(u, v) in front {
        if u == a || v == a || u == b || v == b {
            continue;
        }
        let pu = points[u];
        let pv = points[v];
        if segments_cross(pa, candidate, pu, pv) || segments_cross(pb, candidate, pu, pv) {
            return false;
        }
        let margin = 0.25 * spacing;
        let bbox_dist = |i| {
            let lo = f64::min(pu[i], pv[i]);
            let hi = f64::max(pu[i], pv[i]);
            f64::max(lo - candidate[i], candidate[i] - hi)
        };
        if bbox_dist(0) <= margin
            && bbox_dist(1) <= margin
            && point_segment_distance(candidate, pu, pv) <= margin
        {
            return false;
        }
    }

    triangle_area(pa, pb, candidate) > 0.0 && orient2d(pa, pb, candidate) > 0.0
}

/// Ear-clips the loops that the front left unfilled.
fn close_remaining(
    points: &[[f64; 2]],
    front: &[(usize, usize)],
    triangles: &mut Vec<[usize; 3]>,
) {
    let mut next: HashMap<usize, usize> = HashMap::new();
    for &(u, v) in front {
        next.insert(u, v);
    }

    let mut visited: HashSet<usize> = HashSet::new();
    for &(start, _) in front {
        if !visited.insert(start) {
            continue;
        }

        let mut cycle = vec![start];
        let mut current = start;
        let mut guard = 0;
        while let Some(&following) = next.get(&current) {
            if following == start || !visited.insert(following) {
                break;
            }
            cycle.push(following);
            current = following;
            guard += 1;
            if guard > next.len() + 2 {
                break;
            }
        }

        if cycle.len() >= 3 {
            triangles.extend(ear_clip(points, &cycle));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::advancing_front_2d;
    use crate::kernels::meshgen::util::orient2d;
    use crate::kernels::meshgen::util::triangle_area;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    #[test]
    fn unit_square_is_covered_by_the_front() {
        let boundary = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let mesh = expect_ok(advancing_front_2d(&boundary, 0.25));

        assert_ne!(mesh.triangles.len(), 0);

        let mut total = 0.0;
        for t in &mesh.triangles {
            let a = mesh.points[t[0]];
            let b = mesh.points[t[1]];
            let c = mesh.points[t[2]];
            assert!(orient2d(a, b, c) > 0.0, "inverted triangle {t:?}");
            total += triangle_area(a, b, c);
        }

        assert!((total - 1.0).abs() <= 0.05, "covered area {total}");
    }

    #[test]
    fn clockwise_boundary_is_accepted() {
        let boundary = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        let mesh = expect_ok(advancing_front_2d(&boundary, 0.2));
        let total: f64 = mesh
            .triangles
            .iter()
            .map(|t| triangle_area(mesh.points[t[0]], mesh.points[t[1]], mesh.points[t[2]]))
            .sum();
        assert!((total - 1.0).abs() <= 0.05, "covered area {total}");
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(advancing_front_2d(&[[0.0, 0.0], [1.0, 0.0]], 0.5).is_err());
        assert!(advancing_front_2d(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], 0.0).is_err());
        assert!(advancing_front_2d(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], f64::NAN).is_err());
        assert!(advancing_front_2d(&[[0.0, 0.0], [0.0, 0.0], [1.0, 1.0]], 0.5).is_err());
        assert!(advancing_front_2d(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]], 0.5).is_err());
    }
}
