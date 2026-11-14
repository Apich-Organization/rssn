//! Computational geometry kernels: convex hulls, Delaunay triangulation,
//! Voronoi diagrams, polygon/segment/ray/triangle predicates, k-d trees, and BVH.

/// Error returned by spatial index construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GeometryError {
    /// Input point set was empty or zero-dimensional.
    Empty,
    /// Point dimensions did not match.
    ShapeMismatch(String),
    /// Non-finite coordinates or invalid query parameters.
    Domain(String),
}

impl std::fmt::Display for GeometryError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("geometry input is empty"),
            | Self::ShapeMismatch(msg) => write!(f, "shape mismatch: {msg}"),
            | Self::Domain(msg) => write!(f, "domain error: {msg}"),
        }
    }
}

impl std::error::Error for GeometryError {}

const GEOM_EPS: f64 = 1e-12;

fn cross2(
    o: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

fn sub3(
    a: [f64; 3],
    b: [f64; 3],
) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(
    a: [f64; 3],
    b: [f64; 3],
) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(
    a: [f64; 3],
    b: [f64; 3],
) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The 2D convex hull of `points` in counter-clockwise order, via Andrew's
/// monotone chain algorithm.
///
/// Collinear interior points on hull edges are omitted.
#[must_use]
pub fn convex_hull_2d(points: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a[0].total_cmp(&b[0]).then_with(|| a[1].total_cmp(&b[1])));
    pts.dedup_by(|a, b| (a[0] - b[0]).abs() < GEOM_EPS && (a[1] - b[1]).abs() < GEOM_EPS);
    if pts.len() <= 2 {
        return pts;
    }

    let mut lower: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
    for &p in &pts {
        while lower.len() >= 2
            && cross2(lower[lower.len() - 2], lower[lower.len() - 1], p) <= GEOM_EPS
        {
            lower.pop();
        }
        lower.push(p);
    }

    let mut upper: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
    for &p in pts.iter().rev() {
        while upper.len() >= 2
            && cross2(upper[upper.len() - 2], upper[upper.len() - 1], p) <= GEOM_EPS
        {
            upper.pop();
        }
        upper.push(p);
    }

    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// The 3D convex hull of `points`, returned as outward-oriented triangular faces
/// (indices into `points`).
///
/// Returns an empty vector if fewer than 4 non-coplanar points are given.
#[must_use]
pub fn convex_hull_3d(points: &[[f64; 3]]) -> Vec<[usize; 3]> {
    let n = points.len();
    if n < 4 {
        return Vec::new();
    }

    // Find initial non-degenerate tetrahedron: p0, p1, p2, p3.
    let p0 = 0usize;
    let mut p1 = p0;
    let mut best_dist2 = 0.0;
    for (i, &pt) in points.iter().enumerate().skip(1) {
        let d = sub3(pt, points[p0]);
        let d2 = dot3(d, d);
        if d2 > best_dist2 {
            best_dist2 = d2;
            p1 = i;
        }
    }
    if best_dist2 < GEOM_EPS * GEOM_EPS {
        return Vec::new();
    }

    let e01 = sub3(points[p1], points[p0]);
    let mut p2 = p0;
    let mut best_area2 = 0.0;
    for (i, &pt) in points.iter().enumerate() {
        if i == p0 || i == p1 {
            continue;
        }
        let c = cross3(e01, sub3(pt, points[p0]));
        let a2 = dot3(c, c);
        if a2 > best_area2 {
            best_area2 = a2;
            p2 = i;
        }
    }
    if best_area2 < GEOM_EPS * GEOM_EPS {
        return Vec::new();
    }

    let n012 = cross3(e01, sub3(points[p2], points[p0]));
    let mut p3 = p0;
    let mut best_vol = 0.0;
    for (i, &pt) in points.iter().enumerate() {
        if i == p0 || i == p1 || i == p2 {
            continue;
        }
        let vol = dot3(n012, sub3(pt, points[p0])).abs();
        if vol > best_vol {
            best_vol = vol;
            p3 = i;
        }
    }
    if best_vol < GEOM_EPS * GEOM_EPS {
        return Vec::new();
    }

    let centroid = [
        0.25 * (points[p0][0] + points[p1][0] + points[p2][0] + points[p3][0]),
        0.25 * (points[p0][1] + points[p1][1] + points[p2][1] + points[p3][1]),
        0.25 * (points[p0][2] + points[p1][2] + points[p2][2] + points[p3][2]),
    ];

    let make_outward = |a: usize, b: usize, c: usize| -> [usize; 3] {
        let normal = cross3(sub3(points[b], points[a]), sub3(points[c], points[a]));
        if dot3(normal, sub3(points[a], centroid)) < 0.0 {
            [a, c, b]
        } else {
            [a, b, c]
        }
    };

    let mut faces = vec![
        make_outward(p0, p1, p2),
        make_outward(p0, p1, p3),
        make_outward(p0, p2, p3),
        make_outward(p1, p2, p3),
    ];

    for (i, &pt) in points.iter().enumerate() {
        if i == p0 || i == p1 || i == p2 || i == p3 {
            continue;
        }
        let mut visible = vec![false; faces.len()];
        let mut any_visible = false;
        for (fi, &face) in faces.iter().enumerate() {
            let normal = cross3(
                sub3(points[face[1]], points[face[0]]),
                sub3(points[face[2]], points[face[0]]),
            );
            if dot3(normal, sub3(pt, points[face[0]])) > 1e-11 {
                visible[fi] = true;
                any_visible = true;
            }
        }
        if !any_visible {
            continue;
        }

        let mut visible_edges: Vec<(usize, usize)> = Vec::new();
        let mut next_faces: Vec<[usize; 3]> = Vec::with_capacity(faces.len());
        for (fi, &face) in faces.iter().enumerate() {
            if visible[fi] {
                visible_edges.push((face[0], face[1]));
                visible_edges.push((face[1], face[2]));
                visible_edges.push((face[2], face[0]));
            } else {
                next_faces.push(face);
            }
        }

        for &(u, v) in &visible_edges {
            if !visible_edges.contains(&(v, u)) {
                next_faces.push([u, v, i]);
            }
        }
        faces = next_faces;
    }

    faces
}

fn circumcircle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> Option<([f64; 2], f64)> {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if d.abs() < 1e-15 {
        return None;
    }
    let a2 = a[0] * a[0] + a[1] * a[1];
    let b2 = b[0] * b[0] + b[1] * b[1];
    let c2 = c[0] * c[0] + c[1] * c[1];
    let ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d;
    let uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d;
    let r2 = (ux - a[0]) * (ux - a[0]) + (uy - a[1]) * (uy - a[1]);
    Some(([ux, uy], r2))
}

/// 2D Delaunay triangulation of `points`, returning triangles as `[usize; 3]`
/// vertex indices in counter-clockwise order, via the Bowyer–Watson algorithm.
#[must_use]
pub fn delaunay_2d(points: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = points.len();
    if n < 3 {
        return Vec::new();
    }

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for &p in points {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let delta = (max_x - min_x).max(max_y - min_y).max(1.0);
    let mid_x = f64::midpoint(min_x, max_x);
    let mid_y = f64::midpoint(min_y, max_y);

    let mut all_pts = points.to_vec();
    all_pts.push([mid_x - 20.0 * delta, mid_y - 10.0 * delta]);
    all_pts.push([mid_x + 20.0 * delta, mid_y - 10.0 * delta]);
    all_pts.push([mid_x, mid_y + 20.0 * delta]);

    let make_tri =
        |a: usize, b: usize, c: usize, pts: &[[f64; 2]]| -> Option<([usize; 3], [f64; 2], f64)> {
            let mut tri = [a, b, c];
            if cross2(pts[tri[0]], pts[tri[1]], pts[tri[2]]) < 0.0 {
                tri.swap(1, 2);
            }
            let (center, r2) = circumcircle(pts[tri[0]], pts[tri[1]], pts[tri[2]])?;
            Some((tri, center, r2))
        };

    let mut triangles = Vec::new();
    if let Some(super_tri) = make_tri(n, n + 1, n + 2, &all_pts) {
        triangles.push(super_tri);
    }

    for i in 0..n {
        let p = all_pts[i];
        // Skip near-duplicate points to avoid degenerate zero-area triangles.
        let duplicate = (0..i).any(|j| {
            let dx = p[0] - all_pts[j][0];
            let dy = p[1] - all_pts[j][1];
            dx * dx + dy * dy < GEOM_EPS * GEOM_EPS
        });
        if duplicate {
            continue;
        }

        let mut bad_edges: Vec<(usize, usize)> = Vec::new();
        let mut kept = Vec::with_capacity(triangles.len());
        for &(tri, center, r2) in &triangles {
            let dx = p[0] - center[0];
            let dy = p[1] - center[1];
            if dx * dx + dy * dy <= r2 * (1.0 + 1e-11) + 1e-15 {
                bad_edges.push((tri[0], tri[1]));
                bad_edges.push((tri[1], tri[2]));
                bad_edges.push((tri[2], tri[0]));
            } else {
                kept.push((tri, center, r2));
            }
        }

        for &(u, v) in &bad_edges {
            if !bad_edges.contains(&(v, u))
                && let Some(tri) = make_tri(u, v, i, &all_pts)
            {
                kept.push(tri);
            }
        }
        triangles = kept;
    }

    triangles
        .into_iter()
        .filter(|(tri, _, _)| tri[0] < n && tri[1] < n && tri[2] < n)
        .map(|(tri, _, _)| tri)
        .collect()
}

/// A 2D Voronoi diagram, dual to the Delaunay triangulation.
#[derive(Clone, Debug, PartialEq)]
pub struct VoronoiDiagram {
    /// Circumcenters of the Delaunay triangles.
    pub vertices: Vec<[f64; 2]>,
    /// Finite Voronoi edges, as pairs of indices into `vertices`.
    pub edges: Vec<(usize, usize)>,
    /// Unbounded Voronoi rays `(vertex_index, unit_direction)` for convex-hull edges.
    pub rays: Vec<(usize, [f64; 2])>,
    /// For each input site, the indices into `vertices` forming its cell in CCW order.
    pub cells: Vec<Vec<usize>>,
}

/// Computes the 2D Voronoi diagram of `points` as the dual of [`delaunay_2d`].
#[must_use]
pub fn voronoi_2d(points: &[[f64; 2]]) -> VoronoiDiagram {
    let tris = delaunay_2d(points);
    let mut vertices = Vec::with_capacity(tris.len());
    let mut cells = vec![Vec::new(); points.len()];

    for (ti, &tri) in tris.iter().enumerate() {
        let center = circumcircle(points[tri[0]], points[tri[1]], points[tri[2]])
            .map_or([0.0, 0.0], |(c, _)| c);
        vertices.push(center);
        for &site in &tri {
            cells[site].push(ti);
        }
    }

    for (site, cell) in cells.iter_mut().enumerate() {
        let s = points[site];
        cell.sort_by(|&a, &b| {
            let ang_a = (vertices[a][1] - s[1]).atan2(vertices[a][0] - s[0]);
            let ang_b = (vertices[b][1] - s[1]).atan2(vertices[b][0] - s[0]);
            ang_a.total_cmp(&ang_b)
        });
    }

    let mut edges = Vec::new();
    let mut rays = Vec::new();
    let mut half_edges: Vec<((usize, usize), usize, usize)> = Vec::with_capacity(tris.len() * 3);
    for (ti, &tri) in tris.iter().enumerate() {
        half_edges.push(((tri[0], tri[1]), ti, tri[2]));
        half_edges.push(((tri[1], tri[2]), ti, tri[0]));
        half_edges.push(((tri[2], tri[0]), ti, tri[1]));
    }

    for &((u, v), ti, opp) in &half_edges {
        if let Some(&(_, tj, _)) = half_edges.iter().find(|&&(e, _, _)| e == (v, u)) {
            if ti < tj {
                edges.push((ti, tj));
            }
        } else {
            // Outward normal to edge (u, v) pointing away from opposite vertex `opp`.
            let ex = points[v][0] - points[u][0];
            let ey = points[v][1] - points[u][1];
            let mut nx = ey;
            let mut ny = -ex;
            let to_opp_x = points[opp][0] - points[u][0];
            let to_opp_y = points[opp][1] - points[u][1];
            if nx * to_opp_x + ny * to_opp_y > 0.0 {
                nx = -nx;
                ny = -ny;
            }
            let len = nx.hypot(ny);
            if len > GEOM_EPS {
                rays.push((ti, [nx / len, ny / len]));
            }
        }
    }

    VoronoiDiagram {
        vertices,
        edges,
        rays,
        cells,
    }
}

fn point_on_segment_2d(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
) -> bool {
    if cross2(a, b, p).abs() > GEOM_EPS {
        return false;
    }
    p[0] >= a[0].min(b[0]) - GEOM_EPS
        && p[0] <= a[0].max(b[0]) + GEOM_EPS
        && p[1] >= a[1].min(b[1]) - GEOM_EPS
        && p[1] <= a[1].max(b[1]) + GEOM_EPS
}

/// Returns `true` if `point` lies inside or on the boundary of `polygon`.
#[must_use]
pub fn point_in_polygon(
    point: [f64; 2],
    polygon: &[[f64; 2]],
) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        if point_on_segment_2d(point, a, b) {
            return true;
        }
        let cond = (a[1] > point[1]) != (b[1] > point[1]);
        if cond {
            let x_intersect = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if point[0] < x_intersect {
                inside = !inside;
            }
        }
    }
    inside
}

/// Intersection point of two 2D line segments `[a1, a2]` and `[b1, b2]`, if any.
#[must_use]
pub fn segment_intersection_2d(
    a1: [f64; 2],
    a2: [f64; 2],
    b1: [f64; 2],
    b2: [f64; 2],
) -> Option<[f64; 2]> {
    let da = [a2[0] - a1[0], a2[1] - a1[1]];
    let db = [b2[0] - b1[0], b2[1] - b1[1]];
    let denom = da[0] * db[1] - da[1] * db[0];
    let diff = [b1[0] - a1[0], b1[1] - a1[1]];
    if denom.abs() < GEOM_EPS {
        // Parallel or collinear: check endpoint overlap.
        for candidate in [a1, a2] {
            if point_on_segment_2d(candidate, b1, b2) {
                return Some(candidate);
            }
        }
        for candidate in [b1, b2] {
            if point_on_segment_2d(candidate, a1, a2) {
                return Some(candidate);
            }
        }
        return None;
    }
    let t = (diff[0] * db[1] - diff[1] * db[0]) / denom;
    let u = (diff[0] * da[1] - diff[1] * da[0]) / denom;
    if (-GEOM_EPS..=1.0 + GEOM_EPS).contains(&t) && (-GEOM_EPS..=1.0 + GEOM_EPS).contains(&u) {
        let t = t.clamp(0.0, 1.0);
        Some([a1[0] + t * da[0], a1[1] + t * da[1]])
    } else {
        None
    }
}

/// Intersection of a 2D ray `origin + t * dir` (`t >= 0`) with segment `[b1, b2]`.
///
/// Returns `Some((t, point))` if they intersect.
#[must_use]
pub fn ray_segment_intersection_2d(
    origin: [f64; 2],
    dir: [f64; 2],
    b1: [f64; 2],
    b2: [f64; 2],
) -> Option<(f64, [f64; 2])> {
    let db = [b2[0] - b1[0], b2[1] - b1[1]];
    let denom = dir[0] * db[1] - dir[1] * db[0];
    if denom.abs() < GEOM_EPS {
        return None;
    }
    let diff = [b1[0] - origin[0], b1[1] - origin[1]];
    let t = (diff[0] * db[1] - diff[1] * db[0]) / denom;
    let u = (diff[0] * dir[1] - diff[1] * dir[0]) / denom;
    if t >= -GEOM_EPS && (-GEOM_EPS..=1.0 + GEOM_EPS).contains(&u) {
        let t = t.max(0.0);
        Some((t, [origin[0] + t * dir[0], origin[1] + t * dir[1]]))
    } else {
        None
    }
}

/// Intersection of a 3D ray `origin + t * dir` (`t >= 0`) with triangle `(v0, v1, v2)`,
/// via the Möller–Trumbore algorithm.
///
/// Returns `Some((t, point))` on hit.
#[must_use]
pub fn ray_triangle_intersection_3d(
    origin: [f64; 3],
    dir: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<(f64, [f64; 3])> {
    let e1 = sub3(v1, v0);
    let e2 = sub3(v2, v0);
    let pvec = cross3(dir, e2);
    let det = dot3(e1, pvec);
    if det.abs() < GEOM_EPS {
        return None;
    }
    let inv_det = 1.0 / det;
    let tvec = sub3(origin, v0);
    let u = dot3(tvec, pvec) * inv_det;
    if !(-GEOM_EPS..=1.0 + GEOM_EPS).contains(&u) {
        return None;
    }
    let qvec = cross3(tvec, e1);
    let v = dot3(dir, qvec) * inv_det;
    if v < -GEOM_EPS || u + v > 1.0 + GEOM_EPS {
        return None;
    }
    let t = dot3(e2, qvec) * inv_det;
    if t >= -GEOM_EPS {
        let t = t.max(0.0);
        Some((
            t,
            [
                origin[0] + t * dir[0],
                origin[1] + t * dir[1],
                origin[2] + t * dir[2],
            ],
        ))
    } else {
        None
    }
}

/// Intersection of a 3D line segment `[p0, p1]` with triangle `(v0, v1, v2)`.
#[must_use]
pub fn segment_triangle_intersection_3d(
    p0: [f64; 3],
    p1: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<[f64; 3]> {
    let dir = sub3(p1, p0);
    let (t, hit) = ray_triangle_intersection_3d(p0, dir, v0, v1, v2)?;
    if t <= 1.0 + GEOM_EPS {
        Some(hit)
    } else {
        None
    }
}

#[derive(Clone, Debug)]
struct KdNode {
    point_index: usize,
    axis: usize,
    left: Option<usize>,
    right: Option<usize>,
}

/// A k-d tree for exact nearest-neighbor, k-nearest, and radius queries in `R^d`.
#[derive(Clone, Debug)]
pub struct KdTree {
    points: Vec<Vec<f64>>,
    dim: usize,
    root: Option<usize>,
    nodes: Vec<KdNode>,
}

impl KdTree {
    /// Builds a k-d tree from `points`.
    ///
    /// # Errors
    ///
    /// * [`GeometryError::Empty`] if `points` is empty or 0-dimensional.
    /// * [`GeometryError::ShapeMismatch`] if points have differing dimensions.
    pub fn build(points: Vec<Vec<f64>>) -> Result<Self, GeometryError> {
        if points.is_empty() {
            return Err(GeometryError::Empty);
        }
        let dim = points[0].len();
        if dim == 0 {
            return Err(GeometryError::Empty);
        }
        for pt in &points {
            if pt.len() != dim {
                return Err(GeometryError::ShapeMismatch(format!(
                    "expected dimension {dim}, found {}",
                    pt.len()
                )));
            }
        }
        if points.iter().flatten().any(|v| !v.is_finite()) {
            return Err(GeometryError::Domain(
                "coordinates must be finite".to_string(),
            ));
        }
        let mut indices: Vec<usize> = (0..points.len()).collect();
        let mut nodes = Vec::with_capacity(points.len());
        let root = Self::build_recursive(&points, &mut indices, 0, dim, &mut nodes);
        Ok(Self {
            points,
            dim,
            root,
            nodes,
        })
    }

    fn build_recursive(
        points: &[Vec<f64>],
        indices: &mut [usize],
        depth: usize,
        dim: usize,
        nodes: &mut Vec<KdNode>,
    ) -> Option<usize> {
        if indices.is_empty() {
            return None;
        }
        let axis = depth % dim;
        let mid = indices.len() / 2;
        indices.select_nth_unstable_by(mid, |&a, &b| points[a][axis].total_cmp(&points[b][axis]));
        let point_index = indices[mid];
        let (left_slice, right_with_mid) = indices.split_at_mut(mid);
        let right_slice = &mut right_with_mid[1..];
        let left = Self::build_recursive(points, left_slice, depth + 1, dim, nodes);
        let right = Self::build_recursive(points, right_slice, depth + 1, dim, nodes);
        let node_idx = nodes.len();
        nodes.push(KdNode {
            point_index,
            axis,
            left,
            right,
        });
        Some(node_idx)
    }

    fn dist_sq(
        a: &[f64],
        b: &[f64],
    ) -> f64 {
        a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
    }

    /// Finds the nearest point to `query`, returning `(point_index, euclidean_distance)`.
    ///
    /// # Errors
    ///
    /// [`GeometryError::ShapeMismatch`] if `query.len()` does not match the tree dimension.
    pub fn nearest(
        &self,
        query: &[f64],
    ) -> Result<Option<(usize, f64)>, GeometryError> {
        let mut res = self.k_nearest(query, 1)?;
        Ok(res.pop())
    }

    /// Finds up to `k` nearest points to `query`, sorted by ascending distance.
    ///
    /// # Errors
    ///
    /// [`GeometryError::ShapeMismatch`] if `query.len()` does not match the tree dimension.
    pub fn k_nearest(
        &self,
        query: &[f64],
        k: usize,
    ) -> Result<Vec<(usize, f64)>, GeometryError> {
        if query.iter().any(|v| !v.is_finite()) {
            return Err(GeometryError::Domain(
                "query coordinates must be finite".to_string(),
            ));
        }
        if query.len() != self.dim {
            return Err(GeometryError::ShapeMismatch(format!(
                "query dimension {} does not match tree dimension {}",
                query.len(),
                self.dim
            )));
        }
        if k == 0 {
            return Ok(Vec::new());
        }
        let mut best: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        if let Some(root) = self.root {
            self.k_nearest_rec(root, query, k, &mut best);
        }
        for item in &mut best {
            item.1 = item.1.sqrt();
        }
        Ok(best)
    }

    fn k_nearest_rec(
        &self,
        node_idx: usize,
        query: &[f64],
        k: usize,
        best: &mut Vec<(usize, f64)>,
    ) {
        let node = &self.nodes[node_idx];
        let pt = &self.points[node.point_index];
        let d2 = Self::dist_sq(query, pt);
        if best.len() < k || d2 < best.last().map_or(f64::INFINITY, |x| x.1) {
            let pos = best
                .binary_search_by(|x| x.1.total_cmp(&d2))
                .unwrap_or_else(|e| e);
            best.insert(pos, (node.point_index, d2));
            if best.len() > k {
                best.pop();
            }
        }

        let diff = query[node.axis] - pt[node.axis];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(near_idx) = near {
            self.k_nearest_rec(near_idx, query, k, best);
        }
        let worst_d2 = if best.len() < k {
            f64::INFINITY
        } else {
            best.last().map_or(f64::INFINITY, |x| x.1)
        };
        if diff * diff <= worst_d2
            && let Some(far_idx) = far
        {
            self.k_nearest_rec(far_idx, query, k, best);
        }
    }

    /// Finds all points within Euclidean `radius` of `query`, sorted by ascending distance.
    ///
    /// # Errors
    ///
    /// [`GeometryError::ShapeMismatch`] if `query.len()` does not match the tree dimension.
    pub fn radius_search(
        &self,
        query: &[f64],
        radius: f64,
    ) -> Result<Vec<(usize, f64)>, GeometryError> {
        if query.iter().any(|v| !v.is_finite()) {
            return Err(GeometryError::Domain(
                "query coordinates must be finite".to_string(),
            ));
        }
        if query.len() != self.dim {
            return Err(GeometryError::ShapeMismatch(format!(
                "query dimension {} does not match tree dimension {}",
                query.len(),
                self.dim
            )));
        }
        if radius.is_nan() {
            return Err(GeometryError::Domain("radius must not be NaN".to_string()));
        }
        if radius < 0.0 {
            return Ok(Vec::new());
        }
        let r2 = radius * radius;
        let mut found = Vec::new();
        if let Some(root) = self.root {
            self.radius_rec(root, query, r2, &mut found);
        }
        found.sort_by(|a, b| a.1.total_cmp(&b.1));
        for item in &mut found {
            item.1 = item.1.sqrt();
        }
        Ok(found)
    }

    fn radius_rec(
        &self,
        node_idx: usize,
        query: &[f64],
        r2: f64,
        found: &mut Vec<(usize, f64)>,
    ) {
        let node = &self.nodes[node_idx];
        let pt = &self.points[node.point_index];
        let d2 = Self::dist_sq(query, pt);
        if d2 <= r2 {
            found.push((node.point_index, d2));
        }
        let diff = query[node.axis] - pt[node.axis];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(near_idx) = near {
            self.radius_rec(near_idx, query, r2, found);
        }
        if diff * diff <= r2
            && let Some(far_idx) = far
        {
            self.radius_rec(far_idx, query, r2, found);
        }
    }
}

/// A 3D axis-aligned bounding box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb {
    /// Minimum corner `[x, y, z]`.
    pub min: [f64; 3],
    /// Maximum corner `[x, y, z]`.
    pub max: [f64; 3],
}

impl Aabb {
    fn union(
        self,
        other: Self,
    ) -> Self {
        Self {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    fn centroid(self) -> [f64; 3] {
        [
            f64::midpoint(self.min[0], self.max[0]),
            f64::midpoint(self.min[1], self.max[1]),
            f64::midpoint(self.min[2], self.max[2]),
        ]
    }

    fn dist_sq_to_point(
        self,
        p: [f64; 3],
    ) -> f64 {
        let mut d2 = 0.0;
        for axis in 0..3 {
            if p[axis] < self.min[axis] {
                let d = self.min[axis] - p[axis];
                d2 += d * d;
            } else if p[axis] > self.max[axis] {
                let d = p[axis] - self.max[axis];
                d2 += d * d;
            }
        }
        d2
    }

    fn ray_hits(
        self,
        origin: [f64; 3],
        dir: [f64; 3],
    ) -> bool {
        let mut tmin = 0.0f64;
        let mut tmax = f64::INFINITY;
        for axis in 0..3 {
            if dir[axis].abs() < GEOM_EPS {
                if origin[axis] < self.min[axis] || origin[axis] > self.max[axis] {
                    return false;
                }
            } else {
                let inv = 1.0 / dir[axis];
                let mut t1 = (self.min[axis] - origin[axis]) * inv;
                let mut t2 = (self.max[axis] - origin[axis]) * inv;
                if t1 > t2 {
                    std::mem::swap(&mut t1, &mut t2);
                }
                tmin = tmin.max(t1);
                tmax = tmax.min(t2);
                if tmin > tmax {
                    return false;
                }
            }
        }
        true
    }
}

#[derive(Clone, Debug)]
enum BvhNode {
    Leaf {
        bounds: Aabb,
        item_index: usize,
    },
    Internal {
        bounds: Aabb,
        left: usize,
        right: usize,
    },
}

impl BvhNode {
    fn bounds(&self) -> Aabb {
        match *self {
            | Self::Leaf { bounds, .. } | Self::Internal { bounds, .. } => bounds,
        }
    }
}

/// A 3D Bounding Volume Hierarchy (BVH) over [`Aabb`] primitives for fast
/// nearest-box and ray-intersection queries.
#[derive(Clone, Debug)]
pub struct Bvh {
    nodes: Vec<BvhNode>,
    root: Option<usize>,
}

impl Bvh {
    /// Builds a BVH over `aabbs`.
    #[must_use]
    pub fn build(aabbs: &[Aabb]) -> Self {
        if aabbs.is_empty() {
            return Self {
                nodes: Vec::new(),
                root: None,
            };
        }
        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        let mut nodes = Vec::with_capacity(2 * aabbs.len());
        let root = Some(Self::build_rec(aabbs, &mut indices, &mut nodes));
        Self { nodes, root }
    }

    fn build_rec(
        aabbs: &[Aabb],
        indices: &mut [usize],
        nodes: &mut Vec<BvhNode>,
    ) -> usize {
        if indices.len() == 1 {
            let idx = nodes.len();
            nodes.push(BvhNode::Leaf {
                bounds: aabbs[indices[0]],
                item_index: indices[0],
            });
            return idx;
        }
        let mut bounds = aabbs[indices[0]];
        for &i in &indices[1..] {
            bounds = bounds.union(aabbs[i]);
        }
        let extent = [
            bounds.max[0] - bounds.min[0],
            bounds.max[1] - bounds.min[1],
            bounds.max[2] - bounds.min[2],
        ];
        let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
            0
        } else if extent[1] >= extent[2] {
            1
        } else {
            2
        };
        let mid = indices.len() / 2;
        indices.select_nth_unstable_by(mid, |&a, &b| {
            aabbs[a].centroid()[axis].total_cmp(&aabbs[b].centroid()[axis])
        });
        let (left_slice, right_slice) = indices.split_at_mut(mid);
        let left = Self::build_rec(aabbs, left_slice, nodes);
        let right = Self::build_rec(aabbs, right_slice, nodes);
        let idx = nodes.len();
        nodes.push(BvhNode::Internal { bounds, left, right });
        idx
    }

    /// Returns `(aabb_index, distance)` for the [`Aabb`] closest to `point`.
    #[must_use]
    pub fn nearest_aabb(
        &self,
        point: [f64; 3],
    ) -> Option<(usize, f64)> {
        let root = self.root?;
        let mut best: Option<(usize, f64)> = None;
        self.nearest_rec(root, point, &mut best);
        best.map(|(idx, d2)| (idx, d2.sqrt()))
    }

    fn nearest_rec(
        &self,
        node_idx: usize,
        point: [f64; 3],
        best: &mut Option<(usize, f64)>,
    ) {
        let best_d2 = best.map_or(f64::INFINITY, |x| x.1);
        match self.nodes[node_idx] {
            | BvhNode::Leaf { bounds, item_index } => {
                let d2 = bounds.dist_sq_to_point(point);
                if d2 < best_d2 {
                    *best = Some((item_index, d2));
                }
            },
            | BvhNode::Internal { left, right, .. } => {
                let dl = self.nodes[left].bounds().dist_sq_to_point(point);
                let dr = self.nodes[right].bounds().dist_sq_to_point(point);
                let (first, d_first, second, d_second) = if dl <= dr {
                    (left, dl, right, dr)
                } else {
                    (right, dr, left, dl)
                };
                if d_first < best.map_or(f64::INFINITY, |x| x.1) {
                    self.nearest_rec(first, point, best);
                }
                if d_second < best.map_or(f64::INFINITY, |x| x.1) {
                    self.nearest_rec(second, point, best);
                }
            },
        }
    }

    /// Returns the indices of all [`Aabb`]s intersected by the ray `origin + t * dir` (`t >= 0`).
    #[must_use]
    pub fn ray_intersect(
        &self,
        origin: [f64; 3],
        dir: [f64; 3],
    ) -> Vec<usize> {
        let mut hits = Vec::new();
        let Some(root) = self.root else {
            return hits;
        };
        let mut stack = vec![root];
        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            if !node.bounds().ray_hits(origin, dir) {
                continue;
            }
            match *node {
                | BvhNode::Leaf { item_index, .. } => hits.push(item_index),
                | BvhNode::Internal { left, right, .. } => {
                    stack.push(left);
                    stack.push(right);
                },
            }
        }
        hits
    }
}


/// The area of a polygon given its vertices, by the shoelace formula.
///
/// `points` is `[x0, y0, x1, y1, …]`; the result is non-negative and does not
/// depend on the winding direction.
#[must_use]
pub fn polygon_area(points: &[f64]) -> f64 {
    let count = points.len() / 2;
    if count < 3 {
        return 0.0;
    }
    let mut total = 0.0;
    for index in 0..count {
        let next = (index + 1) % count;
        total +=
            points[index * 2] * points[next * 2 + 1] - points[next * 2] * points[index * 2 + 1];
    }
    (total / 2.0).abs()
}

/// The Euclidean distance between two points of the same dimension.
///
/// # Errors
///
/// [`GeometryError::ShapeMismatch`] when the dimensions differ.
pub fn distance(
    left: &[f64],
    right: &[f64],
) -> Result<f64, GeometryError> {
    if left.len() != right.len() {
        return Err(GeometryError::ShapeMismatch(format!(
            "distance needs two points of one dimension, found {} and {}",
            left.len(),
            right.len()
        )));
    }
    Ok(left
        .iter()
        .zip(right)
        .fold(0.0_f64, |norm, (a, b)| norm.hypot(a - b)))
}

/// The non-negative volume of a tetrahedron, from four vertices in `R^3`.
#[must_use]
pub fn tetrahedron_volume(points: &[f64]) -> f64 {
    if points.len() < 12 {
        return 0.0;
    }
    let a = [
        points[3] - points[0],
        points[4] - points[1],
        points[5] - points[2],
    ];
    let b = [
        points[6] - points[0],
        points[7] - points[1],
        points[8] - points[2],
    ];
    let c = [
        points[9] - points[0],
        points[10] - points[1],
        points[11] - points[2],
    ];
    let cross = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]).abs() / 6.0
}
