//! Computational geometry kernels: convex hulls, Delaunay triangulation,
//! Voronoi diagrams, polygon/segment/ray/triangle predicates, k-d trees, and BVH.

mod kd_tree;
pub use kd_tree::KdTree;

/// Error from spatial index construction, spatial queries or distance evaluation.
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

    fn distance_to_point(
        self,
        p: [f64; 3],
    ) -> f64 {
        let mut distance = 0.0_f64;
        for axis in 0..3 {
            let offset = if p[axis] < self.min[axis] {
                self.min[axis] - p[axis]
            } else if p[axis] > self.max[axis] {
                p[axis] - self.max[axis]
            } else {
                0.0
            };
            distance = distance.hypot(offset);
        }
        distance
    }

    fn ray_hits<const CHECK_RANGE: bool>(
        self,
        origin: [f64; 3],
        dir: [f64; 3],
    ) -> bool {
        let mut tmin = 0.0f64;
        let mut tmax = f64::INFINITY;
        for axis in 0..3 {
            if crate::kernels::exact_eq(dir[axis], 0.0) {
                if origin[axis] < self.min[axis] || origin[axis] > self.max[axis] {
                    return false;
                }
            } else {
                let lo = self.min[axis] - origin[axis];
                let hi = self.max[axis] - origin[axis];
                let mut t1 = lo / dir[axis];
                let mut t2 = hi / dir[axis];
                if CHECK_RANGE
                    && (slab_limit_needs_exact(t1, lo) || slab_limit_needs_exact(t2, hi))
                    && let Some(hit) = self.ray_hits_exact(origin, dir)
                {
                    return hit;
                }
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

    #[cold]
    fn ray_hits_exact(
        self,
        origin: [f64; 3],
        dir: [f64; 3],
    ) -> Option<bool> {
        use num_rational::BigRational;

        let mut enter = BigRational::default();
        let mut leave: Option<BigRational> = None;
        for axis in 0..3 {
            if dir[axis] == 0.0 {
                if origin[axis] < self.min[axis] || origin[axis] > self.max[axis] {
                    return Some(false);
                }
                continue;
            }
            let origin = BigRational::from_float(origin[axis])?;
            let direction = BigRational::from_float(dir[axis])?;
            let mut first = (BigRational::from_float(self.min[axis])? - &origin) / &direction;
            let mut last = (BigRational::from_float(self.max[axis])? - origin) / direction;
            if first > last {
                std::mem::swap(&mut first, &mut last);
            }
            enter = enter.max(first);
            leave = Some(match leave {
                | Some(previous) => previous.min(last),
                | None => last,
            });
            if leave.as_ref().is_some_and(|end| enter > *end) {
                return Some(false);
            }
        }
        Some(true)
    }

    #[cold]
    fn exact_squared_distance(
        self,
        point: [f64; 3],
        exact_point: &[num_bigint::BigInt; 3],
    ) -> ExactSquaredDistance {
        use num_traits::Float;

        let mut squared = num_bigint::BigInt::default();
        for axis in 0..3 {
            let boundary = if point[axis] < self.min[axis] {
                self.min[axis]
            } else if point[axis] > self.max[axis] {
                self.max[axis]
            } else {
                continue;
            };
            if !boundary.is_finite() {
                return ExactSquaredDistance::Unbounded;
            }
            let delta = binary_lattice_coordinate(boundary.integer_decode(), F64_LATTICE_EXPONENT)
                - &exact_point[axis];
            squared += &delta * &delta;
        }
        ExactSquaredDistance::Bounded(squared)
    }
}

fn slab_limit_needs_exact(
    value: f64,
    delta: f64,
) -> bool {
    let magnitude = value.abs();
    magnitude.is_infinite() || (magnitude < f64::MIN_POSITIVE && delta != 0.0)
}

// Float::integer_decode uses exponent -1075 for subnormal f64 values.
const F64_LATTICE_EXPONENT: i32 = -1075;

#[derive(PartialEq, PartialOrd)]
enum ExactSquaredDistance {
    Bounded(num_bigint::BigInt),
    Unbounded,
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
    ordinary_coordinates: bool,
}

impl Bvh {
    /// Builds a BVH over `aabbs`.
    #[must_use]
    pub fn build(aabbs: &[Aabb]) -> Self {
        if aabbs.is_empty() {
            return Self {
                nodes: Vec::new(),
                root: None,
                ordinary_coordinates: true,
            };
        }
        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        let mut nodes = Vec::with_capacity(2 * aabbs.len());
        let root = Some(Self::build_rec(aabbs, &mut indices, &mut nodes));
        let ordinary_coordinates = aabbs.iter().all(|bounds| {
            bounds
                .min
                .into_iter()
                .chain(bounds.max)
                .all(kd_tree::ordinary_coordinate)
        });
        Self {
            nodes,
            root,
            ordinary_coordinates,
        }
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
    /// If its distance overflows `f64`, the index remains the nearest and
    /// the distance is positive infinity.
    #[must_use]
    pub fn nearest_aabb(
        &self,
        point: [f64; 3],
    ) -> Option<(usize, f64)> {
        let root = self.root?;
        let mut best: Option<(usize, f64)> = None;
        self.nearest_rec(root, &mut best, &|bounds| bounds.distance_to_point(point));
        if best.is_some_and(|(_, distance)| distance.is_infinite())
            && point.iter().all(|coordinate| coordinate.is_finite())
        {
            return self.nearest_overflowing(root, point);
        }
        best
    }

    #[cold]
    fn nearest_overflowing(
        &self,
        root: usize,
        point: [f64; 3],
    ) -> Option<(usize, f64)> {
        use num_traits::Float;

        let exact_point = point.map(|coordinate| {
            binary_lattice_coordinate(coordinate.integer_decode(), F64_LATTICE_EXPONENT)
        });
        let mut best = None;
        self.nearest_rec(root, &mut best, &|bounds| {
            bounds.exact_squared_distance(point, &exact_point)
        });
        best.map(|(index, _)| (index, f64::INFINITY))
    }

    fn nearest_rec<D: PartialOrd>(
        &self,
        node_idx: usize,
        best: &mut Option<(usize, D)>,
        distance: &impl Fn(Aabb) -> D,
    ) {
        match self.nodes[node_idx] {
            | BvhNode::Leaf { bounds, item_index } => {
                let key = distance(bounds);
                if best.as_ref().is_none_or(|(_, limit)| key < *limit) {
                    *best = Some((item_index, key));
                }
            },
            | BvhNode::Internal { left, right, .. } => {
                let dl = distance(self.nodes[left].bounds());
                let dr = distance(self.nodes[right].bounds());
                let (first, d_first, second, d_second) = if dl <= dr {
                    (left, dl, right, dr)
                } else {
                    (right, dr, left, dl)
                };
                if best.as_ref().is_none_or(|(_, limit)| d_first < *limit) {
                    self.nearest_rec(first, best, distance);
                }
                if best.as_ref().is_none_or(|(_, limit)| d_second < *limit) {
                    self.nearest_rec(second, best, distance);
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
        let Some(root) = self.root else {
            return Vec::new();
        };
        // If all nonzero operands fall in 2^[-450, 450], the smallest
        // nonzero subtraction is still at least 2^-502. Dividing by any
        // of the directions stays inside the normal f64 range.
        if self.ordinary_coordinates
            && origin
                .into_iter()
                .chain(dir)
                .all(kd_tree::ordinary_coordinate)
        {
            self.ray_intersect_from::<false>(root, origin, dir)
        } else {
            self.ray_intersect_from::<true>(root, origin, dir)
        }
    }

    fn ray_intersect_from<const CHECK_RANGE: bool>(
        &self,
        root: usize,
        origin: [f64; 3],
        dir: [f64; 3],
    ) -> Vec<usize> {
        let mut hits = Vec::new();
        let mut stack = vec![root];
        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            if !node.bounds().ray_hits::<CHECK_RANGE>(origin, dir) {
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

#[cfg(test)]
mod bvh_tests {
    use super::Aabb;
    use super::Bvh;

    #[test]
    fn bvh_ray_is_invariant_when_parameters_leave_f64_range() {
        let boxes = [
            Aabb {
                min: [1.0, 2.0, 0.0],
                max: [1.25, 2.25, 0.0],
            },
            Aabb {
                min: [2.0, 2.0, 0.0],
                max: [3.0, 3.0, 0.0],
            },
        ];
        let tree = Bvh::build(&boxes);
        for direction in [1.0, 1.0e-310, f64::MAX] {
            assert_eq!(
                tree.ray_intersect([0.0; 3], [direction, direction, 0.0]),
                [1]
            );
        }
        let tiny = f64::from_bits(1);
        let scaled_boxes = boxes.map(|aabb| {
            Aabb {
                min: aabb.min.map(|v| v * tiny),
                max: aabb.max.map(|v| v * tiny),
            }
        });
        let scaled = Bvh::build(&scaled_boxes);
        assert_eq!(scaled.ray_intersect([0.0; 3], [1.0, 1.0, 0.0]), [1]);
        assert_eq!(
            scaled.ray_intersect([0.0; 3], [f64::MAX, f64::MAX, 0.0]),
            [1]
        );

        // The subtraction can overflow even when the resulting parameter is small.
        let m = f64::MAX;
        let across = Bvh::build(&[Aabb {
            min: [m, 2.0, 0.0],
            max: [m, 2.5, 0.0],
        }]);
        assert_eq!(across.ray_intersect([-m, 0.0, 0.0], [m, 1.0, 0.0]), [0]);
    }

    #[test]
    fn bvh_nearest_orders_distances_that_overflow_the_result() {
        let m = f64::MAX;
        let boxes = [-m * 0.75, -m * 0.25].map(|x| {
            Aabb {
                min: [x, 0.0, 0.0],
                max: [x, 0.0, 0.0],
            }
        });
        let tree = Bvh::build(&boxes);
        assert_eq!(tree.nearest_aabb([m * 0.125, 0.0, 0.0]).unwrap().0, 1);
        assert_eq!(tree.nearest_aabb([m, 0.0, 0.0]), Some((1, f64::INFINITY)));

        // A small second coordinate can break a tie in the overflowing axis.
        let tied = [
            Aabb {
                min: [-m, 0.0, 0.0],
                max: [-m, 0.0, 0.0],
            },
            Aabb {
                min: [-m, 1.0, 0.0],
                max: [-m, 1.0, 0.0],
            },
        ];
        assert_eq!(
            Bvh::build(&tied).nearest_aabb([m, 10.0, 0.0]),
            Some((1, f64::INFINITY))
        );
        assert_eq!(Bvh::build(&[]).nearest_aabb([m, 0.0, 0.0]), None);
    }
}

/// The area of a polygon given its vertices, by the shoelace formula.
///
/// `points` is `[x0, y0, x1, y1, …]`; the result is non-negative and does not
/// depend on the winding direction. Fewer than three vertices give zero.
/// Ill-conditioned and extreme-scale inputs use the exact represented coordinates.
#[must_use]
pub fn polygon_area(points: &[f64]) -> f64 {
    const MIN_RATIO: f64 = f64::from_bits((1023 - 20) << 52);
    let (vertices, _) = points.as_chunks::<2>();
    if vertices.len() < 3 {
        return 0.0;
    }
    let origin = vertices[0];
    let determinant = |a: &[f64; 2], b: &[f64; 2]| {
        let positive = (a[0] - origin[0]) * (b[1] - origin[1]);
        let negative = (a[1] - origin[1]) * (b[0] - origin[0]);
        (positive - negative, positive.abs() + negative.abs())
    };
    let mut total = 0.0;
    let mut compensation = 0.0;
    let mut magnitude_bound = 0.0;
    let mut accumulate = |term: f64, bound: f64| {
        let term = term - compensation;
        let next = total + term;
        compensation = (next - total) - term;
        total = next;
        magnitude_bound += bound;
    };

    // The fan excludes the two zero edges incident on its origin. Within
    // each block the determinants are independent; only block totals need
    // a compensated accumulation across the polygon.
    let (from, from_tail) = vertices[1..vertices.len() - 1].as_chunks::<8>();
    let (to, to_tail) = vertices[2..].as_chunks::<8>();
    for (from, to) in from.iter().zip(to) {
        let terms = std::array::from_fn(|i| determinant(&from[i], &to[i]));
        let pairwise_sum =
            |[a, b, c, d, e, f, g, h]: [f64; 8]| ((a + b) + (c + d)) + ((e + f) + (g + h));
        accumulate(
            pairwise_sum(terms.map(|term| term.0)),
            pairwise_sum(terms.map(|term| term.1)),
        );
    }
    for (from, to) in from_tail.iter().zip(to_tail) {
        let (value, bound) = determinant(from, to);
        accumulate(value, bound);
    }
    let magnitude = total.abs();
    let underflow_floor = (vertices.len() as f64) * f64::MIN_POSITIVE;
    if magnitude.is_finite()
        && magnitude > MIN_RATIO * magnitude_bound
        && magnitude > underflow_floor
    {
        return magnitude / 2.0;
    }
    if magnitude_bound == 0.0
        && (0..2).any(|axis| {
            vertices
                .iter()
                .all(|point| crate::kernels::exact_eq(point[axis], origin[axis]))
        })
    {
        return 0.0;
    }
    exact_polygon_area(vertices).unwrap_or(magnitude / 2.0)
}

#[cold]
fn exact_polygon_area(vertices: &[[f64; 2]]) -> Option<f64> {
    use num_bigint::BigInt;
    use num_traits::Float;

    if vertices.iter().flatten().any(|point| !point.is_finite()) {
        return None;
    }
    let exponents = [0, 1].map(|axis| {
        binary_lattice_exponent(vertices.iter().map(|point| point[axis].integer_decode()))
    });
    let to_lattice = |point: &[f64; 2]| {
        std::array::from_fn(|axis| {
            binary_lattice_coordinate(point[axis].integer_decode(), exponents[axis])
        })
    };
    let mut vertices = vertices.iter().map(to_lattice);
    let origin: [BigInt; 2] = vertices.next()?;
    let mut previous = origin.clone();
    let mut total = BigInt::ZERO;
    for next in vertices {
        total += &previous[0] * &next[1] - &previous[1] * &next[0];
        previous = next;
    }
    total += &previous[0] * &origin[1] - &previous[1] * &origin[0];
    scaled_lattice_magnitude(total, exponents[0] + exponents[1], 2)
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
///
/// Ill-conditioned and extreme-scale inputs use an exact determinant of the
/// represented coordinates. Fewer than twelve coordinates give zero. Only
/// the first twelve coordinates are used.
#[must_use]
pub fn tetrahedron_volume(points: &[f64]) -> f64 {
    // A ratio above 2^-20 leaves margin for an accurate volume magnitude.
    const MIN_RATIO: f64 = f64::from_bits((1023 - 20) << 52);
    let Some(points) = points.first_chunk::<12>() else {
        return 0.0;
    };
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
    let positive = [b[1] * c[2], b[2] * c[0], b[0] * c[1]];
    let negative = [b[2] * c[1], b[0] * c[2], b[1] * c[0]];
    let cross: [f64; 3] = std::array::from_fn(|i| positive[i] - negative[i]);
    let determinant = a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2];
    let magnitude = determinant.abs();
    let abs_a = a.map(f64::abs);
    // The permanent measures cancellation across the six triple products.
    let permanent = abs_a[0] * (positive[0].abs() + negative[0].abs())
        + abs_a[1] * (positive[1].abs() + negative[1].abs())
        + abs_a[2] * (positive[2].abs() + negative[2].abs());

    // Small product errors can grow when the outer edge is large. A normal
    // increment bounds the combined subnormal rounding loss with ample margin.
    let underflow_floor = (abs_a[0] + abs_a[1] + abs_a[2] + 1.0) * f64::MIN_POSITIVE;
    if magnitude.is_finite() && magnitude > MIN_RATIO * permanent && magnitude > underflow_floor {
        return magnitude / 6.0;
    }
    if permanent == 0.0 && (0..3).any(|axis| a[axis] == 0.0 && b[axis] == 0.0 && c[axis] == 0.0) {
        return 0.0;
    }
    exact_tetrahedron_volume(points).unwrap_or(magnitude / 6.0)
}

#[cold]
fn exact_tetrahedron_volume(points: &[f64; 12]) -> Option<f64> {
    use num_bigint::BigInt;
    use num_traits::Float;

    if points.iter().any(|coordinate| !coordinate.is_finite()) {
        return None;
    }
    let decoded = points.map(Float::integer_decode);
    // Put the coordinates on a common binary lattice so each difference
    // and product is an exact integer operation.
    let exponent = binary_lattice_exponent(decoded);
    let coordinates = decoded.map(|coordinate| binary_lattice_coordinate(coordinate, exponent));
    let [a, b, c]: [[BigInt; 3]; 3] = [3, 6, 9]
        .map(|offset| std::array::from_fn(|axis| &coordinates[offset + axis] - &coordinates[axis]));
    let determinant = &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
        + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
        + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0]);

    // The determinant scales with the cube of the lattice spacing.
    // Restore that factor and divide by six before rounding once.
    scaled_lattice_magnitude(determinant, 3 * exponent, 6)
}

fn binary_lattice_exponent(coordinates: impl IntoIterator<Item = (u64, i16, i8)>) -> i32 {
    coordinates
        .into_iter()
        .filter(|(mantissa, _, _)| *mantissa != 0)
        .map(|(_, power, _)| i32::from(power))
        .min()
        .unwrap_or(0)
}

fn binary_lattice_coordinate(
    (mantissa, power, sign): (u64, i16, i8),
    exponent: i32,
) -> num_bigint::BigInt {
    let magnitude =
        num_bigint::BigInt::from(mantissa) << (i32::from(power) - exponent).max(0).unsigned_abs();
    if sign < 0 {
        -magnitude
    } else {
        magnitude
    }
}

fn scaled_lattice_magnitude(
    value: num_bigint::BigInt,
    power: i32,
    divisor: u8,
) -> Option<f64> {
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::ToPrimitive;

    let mut numerator = BigInt::from(value.into_parts().1);
    let mut denominator = BigInt::from(divisor);
    if power < 0 {
        denominator <<= power.unsigned_abs();
    } else {
        numerator <<= power.unsigned_abs();
    }
    BigRational::new(numerator, denominator).to_f64()
}

#[cfg(test)]
mod tetrahedron_tests {
    use super::tetrahedron_volume;

    const UNIT_AXES: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

    #[test]
    fn tetrahedron_volume_preserves_the_flat_input_contract() {
        assert_eq!(tetrahedron_volume(&UNIT_AXES), 1.0 / 6.0);
        assert_eq!(tetrahedron_volume(&UNIT_AXES[..11]), 0.0);
        let mut extra = UNIT_AXES.to_vec();
        extra.push(f64::NAN);
        assert_eq!(tetrahedron_volume(&extra), 1.0 / 6.0);
        let mut coplanar = UNIT_AXES;
        coplanar[9..].copy_from_slice(&[1.0, 1.0, 0.0]);
        assert_eq!(tetrahedron_volume(&coplanar), 0.0);
        let mut nonfinite = UNIT_AXES;
        nonfinite[3] = f64::INFINITY;
        assert_eq!(tetrahedron_volume(&nonfinite), f64::INFINITY);
        nonfinite[3] = f64::NAN;
        assert!(tetrahedron_volume(&nonfinite).is_nan());
    }

    #[test]
    fn tetrahedron_volume_recovers_representable_extreme_results() {
        // Expected values use rational arithmetic on the binary coordinates.
        let cases = [
            (
                [
                    -1e308, 0.0, 0.0, 1e308, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1e-300,
                ],
                33_333_333.333333336,
            ),
            (
                [
                    0.0, 0.0, 0.0, 1e-300, 0.0, 0.0, 0.0, 1e200, 0.0, 0.0, 0.0, 1e200,
                ],
                1.6666666666666666e99,
            ),
            (
                [0.0, 0.0, 0.0, 1e308, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
                6.6666666666666664e307,
            ),
            (
                [
                    0.0, 0.0, 0.0, 1e300, 0.0, 0.0, 0.0, 1e-300, 0.0, 0.0, 0.0, 1e-300,
                ],
                1.6666666666666669e-301,
            ),
            (
                [
                    0.0, 0.0, 0.0, 1e300, 0.0, 1.0, 1e-200, 1e-200, 0.0, 0.0, 1e100, 1e-200,
                ],
                3.3333333333333336e-101,
            ),
        ];
        for (points, expected) in cases {
            assert_eq!(tetrahedron_volume(&points), expected, "{points:?}");
        }
        for (power, bits) in [(356, 11), (358, 0)] {
            let scale = f64::from_bits((1023 - power) << 52);
            let points = UNIT_AXES.map(|x| x * scale);
            assert_eq!(tetrahedron_volume(&points).to_bits(), bits);
        }
        let huge = f64::from_bits((1023 + 550) << 52);
        let mut mixed = UNIT_AXES;
        mixed[3] = f64::from_bits(1);
        mixed[7] = huge;
        mixed[11] = huge;
        assert_eq!(tetrahedron_volume(&mixed), 67_108_864.0 / 6.0);
    }

    #[test]
    fn tetrahedron_volume_recovers_near_coplanar_magnitudes() {
        let magnitude = 134_217_728.0; // 2^27
        for offset in [1.0, 101.0] {
            let points = [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                magnitude + offset,
                magnitude,
                0.0,
                magnitude + 2.0 * offset,
                magnitude + offset,
                0.0,
            ];
            let mut opposite_winding = points;
            for axis in 0..3 {
                opposite_winding.swap(6 + axis, 9 + axis);
            }
            let expected = offset * offset / 6.0;
            assert_eq!(tetrahedron_volume(&points), expected);
            assert_eq!(tetrahedron_volume(&opposite_winding), expected);
        }
    }
}
