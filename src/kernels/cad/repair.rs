//! Mesh welding, orientation repair and optimal hole filling.
//!
//! [`weld`] merges vertices that coincide within a tolerance and drops the
//! triangles that collapse. [`repair`] then makes the surface consistently
//! oriented and closes every boundary loop using Liepa's (2003) minimum-weight
//! dynamic-programming hole triangulation (with an ear-clipping fallback for
//! very large loops).
//!
//! # Limitations
//!
//! * [`weld`] uses a uniform grid of cells of side `tolerance`, so its cost is
//!   linear in the vertex count for a well-spread mesh and quadratic in the
//!   worst case.
//! * A non-orientable surface (a Möbius strip) cannot be oriented; the pass
//!   stops at the first inconsistency and reports the flips it made.

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use super::types::Mesh;
use super::types::cross;
use super::types::dot;
use super::types::edge_key;
use super::types::length;
use super::types::normalize;
use super::types::sub;
use super::types::triangle_normal;

/// What [`repair`] changed in the mesh it was given.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RepairReport {
    /// Vertices merged away by the weld pass.
    pub merged_vertices: usize,
    /// Triangles dropped because they collapsed to a segment or a point.
    pub removed_degenerate: usize,
    /// Boundary loops closed by optimal hole triangulation.
    pub filled_holes: usize,
    /// Triangles whose winding was reversed to match their neighbours.
    pub flipped: usize,
}

/// Merges vertices closer than `tolerance` and drops degenerate triangles.
///
/// Returns the welded mesh and the number of vertices that disappeared. A
/// non-positive or non-finite `tolerance` falls back to exact matching.
#[must_use]
pub fn weld(
    mesh: &Mesh,
    tolerance: f64,
) -> (Mesh, usize) {
    let (welded, merged, _) = weld_impl(mesh, tolerance);
    (welded, merged)
}

/// Welds, re-orients and closes `mesh`.
///
/// The passes run in that order: coincident vertices are merged, triangles that
/// collapse are dropped, the surface is oriented consistently, and every
/// remaining boundary loop is closed using Liepa's minimum-weight triangulation.
#[must_use]
pub fn repair(
    mesh: &Mesh,
    tolerance: f64,
) -> (Mesh, RepairReport) {
    let (mut welded, merged, removed) = weld_impl(mesh, tolerance);
    let mut flipped = orient_consistently(&mut welded);
    if is_closed(&welded) && closed_signed_volume(&welded) < 0.0 {
        for t in &mut welded.triangles {
            *t = [t[0], t[2], t[1]];
        }
        flipped += welded.triangles.len();
    }
    let filled = fill_holes(&mut welded);
    if filled > 0 && is_closed(&welded) && closed_signed_volume(&welded) < 0.0 {
        for t in &mut welded.triangles {
            *t = [t[0], t[2], t[1]];
        }
        flipped += welded.triangles.len();
    }

    (
        welded,
        RepairReport {
            merged_vertices: merged,
            removed_degenerate: removed,
            filled_holes: filled,
            flipped,
        },
    )
}

/// Whether every undirected edge belongs to exactly two faces of a valid mesh.
#[must_use]
pub fn is_closed(mesh: &Mesh) -> bool {
    topology(mesh).closed
}

/// Whether `mesh` is edge-manifold: no edge is shared by more than two
/// triangles and no triangle repeats a vertex.
#[must_use]
pub fn is_manifold(mesh: &Mesh) -> bool {
    topology(mesh).manifold
}

/// The Euler characteristic `V - E + F` of `mesh`.
#[must_use]
pub fn euler_characteristic(mesh: &Mesh) -> i64 {
    topology(mesh).euler
}

#[derive(Default)]
struct EdgeIncidence {
    faces: usize,
    direction: i64,
}

pub(crate) struct Topology {
    pub(crate) closed: bool,
    pub(crate) manifold: bool,
    pub(crate) oriented: bool,
    pub(crate) euler: i64,
}

/// Reads the shared edge incidence once, so callers can retain the scalar summary.
pub(crate) fn topology(mesh: &Mesh) -> Topology {
    let mut edges: HashMap<_, EdgeIncidence> = HashMap::new();
    let mut valid_faces = true;
    for &[a, b, c] in &mesh.triangles {
        valid_faces &= a != b
            && b != c
            && a != c
            && a < mesh.vertices.len()
            && b < mesh.vertices.len()
            && c < mesh.vertices.len();
        for (start, end) in [(a, b), (b, c), (c, a)] {
            let incidence = edges.entry(edge_key(start, end)).or_default();
            incidence.faces += 1;
            incidence.direction += if start <= end { 1 } else { -1 };
        }
    }
    let mut closed = valid_faces && !edges.is_empty();
    let mut manifold = valid_faces;
    let mut oriented = true;
    for incidence in edges.values() {
        closed &= incidence.faces == 2;
        manifold &= incidence.faces <= 2;
        oriented &= incidence.direction == 0;
    }
    Topology {
        closed,
        manifold,
        oriented: closed && oriented,
        euler: mesh.vertices.len() as i64 - edges.len() as i64 + mesh.triangles.len() as i64,
    }
}

/// The signed volume enclosed by `mesh`, positive for outward-facing normals.
#[must_use]
pub fn signed_volume(mesh: &Mesh) -> f64 {
    if topology(mesh).oriented {
        closed_signed_volume(mesh)
    } else {
        volume_about(mesh, [0.0; 3])
    }
}

pub(crate) fn closed_signed_volume(mesh: &Mesh) -> f64 {
    // Translating a closed surface preserves volume and avoids cancellation
    // between large coordinates when a small solid is far from the origin.
    volume_about(mesh, mesh.vertices.first().copied().unwrap_or([0.0; 3]))
}

fn volume_about(
    mesh: &Mesh,
    origin: [f64; 3],
) -> f64 {
    crate::kernels::par::sum_kahan_fixed(mesh.triangles.len(), 1024, 8192, |index| {
        let Some(&[a, b, c]) = mesh.triangles.get(index) else {
            return 0.0;
        };
        let (Some(a), Some(b), Some(c)) = (
            mesh.vertices.get(a),
            mesh.vertices.get(b),
            mesh.vertices.get(c),
        ) else {
            return 0.0;
        };
        dot(sub(*a, origin), cross(sub(*b, origin), sub(*c, origin)))
    }) / 6.0
}

/// Welds vertices, returning the mesh, the merged count and the dropped count.
fn weld_impl(
    mesh: &Mesh,
    tolerance: f64,
) -> (Mesh, usize, usize) {
    let tol = if tolerance.is_finite() && tolerance > 0.0 {
        tolerance
    } else {
        0.0
    };

    let mut vertices: Vec<[f64; 3]> = Vec::with_capacity(mesh.vertices.len());
    let mut remap: Vec<usize> = Vec::with_capacity(mesh.vertices.len());

    if tol > 0.0 {
        let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        let cell = |v: [f64; 3]| -> (i64, i64, i64) {
            (
                (v[0] / tol).floor() as i64,
                (v[1] / tol).floor() as i64,
                (v[2] / tol).floor() as i64,
            )
        };
        let limit = tol * tol;

        for v in &mesh.vertices {
            let base = cell(*v);
            let mut found: Option<usize> = None;
            'search: for dx in -1_i64..=1 {
                for dy in -1_i64..=1 {
                    for dz in -1_i64..=1 {
                        let (Some(x), Some(y), Some(z)) = (
                            base.0.checked_add(dx),
                            base.1.checked_add(dy),
                            base.2.checked_add(dz),
                        ) else {
                            continue;
                        };
                        let Some(list) = grid.get(&(x, y, z)) else {
                            continue;
                        };
                        for &index in list {
                            let q = vertices[index];
                            let d = [q[0] - v[0], q[1] - v[1], q[2] - v[2]];
                            let in_range = if limit.is_normal() {
                                dot(d, d) <= limit
                            } else {
                                let scaled = d.map(|delta| delta / tol);
                                dot(scaled, scaled) <= 1.0
                            };
                            if in_range {
                                found = Some(index);
                                break 'search;
                            }
                        }
                    }
                }
            }

            match found {
                | Some(index) => remap.push(index),
                | None => {
                    let index = vertices.len();
                    vertices.push(*v);
                    grid.entry(base).or_default().push(index);
                    remap.push(index);
                },
            }
        }
    } else {
        let mut exact: HashMap<[u64; 3], usize> = HashMap::new();
        for v in &mesh.vertices {
            let key = v.map(f64::to_bits);
            match exact.get(&key) {
                | Some(&index) => remap.push(index),
                | None => {
                    let index = vertices.len();
                    vertices.push(*v);
                    exact.insert(key, index);
                    remap.push(index);
                },
            }
        }
    }

    let mut triangles: Vec<[usize; 3]> = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        let (Some(&a), Some(&b), Some(&c)) = (remap.get(t[0]), remap.get(t[1]), remap.get(t[2]))
        else {
            continue;
        };
        if a == b || b == c || c == a {
            continue;
        }
        if degenerate_triangle(vertices[a], vertices[b], vertices[c]) {
            continue;
        }
        triangles.push([a, b, c]);
    }

    let merged = mesh.vertices.len().saturating_sub(vertices.len());
    let removed = mesh.triangles.len().saturating_sub(triangles.len());

    (Mesh { vertices, triangles }, merged, removed)
}

fn degenerate_triangle(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> bool {
    // The sine of the angle between edges is independent of coordinate units.
    // Scale each vector before forming products so tiny triangles stay valid.
    let scaled = |edge: [f64; 3]| {
        let scale = edge.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
        if scale > 0.0 && scale.is_finite() {
            edge.map(|component| component / scale)
        } else {
            [0.0; 3]
        }
    };
    let left = scaled(sub(b, a));
    let right = scaled(sub(c, a));
    length(cross(left, right)) <= 1e-14
}

/// The direction in which `t` traverses the edge `key`.
///
/// Returns `Some(true)` when the triangle runs from the smaller to the larger
/// index of `key`.
fn edge_forward(
    t: [usize; 3],
    key: (usize, usize),
) -> Option<bool> {
    for k in 0..3 {
        let a = t[k];
        let b = t[(k + 1) % 3];
        if edge_key(a, b) == key {
            return Some(a == key.0);
        }
    }
    None
}

/// Reverses triangles so that shared edges are traversed in opposite directions.
///
/// Returns the number of triangles that were flipped.
fn orient_consistently(mesh: &mut Mesh) -> usize {
    let count = mesh.triangles.len();
    if count == 0 {
        return 0;
    }

    let mut adjacency: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (index, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            adjacency
                .entry(edge_key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(index);
        }
    }

    let mut visited = vec![false; count];
    let mut flipped = 0;
    let mut queue: VecDeque<usize> = VecDeque::new();

    for seed in 0..count {
        if visited[seed] {
            continue;
        }
        visited[seed] = true;
        queue.push_back(seed);

        while let Some(i) = queue.pop_front() {
            for k in 0..3 {
                let a = mesh.triangles[i][k];
                let b = mesh.triangles[i][(k + 1) % 3];
                let key = edge_key(a, b);
                let forward_i = a == key.0;

                let Some(list) = adjacency.get(&key) else {
                    continue;
                };
                let neighbours = list.clone();
                for j in neighbours {
                    if j == i || visited[j] {
                        continue;
                    }
                    if edge_forward(mesh.triangles[j], key) == Some(forward_i) {
                        let t = mesh.triangles[j];
                        mesh.triangles[j] = [t[0], t[2], t[1]];
                        flipped += 1;
                    }
                    visited[j] = true;
                    queue.push_back(j);
                }
            }
        }
    }

    flipped
}

/// Maximum boundary loop length solved by Liepa's `O(k^3)` DP before falling
/// back to 2D ear-clipping.
const DP_MAX_LOOP: usize = 128;

/// Closes every boundary loop using Liepa's (2003) minimum-weight dynamic
/// programming triangulation (or 2D ear-clipping for loops longer than 128
/// vertices).
///
/// Returns the number of loops that were filled.
fn fill_holes(mesh: &mut Mesh) -> usize {
    let mut directed: HashMap<(usize, usize), Option<[f64; 3]>> = HashMap::new();
    for &t in &mesh.triangles {
        let normal = triangle_normal(mesh, t);
        for k in 0..3 {
            directed.insert((t[k], t[(k + 1) % 3]), normal);
        }
    }

    // An unmatched edge `a -> b` of the surface is closed by a new triangle
    // that runs `b -> a`, so the hole loop follows `b` to `a`, and the
    // adjacent existing face across `(b, a)` has normal `directed[(a, b)]`.
    let mut next: HashMap<usize, (usize, Option<[f64; 3]>)> = HashMap::new();
    for (&(a, b), &normal) in &directed {
        if !directed.contains_key(&(b, a)) && next.insert(b, (a, normal)).is_some() {
            // Pinched boundary vertices do not determine a unique loop.
            return 0;
        }
    }
    if next.is_empty() {
        return 0;
    }

    let mut starts: Vec<usize> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut visited: HashSet<usize> = HashSet::new();
    let mut filled = 0;

    for start in starts {
        if !visited.insert(start) {
            continue;
        }

        let mut ring = vec![start];
        let mut adj_normals: Vec<Option<[f64; 3]>> = Vec::new();
        let mut current = start;
        let mut closed = false;
        let mut guard = 0;

        while let Some(&(following, adj_normal)) = next.get(&current) {
            adj_normals.push(adj_normal);
            if following == start {
                closed = true;
                break;
            }
            if !visited.insert(following) {
                break;
            }
            ring.push(following);
            current = following;
            guard += 1;
            if guard > next.len() + 2 {
                break;
            }
        }

        if !closed || ring.len() < 3 {
            continue;
        }

        let first_new_face = mesh.triangles.len();
        if ring.len() == 3 {
            mesh.triangles.push([ring[0], ring[1], ring[2]]);
        } else if ring.len() <= DP_MAX_LOOP {
            triangulate_hole_dp(&mesh.vertices, &ring, &adj_normals, &mut mesh.triangles);
        } else {
            triangulate_hole_ear_clipping(&mesh.vertices, &ring, &mut mesh.triangles);
        }
        if mesh.triangles.len() - first_new_face == ring.len() - 2 {
            filled += 1;
        } else {
            mesh.triangles.truncate(first_new_face);
        }
    }

    filled
}

/// Newell's method average outward normal for an ordered 3D polygon loop.
fn newell_normal(
    vertices: &[[f64; 3]],
    ring: &[usize],
) -> [f64; 3] {
    let mut n = [0.0; 3];
    let k = ring.len();
    for i in 0..k {
        let a = vertices[ring[i]];
        let b = vertices[ring[(i + 1) % k]];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    normalize(n)
}

/// Dihedral angle in `[0, pi]` between two unit normals.
fn dihedral_angle(
    n1: [f64; 3],
    n2: [f64; 3],
) -> f64 {
    dot(n1, n2).clamp(-1.0, 1.0).acos()
}

/// Lexicographic DP cost for Liepa's (2003) hole triangulation:
/// 1. Minimize worst dihedral / inversion / degeneracy penalty.
/// 2. Minimize total patch surface area.
/// 3. Break flat ties by minimizing worst triangle aspect ratio.
#[derive(Clone, Copy, Debug)]
struct HoleWeight {
    max_angle: f64,
    area: f64,
    max_aspect: f64,
}

impl HoleWeight {
    const INF: Self = Self {
        max_angle: f64::INFINITY,
        area: f64::INFINITY,
        max_aspect: f64::INFINITY,
    };
    const ZERO: Self = Self {
        max_angle: 0.0,
        area: 0.0,
        max_aspect: 1.0,
    };

    fn is_better(
        &self,
        other: &Self,
        area_eps: f64,
    ) -> bool {
        if self.max_angle + 1e-9 < other.max_angle {
            return true;
        }
        if other.max_angle + 1e-9 < self.max_angle {
            return false;
        }
        if self.area + area_eps < other.area {
            return true;
        }
        if other.area + area_eps < self.area {
            return false;
        }
        self.max_aspect < other.max_aspect
    }
}

/// Liepa's (2003) minimum-weight dynamic-programming hole triangulation.
fn triangulate_hole_dp(
    vertices: &[[f64; 3]],
    ring: &[usize],
    adj_normals: &[Option<[f64; 3]>],
    out: &mut Vec<[usize; 3]>,
) {
    let k = ring.len();
    let loop_normal = newell_normal(vertices, ring);
    let has_loop_normal = length(loop_normal) > 0.5;

    // Estimate perimeter scale for tie-breaking area epsilon.
    let mut perim_sum = 0.0;
    for i in 0..k {
        perim_sum += length(sub(vertices[ring[(i + 1) % k]], vertices[ring[i]]));
    }
    let area_eps = (perim_sum * perim_sum * 1e-12).max(1e-12);

    let mut dp_weight = vec![HoleWeight::INF; k * k];
    let mut dp_split = vec![0_usize; k * k];
    let mut dp_normal: Vec<Option<[f64; 3]>> = vec![None; k * k];

    for i in 0..k - 1 {
        dp_weight[i * k + (i + 1)] = HoleWeight::ZERO;
    }

    for span in 2..k {
        for i in 0..(k - span) {
            let j = i + span;
            let a = vertices[ring[i]];
            let c = vertices[ring[j]];

            let mut best_w = HoleWeight::INF;
            let mut best_m = i + 1;
            let mut best_n: Option<[f64; 3]> = None;

            for m in (i + 1)..j {
                let b = vertices[ring[m]];
                let cr = cross(sub(b, a), sub(c, a));
                let area2 = length(cr);
                let area = 0.5 * area2;

                let (tri_angle, tri_aspect, tri_normal) = if area2 <= 1e-14 {
                    (100.0, 1e9, None)
                } else {
                    let n = [cr[0] / area2, cr[1] / area2, cr[2] / area2];
                    let mut angle = 0.0_f64;

                    if has_loop_normal {
                        let align = dot(n, loop_normal);
                        if align < -1e-6 {
                            angle = angle.max(2.0 * std::f64::consts::PI + (1.0 - align));
                        }
                    }

                    let n_im = if m == i + 1 {
                        adj_normals[i]
                    } else {
                        dp_normal[i * k + m]
                    };
                    if let Some(nb) = n_im {
                        angle = angle.max(dihedral_angle(n, nb));
                    }

                    let n_mj = if j == m + 1 {
                        adj_normals[m]
                    } else {
                        dp_normal[m * k + j]
                    };
                    if let Some(nb) = n_mj {
                        angle = angle.max(dihedral_angle(n, nb));
                    }

                    if i == 0 && j == k - 1 {
                        if let Some(nb) = adj_normals[k - 1] {
                            angle = angle.max(dihedral_angle(n, nb));
                        }
                    }

                    let p = length(sub(b, a)) + length(sub(c, b)) + length(sub(a, c));
                    let aspect = (p * p) / (4.0 * 3.0_f64.sqrt() * area);
                    (angle, aspect, Some(n))
                };

                let w_im = dp_weight[i * k + m];
                let w_mj = dp_weight[m * k + j];
                let cand = HoleWeight {
                    max_angle: tri_angle.max(w_im.max_angle).max(w_mj.max_angle),
                    area: area + w_im.area + w_mj.area,
                    max_aspect: tri_aspect.max(w_im.max_aspect).max(w_mj.max_aspect),
                };

                if cand.is_better(&best_w, area_eps) {
                    best_w = cand;
                    best_m = m;
                    best_n = tri_normal;
                }
            }

            dp_weight[i * k + j] = best_w;
            dp_split[i * k + j] = best_m;
            dp_normal[i * k + j] = best_n;
        }
    }

    reconstruct_dp_triangles(ring, &dp_split, k, 0, k - 1, out);
}

/// Recursively emits triangles from Liepa's DP split table.
fn reconstruct_dp_triangles(
    ring: &[usize],
    split: &[usize],
    k: usize,
    i: usize,
    j: usize,
    out: &mut Vec<[usize; 3]>,
) {
    if j <= i + 1 {
        return;
    }
    let m = split[i * k + j];
    out.push([ring[i], ring[m], ring[j]]);
    reconstruct_dp_triangles(ring, split, k, i, m, out);
    reconstruct_dp_triangles(ring, split, k, m, j, out);
}

/// 2D ear-clipping fallback for very large boundary loops (`k > 128`).
fn triangulate_hole_ear_clipping(
    vertices: &[[f64; 3]],
    ring: &[usize],
    out: &mut Vec<[usize; 3]>,
) {
    let mut n = newell_normal(vertices, ring);
    if length(n) < 0.5 {
        n = [0.0, 0.0, 1.0];
    }
    let ref_vec = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(ref_vec, n));
    let v = cross(n, u);
    let origin = vertices[ring[0]];

    let mut poly: Vec<(usize, [f64; 2])> = ring
        .iter()
        .map(|&idx| {
            let rel = sub(vertices[idx], origin);
            (idx, [dot(rel, u), dot(rel, v)])
        })
        .collect();

    let cross2 = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| -> f64 {
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };

    let point_in_tri = |p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]| -> bool {
        let c1 = cross2(a, b, p);
        let c2 = cross2(b, c, p);
        let c3 = cross2(c, a, p);
        c1 >= -1e-12 && c2 >= -1e-12 && c3 >= -1e-12
    };

    while poly.len() > 3 {
        let len = poly.len();
        let mut best_ear: Option<(usize, f64)> = None;

        for i in 0..len {
            let prev = poly[(i + len - 1) % len].1;
            let curr = poly[i].1;
            let next = poly[(i + 1) % len].1;
            let area2 = cross2(prev, curr, next);
            if area2 <= 1e-12 {
                continue;
            }

            let mut empty = true;
            for j in 0..len {
                if j == i || j == (i + len - 1) % len || j == (i + 1) % len {
                    continue;
                }
                if point_in_tri(poly[j].1, prev, curr, next) {
                    empty = false;
                    break;
                }
            }
            if !empty {
                continue;
            }

            let d0 = ((curr[0] - prev[0]).hypot(curr[1] - prev[1]))
                + ((next[0] - curr[0]).hypot(next[1] - curr[1]))
                + ((prev[0] - next[0]).hypot(prev[1] - next[1]));
            let quality = area2 / (d0 * d0).max(1e-18);
            if best_ear.is_none_or(|(_, q)| quality > q) {
                best_ear = Some((i, quality));
            }
        }

        let ear_idx = match best_ear {
            | Some((idx, _)) => idx,
            | None => {
                (0..len)
                    .max_by(|&a, &b| {
                        let ca = cross2(
                            poly[(a + len - 1) % len].1,
                            poly[a].1,
                            poly[(a + 1) % len].1,
                        );
                        let cb = cross2(
                            poly[(b + len - 1) % len].1,
                            poly[b].1,
                            poly[(b + 1) % len].1,
                        );
                        ca.total_cmp(&cb)
                    })
                    .unwrap_or(0)
            },
        };

        let prev_id = poly[(ear_idx + len - 1) % len].0;
        let curr_id = poly[ear_idx].0;
        let next_id = poly[(ear_idx + 1) % len].0;
        out.push([prev_id, curr_id, next_id]);
        poly.remove(ear_idx);
    }

    if poly.len() == 3 {
        out.push([poly[0].0, poly[1].0, poly[2].0]);
    }
}

#[cfg(test)]
mod tests {
    use super::RepairReport;
    use super::euler_characteristic;
    use super::is_closed;
    use super::is_manifold;
    use super::repair;
    use super::signed_volume;
    use super::weld;
    use crate::kernels::cad::types::Mesh;

    fn cube_faces() -> Vec<[usize; 3]> {
        vec![
            [0, 3, 2],
            [0, 2, 1],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 7, 6],
            [3, 6, 2],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ]
    }

    fn cube_vertices() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    }

    fn cube() -> Mesh {
        Mesh::new(cube_vertices(), cube_faces())
    }

    /// The same cube with a fresh copy of each vertex per triangle.
    fn cube_soup() -> Mesh {
        let vertices = cube_vertices();
        let mut soup = Vec::with_capacity(36);
        let mut triangles = Vec::with_capacity(12);
        for t in cube_faces() {
            let base = soup.len();
            soup.push(vertices[t[0]]);
            soup.push(vertices[t[1]]);
            soup.push(vertices[t[2]]);
            triangles.push([base, base + 1, base + 2]);
        }
        Mesh::new(soup, triangles)
    }

    #[test]
    fn weld_merges_the_cube_soup() {
        let soup = cube_soup();
        assert_eq!(soup.vertices.len(), 36);

        let (welded, merged) = weld(&soup, 1e-9);
        assert_eq!(welded.vertices.len(), 8);
        assert_eq!(merged, 28);
        assert_eq!(welded.triangles.len(), 12);
        assert!((signed_volume(&welded) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn weld_drops_degenerate_triangles() {
        let mesh = Mesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 1, 3]],
        );
        let (welded, merged) = weld(&mesh, 1e-9);
        assert_eq!(welded.vertices.len(), 3);
        assert_eq!(merged, 1);
        assert_eq!(welded.triangles.len(), 1);
    }

    #[test]
    fn repair_fills_a_missing_triangle() {
        let mut mesh = cube();
        mesh.triangles.remove(0);
        assert!(!is_closed(&mesh));

        let (repaired, report) = repair(&mesh, 1e-9);
        assert!(is_closed(&repaired), "hole was not filled");
        assert_eq!(report.filled_holes, 1);
        assert_eq!(repaired.triangles.len(), 12);
        assert!((signed_volume(&repaired) - 1.0).abs() < 1e-12);
        assert_eq!(euler_characteristic(&repaired), 2);
    }

    #[test]
    fn repair_reorients_a_flipped_triangle() {
        let mut mesh = cube();
        let flipped = mesh.triangles[0];
        mesh.triangles[0] = [flipped[0], flipped[2], flipped[1]];

        let (repaired, report) = repair(&mesh, 1e-9);
        assert!(report.flipped >= 1);
        assert!((signed_volume(&repaired) - 1.0).abs() < 1e-12);
        assert!(is_closed(&repaired));
    }

    #[test]
    fn topology_queries_match_the_cube() {
        let mesh = cube();
        assert!(is_closed(&mesh));
        assert!(is_manifold(&mesh));
        assert_eq!(euler_characteristic(&mesh), 2);
        assert!((signed_volume(&mesh) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn open_mesh_is_neither_closed_nor_a_solid() {
        let mut mesh = cube();
        mesh.triangles.remove(3);
        assert!(!is_closed(&mesh));
        assert!(is_manifold(&mesh));
    }

    #[test]
    fn repair_of_an_empty_mesh_reports_nothing() {
        let (repaired, report) = repair(&Mesh::default(), 1e-9);
        assert_eq!(repaired.triangles.len(), 0);
        assert_eq!(report, RepairReport::default());
        assert!(!is_closed(&Mesh::default()));
    }

    /// Builds an extruded prism with a closed bottom and side walls and an open
    /// top face, leaving a single non-convex hole on `z = 1`.
    fn open_top_prism(
        profile: &[[f64; 2]],
        bottom_tris: &[[usize; 3]],
    ) -> Mesh {
        let n = profile.len();
        let mut vertices = Vec::with_capacity(2 * n);
        for &[x, y] in profile {
            vertices.push([x, y, 0.0]);
        }
        for &[x, y] in profile {
            vertices.push([x, y, 1.0]);
        }

        let mut triangles = Vec::with_capacity(bottom_tris.len() + 2 * n);
        triangles.extend_from_slice(bottom_tris);
        for i in 0..n {
            let j = (i + 1) % n;
            let b0 = i;
            let b1 = j;
            let t0 = n + i;
            let t1 = n + j;
            triangles.push([b0, b1, t1]);
            triangles.push([b0, t1, t0]);
        }

        Mesh::new(vertices, triangles)
    }

    #[test]
    fn repair_fills_non_convex_l_shaped_hole_without_inversions() {
        let profile = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let bottom = [[0, 2, 1], [0, 3, 2], [0, 4, 3], [0, 5, 4]];
        let open = open_top_prism(&profile, &bottom);
        let initial_tris = open.triangles.len();
        assert!(!is_closed(&open));

        let (repaired, report) = repair(&open, 1e-9);
        assert_eq!(report.filled_holes, 1);
        assert!(is_closed(&repaired), "L-shaped hole was not closed");
        assert!(is_manifold(&repaired), "repaired L-prism is not manifold");
        assert!(
            (signed_volume(&repaired) - 3.0).abs() < 1e-12,
            "volume was {}",
            signed_volume(&repaired)
        );

        // Every triangle added to close the L-shaped hole must be non-degenerate
        // and face strictly upward (+z).
        for &t in &repaired.triangles[initial_tris..] {
            let n = super::triangle_normal(&repaired, t)
                .expect("hole-fill triangle must not be degenerate");
            assert!(n[2] > 0.999, "inverted hole-fill triangle normal: {n:?}");
        }
    }

    #[test]
    fn repair_fills_non_convex_c_shaped_hole_without_inversions() {
        let profile = [
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [3.0, 2.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ];
        let bottom = [
            [0, 2, 1],
            [0, 3, 2],
            [0, 4, 3],
            [0, 7, 4],
            [7, 5, 4],
            [7, 6, 5],
        ];
        let open = open_top_prism(&profile, &bottom);
        let initial_tris = open.triangles.len();
        assert!(!is_closed(&open));

        let (repaired, report) = repair(&open, 1e-9);
        assert_eq!(report.filled_holes, 1);
        assert!(is_closed(&repaired), "C-shaped hole was not closed");
        assert!(is_manifold(&repaired), "repaired C-prism is not manifold");
        assert!(
            (signed_volume(&repaired) - 7.0).abs() < 1e-12,
            "volume was {}",
            signed_volume(&repaired)
        );

        for &t in &repaired.triangles[initial_tris..] {
            let n = super::triangle_normal(&repaired, t)
                .expect("hole-fill triangle must not be degenerate");
            assert!(n[2] > 0.999, "inverted hole-fill triangle normal: {n:?}");
        }
    }
}
