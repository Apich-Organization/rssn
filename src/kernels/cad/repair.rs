//! Mesh welding, orientation repair and hole filling.
//!
//! [`weld`] merges vertices that coincide within a tolerance and drops the
//! triangles that collapse. [`repair`] then makes the surface consistently
//! oriented and closes every boundary loop with a triangle fan.
//!
//! # Limitations
//!
//! * [`weld`] uses a uniform grid of cells of side `tolerance`, so its cost is
//!   linear in the vertex count for a well-spread mesh and quadratic in the
//!   worst case.
//! * [`repair`] fills each boundary loop with a fan from its first vertex. A
//!   large, highly non-convex hole therefore gets long, thin triangles.
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
use super::types::sub;

/// What [`repair`] changed in the mesh it was given.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RepairReport {
    /// Vertices merged away by the weld pass.
    pub merged_vertices: usize,
    /// Triangles dropped because they collapsed to a segment or a point.
    pub removed_degenerate: usize,
    /// Boundary loops closed with a triangle fan.
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
/// remaining boundary loop is closed with a fan.
#[must_use]
pub fn repair(
    mesh: &Mesh,
    tolerance: f64,
) -> (Mesh, RepairReport) {
    let (mut welded, merged, removed) = weld_impl(mesh, tolerance);
    let mut flipped = orient_consistently(&mut welded);
    if is_closed(&welded) && signed_volume(&welded) < 0.0 {
        for t in &mut welded.triangles {
            *t = [t[0], t[2], t[1]];
        }
        flipped += welded.triangles.len();
    }
    let filled = fill_holes(&mut welded);

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

/// Whether every undirected edge of `mesh` is shared by exactly two triangles.
#[must_use]
pub fn is_closed(mesh: &Mesh) -> bool {
    if mesh.triangles.is_empty() {
        return false;
    }

    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            *counts.entry(edge_key(t[k], t[(k + 1) % 3])).or_insert(0) += 1;
        }
    }

    counts.values().all(|&count| count == 2)
}

/// Whether `mesh` is edge-manifold: no edge is shared by more than two
/// triangles and no triangle repeats a vertex.
#[must_use]
pub fn is_manifold(mesh: &Mesh) -> bool {
    if !super::types::indices_valid(mesh) {
        return false;
    }

    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return false;
        }
        for k in 0..3 {
            *counts.entry(edge_key(t[k], t[(k + 1) % 3])).or_insert(0) += 1;
        }
    }

    counts.values().all(|&count| count <= 2)
}

/// The signed volume enclosed by `mesh`, positive for outward-facing normals.
#[must_use]
pub fn signed_volume(mesh: &Mesh) -> f64 {
    let mut total = 0.0;
    for t in &mesh.triangles {
        let (Some(a), Some(b), Some(c)) = (
            mesh.vertices.get(t[0]),
            mesh.vertices.get(t[1]),
            mesh.vertices.get(t[2]),
        ) else {
            continue;
        };
        total += dot(*a, cross(*b, *c));
    }
    total / 6.0
}

/// The Euler characteristic `V - E + F` of `mesh`.
#[must_use]
pub fn euler_characteristic(mesh: &Mesh) -> i64 {
    let mut edges: HashSet<(usize, usize)> = HashSet::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            edges.insert(edge_key(t[k], t[(k + 1) % 3]));
        }
    }

    mesh.vertices.len() as i64 - edges.len() as i64 + mesh.triangles.len() as i64
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
                        let Some(list) = grid.get(&(base.0 + dx, base.1 + dy, base.2 + dz)) else {
                            continue;
                        };
                        for &index in list {
                            let q = vertices[index];
                            let d = [q[0] - v[0], q[1] - v[1], q[2] - v[2]];
                            if dot(d, d) <= limit {
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
        if length(cross(
            sub(vertices[b], vertices[a]),
            sub(vertices[c], vertices[a]),
        )) <= 1e-14
        {
            continue;
        }
        triangles.push([a, b, c]);
    }

    let merged = mesh.vertices.len().saturating_sub(vertices.len());
    let removed = mesh.triangles.len().saturating_sub(triangles.len());

    (Mesh { vertices, triangles }, merged, removed)
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

/// Closes every boundary loop with a triangle fan.
///
/// Returns the number of loops that were filled.
fn fill_holes(mesh: &mut Mesh) -> usize {
    let mut directed: HashSet<(usize, usize)> = HashSet::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            directed.insert((t[k], t[(k + 1) % 3]));
        }
    }

    // An unmatched edge `a -> b` of the surface is closed by a new triangle
    // that runs `b -> a`, so the hole loop follows `b` to `a`.
    let mut next: HashMap<usize, usize> = HashMap::new();
    for &(a, b) in &directed {
        if !directed.contains(&(b, a)) {
            next.insert(b, a);
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
        let mut current = start;
        let mut closed = false;
        let mut guard = 0;

        while let Some(&following) = next.get(&current) {
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

        for i in 1..ring.len() - 1 {
            mesh.triangles.push([ring[0], ring[i], ring[i + 1]]);
        }
        filled += 1;
    }

    filled
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
}
