//! Three-dimensional Delaunay tetrahedralisation and tetrahedral mesh queries.
//!
//! [`delaunay_3d`] runs an `O(n log n)` expected-time Bowyer–Watson insertion:
//!
//! 1. Duplicate points are rejected in `O(n)` expected time via a 3D spatial
//!    hash grid.
//! 2. Input indices are ordered along a 3D Morton (Z-order) curve so
//!    consecutive insertions are spatially close.
//! 3. Each point is located in `O(1)` average steps by a visibility walk over
//!    the tetrahedron-to-tetrahedron adjacency graph, and its Bowyer–Watson
//!    cavity is discovered by local BFS using exact adaptive-precision
//!    `super::util::insphere_sos` predicates with Simulation of Simplicity
//!    (SoS) tie-breaking for cospherical configurations.

use std::collections::HashMap;
use std::collections::VecDeque;

use super::util::insphere_sos;
use super::util::orient3d;

/// Sentinel value indicating no adjacent tetrahedron across a face.
const NO_NEIGHBOR: usize = usize::MAX;

/// An unstructured tetrahedral mesh in three dimensions.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TetMesh {
    /// Vertex coordinates, one `[x, y, z]` per vertex.
    pub points: Vec<[f64; 3]>,
    /// Positively oriented tetrahedra, as quadruples of vertex indices.
    pub tetrahedra: Vec<[usize; 4]>,
}

/// Delaunay tetrahedralisation of `points`.
///
/// Returns positively oriented tetrahedra as quadruples of indices into
/// `points`.
///
/// # Errors
///
/// Returns `Err` when fewer than four points are supplied, when a coordinate is
/// not finite, or when two points coincide.
pub fn delaunay_3d(points: &[[f64; 3]]) -> Result<Vec<[usize; 4]>, String> {
    if points.len() < 4 {
        return Err("at least four points are required".to_string());
    }
    for p in points {
        if !p[0].is_finite() || !p[1].is_finite() || !p[2].is_finite() {
            return Err("point coordinates must be finite".to_string());
        }
    }

    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in points {
        for axis in 0..3 {
            lo[axis] = lo[axis].min(p[axis]);
            hi[axis] = hi[axis].max(p[axis]);
        }
    }
    let span = (hi[0] - lo[0])
        .max(hi[1] - lo[1])
        .max(hi[2] - lo[2])
        .max(1.0);

    check_duplicates_3d(points, lo, span)?;

    let n = points.len();
    let mut order: Vec<usize> = (0..n).collect();
    let inv_span = 1.0 / span;
    order.sort_unstable_by_key(|&i| morton_key(points[i], lo, inv_span));

    let mut pts = Vec::with_capacity(n + 4);
    pts.extend_from_slice(points);

    let centre = [
        f64::midpoint(lo[0], hi[0]),
        f64::midpoint(lo[1], hi[1]),
        f64::midpoint(lo[2], hi[2]),
    ];

    let s = 100.0 * span;
    pts.push([centre[0] - 3.0 * s, centre[1] - s, centre[2] - s]);
    pts.push([centre[0] + 3.0 * s, centre[1] - s, centre[2] - s]);
    pts.push([centre[0], centre[1] + 3.0 * s, centre[2] - s]);
    pts.push([centre[0], centre[1], centre[2] + 3.0 * s]);

    let super_tet = positive(&pts, [n, n + 1, n + 2, n + 3]);
    let mut tets: Vec<[usize; 4]> = Vec::with_capacity(6 * n + 16);
    let mut adj: Vec<[usize; 4]> = Vec::with_capacity(6 * n + 16);
    let mut alive: Vec<bool> = Vec::with_capacity(6 * n + 16);
    let mut in_cavity: Vec<u32> = Vec::with_capacity(6 * n + 16);
    let mut free_list: Vec<usize> = Vec::new();

    tets.push(super_tet);
    adj.push([NO_NEIGHBOR; 4]);
    alive.push(true);
    in_cavity.push(0);

    let mut stamp: u32 = 0;
    let mut last_tet: usize = 0;
    let mut cavity: Vec<usize> = Vec::new();
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut boundary: Vec<([usize; 3], usize, usize)> = Vec::new();
    let mut edge_map: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

    for &p_idx in &order {
        let p = pts[p_idx];

        let Some(seed) = locate_containing_tet(&pts, &tets, &adj, &alive, last_tet, p_idx) else {
            continue;
        };

        stamp = stamp.wrapping_add(1);
        if stamp == 0 {
            in_cavity.fill(0);
            stamp = 1;
        }

        // BFS over tet-to-tet adjacency to discover the Bowyer-Watson cavity.
        cavity.clear();
        queue.clear();
        in_cavity[seed] = stamp;
        cavity.push(seed);
        queue.push_back(seed);

        while let Some(cur) = queue.pop_front() {
            for f in 0..4 {
                let nb = adj[cur][f];
                if nb == NO_NEIGHBOR || !alive[nb] || in_cavity[nb] == stamp {
                    continue;
                }
                let o = tets[nb];
                if insphere_sos(&pts, o[0], o[1], o[2], o[3], p_idx) {
                    in_cavity[nb] = stamp;
                    cavity.push(nb);
                    queue.push_back(nb);
                }
            }
        }

        // Expand the cavity across any coplanar/non-visible boundary face so
        // every newly coned tetrahedron has strictly positive volume.
        let mut cav_idx = 0_usize;
        while cav_idx < cavity.len() {
            let t = cavity[cav_idx];
            cav_idx += 1;
            let v = tets[t];
            for f in 0..4 {
                let nb = adj[t][f];
                if nb != NO_NEIGHBOR && alive[nb] && in_cavity[nb] == stamp {
                    continue;
                }
                let face = tet_face(v, f);
                if orient3d(pts[face[0]], pts[face[1]], pts[face[2]], p) <= 0.0
                    && nb != NO_NEIGHBOR
                    && alive[nb]
                {
                    in_cavity[nb] = stamp;
                    cavity.push(nb);
                }
            }
        }

        // Collect cavity boundary faces and their exterior neighbor slots.
        boundary.clear();
        for &t in &cavity {
            let v = tets[t];
            for f in 0..4 {
                let nb = adj[t][f];
                if nb == NO_NEIGHBOR || !alive[nb] || in_cavity[nb] != stamp {
                    let ext_f = if nb != NO_NEIGHBOR {
                        adj[nb].iter().position(|&x| x == t).unwrap_or(0)
                    } else {
                        0
                    };
                    boundary.push((tet_face(v, f), nb, ext_f));
                }
            }
        }

        if boundary.is_empty() {
            continue;
        }

        // Retire cavity tetrahedra.
        for &t in &cavity {
            alive[t] = false;
            free_list.push(t);
        }

        // Create new tetrahedra coned to `p_idx` and wire up adjacency.
        edge_map.clear();
        for &(face, ext_nb, ext_f) in &boundary {
            let new_t = if let Some(slot) = free_list.pop() {
                tets[slot] = [face[0], face[1], face[2], p_idx];
                adj[slot] = [NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR, ext_nb];
                alive[slot] = true;
                slot
            } else {
                let slot = tets.len();
                tets.push([face[0], face[1], face[2], p_idx]);
                adj.push([NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR, ext_nb]);
                alive.push(true);
                in_cavity.push(0);
                slot
            };
            last_tet = new_t;

            if ext_nb != NO_NEIGHBOR {
                adj[ext_nb][ext_f] = new_t;
            }

            let internal_edges = [
                (0_usize, face[1], face[2]),
                (1_usize, face[0], face[2]),
                (2_usize, face[0], face[1]),
            ];
            for (f_idx, u, v) in internal_edges {
                let key = if u <= v { (u, v) } else { (v, u) };
                if let Some((other_t, other_f)) = edge_map.remove(&key) {
                    adj[new_t][f_idx] = other_t;
                    adj[other_t][other_f] = new_t;
                } else {
                    edge_map.insert(key, (new_t, f_idx));
                }
            }
        }
    }

    let mut result: Vec<[usize; 4]> = Vec::new();
    for (idx, &t) in tets.iter().enumerate() {
        if alive[idx] && t[0] < n && t[1] < n && t[2] < n && t[3] < n {
            result.push(t);
        }
    }

    Ok(result)
}

/// Checks for duplicate or near-coincident 3D points in `O(n)` expected time
/// using a uniform spatial hash grid.
fn check_duplicates_3d(
    points: &[[f64; 3]],
    lo: [f64; 3],
    span: f64,
) -> Result<(), String> {
    let tol = 1e-12 * span;
    let inv_cell = 1.0 / tol;
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::with_capacity(points.len());

    for (i, p) in points.iter().enumerate() {
        let cx = ((p[0] - lo[0]) * inv_cell).floor() as i64;
        let cy = ((p[1] - lo[1]) * inv_cell).floor() as i64;
        let cz = ((p[2] - lo[2]) * inv_cell).floor() as i64;

        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    if let Some(bucket) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &j in bucket {
                            let q = points[j];
                            if (q[0] - p[0]).abs() <= tol
                                && (q[1] - p[1]).abs() <= tol
                                && (q[2] - p[2]).abs() <= tol
                            {
                                return Err("duplicate points".to_string());
                            }
                        }
                    }
                }
            }
        }
        grid.entry((cx, cy, cz)).or_default().push(i);
    }

    Ok(())
}

/// Spreads the lowest 10 bits of `v` so there are two zero bits between each
/// bit, used for 3D Morton (Z-order) curve keys.
const fn spread_bits_3d(mut v: u32) -> u32 {
    v &= 0x0000_03ff;
    v = (v | (v << 16)) & 0x0300_00ff;
    v = (v | (v << 8)) & 0x0300_f00f;
    v = (v | (v << 4)) & 0x030c_30c3;
    v = (v | (v << 2)) & 0x0924_9249;
    v
}

/// 30-bit Morton (Z-order) key for spatial sorting in `[lo, lo + span]^3`.
fn morton_key(
    p: [f64; 3],
    lo: [f64; 3],
    inv_span: f64,
) -> u32 {
    let q = |x: f64, l: f64| -> u32 { (((x - l) * inv_span).clamp(0.0, 1.0) * 1023.0) as u32 };
    spread_bits_3d(q(p[0], lo[0]))
        | (spread_bits_3d(q(p[1], lo[1])) << 1)
        | (spread_bits_3d(q(p[2], lo[2])) << 2)
}

/// Returns face `f` of the positively oriented tetrahedron `v` (the face
/// opposite `v[f]`), oriented so that `orient3d(face[0], face[1], face[2], v[f]) > 0`.
const fn tet_face(
    v: [usize; 4],
    f: usize,
) -> [usize; 3] {
    match f {
        | 0 => [v[1], v[3], v[2]],
        | 1 => [v[0], v[2], v[3]],
        | 2 => [v[0], v[3], v[1]],
        | _ => [v[0], v[1], v[2]],
    }
}

/// Locates a live tetrahedron whose circumsphere contains `pts[p_idx]` via a
/// visibility walk starting from `start`, falling back to a linear scan if
/// the walk encounters a cycle or boundary.
fn locate_containing_tet(
    pts: &[[f64; 3]],
    tets: &[[usize; 4]],
    adj: &[[usize; 4]],
    alive: &[bool],
    start: usize,
    p_idx: usize,
) -> Option<usize> {
    let p = pts[p_idx];
    let mut cur = if start < alive.len() && alive[start] {
        start
    } else {
        alive.iter().position(|&a| a)?
    };

    let max_steps = tets.len() + 64;
    for step in 0..max_steps {
        let v = tets[cur];
        let mut stepped = false;
        for k in 0..4 {
            let f = (step + k) & 3;
            let face = tet_face(v, f);
            if orient3d(pts[face[0]], pts[face[1]], pts[face[2]], p) < 0.0 {
                let nb = adj[cur][f];
                if nb != NO_NEIGHBOR && alive[nb] {
                    cur = nb;
                    stepped = true;
                    break;
                }
            }
        }
        if !stepped {
            let v = tets[cur];
            if insphere_sos(pts, v[0], v[1], v[2], v[3], p_idx) {
                return Some(cur);
            }
            for &nb in &adj[cur] {
                if nb != NO_NEIGHBOR && alive[nb] {
                    let o = tets[nb];
                    if insphere_sos(pts, o[0], o[1], o[2], o[3], p_idx) {
                        return Some(nb);
                    }
                }
            }
            break;
        }
    }

    // Fallback linear scan if the visibility walk stalled.
    for (idx, &v) in tets.iter().enumerate() {
        if alive[idx] && insphere_sos(pts, v[0], v[1], v[2], v[3], p_idx) {
            return Some(idx);
        }
    }
    None
}

/// The volume of every tetrahedron in `tets`, always non-negative.
///
/// Tetrahedra whose indices do not address `points` are reported as `0.0`.
#[must_use]
pub fn tetra_volumes(
    points: &[[f64; 3]],
    tets: &[[usize; 4]],
) -> Vec<f64> {
    let n = points.len();
    tets.iter()
        .map(|t| {
            if t[0] < n && t[1] < n && t[2] < n && t[3] < n {
                orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]).abs() / 6.0
            } else {
                0.0
            }
        })
        .collect()
}

/// The faces of `tets` that belong to exactly one tetrahedron.
///
/// Each face is returned once, with the winding it has in the tetrahedron that
/// owns it, so the faces of a closed mesh point outwards.
#[must_use]
pub fn boundary_faces(tets: &[[usize; 4]]) -> Vec<[usize; 3]> {
    let mut order: Vec<([usize; 3], [usize; 3], usize)> = Vec::new();
    let mut index: HashMap<[usize; 3], usize> = HashMap::new();

    for t in tets {
        for face in faces(*t) {
            let key = face_key(face);
            match index.get(&key) {
                | Some(&slot) => order[slot].2 += 1,
                | None => {
                    index.insert(key, order.len());
                    order.push((key, face, 1));
                },
            }
        }
    }

    order
        .into_iter()
        .filter(|entry| entry.2 == 1)
        .map(|entry| entry.1)
        .collect()
}

/// The four outward faces of the positively oriented tetrahedron `t`.
fn faces(t: [usize; 4]) -> [[usize; 3]; 4] {
    [
        [t[1], t[2], t[3]],
        [t[0], t[3], t[2]],
        [t[0], t[1], t[3]],
        [t[2], t[1], t[0]],
    ]
}

/// The sorted key of a face, used to pair faces across tetrahedra.
fn face_key(face: [usize; 3]) -> [usize; 3] {
    let mut key = face;
    key.sort_unstable();
    key
}

/// Reorders `t` so that [`orient3d`] reports a positive volume.
fn positive(
    points: &[[f64; 3]],
    mut t: [usize; 4],
) -> [usize; 4] {
    if orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]) < 0.0 {
        t.swap(2, 3);
    }
    t
}

#[cfg(test)]
mod tests {
    use super::boundary_faces;
    use super::delaunay_3d;
    use super::tetra_volumes;
    use crate::kernels::meshgen::util::insphere;
    use crate::kernels::meshgen::util::orient3d;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    fn cube() -> Vec<[f64; 3]> {
        let mut points = Vec::with_capacity(8);
        for &z in &[0.0, 1.0] {
            for &y in &[0.0, 1.0] {
                for &x in &[0.0, 1.0] {
                    points.push([x, y, z]);
                }
            }
        }
        points
    }

    #[test]
    fn unit_cube_fills_a_unit_volume() {
        let points = cube();
        let tets = expect_ok(delaunay_3d(&points));

        assert!(
            tets.len() == 5 || tets.len() == 6,
            "expected 5 or 6 tetrahedra, got {}",
            tets.len()
        );

        for t in &tets {
            assert!(orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]) > 0.0);
        }

        let total: f64 = tetra_volumes(&points, &tets).iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "total volume {total}");

        let faces = boundary_faces(&tets);
        assert_ne!(faces.len(), 0);

        let mut edges: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        for f in &faces {
            for k in 0..3 {
                let a = f[k];
                let b = f[(k + 1) % 3];
                let key = if a <= b { (a, b) } else { (b, a) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
        for (edge, count) in edges {
            assert_eq!(count, 2, "boundary edge {edge:?} appears {count} times");
        }
    }

    #[test]
    fn single_tetrahedron_has_four_boundary_faces() {
        let tets = [[0, 1, 2, 3]];
        let faces = boundary_faces(&tets);
        assert_eq!(faces.len(), 4);
        let mut keys: Vec<[usize; 3]> = faces
            .iter()
            .map(|f| {
                let mut k = *f;
                k.sort_unstable();
                k
            })
            .collect();
        keys.sort_unstable();
        assert_eq!(keys, vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]);
    }

    #[test]
    fn interior_face_is_not_a_boundary_face() {
        let tets = [[0, 1, 2, 3], [1, 2, 3, 4]];
        let faces = boundary_faces(&tets);
        assert_eq!(faces.len(), 6);
    }

    #[test]
    fn tetra_volumes_of_a_reference_simplex() {
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let volumes = tetra_volumes(&points, &[[0, 1, 2, 3]]);
        assert!((volumes[0] - 1.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(delaunay_3d(&[[0.0, 0.0, 0.0]; 3]).is_err());
        assert!(
            delaunay_3d(&[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [f64::NAN, 0.0, 1.0]
            ])
            .is_err()
        );
        assert!(
            delaunay_3d(&[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ])
            .is_err()
        );
    }

    #[test]
    fn regular_cospherical_3x3x3_grid_fills_exact_volume() {
        let mut points = Vec::with_capacity(27);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    points.push([f64::from(x), f64::from(y), f64::from(z)]);
                }
            }
        }
        let tets = expect_ok(delaunay_3d(&points));
        assert_ne!(tets, [] as [[usize; 4]; 0]);

        for t in &tets {
            assert!(orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]) > 0.0);
        }
        let total: f64 = tetra_volumes(&points, &tets).iter().sum();
        assert!((total - 8.0).abs() < 1e-9, "total volume {total}");
    }

    #[test]
    fn large_3d_point_cloud_is_valid_and_delaunay() {
        let mut points = cube();
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next_coord = || -> f64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let u = ((state >> 11) as f64) / ((1_u64 << 53) as f64);
            0.01 + 0.98 * u
        };
        for _ in 0..500 {
            points.push([next_coord(), next_coord(), next_coord()]);
        }

        let tets = expect_ok(delaunay_3d(&points));
        assert!(tets.len() > 500, "expected >500 tets, got {}", tets.len());

        for t in &tets {
            assert!(orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]) > 0.0);
        }

        let total: f64 = tetra_volumes(&points, &tets).iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "total volume {total}");

        // Verify the Delaunay empty-circumsphere property on a broad sample.
        let stride = (tets.len() / 200).max(1);
        for t in tets.iter().step_by(stride) {
            for (idx, &p) in points.iter().enumerate() {
                if idx == t[0] || idx == t[1] || idx == t[2] || idx == t[3] {
                    continue;
                }
                let val = insphere(points[t[0]], points[t[1]], points[t[2]], points[t[3]], p);
                assert!(
                    val <= 0.0,
                    "point {idx} strictly inside circumsphere of {t:?}: {val}"
                );
            }
        }
    }
}
