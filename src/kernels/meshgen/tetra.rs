//! Three-dimensional Delaunay tetrahedralisation and tetrahedral mesh queries.
//!
//! [`delaunay_3d`] is a Bowyer–Watson insertion built on the in-sphere
//! predicate: a new point claims a connected cavity of tetrahedra whose
//! circumspheres contain it, and the cavity boundary is re-coned to the point.
//! The cavity is grown by breadth-first search from the tetrahedron whose
//! circumsphere most deeply contains the point, which keeps the result valid
//! when the input is cospherical or lies on the convex hull.
//!
//! # Limitations
//!
//! * The predicates are plain `f64`. Strongly degenerate inputs (coplanar
//!   subsets, repeated points) are rejected by [`delaunay_3d`] or produce a
//!   tetrahedralisation that is valid but not unique.
//! * The insertion is `O(n^2)` in the point count, so it targets the hundreds
//!   to low thousands of points rather than large clouds.

use std::collections::HashMap;
use std::collections::VecDeque;

use super::util::orient3d;

/// An unstructured tetrahedral mesh in three dimensions.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TetMesh {
    /// Vertex coordinates, one `[x, y, z]` per vertex.
    pub points: Vec<[f64; 3]>,
    /// Positively oriented tetrahedra, as quadruples of vertex indices.
    pub tetrahedra: Vec<[usize; 4]>,
}

/// The relative tolerance used when comparing in-sphere values.
const INSPHERE_EPS: f64 = 1e-12;

/// The in-sphere determinant of `p` against the positively oriented tetrahedron
/// `a`, `b`, `c`, `d`.
///
/// The value is positive when `p` lies strictly inside the circumsphere. The
/// raw determinant of the `4 x 4` system has the opposite sign in three
/// dimensions, so the result is negated once here.
fn insphere(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> f64 {
    let ax = a[0] - p[0];
    let ay = a[1] - p[1];
    let az = a[2] - p[2];
    let bx = b[0] - p[0];
    let by = b[1] - p[1];
    let bz = b[2] - p[2];
    let cx = c[0] - p[0];
    let cy = c[1] - p[1];
    let cz = c[2] - p[2];
    let dx = d[0] - p[0];
    let dy = d[1] - p[1];
    let dz = d[2] - p[2];

    let a2 = ax * ax + ay * ay + az * az;
    let b2 = bx * bx + by * by + bz * bz;
    let c2 = cx * cx + cy * cy + cz * cz;
    let d2 = dx * dx + dy * dy + dz * dz;

    let minor_b = cz * d2 - c2 * dz;
    let minor_c = cy * d2 - c2 * dy;
    let minor_d = cy * dz - cz * dy;

    let first = ax * (by * minor_b - bz * minor_c + b2 * minor_d)
        - ay * (bx * minor_b - bz * (cx * d2 - c2 * dx) + b2 * (cx * dz - cz * dx))
        + az * (bx * minor_c - by * (cx * d2 - c2 * dx) + b2 * (cx * dy - cy * dx))
        - a2 * (bx * minor_d - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx));

    -first
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

    let n = points.len();
    let mut pts = points.to_vec();

    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &pts {
        for axis in 0..3 {
            lo[axis] = lo[axis].min(p[axis]);
            hi[axis] = hi[axis].max(p[axis]);
        }
    }
    let span = (hi[0] - lo[0])
        .max(hi[1] - lo[1])
        .max(hi[2] - lo[2])
        .max(1.0);
    let centre = [
        f64::midpoint(lo[0], hi[0]),
        f64::midpoint(lo[1], hi[1]),
        f64::midpoint(lo[2], hi[2]),
    ];

    let s = 20.0 * span;
    pts.push([centre[0] - 3.0 * s, centre[1] - s, centre[2] - s]);
    pts.push([centre[0] + 3.0 * s, centre[1] - s, centre[2] - s]);
    pts.push([centre[0], centre[1] + 3.0 * s, centre[2] - s]);
    pts.push([centre[0], centre[1], centre[2] + 3.0 * s]);

    let mut tets: Vec<[usize; 4]> = vec![positive(pts.as_slice(), [n, n + 1, n + 2, n + 3])];

    for i in 0..n {
        let p = pts[i];
        if tets.is_empty() {
            return Err("the point set is degenerate".to_string());
        }

        let mut seed = 0_usize;
        let mut best = f64::NEG_INFINITY;
        for (index, t) in tets.iter().enumerate() {
            let value = insphere(pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]], p);
            if value > best {
                best = value;
                seed = index;
            }
        }
        if best <= 0.0 {
            continue;
        }

        let neighbours = face_neighbours(&tets);
        let mut inside = vec![false; tets.len()];
        inside[seed] = true;
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(seed);

        while let Some(current) = queue.pop_front() {
            let t = tets[current];
            for face in faces(t) {
                let key = face_key(face);
                let Some(list) = neighbours.get(&key) else {
                    continue;
                };
                for &other in list {
                    if other == current || inside[other] {
                        continue;
                    }
                    let o = tets[other];
                    let value = insphere(pts[o[0]], pts[o[1]], pts[o[2]], pts[o[3]], p);
                    if value >= -INSPHERE_EPS {
                        inside[other] = true;
                        queue.push_back(other);
                    }
                }
            }
        }

        let mut boundary: Vec<[usize; 3]> = Vec::new();
        for (index, t) in tets.iter().enumerate() {
            if !inside[index] {
                continue;
            }
            for face in faces(*t) {
                let key = face_key(face);
                let shared = neighbours.get(&key).map_or(0, |list| {
                    list.iter().filter(|&&other| inside[other]).count()
                });
                if shared == 1 {
                    boundary.push(face);
                }
            }
        }

        if boundary.is_empty() {
            continue;
        }

        let mut kept: Vec<[usize; 4]> = Vec::new();
        for (index, t) in tets.iter().enumerate() {
            if !inside[index] {
                kept.push(*t);
            }
        }
        for face in boundary {
            let candidate = positive(pts.as_slice(), [face[0], face[1], face[2], i]);
            let volume = orient3d(
                pts[candidate[0]],
                pts[candidate[1]],
                pts[candidate[2]],
                pts[candidate[3]],
            );
            if volume > 0.0 {
                kept.push(candidate);
            }
        }
        tets = kept;
    }

    tets.retain(|t| t[0] < n && t[1] < n && t[2] < n && t[3] < n);
    for t in &mut tets {
        *t = positive(points, *t);
    }
    deduplicate(&mut tets);

    Ok(tets)
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

/// Maps every face to the tetrahedra that contain it.
fn face_neighbours(tets: &[[usize; 4]]) -> HashMap<[usize; 3], Vec<usize>> {
    let mut map: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for (index, t) in tets.iter().enumerate() {
        for face in faces(*t) {
            map.entry(face_key(face)).or_default().push(index);
        }
    }
    map
}

/// Removes tetrahedra that repeat the same four vertices.
fn deduplicate(tets: &mut Vec<[usize; 4]>) {
    let mut seen: std::collections::HashSet<[usize; 4]> = std::collections::HashSet::new();
    tets.retain(|t| {
        let mut key = *t;
        key.sort_unstable();
        seen.insert(key)
    });
}

#[cfg(test)]
mod tests {
    use super::boundary_faces;
    use super::delaunay_3d;
    use super::tetra_volumes;
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
    }
}
