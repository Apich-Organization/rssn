//! Boolean operations on closed triangle meshes by BSP-tree CSG.
//!
//! This is a direct port of the classic `csg.js` algorithm: each operand is
//! turned into a BSP tree, the trees clip each other's polygons, and the
//! surviving polygons are recombined. The result is fan-triangulated, welded,
//! and passed through a conforming edge-splitting pass that eliminates
//! T-junctions introduced by BSP plane splitting so the output [`Mesh`] is
//! watertight and edge-manifold.
//!
//! # Limitations
//!
//! * The operands must be closed and consistently oriented; an open or
//!   self-intersecting input gives an undefined result rather than an error.
//! * Splitting is recursive, so a mesh with many nearly coincident planes can
//!   exhaust the stack.

use std::collections::HashMap;

use super::repair::weld;
use super::types::Mesh;
use super::types::cross;
use super::types::dot;
use super::types::edge_key;
use super::types::length;
use super::types::sub;

/// The distance below which a vertex counts as lying on a splitting plane.
const EPS: f64 = 1e-9;

/// The Boolean operation [`boolean`] performs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolOp {
    /// Everything inside `a` or `b`.
    Union,
    /// Everything inside `a` but outside `b`.
    Difference,
    /// Everything inside both `a` and `b`.
    Intersection,
}

/// An oriented plane `normal . x = w`.
#[derive(Clone, Debug)]
struct Plane {
    normal: [f64; 3],
    w: f64,
}

impl Plane {
    /// Builds the plane through three points, or [`None`] when they are collinear.
    fn from_points(
        a: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
    ) -> Option<Self> {
        let n = cross(sub(b, a), sub(c, a));
        let len = length(n);
        if len <= 1e-14 {
            return None;
        }
        let normal = [n[0] / len, n[1] / len, n[2] / len];
        Some(Self {
            normal,
            w: dot(normal, a),
        })
    }

    /// Reverses the plane's orientation.
    fn flip(&mut self) {
        self.normal = [-self.normal[0], -self.normal[1], -self.normal[2]];
        self.w = -self.w;
    }

    /// Sorts `poly` into the four buckets relative to this plane.
    fn split(
        &self,
        poly: &Polygon,
        coplanar_front: &mut Vec<Polygon>,
        coplanar_back: &mut Vec<Polygon>,
        front: &mut Vec<Polygon>,
        back: &mut Vec<Polygon>,
    ) {
        const COPLANAR: u8 = 0;
        const FRONT: u8 = 1;
        const BACK: u8 = 2;
        const SPANNING: u8 = 3;

        let mut polygon_type: u8 = 0;
        let mut types: Vec<u8> = Vec::with_capacity(poly.vertices.len());
        for v in &poly.vertices {
            let t = dot(self.normal, *v) - self.w;
            let kind = if t < -EPS {
                BACK
            } else if t > EPS {
                FRONT
            } else {
                COPLANAR
            };
            polygon_type |= kind;
            types.push(kind);
        }

        match polygon_type {
            | COPLANAR => {
                if dot(self.normal, poly.plane.normal) > 0.0 {
                    coplanar_front.push(poly.clone());
                } else {
                    coplanar_back.push(poly.clone());
                }
            },
            | FRONT => front.push(poly.clone()),
            | BACK => back.push(poly.clone()),
            | _ => {
                let mut f: Vec<[f64; 3]> = Vec::new();
                let mut b: Vec<[f64; 3]> = Vec::new();
                let n = poly.vertices.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    let ti = types[i];
                    let tj = types[j];
                    let vi = poly.vertices[i];
                    let vj = poly.vertices[j];

                    if ti != BACK {
                        f.push(vi);
                    }
                    if ti != FRONT {
                        b.push(vi);
                    }
                    if (ti | tj) == SPANNING {
                        let denom = dot(self.normal, sub(vj, vi));
                        if denom != 0.0 {
                            let t = (self.w - dot(self.normal, vi)) / denom;
                            let v = [
                                vi[0] + t * (vj[0] - vi[0]),
                                vi[1] + t * (vj[1] - vi[1]),
                                vi[2] + t * (vj[2] - vi[2]),
                            ];
                            f.push(v);
                            b.push(v);
                        }
                    }
                }
                if f.len() >= 3 {
                    if let Some(p) = Polygon::new(f) {
                        front.push(p);
                    }
                }
                if b.len() >= 3 {
                    if let Some(p) = Polygon::new(b) {
                        back.push(p);
                    }
                }
            },
        }
    }
}

/// A planar convex polygon with its supporting plane.
#[derive(Clone, Debug)]
struct Polygon {
    vertices: Vec<[f64; 3]>,
    plane: Plane,
}

impl Polygon {
    /// Builds a polygon from its vertices, or [`None`] when it is degenerate.
    fn new(vertices: Vec<[f64; 3]>) -> Option<Self> {
        let plane = Polygon::plane_of(&vertices)?;
        Some(Self { vertices, plane })
    }

    /// The supporting plane of `vertices`, found from its first non-degenerate
    /// triple.
    fn plane_of(vertices: &[[f64; 3]]) -> Option<Plane> {
        if vertices.len() < 3 {
            return None;
        }
        for i in 1..vertices.len() - 1 {
            if let Some(plane) = Plane::from_points(vertices[0], vertices[i], vertices[i + 1]) {
                return Some(plane);
            }
        }
        None
    }

    /// Reverses the polygon's winding and orientation.
    fn flip(&mut self) {
        self.vertices.reverse();
        self.plane.flip();
    }
}

/// A node of a BSP tree.
#[derive(Clone, Debug, Default)]
struct Node {
    plane: Option<Plane>,
    front: Option<Box<Node>>,
    back: Option<Box<Node>>,
    polygons: Vec<Polygon>,
}

impl Node {
    /// Builds a tree from `polygons`.
    fn from_polygons(polygons: Vec<Polygon>) -> Self {
        let mut node = Node::default();
        node.build(polygons);
        node
    }

    /// Reverses the solid represented by the tree.
    fn invert(&mut self) {
        for poly in &mut self.polygons {
            poly.flip();
        }
        if let Some(plane) = &mut self.plane {
            plane.flip();
        }
        if let Some(front) = &mut self.front {
            front.invert();
        }
        if let Some(back) = &mut self.back {
            back.invert();
        }
        std::mem::swap(&mut self.front, &mut self.back);
    }

    /// Returns the parts of `polygons` that lie outside the solid.
    fn clip_polygons(
        &self,
        polygons: Vec<Polygon>,
    ) -> Vec<Polygon> {
        let Some(plane) = &self.plane else {
            return polygons;
        };

        let mut front: Vec<Polygon> = Vec::new();
        let mut back: Vec<Polygon> = Vec::new();
        for poly in &polygons {
            let mut coplanar_front = Vec::new();
            let mut coplanar_back = Vec::new();
            let mut front_part = Vec::new();
            let mut back_part = Vec::new();
            plane.split(
                poly,
                &mut coplanar_front,
                &mut coplanar_back,
                &mut front_part,
                &mut back_part,
            );
            front.extend(coplanar_front);
            front.extend(front_part);
            back.extend(coplanar_back);
            back.extend(back_part);
        }

        let mut front = match &self.front {
            | Some(node) => node.clip_polygons(front),
            | None => front,
        };
        let back = match &self.back {
            | Some(node) => node.clip_polygons(back),
            | None => Vec::new(),
        };

        front.extend(back);
        front
    }

    /// Removes the parts of this solid that lie inside `other`.
    fn clip_to(
        &mut self,
        other: &Node,
    ) {
        let polygons = std::mem::take(&mut self.polygons);
        self.polygons = other.clip_polygons(polygons);
        if let Some(front) = &mut self.front {
            front.clip_to(other);
        }
        if let Some(back) = &mut self.back {
            back.clip_to(other);
        }
    }

    /// Every polygon held anywhere in the tree.
    fn all_polygons(&self) -> Vec<Polygon> {
        let mut out = self.polygons.clone();
        if let Some(front) = &self.front {
            out.extend(front.all_polygons());
        }
        if let Some(back) = &self.back {
            out.extend(back.all_polygons());
        }
        out
    }

    /// Adds `polygons` to the tree.
    fn build(
        &mut self,
        polygons: Vec<Polygon>,
    ) {
        if polygons.is_empty() {
            return;
        }
        if self.plane.is_none() {
            self.plane = Some(polygons[0].plane.clone());
        }
        let Some(plane) = self.plane.clone() else {
            return;
        };

        let mut front: Vec<Polygon> = Vec::new();
        let mut back: Vec<Polygon> = Vec::new();

        for poly in &polygons {
            let mut coplanar_front = Vec::new();
            let mut coplanar_back = Vec::new();
            let mut front_part = Vec::new();
            let mut back_part = Vec::new();
            plane.split(
                poly,
                &mut coplanar_front,
                &mut coplanar_back,
                &mut front_part,
                &mut back_part,
            );
            self.polygons.extend(coplanar_front);
            self.polygons.extend(coplanar_back);
            front.extend(front_part);
            back.extend(back_part);
        }

        if !front.is_empty() {
            if self.front.is_none() {
                self.front = Some(Box::new(Node::default()));
            }
            if let Some(node) = &mut self.front {
                node.build(front);
            }
        }
        if !back.is_empty() {
            if self.back.is_none() {
                self.back = Some(Box::new(Node::default()));
            }
            if let Some(node) = &mut self.back {
                node.build(back);
            }
        }
    }
}

/// Combines the closed meshes `a` and `b` with the operation `op`.
///
/// # Errors
///
/// Returns `Err` when a mesh has no non-degenerate triangle, so that no solid
/// can be built from it.
pub fn boolean(
    a: &Mesh,
    b: &Mesh,
    op: BoolOp,
) -> Result<Mesh, String> {
    let polygons_a = mesh_to_polygons(a)?;
    let polygons_b = mesh_to_polygons(b)?;

    let mut node_a = Node::from_polygons(polygons_a);
    let mut node_b = Node::from_polygons(polygons_b);

    match op {
        | BoolOp::Union => {
            node_a.clip_to(&node_b);
            node_b.clip_to(&node_a);
            node_b.invert();
            node_b.clip_to(&node_a);
            node_b.invert();
            node_a.build(node_b.all_polygons());
        },
        | BoolOp::Difference => {
            node_a.invert();
            node_a.clip_to(&node_b);
            node_b.clip_to(&node_a);
            node_b.invert();
            node_b.clip_to(&node_a);
            node_b.invert();
            node_a.build(node_b.all_polygons());
            node_a.invert();
        },
        | BoolOp::Intersection => {
            node_a.invert();
            node_b.clip_to(&node_a);
            node_b.invert();
            node_a.clip_to(&node_b);
            node_b.clip_to(&node_a);
            node_a.build(node_b.all_polygons());
            node_a.invert();
        },
    }

    let polygons = node_a.all_polygons();
    if polygons.is_empty() {
        return Ok(Mesh::default());
    }

    Ok(polygons_to_mesh(&polygons))
}

/// Converts the triangles of `mesh` into polygons.
fn mesh_to_polygons(mesh: &Mesh) -> Result<Vec<Polygon>, String> {
    let mut polygons = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        let (Some(a), Some(b), Some(c)) = (
            mesh.vertices.get(t[0]),
            mesh.vertices.get(t[1]),
            mesh.vertices.get(t[2]),
        ) else {
            continue;
        };
        if let Some(poly) = Polygon::new(vec![*a, *b, *c]) {
            polygons.push(poly);
        }
    }

    if polygons.is_empty() {
        return Err("mesh contains no non-degenerate triangles".to_string());
    }
    Ok(polygons)
}

/// Fan-triangulates `polygons`, welds coincident vertices, and eliminates
/// T-junctions so the resulting mesh is watertight and edge-manifold.
fn polygons_to_mesh(polygons: &[Polygon]) -> Mesh {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    for poly in polygons {
        if poly.vertices.len() < 3 {
            continue;
        }
        let base = vertices.len();
        for v in &poly.vertices {
            vertices.push(*v);
        }
        for i in 1..poly.vertices.len() - 1 {
            triangles.push([base, base + i, base + i + 1]);
        }
    }

    let (welded, _) = weld(&Mesh { vertices, triangles }, EPS);
    eliminate_t_junctions(welded)
}

/// Eliminates T-junctions by finding every vertex lying strictly in the
/// interior of any triangle edge and retriangulating affected triangles in a
/// conforming manner that preserves outward winding.
fn eliminate_t_junctions(mut mesh: Mesh) -> Mesh {
    const MAX_PASSES: usize = 8;

    if mesh.triangles.is_empty() || mesh.vertices.len() < 4 {
        return mesh;
    }

    let mut min_pt = [f64::INFINITY; 3];
    let mut max_pt = [f64::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for axis in 0..3 {
            min_pt[axis] = min_pt[axis].min(v[axis]);
            max_pt[axis] = max_pt[axis].max(v[axis]);
        }
    }
    let diag = length(sub(max_pt, min_pt));
    let tol = (diag * 1e-7).clamp(1e-9, 1e-6);
    let tol_sq = tol * tol;

    // Build a uniform spatial hash grid over all vertices once (splitting
    // never adds new vertices).
    let grid_dim = ((mesh.vertices.len() as f64).cbrt().ceil() as usize).clamp(4, 32);
    let cell_size = (diag / grid_dim as f64).max(tol * 16.0);
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (idx, v) in mesh.vertices.iter().enumerate() {
        let cell = (
            (v[0] / cell_size).floor() as i64,
            (v[1] / cell_size).floor() as i64,
            (v[2] / cell_size).floor() as i64,
        );
        grid.entry(cell).or_default().push(idx);
    }

    for _ in 0..MAX_PASSES {
        let mut edge_cache: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        let mut new_triangles: Vec<[usize; 3]> = Vec::with_capacity(mesh.triangles.len());
        let mut any_split = false;

        for &tri in &mesh.triangles {
            let [a, b, c] = tri;
            let s_ab = edge_interior_splits(
                a,
                b,
                &mesh.vertices,
                &grid,
                cell_size,
                tol,
                tol_sq,
                &mut edge_cache,
            );
            let s_bc = edge_interior_splits(
                b,
                c,
                &mesh.vertices,
                &grid,
                cell_size,
                tol,
                tol_sq,
                &mut edge_cache,
            );
            let s_ca = edge_interior_splits(
                c,
                a,
                &mesh.vertices,
                &grid,
                cell_size,
                tol,
                tol_sq,
                &mut edge_cache,
            );

            if s_ab.is_empty() && s_bc.is_empty() && s_ca.is_empty() {
                new_triangles.push(tri);
                continue;
            }

            any_split = true;
            let mut e0 = Vec::with_capacity(s_ab.len() + 2);
            e0.push(a);
            e0.extend_from_slice(&s_ab);
            e0.push(b);

            let mut e1 = Vec::with_capacity(s_bc.len() + 2);
            e1.push(b);
            e1.extend_from_slice(&s_bc);
            e1.push(c);

            let mut e2 = Vec::with_capacity(s_ca.len() + 2);
            e2.push(c);
            e2.extend_from_slice(&s_ca);
            e2.push(a);

            triangulate_subdivided_triangle(&e0, &e1, &e2, &mesh.vertices, &mut new_triangles);
        }

        mesh.triangles = new_triangles;
        if !any_split {
            break;
        }
    }

    mesh
}

/// Returns all vertex indices strictly in the interior of directed edge `u -> v`,
/// ordered by increasing parameter `t` from `u` to `v`.
#[allow(clippy::too_many_arguments)]
fn edge_interior_splits(
    u: usize,
    v: usize,
    vertices: &[[f64; 3]],
    grid: &HashMap<(i64, i64, i64), Vec<usize>>,
    cell_size: f64,
    tol: f64,
    tol_sq: f64,
    cache: &mut HashMap<(usize, usize), Vec<usize>>,
) -> Vec<usize> {
    let key = edge_key(u, v);
    let forward = u == key.0;

    let sorted = cache.entry(key).or_insert_with(|| {
        let p = vertices[key.0];
        let q = vertices[key.1];
        let d = sub(q, p);
        let len_sq = dot(d, d);
        if len_sq <= tol_sq {
            return Vec::new();
        }
        let seg_len = len_sq.sqrt();

        let min_x = p[0].min(q[0]) - tol;
        let max_x = p[0].max(q[0]) + tol;
        let min_y = p[1].min(q[1]) - tol;
        let max_y = p[1].max(q[1]) + tol;
        let min_z = p[2].min(q[2]) - tol;
        let max_z = p[2].max(q[2]) + tol;

        let cx0 = (min_x / cell_size).floor() as i64;
        let cx1 = (max_x / cell_size).floor() as i64;
        let cy0 = (min_y / cell_size).floor() as i64;
        let cy1 = (max_y / cell_size).floor() as i64;
        let cz0 = (min_z / cell_size).floor() as i64;
        let cz1 = (max_z / cell_size).floor() as i64;

        let cell_count =
            i128::from(cx1 - cx0 + 1) * i128::from(cy1 - cy0 + 1) * i128::from(cz1 - cz0 + 1);

        let mut hits: Vec<(f64, usize)> = Vec::new();
        let mut check_vertex = |w: usize| {
            if w == key.0 || w == key.1 {
                return;
            }
            let wp = vertices[w];
            if wp[0] < min_x
                || wp[0] > max_x
                || wp[1] < min_y
                || wp[1] > max_y
                || wp[2] < min_z
                || wp[2] > max_z
            {
                return;
            }
            let t = dot(sub(wp, p), d) / len_sq;
            if t * seg_len <= tol || (1.0 - t) * seg_len <= tol {
                return;
            }
            let proj = [p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]];
            let diff = sub(wp, proj);
            if dot(diff, diff) <= tol_sq {
                hits.push((t, w));
            }
        };

        if cell_count <= vertices.len() as i128 {
            for cx in cx0..=cx1 {
                for cy in cy0..=cy1 {
                    for cz in cz0..=cz1 {
                        if let Some(bucket) = grid.get(&(cx, cy, cz)) {
                            for &w in bucket {
                                check_vertex(w);
                            }
                        }
                    }
                }
            }
        } else {
            for w in 0..vertices.len() {
                check_vertex(w);
            }
        }

        hits.sort_by(|a, b| a.0.total_cmp(&b.0));
        hits.dedup_by_key(|item| item.1);
        hits.into_iter().map(|(_, w)| w).collect()
    });

    if forward {
        sorted.clone()
    } else {
        sorted.iter().copied().rev().collect()
    }
}

/// Conforming triangulation of a triangle `(e0[0], e1[0], e2[0])` whose three
/// directed boundary chains `e0`, `e1`, `e2` may carry interior split vertices.
fn triangulate_subdivided_triangle(
    e0: &[usize],
    e1: &[usize],
    e2: &[usize],
    vertices: &[[f64; 3]],
    out: &mut Vec<[usize; 3]>,
) {
    let n0 = e0.len().saturating_sub(2);
    let n1 = e1.len().saturating_sub(2);
    let n2 = e2.len().saturating_sub(2);
    let a = e0[0];
    let b = e1[0];
    let c = e2[0];

    let push_if_valid = |tri: [usize; 3], out: &mut Vec<[usize; 3]>| {
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0] {
            return;
        }
        let area2 = length(cross(
            sub(vertices[tri[1]], vertices[tri[0]]),
            sub(vertices[tri[2]], vertices[tri[0]]),
        ));
        if area2 > 1e-14 {
            out.push(tri);
        }
    };

    if n0 == 0 && n1 == 0 && n2 == 0 {
        push_if_valid([a, b, c], out);
        return;
    }

    if n1 == 0 && n2 == 0 {
        for i in 0..e0.len() - 1 {
            push_if_valid([e0[i], e0[i + 1], c], out);
        }
        return;
    }
    if n0 == 0 && n2 == 0 {
        for i in 0..e1.len() - 1 {
            push_if_valid([a, e1[i], e1[i + 1]], out);
        }
        return;
    }
    if n0 == 0 && n1 == 0 {
        for i in 0..e2.len() - 1 {
            push_if_valid([e2[i], e2[i + 1], b], out);
        }
        return;
    }

    // Two or three edges are subdivided: split the edge with the most interior
    // vertices at its midpoint and connect to the opposite corner.
    if n0 >= n1 && n0 >= n2 {
        let mid = e0.len() / 2;
        let m = e0[mid];
        let mc = [m, c];
        let cm = [c, m];
        triangulate_subdivided_triangle(&e0[..=mid], &mc, e2, vertices, out);
        triangulate_subdivided_triangle(&e0[mid..], e1, &cm, vertices, out);
    } else if n1 >= n0 && n1 >= n2 {
        let mid = e1.len() / 2;
        let m = e1[mid];
        let ma = [m, a];
        let am = [a, m];
        triangulate_subdivided_triangle(e0, &e1[..=mid], &ma, vertices, out);
        triangulate_subdivided_triangle(&am, &e1[mid..], e2, vertices, out);
    } else {
        let mid = e2.len() / 2;
        let m = e2[mid];
        let mb = [m, b];
        let bm = [b, m];
        triangulate_subdivided_triangle(e0, &bm, &e2[mid..], vertices, out);
        triangulate_subdivided_triangle(&mb, e1, &e2[..=mid], vertices, out);
    }
}

#[cfg(test)]
mod tests {
    use super::BoolOp;
    use super::boolean;
    use crate::kernels::cad::repair::is_closed;
    use crate::kernels::cad::repair::is_manifold;
    use crate::kernels::cad::repair::signed_volume;
    use crate::kernels::cad::types::Mesh;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    fn cube(
        offset: [f64; 3],
        size: f64,
    ) -> Mesh {
        let o = offset;
        let s = size;
        let vertices = vec![
            [o[0], o[1], o[2]],
            [o[0] + s, o[1], o[2]],
            [o[0] + s, o[1] + s, o[2]],
            [o[0], o[1] + s, o[2]],
            [o[0], o[1], o[2] + s],
            [o[0] + s, o[1], o[2] + s],
            [o[0] + s, o[1] + s, o[2] + s],
            [o[0], o[1] + s, o[2] + s],
        ];
        let triangles = vec![
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
        ];
        Mesh::new(vertices, triangles)
    }

    #[test]
    fn boolean_operations_on_offset_cubes() {
        let a = cube([0.0, 0.0, 0.0], 1.0);
        let b = cube([0.5, 0.0, 0.0], 1.0);

        let union = expect_ok(boolean(&a, &b, BoolOp::Union));
        assert!(
            (signed_volume(&union) - 1.5).abs() < 1e-6,
            "{}",
            signed_volume(&union)
        );

        let difference = expect_ok(boolean(&a, &b, BoolOp::Difference));
        assert!(
            (signed_volume(&difference) - 0.5).abs() < 1e-6,
            "{}",
            signed_volume(&difference)
        );

        let intersection = expect_ok(boolean(&a, &b, BoolOp::Intersection));
        assert!(
            (signed_volume(&intersection) - 0.5).abs() < 1e-6,
            "{}",
            signed_volume(&intersection)
        );
    }

    #[test]
    fn union_of_disjoint_cubes_adds_the_volumes() {
        let a = cube([0.0, 0.0, 0.0], 1.0);
        let b = cube([2.0, 0.0, 0.0], 1.0);
        let union = expect_ok(boolean(&a, &b, BoolOp::Union));
        assert!((signed_volume(&union) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn empty_operand_is_rejected() {
        let a = cube([0.0, 0.0, 0.0], 1.0);
        assert!(boolean(&a, &Mesh::default(), BoolOp::Union).is_err());
    }

    #[test]
    fn union_result_encloses_the_combined_box() {
        let a = cube([0.0, 0.0, 0.0], 1.0);
        let b = cube([0.5, 0.0, 0.0], 1.0);
        let union = expect_ok(boolean(&a, &b, BoolOp::Union));

        assert_ne!(union.triangles.len(), 0);
        assert!((signed_volume(&union) - 1.5).abs() < 1e-6);
        assert!(
            is_closed(&union),
            "union should be watertight after T-junction elimination"
        );
        assert!(is_manifold(&union), "union should be edge-manifold");

        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for v in &union.vertices {
            for axis in 0..3 {
                min[axis] = min[axis].min(v[axis]);
                max[axis] = max[axis].max(v[axis]);
            }
        }
        assert!((max[0] - min[0] - 1.5).abs() < 1e-9);
        assert!((max[1] - min[1] - 1.0).abs() < 1e-9);
        assert!((max[2] - min[2] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn boolean_operations_produce_watertight_manifold_meshes() {
        let a = cube([0.0, 0.0, 0.0], 1.0);
        let b = cube([0.5, 0.0, 0.0], 1.0);

        for (op, expected_vol) in [
            (BoolOp::Union, 1.5),
            (BoolOp::Difference, 0.5),
            (BoolOp::Intersection, 0.5),
        ] {
            let res = expect_ok(boolean(&a, &b, op));
            assert!(is_closed(&res), "{op:?} result is not closed");
            assert!(is_manifold(&res), "{op:?} result is not manifold");
            assert!(
                (signed_volume(&res) - expected_vol).abs() < 1e-6,
                "{op:?} volume was {}, expected {expected_vol}",
                signed_volume(&res)
            );
        }
    }

    #[test]
    fn three_axis_offset_cubes_are_watertight_and_accurate() {
        // Offsetting along all three axes creates T-junctions on every
        // intersecting face during BSP splitting.
        let a = cube([0.0, 0.0, 0.0], 1.0);
        let b = cube([0.5, 0.5, 0.5], 1.0);
        // Overlap is a 0.5 x 0.5 x 0.5 cube of volume 0.125.
        for (op, expected_vol) in [
            (BoolOp::Union, 1.875),
            (BoolOp::Difference, 0.875),
            (BoolOp::Intersection, 0.125),
        ] {
            let res = expect_ok(boolean(&a, &b, op));
            assert!(is_closed(&res), "{op:?} 3D-offset result is not closed");
            assert!(is_manifold(&res), "{op:?} 3D-offset result is not manifold");
            assert!(
                (signed_volume(&res) - expected_vol).abs() < 1e-6,
                "{op:?} 3D-offset volume was {}, expected {expected_vol}",
                signed_volume(&res)
            );
        }
    }
}
