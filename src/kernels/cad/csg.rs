//! Boolean operations on closed triangle meshes by BSP-tree CSG.
//!
//! This is a direct port of the classic `csg.js` algorithm: each operand is
//! turned into a BSP tree, the trees clip each other's polygons, and the
//! surviving polygons are recombined. The result is a polygon soup that is
//! welded and fan-triangulated into a [`Mesh`].
//!
//! # Limitations
//!
//! * The operands must be closed and consistently oriented; an open or
//!   self-intersecting input gives an undefined result rather than an error.
//! * Splitting is recursive, so a mesh with many nearly coincident planes can
//!   exhaust the stack.
//! * The output keeps the T-junctions introduced by splitting, so it is closed
//!   but not necessarily edge-manifold until [`super::repair::weld`] is run
//!   with a tolerance.

use super::repair::weld;
use super::types::Mesh;
use super::types::cross;
use super::types::dot;
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

/// Fan-triangulates `polygons` and welds the result.
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

    let (welded, _) = weld(&Mesh { vertices, triangles }, 1e-9);
    welded
}

#[cfg(test)]
mod tests {
    use super::BoolOp;
    use super::boolean;
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

        // Splitting leaves T-junctions on the coplanar faces, so the result is
        // closed geometrically but not edge-matched.
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
}
