//! The triangle-soup container shared by the CAD readers, repair and Boolean
//! operations, plus the small vector helpers they need.

/// A triangle mesh stored as an indexed vertex array and a triangle list.
///
/// The vertices carry no welding guarantee: an STL file, for instance, repeats
/// every vertex once per triangle. Use [`super::repair::weld`] to merge
/// coincident vertices before relying on adjacency.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Vertex coordinates, one `[x, y, z]` per vertex.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles as triples of indices into [`Mesh::vertices`].
    pub triangles: Vec<[usize; 3]>,
}

impl Mesh {
    /// Builds a mesh from raw vertex and triangle arrays.
    #[must_use]
    pub fn new(
        vertices: Vec<[f64; 3]>,
        triangles: Vec<[usize; 3]>,
    ) -> Self {
        Self { vertices, triangles }
    }

    /// The number of vertices.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// The number of triangles.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
}

/// The vector `b - a`.
pub(super) fn sub(
    a: [f64; 3],
    b: [f64; 3],
) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The cross product `a x b`.
pub(super) fn cross(
    a: [f64; 3],
    b: [f64; 3],
) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The dot product `a . b`.
pub(super) fn dot(
    a: [f64; 3],
    b: [f64; 3],
) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The Euclidean length of `a`.
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// `a` scaled to unit length, or the zero vector when `a` is degenerate.
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len <= 0.0 {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / len, a[1] / len, a[2] / len]
    }
}

/// The normalised key of the undirected edge between `a` and `b`.
pub(super) fn edge_key(
    a: usize,
    b: usize,
) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

/// The unit normal of triangle `t`, or [`None`] when it is degenerate.
pub(super) fn triangle_normal(
    mesh: &Mesh,
    t: [usize; 3],
) -> Option<[f64; 3]> {
    if t[0] >= mesh.vertices.len() || t[1] >= mesh.vertices.len() || t[2] >= mesh.vertices.len() {
        return None;
    }

    let a = mesh.vertices[t[0]];
    let b = mesh.vertices[t[1]];
    let c = mesh.vertices[t[2]];
    let n = cross(sub(b, a), sub(c, a));
    if length(n) <= 0.0 {
        None
    } else {
        Some(normalize(n))
    }
}

/// Whether every triangle index addresses a vertex of `mesh`.
pub(crate) fn indices_valid(mesh: &Mesh) -> bool {
    let n = mesh.vertices.len();
    mesh.triangles
        .iter()
        .all(|t| t[0] < n && t[1] < n && t[2] < n)
}

#[cfg(test)]
mod tests {
    use super::Mesh;
    use super::cross;
    use super::dot;
    use super::normalize;
    use super::triangle_normal;

    #[test]
    fn vector_helpers_agree_with_the_hand_computation() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_eq!(cross(a, b), [0.0, 0.0, 1.0]);
        assert!((dot(a, b)).abs() < 1e-15);
        assert!((normalize([3.0, 0.0, 0.0])[0] - 1.0).abs() < 1e-15);
        assert_eq!(normalize([0.0, 0.0, 0.0]), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn triangle_normal_of_a_degenerate_triangle_is_none() {
        let mesh = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0, 1, 2]],
        );
        assert_eq!(triangle_normal(&mesh, [0, 1, 2]), None);
        assert_eq!(mesh.vertex_count(), 3);
        assert_eq!(mesh.triangle_count(), 1);
    }

    #[test]
    fn triangle_normal_of_an_xy_triangle_points_along_z() {
        let mesh = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        assert_eq!(triangle_normal(&mesh, [0, 1, 2]), Some([0.0, 0.0, 1.0]));
    }
}
