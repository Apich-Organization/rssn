use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;
use crate::kernels::matrix::Matrix;
use crate::kernels::solve::LinearSolution;
use crate::kernels::solve::solve_linear_system;

/// A 2D vector.
#[derive(Clone, Copy, Default, Debug, Serialize, Deserialize)]
pub struct Vector2D {
    /// The x component of the vector.
    pub x: f64,
    /// The y component of the vector.
    pub y: f64,
}

impl Vector2D {
    /// Creates a new 2D vector.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
    ) -> Self {
        Self { x, y }
    }

    /// Calculates the norm of the vector.
    #[must_use]
    pub fn norm(&self) -> f64 {
        self.x.hypot(self.y)
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Sub for Vector2D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

/// A 3D vector.
#[derive(Clone, Copy, Default, Debug, Serialize, Deserialize)]
pub struct Vector3D {
    /// The x component of the vector.
    pub x: f64,
    /// The y component of the vector.
    pub y: f64,
    /// The z component of the vector.
    pub z: f64,
}

impl Vector3D {
    /// Creates a new 3D vector.
    #[allow(dead_code)]
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
        z: f64,
    ) -> Self {
        Self { x, y, z }
    }

    /// Calculates the norm of the vector.
    #[allow(dead_code)]
    #[must_use]
    pub fn norm(&self) -> f64 {
        self.z
            .mul_add(self.z, self.x.mul_add(self.x, self.y * self.y))
            .sqrt()
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

/// Specifies the type of boundary condition on an element.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum BoundaryCondition<T> {
    /// A known potential value.
    Potential(T),
    /// A known flux value.
    Flux(T),
}

/// A 2D boundary element.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Element2D {
    /// The first point of the element.
    pub p1: Vector2D,
    /// The second point of the element.
    pub p2: Vector2D,
    /// The midpoint of the element.
    pub midpoint: Vector2D,
    /// The length of the element.
    pub length: f64,
    /// The normal vector of the element.
    pub normal: Vector2D,
}

impl Element2D {
    /// Creates a new 2D boundary element.
    #[must_use]
    pub fn new(
        p1: Vector2D,
        p2: Vector2D,
    ) -> Self {
        let diff = p2 - p1;

        let length = diff.norm();

        let normal = Vector2D::new(diff.y / length, -diff.x / length);

        let midpoint = Vector2D::new(f64::midpoint(p1.x, p2.x), f64::midpoint(p1.y, p2.y));

        Self {
            p1,
            p2,
            midpoint,
            length,
            normal,
        }
    }
}

/// Solves a 2D Laplace problem (e.g., potential flow, steady-state heat conduction)
/// using the Boundary Element Method (BEM) with constant elements.
///
/// This function discretizes the boundary of the domain into elements and applies
/// boundary conditions to solve for unknown potentials or fluxes on the boundary.
///
/// # Arguments
/// * `points` - A `Vec` of `(x, y)` tuples defining the vertices of the boundary polygon.
/// * `bcs` - A `Vec` of `BoundaryCondition` for each element, specifying known potential or flux.
///
/// # Returns
/// A `Result` containing a tuple `(u, q)`, where `u` is a `Vec<f64>` of potentials
/// and `q` is a `Vec<f64>` of normal fluxes on each element. Returns an `Err` string
/// if the system is ill-posed or has no unique solution.
///
/// # Errors
///
/// This function will return an error if the number of points and boundary conditions
/// do not match, or if the BEM system has no unique solution.
pub fn solve_laplace_bem_2d(
    points: &[(f64, f64)],
    bcs: &[BoundaryCondition<f64>],
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = points.len();

    if n != bcs.len() {
        return Err("Number of points and \
             boundary conditions must \
             match."
            .to_string());
    }

    let elements: Vec<_> = (0..n)
        .map(|i| {
            Element2D::new(
                Vector2D::new(points[i].0, points[i].1),
                Vector2D::new(points[(i + 1) % n].0, points[(i + 1) % n].1),
            )
        })
        .collect();

    let mut h_mat = Matrix::zeros(n, n);

    let mut g_mat = Matrix::zeros(n, n);

    // Parallel matrix assembly
    let matrices_data: Vec<Vec<(usize, usize, f64, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(n);

            for j in 0..n {
                if i == j {
                    let g_ii = elements[i].length / (2.0 * std::f64::consts::PI)
                        * (1.0 - (elements[i].length / 2.0).ln());

                    row.push((i, j, 0.0, g_ii)); // Diagonal H will be set later via rigid body motion trick
                } else {
                    let r_vec = elements[j].midpoint - elements[i].midpoint;

                    let r = r_vec.norm();

                    let dot = r_vec
                        .x
                        .mul_add(elements[j].normal.x, r_vec.y * elements[j].normal.y);

                    let h_ij = -dot / (2.0 * std::f64::consts::PI * r * r);

                    let g_ij = -1.0 / (2.0 * std::f64::consts::PI) * r.ln();

                    row.push((i, j, h_ij * elements[j].length, g_ij * elements[j].length));
                }
            }

            row
        })
        .collect();

    for row in matrices_data {
        for (i, j, h, g) in row {
            *h_mat.get_mut(i, j) = h;

            *g_mat.get_mut(i, j) = g;
        }
    }

    // Rigid body motion trick: sum of row H_ij = 0
    // So H_ii = -sum_{j!=i} H_ij
    for i in 0..n {
        let mut row_sum = 0.0;

        for j in 0..n {
            if i != j {
                row_sum += *h_mat.get(i, j);
            }
        }

        *h_mat.get_mut(i, i) = -row_sum;
    }

    let mut a_mat = Matrix::zeros(n, n);

    let mut b_vec = vec![0.0; n];

    for (i, b_val) in b_vec.iter_mut().enumerate() {
        for (j, bc) in bcs.iter().enumerate() {
            match bc {
                // Unknown depends on element j's BC type
                | BoundaryCondition::Potential(u_val) => {
                    // Unknown is flux q_j. Equation side: -G_ij * q_j
                    *a_mat.get_mut(i, j) = -*g_mat.get(i, j);

                    // Known contribution from H_ij * u_j goes to RHS with minus
                    *b_val -= *h_mat.get(i, j) * u_val;
                },
                | BoundaryCondition::Flux(q_val) => {
                    // Unknown is potential u_j. Equation side: H_ij * u_j
                    *a_mat.get_mut(i, j) = *h_mat.get(i, j);

                    // Known contribution from G_ij * q_j goes to RHS
                    *b_val += *g_mat.get(i, j) * q_val;
                },
            }
        }
    }

    let LinearSolution::Unique(solution) = solve_linear_system(&a_mat, &b_vec)? else {
        return Err("BEM system has no unique solution.".to_string());
    };

    let mut u = vec![0.0; n];

    let mut q = vec![0.0; n];

    let mut sol_idx = 0;

    for i in 0..n {
        match bcs[i] {
            | BoundaryCondition::Potential(u_val) => {
                u[i] = u_val;

                q[i] = solution[sol_idx];

                sol_idx += 1;
            },
            | BoundaryCondition::Flux(q_val) => {
                q[i] = q_val;

                u[i] = solution[sol_idx];

                sol_idx += 1;
            },
        }
    }

    Ok((u, q))
}

/// Scenario for 2D BEM: Simulates potential flow around a cylinder.
///
/// This function sets up a circular boundary and applies boundary conditions
/// corresponding to a uniform flow in the x-direction. It then uses the BEM solver
/// to calculate the potential and flux on the cylinder's surface.
///
/// # Returns
/// A `Result` containing a tuple `(u, q)` of potentials and fluxes on the cylinder surface,
/// or an error string if the BEM system cannot be solved.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_laplace_bem_2d`
/// function encounters an error.
pub fn simulate_2d_cylinder_scenario() -> Result<(Vec<f64>, Vec<f64>), String> {
    let n_points = 40;

    let radius = 1.0;

    let mut points = Vec::new();

    let mut bcs = Vec::new();

    for i in 0..n_points {
        let angle = 2.0 * std::f64::consts::PI * (f64::from(i)) / (f64::from(n_points));

        let (x, y) = (radius * angle.cos(), radius * angle.sin());

        points.push((x, y));

        bcs.push(BoundaryCondition::Potential(1.0 * x));
    }

    solve_laplace_bem_2d(&points, &bcs)
}

/// Evaluates the potential at an internal point in the domain after solving the boundary.
///
/// # Arguments
/// * `point` - The `(x, y)` coordinate of the internal point.
/// * `elements` - The boundary elements.
/// * `u` - The solved boundary potentials.
/// * `q` - The solved boundary fluxes.
#[must_use]
pub fn evaluate_potential_2d(
    point: (f64, f64),
    elements: &[Element2D],
    u: &[f64],
    q: &[f64],
) -> f64 {
    let p = Vector2D::new(point.0, point.1);

    let mut result = 0.0;

    for i in 0..elements.len() {
        let r_vec = elements[i].midpoint - p;

        let r = r_vec.norm();

        let dot = r_vec
            .x
            .mul_add(elements[i].normal.x, r_vec.y * elements[i].normal.y);

        let h_ij = -dot / (2.0 * std::f64::consts::PI * r * r);

        let g_ij = -1.0 / (2.0 * std::f64::consts::PI) * r.ln();

        result += (g_ij * elements[i].length).mul_add(q[i], -(h_ij * elements[i].length * u[i]));
    }

    result
}


// ============================================================================
// 3D Laplace BEM on closed triangulated surfaces
// ============================================================================

type V3 = [f64; 3];

fn v_sub(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn v_add(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn v_scale(
    a: V3,
    s: f64,
) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn v_dot(
    a: V3,
    b: V3,
) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn v_cross(
    a: V3,
    b: V3,
) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn v_norm(a: V3) -> f64 {
    v_dot(a, a).sqrt()
}

fn v_mid(
    a: V3,
    b: V3,
) -> V3 {
    v_scale(v_add(a, b), 0.5)
}

/// A closed, consistently oriented triangulated surface.
///
/// Triangles are listed counter-clockwise when seen from outside, so that
/// `(v1 - v0) x (v2 - v0)` points out of the enclosed volume. Use
/// [`SurfaceMesh3D::icosphere`] or [`SurfaceMesh3D::cube`] to generate
/// meshes and [`SurfaceMesh3D::validate`] to check a user-supplied one.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SurfaceMesh3D {
    /// Vertex coordinates.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles as vertex-index triples.
    pub triangles: Vec<[usize; 3]>,
}

impl SurfaceMesh3D {
    /// The three corners of triangle `t`, or `None` if out of range.
    #[must_use]
    pub fn corners(
        &self,
        t: usize,
    ) -> Option<[V3; 3]> {
        let [a, b, c] = *self.triangles.get(t)?;

        Some([*self.vertices.get(a)?, *self.vertices.get(b)?, *self.vertices.get(c)?])
    }

    /// Centroid of triangle `t` (the collocation point), if it exists.
    #[must_use]
    pub fn centroid(
        &self,
        t: usize,
    ) -> Option<[f64; 3]> {
        let [a, b, c] = self.corners(t)?;

        Some(v_scale(v_add(v_add(a, b), c), 1.0 / 3.0))
    }

    /// Outward unit normal of triangle `t`, if it exists and is not degenerate.
    #[must_use]
    pub fn normal(
        &self,
        t: usize,
    ) -> Option<[f64; 3]> {
        let [a, b, c] = self.corners(t)?;

        let n = v_cross(v_sub(b, a), v_sub(c, a));

        let len = v_norm(n);

        (len > 0.0).then(|| v_scale(n, 1.0 / len))
    }

    /// Area of triangle `t`, if it exists.
    #[must_use]
    pub fn area(
        &self,
        t: usize,
    ) -> Option<f64> {
        let [a, b, c] = self.corners(t)?;

        Some(0.5 * v_norm(v_cross(v_sub(b, a), v_sub(c, a))))
    }

    /// Total surface area.
    #[must_use]
    pub fn total_area(&self) -> f64 {
        (0..self.triangles.len()).filter_map(|t| self.area(t)).sum()
    }

    /// Enclosed volume (divergence theorem; positive for outward orientation).
    #[must_use]
    pub fn volume(&self) -> f64 {
        (0..self.triangles.len())
            .filter_map(|t| self.corners(t))
            .map(|[a, b, c]| v_dot(a, v_cross(b, c)) / 6.0)
            .sum()
    }

    /// Checks that the mesh is a closed, consistently oriented surface with
    /// non-degenerate triangles and positive enclosed volume.
    ///
    /// # Errors
    /// Describes the first defect found.
    pub fn validate(&self) -> Result<(), String> {
        if self.triangles.len() < 4 {
            return Err("a closed surface needs at least 4 triangles".to_string());
        }

        let nv = self.vertices.len();

        let mut edges: std::collections::HashMap<(usize, usize), usize> = std::collections::HashMap::new();

        for (t, tri) in self.triangles.iter().enumerate() {
            if tri.iter().any(|&v| v >= nv) {
                return Err(format!("triangle {t} references a missing vertex"));
            }

            if self.area(t).is_none_or(|a| a <= 0.0 || !a.is_finite()) {
                return Err(format!("triangle {t} is degenerate"));
            }

            for k in 0..3 {
                let e = (tri[k], tri[(k + 1) % 3]);

                *edges.entry(e).or_insert(0) += 1;
            }
        }

        for (&(a, b), &count) in &edges {
            if count != 1 || edges.get(&(b, a)) != Some(&1) {
                return Err(format!("edge ({a}, {b}) is not shared by exactly two oppositely oriented triangles"));
            }
        }

        if self.volume() <= 0.0 {
            return Err("surface is inward oriented or encloses no volume".to_string());
        }

        Ok(())
    }

    /// A geodesic sphere: the icosahedron with every triangle split into four
    /// `subdivisions` times (`20 * 4^subdivisions` triangles), vertices
    /// projected onto the sphere.
    #[must_use]
    pub fn icosphere(
        radius: f64,
        center: [f64; 3],
        subdivisions: usize,
    ) -> Self {
        let t = f64::midpoint(1.0, 5.0_f64.sqrt());

        let base: [V3; 12] = [
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ];

        let mut vertices: Vec<V3> = base.iter().map(|&v| v_scale(v, 1.0 / v_norm(v))).collect();

        let mut triangles: Vec<[usize; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];

        for _ in 0..subdivisions {
            let mut cache: std::collections::HashMap<(usize, usize), usize> = std::collections::HashMap::new();

            let mut next = Vec::with_capacity(triangles.len() * 4);

            for &[a, b, c] in &triangles {
                let mut midpoint = |i: usize, j: usize, verts: &mut Vec<V3>| {
                    let key = (i.min(j), i.max(j));

                    *cache.entry(key).or_insert_with(|| {
                        let m = v_mid(verts[i], verts[j]);

                        verts.push(v_scale(m, 1.0 / v_norm(m)));

                        verts.len() - 1
                    })
                };

                let ab = midpoint(a, b, &mut vertices);

                let bc = midpoint(b, c, &mut vertices);

                let ca = midpoint(c, a, &mut vertices);

                next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }

            triangles = next;
        }

        let vertices = vertices.iter().map(|&v| v_add(v_scale(v, radius), center)).collect();

        let mut mesh = Self { vertices, triangles };

        if mesh.volume() < 0.0 {
            for tri in &mut mesh.triangles {
                tri.swap(1, 2);
            }
        }

        mesh
    }

    /// An axis-aligned cube of edge length `side` centred at `center`, each
    /// face split into `n x n` squares of two triangles (`12 n^2`
    /// triangles in total), with outward orientation.
    #[must_use]
    pub fn cube(
        side: f64,
        center: [f64; 3],
        n: usize,
    ) -> Self {
        let n = n.max(1);

        let mut index: std::collections::HashMap<[usize; 3], usize> = std::collections::HashMap::new();

        let mut vertices: Vec<V3> = Vec::new();

        let mut triangles: Vec<[usize; 3]> = Vec::new();

        let mut vertex = |ijk: [usize; 3]| -> usize {
            *index.entry(ijk).or_insert_with(|| {
                let mut p = center;

                for (c, &i) in p.iter_mut().zip(&ijk) {
                    *c += side * (i as f64 / n as f64 - 0.5);
                }

                vertices.push(p);

                vertices.len() - 1
            })
        };

        for axis in 0..3 {
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);

            for (fixed, outward) in [(n, true), (0, false)] {
                for i in 0..n {
                    for j in 0..n {
                        let mut corner = |du: usize, dv: usize| {
                            let mut ijk = [0usize; 3];

                            ijk[axis] = fixed;

                            ijk[u] = i + du;

                            ijk[v] = j + dv;

                            vertex(ijk)
                        };

                        let (p00, p10, p11, p01) = (corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1));

                        if outward {
                            triangles.push([p00, p10, p11]);

                            triangles.push([p00, p11, p01]);
                        } else {
                            triangles.push([p00, p11, p10]);

                            triangles.push([p00, p01, p11]);
                        }
                    }
                }
            }
        }

        Self { vertices, triangles }
    }
}

/// Degree-5, 7-point Gauss rule on the reference triangle: barycentric
/// coordinates and weights (summing to one).
fn triangle_rule() -> [([f64; 3], f64); 7] {
    let s15 = 15.0_f64.sqrt();

    let a1 = (9.0 + 2.0 * s15) / 21.0;

    let b1 = (6.0 - s15) / 21.0;

    let a2 = (9.0 - 2.0 * s15) / 21.0;

    let b2 = (6.0 + s15) / 21.0;

    let w1 = (155.0 - s15) / 1200.0;

    let w2 = (155.0 + s15) / 1200.0;

    [
        ([1.0 / 3.0; 3], 0.225),
        ([a1, b1, b1], w1),
        ([b1, a1, b1], w1),
        ([b1, b1, a1], w1),
        ([a2, b2, b2], w2),
        ([b2, a2, b2], w2),
        ([b2, b2, a2], w2),
    ]
}

/// `integral over triangle(a, b, c) of 1 / |x - y| dS(y)` by adaptive
/// subdivision: triangles close to `x` (relative to their size) are split
/// into four until the 7-point rule is accurate.
fn inv_r_integral(
    x: V3,
    tri: [V3; 3],
    depth: usize,
) -> f64 {
    let [a, b, c] = tri;

    let centroid = v_scale(v_add(v_add(a, b), c), 1.0 / 3.0);

    let diameter = v_norm(v_sub(a, b)).max(v_norm(v_sub(b, c))).max(v_norm(v_sub(c, a)));

    if depth < 7 && v_norm(v_sub(x, centroid)) < 3.0 * diameter {
        let (ab, bc, ca) = (v_mid(a, b), v_mid(b, c), v_mid(c, a));

        return [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
            .iter()
            .map(|&t| inv_r_integral(x, t, depth + 1))
            .sum();
    }

    let area = 0.5 * v_norm(v_cross(v_sub(b, a), v_sub(c, a)));

    let sum: f64 = triangle_rule()
        .iter()
        .map(|&(l, w)| {
            let y = v_add(v_add(v_scale(a, l[0]), v_scale(b, l[1])), v_scale(c, l[2]));

            w / v_norm(v_sub(x, y)).max(f64::MIN_POSITIVE)
        })
        .sum();

    sum * area
}

/// Analytic `integral over triangle of 1 / |x - y| dS(y)` for a point `x` in
/// the plane of the triangle and inside it (the singular self term).
fn inv_r_integral_in_plane(
    x: V3,
    tri: [V3; 3],
) -> f64 {
    let mut total = 0.0;

    for k in 0..3 {
        let (p, q) = (tri[k], tri[(k + 1) % 3]);

        let edge = v_sub(q, p);

        let len = v_norm(edge);

        let t = v_scale(edge, 1.0 / len);

        let (ap, aq) = (v_sub(p, x), v_sub(q, x));

        let (s1, s2) = (v_dot(ap, t), v_dot(aq, t));

        let d = v_norm(v_sub(ap, v_scale(t, s1)));

        let (r1, r2) = (v_norm(ap), v_norm(aq));

        total += d * ((s2 + r2) / (s1 + r1)).ln();
    }

    total
}

/// Solid angle subtended at `x` by a triangle (Van Oosterom-Strackee),
/// positive when `x` lies on the inner side of an outward-oriented triangle.
fn solid_angle(
    x: V3,
    tri: [V3; 3],
) -> f64 {
    let (a, b, c) = (v_sub(tri[0], x), v_sub(tri[1], x), v_sub(tri[2], x));

    let (la, lb, lc) = (v_norm(a), v_norm(b), v_norm(c));

    let num = v_dot(a, v_cross(b, c));

    let den = la * lb * lc + v_dot(a, b) * lc + v_dot(a, c) * lb + v_dot(b, c) * la;

    if num.abs() <= 1e-14 * la * lb * lc {
        return 0.0;
    }

    2.0 * num.atan2(den)
}

/// The solution of a 3D Laplace BEM problem: potential `u` and outward
/// normal derivative `q = du/dn` on every triangle (constant per element).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Bem3dSolution {
    /// Potential on each triangle.
    pub u: Vec<f64>,
    /// Outward normal derivative on each triangle.
    pub q: Vec<f64>,
}

/// Solves the 3D Laplace equation inside a closed triangulated surface by
/// the collocation Boundary Element Method with constant elements.
///
/// The potential `u` and flux `q = du/dn` (outward normal) are piecewise
/// constant on the triangles and collocated at their centroids. Per
/// element, `bcs` prescribes either the potential (Dirichlet) or the flux
/// (Neumann); the other quantity is solved for. The discretised boundary
/// integral equation is
/// `sum_j Hhat_ij u_j = sum_j G_ij q_j`, with the single-layer kernel
/// `G = 1 / (4 pi r)` integrated over each triangle by adaptive
/// subdivision with a 7-point Gauss rule (the singular self term analytically),
/// the double-layer coefficients `H_ij = -Omega_ij / (4 pi)` given exactly by
/// the solid angle subtended by element `j` at the collocation point `i`, and
/// the free-term diagonal fixed by the rigid-body-motion identity
/// `sum_j Hhat_ij = 0`.
///
/// With Neumann data on every element the potential is determined only up
/// to an additive constant; it is then fixed by the condition that the
/// area-weighted mean of `u` vanishes, and the fluxes must satisfy the
/// compatibility condition `sum_j q_j A_j = 0` to within `1e-6` of the flux
/// scale.
///
/// # Errors
/// Fails for an invalid mesh (see [`SurfaceMesh3D::validate`]), a number of
/// boundary conditions different from the number of triangles, non-finite
/// data, an incompatible pure-Neumann problem, or a singular system.
pub fn solve_laplace_bem_3d_mesh(
    mesh: &SurfaceMesh3D,
    bcs: &[BoundaryCondition<f64>],
) -> Result<Bem3dSolution, String> {
    mesh.validate()?;

    let n = mesh.triangles.len();

    if bcs.len() != n {
        return Err(format!("expected {n} boundary conditions, got {}", bcs.len()));
    }

    if bcs
        .iter()
        .any(|bc| matches!(bc, BoundaryCondition::Potential(v) | BoundaryCondition::Flux(v) if !v.is_finite()))
    {
        return Err("boundary data must be finite".to_string());
    }

    let corners: Vec<[V3; 3]> = (0..n).filter_map(|t| mesh.corners(t)).collect();

    let centroids: Vec<V3> = corners
        .iter()
        .map(|&[a, b, c]| v_scale(v_add(v_add(a, b), c), 1.0 / 3.0))
        .collect();

    // Rows of G and of Hhat (free term and rigid-body diagonal included).
    let rows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let x = centroids[i];

            let mut g_row = vec![0.0; n];

            let mut h_row = vec![0.0; n];

            let mut off_diagonal = 0.0;

            for j in 0..n {
                if i == j {
                    g_row[j] = inv_r_integral_in_plane(x, corners[j]) / (4.0 * std::f64::consts::PI);
                } else {
                    g_row[j] = inv_r_integral(x, corners[j], 0) / (4.0 * std::f64::consts::PI);

                    let h = -solid_angle(x, corners[j]) / (4.0 * std::f64::consts::PI);

                    h_row[j] = h;

                    off_diagonal += h;
                }
            }

            h_row[i] = -off_diagonal;

            (g_row, h_row)
        })
        .collect();

    let all_flux = bcs.iter().all(|bc| matches!(bc, BoundaryCondition::Flux(_)));

    let areas: Vec<f64> = (0..n).filter_map(|t| mesh.area(t)).collect();

    if all_flux {
        let scale: f64 = bcs
            .iter()
            .zip(&areas)
            .map(|(bc, a)| if let BoundaryCondition::Flux(q) = bc { (q * a).abs() } else { 0.0 })
            .sum();

        let net: f64 = bcs
            .iter()
            .zip(&areas)
            .map(|(bc, a)| if let BoundaryCondition::Flux(q) = bc { q * a } else { 0.0 })
            .sum();

        if net.abs() > 1e-6 * scale.max(f64::MIN_POSITIVE) {
            return Err("pure Neumann data violates the compatibility condition (net flux must vanish)".to_string());
        }
    }

    // Unknown x_j = q_j on Dirichlet elements, u_j on Neumann elements.
    let dim = if all_flux { n + 1 } else { n };

    let mut a_mat = Mat::zeros(dim, dim);

    let mut rhs = vec![0.0; dim];

    for (i, (g_row, h_row)) in rows.iter().enumerate() {
        for (j, bc) in bcs.iter().enumerate() {
            match *bc {
                BoundaryCondition::Potential(u_val) => {
                    a_mat.set(i, j, -g_row[j]);

                    rhs[i] -= h_row[j] * u_val;
                },
                BoundaryCondition::Flux(q_val) => {
                    a_mat.set(i, j, h_row[j]);

                    rhs[i] += g_row[j] * q_val;
                },
            }
        }

        if all_flux {
            // Lagrange column removes the constant null space.
            a_mat.set(i, n, 1.0);
        }
    }

    if all_flux {
        for (j, &area) in areas.iter().enumerate() {
            a_mat.set(n, j, area);
        }
    }

    let lu = lu_factor(&a_mat).map_err(|e| format!("BEM system is singular: {e:?}"))?;

    let solution = lu.solve(&rhs);

    if solution.iter().any(|v| !v.is_finite()) {
        return Err("BEM system has no finite solution".to_string());
    }

    let mut u = vec![0.0; n];

    let mut q = vec![0.0; n];

    for (j, bc) in bcs.iter().enumerate() {
        match *bc {
            BoundaryCondition::Potential(u_val) => {
                u[j] = u_val;

                q[j] = solution[j];
            },
            BoundaryCondition::Flux(q_val) => {
                q[j] = q_val;

                u[j] = solution[j];
            },
        }
    }

    Ok(Bem3dSolution { u, q })
}

/// Evaluates the potential at an interior point from a solved boundary:
/// `u(x) = sum_j (G_xj q_j + Omega_xj / (4 pi) u_j)`.
///
/// Accuracy degrades as `point` approaches the surface (closer than about
/// one element size).
///
/// # Errors
/// Fails if the solution does not match the mesh.
pub fn evaluate_potential_3d(
    point: [f64; 3],
    mesh: &SurfaceMesh3D,
    solution: &Bem3dSolution,
) -> Result<f64, String> {
    let n = mesh.triangles.len();

    if solution.u.len() != n || solution.q.len() != n {
        return Err("solution does not match the mesh".to_string());
    }

    let mut result = 0.0;

    for t in 0..n {
        let tri = mesh.corners(t).ok_or("invalid triangle index")?;

        let g = inv_r_integral(point, tri, 0) / (4.0 * std::f64::consts::PI);

        let h = solid_angle(point, tri) / (4.0 * std::f64::consts::PI);

        result += g * solution.q[t] + h * solution.u[t];
    }

    Ok(result)
}

/// Self-check of the 3D BEM: solves the Dirichlet problem `u = x` on a
/// 320-triangle unit icosphere with [`solve_laplace_bem_3d_mesh`] and
/// verifies the recovered flux against the exact `q = n_x`.
///
/// This is the former argument-less entry point; the general solver is
/// [`solve_laplace_bem_3d_mesh`].
///
/// # Errors
/// Returns an error if the solve fails or the flux error exceeds 5 %.
pub fn solve_laplace_bem_3d() -> Result<(), String> {
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 2);

    let bcs: Vec<BoundaryCondition<f64>> = (0..mesh.triangles.len())
        .map(|t| BoundaryCondition::Potential(mesh.centroid(t).map_or(0.0, |c| c[0])))
        .collect();

    let sol = solve_laplace_bem_3d_mesh(&mesh, &bcs)?;

    let worst = (0..mesh.triangles.len())
        .map(|t| (sol.q[t] - mesh.normal(t).map_or(0.0, |n| n[0])).abs())
        .fold(0.0, f64::max);

    if worst < 0.05 {
        Ok(())
    } else {
        Err(format!("3D BEM self-check failed: flux error {worst}"))
    }
}
