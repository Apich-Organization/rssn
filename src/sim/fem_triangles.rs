//! # Finite elements on unstructured triangle meshes
//!
//! `-∇·(k ∇u) + c u = f` on a polygonal domain given as a triangle mesh,
//! `u = g` on the boundary:
//!
//! * **P1** (linear, 3 nodes) and **P2** (quadratic, 6 nodes) Lagrange
//!   elements, assembled with a degree-5 seven-point quadrature into a
//!   compressed sparse row matrix;
//! * the symmetric positive definite system solved by **conjugate
//!   gradients with a Jacobi preconditioner** (Dirichlet rows eliminated
//!   symmetrically);
//! * a **residual a-posteriori error estimator** `η_T² = h_T² ‖f - c u_h +
//!   ∇·(k∇u_h)‖²_T + ½ Σ_e h_e ‖[k ∂_n u_h]‖²_e` (P1);
//! * **adaptive refinement**: Dörfler marking and conforming bisection
//!   (every marked edge is split at its midpoint; a triangle is bisected
//!   along its longest marked edge, recursively), so meshes concentrate
//!   where the solution is singular, e.g. at the re-entrant corner of an
//!   L-shaped domain.

use std::collections::HashMap;

/// A triangle mesh: vertex coordinates and counter-clockwise triangles.
#[derive(Clone, Debug, Default)]
pub struct TriMesh {
    /// Vertex coordinates.
    pub nodes: Vec<[f64; 2]>,
    /// Vertex indices of every triangle.
    pub triangles: Vec<[usize; 3]>,
}

impl TriMesh {
    /// The rectangle `[x0, x1] × [y0, y1]` cut into `nx × ny` squares, each
    /// split into two triangles.
    #[must_use]
    pub fn rectangle(
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        nx: usize,
        ny: usize,
    ) -> Self {
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push([x0 + (x1 - x0) * i as f64 / nx as f64, y0 + (y1 - y0) * j as f64 / ny as f64]);
            }
        }
        let id = |i: usize, j: usize| j * (nx + 1) + i;
        let mut triangles = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                triangles.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                triangles.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        Self { nodes, triangles }
    }

    /// The L-shaped domain `[-1, 1]² \ [0, 1] × [-1, 0]` with mesh size
    /// `1/n`.
    #[must_use]
    pub fn l_shape(n: usize) -> Self {
        let square = Self::rectangle(-1.0, 1.0, -1.0, 1.0, 2 * n, 2 * n);
        let keep: Vec<[usize; 3]> = square
            .triangles
            .iter()
            .copied()
            .filter(|t| {
                let cx = t.iter().map(|&v| square.nodes[v][0]).sum::<f64>() / 3.0;
                let cy = t.iter().map(|&v| square.nodes[v][1]).sum::<f64>() / 3.0;
                !(cx > 0.0 && cy < 0.0)
            })
            .collect();
        Self { nodes: square.nodes, triangles: keep }.compacted()
    }

    /// The mesh without unused vertices.
    #[must_use]
    pub fn compacted(self) -> Self {
        let mut map = vec![usize::MAX; self.nodes.len()];
        let mut nodes = Vec::new();
        for t in &self.triangles {
            for &v in t {
                if map[v] == usize::MAX {
                    map[v] = nodes.len();
                    nodes.push(self.nodes[v]);
                }
            }
        }
        let triangles = self.triangles.iter().map(|t| t.map(|v| map[v])).collect();
        Self { nodes, triangles }
    }

    /// Edges `(a, b)`, `a < b`, with the triangles that contain them.
    fn edges(&self) -> HashMap<(usize, usize), Vec<usize>> {
        let mut edges: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (k, t) in self.triangles.iter().enumerate() {
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                edges.entry((a.min(b), a.max(b))).or_default().push(k);
            }
        }
        edges
    }

    /// Vertices on the boundary (on an edge with a single triangle).
    #[must_use]
    pub fn boundary_vertices(&self) -> Vec<bool> {
        let mut on = vec![false; self.nodes.len()];
        for ((a, b), ts) in self.edges() {
            if ts.len() == 1 {
                on[a] = true;
                on[b] = true;
            }
        }
        on
    }

    fn area(
        &self,
        t: &[usize; 3],
    ) -> f64 {
        let [a, b, c] = t.map(|v| self.nodes[v]);
        0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
    }

    /// Refines the triangles flagged in `marked` (and as many neighbours as
    /// conformity needs) by bisection.
    #[must_use]
    pub fn refine(
        &self,
        marked: &[bool],
    ) -> Self {
        let length2 = |a: usize, b: usize| {
            let (p, q) = (self.nodes[a], self.nodes[b]);
            (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)
        };
        let key = |a: usize, b: usize| (a.min(b), a.max(b));
        let longest = |t: &[usize; 3]| -> (usize, usize) {
            (0..3).map(|i| key(t[i], t[(i + 1) % 3])).max_by(|&(a, b), &(c, d)| length2(a, b).total_cmp(&length2(c, d))).unwrap_or((0, 0))
        };
        // Mark the longest edge of every marked triangle, then close: a
        // triangle with any marked edge has its longest edge marked too.
        let mut split: HashMap<(usize, usize), usize> = HashMap::new();
        let mut nodes = self.nodes.clone();
        let mark = |e: (usize, usize), nodes: &mut Vec<[f64; 2]>, split: &mut HashMap<(usize, usize), usize>| -> bool {
            if split.contains_key(&e) {
                return false;
            }
            let (p, q) = (nodes[e.0], nodes[e.1]);
            nodes.push([f64::midpoint(p[0], q[0]), f64::midpoint(p[1], q[1])]);
            split.insert(e, nodes.len() - 1);
            true
        };
        for (t, &m) in self.triangles.iter().zip(marked) {
            if m {
                mark(longest(t), &mut nodes, &mut split);
            }
        }
        loop {
            let mut changed = false;
            for t in &self.triangles {
                let any = (0..3).any(|i| split.contains_key(&key(t[i], t[(i + 1) % 3])));
                if any {
                    changed |= mark(longest(t), &mut nodes, &mut split);
                }
            }
            if !changed {
                break;
            }
        }
        // Bisect recursively along the longest marked edge.
        #[allow(clippy::items_after_statements)] // a recursive helper local to refinement
        fn bisect(
            t: [usize; 3],
            split: &HashMap<(usize, usize), usize>,
            nodes: &[[f64; 2]],
            out: &mut Vec<[usize; 3]>,
        ) {
            let key = |a: usize, b: usize| (a.min(b), a.max(b));
            let length2 = |a: usize, b: usize| {
                let (p, q) = (nodes[a], nodes[b]);
                (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)
            };
            let best = (0..3)
                .filter(|&i| split.contains_key(&key(t[i], t[(i + 1) % 3])))
                .max_by(|&i, &j| length2(t[i], t[(i + 1) % 3]).total_cmp(&length2(t[j], t[(j + 1) % 3])));
            let Some(i) = best else {
                out.push(t);
                return;
            };
            let (a, b, c) = (t[i], t[(i + 1) % 3], t[(i + 2) % 3]);
            let m = split[&key(a, b)];
            bisect([a, m, c], split, nodes, out);
            bisect([m, b, c], split, nodes, out);
        }
        let mut triangles = Vec::with_capacity(self.triangles.len() * 2);
        for &t in &self.triangles {
            bisect(t, &split, &nodes, &mut triangles);
        }
        Self { nodes, triangles }
    }
}

/// Polynomial degree of the elements.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degree {
    /// Linear elements.
    P1,
    /// Quadratic elements.
    P2,
}

/// The problem data.
pub struct Problem<'a> {
    /// Diffusion coefficient `k(x, y) > 0`.
    pub k: &'a dyn Fn(f64, f64) -> f64,
    /// Reaction coefficient `c(x, y) >= 0`.
    pub c: &'a dyn Fn(f64, f64) -> f64,
    /// Source `f(x, y)`.
    pub f: &'a dyn Fn(f64, f64) -> f64,
    /// Dirichlet data `g(x, y)` on the boundary.
    pub g: &'a dyn Fn(f64, f64) -> f64,
}

/// A finite-element solution: degrees of freedom with their coordinates
/// and the element connectivity (3 or 6 indices per triangle).
#[derive(Clone, Debug)]
pub struct Solution {
    /// Coordinates of every degree of freedom.
    pub points: Vec<[f64; 2]>,
    /// Values at the degrees of freedom.
    pub values: Vec<f64>,
    /// Degrees of freedom of each element (vertices, then edge midpoints
    /// opposite... in the order `v0, v1, v2, m01, m12, m20` for P2).
    pub elements: Vec<Vec<usize>>,
    /// The element degree.
    pub degree: Degree,
    /// Conjugate-gradient iterations used.
    pub iterations: usize,
}

/// Strang–Fix seven-point rule on the reference triangle (degree 5):
/// barycentric points and weights summing to 1.
fn quadrature() -> Vec<([f64; 3], f64)> {
    let a1 = 0.059_715_871_789_769_82;
    let b1 = 0.470_142_064_105_115_1;
    let a2 = 0.797_426_985_353_087_3;
    let b2 = 0.101_286_507_323_456_3;
    let w0 = 0.225;
    let w1 = 0.132_394_152_788_506_2;
    let w2 = 0.125_939_180_544_827_1;
    vec![
        ([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], w0),
        ([a1, b1, b1], w1),
        ([b1, a1, b1], w1),
        ([b1, b1, a1], w1),
        ([a2, b2, b2], w2),
        ([b2, a2, b2], w2),
        ([b2, b2, a2], w2),
    ]
}

/// Values and barycentric-gradient coefficients of the basis at `λ`:
/// `(φ_i, ∂φ_i/∂λ_j)`.
fn basis(
    degree: Degree,
    l: [f64; 3],
) -> (Vec<f64>, Vec<[f64; 3]>) {
    match degree {
        | Degree::P1 => (l.to_vec(), vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        | Degree::P2 => {
            let v = |i: usize| l[i] * (2.0 * l[i] - 1.0);
            let dv = |i: usize| {
                let mut g = [0.0; 3];
                g[i] = 4.0 * l[i] - 1.0;
                g
            };
            let e = |i: usize, j: usize| 4.0 * l[i] * l[j];
            let de = |i: usize, j: usize| {
                let mut g = [0.0; 3];
                g[i] = 4.0 * l[j];
                g[j] = 4.0 * l[i];
                g
            };
            (vec![v(0), v(1), v(2), e(0, 1), e(1, 2), e(2, 0)], vec![dv(0), dv(1), dv(2), de(0, 1), de(1, 2), de(2, 0)])
        },
    }
}

/// A sparse matrix in compressed-row form.
struct Csr {
    starts: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
}

impl Csr {
    fn from_triplets(
        n: usize,
        triplets: &HashMap<(usize, usize), f64>,
    ) -> Self {
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for (&(i, j), &v) in triplets {
            rows[i].push((j, v));
        }
        let mut starts = Vec::with_capacity(n + 1);
        let (mut columns, mut values) = (Vec::new(), Vec::new());
        starts.push(0);
        for mut row in rows {
            row.sort_unstable_by_key(|e| e.0);
            for (j, v) in row {
                columns.push(j);
                values.push(v);
            }
            starts.push(columns.len());
        }
        Self { starts, columns, values }
    }

    fn apply(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        (0..self.starts.len() - 1)
            .map(|i| (self.starts[i]..self.starts[i + 1]).map(|k| self.values[k] * x[self.columns[k]]).sum())
            .collect()
    }

    fn diagonal(&self) -> Vec<f64> {
        (0..self.starts.len() - 1)
            .map(|i| (self.starts[i]..self.starts[i + 1]).find(|&k| self.columns[k] == i).map_or(1.0, |k| self.values[k]))
            .collect()
    }
}

/// Jacobi-preconditioned conjugate gradients; returns the solution and the
/// iteration count.
fn pcg(
    a: &Csr,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> (Vec<f64>, usize) {
    let n = b.len();
    let inv_diag: Vec<f64> = a.diagonal().iter().map(|d| if d.abs() > 0.0 { 1.0 / d } else { 1.0 }).collect();
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut z: Vec<f64> = r.iter().zip(&inv_diag).map(|(r, d)| r * d).collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
    for iteration in 0..max_iterations {
        if r.iter().map(|v| v * v).sum::<f64>().sqrt() <= tolerance * b_norm {
            return (x, iteration);
        }
        let ap = a.apply(&p);
        let alpha = rz / p.iter().zip(&ap).map(|(a, b)| a * b).sum::<f64>();
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        z = r.iter().zip(&inv_diag).map(|(r, d)| r * d).collect();
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    (x, max_iterations)
}

/// Solves the problem on `mesh` with elements of `degree`.
///
/// # Errors
/// A degenerate (zero-area) triangle.
pub fn solve(
    mesh: &TriMesh,
    degree: Degree,
    problem: &Problem<'_>,
) -> Result<Solution, String> {
    // Degrees of freedom: vertices, then (P2) one per edge.
    let mut points = mesh.nodes.clone();
    let mut boundary = mesh.boundary_vertices();
    let mut elements: Vec<Vec<usize>> = mesh.triangles.iter().map(|t| t.to_vec()).collect();
    if degree == Degree::P2 {
        let edges = mesh.edges();
        let mut edge_dof: HashMap<(usize, usize), usize> = HashMap::new();
        let mut sorted: Vec<_> = edges.iter().collect();
        sorted.sort_by_key(|(k, _)| **k);
        for (&(a, b), ts) in sorted {
            let (p, q) = (mesh.nodes[a], mesh.nodes[b]);
            edge_dof.insert((a, b), points.len());
            points.push([f64::midpoint(p[0], q[0]), f64::midpoint(p[1], q[1])]);
            boundary.push(ts.len() == 1);
        }
        for (element, t) in elements.iter_mut().zip(&mesh.triangles) {
            for (i, j) in [(0, 1), (1, 2), (2, 0)] {
                element.push(edge_dof[&(t[i].min(t[j]), t[i].max(t[j]))]);
            }
        }
    }
    let n = points.len();
    let rule = quadrature();
    let mut triplets: HashMap<(usize, usize), f64> = HashMap::new();
    let mut load = vec![0.0; n];
    for (t, dofs) in mesh.triangles.iter().zip(&elements) {
        let area = mesh.area(t);
        if area.abs() < 1e-300 {
            return Err("degenerate triangle".to_owned());
        }
        let [a, b, c] = t.map(|v| mesh.nodes[v]);
        // ∇λ_i = (y_j - y_k, x_k - x_j) / (2A)
        let grad_l = [
            [(b[1] - c[1]) / (2.0 * area), (c[0] - b[0]) / (2.0 * area)],
            [(c[1] - a[1]) / (2.0 * area), (a[0] - c[0]) / (2.0 * area)],
            [(a[1] - b[1]) / (2.0 * area), (b[0] - a[0]) / (2.0 * area)],
        ];
        for &(l, w) in &rule {
            let (x, y) = (l[0] * a[0] + l[1] * b[0] + l[2] * c[0], l[0] * a[1] + l[1] * b[1] + l[2] * c[1]);
            let weight = w * area.abs();
            let (phi, dphi) = basis(degree, l);
            let grads: Vec<[f64; 2]> =
                dphi.iter().map(|d| [0, 1].map(|k| d[0] * grad_l[0][k] + d[1] * grad_l[1][k] + d[2] * grad_l[2][k])).collect();
            let (kv, cv, fv) = ((problem.k)(x, y), (problem.c)(x, y), (problem.f)(x, y));
            for (i, &gi) in dofs.iter().enumerate() {
                load[gi] += weight * fv * phi[i];
                for (j, &gj) in dofs.iter().enumerate() {
                    let value = weight * (kv * (grads[i][0] * grads[j][0] + grads[i][1] * grads[j][1]) + cv * phi[i] * phi[j]);
                    *triplets.entry((gi, gj)).or_insert(0.0) += value;
                }
            }
        }
    }
    // Dirichlet rows: u = g, eliminated symmetrically from other rows.
    let fixed: Vec<Option<f64>> =
        points.iter().zip(&boundary).map(|(p, &on)| on.then(|| (problem.g)(p[0], p[1]))).collect();
    let mut reduced: HashMap<(usize, usize), f64> = HashMap::with_capacity(triplets.len());
    for (&(i, j), &v) in &triplets {
        match (fixed[i], fixed[j]) {
            | (Some(_), _) => {},
            | (None, Some(gj)) => load[i] -= v * gj,
            | (None, None) => {
                reduced.insert((i, j), v);
            },
        }
    }
    for (i, g) in fixed.iter().enumerate() {
        if let Some(g) = g {
            reduced.insert((i, i), 1.0);
            load[i] = *g;
        }
    }
    let matrix = Csr::from_triplets(n, &reduced);
    let (values, iterations) = pcg(&matrix, &load, 1e-12, 20 * n + 100);
    Ok(Solution { points, values, elements, degree, iterations })
}

impl Solution {
    /// The maximum nodal error against `exact`.
    #[must_use]
    pub fn max_nodal_error(
        &self,
        exact: &dyn Fn(f64, f64) -> f64,
    ) -> f64 {
        self.points.iter().zip(&self.values).map(|(p, v)| (v - exact(p[0], p[1])).abs()).fold(0.0, f64::max)
    }

    /// The `L²` error against `exact`, by quadrature on every element.
    #[must_use]
    pub fn l2_error(
        &self,
        mesh: &TriMesh,
        exact: &dyn Fn(f64, f64) -> f64,
    ) -> f64 {
        let rule = quadrature();
        let mut total = 0.0;
        for (t, dofs) in mesh.triangles.iter().zip(&self.elements) {
            let [a, b, c] = t.map(|v| mesh.nodes[v]);
            let area = mesh.area(t).abs();
            for &(l, w) in &rule {
                let (x, y) = (l[0] * a[0] + l[1] * b[0] + l[2] * c[0], l[0] * a[1] + l[1] * b[1] + l[2] * c[1]);
                let (phi, _) = basis(self.degree, l);
                let uh: f64 = phi.iter().zip(dofs).map(|(p, &d)| p * self.values[d]).sum();
                total += w * area * (uh - exact(x, y)).powi(2);
            }
        }
        total.sqrt()
    }
}

/// Residual error indicators `η_T` of a P1 solution (one per triangle).
#[must_use]
pub fn error_indicators(
    mesh: &TriMesh,
    solution: &Solution,
    problem: &Problem<'_>,
) -> Vec<f64> {
    let rule = quadrature();
    let gradient = |t: &[usize; 3]| -> [f64; 2] {
        let [a, b, c] = t.map(|v| mesh.nodes[v]);
        let area = mesh.area(t);
        let g = [
            [(b[1] - c[1]) / (2.0 * area), (c[0] - b[0]) / (2.0 * area)],
            [(c[1] - a[1]) / (2.0 * area), (a[0] - c[0]) / (2.0 * area)],
            [(a[1] - b[1]) / (2.0 * area), (b[0] - a[0]) / (2.0 * area)],
        ];
        [0, 1].map(|k| (0..3).map(|i| g[i][k] * solution.values[t[i]]).sum())
    };
    let mut eta2: Vec<f64> = mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| mesh.nodes[v]);
            let area = mesh.area(t).abs();
            let h2 = [(a, b), (b, c), (c, a)].iter().map(|(p, q)| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).fold(0.0, f64::max);
            let mut interior = 0.0;
            for &(l, w) in &rule {
                let (x, y) = (l[0] * a[0] + l[1] * b[0] + l[2] * c[0], l[0] * a[1] + l[1] * b[1] + l[2] * c[1]);
                let uh: f64 = (0..3).map(|i| l[i] * solution.values[t[i]]).sum();
                let r = (problem.f)(x, y) - (problem.c)(x, y) * uh;
                interior += w * area * r * r;
            }
            h2 * interior
        })
        .collect();
    for ((a, b), ts) in mesh.edges() {
        if let [t1, t2] = ts.as_slice() {
            let (p, q) = (mesh.nodes[a], mesh.nodes[b]);
            let he = (p[0] - q[0]).hypot(p[1] - q[1]);
            let normal = [(q[1] - p[1]) / he, (p[0] - q[0]) / he];
            let (g1, g2) = (gradient(&mesh.triangles[*t1]), gradient(&mesh.triangles[*t2]));
            let mid = [f64::midpoint(p[0], q[0]), f64::midpoint(p[1], q[1])];
            let k = (problem.k)(mid[0], mid[1]);
            let jump = k * ((g1[0] - g2[0]) * normal[0] + (g1[1] - g2[1]) * normal[1]);
            let share = 0.5 * he * he * jump * jump;
            eta2[*t1] += 0.5 * share;
            eta2[*t2] += 0.5 * share;
        }
    }
    eta2.into_iter().map(f64::sqrt).collect()
}

/// The final mesh, the solution on it and the estimator history.
pub type Adaptive = (TriMesh, Solution, Vec<(usize, f64)>);

/// Adaptive P1 solution.
///
/// Solve, estimate, mark the smallest set of triangles carrying the
/// fraction `theta` of `Σ η²` (Dörfler), refine; `steps` times. Returns the
/// final mesh, solution and the estimator history `(degrees of freedom,
/// (Σ η²)^(1/2))`.
///
/// # Errors
/// A degenerate triangle.
pub fn solve_adaptive(
    mesh: &TriMesh,
    problem: &Problem<'_>,
    steps: usize,
    theta: f64,
) -> Result<Adaptive, String> {
    let mut mesh = mesh.clone();
    let mut history = Vec::with_capacity(steps + 1);
    loop {
        let solution = solve(&mesh, Degree::P1, problem)?;
        let eta = error_indicators(&mesh, &solution, problem);
        let total: f64 = eta.iter().map(|e| e * e).sum();
        history.push((solution.values.len(), total.sqrt()));
        if history.len() > steps {
            return Ok((mesh, solution, history));
        }
        let mut order: Vec<usize> = (0..eta.len()).collect();
        order.sort_by(|&i, &j| eta[j].total_cmp(&eta[i]));
        let mut marked = vec![false; eta.len()];
        let mut acc = 0.0;
        for i in order {
            if acc >= theta * total {
                break;
            }
            marked[i] = true;
            acc += eta[i] * eta[i];
        }
        mesh = mesh.refine(&marked);
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn sine_problem() -> (impl Fn(f64, f64) -> f64, impl Fn(f64, f64) -> f64) {
        (|x: f64, y: f64| (PI * x).sin() * (PI * y).sin(), |x: f64, y: f64| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin())
    }

    #[test]
    fn convergence_rates_of_p1_and_p2() {
        let (exact, source) = sine_problem();
        let one = |_: f64, _: f64| 1.0;
        let zero = |_: f64, _: f64| 0.0;
        let problem = Problem { k: &one, c: &zero, f: &source, g: &zero };
        let error = |n: usize, degree: Degree| {
            let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, n, n);
            let s = solve(&mesh, degree, &problem).unwrap();
            s.l2_error(&mesh, &exact)
        };
        let (p1a, p1b) = (error(8, Degree::P1), error(16, Degree::P1));
        let (p2a, p2b) = (error(8, Degree::P2), error(16, Degree::P2));
        let rate = |a: f64, b: f64| (a / b).log2();
        assert!((rate(p1a, p1b) - 2.0).abs() < 0.2, "P1 rate {}", rate(p1a, p1b));
        assert!((rate(p2a, p2b) - 3.0).abs() < 0.3, "P2 rate {}", rate(p2a, p2b));
        assert!(p2b < p1b / 20.0, "{p1b} {p2b}");
    }

    #[test]
    fn variable_coefficients_and_boundary_data() {
        // u = x² + x y + y² with k = 1 + x, c = 1: f = -∇·(k∇u) + u.
        let exact = |x: f64, y: f64| x.powi(2) + x * y + y.powi(2);
        let k = |x: f64, _: f64| 1.0 + x;
        let c = |_: f64, _: f64| 1.0;
        // ∇u = (2x + y, x + 2y); ∇·(k∇u) = (2x + y) + (1 + x)(2 + 2)
        let f = move |x: f64, y: f64| -((2.0 * x + y) + 4.0 * (1.0 + x)) + exact(x, y);
        let problem = Problem { k: &k, c: &c, f: &f, g: &exact };
        let mesh = TriMesh::rectangle(-1.0, 2.0, 0.0, 1.0, 6, 3);
        // P2 reproduces a quadratic exactly.
        let s = solve(&mesh, Degree::P2, &problem).unwrap();
        assert!(s.max_nodal_error(&exact) < 1e-9, "{}", s.max_nodal_error(&exact));
    }

    #[test]
    fn adaptive_refinement_on_the_l_shape() {
        // u = r^(2/3) sin(2θ/3) is harmonic with a corner singularity.
        let exact = |x: f64, y: f64| {
            let r = x.hypot(y);
            let mut theta = y.atan2(x);
            if theta < 0.0 {
                theta += 2.0 * PI;
            }
            r.powf(2.0 / 3.0) * (2.0 * theta / 3.0).sin()
        };
        let one = |_: f64, _: f64| 1.0;
        let zero = |_: f64, _: f64| 0.0;
        let problem = Problem { k: &one, c: &zero, f: &zero, g: &exact };
        let start = TriMesh::l_shape(2);
        let (mesh, solution, history) = solve_adaptive(&start, &problem, 8, 0.5).unwrap();
        // The estimator decreases and the mesh is graded towards the corner.
        let (first, last) = (history[0].1, history.last().unwrap().1);
        assert!(last < 0.35 * first, "{history:?}");
        let near: usize = mesh.triangles.iter().filter(|t| t.iter().all(|&v| mesh.nodes[v][0].hypot(mesh.nodes[v][1]) < 0.1)).count();
        assert!(near > 10, "{near} small triangles at the corner");
        assert!(solution.max_nodal_error(&exact) < 0.03, "{}", solution.max_nodal_error(&exact));
        // The refined mesh is conforming: every interior edge has two triangles.
        assert!(mesh.edges().values().all(|ts| ts.len() <= 2));
        let area: f64 = mesh.triangles.iter().map(|t| mesh.area(t)).sum();
        assert!((area - 3.0).abs() < 1e-12 && mesh.triangles.iter().all(|t| mesh.area(t) > 0.0));
    }
}
