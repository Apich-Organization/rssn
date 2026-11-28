//! P1/P2 finite elements on interval, triangle and tetrahedron meshes.
//!
//! # Finite elements on simplex meshes in 1D, 2D and 3D (P1 and P2)
//!
//! Lagrange elements of degree 1 and 2 on intervals, triangles and
//! tetrahedra for `-∇·(k ∇u) + c u = f(t, x)`, with boundary conditions
//! assigned per boundary **marker** ([`MarkerBcs`] of
//! [`BoundaryCondition`]):
//!
//! | Condition | Weak-form treatment |
//! |---|---|
//! | Dirichlet `u = g(t, x)` | rows replaced by the nodal value (P2: vertex and edge nodes) |
//! | Neumann `∂u/∂n = g` (outward) | boundary load `∫ k g v` |
//! | Robin `α u + β ∂u/∂n = g` | boundary mass `∫ (kα/β) u v` and load `∫ (k/β) g v` |
//! | Symmetry | homogeneous Neumann |
//!
//! | Capability | 1D | 2D | 3D |
//! |---|---|---|---|
//! | P1 / P2 Poisson, Helmholtz, reaction–diffusion | yes | yes | yes |
//! | θ-method heat equation `u_t = ...` ([`solve_heat`]) | yes | yes | yes |
//! | built-in meshes ([`SimplexMesh::interval`], [`SimplexMesh::rectangle`], [`SimplexMesh::boxed`]) | yes | yes | yes |
//! | conforming local refinement ([`SimplexMesh::refine_marked`]) | yes | yes | yes |
//! | residual estimator and adaptive solve ([`error_indicators`], [`solve_adaptive`]) | yes | yes | yes |
//!
//! Quadrature is the Grundmann–Möller rule of arbitrary degree
//! ([`simplex_rule`]); linear systems use [`LinearSystem`] (sparse LU or
//! ILU-`BiCGStab`). See [`crate::sim::fem_triangles`] for the adaptive 2D
//! solver on triangle meshes.
//!
//! Coefficients: [`FemProblem`] has constant `k`, `c`; [`FemProblemVar`] with
//! [`solve_static_var`] / [`solve_heat_var`] takes closures `k(x)`, `c(x)` evaluated at the
//! quadrature points of the cells (and of Neumann / Robin facets), bit-identical to the constant
//! case when the closures are constant. The residual estimator
//! ([`error_indicators_var`], [`solve_adaptive_var`]) includes `∇k · ∇u_h` in the element
//! residual and `k(x)` in the facet jump. [`SimplexResidualEstimator`] reuses
//! the geometry and quadrature when estimating several solutions on one space.
//!
//! Limitations: no curved boundaries or higher-order
//! geometry; only refinement, no coarsening; no goal-oriented adaptivity and
//! no multigrid. The residual estimator's face terms need an interior face,
//! so in 1D it reduces to the element residual. In 3D,
//! [`SimplexMesh::refine_marked`] is conforming local red–green refinement (a
//! marked tetrahedron is split into eight by edge midpoints, unmarked
//! neighbours only as far as conformity requires) and [`solve_adaptive`]
//! drives it with the residual indicators.
//!
//! Parallelism: the element matrices/loads of the assembly and the cell residual indicators are
//! computed concurrently and scattered in cell order (the triplet list is the serial one). The
//! integrand `f` of [`FemProblem`] must be `Sync`. The edge-jump sums of [`error_indicators`] are
//! taken in sorted facet order (they used to follow the randomised hash order).

use crate::kernels::krylov::Csr;
use crate::sim::boundary::BoundaryCondition;
use crate::sim::boundary::MarkerBcs;
use crate::sim::structured::LinearSolver;
use crate::sim::structured::LinearSystem;
use std::collections::BTreeMap;
use std::collections::HashMap;

fn compositions(
    total: usize,
    parts: usize,
    cur: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if parts == 1 {
        let mut c = cur.clone();
        c.push(total);
        out.push(c);
        return;
    }
    for k in 0..=total {
        cur.push(k);
        compositions(total - k, parts - 1, cur, out);
        cur.pop();
    }
}

/// Grundmann-Moeller quadrature on the `n`-simplex.
///
/// Grundmann–Möller quadrature on the `n`-simplex, exact for polynomials of
/// degree `2s + 1`. Returns barycentric points (`n + 1` coordinates) and
/// weights that sum to one (multiply by the simplex measure).
#[must_use]
pub fn simplex_rule(
    n: usize,
    s: usize,
) -> (Vec<Vec<f64>>, Vec<f64>) {
    if n == 0 {
        return (vec![vec![1.0]], vec![1.0]);
    }
    let d = 2 * s + 1;
    let fact = |m: usize| (1..=m).map(|x| x as f64).product::<f64>();
    let mut pts = Vec::new();
    let mut w = Vec::new();
    for i in 0..=s {
        let denom = (d + n - 2 * i) as f64;
        let coef = (if i % 2 == 0 { 1.0 } else { -1.0 }) * denom.powi(d as i32)
            / (fact(i) * fact(d + n - i));
        let mut comps = Vec::new();
        compositions(s - i, n + 1, &mut Vec::new(), &mut comps);
        for k in comps {
            pts.push(k.iter().map(|&kj| (2 * kj + 1) as f64 / denom).collect());
            w.push(coef);
        }
    }
    let total: f64 = w.iter().sum();
    for x in &mut w {
        *x /= total;
    }
    (pts, w)
}

/// Simplex mesh with marked boundary facets.
#[derive(Clone, Debug)]
pub struct SimplexMesh {
    /// Topological dimension (1, 2 or 3).
    pub dim: usize,
    /// Vertex coordinates (unused entries zero).
    pub nodes: Vec<[f64; 3]>,
    /// Cells as `dim + 1` vertex indices.
    pub cells: Vec<Vec<usize>>,
    /// Boundary facets: `(cell, local face, marker)`; the local face is the one opposite vertex `local face`.
    pub facets: Vec<(usize, usize, i32)>,
    /// Refinement level of every cell: `0` for the cells of a built-in mesh,
    /// `parent + 1` for every cell created by [`SimplexMesh::refine_marked`].
    pub level: Vec<usize>,
}

/// Marker ids of the box faces produced by the built-in meshes.
pub mod marker {
    /// `x = x_min`.
    pub const XMIN: i32 = 0;
    /// `x = x_max`.
    pub const XMAX: i32 = 1;
    /// `y = y_min`.
    pub const YMIN: i32 = 2;
    /// `y = y_max`.
    pub const YMAX: i32 = 3;
    /// `z = z_min`.
    pub const ZMIN: i32 = 4;
    /// `z = z_max`.
    pub const ZMAX: i32 = 5;
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn inv3(m: [[f64; 3]; 3]) -> ([[f64; 3]; 3], f64) {
    let d = det3(m);
    let c = |i: usize, j: usize| {
        let (r0, r1) = ((i + 1) % 3, (i + 2) % 3);
        let (c0, c1) = ((j + 1) % 3, (j + 2) % 3);
        m[r0][c0] * m[r1][c1] - m[r0][c1] * m[r1][c0]
    };
    let mut inv = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            inv[i][j] = c(j, i) / d;
        }
    }
    (inv, d)
}

impl SimplexMesh {
    fn finish(
        dim: usize,
        nodes: Vec<[f64; 3]>,
        mut cells: Vec<Vec<usize>>,
        lo: [f64; 3],
        hi: [f64; 3],
    ) -> Self {
        // positive orientation
        for c in &mut cells {
            let (_, d) = Self::jacobian_of(dim, &nodes, c);
            if d < 0.0 {
                c.swap(0, 1);
            }
        }
        let mut count: HashMap<Vec<usize>, Vec<(usize, usize)>> = HashMap::new();
        for (ci, c) in cells.iter().enumerate() {
            for f in 0..=dim {
                let mut key: Vec<usize> = c
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != f)
                    .map(|(_, &v)| v)
                    .collect();
                key.sort_unstable();
                count.entry(key).or_default().push((ci, f));
            }
        }
        let mut facets = Vec::new();
        for (key, list) in &count {
            if list.len() == 1 {
                let (ci, f) = list[0];
                let mut cen = [0.0; 3];
                for &v in key {
                    for a in 0..3 {
                        cen[a] += nodes[v][a] / key.len() as f64;
                    }
                }
                let tol = 1e-9;
                let mut m = -1;
                for a in 0..dim {
                    if (cen[a] - lo[a]).abs() < tol {
                        m = 2 * a as i32;
                    } else if (cen[a] - hi[a]).abs() < tol {
                        m = 2 * a as i32 + 1;
                    }
                }
                facets.push((ci, f, m));
            }
        }
        facets.sort_unstable();
        let level = vec![0; cells.len()];
        Self {
            dim,
            nodes,
            cells,
            facets,
            level,
        }
    }

    /// Uniform interval mesh with `n` cells on `[lo, hi]` (markers 0 and 1).
    #[must_use]
    pub fn interval(
        n: usize,
        lo: f64,
        hi: f64,
    ) -> Self {
        let nodes: Vec<[f64; 3]> = (0..=n)
            .map(|i| [lo + (hi - lo) * i as f64 / n as f64, 0.0, 0.0])
            .collect();
        let cells = (0..n).map(|i| vec![i, i + 1]).collect();
        Self::finish(1, nodes, cells, [lo, 0.0, 0.0], [hi, 0.0, 0.0])
    }

    /// Triangulated rectangle (two triangles per quad), markers 0-3.
    #[must_use]
    pub fn rectangle(
        nx: usize,
        ny: usize,
        lo: [f64; 2],
        hi: [f64; 2],
    ) -> Self {
        let mut nodes = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push([
                    lo[0] + (hi[0] - lo[0]) * i as f64 / nx as f64,
                    lo[1] + (hi[1] - lo[1]) * j as f64 / ny as f64,
                    0.0,
                ]);
            }
        }
        let id = |i: usize, j: usize| i + (nx + 1) * j;
        let mut cells = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                if (i + j) % 2 == 0 {
                    cells.push(vec![id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                    cells.push(vec![id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
                } else {
                    cells.push(vec![id(i, j), id(i + 1, j), id(i, j + 1)]);
                    cells.push(vec![id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
                }
            }
        }
        Self::finish(2, nodes, cells, [lo[0], lo[1], 0.0], [hi[0], hi[1], 0.0])
    }

    /// Tetrahedral mesh of a box (six tetrahedra per cube, Kuhn subdivision), markers 0-5.
    #[must_use]
    pub fn boxed(
        n: [usize; 3],
        lo: [f64; 3],
        hi: [f64; 3],
    ) -> Self {
        let mut nodes = Vec::new();
        for k in 0..=n[2] {
            for j in 0..=n[1] {
                for i in 0..=n[0] {
                    nodes.push([
                        lo[0] + (hi[0] - lo[0]) * i as f64 / n[0] as f64,
                        lo[1] + (hi[1] - lo[1]) * j as f64 / n[1] as f64,
                        lo[2] + (hi[2] - lo[2]) * k as f64 / n[2] as f64,
                    ]);
                }
            }
        }
        let id = |i: usize, j: usize, k: usize| i + (n[0] + 1) * (j + (n[1] + 1) * k);
        let perms = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let mut cells = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    for p in perms {
                        let mut cur = [i, j, k];
                        let mut v = vec![id(cur[0], cur[1], cur[2])];
                        for &ax in &p {
                            cur[ax] += 1;
                            v.push(id(cur[0], cur[1], cur[2]));
                        }
                        cells.push(v);
                    }
                }
            }
        }
        Self::finish(3, nodes, cells, lo, hi)
    }

    /// Refines the cells flagged in `marked` by conforming edge bisection.
    ///
    /// Every marked tetrahedron is split at the midpoints of its six edges
    /// into eight tetrahedra. The split is carried out as recursive
    /// bisection: a cell is cut along its longest split edge (ties broken by
    /// vertex indices) and the two halves are cut again until no split edge
    /// is left.
    /// Only midpoints of edges of marked cells are inserted, so an unmarked
    /// neighbour is cut along exactly the edges it shares with marked cells
    /// (**red–green closure**) and a cell without a split edge is kept. Both
    /// cells sharing a face cut that face from the same edge midpoints to the
    /// opposite vertex of the face, hence there are no hanging nodes.
    /// Boundary facets inherit the marker of the parent facet containing
    /// them. The refinement level of every produced cell is one more than the
    /// level of its parent; cells that are not touched keep their level.
    ///
    /// In 1D and 2D the same scheme splits a marked cell into 2 or 4 cells.
    /// Entries of `marked` past the last cell are ignored.
    #[must_use]
    pub fn refine_marked(
        &self,
        marked: &[bool],
    ) -> Self {
        let dim = self.dim;
        let edge = |a: usize, b: usize| (a.min(b), a.max(b));
        let mut nodes = self.nodes.clone();
        let mut split: HashMap<(usize, usize), usize> = HashMap::new();
        let mut mid_edge: HashMap<usize, (usize, usize)> = HashMap::new();
        for (ci, c) in self.cells.iter().enumerate() {
            if !marked.get(ci).copied().unwrap_or(false) {
                continue;
            }
            for i in 0..=dim {
                for j in (i + 1)..=dim {
                    let e = edge(c[i], c[j]);
                    split.entry(e).or_insert_with(|| {
                        let (p, q) = (self.nodes[e.0], self.nodes[e.1]);
                        nodes.push([
                            f64::midpoint(p[0], q[0]),
                            f64::midpoint(p[1], q[1]),
                            f64::midpoint(p[2], q[2]),
                        ]);
                        let m = nodes.len() - 1;
                        mid_edge.insert(m, e);
                        m
                    });
                }
            }
        }
        if split.is_empty() {
            return self.clone();
        }
        // Bisection order: longest split edge first (ties by vertex pair), so
        // the same order is used from both sides of every face.
        #[allow(clippy::items_after_statements)] // recursive helper local to the refinement
        fn bisect(
            cell: &[usize],
            split: &HashMap<(usize, usize), usize>,
            nodes: &[[f64; 3]],
            out: &mut Vec<Vec<usize>>,
        ) {
            let edge = |a: usize, b: usize| (a.min(b), a.max(b));
            let rank = |a: usize, b: usize| {
                let (p, q) = (nodes[a], nodes[b]);
                let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                (d2.to_bits(), a.min(b), a.max(b))
            };
            let mut best: Option<(usize, usize)> = None;
            for i in 0..cell.len() {
                for j in (i + 1)..cell.len() {
                    if split.contains_key(&edge(cell[i], cell[j]))
                        && best.is_none_or(|(bi, bj)| {
                            rank(cell[i], cell[j]) > rank(cell[bi], cell[bj])
                        })
                    {
                        best = Some((i, j));
                    }
                }
            }
            let Some((i, j)) = best else {
                out.push(cell.to_vec());
                return;
            };
            let m = split[&edge(cell[i], cell[j])];
            let mut c1 = cell.to_vec();
            let mut c2 = cell.to_vec();
            c1[j] = m;
            c2[i] = m;
            bisect(&c1, split, nodes, out);
            bisect(&c2, split, nodes, out);
        }
        let mut cells: Vec<Vec<usize>> = Vec::new();
        let mut parent_of: Vec<usize> = Vec::new();
        let mut level: Vec<usize> = Vec::new();
        for (ci, c) in self.cells.iter().enumerate() {
            let before = cells.len();
            bisect(c, &split, &nodes, &mut cells);
            let base = self.level.get(ci).copied().unwrap_or(0);
            let red = cells.len() - before > 1;
            for _ in before..cells.len() {
                parent_of.push(ci);
                level.push(if red { base + 1 } else { base });
            }
        }
        // positive orientation (face indices below refer to the final order)
        for c in &mut cells {
            let (_, d) = Self::jacobian_of(dim, &nodes, c);
            if d < 0.0 {
                c.swap(0, 1);
            }
        }
        let mut faces: HashMap<Vec<usize>, Vec<(usize, usize)>> = HashMap::new();
        for (ci, c) in cells.iter().enumerate() {
            for f in 0..=dim {
                let mut key: Vec<usize> = c
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != f)
                    .map(|(_, &v)| v)
                    .collect();
                key.sort_unstable();
                faces.entry(key).or_default().push((ci, f));
            }
        }
        let mut parent_facets: Vec<Vec<(usize, i32)>> = vec![Vec::new(); self.cells.len()];
        for &(ci, f, mk) in &self.facets {
            parent_facets[ci].push((f, mk));
        }
        let mut facets: Vec<(usize, usize, i32)> = Vec::new();
        for (key, list) in &faces {
            let [(ci, f)] = list.as_slice() else { continue };
            let parent = parent_of[*ci];
            let pcell = &self.cells[parent];
            let mut marker = -1;
            for &(pf, mk) in &parent_facets[parent] {
                let face: Vec<usize> = pcell
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != pf)
                    .map(|(_, &v)| v)
                    .collect();
                let inside = key.iter().all(|v| {
                    face.contains(v)
                        || mid_edge
                            .get(v)
                            .is_some_and(|&(a, b)| face.contains(&a) && face.contains(&b))
                });
                if inside {
                    marker = mk;
                    break;
                }
            }
            facets.push((*ci, *f, marker));
        }
        facets.sort_unstable();
        Self {
            dim,
            nodes,
            cells,
            facets,
            level,
        }
    }

    fn jacobian_of(
        dim: usize,
        nodes: &[[f64; 3]],
        c: &[usize],
    ) -> ([[f64; 3]; 3], f64) {
        let mut j = [[0.0; 3]; 3];
        for a in 0..3 {
            for k in 0..dim {
                j[a][k] = nodes[c[k + 1]][a] - nodes[c[0]][a];
            }
        }
        for k in dim..3 {
            j[k][k] = 1.0;
        }
        (j, det3(j))
    }

    /// Cell measure (length, area or volume).
    #[must_use]
    pub fn measure(
        &self,
        cell: usize,
    ) -> f64 {
        let (_, d) = Self::jacobian_of(self.dim, &self.nodes, &self.cells[cell]);
        let fact: f64 = (1..=self.dim).map(|x| x as f64).product();
        d.abs() / fact
    }

    /// Constant gradients of the barycentric coordinates of a cell.
    #[must_use]
    pub fn bary_gradients(
        &self,
        cell: usize,
    ) -> Vec<[f64; 3]> {
        let (j, _) = Self::jacobian_of(self.dim, &self.nodes, &self.cells[cell]);
        let (inv, _) = inv3(j);
        let mut g = vec![[0.0; 3]; self.dim + 1];
        for k in 0..self.dim {
            for a in 0..3 {
                g[k + 1][a] = inv[k][a];
                g[0][a] -= inv[k][a];
            }
        }
        g
    }

    pub(crate) fn facet_measure(
        &self,
        cell: usize,
        face: usize,
    ) -> f64 {
        let c = &self.cells[cell];
        let v: Vec<[f64; 3]> = c
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != face)
            .map(|(_, &q)| self.nodes[q])
            .collect();
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        match self.dim {
            | 1 => 1.0,
            | 2 => {
                let d = sub(v[1], v[0]);
                d[0].hypot(d[1])
            },
            | _ => {
                let (a, b) = (sub(v[1], v[0]), sub(v[2], v[0]));
                let cr = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
            },
        }
    }
}

/// Finite element space of degree 1 or 2 on a simplex mesh.
#[derive(Clone, Debug)]
pub struct FemSpace {
    /// The mesh.
    pub mesh: SimplexMesh,
    /// Polynomial degree.
    pub degree: usize,
    /// Coordinates of the degrees of freedom.
    pub dofs: Vec<[f64; 3]>,
    /// Local-to-global dof map per cell (vertices, then edges `(i, j)`, `i < j`).
    pub cell_dofs: Vec<Vec<usize>>,
}

impl FemSpace {
    /// Builds the space.
    ///
    /// # Errors
    /// For a degree other than 1 or 2.
    pub fn new(
        mesh: SimplexMesh,
        degree: usize,
    ) -> Result<Self, String> {
        if !(1..=2).contains(&degree) {
            return Err("only P1 and P2 elements are available".into());
        }
        let dim = mesh.dim;
        let mut dofs = mesh.nodes.clone();
        let mut cell_dofs: Vec<Vec<usize>> = mesh.cells.clone();
        if degree == 2 {
            let mut emap: HashMap<(usize, usize), usize> = HashMap::new();
            for (ci, c) in mesh.cells.iter().enumerate() {
                for i in 0..=dim {
                    for j in (i + 1)..=dim {
                        let key = (c[i].min(c[j]), c[i].max(c[j]));
                        let id = *emap.entry(key).or_insert_with(|| {
                            let a = mesh.nodes[key.0];
                            let b = mesh.nodes[key.1];
                            dofs.push([
                                f64::midpoint(a[0], b[0]),
                                f64::midpoint(a[1], b[1]),
                                f64::midpoint(a[2], b[2]),
                            ]);
                            dofs.len() - 1
                        });
                        cell_dofs[ci].push(id);
                    }
                }
            }
        }
        Ok(Self {
            mesh,
            degree,
            dofs,
            cell_dofs,
        })
    }

    /// Number of degrees of freedom.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.dofs.len()
    }

    /// Basis values and barycentric derivatives at barycentric point `lam`.
    pub(crate) fn basis(
        &self,
        lam: &[f64],
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let m = lam.len();
        let mut v = Vec::new();
        let mut dl: Vec<Vec<f64>> = Vec::new();
        let unit = |i: usize| {
            let mut e = vec![0.0; m];
            e[i] = 1.0;
            e
        };
        if self.degree == 1 {
            for i in 0..m {
                v.push(lam[i]);
                dl.push(unit(i));
            }
        } else {
            for i in 0..m {
                v.push(lam[i] * (2.0 * lam[i] - 1.0));
                let mut e = vec![0.0; m];
                e[i] = 4.0 * lam[i] - 1.0;
                dl.push(e);
            }
            for i in 0..m {
                for j in (i + 1)..m {
                    v.push(4.0 * lam[i] * lam[j]);
                    let mut e = vec![0.0; m];
                    e[i] = 4.0 * lam[j];
                    e[j] = 4.0 * lam[i];
                    dl.push(e);
                }
            }
        }
        (v, dl)
    }

    fn local_edge_pairs(&self) -> Vec<(usize, usize)> {
        let m = self.mesh.dim + 1;
        let mut out = Vec::new();
        for i in 0..m {
            for j in (i + 1)..m {
                out.push((i, j));
            }
        }
        out
    }

    /// Dofs lying on a boundary facet `(cell, face)` (global ids).
    pub(crate) fn facet_dofs(
        &self,
        cell: usize,
        face: usize,
    ) -> Vec<usize> {
        let dim = self.mesh.dim;
        let cd = &self.cell_dofs[cell];
        let mut out: Vec<usize> = (0..=dim).filter(|&j| j != face).map(|j| cd[j]).collect();
        if self.degree == 2 {
            for (k, (i, j)) in self.local_edge_pairs().into_iter().enumerate() {
                if i != face && j != face {
                    out.push(cd[dim + 1 + k]);
                }
            }
        }
        out
    }

    /// Interpolates a function at the dofs.
    pub fn interpolate(
        &self,
        f: impl Fn([f64; 3]) -> f64,
    ) -> Vec<f64> {
        self.dofs.iter().map(|&x| f(x)).collect()
    }

    /// Evaluates a finite element function at `x` (brute-force point location).
    #[must_use]
    pub fn eval(
        &self,
        u: &[f64],
        x: [f64; 3],
    ) -> Option<f64> {
        let dim = self.mesh.dim;
        for (ci, c) in self.mesh.cells.iter().enumerate() {
            let g = self.mesh.bary_gradients(ci);
            let mut lam = vec![0.0; dim + 1];
            for m in 0..=dim {
                lam[m] = (0..3)
                    .map(|a| g[m][a] * (x[a] - self.mesh.nodes[c[0]][a]))
                    .sum::<f64>()
                    + if m == 0 { 1.0 } else { 0.0 };
            }
            if lam.iter().all(|&l| l > -1e-10) {
                let (b, _) = self.basis(&lam);
                return Some(
                    b.iter()
                        .zip(&self.cell_dofs[ci])
                        .map(|(bi, &d)| bi * u[d])
                        .sum(),
                );
            }
        }
        None
    }

    /// `L2` error against an exact function.
    pub fn l2_error(
        &self,
        u: &[f64],
        exact: impl Fn([f64; 3]) -> f64,
    ) -> f64 {
        let dim = self.mesh.dim;
        let (pts, w) = simplex_rule(dim, self.degree + 2);
        let mut s = 0.0;
        for ci in 0..self.mesh.cells.len() {
            let vol = self.mesh.measure(ci);
            for (p, wq) in pts.iter().zip(&w) {
                let (b, _) = self.basis(p);
                let uh: f64 = b
                    .iter()
                    .zip(&self.cell_dofs[ci])
                    .map(|(bi, &d)| bi * u[d])
                    .sum();
                let mut x = [0.0; 3];
                for (m, &l) in p.iter().enumerate() {
                    for a in 0..3 {
                        x[a] += l * self.mesh.nodes[self.mesh.cells[ci][m]][a];
                    }
                }
                s += vol * wq * (uh - exact(x)).powi(2);
            }
        }
        s.sqrt()
    }
}

/// Problem data `-∇·(k ∇u) + c u = f(t, x)` with marker boundary conditions.
pub struct FemProblem<'a> {
    /// Diffusivity.
    pub k: f64,
    /// Reaction coefficient.
    pub c: f64,
    /// Source `f(t, x)`.
    pub f: &'a (dyn Fn(f64, [f64; 3]) -> f64 + Sync),
    /// Boundary conditions per marker (unmarked facets use the default).
    pub bcs: &'a MarkerBcs<BoundaryCondition>,
}

/// Problem data `-∇·(k(x) ∇u) + c(x) u = f(t, x)` with coefficients that vary in space.
///
/// `k` and `c` are evaluated at every quadrature point of the cells (and `k` at the quadrature
/// points of Neumann / Robin facets); [`FemProblem`] is the constant-coefficient special case with
/// bit-identical results.
pub struct FemProblemVar<'a> {
    /// Diffusivity `k(x)`; it must be positive.
    pub k: &'a (dyn Fn([f64; 3]) -> f64 + Sync),
    /// Reaction coefficient `c(x)`.
    pub c: &'a (dyn Fn([f64; 3]) -> f64 + Sync),
    /// Source `f(t, x)`.
    pub f: &'a (dyn Fn(f64, [f64; 3]) -> f64 + Sync),
    /// Boundary conditions per marker (unmarked facets use the default).
    pub bcs: &'a MarkerBcs<BoundaryCondition>,
}

/// The coefficients the assembly sees.
#[derive(Clone, Copy)]
enum Coefficients<'a> {
    Constant {
        k: f64,
        c: f64,
    },
    Field {
        k: &'a (dyn Fn([f64; 3]) -> f64 + Sync),
        c: &'a (dyn Fn([f64; 3]) -> f64 + Sync),
    },
}

impl Coefficients<'_> {
    fn k(
        &self,
        x: [f64; 3],
    ) -> f64 {
        match self {
            | Self::Constant { k, .. } => *k,
            | Self::Field { k, .. } => k(x),
        }
    }

    fn c(
        &self,
        x: [f64; 3],
    ) -> f64 {
        match self {
            | Self::Constant { c, .. } => *c,
            | Self::Field { c, .. } => c(x),
        }
    }
}

/// Source, boundary data and coefficients of one assembly.
struct Data<'a> {
    f: &'a (dyn Fn(f64, [f64; 3]) -> f64 + Sync),
    bcs: &'a MarkerBcs<BoundaryCondition>,
    coef: Coefficients<'a>,
}

impl<'a> From<&FemProblem<'a>> for Data<'a> {
    fn from(p: &FemProblem<'a>) -> Self {
        Self {
            f: p.f,
            bcs: p.bcs,
            coef: Coefficients::Constant { k: p.k, c: p.c },
        }
    }
}

impl<'a> From<&FemProblemVar<'a>> for Data<'a> {
    fn from(p: &FemProblemVar<'a>) -> Self {
        Self {
            f: p.f,
            bcs: p.bcs,
            coef: Coefficients::Field { k: p.k, c: p.c },
        }
    }
}

/// Cells from which the element matrices are assembled in parallel.
const ASSEMBLY_PAR_CELLS: usize = 64;

struct Assembled {
    k: Csr,
    m: Csr,
    load: Vec<f64>,
    dirichlet: Vec<Option<f64>>,
}

fn assemble(
    space: &FemSpace,
    prob: &Data<'_>,
    t: f64,
    want_mass: bool,
) -> Result<Assembled, String> {
    let dim = space.mesh.dim;
    let nd = space.size();
    let (pts, w) = simplex_rule(dim, space.degree + 1);
    let mut ktrip: Vec<(usize, usize, f64)> = Vec::new();
    let mut mtrip: Vec<(usize, usize, f64)> = Vec::new();
    let mut load = vec![0.0; nd];
    // Element matrices and loads are independent pure computations: they
    // are evaluated in parallel (index order) and scattered sequentially in
    // cell order, so the triplet list and the load are exactly the serial ones.
    let elements: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> =
        crate::kernels::par::map_indexed(space.mesh.cells.len(), ASSEMBLY_PAR_CELLS, |ci| {
            let vol = space.mesh.measure(ci);
            let g = space.mesh.bary_gradients(ci);
            let cd = &space.cell_dofs[ci];
            let nl = cd.len();
            let mut ke = vec![0.0; nl * nl];
            let mut me = vec![0.0; nl * nl];
            let mut fe = vec![0.0; nl];
            for (p, wq) in pts.iter().zip(&w) {
                let (b, dl) = space.basis(p);
                let grads: Vec<[f64; 3]> = dl
                    .iter()
                    .map(|d| {
                        let mut gr = [0.0; 3];
                        for (m, dm) in d.iter().enumerate() {
                            for a in 0..3 {
                                gr[a] += dm * g[m][a];
                            }
                        }
                        gr
                    })
                    .collect();
                let mut x = [0.0; 3];
                for (m, &l) in p.iter().enumerate() {
                    for a in 0..3 {
                        x[a] += l * space.mesh.nodes[space.mesh.cells[ci][m]][a];
                    }
                }
                let fv = (prob.f)(t, x);
                let (kq, cq) = (prob.coef.k(x), prob.coef.c(x));
                for i in 0..nl {
                    fe[i] += vol * wq * fv * b[i];
                    for j in 0..nl {
                        let gg: f64 = (0..3).map(|a| grads[i][a] * grads[j][a]).sum();
                        ke[i * nl + j] += vol * wq * (kq * gg + cq * b[i] * b[j]);
                        me[i * nl + j] += vol * wq * b[i] * b[j];
                    }
                }
            }
            (ke, me, fe)
        });
    for (ci, (ke, me, fe)) in elements.iter().enumerate() {
        let cd = &space.cell_dofs[ci];
        let nl = cd.len();
        for i in 0..nl {
            load[cd[i]] += fe[i];
            for j in 0..nl {
                ktrip.push((cd[i], cd[j], ke[i * nl + j]));
                if want_mass {
                    mtrip.push((cd[i], cd[j], me[i * nl + j]));
                }
            }
        }
    }
    // boundary facets
    let mut dirichlet: Vec<Option<f64>> = vec![None; nd];
    let (fp, fw) = simplex_rule(dim - 1, space.degree + 1);
    for &(ci, face, mk) in &space.mesh.facets {
        let Some(bc) = prob.bcs.get(mk) else { continue };
        let fdofs = space.facet_dofs(ci, face);
        let meas = space.mesh.facet_measure(ci, face);
        let cd = &space.cell_dofs[ci];
        match bc {
            | BoundaryCondition::Dirichlet(gv) => {
                for d in fdofs {
                    dirichlet[d] = Some(gv.eval(t, space.dofs[d]));
                }
            },
            | BoundaryCondition::Periodic | BoundaryCondition::Absorbing => {
                return Err(
                    "FEM boundary conditions: Dirichlet, Neumann, Robin or symmetry".into(),
                );
            },
            | _ => {
                let (alpha, beta, gv) = match bc {
                    | BoundaryCondition::Neumann(g) => (0.0, 1.0, Some(g)),
                    | BoundaryCondition::Robin { alpha, beta, g } => (*alpha, *beta, Some(g)),
                    | _ => (0.0, 1.0, None),
                };
                if beta == 0.0 {
                    return Err("Robin condition with beta = 0: use Dirichlet".into());
                }
                for (q, wq) in fp.iter().zip(&fw) {
                    // facet barycentric -> cell barycentric
                    let mut lam = vec![0.0; dim + 1];
                    let mut idx = 0;
                    for j in 0..=dim {
                        if j != face {
                            lam[j] = q[idx];
                            idx += 1;
                        }
                    }
                    let (b, _) = space.basis(&lam);
                    let mut x = [0.0; 3];
                    for (m, &l) in lam.iter().enumerate() {
                        for a in 0..3 {
                            x[a] += l * space.mesh.nodes[space.mesh.cells[ci][m]][a];
                        }
                    }
                    let gval = gv.map_or(0.0, |g| g.eval(t, x));
                    let kq = prob.coef.k(x);
                    for i in 0..b.len() {
                        if b[i] == 0.0 {
                            continue;
                        }
                        load[cd[i]] += meas * wq * kq * gval / beta * b[i];
                        if alpha != 0.0 {
                            for j in 0..b.len() {
                                ktrip.push((
                                    cd[i],
                                    cd[j],
                                    meas * wq * kq * alpha / beta * b[i] * b[j],
                                ));
                            }
                        }
                    }
                }
            },
        }
    }
    Ok(Assembled {
        k: Csr::from_triplets(nd, nd, &ktrip),
        m: Csr::from_triplets(nd, nd, &mtrip),
        load,
        dirichlet,
    })
}

/// Replaces rows of Dirichlet dofs by identity rows (`shift * M + scale * K` is formed first).
fn with_dirichlet(
    a: &Csr,
    dir: &[Option<f64>],
) -> Csr {
    let mut trip = Vec::with_capacity(a.values.len());
    for r in 0..a.rows {
        if dir[r].is_some() {
            trip.push((r, r, 1.0));
        } else {
            for k in a.indptr[r]..a.indptr[r + 1] {
                trip.push((r, a.indices[k], a.values[k]));
            }
        }
    }
    Csr::from_triplets(a.rows, a.cols, &trip)
}

fn combine(
    m: &Csr,
    k: &Csr,
    cm: f64,
    ck: f64,
) -> Csr {
    let mut trip = Vec::new();
    for r in 0..k.rows {
        for p in k.indptr[r]..k.indptr[r + 1] {
            trip.push((r, k.indices[p], ck * k.values[p]));
        }
        for p in m.indptr[r]..m.indptr[r + 1] {
            trip.push((r, m.indices[p], cm * m.values[p]));
        }
    }
    Csr::from_triplets(k.rows, k.cols, &trip)
}

/// Solves the static problem at time `t` (the source and boundary data are evaluated at `t`).
///
/// # Errors
/// For unsupported boundary conditions or a singular system.
pub fn solve_static(
    space: &FemSpace,
    prob: &FemProblem<'_>,
    t: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    solve_static_data(space, &prob.into(), t, solver)
}

/// [`solve_static`] with spatially varying `k(x)` and `c(x)`.
///
/// # Errors
/// For unsupported boundary conditions or a singular system.
pub fn solve_static_var(
    space: &FemSpace,
    prob: &FemProblemVar<'_>,
    t: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    solve_static_data(space, &prob.into(), t, solver)
}

fn solve_static_data(
    space: &FemSpace,
    prob: &Data<'_>,
    t: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    let asm = assemble(space, prob, t, false)?;
    let a = with_dirichlet(&asm.k, &asm.dirichlet);
    let mut b = asm.load.clone();
    for (i, d) in asm.dirichlet.iter().enumerate() {
        if let Some(v) = d {
            b[i] = *v;
        }
    }
    LinearSystem::new(a, solver)?.solve(&b, None)
}

/// Integrates `u_t - ∇·(k ∇u) + c u = f` with the θ-method (`θ` in `[1/2, 1]`).
///
/// # Errors
/// For unsupported boundary conditions, a singular system or an invalid `theta`.
pub fn solve_heat(
    space: &FemSpace,
    prob: &FemProblem<'_>,
    init: impl Fn([f64; 3]) -> f64,
    t_end: f64,
    steps: usize,
    theta: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    solve_heat_data(space, &prob.into(), init, t_end, steps, theta, solver)
}

/// [`solve_heat`] with spatially varying `k(x)` and `c(x)`.
///
/// # Errors
/// For unsupported boundary conditions, a singular system or an invalid `theta`.
pub fn solve_heat_var(
    space: &FemSpace,
    prob: &FemProblemVar<'_>,
    init: impl Fn([f64; 3]) -> f64,
    t_end: f64,
    steps: usize,
    theta: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    solve_heat_data(space, &prob.into(), init, t_end, steps, theta, solver)
}

fn solve_heat_data(
    space: &FemSpace,
    prob: &Data<'_>,
    init: impl Fn([f64; 3]) -> f64,
    t_end: f64,
    steps: usize,
    theta: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    if !(0.5..=1.0).contains(&theta) || steps == 0 {
        return Err("theta must be in [1/2, 1] and steps > 0".into());
    }
    let dt = t_end / steps as f64;
    let mut u = space.interpolate(init);
    let a0 = assemble(space, prob, 0.0, true)?;
    let lhs = with_dirichlet(&combine(&a0.m, &a0.k, 1.0, dt * theta), &a0.dirichlet);
    let sys = LinearSystem::new(lhs, solver)?;
    let rhs_mat = combine(&a0.m, &a0.k, 1.0, -dt * (1.0 - theta));
    let mut prev_load = a0.load;
    let mut tmp = vec![0.0; u.len()];
    for s in 0..steps {
        let t1 = (s + 1) as f64 * dt;
        let a1 = assemble(space, prob, t1, false)?;
        rhs_mat.matvec(&u, &mut tmp);
        let mut b: Vec<f64> = (0..u.len())
            .map(|i| tmp[i] + dt * (theta * a1.load[i] + (1.0 - theta) * prev_load[i]))
            .collect();
        for (i, d) in a1.dirichlet.iter().enumerate() {
            if let Some(v) = d {
                b[i] = *v;
            }
        }
        u = sys.solve(&b, Some(&u))?;
        prev_load = a1.load;
    }
    Ok(u)
}

/// Residual error indicators `η_T`, one per cell, of a P1 or P2 solution.
///
/// `η_T² = h_T² ∫_T (f(t,·) - c u_h + ∇·(k ∇u_h))² dx
///       + Σ_{e ⊂ ∂T} ½ h_e ∫_e [k ∂_n u_h]² ds`,
///
/// with `h_T` and `h_e` the diameters of the cell and of its faces. The face
/// sum runs over interior faces only, so the estimator measures the
/// discretisation error of the differential equation, not the approximation
/// of non-homogeneous boundary data. For P1 the element term is `f - c u_h`
/// (`∇·(k ∇u_h) = 0`); for P2 the constant `k Δu_h` is included. In 1D the
/// face term is empty (the jump lives on a set of measure zero).
#[must_use]
pub fn error_indicators(
    space: &FemSpace,
    u: &[f64],
    prob: &FemProblem<'_>,
    t: f64,
) -> Vec<f64> {
    error_indicators_data(space, u, &prob.into(), t)
}

/// [`error_indicators`] for spatially varying `k(x)`, `c(x)`:
///
/// `η_T² = h_T² ∫_T (f - c(x) u_h + ∇·(k ∇u_h))² + Σ_e ½ h_e ∫_e [k(x) ∂_n u_h]²`,
///
/// with `∇·(k ∇u_h) = k Δu_h + ∇k · ∇u_h` (the gradient of `k` by central differences of
/// the closure at the quadrature points, nominal step `10⁻⁶ h_T`) and `k`
/// evaluated at the quadrature points of the face in the jump. At floating-
/// point limits the difference uses distinct finite representable coordinates
/// and their actual spacing.
#[must_use]
pub fn error_indicators_var(
    space: &FemSpace,
    u: &[f64],
    prob: &FemProblemVar<'_>,
    t: f64,
) -> Vec<f64> {
    error_indicators_data(space, u, &prob.into(), t)
}

/// A signed squared-indicator contribution with an independent binary exponent.
/// Keeping the contribution squared permits the signed quadrature weights used
/// by the simplex rule and preserves cancellation before the final square root.
#[derive(Clone, Copy, Debug)]
struct ResidualSquare {
    value: f64,
    exponent: i32,
}

impl ResidualSquare {
    const fn direct(value: f64) -> Self {
        Self { value, exponent: 0 }
    }

    fn times(
        self,
        factor: f64,
    ) -> Self {
        let value = self.value * factor;
        if self.value == 0.0 || factor == 0.0 || !self.value.is_finite() || !factor.is_finite() {
            return Self::direct(value);
        }
        if self.exponent == 0 && value.is_normal() {
            return Self::direct(value);
        }
        let (left, left_exponent) = split_residual_scale(self.value);
        let (right, right_exponent) = split_residual_scale(factor);
        Self {
            value: left * right,
            exponent: self.exponent + left_exponent + right_exponent,
        }
    }

    fn product(
        self,
        other: Self,
    ) -> Self {
        let mut product = self.times(other.value);
        product.exponent += other.exponent;
        product
    }

    fn plus(
        self,
        other: Self,
    ) -> Self {
        if self.value == 0.0 {
            return other;
        }
        if other.value == 0.0 {
            return self;
        }
        let sum = self.value + other.value;
        if !self.value.is_finite()
            || !other.value.is_finite()
            || (self.exponent == 0 && other.exponent == 0 && sum.is_finite())
        {
            return Self::direct(sum);
        }
        let (left, left_exponent) = split_residual_scale(self.value);
        let (right, right_exponent) = split_residual_scale(other.value);
        let left_exponent = self.exponent + left_exponent;
        let right_exponent = other.exponent + right_exponent;
        let exponent = left_exponent.max(right_exponent);
        let rescale = |value, delta| {
            if delta < -1022 {
                0.0
            } else {
                value * residual_power_of_two(delta)
            }
        };
        Self {
            value: rescale(left, left_exponent - exponent)
                + rescale(right, right_exponent - exponent),
            exponent,
        }
    }

    fn sqrt(self) -> f64 {
        if self.exponent == 0 || !(self.value > 0.0 && self.value.is_finite()) {
            return self.value.sqrt();
        }
        let remainder = self.exponent.rem_euclid(2);
        let root = (self.value * f64::from(remainder + 1)).sqrt();
        let (significand, shift) = split_residual_scale(root);
        match self.exponent.div_euclid(2) + shift {
            | 1024.. => f64::INFINITY,
            | exponent @ -1022..=1023 => significand * residual_power_of_two(exponent),
            | exponent @ -1074..=-1023 => significand * f64::from_bits(1_u64 << (exponent + 1074)),
            | -1075 => (significand * 0.5) * f64::from_bits(1),
            | _ => 0.0,
        }
    }
}

const fn residual_power_of_two(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

/// Write a finite nonzero value as a signed significand and binary exponent.
fn split_residual_scale(value: f64) -> (f64, i32) {
    let (normal, offset) = if value.is_subnormal() {
        (value * 18_014_398_509_481_984.0, -54)
    } else {
        (value, 0)
    };
    let bits = normal.to_bits();
    (
        f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000),
        ((bits >> 52) & 0x7ff) as i32 - 1023 + offset,
    )
}

/// The degree-five rule uses at most 15 nodes on a simplex of dimension 3.
const RESIDUAL_RULE_CAPACITY: usize = 15;
const MAX_RESIDUAL_DOFS: usize = 10;

fn integrate_residual_square(
    weights: &[f64],
    mut sample: impl FnMut(usize) -> f64,
) -> ResidualSquare {
    let mut values = [0.0; RESIDUAL_RULE_CAPACITY];
    let mut integral = 0.0;
    let mut maximum = 0.0_f64;
    for (index, &weight) in weights.iter().enumerate() {
        let value = sample(index);
        values[index] = value;
        integral += weight * value * value;
        maximum = maximum.max(value.abs());
    }
    let moderate = (residual_power_of_two(-500)..=residual_power_of_two(500)).contains(&maximum);
    if maximum == 0.0 || !maximum.is_finite() || (moderate && integral.is_finite()) {
        return ResidualSquare::direct(integral);
    }
    let normalized = weights
        .iter()
        .zip(values)
        .map(|(weight, value)| {
            let r = value / maximum;
            weight * r * r
        })
        .sum();
    ResidualSquare::direct(normalized)
        .times(maximum)
        .times(maximum)
}

#[derive(Clone, Debug)]
struct ResidualMap {
    origin: [f64; 3],
    edges: [[f64; 3]; 3],
}

impl ResidualMap {
    fn new(vertices: &[[f64; 3]]) -> Result<Self, String> {
        let origin = vertices[0];
        let mut edges = [[0.0; 3]; 3];
        for (edge, node) in edges.iter_mut().zip(&vertices[1..]) {
            for a in 0..3 {
                edge[a] = node[a] - origin[a];
            }
            if edge.iter().any(|value| !value.is_finite()) {
                return Err("residual estimator requires finite simplex edge differences".into());
            }
        }
        Ok(Self { origin, edges })
    }

    fn point(
        &self,
        barycentric: &[f64],
    ) -> [f64; 3] {
        // Form the local offset before adding the origin. This retains the
        // affine coordinates when the physical mesh is far from zero.
        let mut offset = [0.0; 3];
        for (edge, &weight) in self.edges.iter().zip(&barycentric[1..]) {
            for a in 0..3 {
                offset[a] += weight * edge[a];
            }
        }
        [
            self.origin[0] + offset[0],
            self.origin[1] + offset[1],
            self.origin[2] + offset[2],
        ]
    }
}

fn residual_diameter(vertices: &[[f64; 3]]) -> f64 {
    let mut diameter = 0.0_f64;
    for (i, a) in vertices.iter().enumerate() {
        for b in &vertices[i + 1..] {
            let distance = (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
            diameter = diameter.max(distance);
        }
    }
    diameter
}

#[derive(Clone, Debug)]
struct ResidualCell {
    dofs: [usize; MAX_RESIDUAL_DOFS],
    map: ResidualMap,
    gradients: [[f64; 3]; 4],
    volume_weight: ResidualSquare,
    difference_step: f64,
}

#[derive(Clone, Copy, Debug)]
struct FaceSide {
    cell: usize,
    slots: [usize; 3],
}

impl FaceSide {
    fn barycentric(
        &self,
        weights: &[f64],
    ) -> [f64; 4] {
        let mut lambda = [0.0; 4];
        for (&slot, &weight) in self.slots.iter().zip(weights) {
            lambda[slot] = weight;
        }
        lambda
    }
}

#[derive(Clone, Debug)]
struct ResidualFace {
    left: FaceSide,
    right: FaceSide,
    map: ResidualMap,
    normal: [f64; 3],
    weight: ResidualSquare,
}

#[allow(clippy::single_call_fn)]
fn residual_cell_gradients(
    dim: usize,
    map: &ResidualMap,
) -> Result<([[f64; 3]; 4], ResidualSquare), String> {
    let mut jac = [[0.0; 3]; 3];
    for (a, row) in jac.iter_mut().enumerate() {
        if a < dim {
            for (k, value) in row.iter_mut().enumerate().take(dim) {
                *value = map.edges[k][a];
            }
        } else {
            row[a] = 1.0;
        }
    }
    let (mut inverse, mut determinant) = inv3(jac);
    let mut scale = 1.0;
    if !determinant.is_normal() || inverse.iter().flatten().any(|v| !v.is_finite()) {
        scale = jac[..dim]
            .iter()
            .flat_map(|row| row[..dim].iter())
            .fold(0.0_f64, |max, value| max.max(value.abs()));
        if scale == 0.0 {
            return Err("degenerate cell in residual estimator".into());
        }
        for row in &mut jac[..dim] {
            for value in &mut row[..dim] {
                *value /= scale;
            }
        }
        (inverse, determinant) = inv3(jac);
    }
    if determinant == 0.0 || !determinant.is_finite() {
        return Err("degenerate cell in residual estimator".into());
    }
    let mut gradients = [[0.0; 3]; 4];
    for k in 0..dim {
        for a in 0..dim {
            gradients[k + 1][a] = inverse[k][a] / scale;
            gradients[0][a] -= gradients[k + 1][a];
        }
    }
    if gradients.iter().flatten().any(|value| !value.is_finite()) {
        return Err("cell gradients exceed the residual estimator's numeric range".into());
    }
    let factorial: f64 = (1..=dim).map(|n| n as f64).product();
    let measure = (0..dim).fold(
        ResidualSquare::direct(determinant.abs()).times(factorial.recip()),
        |volume, _| volume.times(scale),
    );
    Ok((gradients, measure))
}

#[allow(clippy::single_call_fn)]
fn residual_face_geometry(
    dim: usize,
    map: &ResidualMap,
) -> Result<([f64; 3], ResidualSquare), String> {
    let a = map.edges[0];
    if dim == 2 {
        let length = a[0].hypot(a[1]);
        if length == 0.0 || !length.is_finite() {
            return Err("degenerate edge in residual estimator".into());
        }
        return Ok((
            [-a[1] / length, a[0] / length, 0.0],
            ResidualSquare::direct(length),
        ));
    }

    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let maximum = |v: [f64; 3]| v.iter().fold(0.0_f64, |max, x| max.max(x.abs()));
    let b = map.edges[1];
    let mut direction = cross(a, b);
    let mut edge_scale = 1.0;
    if !maximum(direction).is_normal() || direction.iter().any(|value| !value.is_finite()) {
        edge_scale = maximum(a).max(maximum(b));
        if edge_scale == 0.0 {
            return Err("degenerate face in residual estimator".into());
        }
        direction = cross(a.map(|x| x / edge_scale), b.map(|x| x / edge_scale));
    }
    let normal_scale = maximum(direction);
    if normal_scale == 0.0 || !normal_scale.is_finite() {
        return Err("degenerate face in residual estimator".into());
    }
    let direction = direction.map(|x| x / normal_scale);
    let length = direction[0].hypot(direction[1]).hypot(direction[2]);
    let measure = ResidualSquare::direct(0.5 * length)
        .times(normal_scale)
        .times(edge_scale)
        .times(edge_scale);
    Ok((direction.map(|x| x / length), measure))
}

impl ResidualCell {
    fn local_values(
        &self,
        u: &[f64],
        count: usize,
    ) -> [f64; MAX_RESIDUAL_DOFS] {
        let mut values = [0.0; MAX_RESIDUAL_DOFS];
        for (value, &dof) in values[..count].iter_mut().zip(&self.dofs) {
            *value = u[dof];
        }
        values
    }

    fn gradient(
        &self,
        degree: usize,
        dim: usize,
        values: &[f64; MAX_RESIDUAL_DOFS],
        lambda: &[f64],
    ) -> [f64; 3] {
        let mut coefficients = [0.0; 4];
        if degree == 1 {
            coefficients[..=dim].copy_from_slice(&values[..=dim]);
        } else {
            for i in 0..=dim {
                coefficients[i] = values[i] * (4.0 * lambda[i] - 1.0);
            }
            let mut edge = dim + 1;
            for i in 0..=dim {
                for j in i + 1..=dim {
                    coefficients[i] += values[edge] * (4.0 * lambda[j]);
                    coefficients[j] += values[edge] * (4.0 * lambda[i]);
                    edge += 1;
                }
            }
        }
        let mut gradient = [0.0; 3];
        for (coefficient, basis_gradient) in coefficients[..=dim].iter().zip(&self.gradients) {
            for a in 0..dim {
                gradient[a] += coefficient * basis_gradient[a];
            }
        }
        gradient
    }

    fn laplacian(
        &self,
        degree: usize,
        dim: usize,
        values: &[f64; MAX_RESIDUAL_DOFS],
    ) -> f64 {
        if degree == 1 {
            return 0.0;
        }
        let dot = |i: usize, j: usize| {
            (0..dim)
                .map(|a| self.gradients[i][a] * self.gradients[j][a])
                .sum::<f64>()
        };
        let mut laplacian = 0.0;
        for (i, &value) in values.iter().enumerate().take(dim + 1) {
            laplacian += 4.0 * value * dot(i, i);
        }
        let mut edge = dim + 1;
        for i in 0..=dim {
            for j in i + 1..=dim {
                laplacian += 8.0 * values[edge] * dot(i, j);
                edge += 1;
            }
        }
        laplacian
    }
}

fn perturb_coordinate(
    value: f64,
    delta: f64,
) -> f64 {
    let candidate = value + delta;
    let distinct = if crate::kernels::exact_eq(candidate, value) {
        if delta.is_sign_negative() {
            value.next_down()
        } else {
            value.next_up()
        }
    } else {
        candidate
    };
    if distinct.is_finite() {
        distinct
    } else {
        value
    }
}

/// Prepared residual indicators for repeated estimates on a fixed P1 or P2 mesh.
///
/// Owns the required geometry, quadrature and degree-of-freedom mapping from
/// the supplied space. Solution values are supplied in the original order.
/// Construct a new estimator after changing or refining the space.
#[derive(Clone, Debug)]
pub struct SimplexResidualEstimator {
    dim: usize,
    degree: usize,
    size: usize,
    local_dofs: usize,
    cells: Vec<ResidualCell>,
    faces: Vec<ResidualFace>,
    cell_points: Vec<Vec<f64>>,
    cell_weights: Vec<f64>,
    basis_values: Vec<Vec<f64>>,
    face_points: Vec<Vec<f64>>,
    face_weights: Vec<f64>,
}

impl SimplexResidualEstimator {
    /// Prepares cell geometry and interior face adjacency.
    ///
    /// # Errors
    /// Returns an error for an invalid space, degenerate cells, or geometry
    /// whose edge differences, diameters or barycentric gradients cannot be
    /// represented in f64.
    pub fn new(space: &FemSpace) -> Result<Self, String> {
        let dim = space.mesh.dim;
        if !(1..=3).contains(&dim) || !(1..=2).contains(&space.degree) {
            return Err(
                "residual estimator requires a P1 or P2 simplex space in dimensions 1 to 3".into(),
            );
        }
        if space.cell_dofs.len() != space.mesh.cells.len() {
            return Err("cell degree-of-freedom map does not match the mesh".into());
        }
        let local_dofs = if space.degree == 1 {
            dim + 1
        } else {
            (dim + 1) * (dim + 2) / 2
        };
        let mut cells = Vec::with_capacity(space.mesh.cells.len());
        for (cell, mapping) in space.mesh.cells.iter().zip(&space.cell_dofs) {
            if cell.len() != dim + 1
                || mapping.len() != local_dofs
                || mapping.iter().any(|&dof| dof >= space.size())
            {
                return Err("invalid residual-estimator cell or degree-of-freedom map".into());
            }
            let mut vertices = [[0.0; 3]; 4];
            for (i, &vertex) in cell.iter().enumerate() {
                vertices[i] = *space
                    .mesh
                    .nodes
                    .get(vertex)
                    .ok_or_else(|| "invalid cell vertex in residual estimator".to_string())?;
                if vertices[i].iter().any(|value| !value.is_finite()) {
                    return Err("nonfinite cell vertex in residual estimator".into());
                }
            }
            let vertices = &vertices[..=dim];
            let map = ResidualMap::new(vertices)?;
            let diameter = residual_diameter(vertices);
            if diameter == 0.0 || !diameter.is_finite() {
                return Err("invalid cell diameter in residual estimator".into());
            }
            let (gradients, measure) = residual_cell_gradients(dim, &map)?;
            let mut dofs = [0; MAX_RESIDUAL_DOFS];
            dofs[..local_dofs].copy_from_slice(mapping);
            cells.push(ResidualCell {
                dofs,
                map,
                gradients,
                volume_weight: measure.times(diameter).times(diameter),
                difference_step: 1e-6 * diameter,
            });
        }
        let faces = residual_faces(&space.mesh)?;
        let (cell_points, cell_weights) = simplex_rule(dim, 2);
        let basis_values = cell_points
            .iter()
            .map(|lambda| space.basis(lambda).0)
            .collect();
        let (face_points, face_weights) = if dim > 1 {
            simplex_rule(dim - 1, 2)
        } else {
            (Vec::new(), Vec::new())
        };
        Ok(Self {
            dim,
            degree: space.degree,
            size: space.size(),
            local_dofs,
            cells,
            faces,
            cell_points,
            cell_weights,
            basis_values,
            face_points,
            face_weights,
        })
    }

    pub(crate) fn heap_bytes(&self) -> usize {
        let rules = [&self.cell_points, &self.basis_values, &self.face_points];
        let nested: usize = rules
            .into_iter()
            .map(|rows| {
                rows.capacity() * size_of::<Vec<f64>>()
                    + rows
                        .iter()
                        .map(|row| row.capacity() * size_of::<f64>())
                        .sum::<usize>()
            })
            .sum();
        self.cells.capacity() * size_of::<ResidualCell>()
            + self.faces.capacity() * size_of::<ResidualFace>()
            + (self.cell_weights.capacity() + self.face_weights.capacity()) * size_of::<f64>()
            + nested
    }

    /// Computes the cell indicators for constant coefficients.
    ///
    /// # Errors
    /// Returns an error when the solution length differs from the space size.
    pub fn indicators(
        &self,
        u: &[f64],
        problem: &FemProblem<'_>,
        time: f64,
    ) -> Result<Vec<f64>, String> {
        self.check_solution(u)?;
        Ok(self.evaluate(u, &problem.into(), time))
    }

    /// Computes the cell indicators for spatially varying coefficients.
    ///
    /// # Errors
    /// Returns an error when the solution length differs from the space size.
    pub fn indicators_var(
        &self,
        u: &[f64],
        problem: &FemProblemVar<'_>,
        time: f64,
    ) -> Result<Vec<f64>, String> {
        self.check_solution(u)?;
        Ok(self.evaluate(u, &problem.into(), time))
    }

    fn check_solution(
        &self,
        u: &[f64],
    ) -> Result<(), String> {
        if u.len() != self.size {
            return Err("residual estimator solution length differs from the space size".into());
        }
        Ok(())
    }

    fn evaluate(
        &self,
        u: &[f64],
        problem: &Data<'_>,
        time: f64,
    ) -> Vec<f64> {
        let mut squares =
            crate::kernels::par::map_indexed(self.cells.len(), ASSEMBLY_PAR_CELLS, |ci| {
                let cell = &self.cells[ci];
                let values = cell.local_values(u, self.local_dofs);
                let laplacian = cell.laplacian(self.degree, self.dim, &values);
                let residual = |q: usize| {
                    let lambda: &[f64] = &self.cell_points[q];
                    let x = cell.map.point(lambda);
                    let value: f64 = self.basis_values[q]
                        .iter()
                        .zip(values)
                        .map(|(b, u)| b * u)
                        .sum();
                    let mut residual = (problem.f)(time, x) - problem.coef.c(x) * value
                        + problem.coef.k(x) * laplacian;
                    if let Coefficients::Field { k, .. } = problem.coef {
                        let gradient = cell.gradient(self.degree, self.dim, &values, lambda);
                        for a in 0..self.dim {
                            let (mut plus, mut minus) = (x, x);
                            plus[a] = perturb_coordinate(x[a], cell.difference_step);
                            minus[a] = perturb_coordinate(x[a], -cell.difference_step);
                            residual += (k(plus) - k(minus)) / (plus[a] - minus[a]) * gradient[a];
                        }
                    }
                    residual
                };
                cell.volume_weight
                    .product(integrate_residual_square(&self.cell_weights, residual))
            });
        for face in &self.faces {
            let left = &self.cells[face.left.cell];
            let right = &self.cells[face.right.cell];
            let left_values = left.local_values(u, self.local_dofs);
            let right_values = right.local_values(u, self.local_dofs);
            let jump = |q: usize| {
                let lambda: &[f64] = &self.face_points[q];
                let x = face.map.point(lambda);
                let left_lambda = face.left.barycentric(lambda);
                let right_lambda = face.right.barycentric(lambda);
                let a = left.gradient(self.degree, self.dim, &left_values, &left_lambda);
                let b = right.gradient(self.degree, self.dim, &right_values, &right_lambda);
                problem.coef.k(x)
                    * (0..self.dim)
                        .map(|i| (b[i] - a[i]) * face.normal[i])
                        .sum::<f64>()
            };
            let share = face
                .weight
                .product(integrate_residual_square(&self.face_weights, jump));
            squares[face.left.cell] = squares[face.left.cell].plus(share);
            squares[face.right.cell] = squares[face.right.cell].plus(share);
        }
        squares.into_iter().map(ResidualSquare::sqrt).collect()
    }
}

#[allow(clippy::single_call_fn)]
fn residual_faces(mesh: &SimplexMesh) -> Result<Vec<ResidualFace>, String> {
    let dim = mesh.dim;
    if dim == 1 {
        return Ok(Vec::new());
    }
    let mut adjacency = BTreeMap::<[usize; 3], [Option<FaceSide>; 3]>::new();
    for (cell, vertices) in mesh.cells.iter().enumerate() {
        for face in 0..=dim {
            let mut indexed = [(0, 0); 3];
            for (entry, (local, &vertex)) in indexed.iter_mut().zip(
                vertices
                    .iter()
                    .enumerate()
                    .filter(|(local, _)| *local != face),
            ) {
                *entry = (vertex, local);
            }
            indexed[..dim].sort_unstable_by_key(|&(vertex, _)| vertex);
            let key = indexed.map(|(vertex, _)| vertex);
            let side = FaceSide {
                cell,
                slots: indexed.map(|(_, local)| local),
            };
            let incident = adjacency.entry(key).or_default();
            if let Some(next) = incident.iter_mut().find(|entry| entry.is_none()) {
                *next = Some(side);
            }
        }
    }

    let mut faces = Vec::new();
    for (key, incident) in adjacency {
        let [Some(left), Some(right), None] = incident else {
            continue;
        };
        let mut vertices = [[0.0; 3]; 3];
        for (position, &vertex) in vertices.iter_mut().zip(&key[..dim]) {
            *position = mesh.nodes[vertex];
        }
        let map = ResidualMap::new(&vertices[..dim])?;
        let diameter = residual_diameter(&vertices[..dim]);
        let (normal, measure) = residual_face_geometry(dim, &map)?;
        faces.push(ResidualFace {
            left,
            right,
            map,
            normal,
            weight: measure.times(diameter).times(0.5),
        });
    }
    Ok(faces)
}

fn error_indicators_data(
    space: &FemSpace,
    u: &[f64],
    prob: &Data<'_>,
    t: f64,
) -> Vec<f64> {
    match SimplexResidualEstimator::new(space) {
        | Ok(estimator) if u.len() >= estimator.size => estimator.evaluate(u, prob, t),
        | _ => vec![f64::NAN; space.mesh.cells.len()],
    }
}

/// The final space, solution and history of [`solve_adaptive`].
///
/// The history holds `(degrees of freedom, (Σ η_T²)^(1/2))` after every solve;
/// the first entry is the initial mesh.
pub type Adaptive = (FemSpace, Vec<f64>, Vec<(usize, f64)>);

/// Adaptive solution of the static problem: solve, estimate, mark and refine.
///
/// Starting from `mesh`, the loop solves `-∇·(k ∇u) + c u = f(t, x)` with
/// [`solve_static`], estimates the error with [`error_indicators`], marks the
/// smallest set of cells whose indicators carry the fraction `theta` of
/// `Σ η_T²` (Dörfler marking), refines them with
/// [`SimplexMesh::refine_marked`] and repeats. It stops when the total
/// estimate is at most `tolerance` or after `max_steps` refinements (so at
/// most `max_steps + 1` solves). `degree` selects P1 or P2.
///
/// # Errors
/// If `degree` is not 1 or 2, `theta` is outside `(0, 1]`, `max_steps` is
/// zero, or one of the solves fails.
#[allow(clippy::too_many_arguments)] // mirrors `solve_static` plus the adaptive controls
pub fn solve_adaptive(
    mesh: &SimplexMesh,
    prob: &FemProblem<'_>,
    degree: usize,
    t: f64,
    theta: f64,
    tolerance: f64,
    max_steps: usize,
    solver: LinearSolver,
) -> Result<Adaptive, String> {
    solve_adaptive_data(
        mesh,
        &prob.into(),
        degree,
        t,
        theta,
        tolerance,
        max_steps,
        solver,
    )
}

/// [`solve_adaptive`] with spatially varying `k(x)`, `c(x)` ([`error_indicators_var`]).
///
/// # Errors
/// As [`solve_adaptive`].
#[allow(clippy::too_many_arguments)] // same controls as `solve_adaptive`
pub fn solve_adaptive_var(
    mesh: &SimplexMesh,
    prob: &FemProblemVar<'_>,
    degree: usize,
    t: f64,
    theta: f64,
    tolerance: f64,
    max_steps: usize,
    solver: LinearSolver,
) -> Result<Adaptive, String> {
    solve_adaptive_data(
        mesh,
        &prob.into(),
        degree,
        t,
        theta,
        tolerance,
        max_steps,
        solver,
    )
}

#[allow(clippy::too_many_arguments)] // same controls as `solve_adaptive`
fn solve_adaptive_data(
    mesh: &SimplexMesh,
    prob: &Data<'_>,
    degree: usize,
    t: f64,
    theta: f64,
    tolerance: f64,
    max_steps: usize,
    solver: LinearSolver,
) -> Result<Adaptive, String> {
    if max_steps == 0 {
        return Err("max_steps must be positive".into());
    }
    if !(theta > 0.0 && theta <= 1.0) {
        return Err("theta must be in (0, 1]".into());
    }
    let mut mesh = mesh.clone();
    let mut history = Vec::with_capacity(max_steps + 1);
    loop {
        let space = FemSpace::new(mesh, degree)?;
        let u = solve_static_data(&space, prob, t, solver)?;
        let eta = error_indicators_data(&space, &u, prob, t);
        let energy: f64 = eta.iter().map(|e| e * e).sum();
        let total = energy.sqrt();
        history.push((space.size(), total));
        if total <= tolerance || history.len() > max_steps {
            return Ok((space, u, history));
        }
        let mut order: Vec<usize> = (0..eta.len()).collect();
        order.sort_by(|&i, &j| eta[j].total_cmp(&eta[i]));
        let mut marked = vec![false; eta.len()];
        let mut acc = 0.0;
        for i in order {
            if acc >= theta * energy {
                break;
            }
            marked[i] = true;
            acc += eta[i] * eta[i];
        }
        mesh = space.mesh.refine_marked(&marked);
    }
}

/// Advance nodal finite-element heat data by one theta-method step.
///
/// Uses the existing consistent mass/stiffness assembly and boundary load
/// machinery at both physical times; no reinterpolation of the state occurs.
///
/// # Errors
/// Invalid state/control data, incompatible boundary conditions, or failed linear solve.
#[allow(clippy::too_many_arguments)] // Explicit restart state and physical time.
pub fn heat_step(
    space: &FemSpace,
    prob: &FemProblem<'_>,
    state: &[f64],
    t: f64,
    dt: f64,
    theta: f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    heat_step_source_data(space, prob, state, t, dt, theta, solver, None)
}

/// Advance FEM heat data with additional nodal source values at both times.
///
/// The source is interpolated in the existing finite-element basis and mapped
/// to a load by the consistent mass matrix, preserving the theta method.
///
/// # Errors
/// Invalid nodal source/state data, controls, boundary data, or linear solve.
#[allow(clippy::too_many_arguments)]
pub fn heat_step_nodal_source(
    space: &FemSpace,
    prob: &FemProblem<'_>,
    state: &[f64],
    t: f64,
    dt: f64,
    theta: f64,
    solver: LinearSolver,
    before: &[f64],
    after: &[f64],
) -> Result<Vec<f64>, String> {
    if before.len() != state.len()
        || after.len() != state.len()
        || before.iter().chain(after).any(|v| !v.is_finite())
    {
        return Err("FEM source must contain one finite rate per nodal degree of freedom".into());
    }
    heat_step_source_data(
        space,
        prob,
        state,
        t,
        dt,
        theta,
        solver,
        Some((before, after)),
    )
}
#[allow(clippy::too_many_arguments)]
fn heat_step_source_data(
    space: &FemSpace,
    prob: &FemProblem<'_>,
    state: &[f64],
    t: f64,
    dt: f64,
    theta: f64,
    solver: LinearSolver,
    nodal: Option<(&[f64], &[f64])>,
) -> Result<Vec<f64>, String> {
    if state.len() != space.dofs.len()
        || state.iter().any(|v| !v.is_finite())
        || !t.is_finite()
        || !dt.is_finite()
        || dt <= 0.0
        || !(0.5..=1.0).contains(&theta)
    {
        return Err("invalid FEM heat step".into());
    }
    let mut before = assemble(space, &prob.into(), t, true)?;
    for (actual, prescribed) in state.iter().zip(&before.dirichlet) {
        if let Some(value) = prescribed
            && (actual - value).abs() > 64.0 * f64::EPSILON * actual.abs().max(value.abs())
        {
            return Err("FEM initial/restart state conflicts with Dirichlet boundary data".into());
        }
    }
    let mut after = assemble(space, &prob.into(), t + dt, false)?;
    if let Some((f0, f1)) = nodal {
        let mut load = vec![0.0; state.len()];
        before.m.matvec(f0, &mut load);
        for (value, source) in before.load.iter_mut().zip(&load) {
            *value += source;
        }
        before.m.matvec(f1, &mut load);
        for (value, source) in after.load.iter_mut().zip(load) {
            *value += source;
        }
    }
    let lhs = with_dirichlet(
        &combine(&before.m, &after.k, 1.0, dt * theta),
        &after.dirichlet,
    );
    let rhs_operator = combine(&before.m, &before.k, 1.0, -dt * (1.0 - theta));
    let mut rhs = vec![0.0; state.len()];
    rhs_operator.matvec(state, &mut rhs);
    for i in 0..state.len() {
        rhs[i] += dt * (theta * after.load[i] + (1.0 - theta) * before.load[i]);
        if let Some(value) = after.dirichlet[i] {
            rhs[i] = value;
        }
    }
    LinearSystem::new(lhs, solver)?.solve(&rhs, Some(state))
}
