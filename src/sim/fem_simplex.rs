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
//!
//! Quadrature is the Grundmann–Möller rule of arbitrary degree
//! ([`simplex_rule`]); linear systems use [`LinearSystem`] (sparse LU or
//! ILU-`BiCGStab`). See [`crate::sim::fem_triangles`] for the adaptive 2D
//! solver with error estimation.
//!
//! Limitations: constant `k` and `c`; no curved boundaries or higher-order
//! geometry; no adaptive refinement in 3D.

use crate::kernels::krylov::Csr;
use crate::sim::boundary::{BoundaryCondition, MarkerBcs};
use crate::sim::structured::{LinearSolver, LinearSystem};
use std::collections::HashMap;

fn compositions(total: usize, parts: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
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
pub fn simplex_rule(n: usize, s: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    if n == 0 {
        return (vec![vec![1.0]], vec![1.0]);
    }
    let d = 2 * s + 1;
    let fact = |m: usize| (1..=m).map(|x| x as f64).product::<f64>();
    let mut pts = Vec::new();
    let mut w = Vec::new();
    for i in 0..=s {
        let denom = (d + n - 2 * i) as f64;
        let coef = (if i % 2 == 0 { 1.0 } else { -1.0 }) * denom.powi(d as i32) / (fact(i) * fact(d + n - i));
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
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
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
    fn finish(dim: usize, nodes: Vec<[f64; 3]>, mut cells: Vec<Vec<usize>>, lo: [f64; 3], hi: [f64; 3]) -> Self {
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
                let mut key: Vec<usize> = c.iter().enumerate().filter(|(j, _)| *j != f).map(|(_, &v)| v).collect();
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
        Self { dim, nodes, cells, facets }
    }

    /// Uniform interval mesh with `n` cells on `[lo, hi]` (markers 0 and 1).
    #[must_use]
    pub fn interval(n: usize, lo: f64, hi: f64) -> Self {
        let nodes: Vec<[f64; 3]> = (0..=n).map(|i| [lo + (hi - lo) * i as f64 / n as f64, 0.0, 0.0]).collect();
        let cells = (0..n).map(|i| vec![i, i + 1]).collect();
        Self::finish(1, nodes, cells, [lo, 0.0, 0.0], [hi, 0.0, 0.0])
    }

    /// Triangulated rectangle (two triangles per quad), markers 0-3.
    #[must_use]
    pub fn rectangle(nx: usize, ny: usize, lo: [f64; 2], hi: [f64; 2]) -> Self {
        let mut nodes = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push([lo[0] + (hi[0] - lo[0]) * i as f64 / nx as f64, lo[1] + (hi[1] - lo[1]) * j as f64 / ny as f64, 0.0]);
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
    pub fn boxed(n: [usize; 3], lo: [f64; 3], hi: [f64; 3]) -> Self {
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
        let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
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

    fn jacobian_of(dim: usize, nodes: &[[f64; 3]], c: &[usize]) -> ([[f64; 3]; 3], f64) {
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
    pub fn measure(&self, cell: usize) -> f64 {
        let (_, d) = Self::jacobian_of(self.dim, &self.nodes, &self.cells[cell]);
        let fact: f64 = (1..=self.dim).map(|x| x as f64).product();
        d.abs() / fact
    }

    /// Constant gradients of the barycentric coordinates of a cell.
    #[must_use]
    pub fn bary_gradients(&self, cell: usize) -> Vec<[f64; 3]> {
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

    pub(crate) fn facet_measure(&self, cell: usize, face: usize) -> f64 {
        let c = &self.cells[cell];
        let v: Vec<[f64; 3]> = c.iter().enumerate().filter(|(j, _)| *j != face).map(|(_, &q)| self.nodes[q]).collect();
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        match self.dim {
            1 => 1.0,
            2 => {
                let d = sub(v[1], v[0]);
                d[0].hypot(d[1])
            },
            _ => {
                let (a, b) = (sub(v[1], v[0]), sub(v[2], v[0]));
                let cr = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
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
    pub fn new(mesh: SimplexMesh, degree: usize) -> Result<Self, String> {
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
                            dofs.push([f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1]), f64::midpoint(a[2], b[2])]);
                            dofs.len() - 1
                        });
                        cell_dofs[ci].push(id);
                    }
                }
            }
        }
        Ok(Self { mesh, degree, dofs, cell_dofs })
    }

    /// Number of degrees of freedom.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.dofs.len()
    }

    /// Basis values and barycentric derivatives at barycentric point `lam`.
    pub(crate) fn basis(&self, lam: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
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
    pub(crate) fn facet_dofs(&self, cell: usize, face: usize) -> Vec<usize> {
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
    pub fn interpolate(&self, f: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
        self.dofs.iter().map(|&x| f(x)).collect()
    }

    /// Evaluates a finite element function at `x` (brute-force point location).
    #[must_use]
    pub fn eval(&self, u: &[f64], x: [f64; 3]) -> Option<f64> {
        let dim = self.mesh.dim;
        for (ci, c) in self.mesh.cells.iter().enumerate() {
            let g = self.mesh.bary_gradients(ci);
            let mut lam = vec![0.0; dim + 1];
            for m in 0..=dim {
                lam[m] = (0..3).map(|a| g[m][a] * (x[a] - self.mesh.nodes[c[0]][a])).sum::<f64>() + if m == 0 { 1.0 } else { 0.0 };
            }
            if lam.iter().all(|&l| l > -1e-10) {
                let (b, _) = self.basis(&lam);
                return Some(b.iter().zip(&self.cell_dofs[ci]).map(|(bi, &d)| bi * u[d]).sum());
            }
        }
        None
    }

    /// `L2` error against an exact function.
    pub fn l2_error(&self, u: &[f64], exact: impl Fn([f64; 3]) -> f64) -> f64 {
        let dim = self.mesh.dim;
        let (pts, w) = simplex_rule(dim, self.degree + 2);
        let mut s = 0.0;
        for ci in 0..self.mesh.cells.len() {
            let vol = self.mesh.measure(ci);
            for (p, wq) in pts.iter().zip(&w) {
                let (b, _) = self.basis(p);
                let uh: f64 = b.iter().zip(&self.cell_dofs[ci]).map(|(bi, &d)| bi * u[d]).sum();
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
    pub f: &'a dyn Fn(f64, [f64; 3]) -> f64,
    /// Boundary conditions per marker (unmarked facets use the default).
    pub bcs: &'a MarkerBcs<BoundaryCondition>,
}

struct Assembled {
    k: Csr,
    m: Csr,
    load: Vec<f64>,
    dirichlet: Vec<Option<f64>>,
}

fn assemble(space: &FemSpace, prob: &FemProblem<'_>, t: f64, want_mass: bool) -> Result<Assembled, String> {
    let dim = space.mesh.dim;
    let nd = space.size();
    let (pts, w) = simplex_rule(dim, space.degree + 1);
    let mut ktrip: Vec<(usize, usize, f64)> = Vec::new();
    let mut mtrip: Vec<(usize, usize, f64)> = Vec::new();
    let mut load = vec![0.0; nd];
    for ci in 0..space.mesh.cells.len() {
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
            for i in 0..nl {
                fe[i] += vol * wq * fv * b[i];
                for j in 0..nl {
                    let gg: f64 = (0..3).map(|a| grads[i][a] * grads[j][a]).sum();
                    ke[i * nl + j] += vol * wq * (prob.k * gg + prob.c * b[i] * b[j]);
                    me[i * nl + j] += vol * wq * b[i] * b[j];
                }
            }
        }
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
            BoundaryCondition::Dirichlet(gv) => {
                for d in fdofs {
                    dirichlet[d] = Some(gv.eval(t, space.dofs[d]));
                }
            },
            BoundaryCondition::Periodic | BoundaryCondition::Absorbing => {
                return Err("FEM boundary conditions: Dirichlet, Neumann, Robin or symmetry".into());
            },
            _ => {
                let (alpha, beta, gv) = match bc {
                    BoundaryCondition::Neumann(g) => (0.0, 1.0, Some(g)),
                    BoundaryCondition::Robin { alpha, beta, g } => (*alpha, *beta, Some(g)),
                    _ => (0.0, 1.0, None),
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
                    for i in 0..b.len() {
                        if b[i] == 0.0 {
                            continue;
                        }
                        load[cd[i]] += meas * wq * prob.k * gval / beta * b[i];
                        if alpha != 0.0 {
                            for j in 0..b.len() {
                                ktrip.push((cd[i], cd[j], meas * wq * prob.k * alpha / beta * b[i] * b[j]));
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
fn with_dirichlet(a: &Csr, dir: &[Option<f64>]) -> Csr {
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

fn combine(m: &Csr, k: &Csr, cm: f64, ck: f64) -> Csr {
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
pub fn solve_static(space: &FemSpace, prob: &FemProblem<'_>, t: f64, solver: LinearSolver) -> Result<Vec<f64>, String> {
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
        let mut b: Vec<f64> = (0..u.len()).map(|i| tmp[i] + dt * (theta * a1.load[i] + (1.0 - theta) * prev_load[i])).collect();
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
