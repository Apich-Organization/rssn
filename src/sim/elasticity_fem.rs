//! Linear elasticity on triangle/tetrahedron meshes (P1/P2).
//!
//! # Linear elasticity on triangle and tetrahedron meshes (P1 and P2)
//!
//! Small-strain isotropic elasticity `-∇·σ(u) = f`, `σ = C ε(u)`, solved with
//! vector Lagrange elements on the simplex meshes of
//! [`crate::sim::fem_simplex`]:
//!
//! | Case | Dimension | Elements | Notes |
//! |---|---|---|---|
//! | plane stress ([`Model::PlaneStress`]) | 2D | P1 (CST), P2 (LST) | thickness factor |
//! | plane strain ([`Model::PlaneStrain`]) | 2D | P1, P2 | |
//! | full 3D ([`Model::Solid3d`]) | 3D | P1, P2 tetrahedra | |
//! | heterogeneous `E(x)`, `ν(x)` ([`ElasticityProblemVar`], [`solve_heterogeneous`]) | 2D, 3D | P1, P2 | all three models; per-cell stress uses the cell material |
//!
//! Boundary conditions per marker ([`ElasticBoundary`]): displacement per
//! component (clamp, roller, prescribed, time dependent), traction vector
//! (`σ n = t`), traction-free. Body force `f(x)`. Stress recovery per cell
//! (exact for P1, centroid value for P2) and von Mises stress are provided.
//!
//! **Heterogeneous material.** [`ElasticityProblemVar`] with
//! [`solve_heterogeneous`] takes the isotropic material as closures `E(x)`,
//! `ν(x)` evaluated at every quadrature point (and at the cell centroid for
//! the per-cell stress recovery, which uses that cell's constitutive matrix),
//! keeping the same P1/P2 spaces, plane stress / plane strain / 3D models,
//! boundary markers and body force. A constant material reproduces [`solve`]
//! to round-off, and a material that is piecewise constant on the element
//! boundaries (layered media) passes the patch test exactly.
//!
//! Large deformation (St. Venant–Kirchhoff, neo-Hookean), elastodynamics
//! (Newmark, generalised-α) and penalty contact are in
//! [`hyperelastic`](crate::sim::hyperelastic).
//!
//! Limitations: the material is isotropic (heterogeneous through
//! [`ElasticityProblemVar`], but no tensorial/anisotropic elasticity, no
//! per-cell damage or plasticity), and the analysis is linearised, static and
//! small strain (no prestress, contact or dynamics).

use crate::kernels::krylov::Csr;
use crate::sim::boundary::ElasticBoundary;
use crate::sim::boundary::MarkerBcs;
use crate::sim::fem_simplex::FemSpace;
use crate::sim::fem_simplex::simplex_rule;
use crate::sim::structured::LinearSolver;
use crate::sim::structured::LinearSystem;

/// Constitutive model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    /// 2D plane stress.
    PlaneStress,
    /// 2D plane strain.
    PlaneStrain,
    /// 3D solid.
    Solid3d,
}

/// Isotropic material.
#[derive(Clone, Copy, Debug)]
pub struct Material {
    /// Young's modulus.
    pub young: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// Thickness (2D only).
    pub thickness: f64,
}

/// Body force closure `f(x)`.
pub type BodyForce<'a> = &'a dyn Fn([f64; 3]) -> [f64; 3];

/// Elasticity problem.
pub struct ElasticityProblem<'a> {
    /// Model.
    pub model: Model,
    /// Material.
    pub material: Material,
    /// Body force.
    pub force: BodyForce<'a>,
    /// Boundary conditions per marker (unlisted: traction free).
    pub bcs: &'a MarkerBcs<ElasticBoundary>,
}

/// Elasticity problem with a spatially varying isotropic material, `E(x)` and
/// `ν(x)` (the heterogeneous counterpart of [`ElasticityProblem`]; a constant
/// material reproduces [`solve`] to round-off).
pub struct ElasticityProblemVar<'a> {
    /// Model.
    pub model: Model,
    /// Young's modulus `E(x)` (positive).
    pub young: &'a dyn Fn([f64; 3]) -> f64,
    /// Poisson ratio `ν(x)`.
    pub poisson: &'a dyn Fn([f64; 3]) -> f64,
    /// Thickness (2D models only).
    pub thickness: f64,
    /// Body force.
    pub force: BodyForce<'a>,
    /// Boundary conditions per marker (unlisted: traction free).
    pub bcs: &'a MarkerBcs<ElasticBoundary>,
}

/// Elasticity matrix in Voigt notation.
#[must_use]
pub fn constitutive_matrix(
    model: Model,
    m: &Material,
) -> Vec<Vec<f64>> {
    let (e, nu) = (m.young, m.poisson);
    match model {
        | Model::PlaneStress => {
            let c = e / (1.0 - nu * nu);
            vec![
                vec![c, c * nu, 0.0],
                vec![c * nu, c, 0.0],
                vec![0.0, 0.0, c * (1.0 - nu) / 2.0],
            ]
        },
        | Model::PlaneStrain => {
            let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            vec![
                vec![c * (1.0 - nu), c * nu, 0.0],
                vec![c * nu, c * (1.0 - nu), 0.0],
                vec![0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0],
            ]
        },
        | Model::Solid3d => {
            let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            let (a, b, g) = (c * (1.0 - nu), c * nu, c * (1.0 - 2.0 * nu) / 2.0);
            let mut d = vec![vec![0.0; 6]; 6];
            for i in 0..3 {
                for j in 0..3 {
                    d[i][j] = if i == j { a } else { b };
                }
                d[3 + i][3 + i] = g;
            }
            d
        },
    }
}

fn strain_column(
    dim: usize,
    comp: usize,
    g: [f64; 3],
) -> Vec<f64> {
    if dim == 2 {
        match comp {
            | 0 => vec![g[0], 0.0, g[1]],
            | _ => vec![0.0, g[1], g[0]],
        }
    } else {
        match comp {
            | 0 => vec![g[0], 0.0, 0.0, g[1], 0.0, g[2]],
            | 1 => vec![0.0, g[1], 0.0, g[0], g[2], 0.0],
            | _ => vec![0.0, 0.0, g[2], 0.0, g[1], g[0]],
        }
    }
}

/// Solution of an elasticity problem.
pub struct ElasticSolution<'a> {
    /// The finite element space.
    pub space: &'a FemSpace,
    /// Displacements, `u[dof * dim + component]`.
    pub u: Vec<f64>,
    /// Constitutive matrix of the homogeneous material (unused when
    /// `cell_c` is present).
    c: Vec<Vec<f64>>,
    /// Constitutive matrix of every cell (heterogeneous material).
    cell_c: Option<Vec<Vec<Vec<f64>>>>,
}

impl ElasticSolution<'_> {
    /// Displacement vector at dof `d`.
    #[must_use]
    pub fn displacement(
        &self,
        d: usize,
    ) -> [f64; 3] {
        let dim = self.space.mesh.dim;
        let mut v = [0.0; 3];
        for a in 0..dim {
            v[a] = self.u[d * dim + a];
        }
        v
    }

    /// Displacement at a point (point location by search).
    #[must_use]
    pub fn eval(
        &self,
        x: [f64; 3],
    ) -> Option<[f64; 3]> {
        let dim = self.space.mesh.dim;
        let mut out = [0.0; 3];
        for a in 0..dim {
            let comp: Vec<f64> = (0..self.space.size())
                .map(|d| self.u[d * dim + a])
                .collect();
            out[a] = self.space.eval(&comp, x)?;
        }
        Some(out)
    }

    /// Stress (Voigt) at the centroid of a cell, using that cell's
    /// constitutive matrix (heterogeneous material) or the homogeneous one.
    #[must_use]
    pub fn stress_at_cell(
        &self,
        cell: usize,
    ) -> Vec<f64> {
        let dim = self.space.mesh.dim;
        let lam = vec![1.0 / (dim + 1) as f64; dim + 1];
        let (_, dl) = self.space.basis(&lam);
        let g = self.space.mesh.bary_gradients(cell);
        let cd = &self.space.cell_dofs[cell];
        let nv = if dim == 2 { 3 } else { 6 };
        let mut eps = vec![0.0; nv];
        for (i, d) in dl.iter().enumerate() {
            let mut gr = [0.0; 3];
            for (m, dm) in d.iter().enumerate() {
                for a in 0..3 {
                    gr[a] += dm * g[m][a];
                }
            }
            for comp in 0..dim {
                let col = strain_column(dim, comp, gr);
                let ui = self.u[cd[i] * dim + comp];
                for r in 0..nv {
                    eps[r] += col[r] * ui;
                }
            }
        }
        let c = match &self.cell_c {
            | Some(cells) => &cells[cell],
            | None => &self.c,
        };
        (0..nv)
            .map(|r| (0..nv).map(|k| c[r][k] * eps[k]).sum())
            .collect()
    }

    /// Von Mises equivalent of a Voigt stress (2D: plane-stress form with `σ_zz = 0`).
    #[must_use]
    pub fn von_mises(s: &[f64]) -> f64 {
        if s.len() == 3 {
            (s[0] * s[0] - s[0] * s[1] + s[1] * s[1] + 3.0 * s[2] * s[2]).sqrt()
        } else {
            let d = [s[0] - s[1], s[1] - s[2], s[2] - s[0]];
            (0.5 * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
                .sqrt()
        }
    }
}

/// Solves the static problem.
///
/// # Errors
/// If the mesh dimension does not match the model, a displacement-free body
/// (singular system) or an unsupported boundary type is met.
pub fn solve<'a>(
    space: &'a FemSpace,
    prob: &ElasticityProblem<'_>,
    solver: LinearSolver,
) -> Result<ElasticSolution<'a>, String> {
    let dim = space.mesh.dim;
    if (prob.model == Model::Solid3d) != (dim == 3) || dim == 1 {
        return Err("model and mesh dimension do not match".into());
    }
    let mat_at = |_: [f64; 3]| prob.material;
    solve_with_material(
        space,
        prob.model,
        prob.material.thickness,
        &mat_at,
        prob.force,
        prob.bcs,
        solver,
        false,
    )
}

/// Solves the static problem with a spatially varying material.
///
/// `E(x)` and `ν(x)` are evaluated at every quadrature point; the per-cell
/// stress recovery ([`ElasticSolution::stress_at_cell`]) uses the
/// constitutive matrix at the cell centroid. A constant material reproduces
/// [`solve`] to round-off.
///
/// # Errors
/// As [`solve`].
pub fn solve_heterogeneous<'a>(
    space: &'a FemSpace,
    prob: &ElasticityProblemVar<'_>,
    solver: LinearSolver,
) -> Result<ElasticSolution<'a>, String> {
    let dim = space.mesh.dim;
    if (prob.model == Model::Solid3d) != (dim == 3) || dim == 1 {
        return Err("model and mesh dimension do not match".into());
    }
    let mat_at = |x: [f64; 3]| {
        Material {
            young: (prob.young)(x),
            poisson: (prob.poisson)(x),
            thickness: prob.thickness,
        }
    };
    solve_with_material(
        space,
        prob.model,
        prob.thickness,
        &mat_at,
        prob.force,
        prob.bcs,
        solver,
        true,
    )
}

fn checked_material(
    model: Model,
    material: Material,
) -> Result<Vec<Vec<f64>>, String> {
    if !material.young.is_finite()
        || material.young <= 0.0
        || !material.poisson.is_finite()
        || material.poisson <= -1.0
        || material.poisson >= 0.5
    {
        return Err("elastic material requires finite E>0 and -1<nu<1/2".into());
    }
    Ok(constitutive_matrix(model, &material))
}

fn stiffness_pattern(space: &FemSpace) -> Csr {
    let dim = space.mesh.dim;
    let rows = space.size() * dim;
    let mut neighbours = vec![Vec::new(); space.size()];
    for nodes in space.cell_dofs.iter().take(space.mesh.cells.len()) {
        for &node in nodes {
            neighbours[node].extend_from_slice(nodes);
        }
    }
    for columns in &mut neighbours {
        columns.sort_unstable();
        columns.dedup();
    }
    let entries = neighbours.iter().map(Vec::len).sum::<usize>() * dim * dim;
    let mut indices = Vec::with_capacity(entries);
    let mut indptr = Vec::with_capacity(rows + 1);
    indptr.push(0);
    for columns in &neighbours {
        for _ in 0..dim {
            for &column in columns {
                indices.extend(column * dim..(column + 1) * dim);
            }
            indptr.push(indices.len());
        }
    }
    // The first addition to -0.0 retains either sign of a zero contribution,
    // just as inserting the first triplet did.
    Csr {
        rows,
        cols: rows,
        indptr,
        indices,
        values: vec![-0.0; entries],
    }
}

fn cell_column_offsets(
    stiffness: &Csr,
    nodes: &[usize],
    dim: usize,
    offsets: &mut Vec<usize>,
) -> Result<(), String> {
    offsets.clear();
    for &node in nodes {
        let row = node * dim;
        let columns = &stiffness.indices[stiffness.indptr[row]..stiffness.indptr[row + 1]];
        for &other in nodes {
            let offset = columns
                .binary_search(&(other * dim))
                .map_err(|_| "elastic element connectivity does not match the stiffness pattern")?;
            offsets.push(offset);
        }
    }
    Ok(())
}

fn replace_displacement_rows(
    stiffness: &mut Csr,
    prescribed: &[Option<f64>],
    load: &mut [f64],
) {
    let mut next = 0;
    for (row, value) in prescribed.iter().enumerate() {
        let (start, end) = (stiffness.indptr[row], stiffness.indptr[row + 1]);
        stiffness.indptr[row] = next;
        if let Some(value) = value {
            // Every prescribed degree of freedom belongs to a cell, so this row
            // has at least its diagonal entry and compaction never grows it.
            stiffness.indices[next] = row;
            stiffness.values[next] = 1.0;
            load[row] = *value;
            next += 1;
        } else {
            if next != start {
                stiffness.indices.copy_within(start..end, next);
                stiffness.values.copy_within(start..end, next);
            }
            next += end - start;
        }
    }
    stiffness.indptr[stiffness.rows] = next;
    stiffness.indices.truncate(next);
    stiffness.values.truncate(next);
}

/// Assembly and solve shared by [`solve`] and [`solve_heterogeneous`]:
/// `mat_at` is sampled at every quadrature point; when `per_cell_c` is set the
/// constitutive matrix at each cell centroid is stored for stress recovery.
#[allow(clippy::too_many_arguments)]
fn solve_with_material<'a>(
    space: &'a FemSpace,
    model: Model,
    thickness: f64,
    mat_at: &dyn Fn([f64; 3]) -> Material,
    force: BodyForce<'_>,
    bcs: &MarkerBcs<ElasticBoundary>,
    solver: LinearSolver,
    per_cell_c: bool,
) -> Result<ElasticSolution<'a>, String> {
    let dim = space.mesh.dim;
    if dim == 2 && (!thickness.is_finite() || thickness <= 0.0) {
        return Err("elastic thickness must be finite positive".into());
    }
    let scale = if dim == 2 { thickness } else { 1.0 };
    let mut stiffness = stiffness_pattern(space);
    let ndof = stiffness.rows;
    let (pts, w) = simplex_rule(dim, space.degree + 1);
    let mut offsets = Vec::new();
    let mut load = vec![0.0; ndof];
    let mut cell_c: Vec<Vec<Vec<f64>>> = if per_cell_c {
        Vec::with_capacity(space.mesh.cells.len())
    } else {
        Vec::new()
    };
    for ci in 0..space.mesh.cells.len() {
        let vol = space.mesh.measure(ci);
        let g = space.mesh.bary_gradients(ci);
        let cd = &space.cell_dofs[ci];
        let nl = cd.len();
        cell_column_offsets(&stiffness, cd, dim, &mut offsets)?;
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
            let cols: Vec<Vec<Vec<f64>>> = grads
                .iter()
                .map(|gr| (0..dim).map(|a| strain_column(dim, a, *gr)).collect())
                .collect();
            let mut x = [0.0; 3];
            for (m, &l) in p.iter().enumerate() {
                for a in 0..3 {
                    x[a] += l * space.mesh.nodes[space.mesh.cells[ci][m]][a];
                }
            }
            let c = checked_material(model, mat_at(x))?;
            let nv = c.len();
            let f = (force)(x);
            if f.iter().any(|v| !v.is_finite()) {
                return Err("elastic body force must be finite".into());
            }
            for i in 0..nl {
                for a in 0..dim {
                    let row = cd[i] * dim + a;
                    let row_start = stiffness.indptr[row];
                    load[row] += scale * vol * wq * f[a] * b[i];
                    // C * B_j,b
                    for j in 0..nl {
                        for bcomp in 0..dim {
                            let mut s = 0.0;
                            for r in 0..nv {
                                let cb: f64 = (0..nv).map(|k| c[r][k] * cols[j][bcomp][k]).sum();
                                s += cols[i][a][r] * cb;
                            }
                            let offset = row_start + offsets[i * nl + j] + bcomp;
                            stiffness.values[offset] += scale * vol * wq * s;
                        }
                    }
                }
            }
        }
        if per_cell_c {
            let mut xc = [0.0; 3];
            for m in 0..=dim {
                for a in 0..3 {
                    xc[a] += space.mesh.nodes[space.mesh.cells[ci][m]][a] / (dim + 1) as f64;
                }
            }
            cell_c.push(checked_material(model, mat_at(xc))?);
        }
    }
    let mut dirichlet: Vec<Option<f64>> = vec![None; ndof];
    let (fp, fw) = simplex_rule(dim - 1, space.degree + 1);
    for &(ci, face, mk) in &space.mesh.facets {
        let Some(bc) = bcs.get(mk) else { continue };
        match bc {
            | ElasticBoundary::Free => {},
            | ElasticBoundary::Displacement(vals) => {
                for d in space.facet_dofs(ci, face) {
                    for a in 0..dim {
                        if let Some(v) = &vals[a] {
                            dirichlet[d * dim + a] = Some(v.eval(0.0, space.dofs[d]));
                        }
                    }
                }
            },
            | ElasticBoundary::Traction(tv) => {
                let meas = space.mesh.facet_measure(ci, face);
                let cd = &space.cell_dofs[ci];
                for (q, wq) in fp.iter().zip(&fw) {
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
                    for i in 0..b.len() {
                        for a in 0..dim {
                            load[cd[i] * dim + a] += scale * meas * wq * tv[a].eval(0.0, x) * b[i];
                        }
                    }
                }
            },
        }
    }
    replace_displacement_rows(&mut stiffness, &dirichlet, &mut load);
    let u = LinearSystem::new(stiffness, solver)?.solve(&load, None)?;
    let c0 = if per_cell_c {
        Vec::new()
    } else {
        checked_material(model, mat_at([0.0; 3]))?
    };
    Ok(ElasticSolution {
        space,
        u,
        c: c0,
        cell_c: per_cell_c.then_some(cell_c),
    })
}

/// Owning adapter snapshot used by session capabilities. This keeps the native
/// constitutive data with its native space without self-referential borrowing.
pub(crate) struct OwnedElasticSolution {
    space: FemSpace,
    u: Vec<f64>,
    c: Vec<Vec<f64>>,
    cell_c: Option<Vec<Vec<Vec<f64>>>>,
}
impl ElasticSolution<'_> {
    pub(crate) fn into_owned(self) -> OwnedElasticSolution {
        OwnedElasticSolution {
            space: self.space.clone(),
            u: self.u,
            c: self.c,
            cell_c: self.cell_c,
        }
    }
}
impl OwnedElasticSolution {
    pub(crate) fn with<R>(
        &self,
        f: impl FnOnce(&ElasticSolution<'_>) -> R,
    ) -> R {
        let borrowed = ElasticSolution {
            space: &self.space,
            u: self.u.clone(),
            c: self.c.clone(),
            cell_c: self.cell_c.clone(),
        };
        f(&borrowed)
    }
}
