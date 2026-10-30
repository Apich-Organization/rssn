//! # Large-deformation elasticity, elastodynamics and penalty contact
//!
//! Total-Lagrangian finite elements on the linear simplex meshes of
//! [`fem_simplex`](crate::sim::fem_simplex) (triangles in **plane strain**, tetrahedra):
//!
//! * hyperelastic models ([`HyperModel`]): **St. Venant–Kirchhoff**
//!   `S = λ tr(E) I + 2 μ E` and compressible **neo-Hookean**
//!   `W = μ/2 (I1 - 3) - μ ln J + λ/2 (ln J)²`, with the exact first
//!   Piola–Kirchhoff stress and the exact material tangent
//!   ([`first_piola`], [`tangent_moduli`]), so Newton converges quadratically;
//! * static solution by **load stepping with Newton iterations** and a
//!   backtracking line search ([`solve_static`]): displacement boundary
//!   conditions per marker (ramped), nominal (dead) tractions, body forces;
//! * **elastodynamics** with consistent mass: the Newmark family and the
//!   **generalised-α** method of Chung and Hulbert with user-chosen
//!   high-frequency dissipation `ρ∞` ([`TimeScheme`], [`solve_dynamic`]);
//! * **penalty contact** with a rigid obstacle given by a gap function and its
//!   normal ([`PenaltyContact`]): node-to-rigid-surface, frictionless, with
//!   area-weighted nodal penalty and the consistent linearisation.
//!
//! Limitations: linear (P1) elements only (volumetric locking in nearly
//! incompressible bending), rigid obstacles only, no friction, no
//! self-contact, no damping other than the algorithmic one of generalised-α,
//! and displacement boundary conditions in dynamics are applied as prescribed
//! values of the new displacement.

use crate::kernels::krylov::Csr;
use crate::sim::boundary::ElasticBoundary;
use crate::sim::boundary::MarkerBcs;
use crate::sim::fem_simplex::SimplexMesh;
use crate::sim::structured::LinearSolver;
use crate::sim::structured::LinearSystem;

type M3 = [[f64; 3]; 3];
type Triplets = Vec<(usize, usize, f64)>;

/// Hyperelastic constitutive model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HyperModel {
    /// St. Venant–Kirchhoff (large displacement, small strain energy).
    StVenantKirchhoff,
    /// Compressible neo-Hookean.
    NeoHookean,
}

/// Isotropic hyperelastic material.
#[derive(Clone, Copy, Debug)]
pub struct HyperMaterial {
    /// Model.
    pub model: HyperModel,
    /// Young's modulus (small-strain).
    pub young: f64,
    /// Poisson ratio (small-strain).
    pub poisson: f64,
    /// Mass density (dynamics).
    pub density: f64,
}

impl HyperMaterial {
    /// Lamé parameters `(λ, μ)`.
    #[must_use]
    pub fn lame(&self) -> (f64, f64) {
        let (e, nu) = (self.young, self.poisson);
        (
            e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
            e / (2.0 * (1.0 + nu)),
        )
    }
}

fn det3(m: &M3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn inv3(m: &M3) -> M3 {
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
    inv
}

fn matmul(
    a: &M3,
    b: &M3,
) -> M3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

fn transpose(a: &M3) -> M3 {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

/// First Piola–Kirchhoff stress `P(F)`.
///
/// # Errors
/// If `det F <= 0` (inverted element).
pub fn first_piola(
    mat: &HyperMaterial,
    f: &M3,
) -> Result<M3, String> {
    let (lam, mu) = mat.lame();
    let j = det3(f);
    if j <= 1e-12 {
        return Err("element inverted during the iteration (det F <= 0)".into());
    }
    let mut p = [[0.0; 3]; 3];
    match mat.model {
        | HyperModel::StVenantKirchhoff => {
            let c = matmul(&transpose(f), f);
            let mut e = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    e[a][b] = 0.5 * (c[a][b] - f64::from(u8::from(a == b)));
                }
            }
            let tr = e[0][0] + e[1][1] + e[2][2];
            let mut s = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] = 2.0 * mu * e[a][b] + if a == b { lam * tr } else { 0.0 };
                }
            }
            p = matmul(f, &s);
        },
        | HyperModel::NeoHookean => {
            let fit = transpose(&inv3(f));
            let lnj = j.ln();
            for a in 0..3 {
                for b in 0..3 {
                    p[a][b] = mu * (f[a][b] - fit[a][b]) + lam * lnj * fit[a][b];
                }
            }
        },
    }
    Ok(p)
}

/// Material tangent `A_{iJkL} = ∂P_{iJ}/∂F_{kL}` stored as `a[i][J][k][L]`.
///
/// # Errors
/// If `det F <= 0`.
pub fn tangent_moduli(
    mat: &HyperMaterial,
    f: &M3,
) -> Result<[[M3; 3]; 3], String> {
    let (lam, mu) = mat.lame();
    let j = det3(f);
    if j <= 1e-12 {
        return Err("element inverted during the iteration (det F <= 0)".into());
    }
    let mut a = [[[[0.0; 3]; 3]; 3]; 3];
    let delta = |x: usize, y: usize| f64::from(u8::from(x == y));
    match mat.model {
        | HyperModel::StVenantKirchhoff => {
            let c = matmul(&transpose(f), f);
            let mut e = [[0.0; 3]; 3];
            for p in 0..3 {
                for q in 0..3 {
                    e[p][q] = 0.5 * (c[p][q] - delta(p, q));
                }
            }
            let tr = e[0][0] + e[1][1] + e[2][2];
            let mut s = [[0.0; 3]; 3];
            for p in 0..3 {
                for q in 0..3 {
                    s[p][q] = 2.0 * mu * e[p][q] + if p == q { lam * tr } else { 0.0 };
                }
            }
            let fft = matmul(f, &transpose(f));
            for i in 0..3 {
                for jj in 0..3 {
                    for k in 0..3 {
                        for l in 0..3 {
                            a[i][jj][k][l] = delta(i, k) * s[jj][l]
                                + lam * f[i][jj] * f[k][l]
                                + mu * (fft[i][k] * delta(jj, l) + f[i][l] * f[k][jj]);
                        }
                    }
                }
            }
        },
        | HyperModel::NeoHookean => {
            let fi = inv3(f);
            let lnj = j.ln();
            for i in 0..3 {
                for jj in 0..3 {
                    for k in 0..3 {
                        for l in 0..3 {
                            a[i][jj][k][l] = mu * delta(i, k) * delta(jj, l)
                                + (mu - lam * lnj) * fi[jj][k] * fi[l][i]
                                + lam * fi[jj][i] * fi[l][k];
                        }
                    }
                }
            }
        },
    }
    Ok(a)
}

/// Frictionless penalty contact of the nodes of a boundary marker with a rigid obstacle.
pub struct PenaltyContact {
    /// Boundary marker whose nodes may touch the obstacle.
    pub marker: i32,
    /// Signed gap `g(x)` (negative = penetration) at a current position.
    pub gap: Box<dyn Fn([f64; 3]) -> f64 + Send + Sync>,
    /// Unit normal of the obstacle surface pointing into the admissible region.
    pub normal: Box<dyn Fn([f64; 3]) -> [f64; 3] + Send + Sync>,
    /// Penalty stiffness per unit area per unit penetration (pressure / length).
    pub penalty: f64,
}

/// Static or dynamic problem on a simplex mesh.
pub struct HyperProblem<'a> {
    /// Mesh of triangles (plane strain) or tetrahedra.
    pub mesh: &'a SimplexMesh,
    /// Material.
    pub material: HyperMaterial,
    /// Boundary conditions per marker: displacement (ramped with the load factor),
    /// nominal traction (dead load per reference area), free.
    pub bcs: &'a MarkerBcs<ElasticBoundary>,
    /// Body force per unit reference volume.
    pub force: &'a dyn Fn([f64; 3]) -> [f64; 3],
    /// Contact surfaces.
    pub contacts: Vec<PenaltyContact>,
}

struct Assembler<'a> {
    prob: &'a HyperProblem<'a>,
    dim: usize,
    grads: Vec<Vec<[f64; 3]>>,
    volumes: Vec<f64>,
}

impl<'a> Assembler<'a> {
    fn new(prob: &'a HyperProblem<'a>) -> Result<Self, String> {
        let mesh = prob.mesh;
        if !(2..=3).contains(&mesh.dim) {
            return Err("the mesh must be 2D (plane strain) or 3D".into());
        }
        let grads = (0..mesh.cells.len())
            .map(|c| mesh.bary_gradients(c))
            .collect();
        let volumes = (0..mesh.cells.len()).map(|c| mesh.measure(c)).collect();
        Ok(Self {
            prob,
            dim: mesh.dim,
            grads,
            volumes,
        })
    }

    const fn ndof(&self) -> usize {
        self.prob.mesh.nodes.len() * self.dim
    }

    fn deformation_gradient(
        &self,
        cell: usize,
        u: &[f64],
    ) -> M3 {
        let mut f = [[0.0; 3]; 3];
        for d in 0..3 {
            f[d][d] = 1.0;
        }
        for (a, &node) in self.prob.mesh.cells[cell].iter().enumerate() {
            for i in 0..self.dim {
                for j in 0..self.dim {
                    f[i][j] += u[node * self.dim + i] * self.grads[cell][a][j];
                }
            }
        }
        f
    }

    /// Internal force vector and (optionally) the tangent triplets.
    fn internal(
        &self,
        u: &[f64],
        want_tangent: bool,
    ) -> Result<(Vec<f64>, Triplets), String> {
        let (dim, mat) = (self.dim, &self.prob.material);
        let mut fint = vec![0.0; self.ndof()];
        let mut trip = Vec::new();
        for (c, cell) in self.prob.mesh.cells.iter().enumerate() {
            let f = self.deformation_gradient(c, u);
            let p = first_piola(mat, &f)?;
            let v = self.volumes[c];
            for (a, &na) in cell.iter().enumerate() {
                for i in 0..dim {
                    let s: f64 = (0..dim).map(|j| p[i][j] * self.grads[c][a][j]).sum();
                    fint[na * dim + i] += v * s;
                }
            }
            if want_tangent {
                let at = tangent_moduli(mat, &f)?;
                for (a, &na) in cell.iter().enumerate() {
                    for (b, &nb) in cell.iter().enumerate() {
                        for i in 0..dim {
                            for k in 0..dim {
                                let mut s = 0.0;
                                for j in 0..dim {
                                    for l in 0..dim {
                                        s += at[i][j][k][l]
                                            * self.grads[c][a][j]
                                            * self.grads[c][b][l];
                                    }
                                }
                                trip.push((na * dim + i, nb * dim + k, v * s));
                            }
                        }
                    }
                }
            }
        }
        Ok((fint, trip))
    }

    /// Consistent mass matrix triplets.
    fn mass(&self) -> Vec<(usize, usize, f64)> {
        let dim = self.dim;
        let mut trip = Vec::new();
        let n1 = (dim + 1) as f64;
        for (c, cell) in self.prob.mesh.cells.iter().enumerate() {
            let m = self.prob.material.density * self.volumes[c] / (n1 * (n1 + 1.0));
            for (a, &na) in cell.iter().enumerate() {
                for (b, &nb) in cell.iter().enumerate() {
                    let w = if a == b { 2.0 * m } else { m };
                    for i in 0..dim {
                        trip.push((na * dim + i, nb * dim + i, w));
                    }
                }
            }
        }
        trip
    }

    /// External load (body force and nominal tractions) scaled by `factor`.
    fn external(
        &self,
        factor: f64,
        t: f64,
    ) -> Vec<f64> {
        let (mesh, dim) = (self.prob.mesh, self.dim);
        let mut load = vec![0.0; self.ndof()];
        let n1 = (dim + 1) as f64;
        for (c, cell) in mesh.cells.iter().enumerate() {
            let mut x = [0.0; 3];
            for &n in cell {
                for a in 0..3 {
                    x[a] += mesh.nodes[n][a] / n1;
                }
            }
            let f = (self.prob.force)(x);
            for &n in cell {
                for i in 0..dim {
                    load[n * dim + i] += factor * self.volumes[c] / n1 * f[i];
                }
            }
        }
        for &(ci, face, mk) in &mesh.facets {
            let Some(ElasticBoundary::Traction(tv)) = self.prob.bcs.get(mk) else {
                continue;
            };
            let meas = mesh.facet_measure(ci, face);
            let cell = &mesh.cells[ci];
            let mut x = [0.0; 3];
            let mut cnt = 0.0;
            for (j, &n) in cell.iter().enumerate() {
                if j != face {
                    cnt += 1.0;
                    for a in 0..3 {
                        x[a] += mesh.nodes[n][a];
                    }
                }
            }
            for xa in &mut x {
                *xa /= cnt;
            }
            for (j, &n) in cell.iter().enumerate() {
                if j != face {
                    for i in 0..dim {
                        load[n * dim + i] += factor * meas / cnt * tv[i].eval(t, x);
                    }
                }
            }
        }
        load
    }

    /// Nodal tributary areas of the facets of a marker.
    fn contact_areas(
        &self,
        marker: i32,
    ) -> Vec<(usize, f64)> {
        let mesh = self.prob.mesh;
        let mut area = vec![0.0; mesh.nodes.len()];
        for &(ci, face, mk) in &mesh.facets {
            if mk != marker {
                continue;
            }
            let meas = mesh.facet_measure(ci, face);
            for (j, &n) in mesh.cells[ci].iter().enumerate() {
                if j != face {
                    area[n] += meas / self.dim as f64;
                }
            }
        }
        area.into_iter()
            .enumerate()
            .filter(|(_, a)| *a > 0.0)
            .collect()
    }

    /// Contact force (acting on the body) and tangent for the current displacement.
    fn contact(
        &self,
        u: &[f64],
        force: &mut [f64],
        trip: Option<&mut Vec<(usize, usize, f64)>>,
    ) {
        let dim = self.dim;
        let mut trip = trip;
        for ct in &self.prob.contacts {
            for (n, a) in self.contact_areas(ct.marker) {
                let mut x = self.prob.mesh.nodes[n];
                for i in 0..dim {
                    x[i] += u[n * dim + i];
                }
                let g = (ct.gap)(x);
                if g >= 0.0 {
                    continue;
                }
                let nrm = (ct.normal)(x);
                let k = ct.penalty * a;
                for i in 0..dim {
                    force[n * dim + i] += -k * g * nrm[i];
                }
                if let Some(t) = trip.as_deref_mut() {
                    for i in 0..dim {
                        for j in 0..dim {
                            t.push((n * dim + i, n * dim + j, k * nrm[i] * nrm[j]));
                        }
                    }
                }
            }
        }
    }

    /// Prescribed displacement values (`None` = free) at the load factor / time.
    fn prescribed(
        &self,
        factor: f64,
        t: f64,
    ) -> Vec<Option<f64>> {
        let (mesh, dim) = (self.prob.mesh, self.dim);
        let mut out = vec![None; self.ndof()];
        for &(ci, face, mk) in &mesh.facets {
            let Some(ElasticBoundary::Displacement(vals)) = self.prob.bcs.get(mk) else {
                continue;
            };
            for (j, &n) in mesh.cells[ci].iter().enumerate() {
                if j == face {
                    continue;
                }
                for a in 0..dim {
                    if let Some(v) = &vals[a] {
                        out[n * dim + a] = Some(factor * v.eval(t, mesh.nodes[n]));
                    }
                }
            }
        }
        out
    }
}

fn solve_system(
    trip: &[(usize, usize, f64)],
    fixed: &[Option<f64>],
    rhs: &[f64],
    u: &[f64],
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    let n = rhs.len();
    let mut t2: Vec<(usize, usize, f64)> = trip
        .iter()
        .filter(|(r, _, _)| fixed[*r].is_none())
        .copied()
        .collect();
    let mut b = rhs.to_vec();
    for r in 0..n {
        if let Some(v) = fixed[r] {
            t2.push((r, r, 1.0));
            b[r] = v - u[r];
        }
    }
    LinearSystem::new(Csr::from_triplets(n, n, &t2), solver)?.solve(&b, None)
}

/// Result of a static computation.
#[derive(Clone, Debug)]
pub struct HyperSolution {
    /// Nodal displacements `u[node * dim + component]`.
    pub u: Vec<f64>,
    /// Newton residual norms of every iteration of the last load step.
    pub residuals: Vec<f64>,
    /// Newton iterations per load step.
    pub iterations: Vec<usize>,
}

/// Static solution by load stepping and Newton iterations.
///
/// The load is applied in `steps` equal increments (the displacement data, tractions and body force
/// are scaled by `s / steps`) with Newton iterations (tolerance on the free-dof residual relative to the
/// largest external or internal force norm).
///
/// # Errors
/// If an element inverts, Newton fails to converge or the system is singular.
pub fn solve_static(
    prob: &HyperProblem<'_>,
    steps: usize,
    tol: f64,
    max_iter: usize,
    solver: LinearSolver,
) -> Result<HyperSolution, String> {
    let asm = Assembler::new(prob)?;
    let n = asm.ndof();
    let mut u = vec![0.0; n];
    let mut iterations = Vec::new();
    let mut residuals = Vec::new();
    for s in 1..=steps.max(1) {
        let factor = s as f64 / steps.max(1) as f64;
        let fixed = asm.prescribed(factor, 0.0);
        for (d, v) in fixed.iter().enumerate() {
            if let Some(v) = v {
                u[d] = *v;
            }
        }
        let ext = asm.external(factor, 0.0);
        residuals.clear();
        let mut converged = false;
        let mut count = 0;
        for _ in 0..max_iter {
            let (fint, mut trip) = asm.internal(&u, true)?;
            let mut fcon = vec![0.0; n];
            asm.contact(&u, &mut fcon, Some(&mut trip));
            let res: Vec<f64> = (0..n)
                .map(|d| {
                    if fixed[d].is_some() {
                        0.0
                    } else {
                        ext[d] + fcon[d] - fint[d]
                    }
                })
                .collect();
            let rn = res.iter().map(|x| x * x).sum::<f64>().sqrt();
            let scale = fint
                .iter()
                .chain(&ext)
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                .max(1e-300);
            residuals.push(rn);
            if rn <= tol.max(1e-10) * scale || rn < 1e-14 {
                converged = true;
                break;
            }
            count += 1;
            let du = solve_system(&trip, &fixed, &res, &u, solver)?;
            // backtracking line search on the residual norm
            let mut alpha = 1.0;
            let mut accepted = false;
            for _ in 0..8 {
                let trial: Vec<f64> = u.iter().zip(&du).map(|(a, b)| a + alpha * b).collect();
                if let Ok((fi, _)) = asm.internal(&trial, false) {
                    let mut fc = vec![0.0; n];
                    asm.contact(&trial, &mut fc, None);
                    let rt: f64 = (0..n)
                        .map(|d| {
                            if fixed[d].is_some() {
                                0.0
                            } else {
                                (ext[d] + fc[d] - fi[d]).powi(2)
                            }
                        })
                        .sum::<f64>()
                        .sqrt();
                    if rt < rn * (1.0 - 1e-4 * alpha) || rt <= tol.max(1e-10) * scale {
                        u = trial;
                        accepted = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if !accepted {
                if rn <= 1e-8 * scale {
                    // stagnation at the round-off floor of an already converged state
                    converged = true;
                    break;
                }
                return Err("Newton line search failed".into());
            }
        }
        if !converged {
            return Err(format!("Newton did not converge in load step {s}"));
        }
        iterations.push(count);
    }
    Ok(HyperSolution {
        u,
        residuals,
        iterations,
    })
}

/// Time integration scheme.
#[derive(Clone, Copy, Debug)]
pub enum TimeScheme {
    /// Newmark with parameters `(β, γ)`; `(1/4, 1/2)` is the energy-conserving average acceleration.
    Newmark {
        /// `β`.
        beta: f64,
        /// `γ`.
        gamma: f64,
    },
    /// Generalised-α (Chung–Hulbert) with spectral radius `ρ∞ ∈ [0, 1]` at infinite frequency.
    GeneralizedAlpha {
        /// `ρ∞` (1 = no dissipation, 0 = asymptotic annihilation).
        rho_inf: f64,
    },
}

impl TimeScheme {
    /// `(α_m, α_f, β, γ)`.
    fn parameters(self) -> (f64, f64, f64, f64) {
        match self {
            | Self::Newmark { beta, gamma } => (0.0, 0.0, beta, gamma),
            | Self::GeneralizedAlpha { rho_inf } => {
                let am = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
                let af = rho_inf / (rho_inf + 1.0);
                (am, af, 0.25 * (1.0 - am + af).powi(2), 0.5 - am + af)
            },
        }
    }
}

/// State of a dynamic run.
#[derive(Clone, Debug)]
pub struct DynamicState {
    /// Displacements.
    pub u: Vec<f64>,
    /// Velocities.
    pub v: Vec<f64>,
    /// Accelerations.
    pub a: Vec<f64>,
    /// Time.
    pub t: f64,
}

/// Kinetic energy `½ vᵀ M v` with the consistent mass matrix.
#[must_use]
pub fn kinetic_energy(
    prob: &HyperProblem<'_>,
    v: &[f64],
) -> f64 {
    let Ok(asm) = Assembler::new(prob) else {
        return 0.0;
    };
    let n = asm.ndof();
    let m = Csr::from_triplets(n, n, &asm.mass());
    let mut mv = vec![0.0; n];
    m.matvec(v, &mut mv);
    0.5 * v.iter().zip(&mv).map(|(a, b)| a * b).sum::<f64>()
}

/// Strain energy `∫ W dV` of the displacement field `u`.
///
/// # Errors
/// If an element is inverted.
pub fn strain_energy(
    prob: &HyperProblem<'_>,
    u: &[f64],
) -> Result<f64, String> {
    let asm = Assembler::new(prob)?;
    let (lam, mu) = prob.material.lame();
    let mut w = 0.0;
    for c in 0..prob.mesh.cells.len() {
        let f = asm.deformation_gradient(c, u);
        let j = det3(&f);
        if j <= 1e-12 {
            return Err("element inverted".into());
        }
        let ct = matmul(&transpose(&f), &f);
        let density = match prob.material.model {
            | HyperModel::StVenantKirchhoff => {
                let mut e2 = 0.0;
                let mut tr = 0.0;
                for a in 0..3 {
                    tr += 0.5 * (ct[a][a] - 1.0);
                    for b in 0..3 {
                        let e = 0.5 * (ct[a][b] - f64::from(u8::from(a == b)));
                        e2 += e * e;
                    }
                }
                0.5 * lam * tr * tr + mu * e2
            },
            | HyperModel::NeoHookean => {
                let i1 = ct[0][0] + ct[1][1] + ct[2][2];
                0.5 * mu * (i1 - 3.0) - mu * j.ln() + 0.5 * lam * j.ln().powi(2)
            },
        };
        w += asm.volumes[c] * density;
    }
    Ok(w)
}

/// Integrates `M a + f_int(u) = f_ext(t)` from the given initial state with fixed step `dt`.
///
/// The external load is the body force and the nominal tractions of `prob`
/// (constant in time, applied suddenly at `t = 0`); displacement conditions are
/// evaluated at the new time. `observer` is called after every step with the new state.
///
/// # Errors
/// If an element inverts, Newton fails or the effective matrix is singular.
#[allow(clippy::too_many_arguments)] // the scheme, state, step and solver are independent controls
pub fn solve_dynamic(
    prob: &HyperProblem<'_>,
    scheme: TimeScheme,
    initial: &DynamicState,
    dt: f64,
    steps: usize,
    tol: f64,
    solver: LinearSolver,
    observer: &mut dyn FnMut(&DynamicState),
) -> Result<DynamicState, String> {
    let asm = Assembler::new(prob)?;
    let n = asm.ndof();
    if initial.u.len() != n || initial.v.len() != n || initial.a.len() != n {
        return Err("initial state has the wrong size".into());
    }
    let (am, af, beta, gamma) = scheme.parameters();
    if beta <= 0.0 || dt <= 0.0 {
        return Err("beta and the time step must be positive".into());
    }
    let mass_trip = asm.mass();
    let mass = Csr::from_triplets(n, n, &mass_trip);
    let ext = asm.external(1.0, 0.0);
    let mut st = initial.clone();
    for _ in 0..steps {
        let t_new = st.t + dt;
        let fixed = asm.prescribed(1.0, t_new);
        let utilde: Vec<f64> = (0..n)
            .map(|i| st.u[i] + dt * st.v[i] + dt * dt * (0.5 - beta) * st.a[i])
            .collect();
        let mut u = st.u.clone();
        for (d, v) in fixed.iter().enumerate() {
            if let Some(v) = v {
                u[d] = *v;
            }
        }
        let mut converged = false;
        for _ in 0..30 {
            let ua: Vec<f64> = (0..n).map(|i| (1.0 - af) * u[i] + af * st.u[i]).collect();
            let (fint, ktrip) = asm.internal(&ua, true)?;
            let mut fcon = vec![0.0; n];
            let mut ctrip = Vec::new();
            asm.contact(&ua, &mut fcon, Some(&mut ctrip));
            let a_new: Vec<f64> = (0..n)
                .map(|i| (u[i] - utilde[i]) / (beta * dt * dt))
                .collect();
            let a_mix: Vec<f64> = (0..n)
                .map(|i| (1.0 - am) * a_new[i] + am * st.a[i])
                .collect();
            let mut ma = vec![0.0; n];
            mass.matvec(&a_mix, &mut ma);
            let res: Vec<f64> = (0..n)
                .map(|i| {
                    if fixed[i].is_some() {
                        0.0
                    } else {
                        ext[i] + fcon[i] - fint[i] - ma[i]
                    }
                })
                .collect();
            let rn = res.iter().map(|x| x * x).sum::<f64>().sqrt();
            let scale = fint
                .iter()
                .chain(&ext)
                .chain(&ma)
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                .max(1e-300);
            if rn <= tol.max(1e-10) * scale || rn < 1e-14 {
                converged = true;
                break;
            }
            let mut trip: Vec<(usize, usize, f64)> = ktrip
                .iter()
                .chain(&ctrip)
                .map(|&(r, c, v)| (r, c, (1.0 - af) * v))
                .collect();
            trip.extend(
                mass_trip
                    .iter()
                    .map(|&(r, c, v)| (r, c, v * (1.0 - am) / (beta * dt * dt))),
            );
            let du = solve_system(&trip, &fixed, &res, &u, solver)?;
            for i in 0..n {
                u[i] += du[i];
            }
        }
        if !converged {
            return Err("Newton did not converge in a time step".into());
        }
        let a_new: Vec<f64> = (0..n)
            .map(|i| (u[i] - utilde[i]) / (beta * dt * dt))
            .collect();
        let v_new: Vec<f64> = (0..n)
            .map(|i| st.v[i] + dt * ((1.0 - gamma) * st.a[i] + gamma * a_new[i]))
            .collect();
        st = DynamicState {
            u,
            v: v_new,
            a: a_new,
            t: t_new,
        };
        observer(&st);
    }
    Ok(st)
}

/// Initial state at rest with the consistent initial acceleration `a = M⁻¹ (f_ext - f_int(u))`.
///
/// # Errors
/// If the mass matrix is singular or an element is inverted.
pub fn rest_state(
    prob: &HyperProblem<'_>,
    u0: &[f64],
    solver: LinearSolver,
) -> Result<DynamicState, String> {
    let asm = Assembler::new(prob)?;
    let n = asm.ndof();
    let (fint, _) = asm.internal(u0, false)?;
    let ext = asm.external(1.0, 0.0);
    let fixed = asm.prescribed(1.0, 0.0);
    let rhs: Vec<f64> = (0..n)
        .map(|i| {
            if fixed[i].is_some() {
                0.0
            } else {
                ext[i] - fint[i]
            }
        })
        .collect();
    let uref: Vec<f64> = fixed.iter().map(|v| v.unwrap_or(0.0)).collect();
    let a = solve_system(&asm.mass(), &fixed, &rhs, &uref, solver)?;
    Ok(DynamicState {
        u: u0.to_vec(),
        v: vec![0.0; n],
        a,
        t: 0.0,
    })
}
