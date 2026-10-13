//! # Iterative 3D BEM with hierarchical (ACA) compression
//!
//! Solves the same collocation boundary integral equation as
//! [`solve_laplace_bem_3d_mesh`](crate::sim::physics_bem::solve_laplace_bem_3d_mesh)
//! (interior Laplace problem on a closed triangulated surface, constant
//! elements, Dirichlet or Neumann data per element) without forming the dense
//! matrices:
//!
//! * the collocation points are organised in a binary cluster tree
//!   (recursive bisection along the widest axis);
//! * block pairs satisfying the admissibility condition
//!   `min(diam_s, diam_t) <= η · dist(s, t)` are compressed with **adaptive
//!   cross approximation** (ACA with partial pivoting) of both the
//!   single-layer matrix `G` and the double-layer matrix `H`; the near field
//!   is stored densely;
//! * the system is solved with restarted **GMRES** (Jacobi preconditioned)
//!   using hierarchical matrix-vector products.
//!
//! The rigid-body diagonal of `H` is the exact value `1/2` (a collocation
//! point at the centre of a flat face of a closed surface), so no row sums
//! are needed.
//!
//! Also provided: pure-Neumann interior problems through a bordered system
//! that removes the constant nullspace, **exterior** problems
//! ([`solve_laplace_bem_3d_exterior`], decaying at infinity) and a dense
//! **linear-element** (P1) collocation solver ([`solve_laplace_bem_3d_p1`]).
//!
//! Limitations: compression helps from a few thousand elements upwards; the
//! P1 solver is dense and not compressed; no cluster-tree rebalancing for
//! strongly anisotropic meshes; exterior problems assume a closed surface.
#![allow(
    clippy::manual_memcpy,
    clippy::too_long_first_doc_paragraph,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::missing_const_for_fn,
    clippy::float_cmp,
    clippy::type_complexity,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::manual_midpoint,
    clippy::items_after_statements,
    clippy::single_match_else
)]

use crate::kernels::krylov::gmres;
use crate::sim::physics_bem::{
    Bem3dSolution, BoundaryCondition, SurfaceMesh3D, inv_r_integral, inv_r_integral_in_plane, solid_angle,
};
use rayon::prelude::*;
use std::f64::consts::PI;

type V3 = [f64; 3];

/// Options of the compressed solver.
#[derive(Clone, Copy, Debug)]
pub struct AcaOptions {
    /// Maximum number of collocation points in a leaf cluster.
    pub leaf_size: usize,
    /// Admissibility parameter `η`.
    pub eta: f64,
    /// Relative ACA tolerance of the low-rank blocks.
    pub aca_tol: f64,
    /// GMRES relative tolerance.
    pub gmres_tol: f64,
    /// GMRES restart length.
    pub restart: usize,
    /// Maximum GMRES iterations.
    pub max_iter: usize,
}

impl Default for AcaOptions {
    fn default() -> Self {
        Self { leaf_size: 32, eta: 1.0, aca_tol: 1e-6, gmres_tol: 1e-10, restart: 60, max_iter: 600 }
    }
}

/// Compression and convergence statistics.
#[derive(Clone, Copy, Debug, Default)]
pub struct AcaStats {
    /// Dense (near-field) blocks.
    pub dense_blocks: usize,
    /// Low-rank (far-field) blocks.
    pub lowrank_blocks: usize,
    /// Stored entries of `G` and `H` together.
    pub stored_entries: usize,
    /// Entries of the two dense matrices (`2 n²`).
    pub dense_entries: usize,
    /// GMRES iterations.
    pub iterations: usize,
    /// Final relative residual.
    pub residual: f64,
}

impl AcaStats {
    /// Stored fraction relative to dense storage.
    #[must_use]
    pub fn compression(&self) -> f64 {
        self.stored_entries as f64 / self.dense_entries.max(1) as f64
    }
}

struct Cluster {
    lo: usize,
    hi: usize,
    bmin: V3,
    bmax: V3,
    kids: Option<(usize, usize)>,
}

enum Block {
    Dense { r0: usize, c0: usize, m: usize, n: usize, data: Vec<f64> },
    LowRank { r0: usize, c0: usize, m: usize, n: usize, k: usize, u: Vec<f64>, v: Vec<f64> },
}

impl Block {
    fn stored(&self) -> usize {
        match self {
            Self::Dense { data, .. } => data.len(),
            Self::LowRank { u, v, .. } => u.len() + v.len(),
        }
    }

    /// `y += B x` in permuted ordering.
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        match self {
            Self::Dense { r0, c0, m, n, data } => {
                for i in 0..*m {
                    let mut s = 0.0;
                    for j in 0..*n {
                        s += data[i * n + j] * x[c0 + j];
                    }
                    y[r0 + i] += s;
                }
            },
            Self::LowRank { r0, c0, m, n, k, u, v } => {
                let mut t = vec![0.0; *k];
                for q in 0..*k {
                    t[q] = (0..*n).map(|j| v[q * n + j] * x[c0 + j]).sum();
                }
                for i in 0..*m {
                    y[r0 + i] += (0..*k).map(|q| u[i * k + q] * t[q]).sum::<f64>();
                }
            },
        }
    }
}

struct HMatrix {
    blocks: Vec<Block>,
}

impl HMatrix {
    fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let n = x.len();
        self.blocks
            .par_chunks(64)
            .map(|chunk| {
                let mut y = vec![0.0; n];
                for b in chunk {
                    b.apply(x, &mut y);
                }
                y
            })
            .reduce(
                || vec![0.0; n],
                |mut a, b| {
                    for (p, q) in a.iter_mut().zip(&b) {
                        *p += q;
                    }
                    a
                },
            )
    }
}

fn build_tree(cent: &[V3], perm: &mut [usize], lo: usize, hi: usize, leaf: usize, out: &mut Vec<Cluster>) -> usize {
    let mut bmin = [f64::INFINITY; 3];
    let mut bmax = [f64::NEG_INFINITY; 3];
    for &p in &perm[lo..hi] {
        for a in 0..3 {
            bmin[a] = bmin[a].min(cent[p][a]);
            bmax[a] = bmax[a].max(cent[p][a]);
        }
    }
    let id = out.len();
    out.push(Cluster { lo, hi, bmin, bmax, kids: None });
    if hi - lo > leaf {
        let axis = (0..3).max_by(|&a, &b| (bmax[a] - bmin[a]).total_cmp(&(bmax[b] - bmin[b]))).unwrap_or(0);
        perm[lo..hi].sort_by(|&p, &q| cent[p][axis].total_cmp(&cent[q][axis]));
        let mid = lo + (hi - lo) / 2;
        let l = build_tree(cent, perm, lo, mid, leaf, out);
        let r = build_tree(cent, perm, mid, hi, leaf, out);
        out[id].kids = Some((l, r));
    }
    id
}

fn diam(c: &Cluster) -> f64 {
    (0..3).map(|a| (c.bmax[a] - c.bmin[a]).powi(2)).sum::<f64>().sqrt()
}

fn dist(a: &Cluster, b: &Cluster) -> f64 {
    (0..3)
        .map(|k| {
            let d = (a.bmin[k] - b.bmax[k]).max(b.bmin[k] - a.bmax[k]).max(0.0);
            d * d
        })
        .sum::<f64>()
        .sqrt()
}

fn partition(tree: &[Cluster], s: usize, t: usize, eta: f64, near: &mut Vec<(usize, usize)>, far: &mut Vec<(usize, usize)>) {
    let (cs, ct) = (&tree[s], &tree[t]);
    if s != t && diam(cs).min(diam(ct)) <= eta * dist(cs, ct) {
        far.push((s, t));
        return;
    }
    match (cs.kids, ct.kids) {
        (None, None) => near.push((s, t)),
        (Some((a, b)), None) => {
            partition(tree, a, t, eta, near, far);
            partition(tree, b, t, eta, near, far);
        },
        (None, Some((a, b))) => {
            partition(tree, s, a, eta, near, far);
            partition(tree, s, b, eta, near, far);
        },
        (Some((a, b)), Some((c, d))) => {
            if cs.hi - cs.lo >= ct.hi - ct.lo {
                partition(tree, a, t, eta, near, far);
                partition(tree, b, t, eta, near, far);
            } else {
                partition(tree, s, c, eta, near, far);
                partition(tree, s, d, eta, near, far);
            }
        },
    }
}

/// ACA with partial pivoting of an `m x n` block given by `entry(i, j)`.
fn aca(m: usize, n: usize, tol: f64, entry: &dyn Fn(usize, usize) -> f64) -> (usize, Vec<f64>, Vec<f64>) {
    let kmax = m.min(n);
    let mut us: Vec<Vec<f64>> = Vec::new();
    let mut vs: Vec<Vec<f64>> = Vec::new();
    let mut used_rows = vec![false; m];
    let mut norm2 = 0.0;
    let mut i = 0;
    for _ in 0..kmax {
        used_rows[i] = true;
        let mut row: Vec<f64> = (0..n).map(|j| entry(i, j)).collect();
        for (u, v) in us.iter().zip(&vs) {
            for j in 0..n {
                row[j] -= u[i] * v[j];
            }
        }
        let (jp, piv) = row.iter().enumerate().fold((0, 0.0_f64), |b, (j, &x)| if x.abs() > b.1.abs() { (j, x) } else { b });
        if piv.abs() < 1e-300 {
            // zero row: pick another unused row
            match (0..m).find(|&r| !used_rows[r]) {
                Some(r) => {
                    i = r;
                    continue;
                },
                None => break,
            }
        }
        let v: Vec<f64> = row.iter().map(|x| x / piv).collect();
        let mut col: Vec<f64> = (0..m).map(|r| entry(r, jp)).collect();
        for (u, vv) in us.iter().zip(&vs) {
            for r in 0..m {
                col[r] -= u[r] * vv[jp];
            }
        }
        let nu: f64 = col.iter().map(|x| x * x).sum::<f64>().sqrt();
        let nv: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        // update Frobenius norm of the approximation
        let mut cross = 0.0;
        for (u, vv) in us.iter().zip(&vs) {
            let du: f64 = (0..m).map(|r| u[r] * col[r]).sum();
            let dv: f64 = (0..n).map(|j| vv[j] * v[j]).sum();
            cross += du * dv;
        }
        norm2 += nu * nu * nv * nv + 2.0 * cross;
        // next row pivot: largest entry of the new column among unused rows
        let next = (0..m).filter(|&r| !used_rows[r]).fold(None, |b: Option<(usize, f64)>, r| match b {
            Some((_, x)) if x >= col[r].abs() => b,
            _ => Some((r, col[r].abs())),
        });
        us.push(col);
        vs.push(v);
        if nu * nv <= tol * norm2.max(0.0).sqrt() {
            break;
        }
        match next {
            Some((r, _)) => i = r,
            None => break,
        }
    }
    let k = us.len();
    let mut u = vec![0.0; m * k];
    let mut v = vec![0.0; k * n];
    for q in 0..k {
        for r in 0..m {
            u[r * k + q] = us[q][r];
        }
        v[q * n..(q + 1) * n].copy_from_slice(&vs[q]);
    }
    (k, u, v)
}

fn build_hmatrix(
    tree: &[Cluster],
    near: &[(usize, usize)],
    far: &[(usize, usize)],
    tol: f64,
    entry: &(dyn Fn(usize, usize) -> f64 + Sync),
    perm: &[usize],
) -> HMatrix {
    let mut jobs: Vec<(usize, usize, bool)> = near.iter().map(|&(s, t)| (s, t, false)).collect();
    jobs.extend(far.iter().map(|&(s, t)| (s, t, true)));
    let blocks: Vec<Block> = jobs
        .par_iter()
        .map(|&(s, t, is_far)| {
            let (cs, ct) = (&tree[s], &tree[t]);
            let (m, n) = (cs.hi - cs.lo, ct.hi - ct.lo);
            let e = |i: usize, j: usize| entry(perm[cs.lo + i], perm[ct.lo + j]);
            if is_far {
                let (k, u, v) = aca(m, n, tol, &e);
                if k * (m + n) < m * n {
                    return Block::LowRank { r0: cs.lo, c0: ct.lo, m, n, k, u, v };
                }
            }
            let mut data = vec![0.0; m * n];
            for i in 0..m {
                for j in 0..n {
                    data[i * n + j] = e(i, j);
                }
            }
            Block::Dense { r0: cs.lo, c0: ct.lo, m, n, data }
        })
        .collect();
    HMatrix { blocks }
}

/// Solves the interior Laplace BEM problem with ACA-compressed operators and GMRES.
///
/// A pure Neumann problem (flux prescribed on every element) is singular up
/// to an additive constant; it is solved as the nonsingular augmented system
/// `[H 1; aᵀ 0] [u; λ] = [G q; 0]` with the area weights `a`, which selects
/// the solution of zero area-weighted mean. The fluxes must satisfy the
/// compatibility condition `Σ q_j a_j = 0`.
///
/// # Errors
/// For an invalid mesh, wrong number or non-finite boundary data, an
/// incompatible pure Neumann problem, or GMRES failure to converge.
pub fn solve_laplace_bem_3d_aca(
    mesh: &SurfaceMesh3D,
    bcs: &[BoundaryCondition<f64>],
    opts: &AcaOptions,
) -> Result<(Bem3dSolution, AcaStats), String> {
    solve_impl(mesh, bcs, opts, false)
}

/// Solves the **exterior** Laplace problem (domain outside the closed surface,
/// solution decaying at infinity) with ACA compression and GMRES.
///
/// The boundary data follow the surface orientation: the flux is the derivative
/// along the outward normal of the body, i.e. into the solution domain. The
/// returned `q` uses the same convention (`u = 1/r` on the unit sphere has
/// `q = -1`). Dirichlet, Neumann or mixed data are all nonsingular.
///
/// # Errors
/// As [`solve_laplace_bem_3d_aca`] (except that pure Neumann data are allowed).
pub fn solve_laplace_bem_3d_exterior(
    mesh: &SurfaceMesh3D,
    bcs: &[BoundaryCondition<f64>],
    opts: &AcaOptions,
) -> Result<(Bem3dSolution, AcaStats), String> {
    solve_impl(mesh, bcs, opts, true)
}

/// Potential at a point **outside** the body from an exterior solution:
/// `u(x) = -Σ_j (G_xj q_j + Ω_xj u_j / 4π)` with `q` along the body's outward normal.
///
/// # Errors
/// If the solution does not match the mesh.
pub fn evaluate_potential_3d_exterior(point: [f64; 3], mesh: &SurfaceMesh3D, solution: &Bem3dSolution) -> Result<f64, String> {
    let n = mesh.triangles.len();
    if solution.u.len() != n || solution.q.len() != n {
        return Err("solution does not match the mesh".into());
    }
    let mut result = 0.0;
    for t in 0..n {
        let tri = mesh.corners(t).ok_or("invalid triangle index")?;
        let g = inv_r_integral(point, tri, 0) / (4.0 * PI);
        let h = solid_angle(point, tri) / (4.0 * PI);
        result -= g * solution.q[t] + h * solution.u[t];
    }
    Ok(result)
}

fn flipped(mesh: &SurfaceMesh3D) -> SurfaceMesh3D {
    SurfaceMesh3D { vertices: mesh.vertices.clone(), triangles: mesh.triangles.iter().map(|&[a, b, c]| [a, c, b]).collect() }
}

fn solve_impl(
    mesh_in: &SurfaceMesh3D,
    bcs_in: &[BoundaryCondition<f64>],
    opts: &AcaOptions,
    exterior: bool,
) -> Result<(Bem3dSolution, AcaStats), String> {
    mesh_in.validate()?;
    let n = mesh_in.triangles.len();
    if bcs_in.len() != n {
        return Err(format!("expected {n} boundary conditions, got {}", bcs_in.len()));
    }
    if bcs_in.iter().any(|bc| matches!(bc, BoundaryCondition::Potential(v) | BoundaryCondition::Flux(v) if !v.is_finite())) {
        return Err("boundary data must be finite".into());
    }
    // exterior: interior formulation on the reversed orientation with flipped flux sign
    let flipped_mesh;
    let (mesh, owned_bcs);
    let bcs: &[BoundaryCondition<f64>] = if exterior {
        flipped_mesh = flipped(mesh_in);
        mesh = &flipped_mesh;
        owned_bcs = bcs_in
            .iter()
            .map(|b| match *b {
                BoundaryCondition::Flux(v) => BoundaryCondition::Flux(-v),
                BoundaryCondition::Potential(v) => BoundaryCondition::Potential(v),
            })
            .collect::<Vec<_>>();
        &owned_bcs
    } else {
        mesh = mesh_in;
        bcs_in
    };
    let all_flux = bcs.iter().all(|bc| matches!(bc, BoundaryCondition::Flux(_)));
    let augmented = all_flux && !exterior;
    let areas: Vec<f64> = (0..n).filter_map(|t| mesh.area(t)).collect();
    if augmented {
        let scale: f64 = bcs.iter().zip(&areas).map(|(b, a)| if let BoundaryCondition::Flux(q) = b { (q * a).abs() } else { 0.0 }).sum();
        let net: f64 = bcs.iter().zip(&areas).map(|(b, a)| if let BoundaryCondition::Flux(q) = b { q * a } else { 0.0 }).sum();
        if net.abs() > 1e-6 * scale.max(f64::MIN_POSITIVE) {
            return Err("pure Neumann data violates the compatibility condition (net flux must vanish)".into());
        }
    }
    let corners: Vec<[V3; 3]> = (0..n).filter_map(|t| mesh.corners(t)).collect();
    let cent: Vec<V3> = corners
        .iter()
        .map(|c| [(c[0][0] + c[1][0] + c[2][0]) / 3.0, (c[0][1] + c[1][1] + c[2][1]) / 3.0, (c[0][2] + c[1][2] + c[2][2]) / 3.0])
        .collect();
    let g_entry = |i: usize, j: usize| -> f64 {
        if i == j {
            inv_r_integral_in_plane(cent[i], corners[j]) / (4.0 * PI)
        } else {
            inv_r_integral(cent[i], corners[j], 0) / (4.0 * PI)
        }
    };
    let h_entry = |i: usize, j: usize| -> f64 {
        if i == j { 0.5 } else { -solid_angle(cent[i], corners[j]) / (4.0 * PI) }
    };
    let mut perm: Vec<usize> = (0..n).collect();
    let mut tree = Vec::new();
    build_tree(&cent, &mut perm, 0, n, opts.leaf_size.max(2), &mut tree);
    let (mut near, mut far) = (Vec::new(), Vec::new());
    partition(&tree, 0, 0, opts.eta, &mut near, &mut far);
    let gm = build_hmatrix(&tree, &near, &far, opts.aca_tol, &g_entry, &perm);
    let hm = build_hmatrix(&tree, &near, &far, opts.aca_tol, &h_entry, &perm);
    let mut inv_perm = vec![0; n];
    for (new, &old) in perm.iter().enumerate() {
        inv_perm[old] = new;
    }
    let is_dir: Vec<bool> = perm.iter().map(|&o| matches!(bcs[o], BoundaryCondition::Potential(_))).collect();
    let u_known: Vec<f64> = perm.iter().map(|&o| if let BoundaryCondition::Potential(v) = bcs[o] { v } else { 0.0 }).collect();
    let q_known: Vec<f64> = perm.iter().map(|&o| if let BoundaryCondition::Flux(v) = bcs[o] { v } else { 0.0 }).collect();
    let hu = hm.matvec(&u_known);
    let gq = gm.matvec(&q_known);
    let mut rhs: Vec<f64> = (0..n).map(|i| gq[i] - hu[i]).collect();
    let dim = if augmented { n + 1 } else { n };
    let area_p: Vec<f64> = perm.iter().map(|&o| areas[o]).collect();
    if augmented {
        rhs.push(0.0);
    }
    let apply = |x: &[f64], y: &mut [f64]| {
        let xd: Vec<f64> = (0..n).map(|i| if is_dir[i] { x[i] } else { 0.0 }).collect();
        let xn: Vec<f64> = (0..n).map(|i| if is_dir[i] { 0.0 } else { x[i] }).collect();
        let gx = gm.matvec(&xd);
        let hx = hm.matvec(&xn);
        for i in 0..n {
            y[i] = hx[i] - gx[i];
        }
        if augmented {
            let lam = x[n];
            for i in 0..n {
                y[i] += lam;
            }
            y[n] = (0..n).map(|i| area_p[i] * x[i]).sum();
        }
    };
    let mut diag: Vec<f64> = (0..n)
        .map(|i| {
            let o = perm[i];
            if is_dir[i] { -g_entry(o, o) } else { 0.5 }
        })
        .collect();
    if augmented {
        diag.push(1.0);
    }
    let precond = |r: &[f64], z: &mut [f64]| {
        for i in 0..dim {
            z[i] = r[i] / diag[i];
        }
    };
    let res = gmres(apply, precond, &rhs, None, opts.restart, opts.gmres_tol, opts.max_iter);
    if !res.converged {
        return Err(format!("GMRES did not converge (residual {})", res.residual));
    }
    let mut u = vec![0.0; n];
    let mut q = vec![0.0; n];
    for o in 0..n {
        let p = inv_perm[o];
        match bcs[o] {
            BoundaryCondition::Potential(v) => {
                u[o] = v;
                q[o] = res.x[p];
            },
            BoundaryCondition::Flux(v) => {
                q[o] = v;
                u[o] = res.x[p];
            },
        }
    }
    if exterior {
        for v in &mut q {
            *v = -*v;
        }
    }
    let stored: usize = gm.blocks.iter().chain(&hm.blocks).map(Block::stored).sum();
    let lowrank = gm.blocks.iter().chain(&hm.blocks).filter(|b| matches!(b, Block::LowRank { .. })).count();
    let stats = AcaStats {
        dense_blocks: gm.blocks.len() + hm.blocks.len() - lowrank,
        lowrank_blocks: lowrank,
        stored_entries: stored,
        dense_entries: 2 * n * n,
        iterations: res.iterations,
        residual: res.residual,
    };
    Ok((Bem3dSolution { u, q }, stats))
}

// ---------------------------------------------------------------------------
// Linear (P1) elements
// ---------------------------------------------------------------------------

/// Boundary data at a mesh vertex for the linear-element BEM.
pub type VertexBc = BoundaryCondition<f64>;

/// Solution of the linear-element BEM: nodal potential and normal derivative.
#[derive(Clone, Debug)]
pub struct BemP1Solution {
    /// Potential at the vertices.
    pub u: Vec<f64>,
    /// Derivative along the outward (body) normal at the vertices.
    pub q: Vec<f64>,
}

const RADON_A: f64 = 0.101_286_507_323_456_3;
const RADON_B: f64 = 0.470_142_064_105_115_1;

fn radon7() -> [([f64; 3], f64); 7] {
    // Radon's weights sum to 1/2 on the reference triangle; scale them to sum to one
    let (wa, wb) = (2.0 * 0.062_969_590_272_413_6, 2.0 * 0.066_197_076_394_253_1);
    [
        ([1.0 / 3.0; 3], 0.225),
        ([RADON_A, RADON_A, 1.0 - 2.0 * RADON_A], wa),
        ([RADON_A, 1.0 - 2.0 * RADON_A, RADON_A], wa),
        ([1.0 - 2.0 * RADON_A, RADON_A, RADON_A], wa),
        ([RADON_B, RADON_B, 1.0 - 2.0 * RADON_B], wb),
        ([RADON_B, 1.0 - 2.0 * RADON_B, RADON_B], wb),
        ([1.0 - 2.0 * RADON_B, RADON_B, RADON_B], wb),
    ]
}

fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn nrm3(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Adds `∫_T φ_k (kernel) dS` for k = 0..3 by adaptive subdivision (`single`
/// selects `1/r`, otherwise `ν·(y - x)/r³`).
fn p1_regular(x: V3, tri: [V3; 3], nu: V3, single: bool, depth: usize, out: &mut [f64; 3]) {
    let centre = [
        (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0,
        (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0,
        (tri[0][2] + tri[1][2] + tri[2][2]) / 3.0,
    ];
    let diam = nrm3(sub3(tri[0], tri[1])).max(nrm3(sub3(tri[1], tri[2]))).max(nrm3(sub3(tri[2], tri[0])));
    let d = nrm3(sub3(centre, x));
    if d < 1.5 * diam && depth < 5 {
        let m = |a: V3, b: V3| [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0];
        let (m01, m12, m20) = (m(tri[0], tri[1]), m(tri[1], tri[2]), m(tri[2], tri[0]));
        // sub-triangle integrals carry the parent's shape functions: evaluate via barycentric lifting
        for (sub, lift) in [
            ([tri[0], m01, m20], [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]]),
            ([m01, tri[1], m12], [[0.5, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.5, 0.5]]),
            ([m20, m12, tri[2]], [[0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [0.0, 0.0, 1.0]]),
            ([m12, m20, m01], [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]),
        ] {
            let mut loc = [0.0; 3];
            p1_regular(x, sub, nu, single, depth + 1, &mut loc);
            // loc[k] integrates the sub-triangle's own vertex function k; parent function j = Σ_k lift[k][j] φ_k
            for j in 0..3 {
                for k in 0..3 {
                    out[j] += lift[k][j] * loc[k];
                }
            }
        }
        return;
    }
    let a2 = {
        let (e1, e2) = (sub3(tri[1], tri[0]), sub3(tri[2], tri[0]));
        nrm3([e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]) * 0.5
    };
    for (bary, w) in radon7() {
        let y = [
            bary[0] * tri[0][0] + bary[1] * tri[1][0] + bary[2] * tri[2][0],
            bary[0] * tri[0][1] + bary[1] * tri[1][1] + bary[2] * tri[2][1],
            bary[0] * tri[0][2] + bary[1] * tri[1][2] + bary[2] * tri[2][2],
        ];
        let r = sub3(y, x);
        let rn = nrm3(r);
        let k = if single { 1.0 / rn } else { (nu[0] * r[0] + nu[1] * r[1] + nu[2] * r[2]) / (rn * rn * rn) };
        for j in 0..3 {
            out[j] += a2 * w * bary[j] * k;
        }
    }
}

/// `∫_T φ_j / |y - x| dS` for `x` equal to vertex `k` of `T` (Duffy transformation).
fn p1_singular(tri: [V3; 3], k: usize, out: &mut [f64; 3]) {
    let order = [k, (k + 1) % 3, (k + 2) % 3];
    let (x, b, c) = (tri[order[0]], tri[order[1]], tri[order[2]]);
    let (e1, e2) = (sub3(b, x), sub3(c, x));
    let area2 = nrm3([e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]);
    // 6-point Gauss-Legendre on [0,1]
    const GX: [f64; 6] = [0.033_765_242_898_423_99, 0.169_395_306_766_867_74, 0.380_690_406_958_401_5, 0.619_309_593_041_598_5, 0.830_604_693_233_132_3, 0.966_234_757_101_576];
    const GW: [f64; 6] = [0.085_662_246_189_585_17, 0.180_380_786_524_069_3, 0.233_956_967_286_345_5, 0.233_956_967_286_345_5, 0.180_380_786_524_069_3, 0.085_662_246_189_585_17];
    for (iu, &u) in GX.iter().enumerate() {
        let _ = u;
        for (iw, &w) in GX.iter().enumerate() {
            let dir = [(1.0 - w) * e1[0] + w * e2[0], (1.0 - w) * e1[1] + w * e2[1], (1.0 - w) * e1[2] + w * e2[2]];
            let f = area2 / nrm3(dir);
            let weight = GW[iu] * GW[iw] * f;
            let xi = GX[iu] * (1.0 - w);
            let eta = GX[iu] * w;
            out[order[0]] += weight * (1.0 - xi - eta);
            out[order[1]] += weight * xi;
            out[order[2]] += weight * eta;
        }
    }
}

/// Solves the interior Laplace problem with **continuous piecewise-linear**
/// potential and flux on the triangulated surface, by collocation at the
/// vertices (dense LU).
///
/// `bcs` gives the data per **vertex**: the potential or the derivative along
/// the outward normal. Single-layer integrals use adaptive subdivision away
/// from the collocation point and a Duffy transformation for the elements
/// that contain it; the double-layer kernel vanishes on those (planar)
/// elements, and the free-term/diagonal is fixed by the rigid-body identity
/// `Σ_j H_ij = 0`. Exterior problems: pass `exterior = true` (orientation
/// reversed and the flux sign convention of
/// [`solve_laplace_bem_3d_exterior`]).
///
/// # Errors
/// For an invalid mesh, wrong data length, or a singular system (pure Neumann
/// interior problems are singular).
pub fn solve_laplace_bem_3d_p1(mesh: &SurfaceMesh3D, bcs: &[VertexBc], exterior: bool) -> Result<BemP1Solution, String> {
    mesh.validate()?;
    let nv = mesh.vertices.len();
    if bcs.len() != nv {
        return Err(format!("expected {nv} vertex conditions, got {}", bcs.len()));
    }
    if bcs.iter().any(|bc| matches!(bc, BoundaryCondition::Potential(v) | BoundaryCondition::Flux(v) if !v.is_finite())) {
        return Err("boundary data must be finite".into());
    }
    let tris: Vec<[V3; 3]> = (0..mesh.triangles.len())
        .filter_map(|t| {
            let c = mesh.corners(t)?;
            Some(if exterior { [c[0], c[2], c[1]] } else { c })
        })
        .collect();
    let idx: Vec<[usize; 3]> = mesh.triangles.iter().map(|&[a, b, c]| if exterior { [a, c, b] } else { [a, b, c] }).collect();
    let normals: Vec<V3> = tris
        .iter()
        .map(|t| {
            let (e1, e2) = (sub3(t[1], t[0]), sub3(t[2], t[0]));
            let c = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
            let l = nrm3(c);
            [c[0] / l, c[1] / l, c[2] / l]
        })
        .collect();
    let rows: Vec<(Vec<f64>, Vec<f64>)> = (0..nv)
        .into_par_iter()
        .map(|i| {
            let x = mesh.vertices[i];
            let mut g = vec![0.0; nv];
            let mut h = vec![0.0; nv];
            for (t, tri) in tris.iter().enumerate() {
                let own = idx[t].iter().position(|&v| v == i);
                let mut gl = [0.0; 3];
                match own {
                    Some(k) => p1_singular(*tri, k, &mut gl),
                    None => {
                        p1_regular(x, *tri, normals[t], true, 0, &mut gl);
                        let mut hl = [0.0; 3];
                        p1_regular(x, *tri, normals[t], false, 0, &mut hl);
                        for j in 0..3 {
                            h[idx[t][j]] -= hl[j] / (4.0 * PI);
                        }
                    },
                }
                for j in 0..3 {
                    g[idx[t][j]] += gl[j] / (4.0 * PI);
                }
            }
            let off: f64 = (0..nv).filter(|&j| j != i).map(|j| h[j]).sum();
            h[i] = -off;
            (g, h)
        })
        .collect();
    let sign = if exterior { -1.0 } else { 1.0 };
    let mut a = crate::kernels::dense::Mat::zeros(nv, nv);
    let mut rhs = vec![0.0; nv];
    for (i, (g, h)) in rows.iter().enumerate() {
        for (j, bc) in bcs.iter().enumerate() {
            match *bc {
                BoundaryCondition::Potential(u) => {
                    a.set(i, j, -g[j]);
                    rhs[i] -= h[j] * u;
                },
                BoundaryCondition::Flux(q) => {
                    a.set(i, j, h[j]);
                    rhs[i] += g[j] * sign * q;
                },
            }
        }
    }
    let lu = crate::kernels::dense::lu_factor(&a).map_err(|e| format!("BEM system is singular: {e:?}"))?;
    let sol = lu.solve(&rhs);
    if sol.iter().any(|v| !v.is_finite()) {
        return Err("BEM system has no finite solution".into());
    }
    let mut u = vec![0.0; nv];
    let mut q = vec![0.0; nv];
    for (j, bc) in bcs.iter().enumerate() {
        match *bc {
            BoundaryCondition::Potential(v) => {
                u[j] = v;
                q[j] = sign * sol[j];
            },
            BoundaryCondition::Flux(v) => {
                q[j] = v;
                u[j] = sol[j];
            },
        }
    }
    Ok(BemP1Solution { u, q })
}
