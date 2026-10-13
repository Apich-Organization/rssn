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
//! Limitations: interior problems with at least one Dirichlet element (a
//! pure Neumann problem is singular); compression helps from a few thousand
//! elements upwards; no cluster-tree rebalancing for strongly anisotropic
//! meshes.
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
    clippy::option_if_let_else
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
/// # Errors
/// For an invalid mesh, wrong number or non-finite boundary data, a pure
/// Neumann problem, or GMRES failure to converge.
pub fn solve_laplace_bem_3d_aca(
    mesh: &SurfaceMesh3D,
    bcs: &[BoundaryCondition<f64>],
    opts: &AcaOptions,
) -> Result<(Bem3dSolution, AcaStats), String> {
    mesh.validate()?;
    let n = mesh.triangles.len();
    if bcs.len() != n {
        return Err(format!("expected {n} boundary conditions, got {}", bcs.len()));
    }
    if bcs.iter().any(|bc| matches!(bc, BoundaryCondition::Potential(v) | BoundaryCondition::Flux(v) if !v.is_finite())) {
        return Err("boundary data must be finite".into());
    }
    if bcs.iter().all(|bc| matches!(bc, BoundaryCondition::Flux(_))) {
        return Err("pure Neumann problems are singular; prescribe the potential on at least one element".into());
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
    // masks in permuted ordering
    let is_dir: Vec<bool> = perm.iter().map(|&o| matches!(bcs[o], BoundaryCondition::Potential(_))).collect();
    let u_known: Vec<f64> = perm.iter().map(|&o| if let BoundaryCondition::Potential(v) = bcs[o] { v } else { 0.0 }).collect();
    let q_known: Vec<f64> = perm.iter().map(|&o| if let BoundaryCondition::Flux(v) = bcs[o] { v } else { 0.0 }).collect();
    // rhs = -H u_D + G q_N
    let hu = hm.matvec(&u_known);
    let gq = gm.matvec(&q_known);
    let rhs: Vec<f64> = (0..n).map(|i| gq[i] - hu[i]).collect();
    // operator: x_j = q_j on Dirichlet elements, u_j on Neumann elements: A x = H x_N - G x_D
    let apply = |x: &[f64], y: &mut [f64]| {
        let xd: Vec<f64> = (0..n).map(|i| if is_dir[i] { x[i] } else { 0.0 }).collect();
        let xn: Vec<f64> = (0..n).map(|i| if is_dir[i] { 0.0 } else { x[i] }).collect();
        let gx = gm.matvec(&xd);
        let hx = hm.matvec(&xn);
        for i in 0..n {
            y[i] = hx[i] - gx[i];
        }
    };
    let diag: Vec<f64> = (0..n)
        .map(|i| {
            let o = perm[i];
            if is_dir[i] { -g_entry(o, o) } else { 0.5 }
        })
        .collect();
    let precond = |r: &[f64], z: &mut [f64]| {
        for i in 0..n {
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
