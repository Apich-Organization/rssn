//! # Geometric multigrid for Poisson / Helmholtz problems in 1D, 2D and 3D
//!
//! Solves `-Δu + σ u = f` on a [`structured::Grid`](crate::sim::structured::Grid)
//! with the boundary conditions of [`crate::sim::boundary`] (Dirichlet,
//! Neumann, Robin, periodic, mixed per face, time-dependent data).
//!
//! * coarse operators are re-discretisations with homogeneous data (the
//!   coarse levels solve residual equations);
//! * smoother: Gauss–Seidel (forward before, backward after the coarse
//!   correction);
//! * transfer: multilinear interpolation and full-weighting restriction
//!   (tensor product, with boundary-consistent weights at non-Dirichlet
//!   faces); an axis is coarsened when it has an odd node count `>= 5`
//!   (periodic: an even count `>= 6`), otherwise it is kept;
//! * cycles: V (`γ = 1`), W (`γ = 2`) and full multigrid ([`MgCycle::Fmg`]);
//! * coarsest level solved with the sparse LU of [`crate::kernels::krylov`].
//!
//! `σ <= 0` Helmholtz problems are solvable while the coarse-grid problems
//! stay nonsingular (σ below the first eigenvalue of the coarsest level).
#![allow(
    clippy::manual_memcpy,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_lossless,
    clippy::single_match_else,
    clippy::manual_contains,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::missing_const_for_fn,
    clippy::float_cmp
)]

use crate::kernels::krylov::{Csr, SparseLu, SparseOrdering};
use crate::sim::boundary::FaceBcs;
use crate::sim::structured::{Grid, LaplacianOp};

/// Multigrid cycle type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MgCycle {
    /// V-cycle (one coarse visit).
    V,
    /// W-cycle (two coarse visits).
    W,
    /// Full multigrid: nested iteration followed by V-cycles.
    Fmg,
}

/// Multigrid options.
#[derive(Clone, Copy, Debug)]
pub struct MgOptions {
    /// Cycle type.
    pub cycle: MgCycle,
    /// Pre-smoothing sweeps.
    pub pre: usize,
    /// Post-smoothing sweeps.
    pub post: usize,
    /// Maximum cycles.
    pub max_cycles: usize,
    /// Relative residual tolerance.
    pub tol: f64,
    /// Stop coarsening at or below this many nodes.
    pub coarsest_nodes: usize,
}

impl Default for MgOptions {
    fn default() -> Self {
        Self { cycle: MgCycle::V, pre: 2, post: 2, max_cycles: 50, tol: 1e-10, coarsest_nodes: 64 }
    }
}

/// Result of a multigrid solve.
#[derive(Clone, Debug)]
pub struct MgResult {
    /// Solution at the nodes.
    pub u: Vec<f64>,
    /// Relative residual after each cycle (index 0 is the initial residual).
    pub residuals: Vec<f64>,
    /// Whether `tol` was reached.
    pub converged: bool,
}

impl MgResult {
    /// Geometric-mean residual reduction per cycle.
    #[must_use]
    pub fn convergence_factor(&self) -> f64 {
        let n = self.residuals.len().saturating_sub(1);
        if n == 0 {
            return 0.0;
        }
        let (a, b) = (self.residuals[0], self.residuals[n]);
        (b / a).powf(1.0 / n as f64)
    }
}

struct Level {
    grid: Grid,
    op: LaplacianOp,
    a: Csr,
    /// Axes coarsened when going to the next level.
    coarsen: [bool; 3],
}

/// Geometric multigrid hierarchy for `-Δu + σ u = f`.
pub struct MultigridSolver {
    levels: Vec<Level>,
    coarse_lu: SparseLu,
    opts: MgOptions,
}

fn coarse_n(n: usize, periodic: bool) -> Option<usize> {
    if periodic {
        (n.is_multiple_of(2) && n >= 6).then_some(n / 2)
    } else {
        (n % 2 == 1 && n >= 5).then_some((n - 1) / 2 + 1)
    }
}

/// Applies a 1D map along `axis` of an array with extents `ns`, producing extent `nd` there.
fn along_axis(src: &[f64], ns: [usize; 3], axis: usize, nd: usize, f: &dyn Fn(&[f64], &mut [f64])) -> Vec<f64> {
    let mut nd_all = ns;
    nd_all[axis] = nd;
    let mut out = vec![0.0; nd_all.iter().product()];
    let stride_s: usize = ns[..axis].iter().product();
    let stride_d: usize = nd_all[..axis].iter().product();
    let (o1, o2) = match axis {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let mut line_in = vec![0.0; ns[axis]];
    let mut line_out = vec![0.0; nd];
    for b in 0..ns[o2] {
        for a in 0..ns[o1] {
            let mut idx = [0usize; 3];
            idx[o1] = a;
            idx[o2] = b;
            let base_s = idx[0] + ns[0] * (idx[1] + ns[1] * idx[2]);
            let base_d = idx[0] + nd_all[0] * (idx[1] + nd_all[1] * idx[2]);
            for i in 0..ns[axis] {
                line_in[i] = src[base_s + i * stride_s];
            }
            f(&line_in, &mut line_out);
            for i in 0..nd {
                out[base_d + i * stride_d] = line_out[i];
            }
        }
    }
    out
}

fn restrict_line(fine: &[f64], coarse: &mut [f64], periodic: bool) {
    let nf = fine.len();
    let nc = coarse.len();
    for j in 0..nc {
        let c = 2 * j;
        if periodic {
            let (m, p) = ((c + nf - 1) % nf, (c + 1) % nf);
            coarse[j] = 0.25 * fine[m] + 0.5 * fine[c] + 0.25 * fine[p];
        } else if j == 0 {
            coarse[j] = (2.0 * fine[0] + fine[1]) / 3.0;
        } else if j == nc - 1 {
            coarse[j] = (2.0 * fine[nf - 1] + fine[nf - 2]) / 3.0;
        } else {
            coarse[j] = 0.25 * fine[c - 1] + 0.5 * fine[c] + 0.25 * fine[c + 1];
        }
    }
}

fn prolong_line(coarse: &[f64], fine: &mut [f64], periodic: bool) {
    let nc = coarse.len();
    for j in 0..nc {
        fine[2 * j] = coarse[j];
        if 2 * j + 1 < fine.len() {
            let nxt = if periodic { coarse[(j + 1) % nc] } else { coarse[j + 1] };
            fine[2 * j + 1] = f64::midpoint(coarse[j], nxt);
        }
    }
}

impl MultigridSolver {
    /// Builds the hierarchy for `-Δu + σ u = f`.
    ///
    /// # Errors
    /// On invalid boundary descriptions or a singular coarsest problem.
    pub fn new(grid: &Grid, bcs: &FaceBcs, sigma: f64, opts: MgOptions) -> Result<Self, String> {
        let hom = bcs.homogeneous();
        let mut levels: Vec<Level> = Vec::new();
        let mut g = grid.clone();
        let mut first = true;
        loop {
            let op = LaplacianOp::new(&g, if first { bcs } else { &hom })?;
            if op.has_absorbing() {
                return Err("absorbing boundaries are only valid for waves".into());
            }
            let a = op.system(sigma, -1.0);
            let mut coarsen = [false; 3];
            let mut nn = g.n;
            if g.len() > opts.coarsest_nodes {
                for ax in 0..g.dim {
                    if let Some(c) = coarse_n(g.n[ax], g.periodic[ax]) {
                        coarsen[ax] = true;
                        nn[ax] = c;
                    }
                }
            }
            let any = coarsen.iter().any(|&c| c);
            levels.push(Level { grid: g.clone(), op, a, coarsen });
            first = false;
            if !any {
                break;
            }
            let mut ng = g.clone();
            ng.n = nn;
            g = ng;
        }
        let last = levels.last().ok_or("empty hierarchy")?;
        let coarse_lu = last
            .a
            .sparse_lu(SparseOrdering::MinDegree)
            .map_err(|e| format!("coarsest-level LU failed: {e:?}"))?;
        Ok(Self { levels, coarse_lu, opts })
    }

    /// Number of levels.
    #[must_use]
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    fn smooth(&self, l: usize, u: &mut [f64], b: &[f64], backward: bool) {
        let a = &self.levels[l].a;
        let n = a.rows;
        let mut sweep = |i: usize| {
            let mut s = b[i];
            let mut d = 1.0;
            for k in a.indptr[i]..a.indptr[i + 1] {
                let j = a.indices[k];
                if j == i {
                    d = a.values[k];
                } else {
                    s -= a.values[k] * u[j];
                }
            }
            u[i] = s / d;
        };
        if backward {
            for i in (0..n).rev() {
                sweep(i);
            }
        } else {
            for i in 0..n {
                sweep(i);
            }
        }
    }

    fn residual(&self, l: usize, u: &[f64], b: &[f64]) -> Vec<f64> {
        let lev = &self.levels[l];
        let mut r = vec![0.0; u.len()];
        lev.a.matvec(u, &mut r);
        for p in 0..r.len() {
            r[p] = if lev.op.is_dirichlet(p) { 0.0 } else { b[p] - r[p] };
        }
        r
    }

    fn restrict(&self, l: usize, r: &[f64]) -> Vec<f64> {
        let lev = &self.levels[l];
        let mut cur = r.to_vec();
        let mut ext = lev.grid.n;
        for ax in 0..lev.grid.dim {
            if lev.coarsen[ax] {
                let per = lev.grid.periodic[ax];
                let nc = coarse_n(ext[ax], per).unwrap_or(ext[ax]);
                cur = along_axis(&cur, ext, ax, nc, &|f, c| restrict_line(f, c, per));
                ext[ax] = nc;
            }
        }
        let coarse = &self.levels[l + 1];
        for (p, v) in cur.iter_mut().enumerate() {
            if coarse.op.is_dirichlet(p) {
                *v = 0.0;
            }
        }
        cur
    }

    fn prolong(&self, l: usize, e: &[f64]) -> Vec<f64> {
        let fine = &self.levels[l];
        let mut ext = self.levels[l + 1].grid.n;
        let mut cur = e.to_vec();
        for ax in 0..fine.grid.dim {
            if fine.coarsen[ax] {
                let per = fine.grid.periodic[ax];
                cur = along_axis(&cur, ext, ax, fine.grid.n[ax], &|c, f| prolong_line(c, f, per));
                ext[ax] = fine.grid.n[ax];
            }
        }
        cur
    }

    fn cycle(&self, l: usize, u: &mut [f64], b: &[f64], gamma: usize) {
        if l + 1 == self.levels.len() {
            let sol = self.coarse_lu.solve(b);
            u.copy_from_slice(&sol);
            return;
        }
        for _ in 0..self.opts.pre {
            self.smooth(l, u, b, false);
        }
        let r = self.residual(l, u, b);
        let rc = self.restrict(l, &r);
        let mut ec = vec![0.0; rc.len()];
        for _ in 0..gamma {
            self.cycle(l + 1, &mut ec, &rc, gamma);
        }
        let corr = self.prolong(l, &ec);
        for p in 0..u.len() {
            if !self.levels[l].op.is_dirichlet(p) {
                u[p] += corr[p];
            }
        }
        for _ in 0..self.opts.post {
            self.smooth(l, u, b, true);
        }
    }

    /// Full-multigrid pass for the error equation `A e = r` (homogeneous data).
    fn fmg(&self, l: usize, r: &[f64]) -> Vec<f64> {
        if l + 1 == self.levels.len() {
            return self.coarse_lu.solve(r);
        }
        let rc = self.restrict(l, r);
        let ec = self.fmg(l + 1, &rc);
        let mut e = self.prolong(l, &ec);
        for p in 0..e.len() {
            if self.levels[l].op.is_dirichlet(p) {
                e[p] = 0.0;
            }
        }
        self.cycle(l, &mut e, r, 1);
        e
    }

    /// Solves with source `f` at the nodes and boundary data evaluated at time `t`.
    ///
    /// # Errors
    /// If `f` has the wrong length.
    pub fn solve(&self, f: &[f64], t: f64) -> Result<MgResult, String> {
        let lev0 = &self.levels[0];
        let nn = lev0.grid.len();
        if f.len() != nn {
            return Err("source has the wrong length".into());
        }
        let force = lev0.op.forcing(t);
        let mut b: Vec<f64> = (0..nn).map(|p| f[p] + force[p]).collect();
        lev0.op.set_dirichlet(&mut b, t);
        let bn = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
        let mut u = vec![0.0; nn];
        lev0.op.set_dirichlet(&mut u, t);
        let rn = |this: &Self, u: &[f64]| this.residual(0, u, &b).iter().map(|v| v * v).sum::<f64>().sqrt() / bn;
        let mut residuals = vec![rn(self, &u)];
        if self.opts.cycle == MgCycle::Fmg {
            let r = self.residual(0, &u, &b);
            let e = self.fmg(0, &r);
            for p in 0..nn {
                u[p] += e[p];
            }
            residuals.push(rn(self, &u));
        }
        let gamma = if self.opts.cycle == MgCycle::W { 2 } else { 1 };
        let mut converged = residuals.last().is_some_and(|&r| r <= self.opts.tol);
        let mut k = 0;
        while !converged && k < self.opts.max_cycles {
            self.cycle(0, &mut u, &b, gamma);
            let r = rn(self, &u);
            residuals.push(r);
            converged = r <= self.opts.tol;
            k += 1;
        }
        Ok(MgResult { u, residuals, converged })
    }
}

/// Convenience wrapper: solves `-Δu + σu = f(x)` with multigrid.
///
/// # Errors
/// See [`MultigridSolver::new`].
pub fn solve_helmholtz_multigrid(
    grid: &Grid,
    bcs: &FaceBcs,
    sigma: f64,
    f: impl Fn([f64; 3]) -> f64,
    opts: MgOptions,
) -> Result<MgResult, String> {
    let mg = MultigridSolver::new(grid, bcs, sigma, opts)?;
    mg.solve(&grid.sample(f), 0.0)
}
