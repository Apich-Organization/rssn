//! Geometric multigrid (V/W/FMG) for Poisson/Helmholtz in 1D/2D/3D.
//!
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
//! Each solve takes a separate lease from the shared scratch pool for tensor
//! transfers, residuals and coarse corrections. Residual reports mask
//! Dirichlet rows and use a scaled norm ratio when the ordinary squared norms
//! are outside the normal floating-point range.
//!
//! `σ <= 0` Helmholtz problems are solvable while the coarse-grid problems
//! stay nonsingular (σ below the first eigenvalue of the coarsest level).

use crate::kernels::krylov::Csr;
use crate::kernels::krylov::SparseLu;
use crate::kernels::krylov::SparseOrdering;
use crate::kernels::krylov::residual_ratio;
use crate::kernels::reuse::f64_scratch;
use crate::sim::boundary::FaceBcs;
use crate::sim::structured::Grid;
use crate::sim::structured::LaplacianOp;

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
        Self {
            cycle: MgCycle::V,
            pre: 2,
            post: 2,
            max_cycles: 50,
            tol: 1e-10,
            coarsest_nodes: 64,
        }
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
    /// Returns zero if there was no cycle or the initial residual was zero.
    #[must_use]
    pub fn convergence_factor(&self) -> f64 {
        let n = self.residuals.len().saturating_sub(1);
        if n == 0 {
            return 0.0;
        }
        let (a, b) = (self.residuals[0], self.residuals[n]);
        if a == 0.0 {
            return 0.0;
        }
        let power = (n as f64).recip();
        let ratio = b / a;
        if ratio.is_normal() || b == 0.0 || n == 1 {
            ratio.powf(power)
        } else {
            b.powf(power) / a.powf(power)
        }
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

fn coarse_n(
    n: usize,
    periodic: bool,
) -> Option<usize> {
    if periodic {
        (n.is_multiple_of(2) && n >= 6).then_some(n / 2)
    } else {
        (n % 2 == 1 && n >= 5).then_some((n - 1) / 2 + 1)
    }
}

struct TransferWork<'a> {
    first: &'a mut [f64],
    second: &'a mut [f64],
    line_in: &'a mut [f64],
    line_out: &'a mut [f64],
}

/// Exact buffer sizes for tensor intermediates and for one gathered line.
fn transfer_sizes(
    fine: &Level,
    coarse: Option<&Level>,
) -> ([usize; 2], usize) {
    let Some(coarse) = coarse else {
        return ([0; 2], 0);
    };
    let (mut f, mut c) = (fine.grid.n, coarse.grid.n);
    let mut axes = (0..fine.grid.dim)
        .filter(|&axis| fine.coarsen[axis])
        .peekable();
    let mut sizes = [0; 2];
    let mut step = 0;
    let mut line = 0;
    while let Some(axis) = axes.next() {
        if axis > 0 {
            line = line.max(f[axis]);
        }
        if axes.peek().is_some() {
            std::mem::swap(&mut f[axis], &mut c[axis]);
            sizes[step] = f.iter().product::<usize>().max(c.iter().product::<usize>());
            step += 1;
        }
    }
    // Every later shape is bounded by the first transition's coarse shape.
    (sizes, line)
}

struct LevelWork<'a> {
    residual: &'a mut [f64],
    coarse_rhs: &'a mut [f64],
    coarse_correction: &'a mut [f64],
    child: &'a mut [f64],
}

impl<'a> LevelWork<'a> {
    const fn split(
        work: &'a mut [f64],
        fine: usize,
        coarse: usize,
    ) -> Self {
        let (residual, rest) = work.split_at_mut(fine);
        let (coarse_rhs, rest) = rest.split_at_mut(coarse);
        let (coarse_correction, child) = rest.split_at_mut(coarse);
        Self {
            residual,
            coarse_rhs,
            coarse_correction,
            child,
        }
    }
}

/// Applies a 1D map along `axis` into the caller's finer or coarser slice.
fn along_axis(
    src: &[f64],
    ns: [usize; 3],
    axis: usize,
    nd: usize,
    out: &mut [f64],
    (line_in, line_out): (&mut [f64], &mut [f64]),
    f: impl Fn(&[f64], &mut [f64]),
) {
    if axis == 0 {
        let lines = ns[1] * ns[2];
        for (source, dest) in src[..ns[0] * lines]
            .chunks_exact(ns[0])
            .zip(out[..nd * lines].chunks_exact_mut(nd))
        {
            f(source, dest);
        }
        return;
    }

    let mut nd_all = ns;
    nd_all[axis] = nd;
    let stride = ns[..axis].iter().product::<usize>();
    let (o1, o2) = if axis == 1 { (0, 2) } else { (0, 1) };
    let line_in = &mut line_in[..ns[axis]];
    let line_out = &mut line_out[..nd];
    for b in 0..ns[o2] {
        for a in 0..ns[o1] {
            let mut idx = [0usize; 3];
            idx[o1] = a;
            idx[o2] = b;
            let base_s = idx[0] + ns[0] * (idx[1] + ns[1] * idx[2]);
            let base_d = idx[0] + nd_all[0] * (idx[1] + nd_all[1] * idx[2]);
            for (i, val) in line_in.iter_mut().enumerate() {
                *val = src[base_s + i * stride];
            }
            f(line_in, line_out);
            for (i, &val) in line_out.iter().enumerate() {
                out[base_d + i * stride] = val;
            }
        }
    }
}

fn tensor_transfer(
    src: &[f64],
    mut ext: [usize; 3],
    dst: &mut [f64],
    target: [usize; 3],
    fine: &Level,
    work: &mut TransferWork<'_>,
    line: impl Fn(&[f64], &mut [f64], bool),
) {
    let TransferWork {
        first,
        second,
        line_in,
        line_out,
    } = work;
    let mut step = |axis: usize, input: &[f64], output: &mut [f64]| {
        along_axis(
            input,
            ext,
            axis,
            target[axis],
            output,
            (line_in, line_out),
            |a, b| {
                line(a, b, fine.grid.periodic[axis]);
            },
        );
        ext[axis] = target[axis];
    };
    let mut axes = (0..fine.grid.dim).filter(|&axis| fine.coarsen[axis]);
    let Some(a) = axes.next() else {
        dst.copy_from_slice(src);
        return;
    };
    let Some(b) = axes.next() else {
        step(a, src, dst);
        return;
    };
    step(a, src, first);
    if let Some(c) = axes.next() {
        step(b, first, second);
        step(c, second, dst);
    } else {
        step(b, first, dst);
    }
}

fn restrict_line(
    fine: &[f64],
    coarse: &mut [f64],
    periodic: bool,
) {
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

fn prolong_line(
    coarse: &[f64],
    fine: &mut [f64],
    periodic: bool,
) {
    let nc = coarse.len();
    for j in 0..nc {
        fine[2 * j] = coarse[j];
        if 2 * j + 1 < fine.len() {
            let nxt = if periodic {
                coarse[(j + 1) % nc]
            } else {
                coarse[j + 1]
            };
            fine[2 * j + 1] = f64::midpoint(coarse[j], nxt);
        }
    }
}

impl MultigridSolver {
    /// Builds the hierarchy for `-Δu + σ u = f`.
    ///
    /// # Errors
    /// On invalid coefficients, boundary data or a singular coarsest problem.
    pub fn new(
        grid: &Grid,
        bcs: &FaceBcs,
        sigma: f64,
        opts: MgOptions,
    ) -> Result<Self, String> {
        if !sigma.is_finite() || !opts.tol.is_finite() || opts.tol < 0.0 {
            return Err(
                "multigrid requires finite sigma and a finite nonnegative tolerance".into(),
            );
        }
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
            levels.push(Level {
                grid: g.clone(),
                op,
                a,
                coarsen,
            });
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
        Ok(Self {
            levels,
            coarse_lu,
            opts,
        })
    }

    /// Number of levels.
    #[must_use]
    pub const fn num_levels(&self) -> usize {
        self.levels.len()
    }

    fn smooth(
        &self,
        l: usize,
        u: &mut [f64],
        b: &[f64],
        backward: bool,
    ) {
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

    fn residual_into(
        &self,
        l: usize,
        u: &[f64],
        b: &[f64],
        out: &mut [f64],
    ) {
        let lev = &self.levels[l];
        lev.a.matvec(u, out);
        for (p, v) in out.iter_mut().enumerate() {
            *v = if lev.op.is_dirichlet(p) {
                0.0
            } else {
                b[p] - *v
            };
        }
    }

    fn relative_residual(
        &self,
        u: &[f64],
        b: &[f64],
        regular_rhs_norm: Option<f64>,
        ax: &mut [f64],
    ) -> Result<f64, String> {
        let lev = &self.levels[0];
        lev.a.matvec(u, ax);
        if let Some(bn) = regular_rhs_norm {
            let mut square = 0.0;
            let mut nonzero = false;
            for (p, (&rhs, &a_u)) in b.iter().zip(ax.iter()).enumerate() {
                if !lev.op.is_dirichlet(p) {
                    let r = rhs - a_u;
                    square += r * r;
                    nonzero |= r != 0.0;
                }
            }
            if square.is_normal() || !nonzero {
                return Ok(square.sqrt() / bn);
            }
        }
        for (p, val) in ax.iter_mut().enumerate() {
            if lev.op.is_dirichlet(p) {
                *val = b[p];
            }
        }
        residual_ratio(b, ax)
            .filter(|r| r.is_finite())
            .ok_or_else(|| "multigrid residual exceeds the finite range".into())
    }

    fn restrict_into(
        &self,
        l: usize,
        r: &[f64],
        out: &mut [f64],
        transfer: &mut TransferWork<'_>,
    ) {
        let fine = &self.levels[l];
        let coarse = &self.levels[l + 1];
        tensor_transfer(
            r,
            fine.grid.n,
            out,
            coarse.grid.n,
            fine,
            transfer,
            restrict_line,
        );
        for (p, val) in out.iter_mut().enumerate() {
            if coarse.op.is_dirichlet(p) {
                *val = 0.0;
            }
        }
    }

    fn prolong_into(
        &self,
        l: usize,
        e: &[f64],
        out: &mut [f64],
        transfer: &mut TransferWork<'_>,
    ) {
        let fine = &self.levels[l];
        tensor_transfer(
            e,
            self.levels[l + 1].grid.n,
            out,
            fine.grid.n,
            fine,
            transfer,
            prolong_line,
        );
    }

    fn cycle(
        &self,
        l: usize,
        u: &mut [f64],
        b: &[f64],
        gamma: usize,
        work: &mut [f64],
        transfer: &mut TransferWork<'_>,
    ) {
        if l + 1 == self.levels.len() {
            let sol = self.coarse_lu.solve(b);
            u.copy_from_slice(&sol);
            return;
        }
        for _ in 0..self.opts.pre {
            self.smooth(l, u, b, false);
        }
        let LevelWork {
            residual,
            coarse_rhs,
            coarse_correction,
            child,
        } = LevelWork::split(work, u.len(), self.levels[l + 1].grid.len());
        self.residual_into(l, u, b, residual);
        self.restrict_into(l, residual, coarse_rhs, transfer);
        coarse_correction.fill(0.0);
        for _ in 0..gamma {
            self.cycle(l + 1, coarse_correction, coarse_rhs, gamma, child, transfer);
        }
        // Restriction has already consumed the fine residual.
        self.prolong_into(l, coarse_correction, residual, transfer);
        for (p, (value, corr)) in u.iter_mut().zip(residual.iter()).enumerate() {
            if !self.levels[l].op.is_dirichlet(p) {
                *value += corr;
            }
        }
        for _ in 0..self.opts.post {
            self.smooth(l, u, b, true);
        }
    }

    /// Full-multigrid pass for the error equation `A e = r` (homogeneous data).
    fn fmg(
        &self,
        l: usize,
        r: &[f64],
        e: &mut [f64],
        work: &mut [f64],
        transfer: &mut TransferWork<'_>,
    ) {
        if l + 1 == self.levels.len() {
            e.copy_from_slice(&self.coarse_lu.solve(r));
            return;
        }
        let LevelWork {
            coarse_rhs,
            coarse_correction,
            child,
            ..
        } = LevelWork::split(work, r.len(), self.levels[l + 1].grid.len());
        self.restrict_into(l, r, coarse_rhs, transfer);
        self.fmg(l + 1, coarse_rhs, coarse_correction, child, transfer);
        self.prolong_into(l, coarse_correction, e, transfer);
        for (p, val) in e.iter_mut().enumerate() {
            if self.levels[l].op.is_dirichlet(p) {
                *val = 0.0;
            }
        }
        self.cycle(l, e, r, 1, work, transfer);
    }

    /// Solves with source `f` at the nodes and boundary data evaluated at time `t`.
    ///
    /// # Errors
    /// If the source length is wrong, the time or active data are nonfinite,
    /// or the residual, solution or workspace cannot be represented.
    pub fn solve(
        &self,
        f: &[f64],
        t: f64,
    ) -> Result<MgResult, String> {
        let lev0 = &self.levels[0];
        let nn = lev0.grid.len();
        if f.len() != nn {
            return Err("source has the wrong length".into());
        }
        if !t.is_finite() {
            return Err("multigrid time must be finite".into());
        }
        let mut b = lev0.op.forcing(t);
        for (rhs, &source) in b.iter_mut().zip(f) {
            *rhs += source;
        }
        lev0.op.set_dirichlet(&mut b, t);
        if b.iter().any(|v| !v.is_finite()) {
            return Err("multigrid source and evaluated boundary data must be finite".into());
        }
        let rhs_square = b.iter().map(|v| v * v).sum::<f64>();
        let regular_rhs_norm = rhs_square.is_normal().then(|| rhs_square.sqrt());
        let mut u = vec![0.0; nn];
        lev0.op.set_dirichlet(&mut u, t);

        let cycle_items = self
            .levels
            .windows(2)
            .try_fold(0usize, |size, pair| {
                size.checked_add(pair[0].grid.len())?
                    .checked_add(pair[1].grid.len().checked_mul(2)?)
            })
            .ok_or("multigrid workspace size overflow")?
            .max(nn);
        let ([axis_first, axis_second], line) = transfer_sizes(lev0, self.levels.get(1));
        let fmg = self.opts.cycle == MgCycle::Fmg;
        let fmg_space = if fmg {
            nn.checked_mul(2)
                .ok_or("multigrid workspace size overflow")?
        } else {
            0
        };
        let items = [axis_first, axis_second, line, line, fmg_space]
            .into_iter()
            .try_fold(cycle_items, usize::checked_add)
            .filter(|&size| size <= isize::MAX as usize / size_of::<f64>())
            .ok_or("multigrid workspace exceeds addressable storage")?;
        let mut buffer = f64_scratch(items);
        buffer.resize(items, 0.0);
        let (work, rest) = buffer.split_at_mut(cycle_items);
        let (first, rest) = rest.split_at_mut(axis_first);
        let (second, rest) = rest.split_at_mut(axis_second);
        let (line_in, rest) = rest.split_at_mut(line);
        let (line_out, extra) = rest.split_at_mut(line);
        let mut transfer = TransferWork {
            first,
            second,
            line_in,
            line_out,
        };

        let rn = |u: &[f64], work: &mut [f64]| {
            self.relative_residual(u, &b, regular_rhs_norm, &mut work[..nn])
        };
        let mut r = rn(&u, work)?;
        let mut residuals = vec![r];
        if fmg {
            let (rhs, e) = extra.split_at_mut(nn);
            self.residual_into(0, &u, &b, rhs);
            self.fmg(0, rhs, e, work, &mut transfer);
            for (value, correction) in u.iter_mut().zip(e) {
                *value += *correction;
            }
            r = rn(&u, work)?;
            residuals.push(r);
        }
        let gamma = if self.opts.cycle == MgCycle::W {
            2
        } else {
            1
        };
        let mut k = 0;
        while r > self.opts.tol && k < self.opts.max_cycles {
            self.cycle(0, &mut u, &b, gamma, work, &mut transfer);
            r = rn(&u, work)?;
            residuals.push(r);
            k += 1;
        }
        if u.iter().any(|v| !v.is_finite()) {
            return Err("multigrid solution exceeds the finite range".into());
        }
        Ok(MgResult {
            u,
            residuals,
            converged: r <= self.opts.tol,
        })
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
