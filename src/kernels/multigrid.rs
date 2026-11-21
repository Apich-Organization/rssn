//! Geometric multigrid (V/W/FMG) and smoothed-aggregation algebraic multigrid.
//!
//! # Multigrid
//!
//! * [`GeometricMg`]: geometric multigrid for the Poisson-type operator
//!   `-Laplacian u + shift u` with homogeneous Dirichlet conditions on the
//!   unit interval, square or cube (1, 2 or 3 dimensions, `n = 2^k - 1`
//!   interior points per direction, second-order 3/5/7-point stencils).
//!   Smoothers: weighted Jacobi or red-black Gauss-Seidel; full-weighting
//!   restriction and (bi/tri)linear interpolation; cycles: V, W and the full
//!   multigrid (FMG, nested iteration) scheme;
//! * [`Amg`]: algebraic multigrid by smoothed aggregation (Vanek, Mandel and
//!   Brezina) for sparse symmetric positive definite [`Csr`] matrices:
//!   strength-of-connection filtering, greedy aggregation, tentative
//!   prolongator built from a (user-extensible) near-null space, prolongator
//!   smoothing by one damped Jacobi step, Galerkin coarse operators and a
//!   symmetric Gauss-Seidel smoother. The V/W cycle is a symmetric positive
//!   definite operator and so an ideal preconditioner for CG
//!   ([`Amg::precondition`]).

use crate::kernels::dense::Lu;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::krylov::Csr;
use crate::kernels::krylov::residual_ratio;

/// Multigrid cycle shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cycle {
    /// One coarse-grid correction per level.
    V,
    /// Two coarse-grid corrections per level.
    W,
}

/// Smoother of the geometric multigrid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Smoother {
    /// Weighted Jacobi (`omega = 2/3`, `4/5`, `6/7` for 1, 2, 3 dimensions).
    Jacobi,
    /// Red-black Gauss-Seidel.
    RedBlackGaussSeidel,
}

/// Parameters of the geometric multigrid cycles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MgOptions {
    /// Cycle shape.
    pub cycle: Cycle,
    /// Pre-smoothing sweeps.
    pub pre: usize,
    /// Post-smoothing sweeps.
    pub post: usize,
    /// Smoother.
    pub smoother: Smoother,
    /// Relative residual tolerance of the solve.
    pub tol: f64,
    /// Maximal number of cycles.
    pub max_cycles: usize,
}

impl Default for MgOptions {
    fn default() -> Self {
        Self {
            cycle: Cycle::V,
            pre: 2,
            post: 2,
            smoother: Smoother::RedBlackGaussSeidel,
            tol: 1e-10,
            max_cycles: 100,
        }
    }
}

/// Outcome of a multigrid solve.
#[derive(Debug, Clone, PartialEq)]
pub struct MgResult {
    /// The solution.
    pub x: Vec<f64>,
    /// Relative residual norm after each cycle (first entry: initial).
    pub history: Vec<f64>,
    /// Whether the tolerance was reached.
    pub converged: bool,
}

impl MgResult {
    /// Geometric mean convergence factor over the recorded cycles.
    #[must_use]
    pub fn mean_factor(&self) -> f64 {
        let n = self.history.len();
        if n < 2 || self.history[0] == 0.0 {
            return 0.0;
        }
        (self.history[n - 1] / self.history[0]).powf(1.0 / (n - 1) as f64)
    }
}

/// Geometric multigrid hierarchy for `-Laplacian + shift` on `(0,1)^dim`.
#[derive(Debug, Clone)]
pub struct GeometricMg {
    dim: usize,
    n: usize,
    shift: f64,
    levels: usize,
}

/// Reusable intermediate storage for adjacent geometric grid transfers.
///
/// A workspace may be reused across hierarchies, levels and directions. It
/// grows only when a larger intermediate grid is needed.
#[derive(Debug, Default)]
pub struct GeometricTransferWorkspace {
    first: Vec<f64>,
    second: Vec<f64>,
}

/// Invalid level or buffer sizes for a geometric grid transfer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometricTransferError {
    /// The requested fine level has no next coarser level.
    Level,
    /// A buffer does not have the expected number of unknowns.
    Shape,
}

impl std::fmt::Display for GeometricTransferError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::Level => "geometric transfer needs an adjacent coarser level",
            | Self::Shape => "geometric transfer buffer length does not match its grid",
        })
    }
}

impl std::error::Error for GeometricTransferError {}

impl GeometricTransferWorkspace {
    pub(crate) const fn heap_bytes(&self) -> usize {
        self.first
            .capacity()
            .saturating_add(self.second.capacity())
            .saturating_mul(std::mem::size_of::<f64>())
    }

    fn transfer(
        &mut self,
        dim: usize,
        from_n: usize,
        to_n: usize,
        from: &[f64],
        to: &mut [f64],
        axis: impl Fn(&[f64], &mut [f64], usize, usize, usize),
    ) {
        if dim == 1 {
            axis(from, to, from_n, to_n, 1);
            return;
        }
        let first_size = to_n * from_n.pow((dim - 1) as u32);
        if self.first.len() < first_size {
            self.first.resize(first_size, 0.0);
        }
        let first = &mut self.first[..first_size];
        axis(from, first, from_n, to_n, 1);
        if dim == 2 {
            axis(first, to, from_n, to_n, to_n);
            return;
        }
        let second_size = to_n * to_n * from_n;
        if self.second.len() < second_size {
            self.second.resize(second_size, 0.0);
        }
        let second = &mut self.second[..second_size];
        axis(first, second, from_n, to_n, to_n);
        axis(second, to, from_n, to_n, to_n * to_n);
    }
}

fn full_weight(
    left: f64,
    center: f64,
    right: f64,
) -> f64 {
    let weighted = 0.5 * center + 0.25 * left + 0.25 * right;
    if weighted.abs() <= f64::MIN_POSITIVE {
        let peak = left.abs().max(center.abs()).max(right.abs());
        if peak > 0.0 && peak <= f64::MIN_POSITIVE {
            // These three values are integral multiples of the smallest
            // subnormal. Form the exact quarter-weighted numerator and round
            // once to the nearest multiple, with ties to even.
            let units = |value: f64| {
                let magnitude = (value.to_bits() & (u64::MAX >> 1)) as i64;
                if value.is_sign_negative() {
                    -magnitude
                } else {
                    magnitude
                }
            };
            let numerator = units(left) + 2 * units(center) + units(right);
            let magnitude = numerator.unsigned_abs();
            let whole = magnitude / 4;
            let remainder = magnitude % 4;
            let rounded = whole + u64::from(remainder > 2 || (remainder == 2 && whole & 1 != 0));
            let sign = if numerator < 0 {
                1_u64 << 63
            } else {
                0
            };
            return f64::from_bits(sign | rounded);
        }
    }
    weighted
}

fn restrict_axis(
    from: &[f64],
    to: &mut [f64],
    nf: usize,
    nc: usize,
    inner: usize,
) {
    for (fine, coarse) in from
        .chunks_exact(nf * inner)
        .zip(to.chunks_exact_mut(nc * inner))
    {
        for (j, column) in coarse.chunks_exact_mut(inner).enumerate() {
            let start = 2 * j * inner;
            for (k, slot) in column.iter_mut().enumerate() {
                *slot = full_weight(
                    fine[start + k],
                    fine[start + inner + k],
                    fine[start + 2 * inner + k],
                );
            }
        }
    }
}

fn prolong_axis(
    from: &[f64],
    to: &mut [f64],
    nc: usize,
    nf: usize,
    inner: usize,
) {
    for (coarse, fine) in from
        .chunks_exact(nc * inner)
        .zip(to.chunks_exact_mut(nf * inner))
    {
        for (j, column) in coarse.chunks_exact(inner).enumerate() {
            let offset = 2 * j * inner;
            let midpoint = &mut fine[offset..offset + inner];
            if j == 0 {
                for (slot, &right) in midpoint.iter_mut().zip(column) {
                    *slot = f64::midpoint(0.0, right);
                }
            } else {
                let left_column = &coarse[(j - 1) * inner..j * inner];
                for ((slot, &left), &right) in midpoint.iter_mut().zip(left_column).zip(column) {
                    *slot = f64::midpoint(left, right);
                }
            }
            fine[offset + inner..offset + 2 * inner].copy_from_slice(column);
        }
        let last = &coarse[(nc - 1) * inner..];
        for (slot, &left) in fine[(nf - 1) * inner..].iter_mut().zip(last) {
            *slot = f64::midpoint(left, 0.0);
        }
    }
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

const fn ordinary_square_sum(sum: f64) -> bool {
    // Leave headroom so underflowing terms cannot affect the ordinary norm.
    sum >= 1e-270 && sum.is_finite()
}

fn fine_residuals<'a>(
    rhs: &'a [f64],
    applied: &'a [f64],
) -> impl Iterator<Item = f64> + Clone + 'a {
    rhs.iter().zip(applied).map(|(rhs, applied)| rhs - applied)
}

fn stable_residual_norm(
    residuals: impl Iterator<Item = f64> + Clone,
    square_sum: f64,
) -> f64 {
    if ordinary_square_sum(square_sum) || square_sum.is_nan() {
        return square_sum.sqrt();
    }
    let peak = residuals.clone().map(f64::abs).fold(0.0, f64::max);
    if peak == 0.0 || !peak.is_finite() {
        return peak;
    }
    let normalized = residuals.map(|r| (r / peak).powi(2)).sum::<f64>().sqrt();
    peak * normalized
}

impl GeometricMg {
    /// Creates the hierarchy for `n = 2^k - 1` interior points per direction.
    ///
    /// Returns `None` when `dim` is not 1, 2 or 3, when `n` is not of the
    /// form `2^k - 1` (`k >= 1`), when the grid size overflows, or when
    /// `shift` is negative or non-finite.
    #[must_use]
    pub fn new(
        dim: usize,
        n: usize,
        shift: f64,
    ) -> Option<Self> {
        let extent = n.checked_add(1)?;
        if !(1..=3).contains(&dim)
            || n == 0
            || !extent.is_power_of_two()
            || !shift.is_finite()
            || shift < 0.0
        {
            return None;
        }
        n.checked_pow(dim as u32)?;
        let levels = extent.trailing_zeros() as usize;
        Some(Self {
            dim,
            n,
            shift,
            levels,
        })
    }

    /// Number of unknowns on the fine grid.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.n.pow(self.dim as u32)
    }

    /// Number of levels (the coarsest has a single unknown).
    #[must_use]
    pub const fn levels(&self) -> usize {
        self.levels
    }

    /// Number of unknowns at a level, or `None` if the level is absent.
    #[must_use]
    pub const fn level_size(
        &self,
        level: usize,
    ) -> Option<usize> {
        if level < self.levels {
            Some(self.n_at(level).pow(self.dim as u32))
        } else {
            None
        }
    }

    const fn n_at(
        &self,
        level: usize,
    ) -> usize {
        ((self.n + 1) >> level) - 1
    }

    fn h_at(
        &self,
        level: usize,
    ) -> f64 {
        1.0 / ((self.n_at(level) + 1) as f64)
    }

    /// Applies the discrete operator on `level` (0 is the finest grid).
    pub fn apply(
        &self,
        level: usize,
        u: &[f64],
        out: &mut [f64],
    ) {
        let n = self.n_at(level);
        let ih2 = 1.0 / (self.h_at(level) * self.h_at(level));
        let diag = 2.0 * self.dim as f64 * ih2 + self.shift;
        let strides: Vec<usize> = (0..self.dim).map(|d| n.pow(d as u32)).collect();
        for idx in 0..u.len() {
            let mut s = diag * u[idx];
            for d in 0..self.dim {
                let c = (idx / strides[d]) % n;
                if c > 0 {
                    s -= ih2 * u[idx - strides[d]];
                }
                if c + 1 < n {
                    s -= ih2 * u[idx + strides[d]];
                }
            }
            out[idx] = s;
        }
    }

    fn smooth(
        &self,
        level: usize,
        u: &mut [f64],
        f: &[f64],
        sweeps: usize,
        sm: Smoother,
    ) {
        if sweeps == 0 {
            return;
        }
        let n = self.n_at(level);
        let ih2 = 1.0 / (self.h_at(level) * self.h_at(level));
        let diag = 2.0 * self.dim as f64 * ih2 + self.shift;
        match sm {
            | Smoother::Jacobi => {
                let omega = match self.dim {
                    | 1 => 2.0 / 3.0,
                    | 2 => 4.0 / 5.0,
                    | _ => 6.0 / 7.0,
                };
                let mut au = vec![0.0; u.len()];
                for _ in 0..sweeps {
                    self.apply(level, u, &mut au);
                    for i in 0..u.len() {
                        u[i] += omega * (f[i] - au[i]) / diag;
                    }
                }
            },
            | Smoother::RedBlackGaussSeidel => {
                let rows = if self.dim >= 2 { n } else { 1 };
                let planes = if self.dim == 3 { n } else { 1 };
                let plane = rows * n;
                for _ in 0..sweeps {
                    for color in 0..2 {
                        for z in 0..planes {
                            for y in 0..rows {
                                let row = z * plane + y * n;
                                // Every grid extent is odd, so a flat index and
                                // the sum of its coordinates have the same parity.
                                let first = color ^ (row & 1);
                                for x in (first..n).step_by(2) {
                                    let idx = row + x;
                                    let mut s = f[idx];
                                    if x > 0 {
                                        s += ih2 * u[idx - 1];
                                    }
                                    if x + 1 < n {
                                        s += ih2 * u[idx + 1];
                                    }
                                    if y > 0 {
                                        s += ih2 * u[idx - n];
                                    }
                                    if y + 1 < rows {
                                        s += ih2 * u[idx + n];
                                    }
                                    if z > 0 {
                                        s += ih2 * u[idx - plane];
                                    }
                                    if z + 1 < planes {
                                        s += ih2 * u[idx + plane];
                                    }
                                    u[idx] = s / diag;
                                }
                            }
                        }
                    }
                }
            },
        }
    }

    /// Restrict a fine grid to the adjacent coarser level by full weighting.
    ///
    /// The stencil is `[1, 2, 1]/4` along each axis. Axis passes run in
    /// x, y, z order on flat buffers with x varying fastest. Each pass
    /// rounds to `f64` before the next. The same workspace may be used
    /// for later transfers in either direction.
    ///
    /// # Errors
    ///
    /// Returns an error for a missing adjacent level or incorrect buffer
    /// lengths. The output is untouched on such errors.
    pub fn restrict_into(
        &self,
        level: usize,
        fine: &[f64],
        coarse: &mut [f64],
        workspace: &mut GeometricTransferWorkspace,
    ) -> Result<(), GeometricTransferError> {
        let (nf, nc) = self.transfer_extents(level, fine.len(), coarse.len())?;
        workspace.transfer(self.dim, nf, nc, fine, coarse, restrict_axis);
        Ok(())
    }

    /// Interpolate from the adjacent coarser level by multilinear interpolation.
    ///
    /// The `level` argument identifies the fine output. Missing boundary
    /// coarse points have value zero, matching homogeneous Dirichlet data.
    ///
    /// # Errors
    ///
    /// Returns an error for a missing adjacent level or incorrect buffer
    /// lengths. The output is untouched on such errors.
    pub fn prolong_into(
        &self,
        level: usize,
        coarse: &[f64],
        fine: &mut [f64],
        workspace: &mut GeometricTransferWorkspace,
    ) -> Result<(), GeometricTransferError> {
        let (nf, nc) = self.transfer_extents(level, fine.len(), coarse.len())?;
        workspace.transfer(self.dim, nc, nf, coarse, fine, prolong_axis);
        Ok(())
    }

    const fn transfer_extents(
        &self,
        level: usize,
        fine_len: usize,
        coarse_len: usize,
    ) -> Result<(usize, usize), GeometricTransferError> {
        if level >= self.levels - 1 {
            return Err(GeometricTransferError::Level);
        }
        let nf = self.n_at(level);
        let nc = self.n_at(level + 1);
        if fine_len != nf.pow(self.dim as u32) || coarse_len != nc.pow(self.dim as u32) {
            return Err(GeometricTransferError::Shape);
        }
        Ok((nf, nc))
    }

    fn restrict(
        &self,
        level: usize,
        fine: &[f64],
        workspace: &mut GeometricTransferWorkspace,
    ) -> Vec<f64> {
        let nc = self.n_at(level + 1);
        let mut coarse = vec![0.0; nc.pow(self.dim as u32)];
        workspace.transfer(
            self.dim,
            self.n_at(level),
            nc,
            fine,
            &mut coarse,
            restrict_axis,
        );
        coarse
    }

    fn prolong(
        &self,
        level: usize,
        coarse: &[f64],
        workspace: &mut GeometricTransferWorkspace,
    ) -> Vec<f64> {
        let nf = self.n_at(level);
        let mut fine = vec![0.0; nf.pow(self.dim as u32)];
        workspace.transfer(
            self.dim,
            self.n_at(level + 1),
            nf,
            coarse,
            &mut fine,
            prolong_axis,
        );
        fine
    }

    fn cycle(
        &self,
        level: usize,
        u: &mut [f64],
        f: &[f64],
        opts: &MgOptions,
        transfers: &mut GeometricTransferWorkspace,
    ) {
        if level + 1 == self.levels {
            // single unknown
            let ih2 = 1.0 / (self.h_at(level) * self.h_at(level));
            u[0] = f[0] / (2.0 * self.dim as f64 * ih2 + self.shift);
            return;
        }
        self.smooth(level, u, f, opts.pre, opts.smoother);
        let mut au = vec![0.0; u.len()];
        self.apply(level, u, &mut au);
        let r: Vec<f64> = f.iter().zip(&au).map(|(a, b)| a - b).collect();
        let rc = self.restrict(level, &r, transfers);
        let mut ec = vec![0.0; rc.len()];
        let gamma = if opts.cycle == Cycle::W {
            2
        } else {
            1
        };
        for _ in 0..gamma {
            self.cycle(level + 1, &mut ec, &rc, opts, transfers);
        }
        let corr = self.prolong(level, &ec, transfers);
        for (ui, c) in u.iter_mut().zip(&corr) {
            *ui += c;
        }
        self.smooth(level, u, f, opts.post, opts.smoother);
    }

    /// One multigrid cycle for `A u = f` starting from the current `u`.
    pub fn cycle_once(
        &self,
        u: &mut [f64],
        f: &[f64],
        opts: &MgOptions,
    ) {
        self.cycle(0, u, f, opts, &mut GeometricTransferWorkspace::default());
    }

    /// Absolute residual 2-norm `||f - A u||₂` on the fine grid.
    ///
    /// Returns infinity if the norm overflows and propagates a NaN stencil
    /// residual. The norm is in the same units as the right-hand side.
    #[must_use]
    pub fn residual_norm(
        &self,
        u: &[f64],
        f: &[f64],
    ) -> f64 {
        let mut au = vec![0.0; u.len()];
        self.apply(0, u, &mut au);
        let residuals = fine_residuals(f, &au);
        let square_sum = residuals.clone().map(|r| r * r).sum::<f64>();
        stable_residual_norm(residuals, square_sum)
    }

    /// Solves `A u = f` by repeated cycles from `u = 0` (or `x0`).
    ///
    /// The history records `||f - A u||₂ / ||f||₂` for nonzero `f`.
    /// For an exactly zero right-hand side, the denominator is
    /// [`f64::MIN_POSITIVE`]. A nonzero relative residual that underflows
    /// is reported as the smallest positive subnormal. Non-finite
    /// right-hand sides or applied iterates give an undefined residual
    /// (NaN). Only a finite residual at or below the tolerance converges.
    #[must_use]
    pub fn solve(
        &self,
        f: &[f64],
        x0: Option<&[f64]>,
        opts: &MgOptions,
    ) -> MgResult {
        let u = x0.map_or_else(|| vec![0.0; f.len()], <[f64]>::to_vec);
        self.finish_cycles(u, f, opts, &mut GeometricTransferWorkspace::default())
    }

    /// Full multigrid: nested iteration from the coarsest grid upwards with
    /// `cycles` cycles per level, then (optionally) further cycles on the fine
    /// grid until the tolerance is met. One FMG sweep with one V-cycle per
    /// level already reaches discretisation accuracy.
    /// Uses the same residual and stopping criterion as [`Self::solve`].
    #[must_use]
    pub fn fmg(
        &self,
        f: &[f64],
        cycles: usize,
        opts: &MgOptions,
    ) -> MgResult {
        let mut transfers = GeometricTransferWorkspace::default();
        let mut fs: Vec<Vec<f64>> = vec![f.to_vec()];
        for l in 0..self.levels - 1 {
            let rc = self.restrict(l, &fs[l], &mut transfers);
            fs.push(rc);
        }
        let mut u = vec![0.0; fs[self.levels - 1].len()];
        self.cycle(
            self.levels - 1,
            &mut u,
            &fs[self.levels - 1],
            opts,
            &mut transfers,
        );
        for l in (0..self.levels - 1).rev() {
            let mut un = self.prolong(l, &u, &mut transfers);
            for _ in 0..cycles.max(1) {
                self.cycle(l, &mut un, &fs[l], opts, &mut transfers);
            }
            u = un;
        }
        self.finish_cycles(u, f, opts, &mut transfers)
    }

    fn finish_cycles(
        &self,
        mut u: Vec<f64>,
        f: &[f64],
        opts: &MgOptions,
        transfers: &mut GeometricTransferWorkspace,
    ) -> MgResult {
        let rhs_square_sum = f.iter().map(|v| v * v).sum::<f64>();
        let rhs_norm = rhs_square_sum.sqrt();
        let ordinary_rhs = ordinary_square_sum(rhs_square_sum);
        let zero_rhs = rhs_square_sum == 0.0 && f.iter().all(|&v| v == 0.0);
        let mut applied = vec![0.0; u.len()];
        let mut relative = |u: &[f64]| {
            self.apply(0, u, &mut applied);
            let residuals = fine_residuals(f, &applied);
            let square_sum = residuals.clone().map(|r| r * r).sum::<f64>();
            if ordinary_rhs && ordinary_square_sum(square_sum) {
                square_sum.sqrt() / rhs_norm
            } else if zero_rhs && applied.iter().all(|v| v.is_finite()) {
                stable_residual_norm(residuals, square_sum) / f64::MIN_POSITIVE
            } else {
                residual_ratio(f, &applied).unwrap_or(f64::NAN)
            }
        };
        let converged_at = |r: f64| r.is_finite() && r <= opts.tol;
        let mut history = vec![relative(&u)];
        let mut converged = converged_at(history[0]);
        while !converged && history.len() <= opts.max_cycles {
            self.cycle(0, &mut u, f, opts, transfers);
            let r = relative(&u);
            history.push(r);
            converged = converged_at(r);
        }
        MgResult {
            x: u,
            history,
            converged,
        }
    }
}

// ----------------------------------------------------------------------
// Algebraic multigrid by smoothed aggregation
// ----------------------------------------------------------------------

/// Parameters of the smoothed-aggregation AMG.
#[derive(Debug, Clone, PartialEq)]
pub struct AmgOptions {
    /// Strength-of-connection threshold `theta`: `j` is a strong neighbour
    /// of `i` when `|a_ij| >= theta sqrt(a_ii a_jj)`.
    pub theta: f64,
    /// Maximal size of the coarsest (directly solved) level.
    pub coarse_size: usize,
    /// Maximal number of levels.
    pub max_levels: usize,
    /// Smoothing sweeps (symmetric Gauss-Seidel) before and after each
    /// coarse-grid correction.
    pub sweeps: usize,
    /// Cycle shape.
    pub cycle: Cycle,
    /// Near-null-space vectors of the matrix; empty selects the constant
    /// vector. Each has the length of the matrix.
    pub near_nullspace: Vec<Vec<f64>>,
    /// Smooth the tentative prolongator (turn off for plain aggregation).
    pub smooth_prolongator: bool,
}

impl Default for AmgOptions {
    fn default() -> Self {
        Self {
            theta: 0.08,
            coarse_size: 40,
            max_levels: 25,
            sweeps: 1,
            cycle: Cycle::V,
            near_nullspace: Vec::new(),
            smooth_prolongator: true,
        }
    }
}

#[derive(Debug, Clone)]
struct AmgLevel {
    a: Csr,
    p: Csr,
    r: Csr,
    diag: Vec<f64>,
}

/// Smoothed-aggregation algebraic multigrid hierarchy.
#[derive(Debug)]
pub struct Amg {
    levels: Vec<AmgLevel>,
    coarse_a: Csr,
    coarse_lu: Option<Lu>,
    sweeps: usize,
    cycle: Cycle,
}

/// Failure of the AMG setup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmgError {
    /// The matrix is not square or has a missing or non-positive diagonal.
    BadMatrix,
}

pub(crate) fn transpose(a: &Csr) -> Csr {
    let mut cnt = vec![0usize; a.cols + 1];
    for &c in &a.indices {
        cnt[c + 1] += 1;
    }
    for i in 0..a.cols {
        cnt[i + 1] += cnt[i];
    }
    let mut pos = cnt.clone();
    let mut indices = vec![0usize; a.indices.len()];
    let mut values = vec![0.0; a.values.len()];
    for r in 0..a.rows {
        for k in a.indptr[r]..a.indptr[r + 1] {
            let c = a.indices[k];
            indices[pos[c]] = r;
            values[pos[c]] = a.values[k];
            pos[c] += 1;
        }
    }
    Csr {
        rows: a.cols,
        cols: a.rows,
        indptr: cnt,
        indices,
        values,
    }
}

pub(crate) fn spgemm(
    a: &Csr,
    b: &Csr,
) -> Csr {
    let mut indptr = vec![0usize; a.rows + 1];
    let mut indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut acc = vec![0.0; b.cols];
    let mut flag = vec![usize::MAX; b.cols];
    for i in 0..a.rows {
        let mut cols: Vec<usize> = Vec::new();
        for ka in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[ka];
            let av = a.values[ka];
            for kb in b.indptr[j]..b.indptr[j + 1] {
                let c = b.indices[kb];
                if flag[c] != i {
                    flag[c] = i;
                    acc[c] = 0.0;
                    cols.push(c);
                }
                acc[c] += av * b.values[kb];
            }
        }
        cols.sort_unstable();
        for c in cols {
            indices.push(c);
            values.push(acc[c]);
        }
        indptr[i + 1] = indices.len();
    }
    Csr {
        rows: a.rows,
        cols: b.cols,
        indptr,
        indices,
        values,
    }
}

/// Greedy aggregation on the strength graph; returns `(aggregate id per
/// node, number of aggregates)`.
pub(crate) fn aggregate(
    a: &Csr,
    theta: f64,
    diag: &[f64],
) -> (Vec<usize>, usize) {
    let n = a.rows;
    let strong = |i: usize, k: usize| -> bool {
        let j = a.indices[k];
        j != i && a.values[k].abs() >= theta * (diag[i].abs() * diag[j].abs()).sqrt()
    };
    let mut agg = vec![usize::MAX; n];
    let mut count = 0;
    // pass 1: seeds whose strong neighbourhood is entirely free
    for i in 0..n {
        if agg[i] != usize::MAX {
            continue;
        }
        let free = (a.indptr[i]..a.indptr[i + 1])
            .filter(|&k| strong(i, k))
            .all(|k| agg[a.indices[k]] == usize::MAX);
        if free {
            agg[i] = count;
            for k in a.indptr[i]..a.indptr[i + 1] {
                if strong(i, k) {
                    agg[a.indices[k]] = count;
                }
            }
            count += 1;
        }
    }
    // pass 2: attach leftovers to the strongest neighbouring aggregate
    let snapshot = agg.clone();
    for i in 0..n {
        if agg[i] != usize::MAX {
            continue;
        }
        let mut best = (0.0, usize::MAX);
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            if j != i && snapshot[j] != usize::MAX && a.values[k].abs() > best.0 {
                best = (a.values[k].abs(), snapshot[j]);
            }
        }
        if best.1 != usize::MAX {
            agg[i] = best.1;
        }
    }
    // pass 3: remaining nodes form new aggregates with their free neighbours
    for i in 0..n {
        if agg[i] != usize::MAX {
            continue;
        }
        agg[i] = count;
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            if agg[j] == usize::MAX {
                agg[j] = count;
            }
        }
        count += 1;
    }
    (agg, count)
}

impl Amg {
    pub(crate) fn heap_bytes(&self) -> usize {
        let levels = self.levels.iter().map(|level| {
            level.a.heap_bytes()
                + level.p.heap_bytes()
                + level.r.heap_bytes()
                + level.diag.capacity() * size_of::<f64>()
        });
        self.levels.capacity() * size_of::<AmgLevel>()
            + levels.sum::<usize>()
            + self.coarse_a.heap_bytes()
            + self.coarse_lu.as_ref().map_or(0, Lu::heap_bytes)
    }

    /// Builds the hierarchy for a symmetric positive definite matrix.
    ///
    /// # Errors
    /// [`AmgError::BadMatrix`] for a non-square matrix or a non-positive
    /// diagonal entry.
    pub fn build(
        a: &Csr,
        opts: &AmgOptions,
    ) -> Result<Self, AmgError> {
        if a.rows != a.cols || a.rows == 0 {
            return Err(AmgError::BadMatrix);
        }
        let mut cur = a.clone();
        let mut b: Vec<Vec<f64>> = if opts.near_nullspace.is_empty() {
            vec![vec![1.0; a.rows]]
        } else {
            opts.near_nullspace.clone()
        };
        if b.iter().any(|v| v.len() != a.rows) {
            return Err(AmgError::BadMatrix);
        }
        let mut levels: Vec<AmgLevel> = Vec::new();
        while cur.rows > opts.coarse_size && levels.len() + 1 < opts.max_levels {
            let diag = cur.diagonal();
            if diag.iter().any(|&d| d <= 0.0) {
                return Err(AmgError::BadMatrix);
            }
            let (agg, nagg) = aggregate(&cur, opts.theta, &diag);
            let nb = b.len();
            if nagg * nb >= cur.rows {
                break; // no coarsening progress
            }
            // tentative prolongator: local QR of the near-null-space block per aggregate
            let mut members: Vec<Vec<usize>> = vec![Vec::new(); nagg];
            for (i, &g) in agg.iter().enumerate() {
                members[g].push(i);
            }
            let mut trip: Vec<(usize, usize, f64)> = Vec::new();
            let mut bc_entries: Vec<(usize, usize, f64)> = Vec::new();
            let mut ncols = 0usize;
            for mem in &members {
                // modified Gram-Schmidt on the columns restricted to the aggregate
                let mut q: Vec<Vec<f64>> = Vec::new();
                let mut rmat = vec![vec![0.0; nb]; nb];
                for c in 0..nb {
                    let mut v: Vec<f64> = mem.iter().map(|&i| b[c][i]).collect();
                    for (qi, qv) in q.iter().enumerate() {
                        let dotp: f64 = v.iter().zip(qv).map(|(x, y)| x * y).sum();
                        rmat[qi][c] = dotp;
                        for (x, y) in v.iter_mut().zip(qv) {
                            *x -= dotp * y;
                        }
                    }
                    let nv = norm2(&v);
                    if nv > 1e-12 && q.len() < mem.len() {
                        rmat[q.len()][c] = nv;
                        q.push(v.iter().map(|x| x / nv).collect());
                    }
                }
                for (qi, qv) in q.iter().enumerate() {
                    for (pos, &i) in mem.iter().enumerate() {
                        trip.push((i, ncols + qi, qv[pos]));
                    }
                }
                for c in 0..nb {
                    for r in 0..q.len() {
                        bc_entries.push((c, ncols + r, rmat[r][c]));
                    }
                }
                ncols += q.len();
            }
            let mut bc: Vec<Vec<f64>> = vec![vec![0.0; ncols]; nb];
            for (c, col, v) in bc_entries {
                bc[c][col] = v;
            }
            let tent = Csr::from_triplets(cur.rows, ncols, &trip);
            let p = if opts.smooth_prolongator {
                // P = (I - omega D^{-1} A) T with omega = 4 / (3 rho), rho by Gershgorin
                let mut rho = 0.0_f64;
                for i in 0..cur.rows {
                    let s: f64 = (cur.indptr[i]..cur.indptr[i + 1])
                        .map(|k| cur.values[k].abs())
                        .sum();
                    rho = rho.max(s / diag[i]);
                }
                let omega = 4.0 / (3.0 * rho);
                let mut scaled = cur.clone();
                for i in 0..scaled.rows {
                    for k in scaled.indptr[i]..scaled.indptr[i + 1] {
                        scaled.values[k] *= -omega / diag[i];
                    }
                }
                // I + scaled
                let mut tr: Vec<(usize, usize, f64)> = Vec::new();
                for i in 0..scaled.rows {
                    tr.push((i, i, 1.0));
                    for k in scaled.indptr[i]..scaled.indptr[i + 1] {
                        tr.push((i, scaled.indices[k], scaled.values[k]));
                    }
                }
                let s = Csr::from_triplets(cur.rows, cur.rows, &tr);
                spgemm(&s, &tent)
            } else {
                tent
            };
            let r = transpose(&p);
            let ap = spgemm(&cur, &p);
            let coarse = spgemm(&r, &ap);
            levels.push(AmgLevel { a: cur, p, r, diag });
            cur = coarse;
            b = bc;
        }
        let dense_n = cur.rows;
        let mut dm = Mat::zeros(dense_n, dense_n);
        for i in 0..dense_n {
            for k in cur.indptr[i]..cur.indptr[i + 1] {
                dm.set(i, cur.indices[k], cur.values[k]);
            }
        }
        let coarse_lu = dense::lu_factor(&dm).ok();
        Ok(Self {
            levels,
            coarse_a: cur,
            coarse_lu,
            sweeps: opts.sweeps.max(1),
            cycle: opts.cycle,
        })
    }

    /// Number of levels including the coarsest.
    #[must_use]
    pub const fn num_levels(&self) -> usize {
        self.levels.len() + 1
    }

    /// Sum of the nonzeros of all level operators divided by the fine one.
    #[must_use]
    pub fn operator_complexity(&self) -> f64 {
        let fine = self
            .levels
            .first()
            .map_or(self.coarse_a.values.len(), |l| l.a.values.len());
        let total: usize = self.levels.iter().map(|l| l.a.values.len()).sum::<usize>()
            + self.coarse_a.values.len();
        total as f64 / fine as f64
    }

    fn gauss_seidel_sym(
        a: &Csr,
        x: &mut [f64],
        b: &[f64],
        diag: &[f64],
        sweeps: usize,
        forward_first: bool,
    ) {
        for _ in 0..sweeps {
            for pass in 0..2 {
                let forward = (pass == 0) == forward_first;
                let order: Box<dyn Iterator<Item = usize>> = if forward {
                    Box::new(0..a.rows)
                } else {
                    Box::new((0..a.rows).rev())
                };
                for i in order {
                    let mut s = b[i];
                    for k in a.indptr[i]..a.indptr[i + 1] {
                        let j = a.indices[k];
                        if j != i {
                            s -= a.values[k] * x[j];
                        }
                    }
                    x[i] = s / diag[i];
                }
            }
        }
    }

    fn cycle(
        &self,
        l: usize,
        x: &mut [f64],
        b: &[f64],
    ) {
        if l == self.levels.len() {
            if let Some(lu) = &self.coarse_lu
                && lu.solve_into(b, x).is_ok()
            {
                return;
            }
            let d = self.coarse_a.diagonal();
            Self::gauss_seidel_sym(&self.coarse_a, x, b, &d, 20, true);
            return;
        }
        let lev = &self.levels[l];
        Self::gauss_seidel_sym(&lev.a, x, b, &lev.diag, self.sweeps, true);
        let mut ax = vec![0.0; x.len()];
        lev.a.matvec(x, &mut ax);
        let r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let mut rc = vec![0.0; lev.r.rows];
        lev.r.matvec(&r, &mut rc);
        let mut ec = vec![0.0; rc.len()];
        let gamma = if self.cycle == Cycle::W {
            2
        } else {
            1
        };
        for _ in 0..gamma {
            self.cycle(l + 1, &mut ec, &rc);
        }
        let mut corr = vec![0.0; x.len()];
        lev.p.matvec(&ec, &mut corr);
        for (xi, c) in x.iter_mut().zip(&corr) {
            *xi += c;
        }
        Self::gauss_seidel_sym(&lev.a, x, b, &lev.diag, self.sweeps, false);
    }

    /// One cycle for `A x = b` starting from `x`.
    pub fn vcycle(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) {
        self.cycle(0, x, b);
    }

    /// Preconditioner action `z = M^{-1} r` (one cycle from zero), symmetric
    /// positive definite for symmetric `A`.
    pub fn precondition(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        z.fill(0.0);
        self.cycle(0, z, r);
    }

    /// Solves `A x = b` by repeated cycles; returns the solution, the
    /// relative residual history and whether `tol` was reached.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_cycles: usize,
    ) -> MgResult {
        let fine = self.levels.first().map_or(&self.coarse_a, |l| &l.a);
        let mut x = x0.map_or_else(|| vec![0.0; b.len()], <[f64]>::to_vec);
        let bn = norm2(b).max(f64::MIN_POSITIVE);
        let res = |x: &[f64]| {
            let mut ax = vec![0.0; b.len()];
            fine.matvec(x, &mut ax);
            norm2(&b.iter().zip(&ax).map(|(p, q)| p - q).collect::<Vec<_>>()) / bn
        };
        let mut history = vec![res(&x)];
        let mut converged = history[0] <= tol;
        while !converged && history.len() <= max_cycles {
            self.cycle(0, &mut x, b);
            let r = res(&x);
            history.push(r);
            converged = r <= tol;
        }
        MgResult {
            x,
            history,
            converged,
        }
    }
}
