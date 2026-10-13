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
#![allow(
    clippy::missing_const_for_fn,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::too_many_arguments,
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::similar_names,
    clippy::float_cmp,
    clippy::suboptimal_flops,
    clippy::type_complexity
)]

use crate::kernels::dense::{self, Lu, Mat};
use crate::kernels::krylov::Csr;

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

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

impl GeometricMg {
    /// Creates the hierarchy for `n = 2^k - 1` interior points per direction.
    ///
    /// Returns `None` when `dim` is not 1, 2 or 3, when `n` is not of the
    /// form `2^k - 1` (`k >= 1`) or when `shift < 0`.
    #[must_use]
    pub fn new(dim: usize, n: usize, shift: f64) -> Option<Self> {
        if !(1..=3).contains(&dim) || n == 0 || (n + 1) & n != 0 || shift < 0.0 {
            return None;
        }
        let levels = (n + 1).trailing_zeros() as usize;
        Some(Self { dim, n, shift, levels })
    }

    /// Number of unknowns on the fine grid.
    #[must_use]
    pub fn size(&self) -> usize {
        self.n.pow(self.dim as u32)
    }

    /// Number of levels (the coarsest has a single unknown).
    #[must_use]
    pub fn levels(&self) -> usize {
        self.levels
    }

    fn n_at(&self, level: usize) -> usize {
        ((self.n + 1) >> level) - 1
    }

    fn h_at(&self, level: usize) -> f64 {
        1.0 / ((self.n_at(level) + 1) as f64)
    }

    /// Applies the discrete operator on `level` (0 is the finest grid).
    pub fn apply(&self, level: usize, u: &[f64], out: &mut [f64]) {
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

    fn smooth(&self, level: usize, u: &mut [f64], f: &[f64], sweeps: usize, sm: Smoother) {
        let n = self.n_at(level);
        let ih2 = 1.0 / (self.h_at(level) * self.h_at(level));
        let diag = 2.0 * self.dim as f64 * ih2 + self.shift;
        let strides: Vec<usize> = (0..self.dim).map(|d| n.pow(d as u32)).collect();
        match sm {
            Smoother::Jacobi => {
                let omega = match self.dim {
                    1 => 2.0 / 3.0,
                    2 => 4.0 / 5.0,
                    _ => 6.0 / 7.0,
                };
                let mut au = vec![0.0; u.len()];
                for _ in 0..sweeps {
                    self.apply(level, u, &mut au);
                    for i in 0..u.len() {
                        u[i] += omega * (f[i] - au[i]) / diag;
                    }
                }
            }
            Smoother::RedBlackGaussSeidel => {
                for _ in 0..sweeps {
                    for color in 0..2 {
                        for idx in 0..u.len() {
                            let mut parity = 0;
                            let mut s = f[idx];
                            for d in 0..self.dim {
                                let c = (idx / strides[d]) % n;
                                parity += c;
                                if c > 0 {
                                    s += ih2 * u[idx - strides[d]];
                                }
                                if c + 1 < n {
                                    s += ih2 * u[idx + strides[d]];
                                }
                            }
                            if parity % 2 == color {
                                u[idx] = s / diag;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Full-weighting restriction from `level` to `level + 1`.
    fn restrict(&self, level: usize, r: &[f64]) -> Vec<f64> {
        let mut cur = r.to_vec();
        let mut dims = vec![self.n_at(level); self.dim];
        let nc = self.n_at(level + 1);
        for axis in 0..self.dim {
            let nf = dims[axis];
            let inner: usize = dims[..axis].iter().product();
            let outer: usize = dims[axis + 1..].iter().product();
            let mut next = vec![0.0; inner * nc * outer];
            for o in 0..outer {
                for j in 0..nc {
                    for i in 0..inner {
                        let c = 2 * j + 1;
                        let mut s = 0.5 * cur[i + inner * (c + nf * o)];
                        if c > 0 {
                            s += 0.25 * cur[i + inner * (c - 1 + nf * o)];
                        }
                        if c + 1 < nf {
                            s += 0.25 * cur[i + inner * (c + 1 + nf * o)];
                        }
                        next[i + inner * (j + nc * o)] = s;
                    }
                }
            }
            cur = next;
            dims[axis] = nc;
        }
        cur
    }

    /// Multilinear interpolation from `level + 1` to `level`.
    fn prolong(&self, level: usize, e: &[f64]) -> Vec<f64> {
        let mut cur = e.to_vec();
        let nc = self.n_at(level + 1);
        let nf = self.n_at(level);
        let mut dims = vec![nc; self.dim];
        for axis in 0..self.dim {
            let inner: usize = dims[..axis].iter().product();
            let outer: usize = dims[axis + 1..].iter().product();
            let mut next = vec![0.0; inner * nf * outer];
            for o in 0..outer {
                for i in 0..nf {
                    for k in 0..inner {
                        let v = if i % 2 == 1 {
                            cur[k + inner * (i / 2 + nc * o)]
                        } else {
                            let left = if i >= 2 { cur[k + inner * (i / 2 - 1 + nc * o)] } else { 0.0 };
                            let right = if i / 2 < nc { cur[k + inner * (i / 2 + nc * o)] } else { 0.0 };
                            0.5 * (left + right)
                        };
                        next[k + inner * (i + nf * o)] = v;
                    }
                }
            }
            cur = next;
            dims[axis] = nf;
        }
        cur
    }

    fn cycle(&self, level: usize, u: &mut [f64], f: &[f64], opts: &MgOptions) {
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
        let rc = self.restrict(level, &r);
        let mut ec = vec![0.0; rc.len()];
        let gamma = if opts.cycle == Cycle::W { 2 } else { 1 };
        for _ in 0..gamma {
            self.cycle(level + 1, &mut ec, &rc, opts);
        }
        let corr = self.prolong(level, &ec);
        for (ui, c) in u.iter_mut().zip(&corr) {
            *ui += c;
        }
        self.smooth(level, u, f, opts.post, opts.smoother);
    }

    /// One multigrid cycle for `A u = f` starting from the current `u`.
    pub fn cycle_once(&self, u: &mut [f64], f: &[f64], opts: &MgOptions) {
        self.cycle(0, u, f, opts);
    }

    /// Residual 2-norm `|f - A u|` on the fine grid.
    #[must_use]
    pub fn residual_norm(&self, u: &[f64], f: &[f64]) -> f64 {
        let mut au = vec![0.0; u.len()];
        self.apply(0, u, &mut au);
        norm2(&f.iter().zip(&au).map(|(a, b)| a - b).collect::<Vec<_>>())
    }

    /// Solves `A u = f` by repeated cycles from `u = 0` (or `x0`).
    #[must_use]
    pub fn solve(&self, f: &[f64], x0: Option<&[f64]>, opts: &MgOptions) -> MgResult {
        let mut u = x0.map_or_else(|| vec![0.0; f.len()], <[f64]>::to_vec);
        let bn = norm2(f).max(f64::MIN_POSITIVE);
        let mut history = vec![self.residual_norm(&u, f) / bn];
        let mut converged = history[0] <= opts.tol;
        while !converged && history.len() <= opts.max_cycles {
            self.cycle(0, &mut u, f, opts);
            let r = self.residual_norm(&u, f) / bn;
            history.push(r);
            converged = r <= opts.tol;
        }
        MgResult { x: u, history, converged }
    }

    /// Full multigrid: nested iteration from the coarsest grid upwards with
    /// `cycles` cycles per level, then (optionally) further cycles on the fine
    /// grid until the tolerance is met. One FMG sweep with one V-cycle per
    /// level already reaches discretisation accuracy.
    #[must_use]
    pub fn fmg(&self, f: &[f64], cycles: usize, opts: &MgOptions) -> MgResult {
        let mut fs: Vec<Vec<f64>> = vec![f.to_vec()];
        for l in 0..self.levels - 1 {
            let rc = self.restrict(l, &fs[l]);
            fs.push(rc);
        }
        let mut u = vec![0.0; fs[self.levels - 1].len()];
        self.cycle(self.levels - 1, &mut u, &fs[self.levels - 1], opts);
        for l in (0..self.levels - 1).rev() {
            let mut un = self.prolong(l, &u);
            for _ in 0..cycles.max(1) {
                self.cycle(l, &mut un, &fs[l], opts);
            }
            u = un;
        }
        let bn = norm2(f).max(f64::MIN_POSITIVE);
        let mut history = vec![self.residual_norm(&u, f) / bn];
        let mut converged = history[0] <= opts.tol;
        while !converged && history.len() <= opts.max_cycles {
            self.cycle(0, &mut u, f, opts);
            let r = self.residual_norm(&u, f) / bn;
            history.push(r);
            converged = r <= opts.tol;
        }
        MgResult { x: u, history, converged }
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

fn transpose(a: &Csr) -> Csr {
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
    Csr { rows: a.cols, cols: a.rows, indptr: cnt, indices, values }
}

fn spgemm(a: &Csr, b: &Csr) -> Csr {
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
    Csr { rows: a.rows, cols: b.cols, indptr, indices, values }
}

/// Greedy aggregation on the strength graph; returns `(aggregate id per
/// node, number of aggregates)`.
fn aggregate(a: &Csr, theta: f64, diag: &[f64]) -> (Vec<usize>, usize) {
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
    /// Builds the hierarchy for a symmetric positive definite matrix.
    ///
    /// # Errors
    /// [`AmgError::BadMatrix`] for a non-square matrix or a non-positive
    /// diagonal entry.
    pub fn build(a: &Csr, opts: &AmgOptions) -> Result<Self, AmgError> {
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
                    let s: f64 = (cur.indptr[i]..cur.indptr[i + 1]).map(|k| cur.values[k].abs()).sum();
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
        Ok(Self { levels, coarse_a: cur, coarse_lu, sweeps: opts.sweeps.max(1), cycle: opts.cycle })
    }

    /// Number of levels including the coarsest.
    #[must_use]
    pub fn num_levels(&self) -> usize {
        self.levels.len() + 1
    }

    /// Sum of the nonzeros of all level operators divided by the fine one.
    #[must_use]
    pub fn operator_complexity(&self) -> f64 {
        let fine = self.levels.first().map_or(self.coarse_a.values.len(), |l| l.a.values.len());
        let total: usize =
            self.levels.iter().map(|l| l.a.values.len()).sum::<usize>() + self.coarse_a.values.len();
        total as f64 / fine as f64
    }

    fn gauss_seidel_sym(a: &Csr, x: &mut [f64], b: &[f64], diag: &[f64], sweeps: usize, forward_first: bool) {
        for _ in 0..sweeps {
            for pass in 0..2 {
                let forward = (pass == 0) == forward_first;
                let order: Box<dyn Iterator<Item = usize>> =
                    if forward { Box::new(0..a.rows) } else { Box::new((0..a.rows).rev()) };
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

    fn cycle(&self, l: usize, x: &mut [f64], b: &[f64]) {
        if l == self.levels.len() {
            if let Some(lu) = &self.coarse_lu {
                x.copy_from_slice(&lu.solve(b));
            } else {
                let d = self.coarse_a.diagonal();
                Self::gauss_seidel_sym(&self.coarse_a, x, b, &d, 20, true);
            }
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
        let gamma = if self.cycle == Cycle::W { 2 } else { 1 };
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
    pub fn vcycle(&self, b: &[f64], x: &mut [f64]) {
        self.cycle(0, x, b);
    }

    /// Preconditioner action `z = M^{-1} r` (one cycle from zero), symmetric
    /// positive definite for symmetric `A`.
    pub fn precondition(&self, r: &[f64], z: &mut [f64]) {
        z.iter_mut().for_each(|v| *v = 0.0);
        self.cycle(0, z, r);
    }

    /// Solves `A x = b` by repeated cycles; returns the solution, the
    /// relative residual history and whether `tol` was reached.
    #[must_use]
    pub fn solve(&self, b: &[f64], x0: Option<&[f64]>, tol: f64, max_cycles: usize) -> MgResult {
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
        MgResult { x, history, converged }
    }
}
