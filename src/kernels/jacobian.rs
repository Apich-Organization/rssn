//! # Structured Jacobians and linear solvers for stiff integrators
//!
//! Infrastructure shared by [`crate::kernels::ode_stiff`]:
//!
//! * [`JacStructure`]: dense, banded (`kl` sub- and `ku` super-diagonals) or
//!   general sparse Jacobian pattern;
//! * [`color_columns`] / [`fd_jacobian_colored`]: greedy column colouring of
//!   a sparsity pattern (Curtis, Powell and Reid) and the resulting
//!   finite-difference Jacobian that needs one right-hand-side evaluation per
//!   colour instead of one per column (`kl + ku + 1` for a band);
//! * [`BandLu`]: banded LU factorisation with partial pivoting
//!   (LAPACK `gbtrf`/`gbtrs` layout);
//! * [`Factor`]: a dense, banded or sparse ([`Csr::sparse_lu`]) factorisation
//!   behind one `solve` interface.

use crate::kernels::dense::{self, Lu, Mat};
use crate::kernels::krylov::{Csr, SparseLu, SparseOrdering};
use crate::kernels::ode_stiff::JacobianFn;

/// User-supplied sparse Jacobian `J(t, y)` of the right-hand side.
pub type SparseJacobianFn<'a> = dyn Fn(f64, &[f64]) -> Csr + 'a;

/// Structure of the Jacobian `df/dy`.
#[derive(Debug, Clone, Copy, Default)]
pub enum JacStructure<'a> {
    /// Full matrix, dense LU.
    #[default]
    Dense,
    /// Entries `J[i][j]` with `-kl <= j - i <= ku` only; banded LU.
    Banded {
        /// Number of sub-diagonals.
        kl: usize,
        /// Number of super-diagonals.
        ku: usize,
    },
    /// Entries within the stored pattern of the given matrix (values are
    /// ignored); sparse LU with approximate-minimum-degree ordering.
    Sparse(&'a Csr),
}

/// Linear-algebra back end chosen from a [`JacStructure`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverKind {
    /// Dense LU.
    Dense,
    /// Banded LU.
    Banded,
    /// Sparse LU.
    Sparse,
}

/// A Jacobian in the representation matching its [`JacStructure`].
#[derive(Debug, Clone)]
pub enum Jac {
    /// Dense matrix.
    Dense(Mat),
    /// Compressed sparse rows (banded and sparse structures).
    Sparse(Csr),
}

impl Jac {
    /// Calls `g(i, j, value)` for every stored entry.
    pub fn for_each(&self, mut g: impl FnMut(usize, usize, f64)) {
        match self {
            Self::Dense(m) => {
                for i in 0..m.rows {
                    for j in 0..m.cols {
                        g(i, j, m.at(i, j));
                    }
                }
            }
            Self::Sparse(c) => {
                for i in 0..c.rows {
                    for k in c.indptr[i]..c.indptr[i + 1] {
                        g(i, c.indices[k], c.values[k]);
                    }
                }
            }
        }
    }
}

/// Mass matrix of `M y' = f`, dense or sparse.
#[derive(Debug, Clone, Copy)]
pub enum Mass<'a> {
    /// Dense mass matrix.
    Dense(&'a Mat),
    /// Sparse mass matrix.
    Sparse(&'a Csr),
}

impl Mass<'_> {
    /// Matrix-vector product `M x`.
    #[must_use]
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        match self {
            Mass::Dense(m) => m.matvec(x),
            Mass::Sparse(c) => {
                let mut y = vec![0.0; c.rows];
                c.matvec(x, &mut y);
                y
            }
        }
    }

    /// Half bandwidths `(kl, ku)` of the stored nonzeros.
    #[must_use]
    pub fn band(&self) -> (usize, usize) {
        let (mut kl, mut ku) = (0, 0);
        self.for_each(|i, j, _| {
            kl = kl.max(i.saturating_sub(j));
            ku = ku.max(j.saturating_sub(i));
        });
        (kl, ku)
    }

    /// Dimension (rows) of the matrix.
    #[must_use]
    pub const fn rows(&self) -> usize {
        match self {
            Mass::Dense(m) => m.rows,
            Mass::Sparse(c) => c.rows,
        }
    }

    /// Dimension (columns) of the matrix.
    #[must_use]
    pub const fn cols(&self) -> usize {
        match self {
            Mass::Dense(m) => m.cols,
            Mass::Sparse(c) => c.cols,
        }
    }

    /// Calls `g(i, j, value)` for every nonzero entry.
    pub fn for_each(&self, mut g: impl FnMut(usize, usize, f64)) {
        match self {
            Mass::Dense(m) => {
                for i in 0..m.rows {
                    for j in 0..m.cols {
                        let v = m.at(i, j);
                        if v != 0.0 {
                            g(i, j, v);
                        }
                    }
                }
            }
            Mass::Sparse(c) => {
                for i in 0..c.rows {
                    for k in c.indptr[i]..c.indptr[i + 1] {
                        g(i, c.indices[k], c.values[k]);
                    }
                }
            }
        }
    }
}

/// Pattern of the band `-kl <= j - i <= ku` of an `n x n` matrix.
#[must_use]
pub fn band_pattern(n: usize, kl: usize, ku: usize) -> Csr {
    let mut trip = Vec::new();
    for i in 0..n {
        let lo = i.saturating_sub(kl);
        let hi = (i + ku).min(n.saturating_sub(1));
        for j in lo..=hi {
            trip.push((i, j, 1.0));
        }
    }
    Csr::from_triplets(n, n, &trip)
}

/// Greedy colouring of the columns of a sparsity pattern: columns that share
/// a row get different colours. Returns the colour of each column and the
/// number of colours.
#[must_use]
pub fn color_columns(pattern: &Csr) -> (Vec<usize>, usize) {
    let nc = pattern.cols;
    // column -> rows
    let mut col_rows: Vec<Vec<usize>> = vec![Vec::new(); nc];
    for i in 0..pattern.rows {
        for k in pattern.indptr[i]..pattern.indptr[i + 1] {
            col_rows[pattern.indices[k]].push(i);
        }
    }
    let mut order: Vec<usize> = (0..nc).collect();
    order.sort_by_key(|&j| std::cmp::Reverse(col_rows[j].len()));
    let mut color = vec![usize::MAX; nc];
    let mut ncolors = 0;
    let mut forbidden = vec![usize::MAX; nc + 1];
    for &j in &order {
        for &i in &col_rows[j] {
            for k in pattern.indptr[i]..pattern.indptr[i + 1] {
                let c = color[pattern.indices[k]];
                if c != usize::MAX {
                    forbidden[c] = j;
                }
            }
        }
        let mut c = 0;
        while forbidden[c] == j {
            c += 1;
        }
        color[j] = c;
        ncolors = ncolors.max(c + 1);
    }
    (color, ncolors)
}

/// Finite-difference Jacobian restricted to `pattern`, computed with one
/// evaluation of `f` per colour of `colors` (from [`color_columns`]).
/// `f0 = f(t, y)`; `nfev` counts the evaluations.
#[must_use]
pub fn fd_jacobian_colored<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t: f64,
    y: &[f64],
    f0: &[f64],
    pattern: &Csr,
    colors: &(Vec<usize>, usize),
    nfev: &mut usize,
) -> Csr {
    let n = y.len();
    let (color, ncolors) = colors;
    let delta: Vec<f64> = y.iter().map(|v| f64::EPSILON.sqrt() * v.abs().max(1e-5)).collect();
    let mut cols_of: Vec<Vec<usize>> = vec![Vec::new(); *ncolors];
    for (j, &c) in color.iter().enumerate() {
        cols_of[c].push(j);
    }
    let mut out = pattern.clone();
    let mut yp = y.to_vec();
    let mut fp = vec![0.0; n];
    // pattern is row-major; fill entry (i, j) from the evaluation of colour(j)
    let mut diffs: Vec<Vec<f64>> = Vec::with_capacity(*ncolors);
    for cols in &cols_of {
        for &j in cols {
            yp[j] = y[j] + delta[j];
        }
        f(t, &yp, &mut fp);
        *nfev += 1;
        for &j in cols {
            yp[j] = y[j];
        }
        diffs.push((0..n).map(|i| fp[i] - f0[i]).collect());
    }
    for i in 0..pattern.rows {
        for k in pattern.indptr[i]..pattern.indptr[i + 1] {
            let j = pattern.indices[k];
            out.values[k] = diffs[color[j]][i] / delta[j];
        }
    }
    out
}

/// Banded LU factorisation with partial pivoting of an `n x n` matrix with
/// `kl` sub- and `ku` super-diagonals.
#[derive(Debug, Clone)]
pub struct BandLu {
    n: usize,
    kl: usize,
    ku: usize,
    /// Row `i` stores column `j` at `a[i][j + kl - i]`, `j - i in [-kl, ku + kl]`.
    a: Vec<Vec<f64>>,
    piv: Vec<usize>,
}

impl BandLu {
    /// Factors the matrix given by triplets `(row, col, value)` (entries
    /// outside the band are an error). Returns `None` when singular or when
    /// an entry lies outside the band.
    #[must_use]
    pub fn factor(n: usize, kl: usize, ku: usize, trip: &[(usize, usize, f64)]) -> Option<Self> {
        let w = 2 * kl + ku + 1;
        let mut a = vec![vec![0.0; w]; n];
        for &(i, j, v) in trip {
            if i >= n || j >= n || j + kl < i || j > i + ku {
                return None;
            }
            a[i][j + kl - i] += v;
        }
        let mut piv = vec![0usize; n];
        for k in 0..n {
            let last_row = (k + kl).min(n - 1);
            let last_col = (k + ku + kl).min(n - 1);
            let mut p = k;
            let mut best = a[k][kl].abs();
            for i in k + 1..=last_row {
                let v = a[i][k + kl - i].abs();
                if v > best {
                    best = v;
                    p = i;
                }
            }
            if best == 0.0 || !best.is_finite() {
                return None;
            }
            piv[k] = p;
            if p != k {
                for j in k..=last_col {
                    let (x, y) = (a[k][j + kl - k], a[p][j + kl - p]);
                    a[k][j + kl - k] = y;
                    a[p][j + kl - p] = x;
                }
            }
            let d = a[k][kl];
            for i in k + 1..=last_row {
                let m = a[i][k + kl - i] / d;
                a[i][k + kl - i] = m;
                if m != 0.0 {
                    for j in k + 1..=last_col {
                        let v = a[k][j + kl - k];
                        a[i][j + kl - i] -= m * v;
                    }
                }
            }
        }
        Some(Self { n, kl, ku, a, piv })
    }

    /// Solves `A x = b`.
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let (n, kl, ku) = (self.n, self.kl, self.ku);
        let mut x = b.to_vec();
        for k in 0..n {
            let p = self.piv[k];
            if p != k {
                x.swap(k, p);
            }
            let last_row = (k + kl).min(n - 1);
            let xk = x[k];
            for i in k + 1..=last_row {
                x[i] -= self.a[i][k + kl - i] * xk;
            }
        }
        for k in (0..n).rev() {
            let last_col = (k + ku + kl).min(n - 1);
            let mut s = x[k];
            for j in k + 1..=last_col {
                s -= self.a[k][j + kl - k] * x[j];
            }
            x[k] = s / self.a[k][kl];
        }
        x
    }
}

/// Factorisation of a (block) system with a dense, banded or sparse solver.
#[derive(Debug, Clone)]
pub enum Factor {
    /// Dense LU.
    Dense(Lu),
    /// Banded LU. When `stages > 1` the unknowns are stored stage-major
    /// (`stage * n + var`) by the caller and permuted to variable-major
    /// (`var * stages + stage`) here, which makes the stage matrix banded.
    Band {
        /// The factorisation of the permuted matrix.
        lu: BandLu,
        /// Variables per stage block.
        n: usize,
        /// Number of stage blocks.
        stages: usize,
    },
    /// Sparse LU.
    Sparse(SparseLu),
}

impl Factor {
    /// Solves the factored system.
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        match self {
            Self::Dense(l) => l.solve(b),
            Self::Sparse(l) => l.solve(b),
            Self::Band { lu, n, stages } => {
                if *stages <= 1 {
                    return lu.solve(b);
                }
                let (n, s) = (*n, *stages);
                let mut pb = vec![0.0; n * s];
                for st in 0..s {
                    for p in 0..n {
                        pb[p * s + st] = b[st * n + p];
                    }
                }
                let px = lu.solve(&pb);
                let mut x = vec![0.0; n * s];
                for st in 0..s {
                    for p in 0..n {
                        x[st * n + p] = px[p * s + st];
                    }
                }
                x
            }
        }
    }

    /// Factors the system with the given triplets in stage-major numbering
    /// (`stages` blocks of size `n`) according to `structure`; `band` is the
    /// `(kl, ku)` of one block (used for banded structure only).
    #[must_use]
    pub fn from_triplets(
        n: usize,
        stages: usize,
        trip: &[(usize, usize, f64)],
        kind: SolverKind,
        band: (usize, usize),
    ) -> Option<Self> {
        match kind {
            SolverKind::Dense => {
                let mut m = Mat::zeros(n * stages, n * stages);
                for &(i, j, v) in trip {
                    let x = m.at(i, j);
                    m.set(i, j, x + v);
                }
                dense::lu_factor(&m).ok().map(Self::Dense)
            }
            SolverKind::Banded => {
                let (kl, ku) = band;
                let s = stages;
                let (kl2, ku2) = if s > 1 {
                    (kl * s + s - 1, ku * s + s - 1)
                } else {
                    (kl, ku)
                };
                let mapped: Vec<(usize, usize, f64)> = if s > 1 {
                    trip.iter()
                        .map(|&(i, j, v)| {
                            let (bi, p) = (i / n, i % n);
                            let (bj, q) = (j / n, j % n);
                            (p * s + bi, q * s + bj, v)
                        })
                        .collect()
                } else {
                    trip.to_vec()
                };
                let lu = BandLu::factor(n * s, kl2, ku2, &mapped)?;
                Some(Self::Band { lu, n, stages })
            }
            SolverKind::Sparse => {
                let c = Csr::from_triplets(n * stages, n * stages, trip);
                c.sparse_lu(SparseOrdering::Amd).ok().map(Self::Sparse)
            }
        }
    }
}

/// Source of Jacobians `J(t, y)`: a user callback (dense or sparse) or
/// finite differences that exploit the declared structure.
pub struct JacobianProvider<'a> {
    dense: Option<&'a JacobianFn<'a>>,
    sparse: Option<&'a SparseJacobianFn<'a>>,
    /// Effective solver kind (a sparse callback forces a sparse solver).
    pub kind: SolverKind,
    pattern: Option<Csr>,
    colors: Option<(Vec<usize>, usize)>,
    /// Effective `(kl, ku)` of the band (banded structure only).
    pub band: (usize, usize),
}

impl<'a> JacobianProvider<'a> {
    /// Prepares the provider for a system of dimension `n`. A sparse
    /// callback with `Dense` structure makes the structure sparse.
    /// Returns `None` if a banded/sparse pattern does not fit `n`.
    #[must_use]
    pub fn new(
        n: usize,
        structure: JacStructure<'a>,
        dense: Option<&'a JacobianFn<'a>>,
        sparse: Option<&'a SparseJacobianFn<'a>>,
    ) -> Option<Self> {
        let (kind, pattern, band) = match structure {
            JacStructure::Dense if sparse.is_some() => (SolverKind::Sparse, None, (0, 0)),
            JacStructure::Dense => (SolverKind::Dense, None, (0, 0)),
            JacStructure::Banded { kl, ku } => {
                (SolverKind::Banded, Some(band_pattern(n, kl, ku)), (kl, ku))
            }
            JacStructure::Sparse(p) => {
                if p.rows != n || p.cols != n {
                    return None;
                }
                (SolverKind::Sparse, Some(p.clone()), (0, 0))
            }
        };
        let colors = pattern.as_ref().map(color_columns);
        Some(Self { dense, sparse, kind, pattern, colors, band })
    }

    /// Evaluates the Jacobian at `(t, y)` with `f0 = f(t, y)`.
    pub fn eval<F: Fn(f64, &[f64], &mut [f64])>(
        &self,
        f: &F,
        t: f64,
        y: &[f64],
        f0: &[f64],
        nfev: &mut usize,
    ) -> Jac {
        if let Some(sj) = self.sparse {
            return Jac::Sparse(sj(t, y));
        }
        if let (Some(pat), Some(col)) = (&self.pattern, &self.colors) {
            self.dense.map_or_else(
                || Jac::Sparse(fd_jacobian_colored(f, t, y, f0, pat, col, nfev)),
                |dj| {
                    let m = dj(t, y);
                    let mut c = pat.clone();
                    for i in 0..c.rows {
                        for k in c.indptr[i]..c.indptr[i + 1] {
                            c.values[k] = m.at(i, c.indices[k]);
                        }
                    }
                    Jac::Sparse(c)
                },
            )
        } else if let Some(dj) = self.dense {
            Jac::Dense(dj(t, y))
        } else {
            *nfev += y.len();
            Jac::Dense(fd_jacobian_dense(f, t, y, f0))
        }
    }
}

/// Dense finite-difference Jacobian (one evaluation per column).
#[must_use]
pub fn fd_jacobian_dense<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t: f64,
    y: &[f64],
    f0: &[f64],
) -> Mat {
    let n = y.len();
    let mut jac = Mat::zeros(n, n);
    let mut yp = y.to_vec();
    let mut fp = vec![0.0; n];
    for j in 0..n {
        let dj = f64::EPSILON.sqrt() * y[j].abs().max(1e-5);
        yp[j] = y[j] + dj;
        f(t, &yp, &mut fp);
        yp[j] = y[j];
        for i in 0..n {
            jac.set(i, j, (fp[i] - f0[i]) / dj);
        }
    }
    jac
}

/// Assembles the triplets of the block system
/// `sum over blocks (bi, bj) of  cm(bi, bj) M + cj(bi, bj) J` in stage-major
/// numbering for `stages` blocks of size `n`.
#[must_use]
pub fn block_triplets(
    n: usize,
    stages: usize,
    cm: &dyn Fn(usize, usize) -> f64,
    cj: &dyn Fn(usize, usize) -> f64,
    mass: Option<Mass<'_>>,
    jac: &Jac,
) -> Vec<(usize, usize, f64)> {
    let mut trip = Vec::new();
    for bi in 0..stages {
        for bj in 0..stages {
            let a = cj(bi, bj);
            if a != 0.0 {
                jac.for_each(|p, q, v| trip.push((bi * n + p, bj * n + q, a * v)));
            }
            let m = cm(bi, bj);
            if m != 0.0 {
                match mass {
                    Some(mm) => mm.for_each(|p, q, v| trip.push((bi * n + p, bj * n + q, m * v))),
                    None => {
                        for p in 0..n {
                            trip.push((bi * n + p, bj * n + p, m));
                        }
                    }
                }
            }
        }
    }
    trip
}
