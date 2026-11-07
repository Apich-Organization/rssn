//! # Sparse matrices and Krylov solvers
//!
//! A compressed-sparse-row matrix ([`Csr`]), Jacobi, ILU(0) and IC(0)
//! preconditioners, and the iterative solvers conjugate gradients (CG),
//! restarted GMRES and `BiCGSTAB`, all matrix-free: they only need a
//! matrix-vector product closure.
//!
//! A sparse direct solver complements them: [`Csr::sparse_lu`] performs a
//! left-looking Gilbert-Peierls LU factorisation with threshold pivoting
//! after a fill-reducing symmetric ordering (reverse Cuthill-McKee via
//! [`Csr::rcm_ordering`] or minimum degree via
//! [`Csr::min_degree_ordering`]), approximate minimum degree
//! ([`Csr::amd_ordering`]) or nested dissection
//! ([`Csr::nested_dissection_ordering`]), see [`crate::kernels::ordering`].
//!
//! Parallelism: [`Csr::matvec`] processes blocks of 512 rows concurrently from 60000 nonzeros on
//! (every row is an independent sequential dot product, so the product is bit-identical).
//!
//! # Limitations
//!
//! * [`Csr::ilut`] offers threshold dropping and bounded per-row fill but is
//!   still unpivoted; it can fail on an unsafe pivot and does not guarantee
//!   convergence or positive definiteness. ILU(0)
//!   ([`Csr::ilu0`]) and IC(0) ([`Csr::ic0`]) retain zero fill. There is no
//!   level-of-fill ILU(k) or incomplete pivoting, so a structurally or numerically
//!   zero pivot makes ILU(0) and IC(0) return `None`; ILUT reports an error.
//!   [`auto_preconditioner`] / [`AutoPreconditioner`] retries the chain
//!   IC(0) (only when the pattern is structurally symmetric) → ILU(0) →
//!   [`Preconditioner::Jacobi`] and reports the kind it used in
//!   [`AutoPreconditioner::kind`], but a failed factorisation is never
//!   repaired: there is still no pivoting or reordering inside a factorisation.
//! * [`Csr::ic0`] reads only the lower triangle and assumes symmetry; on a
//!   nonsymmetric or indefinite matrix it either returns `None` (non-positive
//!   pivot) or is not a valid preconditioner. The automatic helper only tries
//!   IC(0) when the *pattern* is structurally symmetric, which does not prove
//!   that the matrix is SPD. Neither solver checks that the operator is SPD,
//!   so CG can still stagnate on a matrix that is not.
//! * Preconditioners are applied as operators (one forward and one backward
//!   triangular solve); no dense inverse is formed. Those solves are
//!   sequential within a row, so they are not parallelised.
//! * ILU(0)/IC(0) reduce the iteration count but do not make the iteration
//!   count independent of the mesh size: on a 2-D Poisson problem CG still
//!   needs `O(sqrt(h^-1))` iterations, unlike a multigrid cycle.
//! * The sparse direct solver ([`Csr::sparse_lu`]) is the fallback for
//!   matrices on which the incomplete factorisation breaks down, at the cost
//!   of fill-in and `O(nnz(LU))` memory.

#[path = "krylov/ilut.rs"]
mod ilut;
pub use ilut::Ilut;
pub use ilut::IlutError;
pub use ilut::IlutOptions;
pub use ilut::IlutStats;

/// Compressed sparse row matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Row pointers (`rows + 1` entries).
    pub indptr: Vec<usize>,
    /// Column indices, sorted within each row.
    pub indices: Vec<usize>,
    /// Stored values.
    pub values: Vec<f64>,
}

/// Nonzeros from which [`Csr::matvec`] runs in parallel.
const SPMV_PAR_NNZ: usize = 60_000;
/// Rows per parallel block of [`Csr::matvec`].
const SPMV_BLOCK: usize = 512;

impl Csr {
    /// Builds a matrix from `(row, col, value)` triplets; duplicates are summed.
    #[must_use]
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        trip: &[(usize, usize, f64)],
    ) -> Self {
        let mut t: Vec<(usize, usize, f64)> = trip
            .iter()
            .copied()
            .filter(|&(r, c, _)| r < rows && c < cols)
            .collect();
        t.sort_by_key(|a| (a.0, a.1));
        let mut indptr = vec![0; rows + 1];
        let mut indices: Vec<usize> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        let mut last: Option<(usize, usize)> = None;
        for (r, c, v) in t {
            if last == Some((r, c)) {
                if let Some(x) = values.last_mut() {
                    *x += v;
                }
            } else {
                indices.push(c);
                values.push(v);
                indptr[r + 1] += 1;
                last = Some((r, c));
            }
        }
        for i in 0..rows {
            indptr[i + 1] += indptr[i];
        }
        Self {
            rows,
            cols,
            indptr,
            indices,
            values,
        }
    }

    /// Matrix-vector product `y = A x`.
    ///
    /// Every row is an independent sequential dot product, so large matrices
    /// are processed in parallel blocks of rows without changing a bit of the
    /// result.
    pub fn matvec(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) {
        let row = |i: usize| {
            accurate_dot::dot_value(
                (self.indptr[i]..self.indptr[i + 1]).map(|k| (self.values[k], x[self.indices[k]])),
            )
        };
        if self.values.len() >= SPMV_PAR_NNZ {
            let rows = self.rows.min(y.len());
            super::par::for_chunks_mut(&mut y[..rows], SPMV_BLOCK, 2, |b, block| {
                for (r, slot) in block.iter_mut().enumerate() {
                    *slot = row(b * SPMV_BLOCK + r);
                }
            });
        } else {
            for i in 0..self.rows {
                y[i] = row(i);
            }
        }
    }

    /// Diagonal of the matrix.
    #[must_use]
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.rows.min(self.cols))
            .map(|i| {
                (self.indptr[i]..self.indptr[i + 1])
                    .find(|&k| self.indices[k] == i)
                    .map_or(0.0, |k| self.values[k])
            })
            .collect()
    }

    /// Jacobi preconditioner apply function data: the inverse diagonal.
    #[must_use]
    pub fn jacobi_inverse(&self) -> Vec<f64> {
        self.diagonal()
            .iter()
            .map(|&d| if d == 0.0 { 1.0 } else { 1.0 / d })
            .collect()
    }

    /// ILU(0) factorisation (incomplete LU keeping the sparsity pattern).
    ///
    /// The factors are stored in the pattern of `A`: the strictly lower
    /// entries of `L` (unit diagonal) and the entries of `U` in the same
    /// arrays. No fill is created and no pivoting is performed.
    ///
    /// # Returns
    /// `None` when the matrix has no entry on the diagonal of some row or
    /// when a zero pivot is met (structurally or after elimination). The
    /// factorisation never returns `NaN`: on `None` the caller keeps the
    /// unpreconditioned path or falls back to [`Csr::jacobi_inverse`].
    #[must_use]
    pub fn ilu0(&self) -> Option<Ilu0> {
        if !ilut::valid_csr(self) || self.values.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let n = self.rows;
        let mut lu = self.values.clone();
        let mut diag = vec![0usize; n];
        for i in 0..n {
            diag[i] = (self.indptr[i]..self.indptr[i + 1]).find(|&k| self.indices[k] == i)?;
        }
        let mut pos = vec![usize::MAX; n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = k;
            }
            for k in self.indptr[i]..diag[i] {
                let j = self.indices[k];
                let piv = lu[diag[j]];
                if piv == 0.0 {
                    return None;
                }
                let numerator = lu[k];
                lu[k] /= piv;
                let l = lu[k];
                if !l.is_finite() || (l == 0.0 && numerator != 0.0) {
                    return None;
                }
                for kk in diag[j] + 1..self.indptr[j + 1] {
                    let p = pos[self.indices[kk]];
                    if p != usize::MAX {
                        lu[p] -= l * lu[kk];
                        if !lu[p].is_finite() {
                            return None;
                        }
                    }
                }
            }
            if !lu[diag[i]].is_finite() || lu[diag[i]] == 0.0 {
                return None;
            }
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = usize::MAX;
            }
        }
        Some(Ilu0 {
            a: Self {
                values: lu,
                ..self.clone()
            },
            diag,
        })
    }
}

/// ILU(0) factors stored in the pattern of the original matrix.
#[derive(Debug, Clone)]
pub struct Ilu0 {
    a: Csr,
    diag: Vec<usize>,
}

impl Ilu0 {
    /// Solves `L U z = r` (applies the preconditioner `M^-1`).
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let n = self.a.rows;
        for i in 0..n {
            let mut s = r[i];
            for k in self.a.indptr[i]..self.diag[i] {
                s -= self.a.values[k] * z[self.a.indices[k]];
            }
            z[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = z[i];
            for k in self.diag[i] + 1..self.a.indptr[i + 1] {
                s -= self.a.values[k] * z[self.a.indices[k]];
            }
            z[i] = s / self.a.values[self.diag[i]];
        }
    }
}

impl Csr {
    /// IC(0) factorisation `A ≈ L D L^T` (incomplete Cholesky without fill).
    ///
    /// The lower triangular factor uses the strictly lower pattern of `A`;
    /// the upper triangle of `A` is ignored, so the matrix must be symmetric
    /// and positive definite for `M = L D L^T` to be a valid CG
    /// preconditioner. On a symmetric matrix this is the same operator as
    /// [`Csr::ilu0`] written as `L D L^T` instead of `L U`.
    ///
    /// # Returns
    /// `None` when the diagonal is missing in some row or when a pivot `d_i`
    /// is not positive (zero pivot, indefinite matrix, or round-off that
    /// destroys positivity). No `NaN` is produced; the caller falls back to
    /// [`Csr::jacobi_inverse`] or to the unpreconditioned path.
    #[must_use]
    pub fn ic0(&self) -> Option<Ic0> {
        let n = self.rows;
        if self.cols != n {
            return None;
        }
        let mut diag = vec![0usize; n];
        for i in 0..n {
            diag[i] = (self.indptr[i]..self.indptr[i + 1]).find(|&k| self.indices[k] == i)?;
        }
        // `l` starts as a copy of `A`: its strictly lower entries are
        // overwritten by `L`, its upper entries are never read.
        let mut l = self.values.clone();
        let mut d = vec![0.0; n];
        let mut pos = vec![usize::MAX; n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = k;
            }
            let mut acc = 0.0;
            for k in self.indptr[i]..diag[i] {
                let j = self.indices[k];
                let mut s = self.values[k];
                for kk in self.indptr[j]..diag[j] {
                    let p = pos[self.indices[kk]];
                    if p != usize::MAX && p < diag[i] {
                        s -= l[p] * d[self.indices[kk]] * l[kk];
                    }
                }
                let dj = d[j];
                if !dj.is_finite() || dj <= 0.0 {
                    return None;
                }
                l[k] = s / dj;
                acc += l[k] * l[k] * dj;
            }
            let di = self.values[diag[i]] - acc;
            if !di.is_finite() || di <= 0.0 {
                return None;
            }
            d[i] = di;
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = usize::MAX;
            }
        }
        // Transpose of the strictly lower factor, for the backward solve.
        let mut lt = Self {
            rows: n,
            cols: n,
            indptr: vec![0; n + 1],
            indices: Vec::new(),
            values: Vec::new(),
        };
        for i in 0..n {
            for k in self.indptr[i]..diag[i] {
                lt.indptr[self.indices[k] + 1] += 1;
            }
        }
        for i in 0..n {
            lt.indptr[i + 1] += lt.indptr[i];
        }
        lt.indices = vec![0; lt.indptr[n]];
        lt.values = vec![0.0; lt.indptr[n]];
        let mut next = lt.indptr.clone();
        for i in 0..n {
            for k in self.indptr[i]..diag[i] {
                let c = self.indices[k];
                let p = next[c];
                next[c] += 1;
                lt.indices[p] = i;
                lt.values[p] = l[k];
            }
        }
        Some(Ic0 {
            l: Self {
                values: l,
                ..self.clone()
            },
            diag,
            lt,
            d,
        })
    }
}

/// IC(0) factors `M = L D L^T` stored on the pattern of `A`.
#[derive(Debug, Clone)]
pub struct Ic0 {
    l: Csr,
    diag: Vec<usize>,
    lt: Csr,
    d: Vec<f64>,
}

impl Ic0 {
    /// Applies the preconditioner `M^-1 = L^-T D^-1 L^-1` by one forward
    /// substitution, a diagonal scaling and one backward substitution.
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let n = self.l.rows;
        // L u = r (unit diagonal), then u <- D^-1 u.
        for i in 0..n {
            let mut s = r[i];
            for k in self.l.indptr[i]..self.diag[i] {
                s -= self.l.values[k] * z[self.l.indices[k]];
            }
            z[i] = s;
        }
        for (zi, di) in z.iter_mut().zip(&self.d) {
            *zi /= di;
        }
        // L^T z = u.
        for i in (0..n).rev() {
            let mut s = z[i];
            for k in self.lt.indptr[i]..self.lt.indptr[i + 1] {
                s -= self.lt.values[k] * z[self.lt.indices[k]];
            }
            z[i] = s;
        }
    }
}

/// Outcome of an iterative solve.
#[derive(Debug, Clone, PartialEq)]
pub struct IterResult {
    /// Approximate solution.
    pub x: Vec<f64>,
    /// Independently recomputed relative residual `|b - A x| / |b|`.
    /// For a zero RHS, this is zero only for an exactly zero residual;
    /// otherwise it is infinity. Invalid input or operator output also
    /// produces infinity with `converged = false`. A positive ratio that
    /// underflows is reported as the smallest positive subnormal, preserving
    /// the distinction between a tiny residual and an exactly zero one.
    pub residual: f64,
    /// Krylov iterations performed (Arnoldi steps for GMRES).
    /// Operator calibration and residual verification products are additional
    /// and are not counted.
    pub iterations: usize,
    /// Whether the tolerance was reached.
    pub converged: bool,
}

#[path = "krylov/accurate_dot.rs"]
mod accurate_dot;
#[path = "krylov/iterative.rs"]
mod iterative;

/// Correctly rounded signed significand and exponent without underflow.
pub(crate) fn exact_product_sum_scaled<I>(pairs: I) -> (f64, i32)
where
    I: IntoIterator<Item = (f64, f64)>,
    I::IntoIter: Clone,
{
    accurate_dot::product_sum(pairs)
}

/// Shared original-system residual arithmetic for validated internal adapters.
pub(crate) fn residual_ratio(
    b: &[f64],
    ax: &[f64],
) -> Option<f64> {
    iterative::original_residual_ratio(b, ax)
}

/// Identity preconditioner helper.
pub const fn identity_precond(
    r: &[f64],
    z: &mut [f64],
) {
    z.copy_from_slice(r);
}

/// Preconditioner selection for [`cg_with`], [`gmres_with`] and
/// [`bicgstab_with`].
///
/// The variant is only a selector: [`Preconditioner::apply`] runs the
/// operator (a copy for [`Preconditioner::None`], one triangular sweep per
/// level for the incomplete factors), so no dense inverse is ever formed.
#[derive(Debug, Clone, Copy)]
pub enum Preconditioner<'a> {
    /// No preconditioner: the unpreconditioned solver path, unchanged.
    None,
    /// Diagonal scaling; the slice holds the inverse diagonal
    /// ([`Csr::jacobi_inverse`]).
    Jacobi(&'a [f64]),
    /// ILU(0) factors of a general (possibly nonsymmetric) matrix
    /// ([`Csr::ilu0`]).
    Ilu0(&'a Ilu0),
    /// IC(0) factors of a symmetric positive definite matrix
    /// ([`Csr::ic0`]).
    Ic0(&'a Ic0),
}

impl Preconditioner<'_> {
    /// Writes `M^-1 r` into `z`.
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        match *self {
            | Self::None => {
                z.copy_from_slice(r);
            },
            | Self::Jacobi(inv) => {
                for ((zi, ri), di) in z.iter_mut().zip(r).zip(inv) {
                    *zi = ri * di;
                }
            },
            | Self::Ilu0(f) => f.apply(r, z),
            | Self::Ic0(f) => f.apply(r, z),
        }
    }
}

/// Which factorisation [`AutoPreconditioner`] managed to build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AutoPreconditionerKind {
    /// IC(0) factors, on a structurally symmetric pattern.
    Ic0,
    /// ILU(0) factors of the general pattern.
    Ilu0,
    /// Diagonal scaling, used when both incomplete factorisations failed.
    Jacobi,
}

/// A preconditioner chosen automatically for a matrix, owning its factors.
///
/// The selection is the chain IC(0) → ILU(0) → Jacobi: IC(0) is tried only
/// when the sparsity pattern is structurally symmetric (it is the natural
/// preconditioner for CG on an SPD matrix), then ILU(0) on the general
/// pattern, and Jacobi when every factorisation broke down. The factorisation
/// failures never reach the caller: the chosen kind is reported by
/// [`AutoPreconditioner::kind`] and the operator lives in this struct, so it
/// can be applied directly or adapted to a [`Preconditioner`] with
/// [`AutoPreconditioner::as_preconditioner`].
///
/// The choice is a heuristic, not a guarantee of quality: a structurally
/// symmetric pattern does not prove positive definiteness, so IC(0) can be
/// built for an indefinite matrix whose preconditioner is poor (though valid
/// as an operator). A successful factorisation is never rejected; only a
/// `None` from IC(0)/ILU(0) changes the choice.
#[derive(Debug, Clone)]
pub struct AutoPreconditioner {
    kind: AutoPreconditionerKind,
    ic0: Option<Ic0>,
    ilu0: Option<Ilu0>,
    jacobi: Vec<f64>,
}

impl AutoPreconditioner {
    /// Builds the strongest incomplete factorisation that succeeds; the
    /// chain always terminates with Jacobi, so this never fails.
    #[must_use]
    pub fn new(a: &Csr) -> Self {
        let jacobi = a.jacobi_inverse();
        if is_structurally_symmetric(a)
            && let Some(factors) = a.ic0()
        {
            return Self {
                kind: AutoPreconditionerKind::Ic0,
                ic0: Some(factors),
                ilu0: None,
                jacobi,
            };
        }
        if let Some(factors) = a.ilu0() {
            return Self {
                kind: AutoPreconditionerKind::Ilu0,
                ic0: None,
                ilu0: Some(factors),
                jacobi,
            };
        }
        Self {
            kind: AutoPreconditionerKind::Jacobi,
            ic0: None,
            ilu0: None,
            jacobi,
        }
    }

    /// The preconditioner that [`AutoPreconditioner::new`] selected.
    #[must_use]
    pub const fn kind(&self) -> AutoPreconditionerKind {
        self.kind
    }

    /// Writes `M^-1 r` into `z`.
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        self.as_preconditioner().apply(r, z);
    }

    /// Borrows the selected factors as a [`Preconditioner`] for
    /// [`cg_with`], [`gmres_with`] and [`bicgstab_with`].
    #[must_use]
    pub fn as_preconditioner(&self) -> Preconditioner<'_> {
        if let Some(factors) = &self.ic0 {
            return Preconditioner::Ic0(factors);
        }
        if let Some(factors) = &self.ilu0 {
            return Preconditioner::Ilu0(factors);
        }
        Preconditioner::Jacobi(&self.jacobi)
    }
}

/// Builds an [`AutoPreconditioner`] for `a`; see [`AutoPreconditioner::new`].
#[must_use]
pub fn auto_preconditioner(a: &Csr) -> AutoPreconditioner {
    AutoPreconditioner::new(a)
}

/// Whether every stored `(i, j)` has a stored `(j, i)` and the matrix is square.
fn is_structurally_symmetric(a: &Csr) -> bool {
    if a.rows != a.cols {
        return false;
    }
    for i in 0..a.rows {
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            let mirrored = (a.indptr[j]..a.indptr[j + 1]).any(|q| a.indices[q] == i);
            if !mirrored {
                return false;
            }
        }
    }
    true
}

/// Conjugate gradients with an explicit [`Preconditioner`].
///
/// [`Preconditioner::None`] dispatches to [`cg`] with [`identity_precond`],
/// so selection and direct invocation use the same scale-safe algorithm.
pub fn cg_with<A>(
    apply: A,
    precond: &Preconditioner<'_>,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
{
    match *precond {
        | Preconditioner::None => cg(apply, identity_precond, b, x0, tol, max_iter),
        | _ => cg(apply, |r, z| precond.apply(r, z), b, x0, tol, max_iter),
    }
}

/// Restarted GMRES with an explicit [`Preconditioner`]; see [`gmres`].
pub fn gmres_with<A>(
    apply: A,
    precond: &Preconditioner<'_>,
    b: &[f64],
    x0: Option<&[f64]>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
{
    match *precond {
        | Preconditioner::None => gmres(apply, identity_precond, b, x0, restart, tol, max_iter),
        | _ => {
            gmres(
                apply,
                |r, z| precond.apply(r, z),
                b,
                x0,
                restart,
                tol,
                max_iter,
            )
        },
    }
}

/// `BiCGSTAB` with an explicit [`Preconditioner`]; see [`bicgstab`].
pub fn bicgstab_with<A>(
    apply: A,
    precond: &Preconditioner<'_>,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
{
    match *precond {
        | Preconditioner::None => bicgstab(apply, identity_precond, b, x0, tol, max_iter),
        | _ => bicgstab(apply, |r, z| precond.apply(r, z), b, x0, tol, max_iter),
    }
}

/// Preconditioned conjugate gradients for SPD operators.
///
/// `apply(x, y)` computes `y = A x`; `precond(r, z)` computes `z = M^-1 r`.
/// Both operators must be linear and overwrite every output entry. Inputs,
/// output entries, and tolerance must be finite; tolerance is nonnegative and
/// an initial guess must have the RHS length. Invalid inputs or incomplete
/// operator outputs produce a nonconverged result rather than false success.
///
/// The solver normalizes the initial residual correction and search directions;
/// dot products and coefficient ratios retain binary exponents until vector
/// updates. Convergence is checked against the original system, not just the
/// recursive residual. Residual correction is refined under one shared
/// iteration budget when normalization loses tiny components. A zero iteration
/// budget still checks the initial guess.
pub fn cg<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    iterative::solve(apply, precond, b, x0, tol, max_iter, iterative::Method::Cg)
}

/// Restarted GMRES(`restart`) with right preconditioning.
///
/// Uses the same finite-input and complete-output contract as [`cg`]. The
/// preconditioned basis vectors are retained, and Arnoldi uses reorthogonalized
/// modified Gram-Schmidt. A zero restart is clamped to one, as before. Final
/// convergence uses the original system residual.
pub fn gmres<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    iterative::solve(
        apply,
        precond,
        b,
        x0,
        tol,
        max_iter,
        iterative::Method::Gmres(restart),
    )
}

/// Preconditioned `BiCGSTAB` for general nonsymmetric operators.
///
/// Uses the same finite-input, complete linear-operator, zero-RHS, and verified
/// convergence contract as [`cg`]. Coefficient ratios and preconditioned
/// directions retain their scales to avoid spurious arithmetic breakdown.
pub fn bicgstab<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    iterative::solve(
        apply,
        precond,
        b,
        x0,
        tol,
        max_iter,
        iterative::Method::Bicgstab,
    )
}

// ------------------------------------------------------------------
// Sparse direct solver: fill-reducing orderings and left-looking LU
// ------------------------------------------------------------------

/// Fill-reducing symmetric ordering used before sparse factorisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SparseOrdering {
    /// Keep the given ordering.
    Natural,
    /// Reverse Cuthill-McKee (bandwidth / profile reduction).
    Rcm,
    /// Minimum degree on the elimination graph (fill reduction).
    MinDegree,
    /// Approximate minimum degree on the quotient graph
    /// ([`crate::kernels::ordering::amd_order`]); scales to millions of vertices.
    Amd,
    /// Nested dissection by graph bisection
    /// ([`crate::kernels::ordering::nested_dissection_order`]).
    NestedDissection,
}

/// Errors of the sparse direct solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SparseLuError {
    /// The matrix is not square.
    NotSquare,
    /// A structurally or numerically zero pivot was met in this column.
    Singular(usize),
}

/// Breadth-first search from `start`; fills `dist` for the visited
/// vertices (which must be `usize::MAX` on entry) and returns them in
/// visiting order.
fn bfs_order(
    adj: &[Vec<usize>],
    start: usize,
    dist: &mut [usize],
) -> Vec<usize> {
    let mut queue = vec![start];
    dist[start] = 0;
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head];
        head += 1;
        for &w in &adj[v] {
            if dist[w] == usize::MAX {
                dist[w] = dist[v] + 1;
                queue.push(w);
            }
        }
    }
    queue
}

impl Csr {
    /// Adjacency lists of the symmetrised pattern `A + A^T` (no self loops).
    fn symmetric_adjacency(&self) -> Vec<Vec<usize>> {
        let n = self.rows.min(self.cols);
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                let j = self.indices[k];
                if j < n && j != i {
                    adj[i].push(j);
                    adj[j].push(i);
                }
            }
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }
        adj
    }

    /// Half bandwidth `max |i - j|` over the stored entries.
    #[must_use]
    pub fn bandwidth(&self) -> usize {
        let mut bw = 0;
        for i in 0..self.rows {
            for k in self.indptr[i]..self.indptr[i + 1] {
                bw = bw.max(i.abs_diff(self.indices[k]));
            }
        }
        bw
    }

    /// Symmetric permutation `B[i][j] = A[perm[i]][perm[j]]`.
    #[must_use]
    pub fn permute_symmetric(
        &self,
        perm: &[usize],
    ) -> Self {
        let n = self.rows;
        let mut inv = vec![0; n];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }
        let mut trip = Vec::with_capacity(self.values.len());
        for r in 0..n {
            for k in self.indptr[r]..self.indptr[r + 1] {
                trip.push((inv[r], inv[self.indices[k]], self.values[k]));
            }
        }
        Self::from_triplets(n, self.cols, &trip)
    }

    /// Reverse Cuthill-McKee ordering of the symmetrised pattern, with a
    /// George-Liu pseudo-peripheral start vertex per connected component.
    /// Returns `perm` with `perm[new] = old`.
    #[must_use]
    pub fn rcm_ordering(&self) -> Vec<usize> {
        let adj = self.symmetric_adjacency();
        let n = adj.len();
        let deg: Vec<usize> = adj.iter().map(Vec::len).collect();
        let mut visited = vec![false; n];
        let mut order: Vec<usize> = Vec::with_capacity(n);
        let mut dist = vec![usize::MAX; n];
        for seed in 0..n {
            if visited[seed] {
                continue;
            }
            // pseudo-peripheral start vertex of the component
            let mut start = seed;
            let mut queue = bfs_order(&adj, start, &mut dist);
            let mut ecc = dist[*queue.last().unwrap_or(&start)];
            for _ in 0..16 {
                let cand = queue
                    .iter()
                    .copied()
                    .filter(|&v| dist[v] == ecc)
                    .min_by_key(|&v| deg[v])
                    .unwrap_or(start);
                for &v in &queue {
                    dist[v] = usize::MAX;
                }
                let q2 = bfs_order(&adj, cand, &mut dist);
                let e2 = dist[*q2.last().unwrap_or(&cand)];
                if e2 > ecc {
                    start = cand;
                    ecc = e2;
                    for &v in &q2 {
                        dist[v] = usize::MAX;
                    }
                    queue = bfs_order(&adj, start, &mut dist);
                } else {
                    for &v in &q2 {
                        dist[v] = usize::MAX;
                    }
                    queue = bfs_order(&adj, start, &mut dist);
                    break;
                }
            }
            for &v in &queue {
                dist[v] = usize::MAX;
            }
            // Cuthill-McKee numbering from `start`
            let begin = order.len();
            order.push(start);
            visited[start] = true;
            let mut head = begin;
            while head < order.len() {
                let v = order[head];
                head += 1;
                let mut nb: Vec<usize> = adj[v].iter().copied().filter(|&w| !visited[w]).collect();
                nb.sort_by_key(|&w| deg[w]);
                for w in nb {
                    visited[w] = true;
                    order.push(w);
                }
            }
        }
        order.reverse();
        order
    }

    /// Minimum-degree ordering on the elimination graph of the
    /// symmetrised pattern: repeatedly eliminates a vertex of smallest
    /// current degree, turning its neighbourhood into a clique (exact
    /// degrees, lazily maintained heap). Returns `perm[new] = old`.
    #[must_use]
    pub fn min_degree_ordering(&self) -> Vec<usize> {
        use std::cmp::Reverse;
        use std::collections::BTreeSet;
        use std::collections::BinaryHeap;
        let adj0 = self.symmetric_adjacency();
        let n = adj0.len();
        let mut adj: Vec<BTreeSet<usize>> =
            adj0.into_iter().map(|v| v.into_iter().collect()).collect();
        let mut heap: BinaryHeap<Reverse<(usize, usize)>> =
            (0..n).map(|v| Reverse((adj[v].len(), v))).collect();
        let mut done = vec![false; n];
        let mut order = Vec::with_capacity(n);
        while let Some(Reverse((d, v))) = heap.pop() {
            if done[v] || d != adj[v].len() {
                continue;
            }
            done[v] = true;
            order.push(v);
            let nb: Vec<usize> = adj[v].iter().copied().collect();
            for &u in &nb {
                adj[u].remove(&v);
            }
            for (a, &u) in nb.iter().enumerate() {
                for &w in &nb[a + 1..] {
                    if adj[u].insert(w) {
                        adj[w].insert(u);
                    }
                }
            }
            for &u in &nb {
                heap.push(Reverse((adj[u].len(), u)));
            }
            adj[v].clear();
        }
        order
    }

    /// Approximate minimum degree ordering (quotient graph with element
    /// absorption, supervariables and mass elimination); see
    /// [`crate::kernels::ordering::amd_order`]. Returns `perm[new] = old`.
    #[must_use]
    pub fn amd_ordering(&self) -> Vec<usize> {
        crate::kernels::ordering::amd_order(&self.symmetric_adjacency())
    }

    /// Nested dissection ordering; see
    /// [`crate::kernels::ordering::nested_dissection_order`].
    #[must_use]
    pub fn nested_dissection_ordering(&self) -> Vec<usize> {
        crate::kernels::ordering::nested_dissection_order(&self.symmetric_adjacency())
    }

    /// Sparse LU factorisation with a fill-reducing symmetric ordering
    /// and threshold partial pivoting (left-looking Gilbert-Peierls with
    /// a depth-first symbolic reach); see [`SparseLu`].
    ///
    /// # Errors
    /// [`SparseLuError::NotSquare`] or [`SparseLuError::Singular`].
    pub fn sparse_lu(
        &self,
        ordering: SparseOrdering,
    ) -> Result<SparseLu, SparseLuError> {
        self.sparse_lu_with(ordering, 0.1)
    }

    /// As [`Csr::sparse_lu`] with an explicit pivot threshold in `(0, 1]`:
    /// the diagonal entry is kept as pivot when it is at least
    /// `threshold` times the largest candidate (`1.0` is classical
    /// partial pivoting, small values favour sparsity).
    ///
    /// # Errors
    /// [`SparseLuError::NotSquare`] or [`SparseLuError::Singular`].
    pub fn sparse_lu_with(
        &self,
        ordering: SparseOrdering,
        threshold: f64,
    ) -> Result<SparseLu, SparseLuError> {
        if self.rows != self.cols {
            return Err(SparseLuError::NotSquare);
        }
        let n = self.rows;
        let perm: Vec<usize> = match ordering {
            | SparseOrdering::Natural => (0..n).collect(),
            | SparseOrdering::Rcm => self.rcm_ordering(),
            | SparseOrdering::MinDegree => self.min_degree_ordering(),
            | SparseOrdering::Amd => self.amd_ordering(),
            | SparseOrdering::NestedDissection => self.nested_dissection_ordering(),
        };
        let b = self.permute_symmetric(&perm);
        // CSC of B
        let mut cp = vec![0usize; n + 1];
        for &c in &b.indices {
            cp[c + 1] += 1;
        }
        for j in 0..n {
            cp[j + 1] += cp[j];
        }
        let mut ci = vec![0usize; b.indices.len()];
        let mut cx = vec![0.0; b.indices.len()];
        let mut next = cp.clone();
        for r in 0..n {
            for k in b.indptr[r]..b.indptr[r + 1] {
                let c = b.indices[k];
                ci[next[c]] = r;
                cx[next[c]] = b.values[k];
                next[c] += 1;
            }
        }
        let tol = threshold.clamp(1e-12, 1.0);
        let mut lp = vec![0usize; n + 1];
        let mut li: Vec<usize> = Vec::new();
        let mut lx: Vec<f64> = Vec::new();
        let mut up = vec![0usize; n + 1];
        let mut ui: Vec<usize> = Vec::new();
        let mut ux: Vec<f64> = Vec::new();
        let mut pinv = vec![usize::MAX; n];
        let mut x = vec![0.0; n];
        let mut xi = vec![0usize; n];
        let mut pstack = vec![0usize; n];
        let mut mark = vec![usize::MAX; n];
        for k in 0..n {
            lp[k] = li.len();
            up[k] = ui.len();
            // symbolic reach of column k of B in the graph of L
            let mut top = n;
            for p in cp[k]..cp[k + 1] {
                let start = ci[p];
                if mark[start] == k {
                    continue;
                }
                // iterative depth-first search; the stack lives in xi[0..=head]
                let mut head = 0usize;
                xi[0] = start;
                loop {
                    let j = xi[head];
                    let jnew = pinv[j];
                    if mark[j] != k {
                        mark[j] = k;
                        pstack[head] = if jnew == usize::MAX {
                            0
                        } else {
                            lp[jnew]
                        };
                    }
                    let p2 = if jnew == usize::MAX {
                        0
                    } else {
                        lp[jnew + 1]
                    };
                    let mut q = pstack[head];
                    let mut pushed = false;
                    while q < p2 {
                        let i = li[q];
                        q += 1;
                        if mark[i] == k {
                            continue;
                        }
                        pstack[head] = q;
                        head += 1;
                        xi[head] = i;
                        pushed = true;
                        break;
                    }
                    if !pushed {
                        top -= 1;
                        xi[top] = j;
                        if head == 0 {
                            break;
                        }
                        head -= 1;
                    }
                }
            }
            for p in top..n {
                x[xi[p]] = 0.0;
            }
            for p in cp[k]..cp[k + 1] {
                x[ci[p]] = cx[p];
            }
            // sparse triangular solve with the unit lower factor
            for p in top..n {
                let j = xi[p];
                let jn = pinv[j];
                if jn == usize::MAX {
                    continue;
                }
                let xj = x[j];
                for q in lp[jn] + 1..lp[jn + 1] {
                    x[li[q]] -= lx[q] * xj;
                }
            }
            // pivot search
            let mut ipiv = usize::MAX;
            let mut amax = -1.0;
            for p in top..n {
                let i = xi[p];
                if pinv[i] == usize::MAX {
                    let t = x[i].abs();
                    if t > amax {
                        amax = t;
                        ipiv = i;
                    }
                } else {
                    ui.push(pinv[i]);
                    ux.push(x[i]);
                }
            }
            if ipiv == usize::MAX || amax <= 0.0 || !amax.is_finite() {
                return Err(SparseLuError::Singular(k));
            }
            if pinv[k] == usize::MAX && x[k].abs() >= amax * tol {
                ipiv = k;
            }
            let pivot = x[ipiv];
            ui.push(k);
            ux.push(pivot);
            pinv[ipiv] = k;
            li.push(ipiv);
            lx.push(1.0);
            for p in top..n {
                let i = xi[p];
                if pinv[i] == usize::MAX {
                    li.push(i);
                    lx.push(x[i] / pivot);
                }
                x[i] = 0.0;
            }
        }
        lp[n] = li.len();
        up[n] = ui.len();
        for v in &mut li {
            *v = pinv[*v];
        }
        Ok(SparseLu {
            n,
            perm,
            pinv,
            lp,
            li,
            lx,
            up,
            ui,
            ux,
        })
    }
}

/// Sparse LU factors `P B = L U` of a symmetrically permuted matrix
/// `B = A[perm][:, perm]` (CSC storage of `L` with unit diagonal and `U`).
#[derive(Debug, Clone)]
pub struct SparseLu {
    n: usize,
    perm: Vec<usize>,
    pinv: Vec<usize>,
    lp: Vec<usize>,
    li: Vec<usize>,
    lx: Vec<f64>,
    up: Vec<usize>,
    ui: Vec<usize>,
    ux: Vec<f64>,
}

impl SparseLu {
    /// Number of stored entries in `L` and `U` (a measure of fill-in).
    #[must_use]
    pub const fn nnz(&self) -> usize {
        self.lx.len() + self.ux.len()
    }

    /// Dimension of the factored matrix.
    #[must_use]
    pub const fn dim(&self) -> usize {
        self.n
    }

    /// The symmetric fill-reducing permutation used (`perm[new] = old`).
    #[must_use]
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// Solves `A x = b`.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
    ) -> Vec<f64> {
        let n = self.n;
        let mut w = vec![0.0; n];
        for i in 0..n {
            w[self.pinv[i]] = b[self.perm[i]];
        }
        // forward: L w = w (unit diagonal stored first in each column)
        for j in 0..n {
            let wj = w[j];
            if wj != 0.0 {
                for q in self.lp[j] + 1..self.lp[j + 1] {
                    w[self.li[q]] -= self.lx[q] * wj;
                }
            }
        }
        // backward: U z = w (diagonal stored last in each column)
        for j in (0..n).rev() {
            let last = self.up[j + 1] - 1;
            w[j] /= self.ux[last];
            let wj = w[j];
            if wj != 0.0 {
                for q in self.up[j]..last {
                    w[self.ui[q]] -= self.ux[q] * wj;
                }
            }
        }
        let mut x = vec![0.0; n];
        for i in 0..n {
            x[self.perm[i]] = w[i];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 5-point Laplacian on an `m x m` grid with homogeneous Dirichlet data:
    /// `4` on the diagonal and `-1` on the four grid neighbours.
    fn poisson2d(m: usize) -> Csr {
        let mut t = Vec::with_capacity(5 * m * m);
        for i in 0..m {
            for j in 0..m {
                let r = i * m + j;
                t.push((r, r, 4.0));
                if i > 0 {
                    t.push((r, r - m, -1.0));
                }
                if i + 1 < m {
                    t.push((r, r + m, -1.0));
                }
                if j > 0 {
                    t.push((r, r - 1, -1.0));
                }
                if j + 1 < m {
                    t.push((r, r + 1, -1.0));
                }
            }
        }
        Csr::from_triplets(m * m, m * m, &t)
    }

    fn rhs(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| ((i * 7919 % 101) as f64) / 50.0 - 1.0)
            .collect()
    }

    fn dense(a: &Csr) -> Vec<Vec<f64>> {
        let n = a.rows;
        let mut out = vec![vec![0.0; n]; n];
        for j in 0..n {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let mut col = vec![0.0; n];
            a.matvec(&e, &mut col);
            for i in 0..n {
                out[i][j] = col[i];
            }
        }
        out
    }

    fn matmul(
        l: &[Vec<f64>],
        u: &[Vec<f64>],
    ) -> Vec<Vec<f64>> {
        let n = l.len();
        let mut p = vec![vec![0.0; n]; n];
        for i in 0..n {
            for k in 0..n {
                let lik = l[i][k];
                if lik != 0.0 {
                    for j in 0..n {
                        p[i][j] += lik * u[k][j];
                    }
                }
            }
        }
        p
    }

    /// Exact eigenvalues of the 5-point Laplacian on an `m x m` grid, used as
    /// the condition-number proxy `lambda_max / lambda_min`.
    fn poisson_condition(m: usize) -> f64 {
        let theta = std::f64::consts::PI / (m as f64 + 1.0);
        let lam_min = 4.0 - 4.0 * theta.cos();
        let lam_max = 4.0 + 4.0 * theta.cos();
        lam_max / lam_min
    }

    fn rel_residual(
        a: &Csr,
        x: &[f64],
        b: &[f64],
    ) -> f64 {
        let mut ax = vec![0.0; b.len()];
        a.matvec(x, &mut ax);
        let num: f64 = ax
            .iter()
            .zip(b)
            .map(|(p, q)| (p - q) * (p - q))
            .sum::<f64>()
            .sqrt();
        let den: f64 = b.iter().map(|q| q * q).sum::<f64>().sqrt();
        num / den.max(f64::MIN_POSITIVE)
    }

    /// Textbook unpreconditioned CG, kept independent of [`cg`] so the
    /// identity-preconditioner path can be compared bit for bit.
    fn reference_cg(
        a: &Csr,
        b: &[f64],
        tol: f64,
        max_iter: usize,
    ) -> (Vec<f64>, usize) {
        let n = b.len();
        let mut x = vec![0.0; n];
        let mut ax = vec![0.0; n];
        a.matvec(&x, &mut ax);
        let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
        let mut p = r.clone();
        let mut rs = r.iter().map(|v| v * v).sum::<f64>();
        let bn = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);
        let mut it = 0;
        while it < max_iter && rs.sqrt() / bn > tol {
            let mut ap = vec![0.0; n];
            a.matvec(&p, &mut ap);
            let alpha = rs / p.iter().zip(&ap).map(|(u, v)| u * v).sum::<f64>();
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let rs_new = r.iter().map(|v| v * v).sum::<f64>();
            let beta = rs_new / rs;
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }
            rs = rs_new;
            it += 1;
        }
        (x, it)
    }

    /// `L U` of an ILU(0) factor, as a dense matrix.
    fn ilu_dense(f: &Ilu0) -> Vec<Vec<f64>> {
        let n = f.a.rows;
        let mut l = vec![vec![0.0; n]; n];
        let mut u = vec![vec![0.0; n]; n];
        for i in 0..n {
            l[i][i] = 1.0;
            for k in f.a.indptr[i]..f.a.indptr[i + 1] {
                let j = f.a.indices[k];
                if j < i {
                    l[i][j] = f.a.values[k];
                } else {
                    u[i][j] = f.a.values[k];
                }
            }
        }
        matmul(&l, &u)
    }

    /// `L D L^T` of an IC(0) factor, as a dense matrix.
    fn ic0_dense(f: &Ic0) -> Vec<Vec<f64>> {
        let n = f.l.rows;
        let mut l = vec![vec![0.0; n]; n];
        for i in 0..n {
            l[i][i] = 1.0;
            for k in f.l.indptr[i]..f.diag[i] {
                l[i][f.l.indices[k]] = f.l.values[k];
            }
        }
        let mut p = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                p[i][j] = (0..n).map(|k| l[i][k] * f.d[k] * l[j][k]).sum();
            }
        }
        p
    }

    /// On-pattern error and off-pattern fill of `P = M` against `A`.
    fn pattern_error(
        a: &[Vec<f64>],
        p: &[Vec<f64>],
    ) -> (f64, f64, f64) {
        let n = a.len();
        let mut on = 0.0_f64;
        let mut off = 0.0_f64;
        let mut off_fro = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                if a[i][j] == 0.0 {
                    off = off.max(p[i][j].abs());
                    off_fro += p[i][j] * p[i][j];
                } else {
                    on = on.max((p[i][j] - a[i][j]).abs());
                }
            }
        }
        (on, off, off_fro.sqrt())
    }

    #[test]
    fn sparse_poisson_cg_ic0_cuts_iterations() {
        for &m in &[64_usize, 128] {
            let a = poisson2d(m);
            let n = a.rows;
            let b = rhs(n);
            let dinv = a.jacobi_inverse();
            let ic = a.ic0().unwrap();
            let ilu = a.ilu0().unwrap();
            let plain = cg(
                |x, y| a.matvec(x, y),
                identity_precond,
                &b,
                None,
                1e-10,
                5000,
            );
            let jac = cg_with(
                |x, y| a.matvec(x, y),
                &Preconditioner::Jacobi(&dinv),
                &b,
                None,
                1e-10,
                5000,
            );
            let pre = cg_with(
                |x, y| a.matvec(x, y),
                &Preconditioner::Ic0(&ic),
                &b,
                None,
                1e-10,
                5000,
            );
            let pilu = cg_with(
                |x, y| a.matvec(x, y),
                &Preconditioner::Ilu0(&ilu),
                &b,
                None,
                1e-10,
                5000,
            );
            println!(
                "poisson m={m} n={n} kappa~{:.0}: none={} jacobi={} ic0={} ilu0={}",
                poisson_condition(m),
                plain.iterations,
                jac.iterations,
                pre.iterations,
                pilu.iterations
            );
            assert!(plain.converged && jac.converged && pre.converged && pilu.converged);
            assert!(
                pre.iterations * 2 < plain.iterations,
                "IC(0) did not cut the count: {} vs {}",
                pre.iterations,
                plain.iterations
            );
            assert!(pilu.iterations * 2 < plain.iterations);
            assert!(rel_residual(&a, &pre.x, &b) < 1e-8);
            assert!(rel_residual(&a, &pilu.x, &b) < 1e-8);
        }
    }

    #[test]
    fn sparse_ilu0_factor_matches_a_on_pattern() {
        let a = poisson2d(8);
        let f = a.ilu0().unwrap();
        let (on, off, off_fro) = pattern_error(&dense(&a), &ilu_dense(&f));
        println!(
            "ILU(0) m=8: max |LU-A| on pattern={on:.3e}, max off-pattern |LU|={off:.3e}, \
             off-pattern Frobenius={off_fro:.3e}"
        );
        assert!(on < 1e-12, "on-pattern error {on}");
        assert!(off > 1e-6, "no fill was dropped: {off}");
        assert!(off_fro > 1e-6);
        // The apply really is (LU)^-1: compare with a dense solve.
        let r = rhs(64);
        let mut z = vec![0.0; 64];
        f.apply(&r, &mut z);
        let mat = crate::kernels::dense::Mat::from_rows(&ilu_dense(&f)).unwrap();
        let want = crate::kernels::dense::solve(&mat, &r).unwrap();
        let d = z
            .iter()
            .zip(&want)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(d < 1e-12, "apply is not (LU)^-1: {d}");
    }

    #[test]
    fn sparse_ic0_factor_matches_a_on_pattern() {
        let a = poisson2d(8);
        let f = a.ic0().unwrap();
        let (on, off, off_fro) = pattern_error(&dense(&a), &ic0_dense(&f));
        println!(
            "IC(0) m=8: max |LDL^T-A| on pattern={on:.3e}, max off-pattern |LDL^T|={off:.3e}, \
             off-pattern Frobenius={off_fro:.3e}"
        );
        assert!(on < 1e-12, "on-pattern error {on}");
        assert!(off > 1e-6, "no fill was dropped: {off}");
        assert!(off_fro > 1e-6);
        // The apply really is (L D L^T)^-1: compare with a dense solve.
        let r = rhs(64);
        let mut z = vec![0.0; 64];
        f.apply(&r, &mut z);
        let mat = crate::kernels::dense::Mat::from_rows(&ic0_dense(&f)).unwrap();
        let want = crate::kernels::dense::cholesky_solve(&mat, &r).unwrap();
        let d = z
            .iter()
            .zip(&want)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(d < 1e-12, "apply is not (LDL^T)^-1: {d}");
    }

    #[test]
    fn sparse_ilu0_breakdown_is_reported() {
        // Zero pivot after the first elimination step.
        let singular =
            Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]);
        assert!(singular.ilu0().is_none());
        // Structurally missing diagonal.
        let off_diag = Csr::from_triplets(2, 2, &[(0, 1, 1.0), (1, 0, 1.0)]);
        assert!(off_diag.ilu0().is_none());
        assert!(off_diag.ic0().is_none());
        // Indefinite: d_1 = 1 - 2^2 / 1 < 0.
        let indefinite =
            Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0)]);
        assert!(indefinite.ic0().is_none());
        // A valid factor never produces NaN on apply.
        let f = poisson2d(8).ic0().unwrap();
        let r = vec![1.0; 64];
        let mut z = vec![0.0; 64];
        f.apply(&r, &mut z);
        assert!(z.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn sparse_preconditioned_cg_matches_direct_solve() {
        let a = poisson2d(16);
        let n = a.rows;
        let b = rhs(n);
        let ic = a.ic0().unwrap();
        let got = cg_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::Ic0(&ic),
            &b,
            None,
            1e-14,
            2000,
        );
        assert!(got.converged);
        let mat = crate::kernels::dense::Mat::from_rows(&dense(&a)).unwrap();
        let direct = crate::kernels::dense::cholesky_solve(&mat, &b).unwrap();
        let err = got
            .x
            .iter()
            .zip(&direct)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        println!("max |x_ic0-cg - x_direct| = {err:.3e}");
        assert!(err < 1e-10, "solution mismatch {err}");
    }

    #[test]
    fn sparse_unpreconditioned_selector_matches_scale_safe_cg() {
        let a = poisson2d(32);
        let n = a.rows;
        let b = rhs(n);
        let (xr, it) = reference_cg(&a, &b, 1e-10, 2000);
        let got = cg(
            |x, y| a.matvec(x, y),
            identity_precond,
            &b,
            None,
            1e-10,
            2000,
        );
        assert_eq!(got.iterations, it);
        for i in 0..n {
            assert!(
                (got.x[i] - xr[i]).abs() <= 1e-10 * xr[i].abs().max(1.0),
                "row {i}: {} vs {}",
                got.x[i],
                xr[i]
            );
        }
        let sel = cg_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::None,
            &b,
            None,
            1e-10,
            2000,
        );
        assert_eq!(sel, got);
    }

    #[test]
    fn sparse_gmres_bicgstab_preconditioner_selection() {
        let mut t = Vec::new();
        for i in 0..400 {
            t.push((i, i, 3.0));
            if i > 0 {
                t.push((i, i - 1, -2.0));
            }
            if i + 1 < 400 {
                t.push((i, i + 1, -0.5));
            }
        }
        let a = Csr::from_triplets(400, 400, &t);
        let b = vec![1.0; 400];
        let ilu = a.ilu0().unwrap();
        let g0 = gmres(
            |x, y| a.matvec(x, y),
            identity_precond,
            &b,
            None,
            30,
            1e-10,
            2000,
        );
        let g1 = gmres_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::None,
            &b,
            None,
            30,
            1e-10,
            2000,
        );
        assert_eq!(g0, g1);
        let gp = gmres_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::Ilu0(&ilu),
            &b,
            None,
            10,
            1e-12,
            100,
        );
        assert!(gp.converged && gp.iterations < g0.iterations);
        let b0 = bicgstab(
            |x, y| a.matvec(x, y),
            identity_precond,
            &b,
            None,
            1e-10,
            1000,
        );
        let b1 = bicgstab_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::None,
            &b,
            None,
            1e-10,
            1000,
        );
        assert_eq!(b0, b1);
        let bp = bicgstab_with(
            |x, y| a.matvec(x, y),
            &Preconditioner::Ilu0(&ilu),
            &b,
            None,
            1e-12,
            100,
        );
        assert!(bp.converged && bp.iterations < b0.iterations);
        println!(
            "gmres none={} ilu0={}, bicgstab none={} ilu0={}",
            g0.iterations, gp.iterations, b0.iterations, bp.iterations
        );
    }
}
