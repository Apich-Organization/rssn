//! Dense linear algebra: LU, Cholesky, QR, SVD, symmetric and general eigenvalues.
//!
//! # Dense linear algebra on row-major `f64` matrices
//!
//! A small self-contained dense toolkit: LU with partial pivoting,
//! Cholesky, Householder QR, least squares, the Jacobi eigenvalue algorithm
//! for symmetric matrices, a thin bidiagonal SVD (singular values,
//! pseudo-inverse, rank, 2-norm condition number) and the eigenvalues of a
//! general real matrix by Hessenberg reduction and the shifted QR
//! algorithm. [`eigen_generalized_symmetric`] solves `A x = lambda B x` (Cholesky reduction for
//! positive definite `B`; `B^{-1} A` with possibly complex eigenvalues for indefinite `B`).
//!
//! # Limitations
//!
//! * [`eigen_symmetric`] is the symmetric path: it returns real eigenvalues and
//!   orthonormal real eigenvectors, and it must not be fed a non-symmetric
//!   matrix (it does not symmetrise the input). Use [`eigen`] for a general matrix.
//! * [`eigen`] returns complex eigenpairs, so a real matrix with a
//!   complex-conjugate pair gets genuinely complex vectors. It computes a real
//!   Schur form, which is quasi-upper-triangular: a conjugate pair stays as a
//!   2x2 diagonal block (or as two 1x1 blocks when the pair is real).
//! * A defective matrix (for example a Jordan block) has fewer independent
//!   eigenvectors than eigenvalues. [`eigen`] still returns one vector per
//!   eigenvalue, but the extra vectors are parallel to an earlier one and the
//!   result is flagged with [`Eigen::defective`]; the returned vectors are then
//!   not a basis. The residual of a genuine eigenvector stays small even when
//!   the matrix is defective, so the flag is the reliable signal.
//! * [`eigen`] is a dense method with cubic work per QR sweep; for large sparse
//!   operators use [`crate::kernels::eigs`] / [`crate::kernels::eigs_ext`].
//! * [`eigen_generalized_symmetric`] returns no eigenvectors when `B` is
//!   indefinite (the reduction is to `B^{-1} A`, which is not symmetric).
//!
//! # Parallelism and determinism
//!
//! The parallel kernels keep the arithmetic order of the sequential reference: only
//! the independent outer loop (a row, a column, a row block) is handed to
//! `kernels::par`, and each entry is still accumulated over the same index in
//! the same direction. A reduction never changes association, so the results
//! are bit-identical for one thread, four threads or the sequential path.
//!
//! [`Mat::matmul`] updates blocks of output rows so one panel of the right
//! operand is streamed once per block instead of once per row; [`lu_factor`]
//! splits the rows of its trailing update; [`cholesky`] splits the rows of a
//! column; [`qr`] splits the rows of the `R` and `Q` updates; and the
//! Hessenberg reduction and [`Lu::inverse`] split rows and columns. Work below
//! `DENSE_PAR_WORK` element updates stays on the sequential path, so small
//! matrices never pay a fork/join.

use std::fmt;

use num_complex::Complex64;

/// Errors of the dense routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenseError {
    /// The embedding host cancelled the computation.
    Cancelled,
    /// Shapes of the operands do not match.
    Shape,
    /// The matrix is singular (or not positive definite for Cholesky).
    Singular,
    /// An iteration did not converge.
    NoConvergence,
    /// The input holds NaN or an infinity.
    NonFinite,
}

impl fmt::Display for DenseError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Cancelled => write!(f, "computation cancelled"),
            | Self::Shape => write!(f, "incompatible shapes"),
            | Self::Singular => write!(f, "matrix is singular or not positive definite"),
            | Self::NoConvergence => write!(f, "iteration did not converge"),
            | Self::NonFinite => write!(f, "input contains NaN or infinity"),
        }
    }
}

impl std::error::Error for DenseError {}
impl From<crate::cancellation::Cancelled> for DenseError {
    fn from(_: crate::cancellation::Cancelled) -> Self {
        Self::Cancelled
    }
}

/// Multiply-adds from which dense products run in parallel.
const DENSE_PAR_FLOPS: usize = 4_000_000;

/// Element updates from which a single LU trailing update is split across
/// threads.
///
/// Below it the fork/join costs more than the step; the value keeps every
/// 200 x 200 factorization on the sequential path.
const DENSE_PAR_WORK: usize = 250_000;

/// Element updates from which one Cholesky column is split across threads.
///
/// The columns of a factorisation are a serial chain, so it pays to split a
/// column earlier than a rank-one update. The largest column of a 256 x 256
/// matrix is `128 * 128` and stays sequential.
const DENSE_PAR_WORK_COLUMN: usize = 20_000;

/// Element updates from which one Householder step's `R` and `Q` updates are
/// split across threads.
const DENSE_PAR_WORK_QR: usize = 50_000;

/// Rows of a product that are updated together, so one panel of the right
/// operand is streamed once per block instead of once per output row.
const DENSE_MATMUL_BLOCK_F64: usize = 16_384;

/// Row-major dense matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Mat {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Entries, `data[i * cols + j]`.
    pub data: Vec<f64>,
}

impl Mat {
    /// Zero matrix.
    #[must_use]
    pub fn zeros(
        rows: usize,
        cols: usize,
    ) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Identity matrix.
    #[must_use]
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.data[i * n + i] = 1.0;
        }
        m
    }

    /// Builds a matrix from row slices.
    ///
    /// # Errors
    /// Returns [`DenseError::Shape`] when rows have different lengths.
    pub fn from_rows(rows: &[Vec<f64>]) -> Result<Self, DenseError> {
        let c = rows.first().map_or(0, Vec::len);
        if rows.iter().any(|r| r.len() != c) {
            return Err(DenseError::Shape);
        }
        Ok(Self {
            rows: rows.len(),
            cols: c,
            data: rows.concat(),
        })
    }

    fn has_valid_shape(&self) -> bool {
        self.rows.checked_mul(self.cols) == Some(self.data.len())
    }

    /// Whether every entry is finite (no NaN, no infinity).
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.data.iter().all(|x| x.is_finite())
    }

    /// Entry `(i, j)`.
    #[must_use]
    pub fn at(
        &self,
        i: usize,
        j: usize,
    ) -> f64 {
        self.data[i * self.cols + j]
    }

    /// Sets entry `(i, j)`.
    pub fn set(
        &mut self,
        i: usize,
        j: usize,
        v: f64,
    ) {
        self.data[i * self.cols + j] = v;
    }

    /// Transpose (cache-blocked in 32x32 tiles and parallelised for large matrices).
    #[must_use]
    pub fn transpose(&self) -> Self {
        const TILE: usize = 32;
        let rows = self.rows;
        let cols = self.cols;
        let mut t = Self::zeros(cols, rows);
        if rows == 0 || cols == 0 {
            return t;
        }
        if rows.saturating_mul(cols) < 4_096 {
            for i in 0..rows {
                let src_row = i * cols;
                for j in 0..cols {
                    t.data[j * rows + i] = self.data[src_row + j];
                }
            }
            return t;
        }
        // Output matrix `t` has shape (cols, rows). Partition its rows (j in 0..cols)
        // in blocks of TILE so each parallel chunk owns disjoint rows of `t`.
        let grain = if rows.saturating_mul(cols) >= DENSE_PAR_WORK {
            2
        } else {
            usize::MAX
        };
        let chunk_elems = TILE * rows;
        super::par::for_chunks_mut(&mut t.data, chunk_elems, grain, |jb, dst_band| {
            let j0 = jb * TILE;
            let j_len = dst_band.len() / rows;
            let mut i0 = 0;
            while i0 < rows {
                let i_end = (i0 + TILE).min(rows);
                for j_off in 0..j_len {
                    let j = j0 + j_off;
                    let dst_row = j_off * rows;
                    for i in i0..i_end {
                        dst_band[dst_row + i] = self.data[i * cols + j];
                    }
                }
                i0 += TILE;
            }
        });
        t
    }

    /// Matrix-vector product.
    #[must_use]
    pub fn matvec(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        let row = |i| self.row_product(i, x);
        if self.rows * self.cols >= DENSE_PAR_FLOPS / 8 {
            super::par::map_indexed(self.rows, 2, row)
        } else {
            (0..self.rows).map(row).collect()
        }
    }

    /// Matrix-vector product into an existing buffer with one entry per row.
    pub(crate) fn matvec_into(
        &self,
        x: &[f64],
        output: &mut [f64],
    ) {
        assert_eq!(output.len(), self.rows);
        let grain = if self.rows * self.cols >= DENSE_PAR_FLOPS / 8 {
            2
        } else {
            usize::MAX
        };
        super::par::for_each_mut(output, grain, |i, value| *value = self.row_product(i, x));
    }

    fn row_product(
        &self,
        row: usize,
        x: &[f64],
    ) -> f64 {
        self.data[row * self.cols..(row + 1) * self.cols]
            .iter()
            .zip(x)
            .map(|(a, b)| a * b)
            .sum()
    }

    /// Matrix product.
    ///
    /// # Errors
    /// Returns [`DenseError::Shape`] for incompatible matrix dimensions or storage.
    pub fn matmul(
        &self,
        other: &Self,
    ) -> Result<Self, DenseError> {
        if self.cols != other.rows || !self.has_valid_shape() || !other.has_valid_shape() {
            return Err(DenseError::Shape);
        }
        let size = self.rows.checked_mul(other.cols).ok_or(DenseError::Shape)?;
        let mut data = Vec::new();
        data.try_reserve_exact(size)
            .map_err(|_| DenseError::Shape)?;
        data.resize(size, 0.0);
        let mut out = Self {
            rows: self.rows,
            cols: other.cols,
            data,
        };
        let (m, oc) = (self.cols, other.cols);
        if oc == 0 {
            return Ok(out);
        }
        // A block of output rows shares one pass over `other`, so a panel of the
        // right operand is streamed once per block instead of once per row.
        // Inside a block `k` stays the outer loop and every entry accumulates
        // in ascending `k`, which is exactly the serial order.
        let block_rows = (DENSE_MATMUL_BLOCK_F64 / oc).clamp(1, 64);
        let grain = if size.saturating_mul(m) >= DENSE_PAR_FLOPS {
            2
        } else {
            usize::MAX
        };
        super::par::for_chunks_mut(
            &mut out.data,
            block_rows.saturating_mul(oc),
            grain,
            |block, rows| {
                let row0 = block.saturating_mul(block_rows);
                for k in 0..m {
                    let right = &other.data[k * oc..(k + 1) * oc];
                    for (local, out_row) in rows.chunks_exact_mut(oc).enumerate() {
                        let a = self.data[(row0 + local) * m + k];
                        if a != 0.0 {
                            for (slot, &factor) in out_row.iter_mut().zip(right) {
                                *slot += a * factor;
                            }
                        }
                    }
                }
            },
        );
        Ok(out)
    }

    /// Frobenius norm with a serial SIMD reduction and scaling at extreme magnitudes.
    #[must_use]
    pub fn norm_fro(&self) -> f64 {
        super::vector::norm_serial(&self.data)
    }
}

/// LU factorisation `P A = L U` with partial pivoting.
#[derive(Debug, Clone)]
pub struct Lu {
    factors: Mat,
    piv: Vec<usize>,
    sign: f64,
}

/// Computes the LU factorisation of a square matrix.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input or incompatible matrix storage,
/// [`DenseError::Singular`] when a pivot is exactly zero, or
/// [`DenseError::NonFinite`] for NaN or infinite entries.
pub fn lu_factor(a: &Mat) -> Result<Lu, DenseError> {
    crate::cancellation::checkpoint()?;
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    let mut lu = a.clone();
    let mut piv: Vec<usize> = (0..n).collect();
    let mut sign = 1.0;
    for k in 0..n {
        crate::cancellation::checkpoint()?;
        let mut p = k;
        for i in k + 1..n {
            if lu.at(i, k).abs() > lu.at(p, k).abs() {
                p = i;
            }
        }
        if lu.at(p, k) == 0.0 {
            return Err(DenseError::Singular);
        }
        if p != k {
            for j in 0..n {
                lu.data.swap(k * n + j, p * n + j);
            }
            piv.swap(k, p);
            sign = -sign;
        }
        let d = lu.data[k * n + k];
        let trailing = n - k - 1;
        if trailing == 0 {
            break;
        }
        // Split off the pivot row so the rows below can be updated concurrently
        // while it is read. Each trailing entry is updated exactly once per `k`
        // and the row update runs left to right, so the bits are the serial
        // ones; only the order in which independent rows are visited changes.
        let (head, tail) = lu.data.split_at_mut((k + 1) * n);
        let pivot_row = &head[k * n + k + 1..(k + 1) * n];
        let grain = if trailing.saturating_mul(trailing) >= DENSE_PAR_WORK {
            2
        } else {
            usize::MAX
        };
        super::par::for_chunks_mut(tail, n, grain, |_, row| {
            let factor = row[k] / d;
            row[k] = factor;
            if factor != 0.0 {
                for (slot, &pivot) in row[k + 1..].iter_mut().zip(pivot_row) {
                    *slot -= factor * pivot;
                }
            }
        });
    }
    Ok(Lu {
        factors: lu,
        piv,
        sign,
    })
}

impl Lu {
    pub(crate) const fn heap_bytes(&self) -> usize {
        self.factors.data.capacity() * size_of::<f64>() + self.piv.capacity() * size_of::<usize>()
    }

    /// Solves `A x = b`, returning a new solution vector.
    ///
    /// Returns an empty vector for a short right-hand side and ignores any
    /// trailing entries. Use [`Lu::solve_into`] to check shape and finiteness.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
    ) -> Vec<f64> {
        if b.len() < self.factors.rows {
            return Vec::new();
        }
        let mut x = vec![0.0; self.factors.rows];
        self.substitute(b, &mut x);
        x
    }

    /// Solves `A x = b` into caller-provided storage without allocating.
    ///
    /// # Errors
    /// [`DenseError::Shape`] if either slice has the wrong length, or
    /// [`DenseError::NonFinite`] for a NaN or infinite right-hand side.
    /// On an input error, `x` is left untouched.
    pub fn solve_into(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) -> Result<(), DenseError> {
        let n = self.factors.rows;
        if b.len() != n || x.len() != n {
            return Err(DenseError::Shape);
        }
        if b.iter().any(|v| !v.is_finite()) {
            return Err(DenseError::NonFinite);
        }
        self.substitute(b, x);
        Ok(())
    }

    fn substitute(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) {
        for (value, &pivot) in x.iter_mut().zip(&self.piv) {
            *value = b[pivot];
        }
        self.forward_backward(x);
    }

    /// Solves `L U x = b` for an already permuted right-hand side `x`.
    fn forward_backward(
        &self,
        x: &mut [f64],
    ) {
        let n = self.factors.rows;
        for i in 0..n {
            for j in 0..i {
                x[i] -= self.factors.at(i, j) * x[j];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= self.factors.at(i, j) * x[j];
            }
            x[i] /= self.factors.at(i, i);
        }
    }

    /// Solves `A x = e_j` for the `j`-th unit vector `e_j`.
    fn substitute_unit(
        &self,
        j: usize,
        x: &mut [f64],
    ) {
        for (value, &pivot) in x.iter_mut().zip(&self.piv) {
            *value = if pivot == j { 1.0 } else { 0.0 };
        }
        self.forward_backward(x);
    }

    /// Determinant of the factored matrix.
    #[must_use]
    pub fn det(&self) -> f64 {
        (0..self.factors.rows)
            .map(|i| self.factors.at(i, i))
            .product::<f64>()
            * self.sign
    }

    /// Inverse of the factored matrix.
    #[must_use]
    pub fn inverse(&self) -> Mat {
        let n = self.factors.rows;
        let mut inv = Mat::zeros(n, n);
        if n == 0 {
            return inv;
        }
        // One column per task: the substitutions of distinct columns are
        // independent and each keeps the serial left-to-right order.
        let grain = if n.saturating_mul(n).saturating_mul(n) >= DENSE_PAR_FLOPS {
            2
        } else {
            usize::MAX
        };
        let mut columns = crate::kernels::reuse::global_pool().scratch_f64(n * n);
        columns.resize(n * n, 0.0);
        super::par::for_chunks_mut(&mut columns[..], n, grain, |j, column| {
            self.substitute_unit(j, column);
        });
        for (j, column) in columns.chunks_exact(n).enumerate() {
            for (i, &value) in column.iter().enumerate() {
                inv.data[i * n + j] = value;
            }
        }
        inv
    }
}

/// Solves the square system `A x = b` by LU.
///
/// # Errors
/// Propagates the errors of [`lu_factor`] and [`Lu::solve_into`].
pub fn solve(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if b.len() != a.rows {
        return Err(DenseError::Shape);
    }
    let factors = lu_factor(a)?;
    let mut x = vec![0.0; b.len()];
    factors.solve_into(b, &mut x)?;
    Ok(x)
}

/// Cholesky factor `L` (lower triangular) with `A = L L^T`.
///
/// # Errors
/// [`DenseError::Singular`] if `A` is not symmetric positive definite.
pub fn cholesky(a: &Mat) -> Result<Mat, DenseError> {
    crate::cancellation::checkpoint()?;
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    let mut l = Mat::zeros(n, n);
    // Scratch column for the entries below the diagonal of one column. The
    // entries of a column only need the columns to their left, so they are
    // independent and can be computed concurrently; each one subtracts
    // `k = 0..j` in the same order as the sequential sweep, so the bits match.
    let mut column = crate::kernels::reuse::global_pool().scratch_f64(n);
    for j in 0..n {
        crate::cancellation::checkpoint()?;
        let mut s = a.data[j * n + j];
        for k in 0..j {
            s -= l.data[j * n + k] * l.data[j * n + k];
        }
        if s <= 0.0 {
            return Err(DenseError::Singular);
        }
        let diagonal = s.sqrt();
        l.data[j * n + j] = diagonal;
        let rows = n - j - 1;
        if rows == 0 {
            continue;
        }
        let grain = if rows.saturating_mul(j) >= DENSE_PAR_WORK_COLUMN {
            2
        } else {
            usize::MAX
        };
        super::par::map_indexed_into(rows, grain, &mut column, |offset| {
            let i = j + 1 + offset;
            let mut s = a.data[i * n + j];
            for k in 0..j {
                s -= l.data[i * n + k] * l.data[j * n + k];
            }
            s / diagonal
        });
        for (offset, &value) in column.iter().enumerate() {
            l.data[(j + 1 + offset) * n + j] = value;
        }
    }
    Ok(l)
}

/// Solves `A x = b` for symmetric positive definite `A` using Cholesky.
///
/// # Errors
/// As [`cholesky`].
pub fn cholesky_solve(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if b.len() != a.rows {
        return Err(DenseError::Shape);
    }
    let l = cholesky(a)?;
    let mut solution = b.to_vec();
    solve_cholesky_factor_in_place(&l, &mut solution);
    Ok(solution)
}

/// Solves `L x = rhs` with a Cholesky factor from a shape-checked system.
pub(crate) fn solve_cholesky_lower_in_place(
    l: &Mat,
    rhs: &mut [f64],
) {
    for i in 0..l.rows {
        for k in 0..i {
            rhs[i] -= l.at(i, k) * rhs[k];
        }
        rhs[i] /= l.at(i, i);
    }
}

/// Solves `L L^T x = rhs` in caller-owned storage using a retained factor.
pub(crate) fn solve_cholesky_factor_in_place(
    l: &Mat,
    rhs: &mut [f64],
) {
    solve_cholesky_lower_in_place(l, rhs);
    for i in (0..l.rows).rev() {
        for k in i + 1..l.rows {
            rhs[i] -= l.at(k, i) * rhs[k];
        }
        rhs[i] /= l.at(i, i);
    }
}

/// Householder QR: returns `(Q, R)` with `Q` of size `m x m`, `R` of size `m x n`.
#[must_use]
pub fn qr(a: &Mat) -> (Mat, Mat) {
    let (m, n) = (a.rows, a.cols);
    let mut r = a.clone();
    let mut q = Mat::identity(m);
    let pool = crate::kernels::reuse::global_pool();
    let mut v = pool.scratch_f64(m);
    v.resize(m, 0.0);
    let mut weights = pool.scratch_f64(n.max(1));
    weights.resize(n.max(1), 0.0);
    for k in 0..n.min(m.saturating_sub(1)) {
        let mut norm = 0.0;
        for i in k..m {
            norm += r.data[i * n + k] * r.data[i * n + k];
        }
        let norm = norm.sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if r.data[k * n + k] > 0.0 {
            -norm
        } else {
            norm
        };
        // The reflector is zero above row `k`; clear the reused buffer so the
        // full-length sums below see the same terms as a fresh allocation.
        v[..k].fill(0.0);
        for i in k..m {
            v[i] = r.data[i * n + k];
        }
        v[k] -= alpha;
        let vn: f64 = v.iter().map(|x| x * x).sum();
        if vn == 0.0 {
            continue;
        }
        let grain = if (m - k).saturating_mul(n) >= DENSE_PAR_WORK_QR {
            2
        } else {
            usize::MAX
        };
        // Every column of `R` is updated by its own scalar `2 v^T r_j / v^T v`,
        // computed from the pre-update matrix exactly as the column-by-column
        // loop did. Materialising the scalars first lets the rows be updated in
        // parallel while each entry sees the identical single subtraction.
        super::par::map_indexed_into(n, grain, &mut weights, |j| {
            let mut sum = 0.0;
            for i in k..m {
                sum += v[i] * r.data[i * n + j];
            }
            sum * 2.0 / vn
        });
        let (_upper, lower) = r.data.split_at_mut(k * n);
        super::par::for_chunks_mut(lower, n, grain, |offset, row| {
            let factor = v[k + offset];
            for (slot, &weight) in row.iter_mut().zip(weights.iter()) {
                *slot -= weight * factor;
            }
        });
        // Each row of `Q` projects onto the reflector independently, and the
        // projection keeps the serial `j = k..m` order.
        super::par::for_chunks_mut(&mut q.data, m, grain, |_, row| {
            let mut sum = 0.0;
            for j in k..m {
                sum += row[j] * v[j];
            }
            let factor = sum * 2.0 / vn;
            for (slot, &reflector) in row[k..].iter_mut().zip(&v[k..]) {
                *slot -= factor * reflector;
            }
        });
    }
    (q, r)
}

/// Economy (thin) Householder QR: returns `(Q, R)` with `Q` of size `m x k` and
/// `R` of size `k x n`, where `k = min(m, n)`.
///
/// Avoids materialising an `m x m` orthogonal matrix when `m >> n`.
///
/// # Errors
/// Returns [`DenseError::Shape`] if `a` has inconsistent internal storage.
pub fn qr_thin(a: &Mat) -> Result<(Mat, Mat), DenseError> {
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    let (m, n) = (a.rows, a.cols);
    let k_dim = m.min(n);
    if k_dim == 0 {
        return Ok((Mat::zeros(m, 0), Mat::zeros(0, n)));
    }
    let mut r_full = a.clone();
    let mut reflectors: Vec<(usize, Vec<f64>, f64)> = Vec::with_capacity(k_dim);
    let pool = crate::kernels::reuse::global_pool();
    let mut v = pool.scratch_f64(m);
    v.resize(m, 0.0);
    let mut weights = pool.scratch_f64(n.max(1));
    weights.resize(n.max(1), 0.0);

    for k in 0..n.min(m.saturating_sub(1)) {
        let mut norm_sq = 0.0;
        for i in k..m {
            let val = r_full.data[i * n + k];
            norm_sq += val * val;
        }
        let norm = norm_sq.sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if r_full.data[k * n + k] > 0.0 {
            -norm
        } else {
            norm
        };
        v[..k].fill(0.0);
        for i in k..m {
            v[i] = r_full.data[i * n + k];
        }
        v[k] -= alpha;
        let vn: f64 = v[k..m].iter().map(|x| x * x).sum();
        if vn == 0.0 {
            continue;
        }
        let grain = if (m - k).saturating_mul(n) >= DENSE_PAR_WORK_QR {
            2
        } else {
            usize::MAX
        };
        super::par::map_indexed_into(n, grain, &mut weights, |j| {
            let mut sum = 0.0;
            for i in k..m {
                sum += v[i] * r_full.data[i * n + j];
            }
            sum * 2.0 / vn
        });
        let (_upper, lower) = r_full.data.split_at_mut(k * n);
        super::par::for_chunks_mut(lower, n, grain, |offset, row| {
            let factor = v[k + offset];
            for (slot, &weight) in row.iter_mut().zip(weights.iter()) {
                *slot -= weight * factor;
            }
        });
        reflectors.push((k, v[k..m].to_vec(), 2.0 / vn));
    }

    // Build thin Q (m x k_dim) by applying reflectors in reverse to I_{m x k_dim}.
    let mut q = Mat::zeros(m, k_dim);
    for i in 0..k_dim {
        q.data[i * k_dim + i] = 1.0;
    }
    let mut col_dots = pool.scratch_f64(k_dim);
    col_dots.resize(k_dim, 0.0);
    for (k, ref_v, beta) in reflectors.into_iter().rev() {
        let grain = if (m - k).saturating_mul(k_dim) >= DENSE_PAR_WORK_QR {
            2
        } else {
            usize::MAX
        };
        super::par::map_indexed_into(k_dim, grain, &mut col_dots, |j| {
            let mut sum = 0.0;
            for (off, &vi) in ref_v.iter().enumerate() {
                sum += vi * q.data[(k + off) * k_dim + j];
            }
            sum * beta
        });
        let (_upper, lower) = q.data.split_at_mut(k * k_dim);
        super::par::for_chunks_mut(lower, k_dim, grain, |off, row| {
            let vi = ref_v[off];
            for (slot, &dot) in row.iter_mut().zip(col_dots.iter()) {
                *slot -= vi * dot;
            }
        });
    }

    let r_thin = Mat {
        rows: k_dim,
        cols: n,
        data: r_full.data[..k_dim * n].to_vec(),
    };
    Ok((q, r_thin))
}

/// Solves a dense square linear system `A x = b` given as row slices, checking
/// that the solution is finite.
///
/// # Errors
/// Returns [`DenseError::Shape`] for mismatched dimensions,
/// [`DenseError::Singular`] if the matrix is singular or produces a non-finite
/// solution, or [`DenseError::NonFinite`] if inputs contain NaN/Inf.
pub fn solve_linear_rows(
    rows: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let mat = Mat::from_rows(rows)?;
    if mat.rows != mat.cols || mat.rows != rhs.len() {
        return Err(DenseError::Shape);
    }
    if mat.rows == 0 {
        return Ok(Vec::new());
    }
    if !mat.is_finite() || rhs.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let lu = lu_factor(&mat)?;
    let x = lu.solve(rhs);
    if x.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::Singular);
    }
    Ok(x)
}

/// Solves a general tridiagonal linear system in `O(n)` time via the Thomas algorithm:
/// `sub[i - 1] * x[i - 1] + diag[i] * x[i] + sup[i] * x[i + 1] = rhs[i]`.
///
/// `diag` and `rhs` have length `n`; `sub` and `sup` have length `n.saturating_sub(1)`.
///
/// # Errors
/// Returns [`DenseError::Shape`] for mismatched lengths,
/// [`DenseError::NonFinite`] if any input is non-finite, or
/// [`DenseError::Singular`] on a zero/non-finite pivot.
pub fn solve_tridiagonal(
    sub: &[f64],
    diag: &[f64],
    sup: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let n = diag.len();
    let expected_off = n.saturating_sub(1);
    if rhs.len() != n || sub.len() != expected_off || sup.len() != expected_off {
        return Err(DenseError::Shape);
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    if diag
        .iter()
        .chain(sub)
        .chain(sup)
        .chain(rhs)
        .any(|v| !v.is_finite())
    {
        return Err(DenseError::NonFinite);
    }
    if diag[0] == 0.0 {
        return Err(DenseError::Singular);
    }
    let mut cp = super::reuse::f64_scratch(n);
    cp.resize(n, 0.0);
    let mut dp = super::reuse::f64_scratch(n);
    dp.resize(n, 0.0);

    cp[0] = if n > 1 {
        sup[0] / diag[0]
    } else {
        0.0
    };
    dp[0] = rhs[0] / diag[0];
    if !cp[0].is_finite() || !dp[0].is_finite() {
        return Err(DenseError::Singular);
    }
    for i in 1..n {
        let denom = diag[i] - sub[i - 1] * cp[i - 1];
        if denom == 0.0 || !denom.is_finite() {
            return Err(DenseError::Singular);
        }
        cp[i] = if i + 1 < n {
            sup[i] / denom
        } else {
            0.0
        };
        dp[i] = (rhs[i] - sub[i - 1] * dp[i - 1]) / denom;
        if !cp[i].is_finite() || !dp[i].is_finite() {
            return Err(DenseError::Singular);
        }
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::Singular);
    }
    Ok(x)
}

/// Solves a symmetric positive-definite tridiagonal linear system in `O(n)` time
/// via `LDL^T` factorisation.
///
/// `diag` and `rhs` have length `n`; `off` has length `n.saturating_sub(1)`.
///
/// # Errors
/// Returns [`DenseError::Shape`] for mismatched lengths,
/// [`DenseError::NonFinite`] if any input is non-finite, or
/// [`DenseError::Singular`] if the matrix is not strictly positive definite.
pub fn solve_spd_tridiagonal(
    diag: &[f64],
    off: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let n = diag.len();
    let expected_off = n.saturating_sub(1);
    if rhs.len() != n || off.len() != expected_off {
        return Err(DenseError::Shape);
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    if diag.iter().chain(off).chain(rhs).any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let mut pivot = diag[0];
    let mut forward_rhs = rhs[0];
    if pivot <= 0.0 {
        return Err(DenseError::Singular);
    }
    let mut l = super::reuse::f64_scratch(expected_off);
    l.clear();
    let mut x = Vec::with_capacity(n);
    x.push(forward_rhs / pivot);

    // Only the latest diagonal and unscaled right-hand side feed the
    // next row. Store the scaled solution and the lower factor for the
    // backward pass without keeping a separate diagonal buffer.
    for (i, &coefficient) in off.iter().enumerate() {
        let factor = coefficient / pivot;
        pivot = diag[i + 1] - factor * coefficient;
        forward_rhs = rhs[i + 1] - factor * forward_rhs;
        if pivot <= 0.0 || !pivot.is_finite() || !forward_rhs.is_finite() || !factor.is_finite() {
            return Err(DenseError::Singular);
        }
        l.push(factor);
        x.push(forward_rhs / pivot);
    }
    for i in (0..n - 1).rev() {
        x[i] -= l[i] * x[i + 1];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::Singular);
    }
    Ok(x)
}

/// Least-squares solution of `min |A x - b|_2` (`m >= n`, full column rank) by QR.
///
/// # Errors
/// [`DenseError::Shape`] if `m < n` or `b` has the wrong length,
/// [`DenseError::Singular`] if `A` is rank deficient.
pub fn lstsq(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let (m, n) = (a.rows, a.cols);
    if m < n || b.len() != m {
        return Err(DenseError::Shape);
    }
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() || !b.iter().all(|x| x.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let (q, r) = qr(a);
    let qtb = q.transpose().matvec(b);
    let mut x = vec![0.0; n];
    let scale = (0..n).map(|i| r.at(i, i).abs()).fold(0.0, f64::max);
    for i in (0..n).rev() {
        let mut s = qtb[i];
        for j in i + 1..n {
            s -= r.at(i, j) * x[j];
        }
        if r.at(i, i).abs() <= 1e-14 * scale {
            return Err(DenseError::Singular);
        }
        x[i] = s / r.at(i, i);
    }
    Ok(x)
}

/// Power-of-two factor that brings the largest entry of `a` near 1 when it
/// is far outside `[1e-100, 1e100]` (so that squares in the QR/Jacobi
/// iterations neither overflow nor underflow), else `None`.
fn eigen_scale(a: &Mat) -> Option<f64> {
    let big = a.data.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if big.is_finite() && big > 0.0 && !(1e-100..=1e100).contains(&big) {
        // Use the exact binary exponent, including for subnormal inputs.
        // Rounding log2(MAX) would give 1024, an overflowing scale.
        let power = if big.is_normal() {
            big.to_bits() & 0x7ff0_0000_0000_0000
        } else {
            1 << big.to_bits().ilog2()
        };
        Some(f64::from_bits(power))
    } else {
        None
    }
}

fn scaled(
    a: &Mat,
    factor: f64,
) -> Mat {
    let mut out = a.clone();
    for x in &mut out.data {
        *x /= factor;
    }
    out
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
///
/// Returns eigenvalues in ascending order and the matrix whose columns are
/// the corresponding orthonormal eigenvectors.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input, [`DenseError::NoConvergence`]
/// after 100 sweeps.
pub fn eigen_symmetric(a: &Mat) -> Result<(Vec<f64>, Mat), DenseError> {
    crate::cancellation::checkpoint()?;
    if let Some(f) = eigen_scale(a)
        && a.rows == a.cols
        && a.is_finite()
    {
        let (vals, vecs) = eigen_symmetric(&scaled(a, f))?;
        return Ok((vals.into_iter().map(|v| v * f).collect(), vecs));
    }
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if n == 0 {
        return Ok((Vec::new(), Mat::zeros(0, 0)));
    }
    if n == 1 {
        return Ok((vec![a.at(0, 0)], Mat::identity(1)));
    }
    let mut m = a.clone();
    let mut v = Mat::identity(n);
    for _ in 0..100 {
        crate::cancellation::checkpoint()?;
        let mut off = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    off += m.at(i, j) * m.at(i, j);
                }
            }
        }
        if off.sqrt() <= 1e-15 * m.norm_fro().max(f64::MIN_POSITIVE) {
            let mut idx: Vec<usize> = (0..n).collect();
            idx.sort_by(|&x, &y| m.at(x, x).total_cmp(&m.at(y, y)));
            let vals = idx.iter().map(|&i| m.at(i, i)).collect();
            let mut vs = Mat::zeros(n, n);
            for (c, &i) in idx.iter().enumerate() {
                for r in 0..n {
                    vs.set(r, c, v.at(r, i));
                }
            }
            return Ok((vals, vs));
        }
        for p in 0..n {
            crate::cancellation::checkpoint()?;
            for q in p + 1..n {
                let apq = m.at(p, q);
                if apq == 0.0 {
                    continue;
                }
                let theta = (m.at(q, q) - m.at(p, p)) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + theta.mul_add(theta, 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / t.mul_add(t, 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (kp, kq) = (m.at(k, p), m.at(k, q));
                    m.set(k, p, c * kp - s * kq);
                    m.set(k, q, s * kp + c * kq);
                }
                for k in 0..n {
                    let (pk, qk) = (m.at(p, k), m.at(q, k));
                    m.set(p, k, c * pk - s * qk);
                    m.set(q, k, s * pk + c * qk);
                }
                for k in 0..n {
                    let (kp, kq) = (v.at(k, p), v.at(k, q));
                    v.set(k, p, c * kp - s * kq);
                    v.set(k, q, s * kp + c * kq);
                }
            }
        }
    }
    Err(DenseError::NoConvergence)
}

/// Singular value decomposition result `A = U diag(s) V^T`.
#[derive(Debug, Clone)]
pub struct Svd {
    /// Left singular vectors (`m x k`, `k = min(m, n)`).
    pub u: Mat,
    /// Singular values in descending order.
    pub s: Vec<f64>,
    /// Right singular vectors (`n x k`).
    pub v: Mat,
}

/// Thin singular-value decomposition using the bidiagonal SVD backend.
///
/// Rank-deficient matrices retain orthonormal null-space columns. No Gram
/// matrix is formed, so small singular values do not lose accuracy through
/// squaring the condition number.
///
/// # Errors
/// Returns [`DenseError::Shape`] for inconsistent storage or non-finite input,
/// and [`DenseError::NoConvergence`] if the backend does not converge.
pub fn svd(a: &Mat) -> Result<Svd, DenseError> {
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let k = a.rows.min(a.cols);
    if k == 0 {
        return Ok(Svd {
            u: Mat::zeros(a.rows, k),
            s: Vec::new(),
            v: Mat::zeros(a.cols, k),
        });
    }
    // Bidiagonalization still forms squared local quantities internally.
    // Normalize first so uniformly tiny/huge matrices retain their rank and
    // orthogonality instead of underflowing or overflowing those quantities.
    let scale = a.data.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let scale = if scale == 0.0 { 1.0 } else { scale };
    let input = faer::Mat::from_fn(a.rows, a.cols, |i, j| a.at(i, j) / scale);
    let decomposition = input.thin_svd().map_err(|_| DenseError::NoConvergence)?;
    let singular = decomposition.S().column_vector();
    let s: Vec<f64> = (0..k).map(|i| *singular.get(i) * scale).collect();
    if !s.iter().all(|x| x.is_finite()) {
        return Err(DenseError::NoConvergence);
    }
    let mut u = Mat::zeros(a.rows, k);
    let mut v = Mat::zeros(a.cols, k);
    for i in 0..a.rows {
        for j in 0..k {
            u.set(i, j, *decomposition.U().get(i, j));
        }
    }
    for i in 0..a.cols {
        for j in 0..k {
            v.set(i, j, *decomposition.V().get(i, j));
        }
    }
    Ok(Svd { u, s, v })
}

/// 2-norm condition number `sigma_max / sigma_min` (infinite when singular).
///
/// # Errors
/// As [`svd`].
pub fn cond(a: &Mat) -> Result<f64, DenseError> {
    let s = svd(a)?.s;
    let (hi, lo) = (
        s.first().copied().unwrap_or(0.0),
        s.last().copied().unwrap_or(0.0),
    );
    Ok(if lo == 0.0 {
        f64::INFINITY
    } else {
        hi / lo
    })
}

/// Numerical rank: singular values above `tol` (or the default
/// `max(m, n) * eps * sigma_max` when `tol` is `None`).
///
/// # Errors
/// As [`svd`].
pub fn rank(
    a: &Mat,
    tol: Option<f64>,
) -> Result<usize, DenseError> {
    let s = svd(a)?.s;
    let smax = s.first().copied().unwrap_or(0.0);
    let t = tol.unwrap_or_else(|| a.rows.max(a.cols) as f64 * f64::EPSILON * smax);
    Ok(s.iter().filter(|&&v| v > t).count())
}

/// Orthonormal basis of the null space `{x : A x = 0}` as the columns of the
/// result (`cols x (cols - rank)`), from the SVD.
///
/// `tol` is the singular-value threshold (default `max(m, n) * eps * sigma_max`).
/// An empty matrix has the whole space as null space; a zero matrix too.
///
/// # Errors
/// As [`svd`].
pub fn null_space(
    a: &Mat,
    tol: Option<f64>,
) -> Result<Mat, DenseError> {
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    let n = a.cols;
    if n == 0 {
        return Ok(Mat::zeros(0, 0));
    }
    if a.rows == 0 {
        return Ok(Mat::identity(n));
    }
    // Pad a wide matrix with zero rows so that the thin SVD returns all of V.
    let mut padded = a.clone();
    if a.rows < n {
        padded.data.resize(n * n, 0.0);
        padded.rows = n;
    }
    let d = svd(&padded)?;
    let smax = d.s.first().copied().unwrap_or(0.0);
    let t = tol.unwrap_or_else(|| a.rows.max(n) as f64 * f64::EPSILON * smax);
    let rank = d.s.iter().filter(|&&v| v > t).count();
    let mut out = Mat::zeros(n, n - rank);
    for j in rank..n {
        for i in 0..n {
            out.set(i, j - rank, d.v.at(i, j));
        }
    }
    Ok(out)
}

/// Orthonormal basis of the column space (range) of `A` as the columns of
/// the result (`rows x rank`), from the SVD; `tol` as in [`null_space`].
///
/// # Errors
/// As [`svd`].
pub fn orth(
    a: &Mat,
    tol: Option<f64>,
) -> Result<Mat, DenseError> {
    let d = svd(a)?;
    let smax = d.s.first().copied().unwrap_or(0.0);
    let t = tol.unwrap_or_else(|| a.rows.max(a.cols) as f64 * f64::EPSILON * smax);
    let rank = d.s.iter().filter(|&&v| v > t).count();
    let mut out = Mat::zeros(a.rows, rank);
    for j in 0..rank {
        for i in 0..a.rows {
            out.set(i, j, d.u.at(i, j));
        }
    }
    Ok(out)
}

/// Moore-Penrose pseudo-inverse through the SVD.
///
/// # Errors
/// As [`svd`].
pub fn pinv(a: &Mat) -> Result<Mat, DenseError> {
    let d = svd(a)?;
    let smax = d.s.first().copied().unwrap_or(0.0);
    let tol = a.rows.max(a.cols) as f64 * f64::EPSILON * smax;
    let k = d.s.len();
    let mut out = Mat::zeros(a.cols, a.rows);
    for l in 0..k {
        if d.s[l] <= tol {
            continue;
        }
        for i in 0..a.cols {
            for j in 0..a.rows {
                out.data[i * a.rows + j] += d.v.at(i, l) * d.u.at(j, l) / d.s[l];
            }
        }
    }
    Ok(out)
}

/// Minimum-norm least-squares solution through the SVD (handles rank
/// deficiency).
///
/// # Errors
/// Propagates errors from [`svd`]. Returns [`DenseError::Shape`] for a wrong
/// right-hand side length or [`DenseError::NonFinite`] for non-finite entries.
pub fn lstsq_svd(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if b.len() != a.rows {
        return Err(DenseError::Shape);
    }
    if b.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let d = svd(a)?;
    let k = d.s.len();
    if k == 0 {
        return Ok(vec![0.0; a.cols]);
    }
    let smax = d.s.first().copied().unwrap_or(0.0);
    let tol = a.rows.max(a.cols) as f64 * f64::EPSILON * smax;
    let active = d.s.partition_point(|&sigma| sigma > tol);
    let max_b = b.iter().fold(0.0_f64, |max, v| max.max(v.abs()));
    let scale = if max_b > 0.0 && !(1e-100..=1e100).contains(&max_b) {
        max_b
    } else {
        1.0
    };
    // Apply V Σ⁺ Uᵀ to b directly. A row-major pass projects onto the
    // retained columns of U without materialising either Uᵀ or A⁺.
    let mut coefficients = vec![0.0; k];
    for (row, &rhs) in d.u.data.chunks_exact(k).zip(b) {
        let rhs = rhs / scale;
        for (coefficient, u) in coefficients.iter_mut().take(active).zip(row) {
            *coefficient += u * rhs;
        }
    }
    for (coefficient, &sigma) in coefficients.iter_mut().take(active).zip(&d.s) {
        *coefficient = if crate::kernels::exact_eq(scale, 1.0) {
            *coefficient / sigma
        } else {
            crate::kernels::numeric_common::scale_product(*coefficient, &[scale], &[sigma])
                .map_err(|_| DenseError::NoConvergence)?
        };
    }
    Ok(d.v.matvec(&coefficients))
}

/// Result of [`eigen_generalized_symmetric`].
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedEigen {
    /// Eigenvalues `(re, im)`: real and ascending when `b_positive_definite`; otherwise in the
    /// order of [`eigenvalues`] of `B^{-1} A` (complex conjugate pairs may occur).
    pub values: Vec<(f64, f64)>,
    /// Eigenvectors as columns (`B`-orthonormal, `X^T B X = I`), only when `B` is positive
    /// definite.
    pub vectors: Option<Mat>,
    /// Whether the Cholesky factorisation of `B` succeeded (the symmetric-definite reduction was
    /// used).
    pub b_positive_definite: bool,
}

/// Dense symmetric generalised eigenproblem `A x = lambda B x` (`A`, `B` symmetric).
///
/// For positive definite `B` (Cholesky `B = L L^T` succeeds) the problem is reduced to the
/// symmetric `L^{-1} A L^{-T}`; the eigenvalues are real and the eigenvectors `B`-orthonormal.
/// For **indefinite** `B` the Cholesky reduction does not exist, and the eigenvalues are those of
/// `B^{-1} A` (LU solve, then the general QR algorithm), which may be complex pairs; no
/// eigenvectors are returned then.
///
/// # Errors
/// [`DenseError::Shape`] for mismatched or non-square input, [`DenseError::Singular`] when `B`
/// is indefinite and singular, [`DenseError::NoConvergence`] from the QR iteration.
pub fn eigen_generalized_symmetric(
    a: &Mat,
    b: &Mat,
) -> Result<GeneralizedEigen, DenseError> {
    if a.rows != a.cols
        || b.rows != b.cols
        || a.rows != b.rows
        || !a.has_valid_shape()
        || !b.has_valid_shape()
    {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() || !b.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if let Ok(l) = cholesky(b) {
        // C = L^{-1} A L^{-T}: solve L W = A, then C^T = L^{-1} W^T.
        let mut w = a.clone();
        for c in 0..n {
            for i in 0..n {
                let mut s = w.at(i, c);
                for k in 0..i {
                    s -= l.at(i, k) * w.at(k, c);
                }
                w.set(i, c, s / l.at(i, i));
            }
        }
        let mut cmat = w.transpose();
        for c in 0..n {
            for i in 0..n {
                let mut s = cmat.at(i, c);
                for k in 0..i {
                    s -= l.at(i, k) * cmat.at(k, c);
                }
                cmat.set(i, c, s / l.at(i, i));
            }
        }
        for i in 0..n {
            for j in 0..i {
                let s = f64::midpoint(cmat.at(i, j), cmat.at(j, i));
                cmat.set(i, j, s);
                cmat.set(j, i, s);
            }
        }
        let (theta, y) = eigen_symmetric(&cmat)?;
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&p, &q| theta[p].total_cmp(&theta[q]));
        // x = L^{-T} y
        let mut x = Mat::zeros(n, n);
        for (col, &src) in order.iter().enumerate() {
            let mut v: Vec<f64> = (0..n).map(|i| y.at(i, src)).collect();
            for i in (0..n).rev() {
                for k in i + 1..n {
                    v[i] -= l.at(k, i) * v[k];
                }
                v[i] /= l.at(i, i);
            }
            for (i, e) in v.iter().enumerate() {
                x.set(i, col, *e);
            }
        }
        return Ok(GeneralizedEigen {
            values: order.iter().map(|&i| (theta[i], 0.0)).collect(),
            vectors: Some(x),
            b_positive_definite: true,
        });
    }
    // Indefinite B: B^{-1} A.
    let lu = lu_factor(b)?;
    let mut m = Mat::zeros(n, n);
    for c in 0..n {
        let col: Vec<f64> = (0..n).map(|i| a.at(i, c)).collect();
        let sol = lu.solve(&col);
        for (i, e) in sol.iter().enumerate() {
            m.set(i, c, *e);
        }
    }
    Ok(GeneralizedEigen {
        values: eigenvalues(&m)?,
        vectors: None,
        b_positive_definite: false,
    })
}

/// Eigenvalues `(re, im)` of a general real square matrix (Hessenberg
/// reduction by elimination and the shifted QR algorithm of EISPACK `hqr`).
///
/// # Errors
/// [`DenseError::NoConvergence`] when a root needs more than 60 iterations.
pub fn eigenvalues(a: &Mat) -> Result<Vec<(f64, f64)>, DenseError> {
    if let Some(f) = eigen_scale(a)
        && a.rows == a.cols
        && a.is_finite()
    {
        let v = eigenvalues(&scaled(a, f))?;
        return Ok(v.into_iter().map(|(re, im)| (re * f, im * f)).collect());
    }
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![(a.at(0, 0), 0.0)]);
    }
    // 1-based working copy in one flat buffer, `h[i * ld + j]` for
    // `1 <= i, j <= n`; the two guard rows and columns stay zero.
    let ld = n + 2;
    let mut h = vec![0.0; ld * ld];
    for i in 0..n {
        for j in 0..n {
            h[(i + 1) * ld + j + 1] = a.at(i, j);
        }
    }
    // Balancing (power-of-two diagonal similarity) to equalise row/column norms.
    let mut balanced = false;
    while !balanced {
        balanced = true;
        for i in 1..=n {
            let (mut c, mut r) = (0.0, 0.0);
            for j in 1..=n {
                if j != i {
                    c += h[j * ld + i].abs();
                    r += h[i * ld + j].abs();
                }
            }
            if c != 0.0 && r != 0.0 {
                let mut f = 1.0;
                let s = c + r;
                let mut g = r / 2.0;
                while c < g {
                    f *= 2.0;
                    c *= 4.0;
                }
                g = r * 2.0;
                while c > g {
                    f /= 2.0;
                    c /= 4.0;
                }
                if (c + r) / f < 0.95 * s {
                    balanced = false;
                    for j in 1..=n {
                        h[i * ld + j] /= f;
                    }
                    for j in 1..=n {
                        h[j * ld + i] *= f;
                    }
                }
            }
        }
    }
    // Reduction to Hessenberg form by stabilised elementary similarity.
    for m in 2..n {
        let mut x: f64 = 0.0;
        let mut i = m;
        for j in m..=n {
            if h[j * ld + m - 1].abs() > x.abs() {
                x = h[j * ld + m - 1];
                i = j;
            }
        }
        if i != m {
            for j in m - 1..=n {
                h.swap(i * ld + j, m * ld + j);
            }
            for j in 1..=n {
                h.swap(j * ld + i, j * ld + m);
            }
        }
        if x != 0.0 {
            for i in m + 1..=n {
                let mut y = h[i * ld + m - 1];
                if y != 0.0 {
                    y /= x;
                    h[i * ld + m - 1] = y;
                    for j in m..=n {
                        h[i * ld + j] -= y * h[m * ld + j];
                    }
                    for j in 1..=n {
                        h[j * ld + m] += y * h[j * ld + i];
                    }
                }
            }
        }
    }
    for i in 3..=n {
        for j in 1..i - 1 {
            h[i * ld + j] = 0.0;
        }
    }
    let mut anorm = 0.0;
    for i in 1..=n {
        for j in i.saturating_sub(1).max(1)..=n {
            anorm += h[i * ld + j].abs();
        }
    }
    let mut wr = vec![0.0; n + 1];
    let mut wi = vec![0.0; n + 1];
    let mut nn = n;
    let mut t = 0.0;
    let (mut p, mut q, mut r, mut s, mut x, mut y, mut z, mut w);
    while nn >= 1 {
        let mut its = 0;
        loop {
            let mut l = nn;
            while l >= 2 {
                s = h[(l - 1) * ld + l - 1].abs() + h[l * ld + l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if crate::kernels::exact_eq(h[l * ld + l - 1].abs() + s, s) {
                    h[l * ld + l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            x = h[nn * ld + nn];
            if l == nn {
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            y = h[(nn - 1) * ld + nn - 1];
            w = h[nn * ld + nn - 1] * h[(nn - 1) * ld + nn];
            if l == nn - 1 {
                p = 0.5 * (y - x);
                q = p * p + w;
                z = q.abs().sqrt();
                x += t;
                if q >= 0.0 {
                    z = p + z.copysign(p);
                    wr[nn - 1] = x + z;
                    wr[nn] = wr[nn - 1];
                    if z != 0.0 {
                        wr[nn] = x - w / z;
                    }
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = z;
                    wi[nn] = -z;
                }
                nn = nn.saturating_sub(2);
                break;
            }
            if its >= 60 {
                return Err(DenseError::NoConvergence);
            }
            if its == 10 || its == 20 {
                t += x;
                for i in 1..=nn {
                    h[i * ld + i] -= x;
                }
                s = h[nn * ld + nn - 1].abs() + h[(nn - 1) * ld + nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;
            let mut m = nn - 2;
            loop {
                z = h[m * ld + m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[(m + 1) * ld + m] + h[m * ld + m + 1];
                q = h[(m + 1) * ld + m + 1] - z - r - s;
                r = h[(m + 2) * ld + m + 1];
                s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == l {
                    break;
                }
                let u = h[m * ld + m - 1].abs() * (q.abs() + r.abs());
                let v = p.abs()
                    * (h[(m - 1) * ld + m - 1].abs() + z.abs() + h[(m + 1) * ld + m + 1].abs());
                if crate::kernels::exact_eq(u + v, v) {
                    break;
                }
                m -= 1;
            }
            for i in m + 2..=nn {
                h[i * ld + i - 2] = 0.0;
                if i != m + 2 {
                    h[i * ld + i - 3] = 0.0;
                }
            }
            for k in m..nn {
                if k != m {
                    p = h[k * ld + k - 1];
                    q = h[(k + 1) * ld + k - 1];
                    r = if k == nn - 1 {
                        0.0
                    } else {
                        h[(k + 2) * ld + k - 1]
                    };
                    x = p.abs() + q.abs() + r.abs();
                    if x != 0.0 {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }
                s = (p * p + q * q + r * r).sqrt().copysign(p);
                if s != 0.0 {
                    if k == m {
                        if l != m {
                            h[k * ld + k - 1] = -h[k * ld + k - 1];
                        }
                    } else {
                        h[k * ld + k - 1] = -s * x;
                    }
                    p += s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q /= p;
                    r /= p;
                    for j in k..=nn {
                        p = h[k * ld + j] + q * h[(k + 1) * ld + j];
                        if k != nn - 1 {
                            p += r * h[(k + 2) * ld + j];
                            h[(k + 2) * ld + j] -= p * z;
                        }
                        h[(k + 1) * ld + j] -= p * y;
                        h[k * ld + j] -= p * x;
                    }
                    let mmin = if nn < k + 3 { nn } else { k + 3 };
                    for i in l..=mmin {
                        p = x * h[i * ld + k] + y * h[i * ld + k + 1];
                        if k != nn - 1 {
                            p += z * h[i * ld + k + 2];
                            h[i * ld + k + 2] -= p * r;
                        }
                        h[i * ld + k + 1] -= p * q;
                        h[i * ld + k] -= p;
                    }
                }
            }
        }
    }
    Ok((1..=n).map(|i| (wr[i], wi[i])).collect())
}

/// Eigen-decomposition of a general real matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Eigen {
    /// Eigenvalues `(re, im)` in real-Schur order.
    pub values: Vec<(f64, f64)>,
    /// Unit 2-norm eigenvectors, one per entry of `values`. A
    /// complex-conjugate pair gets genuinely complex vectors (the imaginary
    /// parts of the pair have opposite signs).
    pub vectors: Vec<Vec<Complex64>>,
    /// Relative residuals `||A x - lambda x|| / (||A||_F ||x||)`.
    pub residuals: Vec<f64>,
    /// `true` when the matrix is numerically defective: an eigenvalue is
    /// repeated but the computed eigenvectors are parallel, or a pair misses
    /// the residual tolerance. The vectors are then not a full eigenbasis.
    pub defective: bool,
}

/// Eigenvalues and eigenvectors of a general (non-symmetric) real matrix.
///
/// The matrix is first balanced by an exact power-of-two diagonal similarity,
/// reduced to upper Hessenberg form by Householder similarities and driven to
/// real Schur form by the implicitly shifted QR algorithm (Francis double
/// shift). A complex-conjugate pair of a real matrix therefore appears as a
/// 2x2 diagonal block and is returned as two complex eigenvalues with complex
/// eigenvectors; the transformation is accumulated, so the vectors belong to
/// the original matrix. Eigenvectors of the quasi-triangular factor are
/// obtained by block back substitution.
///
/// For a symmetric matrix prefer [`eigen_symmetric`]: it is cheaper and
/// returns real orthonormal vectors.
///
/// A defective matrix (a Jordan block, for instance) has fewer independent
/// eigenvectors than eigenvalues. This routine returns one vector per
/// eigenvalue, flags the result through [`Eigen::defective`] and leaves the
/// dependent vectors parallel; see the module-level limitations.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input, [`DenseError::NoConvergence`]
/// when a Schur eigenvalue needs more than `200 n` QR sweeps.
pub fn eigen(a: &Mat) -> Result<Eigen, DenseError> {
    if let Some(f) = eigen_scale(a)
        && a.rows == a.cols
        && a.is_finite()
    {
        let mut e = eigen(&scaled(a, f))?;
        for v in &mut e.values {
            *v = (v.0 * f, v.1 * f);
        }
        return Ok(e);
    }
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Eigen {
            values: Vec::new(),
            vectors: Vec::new(),
            residuals: Vec::new(),
            defective: false,
        });
    }
    if n == 1 {
        return Ok(Eigen {
            values: vec![(a.at(0, 0), 0.0)],
            vectors: vec![vec![Complex64::new(1.0, 0.0)]],
            residuals: vec![0.0],
            defective: false,
        });
    }
    let (balanced, scale) = balance(a);
    let (h, q0) = hessenberg(&balanced);
    let (t, q) = real_schur(h, q0)?;
    let (values, mut vectors) = schur_eigenpairs(&t, &q);
    // Undo the balancing: `B = D^{-1} A D`, so an eigenvector `y` of `B`
    // gives the eigenvector `D y` of `A`; the eigenvalues are unchanged.
    for x in &mut vectors {
        for (c, s) in x.iter_mut().zip(&scale) {
            *c *= *s;
        }
        let nv = x
            .iter()
            .map(num_complex::Complex::norm_sqr)
            .sum::<f64>()
            .sqrt();
        if nv > 0.0 {
            for c in &mut *x {
                *c /= nv;
            }
        }
    }
    let anorm = a.norm_fro().max(f64::MIN_POSITIVE);
    let residuals: Vec<f64> = values
        .iter()
        .zip(&vectors)
        .map(|(lam, x)| relative_residual(a, *lam, x, anorm))
        .collect();
    let defective = is_defective(&values, &vectors, &residuals);
    Ok(Eigen {
        values,
        vectors,
        residuals,
        defective,
    })
}

/// Balances `a` by an exact power-of-two diagonal similarity
/// `B = D^{-1} A D` and returns `(B, diag(D))`.
fn balance(a: &Mat) -> (Mat, Vec<f64>) {
    let n = a.rows;
    let mut m = a.clone();
    let mut scale = vec![1.0; n];
    for _ in 0..100 {
        let mut changed = false;
        for i in 0..n {
            let mut c = 0.0;
            let mut r = 0.0;
            for j in 0..n {
                if j != i {
                    c += m.at(j, i).abs();
                    r += m.at(i, j).abs();
                }
            }
            if c == 0.0 || r == 0.0 {
                continue;
            }
            let s = c + r;
            let mut f: f64 = 1.0;
            let mut cc = c;
            let mut g = r / 2.0;
            while cc < g {
                f *= 2.0;
                cc *= 4.0;
            }
            g = r * 2.0;
            while cc > g {
                f /= 2.0;
                cc /= 4.0;
            }
            if f.is_finite() && f > 0.0 && (c + r) / f < 0.95 * s {
                changed = true;
                for j in 0..n {
                    m.set(i, j, m.at(i, j) / f);
                }
                for j in 0..n {
                    m.set(j, i, m.at(j, i) * f);
                }
                if let Some(d) = scale.get_mut(i) {
                    *d *= f;
                }
            }
        }
        if !changed {
            break;
        }
    }
    (m, scale)
}

/// Householder reduction to upper Hessenberg form: returns `(H, Q)` with
/// `a = Q H Q^T` and `Q` orthogonal.
fn hessenberg(a: &Mat) -> (Mat, Mat) {
    let n = a.rows;
    let mut h = a.clone();
    let mut q = Mat::identity(n);
    if n < 3 {
        return (h, q);
    }
    let pool = crate::kernels::reuse::global_pool();
    let mut reflector = pool.scratch_f64(n);
    let mut weights = pool.scratch_f64(n);
    reflector.resize(n, 0.0);
    weights.resize(n, 0.0);
    for k in 0..n - 2 {
        let rows = n - k - 1;
        let mut norm = 0.0;
        for i in k + 1..n {
            norm += h.data[i * n + k] * h.data[i * n + k];
        }
        let norm = norm.sqrt();
        if norm == 0.0 {
            continue;
        }
        let x0 = h.data[(k + 1) * n + k];
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        reflector[0] = x0 - alpha;
        for t in 1..rows {
            reflector[t] = h.data[(k + 1 + t) * n + k];
        }
        let v = &reflector[..rows];
        let mut vn2 = 0.0;
        for &value in v {
            vn2 += value * value;
        }
        if vn2 == 0.0 {
            continue;
        }
        let beta = 2.0 / vn2;
        let grain = if rows.saturating_mul(n - k) >= DENSE_PAR_WORK {
            2
        } else {
            usize::MAX
        };
        // Left apply: one scalar per column, `beta * v^T H[k+1.., j]`, from the
        // pre-update matrix; then each row of the block is updated in place.
        // Every entry sees the same single subtraction as the serial loop.
        for j in k..n {
            let mut sum = 0.0;
            for (t, &value) in v.iter().enumerate() {
                sum += value * h.data[(k + 1 + t) * n + j];
            }
            weights[j - k] = sum * beta;
        }
        let column_weights = &weights[..n - k];
        let (_upper, lower) = h.data.split_at_mut((k + 1) * n);
        super::par::for_chunks_mut(lower, n, grain, |offset, row| {
            let value = v[offset];
            for (slot, &weight) in row[k..].iter_mut().zip(column_weights) {
                *slot -= weight * value;
            }
        });
        // Right apply on `H`, then on `Q`: each row's projection is a fixed
        // `t = 0..rows` sum, so the rows are independent.
        super::par::for_chunks_mut(&mut h.data, n, grain, |_, row| {
            let mut sum = 0.0;
            for (t, &value) in v.iter().enumerate() {
                sum += row[k + 1 + t] * value;
            }
            let weight = sum * beta;
            for (slot, &value) in row[k + 1..].iter_mut().zip(v) {
                *slot -= weight * value;
            }
        });
        super::par::for_chunks_mut(&mut q.data, n, grain, |_, row| {
            let mut sum = 0.0;
            for (t, &value) in v.iter().enumerate() {
                sum += row[k + 1 + t] * value;
            }
            let weight = sum * beta;
            for (slot, &value) in row[k + 1..].iter_mut().zip(v) {
                *slot -= weight * value;
            }
        });
        h.data[(k + 1) * n + k] = alpha;
        for i in k + 2..n {
            h.data[i * n + k] = 0.0;
        }
    }
    (h, q)
}

/// Reflector `I - beta v v^T` mapping the vector `x` to `alpha e_1`.
/// Returns `(v, beta)`; `None` when `x` is already a multiple of `e_1`.
fn householder3(
    x: [f64; 3],
    rows: usize,
) -> Option<([f64; 3], f64)> {
    let [a, b, c] = x;
    let c = if rows >= 3 { c } else { 0.0 };
    let norm = (a * a + b * b + c * c).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    let alpha = if a >= 0.0 { -norm } else { norm };
    let [v0, v1, v2] = [a - alpha, b, c];
    let vn2 = v0 * v0 + v1 * v1 + v2 * v2;
    if vn2 == 0.0 {
        return None;
    }
    Some(([v0, v1, v2], 2.0 / vn2))
}

/// `rows k..k+rows` of `h` (columns `j0..`) are multiplied by the reflector.
fn apply_left(
    h: &mut Mat,
    k: usize,
    rows: usize,
    v: &[f64; 3],
    beta: f64,
    j0: usize,
) {
    let ncols = h.cols;
    for j in j0..ncols {
        let s: f64 = v
            .iter()
            .take(rows)
            .enumerate()
            .map(|(t, vt)| vt * h.at(k + t, j))
            .sum::<f64>()
            * beta;
        for (t, vt) in v.iter().take(rows).enumerate() {
            h.set(k + t, j, h.at(k + t, j) - s * vt);
        }
    }
}

/// Columns `k..k+rows` of `h` (rows below `row_limit`) are multiplied by the
/// reflector.
fn apply_right(
    h: &mut Mat,
    k: usize,
    rows: usize,
    v: &[f64; 3],
    beta: f64,
    row_limit: usize,
) {
    let nrows = row_limit.min(h.rows);
    for i in 0..nrows {
        let s: f64 = v
            .iter()
            .take(rows)
            .enumerate()
            .map(|(t, vt)| h.at(i, k + t) * vt)
            .sum::<f64>()
            * beta;
        for (t, vt) in v.iter().take(rows).enumerate() {
            h.set(i, k + t, h.at(i, k + t) - s * vt);
        }
    }
}

/// One implicitly shifted (Francis double-shift) QR sweep on the active block
/// `h[l..m][l..m]`, updating the accumulated orthogonal factor `q`.
fn francis_step(
    h: &mut Mat,
    q: &mut Mat,
    l: usize,
    m: usize,
) {
    let a1 = m - 2;
    let a2 = m - 1;
    let s = h.at(a1, a1) + h.at(a2, a2);
    let t = h.at(a1, a1) * h.at(a2, a2) - h.at(a1, a2) * h.at(a2, a1);
    let mut x = h.at(l, l) * (h.at(l, l) - s) + h.at(l, l + 1) * h.at(l + 1, l) + t;
    let mut y = h.at(l + 1, l) * (h.at(l, l) + h.at(l + 1, l + 1) - s);
    let mut z = h.at(l + 1, l) * h.at(l + 2, l + 1);
    for k in l..m - 1 {
        let rows = (m - k).min(3);
        if let Some((v, beta)) = householder3([x, y, z], rows) {
            let j0 = if k == l { l } else { k - 1 };
            apply_left(h, k, rows, &v, beta, j0);
            apply_right(h, k, rows, &v, beta, m);
            let qrows = q.rows;
            apply_right(q, k, rows, &v, beta, qrows);
            if k > l {
                // The bulge below the subdiagonal of column k-1 is now zero.
                h.set(k + 1, k - 1, 0.0);
                if k + 2 < m {
                    h.set(k + 2, k - 1, 0.0);
                }
            }
        }
        x = h.at(k + 1, k);
        y = if k + 2 < m {
            h.at(k + 2, k)
        } else {
            0.0
        };
        z = if k + 3 < m {
            h.at(k + 3, k)
        } else {
            0.0
        };
    }
}

/// Real Schur form of an upper Hessenberg matrix: returns `(T, Q)` with `T`
/// quasi-upper-triangular and `Q` orthogonal, `H = Q T Q^T`.
fn real_schur(
    mut h: Mat,
    mut q: Mat,
) -> Result<(Mat, Mat), DenseError> {
    let n = h.rows;
    if n < 2 {
        return Ok((h, q));
    }
    let eps = 4.0 * f64::EPSILON;
    let mut m = n;
    let mut sweeps = 0usize;
    let limit = 200 * n;
    while m > 1 {
        let d1 = h.at(m - 1, m - 1).abs() + h.at(m - 2, m - 2).abs();
        if h.at(m - 1, m - 2).abs() <= eps * d1.max(f64::MIN_POSITIVE) {
            h.set(m - 1, m - 2, 0.0);
            m -= 1;
            continue;
        }
        if m >= 3 {
            let d2 = h.at(m - 2, m - 2).abs() + h.at(m - 3, m - 3).abs();
            if h.at(m - 2, m - 3).abs() <= eps * d2.max(f64::MIN_POSITIVE) {
                h.set(m - 2, m - 3, 0.0);
                m -= 2;
                continue;
            }
        }
        if m == 2 {
            break;
        }
        let mut l = m - 1;
        while l > 0 {
            let d = h.at(l, l).abs() + h.at(l - 1, l - 1).abs();
            if h.at(l, l - 1).abs() <= eps * d.max(f64::MIN_POSITIVE) {
                break;
            }
            l -= 1;
        }
        if sweeps >= limit {
            return Err(DenseError::NoConvergence);
        }
        francis_step(&mut h, &mut q, l, m);
        sweeps += 1;
    }
    Ok((h, q))
}

/// Eigenvalues and eigenvectors of a quasi-upper-triangular `t` (with
/// orthogonal `q`). Returns one vector per eigenvalue.
fn schur_eigenpairs(
    t: &Mat,
    q: &Mat,
) -> (Vec<(f64, f64)>, Vec<Vec<Complex64>>) {
    let n = t.rows;
    let mut values = Vec::with_capacity(n);
    let mut vectors = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        let block = if i + 1 < n && t.at(i + 1, i) != 0.0 {
            2
        } else {
            1
        };
        let block_values = if block == 2 {
            let (a, b, c, d) = (
                t.at(i, i),
                t.at(i, i + 1),
                t.at(i + 1, i),
                t.at(i + 1, i + 1),
            );
            let p = f64::midpoint(a, d);
            let det = a * d - b * c;
            let disc = p * p - det;
            if disc >= 0.0 {
                let s = disc.sqrt();
                vec![(p + s, 0.0), (p - s, 0.0)]
            } else {
                let s = (-disc).sqrt();
                vec![(p, s), (p, -s)]
            }
        } else {
            vec![(t.at(i, i), 0.0)]
        };
        for lam in block_values {
            let y = schur_back_substitute(t, i, block, lam);
            let mut x = vec![Complex64::new(0.0, 0.0); n];
            for (row, xi) in x.iter_mut().enumerate() {
                let mut acc = Complex64::new(0.0, 0.0);
                for (j, yj) in y.iter().enumerate() {
                    let qij = q.at(row, j);
                    if qij != 0.0 {
                        acc += *yj * qij;
                    }
                }
                *xi = acc;
            }
            let nv = x
                .iter()
                .map(num_complex::Complex::norm_sqr)
                .sum::<f64>()
                .sqrt();
            if nv > 0.0 {
                for c in &mut *x {
                    *c /= nv;
                }
            }
            values.push(lam);
            vectors.push(x);
        }
        i += block;
    }
    (values, vectors)
}

/// Eigenvector of the quasi-triangular `t` for the eigenvalue `lam` of the
/// diagonal block `(s, s + block - 1)`, by block back substitution.
fn schur_back_substitute(
    t: &Mat,
    s: usize,
    block: usize,
    lam: (f64, f64),
) -> Vec<Complex64> {
    let n = t.rows;
    let l = Complex64::new(lam.0, lam.1);
    let zero = Complex64::new(0.0, 0.0);
    let mut y = vec![zero; n];
    if block == 1 {
        if let Some(slot) = y.get_mut(s) {
            *slot = Complex64::new(1.0, 0.0);
        }
    } else {
        let a = t.at(s, s);
        let b = t.at(s, s + 1);
        let c = t.at(s + 1, s);
        let d = t.at(s + 1, s + 1);
        let mut v0 = Complex64::new(b, 0.0);
        let mut v1 = l - Complex64::new(a, 0.0);
        if v0.norm() + v1.norm() == 0.0 {
            v0 = l - Complex64::new(d, 0.0);
            v1 = Complex64::new(c, 0.0);
        }
        let nv = (v0.norm_sqr() + v1.norm_sqr()).sqrt();
        if nv > 0.0 {
            if let Some(slot) = y.get_mut(s) {
                *slot = v0 / nv;
            }
            if let Some(slot) = y.get_mut(s + 1) {
                *slot = v1 / nv;
            }
        }
    }
    let mut i = s;
    while i > 0 {
        let row = i - 1;
        if row >= 1 && t.at(row, row - 1) != 0.0 {
            let top = row - 1;
            let a11 = Complex64::new(t.at(top, top), 0.0) - l;
            let a12 = Complex64::new(t.at(top, row), 0.0);
            let a21 = Complex64::new(t.at(row, top), 0.0);
            let a22 = Complex64::new(t.at(row, row), 0.0) - l;
            let mut r1 = zero;
            let mut r2 = zero;
            for j in row + 1..n {
                if let Some(yj) = y.get(j) {
                    r1 -= Complex64::new(t.at(top, j), 0.0) * *yj;
                    r2 -= Complex64::new(t.at(row, j), 0.0) * *yj;
                }
            }
            let det = a11 * a22 - a12 * a21;
            let sc = (1.0 + a11.norm() + a12.norm()) * (1.0 + a21.norm() + a22.norm());
            if det.norm() > 1e-14 * sc {
                if let Some(slot) = y.get_mut(top) {
                    *slot = (r1 * a22 - a12 * r2) / det;
                }
                if let Some(slot) = y.get_mut(row) {
                    *slot = (a11 * r2 - r1 * a21) / det;
                }
            }
            i = top;
        } else {
            let mut rhs = zero;
            for j in row + 1..n {
                if let Some(yj) = y.get(j) {
                    rhs -= Complex64::new(t.at(row, j), 0.0) * *yj;
                }
            }
            let piv = Complex64::new(t.at(row, row), 0.0) - l;
            let sc = 1.0 + t.at(row, row).abs() + lam.0.abs() + lam.1.abs();
            if piv.norm() > 1e-14 * sc
                && let Some(slot) = y.get_mut(row)
            {
                *slot = rhs / piv;
            }
            i = row;
        }
    }
    y
}

/// `||A x - lambda x|| / (||A||_F ||x||)` for a unit vector `x`.
fn relative_residual(
    a: &Mat,
    lam: (f64, f64),
    x: &[Complex64],
    anorm: f64,
) -> f64 {
    let n = a.rows;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..n {
        let mut sr = 0.0;
        let mut si = 0.0;
        for j in 0..n {
            let aij = a.at(i, j);
            if let Some(xj) = x.get(j) {
                sr += aij * xj.re;
                si += aij * xj.im;
            }
        }
        if let Some(xi) = x.get(i) {
            let er = sr - (lam.0 * xi.re - lam.1 * xi.im);
            let ei = si - (lam.0 * xi.im + lam.1 * xi.re);
            num += er * er + ei * ei;
            den += xi.norm_sqr();
        }
    }
    let d = anorm * den.sqrt();
    if d > 0.0 {
        num.sqrt() / d
    } else {
        0.0
    }
}

/// Defectiveness: a repeated eigenvalue with parallel eigenvectors, or a pair
/// whose relative residual is not small.
fn is_defective(
    values: &[(f64, f64)],
    vectors: &[Vec<Complex64>],
    residuals: &[f64],
) -> bool {
    if residuals.iter().any(|r| *r > 1e-8) {
        return true;
    }
    for (i, (li, vi)) in values.iter().zip(vectors).enumerate() {
        for (lj, vj) in values.iter().zip(vectors).skip(i + 1) {
            let close = (li.0 - lj.0).abs() <= 1e-8 * (1.0 + li.0.abs())
                && (li.1 - lj.1).abs() <= 1e-8 * (1.0 + li.1.abs());
            if !close {
                continue;
            }
            let ip: Complex64 = vi.iter().zip(vj).map(|(a, b)| a.conj() * *b).sum();
            if ip.norm() > 1.0 - 1e-8 {
                return true;
            }
        }
    }
    false
}

/// Real Schur decomposition `A = Q T Q^T` with orthogonal `Q` and
/// quasi-upper-triangular `T`.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input, [`DenseError::NonFinite`] for
/// non-finite entries, [`DenseError::NoConvergence`] if QR sweeps fail.
pub fn schur(a: &Mat) -> Result<(Mat, Mat), DenseError> {
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if n == 0 {
        return Ok((Mat::zeros(0, 0), Mat::zeros(0, 0)));
    }
    if n == 1 {
        return Ok((Mat::identity(1), a.clone()));
    }
    let (h, q0) = hessenberg(a);
    let (t, q) = real_schur(h, q0)?;
    Ok((q, t))
}

/// Spectral radius estimate through [`eigenvalues`].
///
/// # Errors
/// As [`eigenvalues`].
pub fn spectral_radius(a: &Mat) -> Result<f64, DenseError> {
    Ok(eigenvalues(a)?
        .iter()
        .map(|&(re, im)| re.hypot(im))
        .fold(0.0, f64::max))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat(rows: &[&[f64]]) -> Mat {
        Mat::from_rows(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap()
    }

    /// Independent residual `||A x - lambda x|| / (||A||_F ||x||)`.
    fn residual(
        a: &Mat,
        lam: (f64, f64),
        x: &[Complex64],
    ) -> f64 {
        let n = a.rows;
        let mut xr = vec![0.0; n];
        let mut xi = vec![0.0; n];
        for (i, c) in x.iter().enumerate() {
            xr[i] = c.re;
            xi[i] = c.im;
        }
        let ar = a.matvec(&xr);
        let ai = a.matvec(&xi);
        let mut num = 0.0;
        for i in 0..n {
            let er = ar[i] - (lam.0 * xr[i] - lam.1 * xi[i]);
            let ei = ai[i] - (lam.0 * xi[i] + lam.1 * xr[i]);
            num += er * er + ei * ei;
        }
        let vn = x
            .iter()
            .map(num_complex::Complex::norm_sqr)
            .sum::<f64>()
            .sqrt();
        num.sqrt() / (a.norm_fro() * vn)
    }

    fn assert_residuals(
        a: &Mat,
        e: &Eigen,
    ) {
        assert_eq!(e.values.len(), a.rows);
        assert_eq!(e.vectors.len(), a.rows);
        for ((lam, x), r) in e.values.iter().zip(&e.vectors).zip(&e.residuals) {
            let vn = x
                .iter()
                .map(num_complex::Complex::norm_sqr)
                .sum::<f64>()
                .sqrt();
            assert!((vn - 1.0).abs() <= 1e-12, "vector not unit: {vn}");
            assert!(*r <= 1e-10, "reported residual {r} for {lam:?}");
            let res = residual(a, *lam, x);
            assert!(res <= 1e-10, "actual residual {res} for {lam:?}");
        }
    }

    #[test]
    fn eigen_edge_cases_empty_and_one_by_one() {
        // n == 0
        let a0 = Mat::zeros(0, 0);
        let (vals0, vecs0) = eigen_symmetric(&a0).unwrap();
        assert_eq!(vals0, [] as [f64; 0]);
        assert_eq!(vecs0.rows, 0);
        assert_eq!(vecs0.cols, 0);
        let e0 = eigen(&a0).unwrap();
        assert_eq!(e0.values, [] as [(f64, f64); 0]);
        assert_eq!(
            e0.vectors,
            [] as [std::vec::Vec<num_complex::Complex<f64>>; 0]
        );
        assert_eq!(e0.residuals, [] as [f64; 0]);
        assert!(!e0.defective);

        // n == 1
        let a1 = Mat {
            rows: 1,
            cols: 1,
            data: vec![42.5],
        };
        let (vals1, vecs1) = eigen_symmetric(&a1).unwrap();
        assert_eq!(vals1, vec![42.5]);
        assert_eq!(vecs1.rows, 1);
        assert_eq!(vecs1.cols, 1);
        assert_eq!(vecs1.at(0, 0), 1.0);

        let e1 = eigen(&a1).unwrap();
        assert_eq!(e1.values, vec![(42.5, 0.0)]);
        assert_eq!(e1.vectors, vec![vec![Complex64::new(1.0, 0.0)]]);
        assert_eq!(e1.residuals, vec![0.0]);
        assert!(!e1.defective);

        let ev1 = eigenvalues(&a1).unwrap();
        assert_eq!(ev1, vec![(42.5, 0.0)]);
    }

    #[test]
    fn eigen_symmetric_results_unchanged() {
        let a = mat(&[&[2.0, 1.0], &[1.0, 2.0]]);
        let (vals, vecs) = eigen_symmetric(&a).unwrap();
        assert!((vals[0] - 1.0).abs() < 1e-14);
        assert!((vals[1] - 3.0).abs() < 1e-14);
        // Orthonormal eigenvectors of the Jacobi path, unchanged.
        let c0: Vec<f64> = (0..2).map(|i| vecs.at(i, 0)).collect();
        let c1: Vec<f64> = (0..2).map(|i| vecs.at(i, 1)).collect();
        assert!((c0[0] * c0[0] + c0[1] * c0[1] - 1.0).abs() < 1e-14);
        assert!((c0[0] * c1[0] + c0[1] * c1[1]).abs() < 1e-14);
        // The general solver agrees on a symmetric input.
        let e = eigen(&a).unwrap();
        let mut v: Vec<f64> = e.values.iter().map(|x| x.0).collect();
        v.sort_by(f64::total_cmp);
        assert!((v[0] - 1.0).abs() < 1e-12 && (v[1] - 3.0).abs() < 1e-12);
        for lam in &e.values {
            assert!(lam.1.abs() < 1e-12);
        }
        assert_residuals(&a, &e);
        assert!(!e.defective);
    }

    #[test]
    #[allow(clippy::assert_is_empty)]
    fn eigen_general_residuals() {
        let a = mat(&[
            &[4.0, 1.0, 2.0, 0.0, 1.0],
            &[0.0, 3.0, 1.0, 1.0, 0.0],
            &[2.0, 0.0, 2.0, 1.0, 1.0],
            &[1.0, 0.0, 0.0, 1.0, 2.0],
            &[0.0, 1.0, 0.0, 1.0, 3.0],
        ]);
        let e = eigen(&a).unwrap();
        assert_residuals(&a, &e);
        // The spectrum must agree with the (vector-free) Hessenberg/QR path.
        let mut other = eigenvalues(&a).unwrap();
        for lam in &e.values {
            let best = other
                .iter()
                .map(|o| (lam.0 - o.0).abs().hypot(lam.1 - o.1))
                .fold(f64::INFINITY, f64::min);
            assert!(best <= 1e-8, "spectrum mismatch for {lam:?}: {best}");
            if let Some(pos) = other
                .iter()
                .position(|o| (lam.0 - o.0).abs().hypot(lam.1 - o.1) <= 1e-8)
            {
                other.swap_remove(pos);
            }
        }
        assert!(other.is_empty());
    }

    #[test]
    fn eigen_rotation_complex_pair() {
        let a = mat(&[&[0.0, -1.0], &[1.0, 0.0]]);
        let e = eigen(&a).unwrap();
        assert_eq!(e.values.len(), 2);
        let pos = e.values.iter().position(|v| v.1 > 0.5).unwrap();
        let neg = e.values.iter().position(|v| v.1 < -0.5).unwrap();
        assert!(e.values[pos].0.abs() < 1e-12);
        assert!((e.values[pos].1 - 1.0).abs() < 1e-12);
        assert!(e.values[neg].0.abs() < 1e-12);
        assert!((e.values[neg].1 + 1.0).abs() < 1e-12);
        // Genuinely complex eigenvectors, not real ones relabelled.
        let im_pos = e.vectors[pos]
            .iter()
            .map(|c| c.im.abs())
            .fold(0.0, f64::max);
        let im_neg = e.vectors[neg]
            .iter()
            .map(|c| c.im.abs())
            .fold(0.0, f64::max);
        assert!(im_pos > 1e-3, "real eigenvector for i: {im_pos}");
        assert!(im_neg > 1e-3, "real eigenvector for -i: {im_neg}");
        assert_residuals(&a, &e);
    }

    #[test]
    fn eigen_two_complex_pairs() {
        let a = mat(&[
            &[0.0, -1.0, 0.0, 0.0],
            &[1.0, 0.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0, -1.0],
            &[0.0, 0.0, 1.0, 0.0],
        ]);
        let e = eigen(&a).unwrap();
        assert_eq!(e.values.len(), 4);
        let up = e.values.iter().filter(|v| v.1 > 0.5).count();
        let down = e.values.iter().filter(|v| v.1 < -0.5).count();
        assert_eq!((up, down), (2, 2));
        for v in &e.values {
            assert!(
                v.0.abs() < 1e-12 && (v.1.abs() - 1.0).abs() < 1e-12,
                "{v:?}"
            );
        }
        for x in &e.vectors {
            let im = x.iter().map(|c| c.im.abs()).fold(0.0, f64::max);
            assert!(im > 1e-3, "pair eigenvector is real");
        }
        assert_residuals(&a, &e);
    }

    #[test]
    fn eigen_jordan_block_reports_deficiency() {
        let a = mat(&[&[1.0, 1.0], &[0.0, 1.0]]);
        let e = eigen(&a).unwrap();
        assert_eq!(e.values.len(), 2);
        for v in &e.values {
            assert!((v.0 - 1.0).abs() < 1e-12 && v.1.abs() < 1e-12);
        }
        assert!(e.defective, "Jordan block must be reported as defective");
        // Only one independent eigenvector exists; the second pair cannot
        // satisfy the eigen-equation, and the residual says so.
        let worst = e.residuals.iter().fold(0.0_f64, |m, r| m.max(*r));
        assert!(worst > 1e-6, "residuals hide the deficiency: {worst}");
    }

    #[test]
    fn eigen_scaling_invariance() {
        let a = mat(&[
            &[1.0, 2.0, 3.0, 4.0],
            &[0.5, 2.0, 1.0, 0.0],
            &[1.0, 0.0, 3.0, 1.0],
            &[0.0, 1.0, 2.0, 4.0],
        ]);
        let d = [1e6, 1e-6, 1e3, 1e-3];
        let rows: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| a.at(i, j) * d[i] / d[j]).collect())
            .collect();
        let b = Mat::from_rows(&rows).unwrap();
        let ea = eigen(&a).unwrap();
        let eb = eigen(&b).unwrap();
        assert_eq!(ea.values.len(), eb.values.len());
        for lam in &ea.values {
            let best = eb
                .values
                .iter()
                .map(|o| (lam.0 - o.0).abs().hypot(lam.1 - o.1))
                .fold(f64::INFINITY, f64::min);
            assert!(
                best <= 1e-10 * (1.0 + lam.0.hypot(lam.1)),
                "balancing changed the spectrum: {lam:?} -> {best}"
            );
        }
        assert_residuals(&b, &eb);
    }

    /// Raw bits of a slice, for exact comparisons.
    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|value| value.to_bits()).collect()
    }

    /// Deterministic xorshift values in `[-1, 1)`, so the kernels are fed
    /// non-trivial data without a dependency or a fixed literal table.
    fn sample(
        seed: u64,
        len: usize,
    ) -> Vec<f64> {
        let mut state = seed | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
            })
            .collect()
    }

    /// Strictly diagonally dominant matrix: every Gershgorin disk is centred on
    /// a positive diagonal entry with a smaller radius, so `LU` never meets a
    /// zero pivot and the spectrum is real and positive.
    fn well_conditioned(
        n: usize,
        seed: u64,
    ) -> Mat {
        let mut data = sample(seed, n * n);
        for i in 0..n {
            let mut radius = 0.0;
            for j in 0..n {
                if i != j {
                    radius += data[i * n + j].abs();
                }
            }
            data[i * n + i] = radius + 1.0;
        }
        Mat {
            rows: n,
            cols: n,
            data,
        }
    }

    /// Symmetric positive definite `B^T B + n I`, formed without the kernel
    /// under test.
    fn spd(
        n: usize,
        seed: u64,
    ) -> Mat {
        let b = sample(seed, n * n);
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += b[k * n + i] * b[k * n + j];
                }
                data[i * n + j] = sum;
            }
            data[i * n + i] += n as f64;
        }
        Mat {
            rows: n,
            cols: n,
            data,
        }
    }

    /// The plain triple loop `C = A B`, with the serial accumulation order.
    fn scalar_matmul(
        rows: usize,
        inner: usize,
        cols: usize,
        a: &[f64],
        b: &[f64],
    ) -> Vec<f64> {
        let mut out = vec![0.0; rows * cols];
        for i in 0..rows {
            for k in 0..inner {
                let factor = a[i * inner + k];
                if factor != 0.0 {
                    for j in 0..cols {
                        out[i * cols + j] += factor * b[k * cols + j];
                    }
                }
            }
        }
        out
    }

    /// The pre-optimisation sequential `lu_factor`, kept as the oracle.
    fn reference_lu(a: &Mat) -> (Mat, Vec<usize>, f64) {
        let n = a.rows;
        let mut lu = a.clone();
        let mut piv: Vec<usize> = (0..n).collect();
        let mut sign = 1.0;
        for k in 0..n {
            let mut p = k;
            for i in k + 1..n {
                if lu.at(i, k).abs() > lu.at(p, k).abs() {
                    p = i;
                }
            }
            assert_ne!(lu.at(p, k), 0.0, "oracle fed a singular matrix");
            if p != k {
                for j in 0..n {
                    lu.data.swap(k * n + j, p * n + j);
                }
                piv.swap(k, p);
                sign = -sign;
            }
            let d = lu.at(k, k);
            for i in k + 1..n {
                let m = lu.at(i, k) / d;
                lu.data[i * n + k] = m;
                if m != 0.0 {
                    for j in k + 1..n {
                        lu.data[i * n + j] -= m * lu.data[k * n + j];
                    }
                }
            }
        }
        (lu, piv, sign)
    }

    /// The pre-optimisation sequential `cholesky`, kept as the oracle.
    fn reference_cholesky(a: &Mat) -> Result<Mat, DenseError> {
        let n = a.rows;
        let mut l = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let mut s = a.at(i, j);
                for k in 0..j {
                    s -= l.at(i, k) * l.at(j, k);
                }
                if i == j {
                    if s <= 0.0 {
                        return Err(DenseError::Singular);
                    }
                    l.set(i, i, s.sqrt());
                } else {
                    l.set(i, j, s / l.at(j, j));
                }
            }
        }
        Ok(l)
    }

    /// The pre-optimisation sequential `qr`, kept as the oracle.
    fn reference_qr(a: &Mat) -> (Mat, Mat) {
        let (m, n) = (a.rows, a.cols);
        let mut r = a.clone();
        let mut q = Mat::identity(m);
        for k in 0..n.min(m.saturating_sub(1)) {
            let mut norm = 0.0;
            for i in k..m {
                norm += r.at(i, k) * r.at(i, k);
            }
            let norm = norm.sqrt();
            if norm == 0.0 {
                continue;
            }
            let alpha = if r.at(k, k) > 0.0 {
                -norm
            } else {
                norm
            };
            let mut v = vec![0.0; m];
            for i in k..m {
                v[i] = r.at(i, k);
            }
            v[k] -= alpha;
            let vn: f64 = v.iter().map(|x| x * x).sum();
            if vn == 0.0 {
                continue;
            }
            for j in 0..n {
                let d: f64 = (k..m).map(|i| v[i] * r.at(i, j)).sum::<f64>() * 2.0 / vn;
                for i in k..m {
                    r.data[i * n + j] -= d * v[i];
                }
            }
            for i in 0..m {
                let d: f64 = (k..m).map(|j| q.at(i, j) * v[j]).sum::<f64>() * 2.0 / vn;
                for j in k..m {
                    q.data[i * m + j] -= d * v[j];
                }
            }
        }
        (q, r)
    }

    /// Serial forward/back substitution with a known pivot order.
    fn reference_solve(
        factors: &Mat,
        piv: &[usize],
        b: &[f64],
    ) -> Vec<f64> {
        let n = factors.rows;
        let mut x = vec![0.0; n];
        for (value, &pivot) in x.iter_mut().zip(piv) {
            *value = b[pivot];
        }
        for i in 0..n {
            for j in 0..i {
                x[i] -= factors.at(i, j) * x[j];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= factors.at(i, j) * x[j];
            }
            x[i] /= factors.at(i, i);
        }
        x
    }

    /// Every kernel result that the blocked/parallel paths touch, as raw bits.
    fn kernel_bits(
        a: &Mat,
        s: &Mat,
        rhs: &[f64],
    ) -> Vec<u64> {
        let mut out = Vec::new();
        out.extend(bits(&a.matmul(s).expect("shapes agree").data));
        let lu = lu_factor(a).expect("nonsingular");
        out.extend(bits(&lu.inverse().data));
        out.extend(bits(&lu.solve(rhs)));
        out.push(lu.det().to_bits());
        out.extend(bits(&cholesky(s).expect("spd").data));
        let (q, r) = qr(a);
        out.extend(bits(&q.data));
        out.extend(bits(&r.data));
        out.extend(
            eigenvalues(a)
                .expect("converges")
                .iter()
                .flat_map(|value| [value.0.to_bits(), value.1.to_bits()]),
        );
        out
    }

    /// The two routines whose parallel path needs a sizeable matrix.
    fn factorization_bits(a: &Mat) -> Vec<u64> {
        let lu = lu_factor(a).expect("nonsingular");
        let (q, r) = qr(a);
        let mut out = bits(&lu.factors.data);
        out.extend(bits(&q.data));
        out.extend(bits(&r.data));
        out
    }

    #[test]
    fn blocked_kernels_match_the_sequential_oracle_bit_for_bit() {
        /// Runs every kernel and compares against the sequential oracle.
        fn check_oracle(
            n: usize,
            a: &Mat,
            s: &Mat,
            rhs: &[f64],
        ) {
            for (rows, inner, cols) in [(n, n, n), (n, n, 3), (n, 3, n), (n, n, 1)] {
                let left = Mat {
                    rows,
                    cols: inner,
                    data: sample(0x5eed_3000 + rows as u64, rows * inner),
                };
                let right = Mat {
                    rows: inner,
                    cols,
                    data: sample(0x5eed_4000 + cols as u64, inner * cols),
                };
                let product = left.matmul(&right).expect("shapes agree");
                let naive = scalar_matmul(rows, inner, cols, &left.data, &right.data);
                assert_eq!(
                    bits(&product.data),
                    bits(&naive),
                    "matmul {rows}x{inner}x{cols}"
                );
            }

            let lu = lu_factor(a).expect("nonsingular");
            let (factors, piv, sign) = reference_lu(a);
            assert_eq!(bits(&lu.factors.data), bits(&factors.data), "lu n={n}");
            assert_eq!(lu.piv, piv, "piv n={n}");
            assert_eq!(lu.sign.to_bits(), sign.to_bits(), "sign n={n}");
            assert_eq!(
                bits(&lu.solve(rhs)),
                bits(&reference_solve(&factors, &piv, rhs))
            );

            let mut reference_inverse = vec![0.0; n * n];
            for j in 0..n {
                let mut unit = vec![0.0; n];
                unit[j] = 1.0;
                for (i, value) in reference_solve(&factors, &piv, &unit).iter().enumerate() {
                    reference_inverse[i * n + j] = *value;
                }
            }
            assert_eq!(
                bits(&lu.inverse().data),
                bits(&reference_inverse),
                "inverse n={n}"
            );

            let l = cholesky(s).expect("spd");
            let oracle = reference_cholesky(s).expect("spd");
            assert_eq!(bits(&l.data), bits(&oracle.data), "cholesky n={n}");

            let (q, r) = qr(a);
            let (oracle_q, oracle_r) = reference_qr(a);
            assert_eq!(bits(&q.data), bits(&oracle_q.data), "qr q n={n}");
            assert_eq!(bits(&r.data), bits(&oracle_r.data), "qr r n={n}");
        }

        for n in [1usize, 3, 8, 33, 200] {
            let a = well_conditioned(n, 0x5eed_0000 + n as u64);
            let s = spd(n, 0x5eed_1000 + n as u64);
            let rhs = sample(0x5eed_2000 + n as u64, n);
            crate::execution::with_scoped_pool(Some(1), || check_oracle(n, &a, &s, &rhs));
            crate::execution::with_scoped_pool(Some(4), || check_oracle(n, &a, &s, &rhs));
        }
    }

    #[test]
    fn blocked_kernels_are_bit_identical_for_one_and_four_threads() {
        for n in [1usize, 3, 8, 33, 200] {
            let a = well_conditioned(n, 0xdead_0000 + n as u64);
            let s = spd(n, 0xdead_1000 + n as u64);
            let rhs = sample(0xdead_2000 + n as u64, n);
            let serial = crate::execution::with_scoped_pool(Some(1), || kernel_bits(&a, &s, &rhs));
            let parallel =
                crate::execution::with_scoped_pool(Some(4), || kernel_bits(&a, &s, &rhs));
            assert_eq!(serial, parallel, "thread count changed the bits at n={n}");
        }
        // `504^2 > DENSE_PAR_WORK`, so the LU trailing update and the QR row
        // update really do split at four threads.
        let n = 504;
        let a = well_conditioned(n, 0xdead_5000);
        let serial = crate::execution::with_scoped_pool(Some(1), || factorization_bits(&a));
        let parallel = crate::execution::with_scoped_pool(Some(4), || factorization_bits(&a));
        assert_eq!(serial, parallel, "thread count changed the bits at n={n}");
    }

    #[test]
    fn hessenberg_qr_spectrum_matches_the_schur_path() {
        for n in [3usize, 8, 17, 33] {
            let a = well_conditioned(n, 0x51ed_0000 + n as u64);
            let mut direct = eigenvalues(&a).expect("converges");
            let schur = eigen(&a).expect("converges");
            assert_eq!(direct.len(), n);
            for lam in &schur.values {
                let tolerance = 1e-8 * (1.0 + lam.0.hypot(lam.1));
                let best = direct
                    .iter()
                    .map(|other| (lam.0 - other.0).abs().hypot(lam.1 - other.1))
                    .fold(f64::INFINITY, f64::min);
                assert!(best <= tolerance, "n={n}: {lam:?} is {best} away");
                if let Some(position) = direct
                    .iter()
                    .position(|other| (lam.0 - other.0).abs().hypot(lam.1 - other.1) <= tolerance)
                {
                    direct.swap_remove(position);
                }
            }
            assert!(direct.is_empty(), "n={n}: unmatched eigenvalues left over");
        }
    }
}

/// Matrix norm kind for [`norm`] and [`cond_norm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixNorm {
    /// Maximum absolute column sum: `max_j sum_i |A_{i,j}|`.
    One,
    /// Maximum absolute row sum: `max_i sum_j |A_{i,j}|`.
    Inf,
    /// Frobenius norm: `sqrt(sum_{i,j} |A_{i,j}|^2)`.
    Frobenius,
    /// Spectral (2-norm): largest singular value `sigma_1(A)`.
    Spectral,
}

/// Matrix norm with checked shape and finite entries.
///
/// # Errors
/// Returns shape, non-finite input, or SVD convergence errors.
pub fn norm(
    a: &Mat,
    kind: MatrixNorm,
) -> Result<f64, DenseError> {
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    if a.rows == 0 || a.cols == 0 {
        return Ok(0.0);
    }
    Ok(match kind {
        | MatrixNorm::One => crate::kernels::matfun::norm1(a),
        | MatrixNorm::Inf => crate::kernels::schur::norm_inf(a),
        | MatrixNorm::Frobenius => a.norm_fro(),
        | MatrixNorm::Spectral => svd(a)?.s.first().copied().unwrap_or(0.0),
    })
}

/// Condition number in the selected norm; singular matrices have infinite condition number.
///
/// # Errors
/// Returns shape, non-finite input, or factorization errors.
pub fn cond_norm(
    a: &Mat,
    kind: MatrixNorm,
) -> Result<f64, DenseError> {
    let norm_a = norm(a, kind)?;
    if kind == MatrixNorm::Spectral {
        return cond(a);
    }
    if a.rows != a.cols || !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if a.rows == 0 {
        return Ok(1.0);
    }
    let lu = match lu_factor(a) {
        | Ok(lu) => lu,
        | Err(DenseError::Singular) => return Ok(f64::INFINITY),
        | Err(error) => return Err(error),
    };
    Ok(norm_a * norm(&lu.inverse(), kind)?)
}

/// Moore-Penrose pseudoinverse with a relative singular-value threshold.
///
/// # Errors
/// Returns non-finite input for invalid relative tolerance, or SVD errors.
pub fn pinv_tol(
    a: &Mat,
    rcond: f64,
) -> Result<Mat, DenseError> {
    if !rcond.is_finite() || rcond < 0.0 {
        return Err(DenseError::NonFinite);
    }
    if !a.has_valid_shape() {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    if a.rows == 0 || a.cols == 0 {
        return Ok(Mat::zeros(a.cols, a.rows));
    }
    let d = svd(a)?;
    let cutoff = rcond * d.s.first().copied().unwrap_or(0.0);
    let mut out = Mat::zeros(a.cols, a.rows);
    for r in 0..d.s.len() {
        if d.s[r] <= cutoff {
            continue;
        }
        for i in 0..a.cols {
            for j in 0..a.rows {
                out.data[i * a.rows + j] += d.v.at(i, r) * d.u.at(j, r) / d.s[r];
            }
        }
    }
    Ok(out)
}

/// Optimized host-compatible dense operations preserving blocked faer arithmetic.
pub mod optimized;

/// Exact rational dense factorizations and perfect-square arithmetic.
pub mod exact;
/// Sequential reference operations with stable row-major layouts.
pub mod reference;

#[cfg(test)]
mod extra_tests {
    use super::*;

    #[test]
    fn tiled_transpose_matches_naive_large() {
        let rows = 97;
        let cols = 131;
        let data: Vec<f64> = (0..rows * cols).map(|i| (i as f64) * 0.125 - 7.0).collect();
        let m = Mat { rows, cols, data };
        let t = m.transpose();
        assert_eq!(t.rows, cols);
        assert_eq!(t.cols, rows);
        for i in 0..rows {
            for j in 0..cols {
                assert_eq!(t.at(j, i).to_bits(), m.at(i, j).to_bits());
            }
        }
    }

    #[test]
    fn qr_thin_reconstructs_tall_matrix() {
        let m = 20;
        let n = 5;
        let data: Vec<f64> = (0..m * n)
            .map(|k| ((k as f64) * 0.7).sin() + 0.1 * (k as f64))
            .collect();
        let a = Mat {
            rows: m,
            cols: n,
            data,
        };
        let (q, r) = qr_thin(&a).unwrap();
        assert_eq!((q.rows, q.cols), (m, n));
        assert_eq!((r.rows, r.cols), (n, n));
        let recon = q.matmul(&r).unwrap();
        for (x, y) in recon.data.iter().zip(&a.data) {
            assert!((x - y).abs() < 1e-11);
        }
        let qtq = q.transpose().matmul(&q).unwrap();
        for i in 0..n {
            for j in 0..n {
                let target = if i == j { 1.0 } else { 0.0 };
                assert!((qtq.at(i, j) - target).abs() < 1e-11);
            }
        }
    }

    #[test]
    fn tridiagonal_solvers_agree_with_exact_solution() {
        let sub = vec![-1.0, -1.0, -1.0, -1.0];
        let diag = vec![4.0, 4.0, 4.0, 4.0, 4.0];
        let sup = vec![-1.0, -1.0, -1.0, -1.0];
        let x_true = [1.0, -2.0, 3.0, -4.0, 5.0];
        let mut rhs = vec![0.0; 5];
        for i in 0..5 {
            rhs[i] = diag[i] * x_true[i];
            if i > 0 {
                rhs[i] += sub[i - 1] * x_true[i - 1];
            }
            if i + 1 < 5 {
                rhs[i] += sup[i] * x_true[i + 1];
            }
        }
        let x_gen = solve_tridiagonal(&sub, &diag, &sup, &rhs).unwrap();
        let x_spd = solve_spd_tridiagonal(&diag, &sub, &rhs).unwrap();
        for i in 0..5 {
            assert!((x_gen[i] - x_true[i]).abs() < 1e-13);
            assert!((x_spd[i] - x_true[i]).abs() < 1e-13);
        }
    }
}
