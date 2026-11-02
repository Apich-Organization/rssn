//! # Dense linear algebra on row-major `f64` matrices
//!
//! A small self-contained dense toolkit: LU with partial pivoting,
//! Cholesky, Householder QR, least squares, the Jacobi eigenvalue algorithm
//! for symmetric matrices, a one-sided Jacobi SVD (singular values,
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
//! Parallelism: [`Mat::matmul`] and [`Mat::matvec`] run row by row concurrently above a size
//! threshold; every row keeps its sequential accumulation order, so the results are bit-identical
//! to a one-thread run. [`lu_factor`] stays sequential: the per-step dispatch cost of a parallel
//! trailing update did not pay at 700 x 700.

use std::fmt;

use num_complex::Complex64;

/// Errors of the dense routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenseError {
    /// Shapes of the operands do not match.
    Shape,
    /// The matrix is singular (or not positive definite for Cholesky).
    Singular,
    /// An iteration did not converge.
    NoConvergence,
}

impl fmt::Display for DenseError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::Shape => write!(f, "incompatible shapes"),
            | Self::Singular => write!(f, "matrix is singular or not positive definite"),
            | Self::NoConvergence => write!(f, "iteration did not converge"),
        }
    }
}

impl std::error::Error for DenseError {}

/// Multiply-adds from which dense products run in parallel.
const DENSE_PAR_FLOPS: usize = 4_000_000;

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

    /// Transpose.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let mut t = Self::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.data[j * self.rows + i] = self.data[i * self.cols + j];
            }
        }
        t
    }

    /// Matrix-vector product.
    #[must_use]
    pub fn matvec(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        let row = |i: usize| -> f64 {
            self.data[i * self.cols..(i + 1) * self.cols]
                .iter()
                .zip(x)
                .map(|(a, b)| a * b)
                .sum()
        };
        if self.rows * self.cols >= DENSE_PAR_FLOPS / 8 {
            super::par::map_indexed(self.rows, 2, row)
        } else {
            (0..self.rows).map(row).collect()
        }
    }

    /// Matrix product.
    ///
    /// # Errors
    /// Returns [`DenseError::Shape`] on mismatched inner dimensions.
    pub fn matmul(
        &self,
        other: &Self,
    ) -> Result<Self, DenseError> {
        if self.cols != other.rows {
            return Err(DenseError::Shape);
        }
        let mut out = Self::zeros(self.rows, other.cols);
        let (m, oc) = (self.cols, other.cols);
        if oc == 0 {
            return Ok(out);
        }
        // one output row per task; the row's own accumulation order is the
        // serial one, so the product is bit-identical
        let grain = if self.rows * m * oc >= DENSE_PAR_FLOPS {
            2
        } else {
            usize::MAX
        };
        super::par::for_chunks_mut(&mut out.data, oc, grain, |i, out_row| {
            for k in 0..m {
                let a = self.data[i * m + k];
                if a != 0.0 {
                    for j in 0..oc {
                        out_row[j] += a * other.data[k * oc + j];
                    }
                }
            }
        });
        Ok(out)
    }

    /// Frobenius norm.
    #[must_use]
    pub fn norm_fro(&self) -> f64 {
        self.data.iter().map(|v| v * v).sum::<f64>().sqrt()
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
/// [`DenseError::Shape`] for non-square input, [`DenseError::Singular`]
/// when a pivot is exactly zero.
pub fn lu_factor(a: &Mat) -> Result<Lu, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
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
    Ok(Lu {
        factors: lu,
        piv,
        sign,
    })
}

impl Lu {
    /// Solves `A x = b`.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
    ) -> Vec<f64> {
        let n = self.factors.rows;
        let mut x: Vec<f64> = self.piv.iter().map(|&p| b[p]).collect();
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
        x
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
        for j in 0..n {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let c = self.solve(&e);
            for i in 0..n {
                inv.data[i * n + j] = c[i];
            }
        }
        inv
    }
}

/// Solves the square system `A x = b` by LU.
///
/// # Errors
/// Propagates the errors of [`lu_factor`]; [`DenseError::Shape`] if `b` has the wrong length.
pub fn solve(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if b.len() != a.rows {
        return Err(DenseError::Shape);
    }
    Ok(lu_factor(a)?.solve(b))
}

/// Cholesky factor `L` (lower triangular) with `A = L L^T`.
///
/// # Errors
/// [`DenseError::Singular`] if `A` is not symmetric positive definite.
pub fn cholesky(a: &Mat) -> Result<Mat, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
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

/// Solves `A x = b` for symmetric positive definite `A` using Cholesky.
///
/// # Errors
/// As [`cholesky`].
pub fn cholesky_solve(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let l = cholesky(a)?;
    let n = a.rows;
    let mut y = b.to_vec();
    for i in 0..n {
        for k in 0..i {
            y[i] -= l.at(i, k) * y[k];
        }
        y[i] /= l.at(i, i);
    }
    for i in (0..n).rev() {
        for k in i + 1..n {
            y[i] -= l.at(k, i) * y[k];
        }
        y[i] /= l.at(i, i);
    }
    Ok(y)
}

/// Householder QR: returns `(Q, R)` with `Q` of size `m x m`, `R` of size `m x n`.
#[must_use]
pub fn qr(a: &Mat) -> (Mat, Mat) {
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

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
///
/// Returns eigenvalues in ascending order and the matrix whose columns are
/// the corresponding orthonormal eigenvectors.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input, [`DenseError::NoConvergence`]
/// after 100 sweeps.
pub fn eigen_symmetric(a: &Mat) -> Result<(Vec<f64>, Mat), DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    let mut m = a.clone();
    let mut v = Mat::identity(n);
    for _ in 0..100 {
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

/// One-sided Jacobi SVD.
///
/// # Errors
/// [`DenseError::NoConvergence`] after 60 sweeps.
pub fn svd(a: &Mat) -> Result<Svd, DenseError> {
    if a.rows < a.cols {
        let t = svd(&a.transpose())?;
        return Ok(Svd {
            u: t.v,
            s: t.s,
            v: t.u,
        });
    }
    let (m, n) = (a.rows, a.cols);
    let mut u = a.clone();
    let mut v = Mat::identity(n);
    let mut converged = false;
    for _ in 0..60 {
        let mut rotated = false;
        for p in 0..n {
            for q in p + 1..n {
                let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
                for i in 0..m {
                    let (x, y) = (u.at(i, p), u.at(i, q));
                    alpha += x * x;
                    beta += y * y;
                    gamma += x * y;
                }
                if gamma == 0.0 || gamma.abs() <= 1e-15 * (alpha * beta).sqrt() {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = if zeta == 0.0 {
                    1.0
                } else {
                    zeta.signum() / (zeta.abs() + zeta.mul_add(zeta, 1.0).sqrt())
                };
                let c = 1.0 / t.mul_add(t, 1.0).sqrt();
                let s = c * t;
                for i in 0..m {
                    let (x, y) = (u.at(i, p), u.at(i, q));
                    u.set(i, p, c * x - s * y);
                    u.set(i, q, s * x + c * y);
                }
                for i in 0..n {
                    let (x, y) = (v.at(i, p), v.at(i, q));
                    v.set(i, p, c * x - s * y);
                    v.set(i, q, s * x + c * y);
                }
            }
        }
        if !rotated {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(DenseError::NoConvergence);
    }
    let mut sv: Vec<f64> = (0..n)
        .map(|j| (0..m).map(|i| u.at(i, j) * u.at(i, j)).sum::<f64>().sqrt())
        .collect();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&x, &y| sv[y].total_cmp(&sv[x]));
    let mut uo = Mat::zeros(m, n);
    let mut vo = Mat::zeros(n, n);
    for (c, &j) in idx.iter().enumerate() {
        for i in 0..m {
            uo.set(
                i,
                c,
                if sv[j] > 0.0 {
                    u.at(i, j) / sv[j]
                } else {
                    0.0
                },
            );
        }
        for i in 0..n {
            vo.set(i, c, v.at(i, j));
        }
    }
    sv = idx.iter().map(|&j| sv[j]).collect();
    Ok(Svd { u: uo, s: sv, v: vo })
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
/// As [`svd`].
pub fn lstsq_svd(
    a: &Mat,
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    Ok(pinv(a)?.matvec(b))
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
    if a.rows != a.cols || b.rows != b.cols || a.rows != b.rows {
        return Err(DenseError::Shape);
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
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Vec::new());
    }
    // 1-based working copy.
    let mut h = vec![vec![0.0; n + 2]; n + 2];
    for i in 0..n {
        for j in 0..n {
            h[i + 1][j + 1] = a.at(i, j);
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
                    c += h[j][i].abs();
                    r += h[i][j].abs();
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
                        h[i][j] /= f;
                    }
                    for j in 1..=n {
                        h[j][i] *= f;
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
            if h[j][m - 1].abs() > x.abs() {
                x = h[j][m - 1];
                i = j;
            }
        }
        if i != m {
            for j in m - 1..=n {
                let t = h[i][j];
                h[i][j] = h[m][j];
                h[m][j] = t;
            }
            for j in 1..=n {
                h[j].swap(i, m);
            }
        }
        if x != 0.0 {
            for i in m + 1..=n {
                let mut y = h[i][m - 1];
                if y != 0.0 {
                    y /= x;
                    h[i][m - 1] = y;
                    for j in m..=n {
                        h[i][j] -= y * h[m][j];
                    }
                    for j in 1..=n {
                        h[j][m] += y * h[j][i];
                    }
                }
            }
        }
    }
    for i in 3..=n {
        for j in 1..i - 1 {
            h[i][j] = 0.0;
        }
    }
    let mut anorm = 0.0;
    for i in 1..=n {
        for j in i.saturating_sub(1).max(1)..=n {
            anorm += h[i][j].abs();
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
                s = h[l - 1][l - 1].abs() + h[l][l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if crate::kernels::exact_eq(h[l][l - 1].abs() + s, s) {
                    h[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            x = h[nn][nn];
            if l == nn {
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            y = h[nn - 1][nn - 1];
            w = h[nn][nn - 1] * h[nn - 1][nn];
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
                    h[i][i] -= x;
                }
                s = h[nn][nn - 1].abs() + h[nn - 1][nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;
            let mut m = nn - 2;
            loop {
                z = h[m][m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[m + 1][m] + h[m][m + 1];
                q = h[m + 1][m + 1] - z - r - s;
                r = h[m + 2][m + 1];
                s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == l {
                    break;
                }
                let u = h[m][m - 1].abs() * (q.abs() + r.abs());
                let v = p.abs() * (h[m - 1][m - 1].abs() + z.abs() + h[m + 1][m + 1].abs());
                if crate::kernels::exact_eq(u + v, v) {
                    break;
                }
                m -= 1;
            }
            for i in m + 2..=nn {
                h[i][i - 2] = 0.0;
                if i != m + 2 {
                    h[i][i - 3] = 0.0;
                }
            }
            for k in m..nn {
                if k != m {
                    p = h[k][k - 1];
                    q = h[k + 1][k - 1];
                    r = if k == nn - 1 {
                        0.0
                    } else {
                        h[k + 2][k - 1]
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
                            h[k][k - 1] = -h[k][k - 1];
                        }
                    } else {
                        h[k][k - 1] = -s * x;
                    }
                    p += s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q /= p;
                    r /= p;
                    for j in k..=nn {
                        p = h[k][j] + q * h[k + 1][j];
                        if k != nn - 1 {
                            p += r * h[k + 2][j];
                            h[k + 2][j] -= p * z;
                        }
                        h[k + 1][j] -= p * y;
                        h[k][j] -= p * x;
                    }
                    let mmin = if nn < k + 3 { nn } else { k + 3 };
                    for i in l..=mmin {
                        p = x * h[i][k] + y * h[i][k + 1];
                        if k != nn - 1 {
                            p += z * h[i][k + 2];
                            h[i][k + 2] -= p * r;
                        }
                        h[i][k + 1] -= p * q;
                        h[i][k] -= p;
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
    if a.rows != a.cols {
        return Err(DenseError::Shape);
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
    for k in 0..n - 2 {
        let mut norm = 0.0;
        for i in k + 1..n {
            norm += h.at(i, k) * h.at(i, k);
        }
        let norm = norm.sqrt();
        if norm == 0.0 {
            continue;
        }
        let x0 = h.at(k + 1, k);
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        let mut v = Vec::with_capacity(n - k - 1);
        v.push(x0 - alpha);
        for i in k + 2..n {
            v.push(h.at(i, k));
        }
        let vn2: f64 = v.iter().map(|e| e * e).sum();
        if vn2 == 0.0 {
            continue;
        }
        let beta = 2.0 / vn2;
        for j in k..n {
            let s: f64 = v
                .iter()
                .enumerate()
                .map(|(t, vt)| vt * h.at(k + 1 + t, j))
                .sum::<f64>()
                * beta;
            for (t, vt) in v.iter().enumerate() {
                h.set(k + 1 + t, j, h.at(k + 1 + t, j) - s * vt);
            }
        }
        for i in 0..n {
            let s: f64 = v
                .iter()
                .enumerate()
                .map(|(t, vt)| h.at(i, k + 1 + t) * vt)
                .sum::<f64>()
                * beta;
            for (t, vt) in v.iter().enumerate() {
                h.set(i, k + 1 + t, h.at(i, k + 1 + t) - s * vt);
            }
        }
        for i in 0..n {
            let s: f64 = v
                .iter()
                .enumerate()
                .map(|(t, vt)| q.at(i, k + 1 + t) * vt)
                .sum::<f64>()
                * beta;
            for (t, vt) in v.iter().enumerate() {
                q.set(i, k + 1 + t, q.at(i, k + 1 + t) - s * vt);
            }
        }
        h.set(k + 1, k, alpha);
        for i in k + 2..n {
            h.set(i, k, 0.0);
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
}
