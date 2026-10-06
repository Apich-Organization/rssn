//! # Dense linear algebra on row-major `f64` matrices
//!
//! A small self-contained dense toolkit: LU with partial pivoting,
//! Cholesky, Householder QR, least squares, the Jacobi eigenvalue algorithm
//! for symmetric matrices, a one-sided Jacobi SVD (singular values,
//! pseudo-inverse, rank, 2-norm condition number) and the eigenvalues of a
//! general real matrix by Hessenberg reduction and the shifted QR
//! algorithm.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use std::fmt;

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
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Shape => write!(f, "incompatible shapes"),
            Self::Singular => write!(f, "matrix is singular or not positive definite"),
            Self::NoConvergence => write!(f, "iteration did not converge"),
        }
    }
}

impl std::error::Error for DenseError {}

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
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![0.0; rows * cols] }
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
        Ok(Self { rows: rows.len(), cols: c, data: rows.concat() })
    }

    /// Entry `(i, j)`.
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }

    /// Sets entry `(i, j)`.
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
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
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        (0..self.rows)
            .map(|i| {
                self.data[i * self.cols..(i + 1) * self.cols]
                    .iter()
                    .zip(x)
                    .map(|(a, b)| a * b)
                    .sum()
            })
            .collect()
    }

    /// Matrix product.
    ///
    /// # Errors
    /// Returns [`DenseError::Shape`] on mismatched inner dimensions.
    pub fn matmul(&self, other: &Self) -> Result<Self, DenseError> {
        if self.cols != other.rows {
            return Err(DenseError::Shape);
        }
        let mut out = Self::zeros(self.rows, other.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let a = self.data[i * self.cols + k];
                if a != 0.0 {
                    for j in 0..other.cols {
                        out.data[i * other.cols + j] += a * other.data[k * other.cols + j];
                    }
                }
            }
        }
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
    lu: Mat,
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
    Ok(Lu { lu, piv, sign })
}

impl Lu {
    /// Solves `A x = b`.
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.lu.rows;
        let mut x: Vec<f64> = self.piv.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            for j in 0..i {
                x[i] -= self.lu.at(i, j) * x[j];
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                x[i] -= self.lu.at(i, j) * x[j];
            }
            x[i] /= self.lu.at(i, i);
        }
        x
    }

    /// Determinant of the factored matrix.
    #[must_use]
    pub fn det(&self) -> f64 {
        (0..self.lu.rows).map(|i| self.lu.at(i, i)).product::<f64>() * self.sign
    }

    /// Inverse of the factored matrix.
    #[must_use]
    pub fn inverse(&self) -> Mat {
        let n = self.lu.rows;
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
pub fn solve(a: &Mat, b: &[f64]) -> Result<Vec<f64>, DenseError> {
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
pub fn cholesky_solve(a: &Mat, b: &[f64]) -> Result<Vec<f64>, DenseError> {
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
        let alpha = if r.at(k, k) > 0.0 { -norm } else { norm };
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
pub fn lstsq(a: &Mat, b: &[f64]) -> Result<Vec<f64>, DenseError> {
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
        return Ok(Svd { u: t.v, s: t.s, v: t.u });
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
            uo.set(i, c, if sv[j] > 0.0 { u.at(i, j) / sv[j] } else { 0.0 });
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
    let (hi, lo) = (s.first().copied().unwrap_or(0.0), s.last().copied().unwrap_or(0.0));
    Ok(if lo == 0.0 { f64::INFINITY } else { hi / lo })
}

/// Numerical rank: singular values above `tol` (or the default
/// `max(m, n) * eps * sigma_max` when `tol` is `None`).
///
/// # Errors
/// As [`svd`].
pub fn rank(a: &Mat, tol: Option<f64>) -> Result<usize, DenseError> {
    let s = svd(a)?.s;
    let smax = s.first().copied().unwrap_or(0.0);
    let t = tol.unwrap_or(a.rows.max(a.cols) as f64 * f64::EPSILON * smax);
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
pub fn lstsq_svd(a: &Mat, b: &[f64]) -> Result<Vec<f64>, DenseError> {
    Ok(pinv(a)?.matvec(b))
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
                let t = h[j][i];
                h[j][i] = h[j][m];
                h[j][m] = t;
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
                if h[l][l - 1].abs() + s == s {
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
                if u + v == v {
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
                    r = if k != nn - 1 { h[k + 2][k - 1] } else { 0.0 };
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

/// Spectral radius estimate through [`eigenvalues`].
///
/// # Errors
/// As [`eigenvalues`].
pub fn spectral_radius(a: &Mat) -> Result<f64, DenseError> {
    Ok(eigenvalues(a)?.iter().map(|&(re, im)| re.hypot(im)).fold(0.0, f64::max))
}
