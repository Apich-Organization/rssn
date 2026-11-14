//! The `faer` backend for the `f64` linear-algebra kernels.
//!
//! `faer` is column-major and takes its own matrix type, so every entry point
//! copies the caller's row-major slice into a `Mat`. That copy is O(n²) against
//! an O(n³) kernel and is not the thing worth optimising.
//!
//! Everything here is `f64`. There is no generic element type and no strided
//! entry point: a non-contiguous operand cannot be expressed, so it cannot be
//! read with the wrong stride. A view is materialised by the caller before it
//! reaches this module.

use super::DenseError;
use faer::Accum;
use faer::Mat;
use faer::get_global_parallelism;
use faer::linalg::solvers::PartialPivLu;
use faer::linalg::solvers::Qr;
use faer::linalg::solvers::Solve;

/// A row-major slice as a `faer` matrix, checking the shape.
fn matrix(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<Mat<f64>, DenseError> {
    if rows.checked_mul(cols) != Some(a.len()) {
        return Err(DenseError::Shape);
    }
    if a.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let mut out = Mat::<f64>::zeros(rows, cols);
    for row in 0..rows {
        for column in 0..cols {
            out[(row, column)] = a[row * cols + column];
        }
    }
    Ok(out)
}

/// A `faer` matrix as a row-major `Vec`.
fn to_row_major(matrix: faer::MatRef<'_, f64>) -> Vec<f64> {
    let (rows, cols) = matrix.shape();
    let mut out = vec![0.0; rows * cols];
    for row in 0..rows {
        for column in 0..cols {
            out[row * cols + column] = matrix[(row, column)];
        }
    }
    out
}

/// `true` when the `U` factor has a zero on its diagonal, i.e. the matrix the
/// factorisation was taken of is exactly singular.
///
/// An exactly singular matrix is the case the old path reported
/// [`DenseError::Singular`] for, and a nearly singular one is not: the
/// answer is large, not absent, and refusing it would be a semantic change.
fn is_singular(
    u: faer::MatRef<'_, f64>,
    n: usize,
) -> bool {
    (0..n).any(|index| u[(index, index)] == 0.0)
}

/// The determinant of a square matrix.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn det(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<f64, DenseError> {
    if rows != cols {
        return Err(DenseError::Shape);
    }
    let matrix = matrix(rows, cols, a)?;
    if rows == 0 {
        return Ok(1.0);
    }
    Ok(matrix.as_ref().determinant())
}

/// The inverse of a square matrix.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn inv(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if rows != cols {
        return Err(DenseError::Shape);
    }
    let matrix = matrix(rows, cols, a)?;
    if rows == 0 {
        return Ok(Vec::new());
    }
    let lu = PartialPivLu::new(matrix.as_ref());
    if is_singular(lu.U(), rows) {
        return Err(DenseError::Singular);
    }
    let mut inverse = Mat::<f64>::identity(rows, rows);
    lu.solve_in_place(inverse.as_mut());
    Ok(to_row_major(inverse.as_ref()))
}

/// `A x = b` for a square `A`.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn solve(
    rows: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if rows != cols {
        return Err(DenseError::Shape);
    }
    let matrix = matrix(rows, cols, a)?;
    if b.len() != rows {
        return Err(DenseError::Shape);
    }
    if b.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    if rows == 0 {
        return Ok(Vec::new());
    }
    let lu = PartialPivLu::new(matrix.as_ref());
    if is_singular(lu.U(), rows) {
        return Err(DenseError::Singular);
    }
    let mut x = Mat::<f64>::zeros(rows, 1);
    for (row, value) in b.iter().enumerate() {
        x[(row, 0)] = *value;
    }
    lu.solve_in_place(x.as_mut());
    Ok((0..rows).map(|row| x[(row, 0)]).collect())
}

/// `P A = L U` with partial pivoting.
///
/// `faer` reports the row permutation as a pair of index arrays; the `P` this
/// crate promises is the permutation *matrix*, so it is materialised from the
/// forward array as `P[i][forward[i]] = 1`, which is `(P A)[i] = A[forward[i]]`.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
// The established decomposition API returns the three P, L, U factors together.
#[allow(clippy::type_complexity)]
pub fn lu(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), DenseError> {
    if rows != cols {
        return Err(DenseError::Shape);
    }
    let matrix = matrix(rows, cols, a)?;
    if rows == 0 {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    }
    let lu = PartialPivLu::new(matrix.as_ref());
    let (forward, _) = lu.P().arrays();
    let mut p = vec![0.0; rows * rows];
    for (row, &source) in forward.iter().enumerate() {
        p[row * rows + source] = 1.0;
    }
    Ok((p, to_row_major(lu.L()), to_row_major(lu.U())))
}

/// `A = Q R` with `Q` orthonormal and `rows × rows`.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn qr(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), DenseError> {
    let q_len = rows.checked_mul(rows).ok_or(DenseError::Shape)?;
    let matrix = matrix(rows, cols, a)?;
    if rows == 0 || cols == 0 {
        let mut q = vec![0.0; q_len];
        for i in 0..rows {
            q[i * rows + i] = 1.0;
        }
        return Ok((q, vec![0.0; rows * cols]));
    }
    let qr = Qr::new(matrix.as_ref());
    let q = qr.compute_Q();
    Ok((to_row_major(q.as_ref()), to_row_major(qr.R())))
}

/// `A B` for `A` of shape `rows × inner` and `B` of `inner × cols`.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn matmul(
    rows: usize,
    inner: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    let output_len = rows.checked_mul(cols).ok_or(DenseError::Shape)?;
    let lhs = matrix(rows, inner, a)?;
    let rhs = matrix(inner, cols, b)?;
    if rows == 0 || inner == 0 || cols == 0 {
        return Ok(vec![0.0; output_len]);
    }
    let mut out = Mat::<f64>::zeros(rows, cols);
    faer::linalg::matmul::matmul(
        out.as_mut(),
        Accum::Replace,
        lhs.as_ref(),
        rhs.as_ref(),
        1.0,
        get_global_parallelism(),
    );
    Ok(to_row_major(out.as_ref()))
}

/// Minimum-norm least-squares solve using the shared thin SVD and an absolute
/// singular-value threshold of `1e-12`, retained for host API compatibility.
///
/// # Errors
/// Inconsistent shapes, non-finite input, or failed SVD convergence.
pub fn lstsq(
    rows: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if rows.checked_mul(cols) != Some(a.len()) || b.len() != rows {
        return Err(DenseError::Shape);
    }
    if a.iter().chain(b).any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let a = super::Mat {
        rows,
        cols,
        data: a.to_vec(),
    };
    let factors = super::svd(&a)?;
    let mut solution = vec![0.0; cols];
    for k in 0..factors.s.len() {
        if factors.s[k] <= 1e-12 {
            continue;
        }
        let projection = (0..rows).map(|i| factors.u.at(i, k) * b[i]).sum::<f64>() / factors.s[k];
        for j in 0..cols {
            solution[j] += factors.v.at(j, k) * projection;
        }
    }
    Ok(solution)
}
