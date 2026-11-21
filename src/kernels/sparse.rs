//! # Numerical Sparse Matrix Operations
//!
//! This module provides utilities for working with sparse matrices. The storage
//! is [`faer::sparse::SparseColMat`], which is faer's compressed sparse type:
//! **compressed sparse column** (CSC), not CSR. It includes functions for
//! building a matrix from triplets, performing sparse matrix-vector
//! multiplication, and solving sparse linear systems using iterative methods
//! like Conjugate Gradient.
//!
//! ## Row/column convention
//!
//! CSC stores the row indices and values of each column `j` in matching slots.
//! In compressed storage these slots span `col_ptr[j]..col_ptr[j + 1]`.
//! Uncompressed storage can reserve additional capacity in each column, so
//! the matrix-vector kernels use faer's range of live entries. The ordinary
//! product accumulates `y[row] += val * x[col]`, and the transposed product
//! accumulates `y[col] += val * x[row]`.

use faer::sparse::SparseColMat;
use faer::sparse::SparseRowMat;
use faer::sparse::SymbolicSparseRowMat;
use faer::sparse::Triplet;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::ArrayD;
use serde::Deserialize;
use serde::Serialize;
use std::ops::Range;

use crate::kernels::matrix::Matrix;

/// Alias for a dynamic-dimensional array of f64.
pub type Array = ArrayD<f64>;

/// The sparse matrix type used throughout this module.
///
/// This is faer's compressed sparse column matrix with `usize` indices. The
/// constructors keep the historical `csr_*`/`to_csr` names (see
/// [`csr_from_triplets`]), but the storage is CSC; no caller depends on the
/// physical layout, only on the logical `(row, column)` indexing.
pub type SparseMat = SparseColMat<usize, f64>;

/// Creates a sparse matrix from a list of triplets.
///
/// The historical name is kept for API stability; faer stores the result in CSC
/// form regardless of the name. Duplicate `(row, col)` entries are summed, as
/// they were before.
///
/// # Arguments
/// * `rows` - The number of rows in the matrix.
/// * `cols` - The number of columns in the matrix.
/// * `triplets` - A slice of `(row_index, col_index, value)` tuples representing the non-zero entries.
///
/// # Returns
/// A [`SparseMat`] representing the sparse matrix.
///
/// # Panics
/// Panics if any triplet index is out of bounds for the given shape.
#[must_use]
pub fn csr_from_triplets(
    rows: usize,
    cols: usize,
    triplets: &[(usize, usize, f64)],
) -> SparseMat {
    let entries: Vec<Triplet<usize, usize, f64>> = triplets
        .iter()
        .map(|&(r, c, v)| Triplet::new(r, c, v))
        .collect();

    SparseMat::try_new_from_triplets(rows, cols, &entries)
        .unwrap_or_else(|e| panic!("invalid sparse triplets for a {rows}x{cols} matrix: {e}"))
}

/// Performs sparse matrix-vector multiplication for a sparse matrix and a standard `Vec`.
///
/// The matrix is stored in CSC, so the loop walks columns: every stored entry
/// contributes `val * vector[column]` to `result[row]`. This is the same
/// arithmetic as `y[i] = sum_j A[i][j] * x[j]`, but it never needs the transpose
/// of the matrix.
///
/// # Arguments
/// * `matrix` - The sparse matrix.
/// * `vector` - The dense vector.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` representing the result of the multiplication.
///
/// # Errors
/// Returns an error if matrix and vector dimensions are not compatible.
pub fn sp_mat_vec_mul(
    matrix: &SparseMat,
    vector: &[f64],
) -> Result<Vec<f64>, String> {
    if matrix.ncols() != vector.len() {
        return Err(MATVEC_SHAPE_ERROR.into());
    }
    let mut result = vec![0.0; matrix.nrows()];
    accumulate_matvec(matrix, vector, &mut result);
    Ok(result)
}

/// Computes `A x` into caller-owned storage, replacing every output element.
///
/// # Errors
/// The input needs one value per column and the output one per row.
/// On a shape error the output is untouched.
pub fn sp_mat_vec_mul_into(
    matrix: &SparseMat,
    vector: &[f64],
    result: &mut [f64],
) -> Result<(), String> {
    if matrix.ncols() != vector.len() || matrix.nrows() != result.len() {
        return Err(MATVEC_SHAPE_ERROR.into());
    }
    result.fill(0.0);
    accumulate_matvec(matrix, vector, result);
    Ok(())
}

/// Computes `Aᵀ x` by taking a dot product over each original column.
///
/// Each dot product uses the stored order of the live entries.
///
/// # Errors
/// The vector length differs from the number of rows of `A`.
pub fn sp_mat_transpose_vec_mul(
    matrix: &SparseMat,
    vector: &[f64],
) -> Result<Vec<f64>, String> {
    if matrix.nrows() != vector.len() {
        return Err(MATVEC_SHAPE_ERROR.into());
    }
    Ok(transpose_products(matrix, vector).collect())
}

/// Computes `Aᵀ x` into caller-owned storage, replacing every output element.
///
/// # Errors
/// The input needs one value per row and the output one per column.
/// On a shape error the output is untouched.
pub fn sp_mat_transpose_vec_mul_into(
    matrix: &SparseMat,
    vector: &[f64],
    result: &mut [f64],
) -> Result<(), String> {
    if matrix.nrows() != vector.len() || matrix.ncols() != result.len() {
        return Err(MATVEC_SHAPE_ERROR.into());
    }
    for (slot, value) in result.iter_mut().zip(transpose_products(matrix, vector)) {
        *slot = value;
    }
    Ok(())
}

const MATVEC_SHAPE_ERROR: &str =
    "Matrix and vector dimensions are not compatible for multiplication.";

fn live_columns(matrix: &SparseMat) -> impl Iterator<Item = (&[usize], &[f64])> {
    entries_at(matrix, live_column_spans(matrix))
}

fn live_column_spans(matrix: &SparseMat) -> impl Iterator<Item = Range<usize>> {
    let symbolic = matrix.symbolic();
    (0..matrix.ncols()).map(move |column| symbolic.col_range(column))
}

fn entries_at<'a>(
    matrix: &'a SparseMat,
    spans: impl Iterator<Item = Range<usize>> + 'a,
) -> impl Iterator<Item = (&'a [usize], &'a [f64])> + 'a {
    let rows = matrix.row_idx();
    let values = matrix.val();
    spans.map(move |span| (&rows[span.clone()], &values[span]))
}

fn accumulate_matvec(
    matrix: &SparseMat,
    vector: &[f64],
    result: &mut [f64],
) {
    for ((rows, values), &x) in live_columns(matrix).zip(vector) {
        for (&row, &value) in rows.iter().zip(values) {
            result[row] += value * x;
        }
    }
}

fn transpose_products<'a>(
    matrix: &'a SparseMat,
    vector: &'a [f64],
) -> impl Iterator<Item = f64> + 'a {
    live_columns(matrix).map(|(rows, values)| {
        rows.iter()
            .zip(values)
            .fold(0.0, |sum, (&row, &value)| sum + value * vector[row])
    })
}

/// Converts a dense `ndarray::Array` to a sparse matrix.
///
/// Zero entries are dropped, so the result holds one stored entry per non-zero.
///
/// # Arguments
/// * `arr` - The dense `ndarray::Array` to convert.
///
/// # Returns
/// A [`SparseMat`] representing the sparse matrix.
///
/// # Panics
/// Panics if the input array is not 2D.
#[must_use]
pub fn to_csr(arr: &Array) -> SparseMat {
    assert_eq!(
        arr.ndim(),
        2,
        "Input array must be 2D for \
         CSR conversion."
    );

    let rows = arr.shape()[0];

    let cols = arr.shape()[1];

    let mut triplets = Vec::new();

    for (i, row) in arr.outer_iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if val != 0.0 {
                triplets.push(Triplet::new(i, j, val));
            }
        }
    }

    SparseMat::try_new_from_triplets(rows, cols, &triplets)
        .unwrap_or_else(|e| panic!("invalid dense to sparse conversion: {e}"))
}

/// Converts a sparse matrix to a dense `ndarray::Array2`.
///
/// # Arguments
/// * `matrix` - The sparse matrix.
///
/// # Returns
/// An `ndarray::Array2<f64>` representing the dense matrix.
#[must_use]
pub fn to_dense(matrix: &SparseMat) -> Array2<f64> {
    let dense = matrix.to_dense();

    let rows = matrix.nrows();

    let cols = matrix.ncols();

    let mut out = Array2::<f64>::zeros((rows, cols));

    for j in 0..cols {
        let column = dense.col_as_slice(j);

        for (i, &val) in column.iter().enumerate() {
            out[[i, j]] = val;
        }
    }

    out
}

/// Computes the rank of a sparse matrix by converting to dense and performing RREF.
///
/// Note: This is inefficient for large sparse matrices and should only be used for small matrices.
///
/// # Arguments
/// * `matrix` - The sparse matrix.
///
/// # Returns
/// The rank of the matrix as a `usize`.
#[must_use]
pub fn rank(matrix: &SparseMat) -> usize {
    let dense_array2: Array2<f64> = to_dense(matrix);

    let rows = dense_array2.nrows();

    let cols = dense_array2.ncols();

    let mut dense_matrix = Matrix::new(rows, cols, dense_array2.into_raw_vec_and_offset().0);

    dense_matrix.rref().unwrap_or_default()
}

/// Transposes a sparse matrix.
///
/// The transpose is a column-major view of a row-major matrix, so faer converts
/// it without re-sorting the entries.
///
/// # Panics
/// Panics if the matrix cannot be converted to column-major format.
#[must_use]
pub fn transpose(matrix: &SparseMat) -> SparseMat {
    matrix
        .as_ref()
        .transpose()
        .to_col_major()
        .unwrap_or_else(|e| panic!("failed to transpose sparse matrix: {e}"))
}

/// Computes the trace of a square sparse matrix.
///
/// # Errors
/// Returns an error if the matrix is not square.
pub fn trace(matrix: &SparseMat) -> Result<f64, String> {
    if matrix.nrows() != matrix.ncols() {
        return Err("Matrix must be \
                    square to compute \
                    trace."
            .to_string());
    }

    let mut sum = 0.0;

    for i in 0..matrix.nrows() {
        if let Some(&val) = matrix.get(i, i) {
            sum += val;
        }
    }

    Ok(sum)
}

/// Checks if the sparse matrix is symmetric ($A = A^T$).
#[must_use]
pub fn is_symmetric(
    matrix: &SparseMat,
    epsilon: f64,
) -> bool {
    if matrix.nrows() != matrix.ncols() {
        return false;
    }

    for entry in matrix.triplet_iter() {
        let other = matrix.get(entry.col, entry.row).copied().unwrap_or(0.0);

        if (entry.val - other).abs() > epsilon {
            return false;
        }
    }

    true
}

/// Checks if the sparse matrix is diagonal.
#[must_use]
pub fn is_diagonal(matrix: &SparseMat) -> bool {
    for entry in matrix.triplet_iter() {
        if entry.row != entry.col && *entry.val != 0.0 {
            return false;
        }
    }

    true
}

/// Computes the Frobenius norm from the live CSC entries.
///
/// Extreme magnitudes are rescaled when the direct sum of squares leaves the
/// normal range. An empty matrix has norm zero.
#[must_use]
pub fn frobenius_norm(matrix: &SparseMat) -> f64 {
    let squares = fold_live_values(matrix, |sum, &value| sum + value * value);
    finish_frobenius(matrix, squares)
}

/// Computes the maximum absolute column sum from the live CSC entries.
///
/// An empty matrix has norm zero. Column sums containing `NaN` are skipped.
#[must_use]
pub fn l1_norm(matrix: &SparseMat) -> f64 {
    finish_l1(matrix, max_column_sum(matrix, |_, sum| sum))
}

/// Computes the maximum absolute row sum from the live CSC entries.
///
/// An empty matrix has norm zero. Row sums containing `NaN` are skipped.
#[must_use]
pub fn linf_norm(matrix: &SparseMat) -> f64 {
    if matrix.nrows() == 0 || matrix.ncols() == 0 || matrix.val().is_empty() {
        return 0.0;
    }
    let mut row_sums = vec![0.0; matrix.nrows()];
    for_each_live_chunk(matrix, |rows, values| {
        for (&row, &value) in rows.iter().zip(values) {
            row_sums[row] += value.abs();
        }
    });
    finish_row_sums(matrix, &mut row_sums)
}

/// Frobenius and induced norms of a real sparse matrix.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SparseNorms {
    /// Square root of the sum of squared entries.
    pub frobenius: f64,
    /// Maximum absolute column sum.
    pub l1: f64,
    /// Maximum absolute row sum.
    pub linf: f64,
}

/// Computes the Frobenius, induced 1- and induced infinity norms together.
///
/// Ordinary inputs use one traversal of the live entries. Extreme scales
/// may need additional passes. Empty matrices produce three zeros.
#[must_use]
pub fn matrix_norms(matrix: &SparseMat) -> SparseNorms {
    if matrix.nrows() == 0 || matrix.ncols() == 0 || matrix.val().is_empty() {
        return SparseNorms::default();
    }
    let mut row_sums = vec![0.0; matrix.nrows()];
    let mut squares = 0.0;
    let mut l1 = 0.0_f64;
    let populated = live_column_spans(matrix).filter(|span| !span.is_empty());
    for (rows, values) in entries_at(matrix, populated) {
        let mut col_sum = 0.0;
        for (&row, &value) in rows.iter().zip(values) {
            let magnitude = value.abs();
            squares += value * value;
            col_sum += magnitude;
            row_sums[row] += magnitude;
        }
        l1 = l1.max(col_sum);
    }
    SparseNorms {
        frobenius: finish_frobenius(matrix, squares),
        l1: finish_l1(matrix, l1),
        linf: finish_row_sums(matrix, &mut row_sums),
    }
}

fn for_each_live_chunk(
    matrix: &SparseMat,
    mut visit: impl FnMut(&[usize], &[f64]),
) {
    let pointers = matrix.col_ptr();
    // Compressed columns occupy one contiguous live span, even when the
    // backing arrays contain unused entries before or after it.
    if matrix.col_nnz().is_none()
        && let (Some(&start), Some(&end)) = (pointers.first(), pointers.last())
    {
        visit(&matrix.row_idx()[start..end], &matrix.val()[start..end]);
    } else {
        for (rows, values) in live_columns(matrix) {
            visit(rows, values);
        }
    }
}

fn fold_live_values(
    matrix: &SparseMat,
    fold: impl Fn(f64, &f64) -> f64,
) -> f64 {
    let mut result = 0.0;
    for_each_live_chunk(matrix, |_, values| {
        result = values.iter().fold(result, &fold);
    });
    result
}

fn finish_frobenius(
    matrix: &SparseMat,
    squares: f64,
) -> f64 {
    if squares.is_normal() || squares.is_nan() {
        return squares.sqrt();
    }
    let scale = fold_live_values(matrix, |largest, value| largest.max(value.abs()));
    if scale == 0.0 || scale.is_infinite() {
        return scale;
    }
    let scaled = fold_live_values(matrix, |sum, &value| {
        let relative = value / scale;
        sum + relative * relative
    });
    scale * scaled.sqrt()
}

const LARGE_TOTAL: f64 = 0.5 * f64::MAX;

fn max_column_sum(
    matrix: &SparseMat,
    finish: impl Fn(&[f64], f64) -> f64,
) -> f64 {
    let values = matrix.val();
    live_column_spans(matrix)
        .filter(|span| !span.is_empty())
        .fold(0.0, |largest, span| {
            let column = &values[span];
            let sum = column.iter().fold(0.0, |sum, value| sum + value.abs());
            largest.max(finish(column, sum))
        })
}

fn finish_l1(
    matrix: &SparseMat,
    largest: f64,
) -> f64 {
    if largest >= LARGE_TOTAL {
        max_column_sum(matrix, finish_column_sum)
    } else {
        largest
    }
}

// Keeping half of each positive term gives room for the compensated sum
// before rounding the result at the original scale.
#[derive(Clone, Copy, Default)]
struct HalfSum {
    sum: f64,
    correction: f64,
}

impl HalfSum {
    fn add(
        &mut self,
        value: f64,
    ) {
        let half = 0.5 * value.abs();
        let total = self.sum + half;
        if total.is_finite() {
            self.correction += if self.sum >= half {
                (self.sum - total) + half
            } else {
                (half - total) + self.sum
            };
        }
        self.sum = total;
    }

    fn total(self) -> f64 {
        2.0 * (self.sum + self.correction)
    }
}

fn finish_column_sum(
    values: &[f64],
    sum: f64,
) -> f64 {
    if sum >= LARGE_TOTAL {
        let mut recovered = HalfSum::default();
        for &value in values {
            recovered.add(value);
        }
        recovered.total()
    } else {
        sum
    }
}

fn finish_row_sums(
    matrix: &SparseMat,
    sums: &mut [f64],
) -> f64 {
    let maximum = |sums: &[f64]| sums.iter().copied().fold(0.0, f64::max);
    let largest = maximum(sums);
    if largest < LARGE_TOTAL {
        return largest;
    }
    let mut recovered = vec![HalfSum::default(); sums.len()];
    for_each_live_chunk(matrix, |rows, values| {
        for (&row, &value) in rows.iter().zip(values) {
            if sums[row] >= LARGE_TOTAL {
                recovered[row].add(value);
            }
        }
    });
    for (sum, recovered) in sums.iter_mut().zip(recovered) {
        if *sum >= LARGE_TOTAL {
            *sum = recovered.total();
        }
    }
    maximum(sums)
}

/// A serializable representation of a sparse matrix in CSR format.
///
/// The field layout is unchanged from the `sprs` days and stays genuinely
/// row-major: [`From<&SparseMat>`] converts the CSC storage to CSR before
/// copying, and [`SparseMatrixData::to_csmat`] converts it back, so serialized
/// data keeps the meaning its field names promise.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SparseMatrixData {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// CSR row pointer.
    pub indptr: Vec<usize>,
    /// CSR column indices.
    pub indices: Vec<usize>,
    /// Non-zero data values.
    pub data: Vec<f64>,
}

impl From<&SparseMat> for SparseMatrixData {
    fn from(mat: &SparseMat) -> Self {
        let row_major = mat
            .to_row_major()
            .unwrap_or_else(|e| panic!("failed to convert sparse matrix to CSR: {e}"));

        Self {
            rows: mat.nrows(),
            cols: mat.ncols(),
            indptr: row_major.row_ptr().to_vec(),
            indices: row_major.col_idx().to_vec(),
            data: row_major.val().to_vec(),
        }
    }
}

impl SparseMatrixData {
    /// Converts back to a [`SparseMat`] sparse matrix.
    ///
    /// # Panics
    /// Panics if `indptr`, `indices` and `data` are not a consistent CSR
    /// representation of a `rows` by `cols` matrix.
    #[must_use]
    pub fn to_csmat(&self) -> SparseMat {
        let symbolic = SymbolicSparseRowMat::<usize>::new_checked(
            self.rows,
            self.cols,
            self.indptr.clone(),
            None,
            self.indices.clone(),
        );

        let row_major = SparseRowMat::<usize, f64>::new(symbolic, self.data.clone());

        row_major
            .to_col_major()
            .unwrap_or_else(|e| panic!("failed to convert CSR data to a sparse matrix: {e}"))
    }
}

/// Solves a sparse linear system `Ax=b` using the Conjugate Gradient method.
///
/// This method is suitable for symmetric, positive-definite matrices. It is an iterative
/// algorithm that converges to the exact solution in at most `n` iterations (where `n` is
/// the matrix dimension) in exact arithmetic, but is typically stopped earlier based on a tolerance.
///
/// The tolerance is an absolute bound on the Euclidean norm of the residual
/// `b - A x`, checked before the first iteration and after every update. If the
/// iteration limit is reached without meeting it, the current iterate is
/// returned anyway. faer 0.24 does ship a conjugate gradient
/// (`faer::matrix_free::conjugate_gradient`), but it tests a *relative*
/// residual, reports `NonConvergence` as an error and requires a positive
/// definite operator, so using it would change all three behaviours; the loop is
/// kept, iterating over the faer storage.
///
/// # Arguments
/// * `a` - The sparse matrix `A` ([`SparseMat`]).
/// * `b` - The vector `b` (`Array1<f64>`).
/// * `x0` - An initial guess for the solution `x`.
/// * `max_iter` - The maximum number of iterations.
/// * `tolerance` - The desired tolerance for the residual norm.
///
/// # Returns
/// A `Result` containing the solution vector `x`, or an error string.
///
/// # Errors
/// Returns an error if the matrix and vector dimensions are incompatible.
pub fn solve_conjugate_gradient(
    a: &SparseMat,
    b: &Array1<f64>,
    x0: Option<&Array1<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> Result<Array1<f64>, String> {
    let n = a.ncols();

    if a.nrows() != n || b.len() != n {
        return Err("Matrix and \
                    vector dimensions \
                    are incompatible.\
                    "
        .to_string());
    }

    let mut x = x0.cloned().unwrap_or_else(|| Array1::zeros(n));

    let mut r = b - &apply(a, &x)?;

    let mut p = r.clone();

    let mut rs_old = r.dot(&r);

    if rs_old.sqrt() < tolerance {
        return Ok(x);
    }

    for _ in 0..max_iter {
        let ap = apply(a, &p)?;

        let alpha = rs_old / p.dot(&ap);

        x = &x + &(&p * alpha);

        r = &r - &(&ap * alpha);

        let rs_new = r.dot(&r);

        if rs_new.sqrt() < tolerance {
            break;
        }

        p = &r + &(&p * (rs_new / rs_old));

        rs_old = rs_new;
    }

    Ok(x)
}

/// Applies the sparse matrix to a dense vector.
fn apply(
    a: &SparseMat,
    x: &Array1<f64>,
) -> Result<Array1<f64>, String> {
    let dense = x
        .as_slice()
        .ok_or_else(|| "initial guess must be contiguous".to_string())?;

    Ok(Array1::from_vec(sp_mat_vec_mul(a, dense)?))
}

#[cfg(test)]
mod tests {

    use ndarray::array;

    use super::*;

    #[test]
    pub fn test_csr_from_triplets() {
        let triplets = vec![(0, 0, 1.0), (1, 2, 2.0), (2, 1, 3.0)];

        let mat = csr_from_triplets(3, 3, &triplets);

        assert_eq!(mat.nrows(), 3);

        assert_eq!(mat.ncols(), 3);

        assert_eq!(mat.compute_nnz(), 3);

        assert_eq!(mat.get(0, 0), Some(&1.0));

        assert_eq!(mat.get(1, 2), Some(&2.0));

        assert_eq!(mat.get(2, 1), Some(&3.0));

        assert_eq!(mat.get(0, 1), None);
    }

    #[test]
    pub fn test_sp_mat_vec_mul() {
        let triplets = vec![(0, 0, 1.0), (0, 2, 2.0), (2, 1, 3.0)];

        let mat = csr_from_triplets(3, 3, &triplets);

        let vec = vec![10.0, 20.0, 30.0];

        let result = sp_mat_vec_mul(&mat, &vec);

        match result {
            | Ok(res) => {
                assert_eq!(res, vec![70.0, 0.0, 60.0]);
            },
            | Err(e) => {
                panic!(
                    "sp_mat_vec_mul \
                     failed with: {e}"
                )
            },
        }
    }

    /// An asymmetric, non-square matrix: a row/column mix-up cannot pass this.
    #[test]
    pub fn test_sp_mat_vec_mul_asymmetric_non_square() {
        // [[0, 2, 3],
        //  [5, 0, 0]]
        let mat = csr_from_triplets(2, 3, &[(0, 1, 2.0), (0, 2, 3.0), (1, 0, 5.0)]);

        assert_eq!(to_dense(&mat), array![[0.0, 2.0, 3.0], [5.0, 0.0, 0.0]]);

        // 2*10 + 3*100 = 320, 5*1 = 5.
        assert_eq!(
            sp_mat_vec_mul(&mat, &[1.0, 10.0, 100.0]),
            Ok(vec![320.0, 5.0])
        );

        // The transpose is 3x2 and must multiply a 2-vector.
        let transposed = transpose(&mat);

        assert_eq!(
            to_dense(&transposed),
            array![[0.0, 5.0], [2.0, 0.0], [3.0, 0.0]]
        );

        assert_eq!(
            sp_mat_vec_mul(&transposed, &[1.0, 10.0]),
            Ok(vec![50.0, 2.0, 3.0])
        );
    }

    /// Duplicate triplets must sum, as `TriMat::to_csr` did.
    #[test]
    pub fn test_csr_from_triplets_sums_duplicates() {
        let mat = csr_from_triplets(2, 2, &[(0, 1, 1.5), (0, 1, 2.5)]);

        assert_eq!(mat.compute_nnz(), 1);

        assert_eq!(mat.get(0, 1), Some(&4.0));
    }

    #[test]
    pub fn test_to_csr() {
        let dense_arr = array![[1.0, 0.0, 2.0], [0.0, 0.0, 0.0], [3.0, 0.0, 4.0]].into_dyn();

        let csr_mat = to_csr(&dense_arr);

        assert_eq!(csr_mat.nrows(), 3);

        assert_eq!(csr_mat.ncols(), 3);

        assert_eq!(csr_mat.compute_nnz(), 4);

        assert_eq!(csr_mat.get(0, 0), Some(&1.0));

        assert_eq!(csr_mat.get(0, 2), Some(&2.0));

        assert_eq!(csr_mat.get(2, 0), Some(&3.0));

        assert_eq!(csr_mat.get(2, 2), Some(&4.0));

        assert_eq!(csr_mat.get(1, 1), None);
    }

    /// `SparseMatrixData` is documented as CSR, so the round trip must preserve
    /// the orientation even though the in-memory storage is CSC.
    #[test]
    pub fn test_sparse_matrix_data_is_row_major() {
        // [[0, 1, 0],
        //  [2, 0, 3]]
        let mat = csr_from_triplets(2, 3, &[(0, 1, 1.0), (1, 0, 2.0), (1, 2, 3.0)]);

        let data = SparseMatrixData::from(&mat);

        assert_eq!(data.rows, 2);

        assert_eq!(data.cols, 3);

        assert_eq!(data.indptr, vec![0, 1, 3]);

        assert_eq!(data.indices, vec![1, 0, 2]);

        assert_eq!(data.data, vec![1.0, 2.0, 3.0]);

        let back = data.to_csmat();

        assert_eq!(to_dense(&back), to_dense(&mat));
    }
}
