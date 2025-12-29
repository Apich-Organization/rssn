//! # Numerical Sparse Matrix Operations
//!
//! This module provides utilities for working with sparse matrices, particularly
//! in the Compressed Sparse Row (CSR) format. It includes functions for creating
//! CSR matrices, performing sparse matrix-vector multiplication, and solving
//! sparse linear systems using iterative methods like Conjugate Gradient.

use ndarray::ArrayD;
use sprs_rssn::CsMat;
use sprs_rssn::TriMat;

/// Alias for a dynamic-dimensional array of f64.

pub type Array = ArrayD<f64>;

/// Creates a new CSR matrix from a triplet matrix.
///
/// # Arguments
/// * `rows` - The number of rows in the matrix.
/// * `cols` - The number of columns in the matrix.
/// * `triplets` - A slice of `(row_index, col_index, value)` tuples representing the non-zero entries.
///
/// # Returns
/// A `CsMat<f64>` representing the sparse matrix.
#[must_use]

pub fn csr_from_triplets(
    rows: usize,
    cols: usize,
    triplets: &[(usize, usize, f64)],
) -> CsMat<f64> {

    let mut mat =
        TriMat::new((rows, cols));

    for &(r, c, v) in triplets {

        mat.add_triplet(r, c, v);
    }

    mat.to_csr()
}

/// Performs sparse matrix-vector multiplication for a CSR matrix and a standard `Vec`.
///
/// # Arguments
/// * `matrix` - The sparse matrix in CSR format.
/// * `vector` - The dense vector.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` representing the result of the multiplication.
///
/// # Errors
/// Returns an error if matrix and vector dimensions are not compatible.

pub fn sp_mat_vec_mul(
    matrix: &CsMat<f64>,
    vector: &[f64],
) -> Result<Vec<f64>, String> {

    if matrix.cols() != vector.len() {

        return Err("Matrix and \
                    vector dimensions \
                    are not compatible \
                    for multiplication.\
                    "
        .to_string());
    }

    let mut result =
        vec![0.0; matrix.rows()];

    for (i, row) in matrix
        .outer_iterator()
        .enumerate()
    {

        let mut row_sum = 0.0;

        for (j, &val) in row.iter() {

            row_sum += val * vector[j];
        }

        result[i] = row_sum;
    }

    Ok(result)
}

/// Converts a dense `ndarray::Array` to a Compressed Sparse Row (CSR) matrix.
///
/// This implementation directly constructs the CSR vectors for efficiency.
///
/// # Arguments
/// * `arr` - The dense `ndarray::Array` to convert.
///
/// # Returns
/// A `CsMat<f64>` representing the sparse matrix.
///
/// # Panics
/// Panics if the input array is not 2D.
#[must_use]

pub fn to_csr(
    arr: &Array
) -> CsMat<f64> {

    assert_eq!(
        arr.ndim(),
        2,
        "Input array must be 2D for \
         CSR conversion."
    );

    let rows = arr.shape()[0];

    let cols = arr.shape()[1];

    let mut indptr =
        Vec::with_capacity(rows + 1);

    let mut indices = Vec::new();

    let mut data = Vec::new();

    indptr.push(0);

    let mut non_zeros = 0;

    for row in arr.outer_iter() {

        for (j, &val) in row
            .iter()
            .enumerate()
        {

            if val != 0.0 {

                non_zeros += 1;

                indices.push(j);

                data.push(val);
            }
        }

        indptr.push(non_zeros);
    }

    CsMat::new(
        (rows, cols),
        indptr,
        indices,
        data,
    )
}

/// Converts a Compressed Sparse Row (CSR) matrix to a dense `ndarray::Array`.
///
/// # Arguments
/// * `matrix` - The sparse matrix in CSR format.
///
/// # Returns
/// An `ndarray::Array2<f64>` representing the dense matrix.
#[must_use]

pub fn to_dense(
    matrix: &CsMat<f64>
) -> Array2<f64> {

    matrix.to_dense()
}

use ndarray::Array2;

use crate::numerical::matrix::Matrix;

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

pub fn rank(
    matrix: &CsMat<f64>
) -> usize {

    let dense_array2: Array2<f64> =
        matrix.to_dense();

    let rows = dense_array2.nrows();

    let cols = dense_array2.ncols();

    let mut dense_matrix = Matrix::new(
        rows,
        cols,
        dense_array2
            .into_raw_vec_and_offset()
            .0,
    );

    dense_matrix
        .rref()
        .unwrap_or_default()
}

/// Transposes a sparse matrix.
#[must_use]

pub fn transpose(
    matrix: &CsMat<f64>
) -> CsMat<f64> {

    matrix
        .clone()
        .transpose_into()
}

/// Computes the trace of a square sparse matrix.
///
/// # Errors
/// Returns an error if the matrix is not square.

pub fn trace(
    matrix: &CsMat<f64>
) -> Result<f64, String> {

    if matrix.rows() != matrix.cols() {

        return Err("Matrix must be \
                    square to compute \
                    trace."
            .to_string());
    }

    let mut sum = 0.0;

    for i in 0 .. matrix.rows() {

        if let Some(&val) =
            matrix.get(i, i)
        {

            sum += val;
        }
    }

    Ok(sum)
}

/// Checks if the sparse matrix is symmetric ($A = A^T$).
#[must_use]

pub fn is_symmetric(
    matrix: &CsMat<f64>,
    epsilon: f64,
) -> bool {

    if matrix.rows() != matrix.cols() {

        return false;
    }

    for (val, (r, c)) in matrix {

        let other = matrix
            .get(c, r)
            .copied()
            .unwrap_or(0.0);

        if (val - other).abs() > epsilon
        {

            return false;
        }
    }

    true
}

/// Checks if the sparse matrix is diagonal.
#[must_use]

pub fn is_diagonal(
    matrix: &CsMat<f64>
) -> bool {

    for (&val, (r, c)) in matrix {

        if r != c && val != 0.0 {

            return false;
        }
    }

    true
}

/// Computes the Frobenius norm of a sparse matrix.
#[must_use]

pub fn frobenius_norm(
    matrix: &CsMat<f64>
) -> f64 {

    let mut sum = 0.0;

    for &val in matrix.data() {

        sum += val * val;
    }

    sum.sqrt()
}

/// Computes the L1 norm of a sparse matrix (max column sum).
#[must_use]

pub fn l1_norm(
    matrix: &CsMat<f64>
) -> f64 {

    let mut col_sums =
        vec![0.0; matrix.cols()];

    for (val, (_, c)) in matrix {

        col_sums[c] += val.abs();
    }

    col_sums
        .into_iter()
        .fold(0.0, |max, s| {
            if s > max {

                s
            } else {

                max
            }
        })
}

/// Computes the Linf norm of a sparse matrix (max row sum).
#[must_use]

pub fn linf_norm(
    matrix: &CsMat<f64>
) -> f64 {

    let mut max_sum = 0.0;

    for row in matrix.outer_iterator() {

        let mut row_sum = 0.0;

        for (_, &val) in row.iter() {

            row_sum += val.abs();
        }

        if row_sum > max_sum {

            max_sum = row_sum;
        }
    }

    max_sum
}

use serde::Deserialize;
use serde::Serialize;

/// A serializable representation of a sparse matrix in CSR format.
#[derive(
    Serialize, Deserialize, Debug, Clone,
)]

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

impl From<&CsMat<f64>>
    for SparseMatrixData
{
    fn from(mat: &CsMat<f64>) -> Self {

        Self {
            rows: mat.rows(),
            cols: mat.cols(),
            indptr: mat
                .indptr()
                .as_slice()
                .unwrap_or(&[])
                .to_vec(),
            indices: mat
                .indices()
                .to_vec(),
            data: mat.data().to_vec(),
        }
    }
}

impl SparseMatrixData {
    /// Converts back to a `CsMat` sparse matrix.
    #[must_use]

    pub fn to_csmat(
        &self
    ) -> CsMat<f64> {

        CsMat::new(
            (self.rows, self.cols),
            self.indptr.clone(),
            self.indices.clone(),
            self.data.clone(),
        )
    }
}

#[cfg(test)]

mod tests {

    use ndarray::array;

    use super::*;

    #[test]

    pub(crate) fn test_csr_from_triplets(
    ) {

        let triplets = vec![
            (0, 0, 1.0),
            (1, 2, 2.0),
            (2, 1, 3.0),
        ];

        let mat = csr_from_triplets(
            3,
            3,
            &triplets,
        );

        assert_eq!(mat.rows(), 3);

        assert_eq!(mat.cols(), 3);

        assert_eq!(mat.nnz(), 3);

        assert_eq!(
            mat.get(0, 0),
            Some(&1.0)
        );

        assert_eq!(
            mat.get(1, 2),
            Some(&2.0)
        );

        assert_eq!(
            mat.get(2, 1),
            Some(&3.0)
        );

        assert_eq!(mat.get(0, 1), None);
    }

    #[test]

    pub(crate) fn test_sp_mat_vec_mul()
    {

        let triplets = vec![
            (0, 0, 1.0),
            (0, 2, 2.0),
            (2, 1, 3.0),
        ];

        let mat = csr_from_triplets(
            3,
            3,
            &triplets,
        );

        let vec =
            vec![10.0, 20.0, 30.0];

        let result =
            sp_mat_vec_mul(&mat, &vec);

        match result {
            | Ok(res) => {

                assert_eq!(
                    res,
                    vec![
                        70.0, 0.0, 60.0
                    ]
                )
            },
            | Err(e) => {

                panic!(
                    "sp_mat_vec_mul \
                     failed with: {}",
                    e
                )
            },
        }
    }

    #[test]

    pub(crate) fn test_to_csr() {

        let dense_arr = array![
            [1.0, 0.0, 2.0],
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 4.0]
        ]
        .into_dyn();

        let csr_mat =
            to_csr(&dense_arr);

        assert_eq!(csr_mat.rows(), 3);

        assert_eq!(csr_mat.cols(), 3);

        assert_eq!(csr_mat.nnz(), 4);

        assert_eq!(
            csr_mat.get(0, 0),
            Some(&1.0)
        );

        assert_eq!(
            csr_mat.get(0, 2),
            Some(&2.0)
        );

        assert_eq!(
            csr_mat.get(2, 0),
            Some(&3.0)
        );

        assert_eq!(
            csr_mat.get(2, 2),
            Some(&4.0)
        );

        assert_eq!(
            csr_mat.get(1, 1),
            None
        );
    }
}

use ndarray::Array1;

/// Solves a sparse linear system `Ax=b` using the Conjugate Gradient method.
///
/// This method is suitable for symmetric, positive-definite matrices. It is an iterative
/// algorithm that converges to the exact solution in at most `n` iterations (where `n` is
/// the matrix dimension) in exact arithmetic, but is typically stopped earlier based on a tolerance.
///
/// # Arguments
/// * `a` - The sparse matrix `A` (`CsMat<f64>`).
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
    a: &CsMat<f64>,
    b: &Array1<f64>,
    x0: Option<&Array1<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> Result<Array1<f64>, String> {

    let n = a.cols();

    if a.rows() != n || b.len() != n {

        return Err("Matrix and \
                    vector dimensions \
                    are incompatible.\
                    "
        .to_string());
    }

    let mut x = x0
        .cloned()
        .unwrap_or_else(|| {

            Array1::zeros(n)
        });

    let mut r = b - &(a * &x);

    let mut p = r.clone();

    let mut rs_old = r.dot(&r);

    if rs_old.sqrt() < tolerance {

        return Ok(x);
    }

    for _ in 0 .. max_iter {

        let ap = a * &p;

        let alpha = rs_old / p.dot(&ap);

        x = &x + &(&p * alpha);

        r = &r - &(&ap * alpha);

        let rs_new = r.dot(&r);

        if rs_new.sqrt() < tolerance {

            break;
        }

        p = &r
            + &(&p * (rs_new / rs_old));

        rs_old = rs_new;
    }

    Ok(x)
}
