//! # Numerical Sparse Matrix Operations
//!
//! This module provides utilities for working with sparse matrices, particularly
//! in the Compressed Sparse Row (CSR) format. It includes functions for creating
//! CSR matrices, performing sparse matrix-vector multiplication, and solving
//! sparse linear systems using iterative methods like Conjugate Gradient.

use ndarray::ArrayD;
use sprs::{CsMat, TriMat};

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
pub fn csr_from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, f64)]) -> CsMat<f64> {
    let mut mat = TriMat::new((rows, cols));
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
/// A `Vec<f64>` representing the result of the multiplication.
///
/// # Panics
/// Panics if matrix and vector dimensions are not compatible.
pub fn sp_mat_vec_mul(matrix: &CsMat<f64>, vector: &[f64]) -> Vec<f64> {
    if matrix.cols() != vector.len() {
        panic!("Matrix and vector dimensions are not compatible for multiplication.");
    }
    let mut result = vec![0.0; matrix.rows()];
    for (i, row) in matrix.outer_iterator().enumerate() {
        let mut row_sum = 0.0;
        for (j, &val) in row.iter() {
            row_sum += val * vector[j];
        }
        result[i] = row_sum;
    }
    result
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
pub fn to_csr(arr: &Array) -> CsMat<f64> {
    assert_eq!(arr.ndim(), 2, "Input array must be 2D for CSR conversion.");
    let rows = arr.shape()[0];
    let cols = arr.shape()[1];

    let mut indptr = Vec::with_capacity(rows + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();

    indptr.push(0);
    let mut non_zeros = 0;

    for row in arr.outer_iter() {
        for (j, &val) in row.iter().enumerate() {
            if val != 0.0 {
                non_zeros += 1;
                indices.push(j);
                data.push(val);
            }
        }
        indptr.push(non_zeros);
    }

    CsMat::new((rows, cols), indptr, indices, data)
}

/// Converts a Compressed Sparse Row (CSR) matrix to a dense `ndarray::Array`.
///
/// # Arguments
/// * `matrix` - The sparse matrix in CSR format.
///
/// # Returns
/// An `ndarray::Array2<f64>` representing the dense matrix.
pub fn to_dense(matrix: &CsMat<f64>) -> Array2<f64> {
    matrix.to_dense()
}

use crate::numerical::matrix::Matrix;
use ndarray::Array2;

/// Computes the rank of a sparse matrix by converting to dense and performing RREF.
///
/// Note: This is inefficient for large sparse matrices and should only be used for small matrices.
///
/// # Arguments
/// * `matrix` - The sparse matrix.
///
/// # Returns
/// The rank of the matrix as a `usize`.
pub fn rank(matrix: &CsMat<f64>) -> usize {
    let dense_array2: Array2<f64> = matrix.to_dense();
    let rows = dense_array2.nrows();
    let cols = dense_array2.ncols();
    let mut dense_matrix = Matrix::new(rows, cols, dense_array2.into_raw_vec_and_offset().0);
    dense_matrix.rref()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    #[test]
    pub(crate) fn test_csr_from_triplets() {
        let triplets = vec![(0, 0, 1.0), (1, 2, 2.0), (2, 1, 3.0)];
        let mat = csr_from_triplets(3, 3, &triplets);
        assert_eq!(mat.rows(), 3);
        assert_eq!(mat.cols(), 3);
        assert_eq!(mat.nnz(), 3);
        assert_eq!(mat.get(0, 0), Some(&1.0));
        assert_eq!(mat.get(1, 2), Some(&2.0));
        assert_eq!(mat.get(2, 1), Some(&3.0));
        assert_eq!(mat.get(0, 1), None);
    }
    #[test]
    pub(crate) fn test_sp_mat_vec_mul() {
        let triplets = vec![(0, 0, 1.0), (0, 2, 2.0), (2, 1, 3.0)];
        let mat = csr_from_triplets(3, 3, &triplets);
        let vec = vec![10.0, 20.0, 30.0];
        let result = sp_mat_vec_mul(&mat, &vec);
        assert_eq!(result, vec![70.0, 0.0, 60.0]); // 1*10+2*30=70, 0, 3*20=60
    }
    #[test]
    pub(crate) fn test_to_csr() {
        let dense_arr = array![[1.0, 0.0, 2.0], [0.0, 0.0, 0.0], [3.0, 0.0, 4.0]].into_dyn();
        let csr_mat = to_csr(&dense_arr);
        assert_eq!(csr_mat.rows(), 3);
        assert_eq!(csr_mat.cols(), 3);
        assert_eq!(csr_mat.nnz(), 4);
        assert_eq!(csr_mat.get(0, 0), Some(&1.0));
        assert_eq!(csr_mat.get(0, 2), Some(&2.0));
        assert_eq!(csr_mat.get(2, 0), Some(&3.0));
        assert_eq!(csr_mat.get(2, 2), Some(&4.0));
        assert_eq!(csr_mat.get(1, 1), None);
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
pub fn solve_conjugate_gradient(
    a: &CsMat<f64>,
    b: &Array1<f64>,
    x0: Option<&Array1<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> Result<Array1<f64>, String> {
    let n = a.cols();
    if a.rows() != n || b.len() != n {
        return Err("Matrix and vector dimensions are incompatible.".to_string());
    }

    let mut x = x0.cloned().unwrap_or_else(|| Array1::zeros(n));
    let mut r = b - &(&*a * &x);
    let mut p = r.clone();
    let mut rs_old = r.dot(&r);

    if rs_old.sqrt() < tolerance {
        return Ok(x);
    }

    for _ in 0..max_iter {
        let ap = &*a * &p;
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
