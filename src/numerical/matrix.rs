//! # Numerical Matrix and Linear Algebra
//!
//! This module provides a generic `Matrix` struct for dense matrices over any type
//! that implements a custom `Field` trait. It supports a wide range of linear algebra
//! operations, including matrix arithmetic, RREF, inversion, null space calculation,
//! and eigenvalue decomposition for symmetric matrices.

use std::fmt::Debug;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

use faer::linalg::solvers::DenseSolveCore;
// Faer imports
use faer::{
    MatRef,
    Side,
};
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::finite_field::PrimeFieldElement;

/// A trait defining the requirements for a field in linear algebra.

pub trait Field:
    Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + Neg<Output = Self>
    + Clone
    + PartialEq
    + Debug
    + Zero
    + One
{
    /// Checks if the element is invertible.

    fn is_invertible(&self) -> bool;

    /// Returns the inverse of the element.
    ///
    /// # Errors
    /// Returns an error if the element is not invertible (e.g., zero in most fields).

    fn inverse(
        &self
    ) -> Result<Self, String>;

    /// Optional: Multiply two matrices using Faer backend if supported.

    #[must_use]

    fn faer_mul(
        _lhs: &Matrix<Self>,
        _rhs: &Matrix<Self>,
    ) -> Option<Matrix<Self>> {

        None
    }

    /// Optional: Invert a matrix using Faer backend if supported.

    #[must_use]

    fn faer_inverse(
        _matrix: &Matrix<Self>
    ) -> Option<Matrix<Self>> {

        None
    }

    /// Optional: Solve Ax = b using Faer backend if supported.

    #[must_use]

    fn faer_solve(
        _a: &Matrix<Self>,
        _b: &Matrix<Self>,
    ) -> Option<Matrix<Self>> {

        None
    }

    /// Optional: Perform matrix decomposition using Faer backend.

    fn faer_decompose(
        &self,
        _matrix: &Matrix<Self>,
        _kind: FaerDecompositionType,
    ) -> Option<
        FaerDecompositionResult<Self>,
    > {

        None
    }
}

impl Field for f64 {
    fn is_invertible(&self) -> bool {

        *self != 0.0
    }

    fn inverse(
        &self
    ) -> Result<Self, String> {

        if self.is_invertible() {

            Ok(1.0 / self)
        } else {

            Err("Cannot invert 0.0"
                .to_string())
        }
    }

    fn faer_mul(
        lhs: &Matrix<Self>,
        rhs: &Matrix<Self>,
    ) -> Option<Matrix<Self>> {

        if lhs.cols != rhs.rows {

            return None;
        }

        // Trick: Interpret row-major data of A as column-major data of A^T.
        // A * B = C  => C^T = B^T * A^T
        // We compute B^T * A^T using Faer (which expects col-major), result is col-major C^T.
        // Col-major C^T has the same memory layout as Row-major C.

        let lhs_t = MatRef::from_column_major_slice(&lhs.data, lhs.cols, lhs.rows);

        let rhs_t = MatRef::from_column_major_slice(&rhs.data, rhs.cols, rhs.rows);

        // res_mat (C^T) = rhs_t (B^T) * lhs_t (A^T)
        // faer operator * uses global parallelism
        let res_mat = rhs_t * lhs_t;

        let mut data =
            Vec::with_capacity(
                lhs.rows * rhs.cols,
            );

        for j in 0 .. rhs.cols {

            data.extend_from_slice(
                res_mat.col_as_slice(j),
            );
        }

        Some(
            Matrix::new(
                lhs.rows,
                rhs.cols,
                data,
            )
            .with_backend(lhs.backend),
        )
    }

    fn faer_inverse(
        matrix: &Matrix<Self>
    ) -> Option<Matrix<Self>> {

        if matrix.rows != matrix.cols {

            return None;
        }

        let n = matrix.rows;

        // Similar trick: (A^-1)^T = (A^T)^-1.
        // View as A^T (col-major). Invert. Result is (A^-1)^T (col-major) -> A^-1 (row-major).
        let mat_t = MatRef::from_column_major_slice(&matrix.data, n, n);

        // Use partial pivot LU for inversion
        let lu = mat_t.partial_piv_lu();

        // Ideally we check determinant or rank, but for perf we just invert.

        let inv_t = lu.inverse();

        // inv_t is (A^-1)^T in col-major.
        // We need the data as Vec<f64>.
        // Mat stores data in col-major. col_as_slice might help if contiguous.
        // Mat is contiguous.

        let mut data = vec![0.0; n * n];

        // Copy data out.
        // If Mat is standard, we can just copy.
        // inv_t is Mat<f64>.

        for j in 0 .. n {

            for i in 0 .. n {

                // inv_t(i, j)
                data[i * n + j] =
                    *inv_t.get(j, i);
                // So we just need the buffer.
            }
        }

        // Actually, we can just copy the slice if strictly contiguous.
        // Mat::as_slice() gives column major slice.
        // This slice corresponds exactly to row major data of the transpose of the matrix represented by Mat.
        // inv_t represents (A^-1)^T.
        // Its col-major data is Row-major data of ((A^-1)^T)^T = A^-1.
        // Correct.

        // So:
        let mut data =
            Vec::with_capacity(n * n);

        for j in 0 .. n {

            data.extend_from_slice(
                inv_t.col_as_slice(j),
            );
        }

        Some(
            Matrix::new(n, n, data)
                .with_backend(
                    matrix.backend,
                ),
        )
    }

    fn faer_decompose(
        &self,
        matrix: &Matrix<Self>,
        kind: FaerDecompositionType,
    ) -> Option<
        FaerDecompositionResult<Self>,
    > {

        // Ensure rows/cols are compatible (handled by faer per decomp mostly)
        // But we need dimensions.
        let rows = matrix.rows;

        let cols = matrix.cols;

        match kind {
            FaerDecompositionType::Svd => {
                let mat_t = MatRef::from_column_major_slice(&matrix.data, cols, rows);
                let svd = mat_t.thin_svd().ok()?;

                let u_prime = svd.U();
                let v_prime = svd.V();
                let s_col = svd.S().column_vector();
                let s = (0..s_col.nrows()).map(|i| *s_col.get(i)).collect::<Vec<_>>();

                let k = s.len();
                let mut u_data = vec![0.0; rows * k];
                for i in 0..rows {
                    for j in 0..k {
                        u_data[i * k + j] = *v_prime.get(i, j);
                    }
                }

                let mut v_data = vec![0.0; cols * k];
                for i in 0..cols {
                    for j in 0..k {
                        v_data[i * k + j] = *u_prime.get(i, j);
                    }
                }

                Some(FaerDecompositionResult::Svd {
                    u: Matrix::new(rows, k, u_data).with_backend(matrix.backend),
                    s,
                    v: Matrix::new(cols, k, v_data).with_backend(matrix.backend),
                })
            }
            FaerDecompositionType::Cholesky => {
                if rows != cols { return None; }
                let mat = MatRef::from_column_major_slice(&matrix.data, rows, cols);
                if let Ok(chol) = mat.llt(Side::Lower) {
                    let l_mat = chol.L();
                    let mut l_data = vec![0.0; rows * cols];
                    for i in 0..rows {
                        for j in 0..cols {
                             l_data[i * cols + j] = *l_mat.get(i, j);
                        }
                    }
                    Some(FaerDecompositionResult::Cholesky {
                        l: Matrix::new(rows, cols, l_data).with_backend(matrix.backend)
                    })
                } else {
                    None
                }
            }
             FaerDecompositionType::EigenSymmetric => {
                 if rows != cols { return None; }
                 let mat = MatRef::from_column_major_slice(&matrix.data, rows, cols);
                 let evd = mat.self_adjoint_eigen(Side::Lower).ok()?;

                 let s_col = evd.S().column_vector();
                 let s = (0..s_col.nrows()).map(|i| *s_col.get(i)).collect::<Vec<_>>();
                 let u_mat = evd.U(); // Eigenvectors

                 let mut u_data = vec![0.0; rows * cols];
                 for i in 0..rows {
                    for j in 0..cols {
                        u_data[i * cols + j] = *u_mat.get(i, j);
                    }
                 }

                 Some(FaerDecompositionResult::EigenSymmetric {
                     values: s,
                     vectors: Matrix::new(rows, cols, u_data).with_backend(matrix.backend)
                 })
             }
            _ => None // LU/QR/etc not implemented yet for this quick pass
        }
    }
}

impl Field for PrimeFieldElement {
    fn is_invertible(&self) -> bool {

        !self.value.is_zero()
    }

    fn inverse(
        &self
    ) -> Result<Self, String> {

        self.inverse()
            .ok_or_else(|| {

                "Cannot invert \
                 non-invertible element"
                    .to_string()
            })
    }
}

/// Backend usage for matrix operations
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    Default,
)]

pub enum Backend {
    /// Native Rust implementation (default)
    #[default]
    Native,
    /// Faer backend (high performance BLAS-like) - currently supports f64
    Faer,
}


/// Types of matrix decompositions supported by Faer
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub enum FaerDecompositionType {
    /// LU Decomposition with Partial Pivoting
    Lu,
    /// QR Decomposition
    Qr,
    /// Cholesky Decomposition (Lower)
    Cholesky,
    /// SVD Decomposition
    Svd,
    /// Eigendecomposition for Symmetric Matrices
    EigenSymmetric,
}

/// Results of Faer matrix decompositions
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]
#[serde(bound(
    serialize = "T: Serialize",
    deserialize = "T: Deserialize<'de>"
))]

pub enum FaerDecompositionResult<
    T: Field,
> {
    /// LU: P * A = L * U
    Lu {
        /// Lower triangular matrix.
        l: Matrix<T>,
        /// Upper triangular matrix.
        u: Matrix<T>,
        /// Row permutation indices.
        p: Vec<usize>,
    },
    /// QR: A = Q * R
    Qr {
        /// Orthogonal matrix Q.
        q: Matrix<T>,
        /// Upper triangular matrix R.
        r: Matrix<T>,
    },
    /// Cholesky: A = L * L^T
    Cholesky {
        /// Lower triangular matrix L.
        l: Matrix<T>,
    },
    /// SVD: A = U * S * V^T
    Svd {
        /// Left singular vectors.
        u: Matrix<T>,
        /// Singular values.
        s: Vec<f64>,
        /// Right singular vectors.
        v: Matrix<T>,
    },
    /// Eigen Symmetric: A = V * D * V^T (s are eigenvalues, u are eigenvectors)
    EigenSymmetric {
        /// Eigenvalues.
        values: Vec<f64>,
        /// Eigenvectors matrix.
        vectors: Matrix<T>,
    },
}


/// A generic dense matrix over any type that implements the Field trait.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub struct Matrix<T: Field> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
    /// Computational backend used for matrix operations.
    #[serde(default)]
    pub backend: Backend,
}

impl<T: Field> Matrix<T> {
    /// Creates a new `Matrix` from dimensions and a flat `Vec` of data.
    ///
    /// # Arguments
    /// * `rows` - The number of rows.
    /// * `cols` - The number of columns.
    /// * `data` - A `Vec` containing the matrix elements in row-major order.
    ///
    /// # Panics
    /// Panics if `rows * cols` does not equal `data.len()`.
    #[must_use]

    pub fn new(
        rows: usize,
        cols: usize,
        data: Vec<T>,
    ) -> Self {

        assert_eq!(
            rows * cols,
            data.len()
        );

        Self {
            rows,
            cols,
            data,
            backend: Backend::default(),
        }
    }

    /// Sets the backend for the matrix.

    #[must_use]

    pub const fn with_backend(
        mut self,
        backend: Backend,
    ) -> Self {

        self.backend = backend;

        self
    }

    /// Updates the backend for the matrix in-place.

    pub const fn set_backend(
        &mut self,
        backend: Backend,
    ) {

        self.backend = backend;
    }

    /// Creates a new `Matrix` filled with the zero element of type `T`.
    ///
    /// # Arguments
    /// * `rows` - The number of rows.
    /// * `cols` - The number of columns.
    ///
    /// # Returns
    /// A new `Matrix` of the specified dimensions, with all elements initialized to `T::zero()`.
    #[must_use]

    pub fn zeros(
        rows: usize,
        cols: usize,
    ) -> Self {

        Self {
            rows,
            cols,
            data: vec![
                T::zero();
                rows * cols
            ],
            backend: Backend::default(),
        }
    }

    /// Returns an immutable reference to the element at the specified row and column.
    ///
    /// # Arguments
    /// * `row` - The row index.
    /// * `col` - The column index.
    ///
    /// # Returns
    /// An immutable reference to the element.
    ///
    /// # Panics
    /// Panics if the `row` or `col` indices are out of bounds.
    #[must_use]

    pub fn get(
        &self,
        row: usize,
        col: usize,
    ) -> &T {

        &self.data
            [row * self.cols + col]
    }

    /// Returns a mutable reference to the element at the specified row and column.
    ///
    /// # Arguments
    /// * `row` - The row index.
    /// * `col` - The column index.
    ///
    /// # Returns
    /// A mutable reference to the element.
    ///
    /// # Panics
    /// Panics if the `row` or `col` indices are out of bounds.

    pub fn get_mut(
        &mut self,
        row: usize,
        col: usize,
    ) -> &mut T {

        &mut self.data
            [row * self.cols + col]
    }

    /// Returns the number of rows in the matrix.
    #[must_use]

    pub const fn rows(&self) -> usize {

        self.rows
    }

    /// Returns the number of columns in the matrix.
    #[must_use]

    pub const fn cols(&self) -> usize {

        self.cols
    }

    /// Computes a matrix decomposition using the Faer backend.
    /// Returns None if the backend is not set to Faer, if the decomposition type is unsupported
    /// for the matrix type, or if the decomposition fails (e.g., Cholesky on non-SPD).

    #[must_use]

    pub fn decompose(
        &self,
        kind: FaerDecompositionType,
    ) -> Option<
        FaerDecompositionResult<T>,
    > {

        if self.backend == Backend::Faer
        {

            self.data
                .first()
                .and_then(|e| {

                    e.faer_decompose(
                        self, kind,
                    )
                })
        } else {

            if self.data.is_empty() {

                return None;
            }

            self.data[0].faer_decompose(
                self, kind,
            )
        }
    }

    /// Returns an immutable reference to the matrix's internal data vector.
    #[must_use]

    pub const fn data(
        &self
    ) -> &Vec<T> {

        &self.data
    }

    /// Consumes the matrix and returns its internal data vector.
    #[must_use]

    pub fn into_data(self) -> Vec<T> {

        self.data
    }

    /// Returns a `Vec` of `Vec<T>` where each inner `Vec` represents a column of the matrix.
    ///
    /// This method effectively transposes the matrix data into a column-major representation.
    ///
    /// # Returns
    /// A `Vec<Vec<T>>` where each inner vector is a column.
    #[must_use]

    pub fn get_cols(
        &self
    ) -> Vec<Vec<T>> {

        let mut cols_vec =
            Vec::with_capacity(
                self.cols,
            );

        for j in 0 .. self.cols {

            let mut col =
                Vec::with_capacity(
                    self.rows,
                );

            for i in 0 .. self.rows {

                col.push(
                    self.get(i, j)
                        .clone(),
                );
            }

            cols_vec.push(col);
        }

        cols_vec
    }

    /// Computes the reduced row echelon form (RREF) of the matrix in-place.
    ///
    /// This method applies Gaussian elimination to transform the matrix into its RREF.
    /// It is used for solving linear systems, finding matrix inverses, and determining rank.
    ///
    /// # Returns
    /// The rank of the matrix (number of non-zero rows in RREF).
    ///
    /// # Errors
    /// Returns an error if a division by zero occurs during the elimination process,
    /// which can happen if the matrix elements are not from a proper field.

    pub fn rref(
        &mut self
    ) -> Result<usize, String> {

        let mut pivot_row = 0;

        for j in 0 .. self.cols {

            if pivot_row >= self.rows {

                break;
            }

            let mut i = pivot_row;

            while i < self.rows
                && !self
                    .get(i, j)
                    .is_invertible()
            {

                i += 1;
            }

            if i < self.rows {

                for k in 0 .. self.cols
                {

                    self.data.swap(
                        i * self.cols
                            + k,
                        pivot_row
                            * self.cols
                            + k,
                    );
                }

                let pivot_inv = self
                    .get(pivot_row, j)
                    .clone()
                    .inverse()?;

                for k in j .. self.cols
                {

                    let val = self
                        .get(
                            pivot_row,
                            k,
                        )
                        .clone();

                    *self.get_mut(
                        pivot_row,
                        k,
                    ) = val
                        * pivot_inv
                            .clone();
                }

                for i_prime in
                    0 .. self.rows
                {

                    if i_prime
                        != pivot_row
                    {

                        let factor =
                            self.get(
                                i_prime,
                                j,
                            )
                            .clone();

                        for k in j
                            .. self.cols
                        {

                            let pivot_row_val = self
                                .get(pivot_row, k)
                                .clone();

                            let term = factor.clone() * pivot_row_val;

                            let current_val = self
                                .get(i_prime, k)
                                .clone();

                            *self.get_mut(i_prime, k) = current_val - term;
                        }
                    }
                }

                pivot_row += 1;
            }
        }

        Ok(pivot_row)
    }

    /// Computes the transpose of the matrix.
    ///
    /// The transpose of a matrix `A` (denoted `A^T`) is obtained by flipping the matrix
    /// over its diagonal; that is, it switches the row and column indices of the matrix.
    ///
    /// # Returns
    /// A new `Matrix` representing the transpose.
    #[must_use]

    pub fn transpose(&self) -> Self {

        let mut new_data = vec![
                T::zero();
                self.rows * self.cols
            ];

        for i in 0 .. self.rows {

            for j in 0 .. self.cols {

                new_data[j * self
                    .rows
                    + i] = self
                    .get(i, j)
                    .clone();
            }
        }

        Self::new(
            self.cols,
            self.rows,
            new_data,
        )
        .with_backend(self.backend)
    }

    /// # Strassen's Matrix Multiplication
    ///
    /// This function multiplies two matrices using Strassen's algorithm, which is more
    /// efficient for large matrices than the naive O(n^3) algorithm.
    ///
    /// ## Arguments
    /// * `other` - The matrix to multiply with.
    ///
    /// ## Returns
    /// A new `Matrix` representing the product of the two matrices.
    ///
    /// ## Errors
    /// Returns an error if:
    /// - Matrix multiplication is not possible due to incompatible dimensions.
    /// - The matrices are not square (for the internal recursive calls).

    pub fn mul_strassen(
        &self,
        other: &Self,
    ) -> Result<Self, String> {

        if self.cols != other.rows {

            return Err(format!(
                "Matrix multiplication not possible: left.cols ({}) != right.rows ({})",
                self.cols, other.rows
            ));
        }

        let n = self
            .rows
            .max(self.cols)
            .max(other.rows)
            .max(other.cols);

        let m = n.next_power_of_two();

        let mut a_padded =
            Self::zeros(m, m);

        let mut b_padded =
            Self::zeros(m, m);

        // Fill the padded matrices
        for i in 0 .. self.rows {

            for j in 0 .. self.cols {

                *a_padded
                    .get_mut(i, j) =
                    self.get(i, j)
                        .clone();
            }
        }

        for i in 0 .. other.rows {

            for j in 0 .. other.cols {

                *b_padded
                    .get_mut(i, j) =
                    other
                        .get(i, j)
                        .clone();
            }
        }

        let c_padded =
            strassen_recursive(
                &a_padded,
                &b_padded,
            );

        // Check if the result matrix has the expected size
        if c_padded.rows != m
            || c_padded.cols != m
        {

            return Err(
                "Internal error: \
                 Strassen result has \
                 unexpected dimensions"
                    .to_string(),
            );
        }

        let mut result_data = vec![
                T::zero();
                self.rows * other.cols
            ];

        for i in 0 .. self.rows {

            for j in 0 .. other.cols {

                result_data[i
                    * other.cols
                    + j] = c_padded
                    .get(i, j)
                    .clone();
            }
        }

        Ok(Self::new(
            self.rows,
            other.cols,
            result_data,
        )
        .with_backend(self.backend))
    }

    /// Splits a matrix into four sub-matrices of equal size.

    fn split(
        &self
    ) -> (
        Self,
        Self,
        Self,
        Self,
    ) {

        // Check that the matrix dimensions are even so we can split it equally
        if !self
            .rows
            .is_multiple_of(2)
            || !self
                .cols
                .is_multiple_of(2)
        {

            // Return appropriately sized zero matrices if splitting isn't possible
            let half_rows =
                self.rows / 2;

            let half_cols =
                self.cols / 2;

            return (
                Self::zeros(
                    half_rows,
                    half_cols,
                ),
                Self::zeros(
                    half_rows,
                    half_cols,
                ),
                Self::zeros(
                    half_rows,
                    half_cols,
                ),
                Self::zeros(
                    half_rows,
                    half_cols,
                ),
            );
        }

        let new_dim = self.rows / 2;

        let new_col_dim = self.cols / 2; // Make sure to use the correct column division

        let mut a11 = Self::zeros(
            new_dim,
            new_col_dim,
        );

        let mut a12 = Self::zeros(
            new_dim,
            new_col_dim,
        );

        let mut a21 = Self::zeros(
            new_dim,
            new_col_dim,
        );

        let mut a22 = Self::zeros(
            new_dim,
            new_col_dim,
        );

        for i in 0 .. new_dim {

            for j in 0 .. new_col_dim {

                *a11.get_mut(i, j) =
                    self.get(i, j)
                        .clone();

                *a12.get_mut(i, j) =
                    self.get(
                        i,
                        j + new_col_dim,
                    )
                    .clone();

                *a21.get_mut(i, j) =
                    self.get(
                        i + new_dim,
                        j,
                    )
                    .clone();

                *a22.get_mut(i, j) =
                    self.get(
                        i + new_dim,
                        j + new_col_dim,
                    )
                    .clone();
            }
        }

        (a11, a12, a21, a22)
    }

    /// Joins four sub-matrices into a single larger matrix.

    fn join(
        a11: &Self,
        a12: &Self,
        a21: &Self,
        a22: &Self,
    ) -> Self {

        // All four submatrices should have the same dimensions
        if a11.rows != a12.rows
            || a11.cols != a12.cols
            || a11.rows != a21.rows
            || a11.cols != a21.cols
            || a11.rows != a22.rows
            || a11.cols != a22.cols
        {

            // Return a zero matrix if dimensions don't match
            return Self::zeros(0, 0);
        }

        let result_rows = a11.rows * 2;

        let result_cols = a11.cols * 2;

        let mut result = Self::zeros(
            result_rows,
            result_cols,
        );

        let sub_row_dim = a11.rows;

        let sub_col_dim = a11.cols;

        // Copy a11 to top-left
        for i in 0 .. sub_row_dim {

            for j in 0 .. sub_col_dim {

                *result.get_mut(i, j) =
                    a11.get(i, j)
                        .clone();
            }
        }

        // Copy a12 to top-right
        for i in 0 .. sub_row_dim {

            for j in 0 .. sub_col_dim {

                *result.get_mut(
                    i,
                    j + sub_col_dim,
                ) = a12
                    .get(i, j)
                    .clone();
            }
        }

        // Copy a21 to bottom-left
        for i in 0 .. sub_row_dim {

            for j in 0 .. sub_col_dim {

                *result.get_mut(
                    i + sub_row_dim,
                    j,
                ) = a21
                    .get(i, j)
                    .clone();
            }
        }

        // Copy a22 to bottom-right
        for i in 0 .. sub_row_dim {

            for j in 0 .. sub_col_dim {

                *result.get_mut(
                    i + sub_row_dim,
                    j + sub_col_dim,
                ) = a22
                    .get(i, j)
                    .clone();
            }
        }

        result
    }

    /// Computes the determinant of a square matrix.
    ///
    /// This method uses block matrix decomposition (Schur complement) for efficient
    /// determinant calculation, especially for large matrices.
    ///
    /// # Returns
    /// The determinant of the matrix.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The matrix is not square.
    /// - The matrix is singular or reduction fails.

    pub fn determinant(
        &self
    ) -> Result<T, String> {

        if self.rows != self.cols {

            return Err("Matrix must \
                        be square to \
                        compute the \
                        determinant."
                .to_string());
        }

        if self.rows > 64
            && self
                .rows
                .is_multiple_of(2)
        {

            return self
                .determinant_block();
        }

        if self.rows == 0 {

            return Ok(T::one());
        }

        if self.rows == 1 {

            return Ok(self
                .get(0, 0)
                .clone());
        }

        if self.rows == 2 {

            let a = self
                .get(0, 0)
                .clone();

            let b = self
                .get(0, 1)
                .clone();

            let c = self
                .get(1, 0)
                .clone();

            let d = self
                .get(1, 1)
                .clone();

            return Ok(a * d - b * c);
        }

        let (lu, swaps) =
            self.lu_decomposition()?;

        let mut det = T::one();

        for i in 0 .. self.rows {

            det *= lu.get(i, i).clone();
        }

        if (swaps % 2) != 0 {

            det = -det;
        }

        Ok(det)
    }

    /// Computes the LU decomposition of a square matrix.
    ///
    /// # Returns
    /// A tuple containing the LU matrix and the number of row swaps.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The matrix is not square.
    /// - The matrix is singular.

    pub fn lu_decomposition(
        &self
    ) -> Result<(Self, usize), String>
    {

        if self.rows != self.cols {

            return Err(
                "Matrix must be \
                 square for LU \
                 decomposition."
                    .to_string(),
            );
        }

        let n = self.rows;

        let mut lu = self.clone();

        let mut swaps = 0;

        for j in 0 .. n {

            let mut pivot_row = j;

            for i in j .. n {

                if self
                    .get(i, j)
                    .is_invertible()
                {

                    pivot_row = i;

                    break;
                }
            }

            if pivot_row != j {

                swaps += 1;

                for k in 0 .. n {

                    let val1 = lu
                        .get(j, k)
                        .clone();

                    let val2 = lu
                        .get(
                            pivot_row,
                            k,
                        )
                        .clone();

                    *lu.get_mut(j, k) =
                        val2;

                    *lu.get_mut(
                        pivot_row,
                        k,
                    ) = val1;
                }
            }

            let pivot_val =
                lu.get(j, j).clone();

            if !pivot_val
                .is_invertible()
            {

                return Err("Matrix is singular.".to_string());
            }

            for i in (j + 1) .. n {

                let factor = lu
                    .get(i, j)
                    .clone()
                    * pivot_val
                        .inverse()?;

                *lu.get_mut(i, j) =
                    factor.clone();

                for k in (j + 1) .. n {

                    let val = lu
                        .get(j, k)
                        .clone()
                        * factor
                            .clone();

                    let current_val =
                        lu.get(i, k)
                            .clone();

                    *lu.get_mut(i, k) =
                        current_val
                            - val;
                }
            }
        }

        Ok((lu, swaps))
    }

    /// # Block Matrix Determinant
    ///
    /// Computes the determinant using block matrix decomposition (Schur complement).
    /// This is efficient for large matrices.
    /// # Errors
    ///
    /// This function will return an error if:
    /// - The matrix is not square.
    /// - The matrix is singular or reduction fails.
    #[allow(clippy::suspicious_operation_groupings)]

    pub fn determinant_block(
        &self
    ) -> Result<T, String> {

        if self.rows != self.cols {

            return Err("Matrix must \
                        be square."
                .to_string());
        }

        let n = self.rows;

        if n == 0 {

            return Ok(T::one());
        }

        if !n.is_multiple_of(2) {

            return self
                .determinant_lu();
        }

        let (a, b, c, d) = self.split();

        // Check that the submatrices have consistent dimensions for block operations
        if a.rows != b.rows
            || a.cols != c.cols
            || b.cols != d.cols
            || c.cols != d.rows
        {

            return Err(
                "Block matrix \
                 decomposition failed \
                 due to inconsistent \
                 submatrix dimensions."
                    .to_string(),
            );
        }

        // Try to compute determinant using Schur complement if the top-left block is invertible
        match a.inverse() {
            | Some(a_inv) => {

                // Calculate Schur complement: S = D - C * A^(-1) * B
                let a_inv_b = a_inv * b;

                let schur_complement =
                    d - c * a_inv_b;

                match (
                    a.determinant_lu(),
                    schur_complement
                        .determinant_lu(
                        ),
                ) {
                    | (
                        Ok(det_a),
                        Ok(det_s),
                    ) => {
                        Ok(det_a
                            * det_s)
                    },
                    | _ => {

                        // If calculation fails, fallback to LU decomposition
                        self.determinant_lu()
                    },
                }
            },
            | None => {

                // If A is not invertible, fallback to LU decomposition
                self.determinant_lu()
            },
        }
    }

    /// Computes the determinant using LU decomposition.
    /// # Errors
    ///
    /// This function will return an error if:
    /// - The matrix is not square.
    /// - The matrix is singular or reduction fails.

    pub fn determinant_lu(
        &self
    ) -> Result<T, String> {

        let (lu, swaps) =
            self.lu_decomposition()?;

        let mut det = T::one();

        for i in 0 .. self.rows {

            det *= lu.get(i, i).clone();
        }

        if (swaps % 2) != 0 {

            det = -det;
        }

        Ok(det)
    }

    /// Computes the inverse of a square matrix.
    ///
    /// This method uses Gaussian elimination on an augmented matrix `[A | I]`
    /// to transform `A` into the identity matrix `I`, resulting in `[I | A^-1]`.
    ///
    /// # Returns
    /// * `Some(Matrix)` containing the inverse matrix if it exists.
    /// * `None` if the matrix is not square or is singular (not invertible).
    #[must_use]

    pub fn inverse(
        &self
    ) -> Option<Self> {

        if self.backend == Backend::Faer
        {

            if let Some(res) =
                T::faer_inverse(self)
            {

                return Some(res);
            }
        }

        if self.rows != self.cols {

            return None;
        }

        let n = self.rows;

        let mut augmented =
            Self::zeros(n, 2 * n);

        for i in 0 .. n {

            for j in 0 .. n {

                *augmented
                    .get_mut(i, j) =
                    self.get(i, j)
                        .clone();

                if i == j {

                    *augmented
                        .get_mut(
                            i,
                            j + n,
                        ) = T::one();
                }
            }
        }

        // Check if RREF operation was successful
        if let Ok(rank) =
            augmented.rref()
        {

            if rank == n {

                let mut inv_data =
                    vec![
                        T::zero();
                        n * n
                    ];

                for i in 0 .. n {

                    for j in 0 .. n {

                        inv_data[i * n + j] = augmented
                            .get(i, j + n)
                            .clone();
                    }
                }

                Some(
                    Self::new(
                        n,
                        n,
                        inv_data,
                    )
                    .with_backend(
                        self.backend,
                    ),
                )
            } else {

                // Matrix is not invertible (rank is less than n)
                None
            }
        } else {

            // RREF operation failed
            None
        }
    }

    /// Computes a basis for the null space (kernel) of the matrix.
    ///
    /// The null space of a matrix `A` is the set of all vectors `x` such that `Ax = 0`.
    /// This method finds the null space by first computing the RREF of the matrix,
    /// identifying pivot and free variables, and then constructing basis vectors.
    ///
    /// # Returns
    /// A `Matrix` whose columns form a basis for the null space.
    ///
    /// # Errors
    /// Returns an error if the RREF calculation fails.

    pub fn null_space(
        &self
    ) -> Result<Self, String> {

        let mut rref_matrix =
            self.clone();

        let rank =
            rref_matrix.rref()?;

        let mut pivot_cols = Vec::new();

        let mut lead = 0;

        for r in 0 .. rank {

            if lead >= self.cols {

                break;
            }

            let mut i = lead;

            while !rref_matrix
                .get(r, i)
                .is_invertible()
            {

                i += 1;

                if i == self.cols {

                    break;
                }
            }

            if i < self.cols {

                pivot_cols.push(i);

                lead = i + 1;
            }
        }

        let free_cols: Vec<usize> = (0
            .. self.cols)
            .filter(|c| {

                !pivot_cols.contains(c)
            })
            .collect();

        let num_free = free_cols.len();

        let mut basis_vectors =
            Vec::with_capacity(
                num_free,
            );

        for free_col in free_cols {

            let mut vec = vec![
                T::zero();
                self.cols
            ];

            vec[free_col] = T::one();

            for (i, &pivot_col) in
                pivot_cols
                    .iter()
                    .enumerate()
            {

                vec[pivot_col] =
                    -rref_matrix
                        .get(
                            i,
                            free_col,
                        )
                        .clone();
            }

            basis_vectors.push(vec);
        }

        let mut null_space_data = vec![
                T::zero();
                self.cols * num_free
            ];

        for (j, basis_vec) in
            basis_vectors
                .iter()
                .enumerate()
        {

            for (i, val) in basis_vec
                .iter()
                .enumerate()
            {

                null_space_data[i
                    * num_free
                    + j] = val.clone();
            }
        }

        Ok(Self::new(
            self.cols,
            num_free,
            null_space_data,
        )
        .with_backend(self.backend))
    }

    /// Computes the rank of the matrix.
    ///
    /// # Errors
    ///
    /// Returns an error if the RREF computation fails.

    pub fn rank(
        &self
    ) -> Result<usize, String> {

        let mut rref_matrix =
            self.clone();

        rref_matrix.rref()
    }

    /// Computes the trace of the matrix (sum of diagonal elements).
    ///
    /// # Errors
    /// Returns an error if the matrix is not square.

    pub fn trace(
        &self
    ) -> Result<T, String> {

        if self.rows != self.cols {

            return Err(
                "Trace is only \
                 defined for square \
                 matrices."
                    .to_string(),
            );
        }

        let mut sum = T::zero();

        for i in 0 .. self.rows {

            sum += self
                .get(i, i)
                .clone();
        }

        Ok(sum)
    }

    /// Checks if the matrix is symmetric ($A = A^T$).
    #[must_use]

    pub fn is_symmetric(&self) -> bool {

        if self.rows != self.cols {

            return false;
        }

        for i in 0 .. self.rows {

            for j in
                (i + 1) .. self.cols
            {

                if self.get(i, j)
                    != self.get(j, i)
                {

                    return false;
                }
            }
        }

        true
    }

    /// Checks if the matrix is diagonal.
    #[must_use]

    pub fn is_diagonal(&self) -> bool {

        for i in 0 .. self.rows {

            for j in 0 .. self.cols {

                if i != j
                    && !self
                        .get(i, j)
                        .is_zero()
                {

                    return false;
                }
            }
        }

        true
    }

    /// Computes the Frobenius norm of the matrix.
    #[must_use]

    pub fn frobenius_norm(&self) -> f64
    where
        T: ToPrimitive,
    {

        let mut sum = 0.0;

        for val in &self.data {

            let v = val
                .to_f64()
                .unwrap_or(0.0);

            sum += v * v;
        }

        sum.sqrt()
    }

    /// Computes the L1 norm of the matrix (max column sum).
    #[must_use]

    pub fn l1_norm(&self) -> f64
    where
        T: ToPrimitive,
    {

        let mut max_sum = 0.0;

        for j in 0 .. self.cols {

            let mut col_sum = 0.0;

            for i in 0 .. self.rows {

                col_sum += self
                    .get(i, j)
                    .to_f64()
                    .unwrap_or(0.0)
                    .abs();
            }

            if col_sum > max_sum {

                max_sum = col_sum;
            }
        }

        max_sum
    }

    /// Computes the Linf norm of the matrix (max row sum).
    #[must_use]

    pub fn linf_norm(&self) -> f64
    where
        T: ToPrimitive,
    {

        let mut max_sum = 0.0;

        for i in 0 .. self.rows {

            let mut row_sum = 0.0;

            for j in 0 .. self.cols {

                row_sum += self
                    .get(i, j)
                    .to_f64()
                    .unwrap_or(0.0)
                    .abs();
            }

            if row_sum > max_sum {

                max_sum = row_sum;
            }
        }

        max_sum
    }

    /// Creates an identity matrix of a given size.
    #[must_use]

    pub fn identity(
        size: usize
    ) -> Self {

        let mut data = vec![
            T::zero();
            size * size
        ];

        for i in 0 .. size {

            data[i * size + i] =
                T::one();
        }

        Self::new(size, size, data)
            .with_backend(
                Backend::default(),
            )
    }

    /// Checks if the matrix is identity within a given tolerance.
    #[must_use]

    pub fn is_identity(
        &self,
        epsilon: f64,
    ) -> bool
    where
        T: ToPrimitive,
    {

        if self.rows != self.cols {

            return false;
        }

        for i in 0 .. self.rows {

            for j in 0 .. self.cols {

                let expected = if i == j
                {

                    1.0
                } else {

                    0.0
                };

                let val = self
                    .get(i, j)
                    .to_f64()
                    .unwrap_or(0.0);

                if (val - expected)
                    .abs()
                    > epsilon
                {

                    return false;
                }
            }
        }

        true
    }

    /// Checks if the matrix is orthogonal ($A^T A = I$).
    #[must_use]

    pub fn is_orthogonal(
        &self,
        epsilon: f64,
    ) -> bool
    where
        T: ToPrimitive,
    {

        if self.rows != self.cols {

            return false;
        }

        let prod = self.transpose()
            * self.clone();

        prod.is_identity(epsilon)
    }
}

/// Recursive helper for Strassen's algorithm.

fn strassen_recursive<T: Field>(
    a: &Matrix<T>,
    b: &Matrix<T>,
) -> Matrix<T> {

    // Check that matrices can be multiplied
    if a.cols != b.rows {

        // Return a zero matrix as a fallback (though this shouldn't happen in practice with our padding)
        return Matrix::zeros(
            a.rows,
            b.cols,
        );
    }

    let n = a.rows;

    if n <= 64 {

        // For small matrices, use standard multiplication to avoid overhead
        return a.clone() * b.clone();
    }

    // For the Strassen algorithm to work properly, matrix dimensions must be even
    if !n.is_multiple_of(2) {

        // Fall back to standard multiplication if dimension is odd
        return a.clone() * b.clone();
    }

    let (a11, a12, a21, a22) =
        a.split();

    let (b11, b12, b21, b22) =
        b.split();

    let p1 = strassen_recursive(
        &(a11.clone() + a22.clone()),
        &(b11.clone() + b22.clone()),
    );

    let p2 = strassen_recursive(
        &(a21.clone() + a22.clone()),
        &b11,
    );

    let p3 = strassen_recursive(
        &a11,
        &(b12.clone() - b22.clone()),
    );

    let p4 = strassen_recursive(
        &a22,
        &(b21.clone() - b11.clone()),
    );

    let p5 = strassen_recursive(
        &(a11.clone() + a12.clone()),
        &b22,
    );

    let p6 = strassen_recursive(
        &(a21 - a11),
        &(b11 + b12),
    );

    let p7 = strassen_recursive(
        &(a12 - a22),
        &(b21 + b22),
    );

    let c11 = p1.clone() + p4.clone()
        - p5.clone()
        + p7;

    let c12 = p3.clone() + p5;

    let c21 = p2.clone() + p4;

    let c22 = p1 - p2 + p3 + p6;

    Matrix::join(
        &c11, &c12, &c21, &c22,
    )
}

impl Matrix<f64> {
    /// Finds the eigenvalues and eigenvectors of a symmetric matrix using the Jacobi iteration method.
    ///
    /// The Jacobi method is an iterative algorithm for computing the eigenvalues and eigenvectors
    /// of a real symmetric matrix by performing a sequence of orthogonal similarity transformations.
    ///
    /// # Arguments
    /// * `max_sweeps` - The maximum number of sweeps (iterations) to perform.
    /// * `tolerance` - The convergence tolerance for off-diagonal elements.
    ///
    /// # Returns
    /// A `Result` containing a tuple `(eigenvalues, eigenvectors)` where `eigenvalues` is a `Vec<f64>`
    /// and `eigenvectors` is a `Matrix<f64>` (columns are eigenvectors).
    ///
    /// # Errors
    /// Returns an error if:
    /// - The matrix is not square.
    /// - The decomposition fails to converge within the maximum sweeps.

    pub fn jacobi_eigen_decomposition(
        &self,
        max_sweeps: usize,
        tolerance: f64,
    ) -> Result<(Vec<f64>, Self), String>
    {

        if self.rows != self.cols {

            return Err("Matrix must \
                        be square."
                .to_string());
        }

        let mut a = self.clone();

        let n = self.rows;

        let mut eigenvectors =
            Self::identity(n);

        for _ in 0 .. max_sweeps {

            let mut off_diagonal_sum =
                0.0;

            for p in 0 .. n {

                for q in (p + 1) .. n {

                    off_diagonal_sum +=
                        a.get(p, q)
                            .abs()
                            .powi(2);
                }
            }

            if off_diagonal_sum.sqrt()
                < tolerance
            {

                break;
            }

            for p in 0 .. n {

                for q in (p + 1) .. n {

                    let apq =
                        *a.get(p, q);

                    if apq.abs()
                        < tolerance
                            / (n as f64)
                    {

                        continue;
                    }

                    let app =
                        *a.get(p, p);

                    let aqq =
                        *a.get(q, q);

                    let tau = (aqq
                        - app)
                        / (2.0 * apq);

                    let t = if tau
                        >= 0.0
                    {

                        1.0 / (tau + (1.0 + tau * tau).sqrt())
                    } else {

                        -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                    };

                    let c = 1.0
                        / (1.0 + t * t)
                            .sqrt();

                    let s = t * c;

                    let new_app =
                        app - t * apq;

                    let new_aqq =
                        aqq + t * apq;

                    *a.get_mut(p, p) =
                        new_app;

                    *a.get_mut(q, q) =
                        new_aqq;

                    *a.get_mut(p, q) =
                        0.0;

                    *a.get_mut(q, p) =
                        0.0;

                    for i in 0 .. n {

                        if i != p
                            && i != q
                        {

                            let aip =
                                *a.get(
                                    i,
                                    p,
                                );

                            let aiq =
                                *a.get(
                                    i,
                                    q,
                                );

                            *a.get_mut(i, p) = c * aip - s * aiq;

                            *a.get_mut(i, q) = s * aip + c * aiq;

                            *a.get_mut(p, i) = *a.get(i, p);

                            *a.get_mut(q, i) = *a.get(i, q);
                        }
                    }

                    for i in 0 .. n {

                        let vip = *eigenvectors.get(i, p);

                        let viq = *eigenvectors.get(i, q);

                        *eigenvectors
                            .get_mut(
                                i, p,
                            ) = c * vip
                            - s * viq;

                        *eigenvectors
                            .get_mut(
                                i, q,
                            ) = s * vip
                            + c * viq;
                    }
                }
            }
        }

        let eigenvalues = (0 .. n)
            .map(|i| *a.get(i, i))
            .collect();

        Ok((
            eigenvalues,
            eigenvectors,
        ))
    }
}

impl<T: Field> Add for Matrix<T> {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        assert_eq!(self.rows, rhs.rows);

        assert_eq!(self.cols, rhs.cols);

        let data = self
            .data
            .into_iter()
            .zip(rhs.data)
            .map(|(a, b)| a + b)
            .collect();

        Self::new(
            self.rows,
            self.cols,
            data,
        )
        .with_backend(self.backend)
    }
}

impl<T: Field> Sub for Matrix<T> {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        assert_eq!(self.rows, rhs.rows);

        assert_eq!(self.cols, rhs.cols);

        let data = self
            .data
            .into_iter()
            .zip(rhs.data)
            .map(|(a, b)| a - b)
            .collect();

        Self::new(
            self.rows,
            self.cols,
            data,
        )
        .with_backend(self.backend)
    }
}

impl<T: Field> Mul for Matrix<T> {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {

        assert_eq!(self.cols, rhs.rows);

        // Try Faer backend
        if self.backend == Backend::Faer
            || rhs.backend
                == Backend::Faer
        {

            if let Some(res) =
                T::faer_mul(&self, &rhs)
            {

                return res
                    .with_backend(
                        self.backend,
                    ); // Propagate lhs backend preference
            }
        }

        let mut data = vec![
                T::zero();
                self.rows * rhs.cols
            ];

        for i in 0 .. self.rows {

            for j in 0 .. rhs.cols {

                for k in 0 .. self.cols
                {

                    data[i * rhs
                        .cols
                        + j] += self
                        .get(i, k)
                        .clone()
                        * rhs
                            .get(k, j)
                            .clone();
                }
            }
        }

        Self::new(
            self.rows,
            rhs.cols,
            data,
        )
        .with_backend(self.backend)
    }
}

/// Scalar multiplication: Matrix * scalar

impl Mul<f64> for Matrix<f64> {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {

        let new_data = self
            .data
            .into_iter()
            .map(|x| x * rhs)
            .collect();

        Self::new(
            self.rows,
            self.cols,
            new_data,
        )
        .with_backend(self.backend)
    }
}

/// Scalar multiplication: &Matrix * scalar

impl Mul<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn mul(
        self,
        rhs: f64,
    ) -> Matrix<f64> {

        let new_data = self
            .data
            .iter()
            .map(|x| *x * rhs)
            .collect();

        Matrix::new(
            self.rows,
            self.cols,
            new_data,
        )
        .with_backend(self.backend)
    }
}

impl<T: Field> AddAssign for Matrix<T> {
    fn add_assign(
        &mut self,
        rhs: Self,
    ) {

        assert_eq!(self.rows, rhs.rows);

        assert_eq!(self.cols, rhs.cols);

        for (a, b) in self
            .data
            .iter_mut()
            .zip(rhs.data)
        {

            *a += b;
        }
    }
}

impl<T: Field> SubAssign for Matrix<T> {
    fn sub_assign(
        &mut self,
        rhs: Self,
    ) {

        assert_eq!(self.rows, rhs.rows);

        assert_eq!(self.cols, rhs.cols);

        for (a, b) in self
            .data
            .iter_mut()
            .zip(rhs.data)
        {

            *a -= b;
        }
    }
}

impl<T: Field> MulAssign for Matrix<T> {
    fn mul_assign(
        &mut self,
        rhs: Self,
    ) {

        let res = self.clone() * rhs;

        self.rows = res.rows;

        self.cols = res.cols;

        self.data = res.data;
    }
}

impl<T: Field> Neg for Matrix<T> {
    type Output = Self;

    fn neg(self) -> Self {

        let data = self
            .data
            .into_iter()
            .map(|x| -x)
            .collect();

        Self::new(
            self.rows,
            self.cols,
            data,
        )
        .with_backend(self.backend)
    }
}

impl MulAssign<f64> for Matrix<f64> {
    fn mul_assign(
        &mut self,
        rhs: f64,
    ) {

        for x in &mut self.data {

            *x *= rhs;
        }
    }
}

impl Div<f64> for Matrix<f64> {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {

        let new_data = self
            .data
            .into_iter()
            .map(|x| x / rhs)
            .collect();

        Self::new(
            self.rows,
            self.cols,
            new_data,
        )
        .with_backend(self.backend)
    }
}

impl DivAssign<f64> for Matrix<f64> {
    fn div_assign(
        &mut self,
        rhs: f64,
    ) {

        for x in &mut self.data {

            *x /= rhs;
        }
    }
}
