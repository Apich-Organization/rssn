//! # Numerical Matrix and Linear Algebra
//!
//! This module provides a generic `Matrix` struct for dense matrices over any type
//! that implements a custom `Field` trait. It supports a wide range of linear algebra
//! operations, including matrix arithmetic, RREF, inversion, null space calculation,
//! and eigenvalue decomposition for symmetric matrices.
use crate::symbolic::finite_field::PrimeFieldElement;
use num_traits::{One, Zero};
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
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
    fn is_invertible(&self) -> bool;
    fn inverse(&self) -> Result<Self, String>;
}
impl Field for f64 {
    fn is_invertible(&self) -> bool {
        *self != 0.0
    }
    fn inverse(&self) -> Result<Self, String> {
        if self.is_invertible() {
            Ok(1.0 / self)
        } else {
            Err("Cannot invert 0.0".to_string())
        }
    }
}
impl Field for PrimeFieldElement {
    fn is_invertible(&self) -> bool {
        !self.value.is_zero()
    }
    fn inverse(&self) -> Result<Self, String> {
        self.inverse()
            .ok_or_else(|| "Cannot invert non-invertible element".to_string())
    }
}
/// A generic dense matrix over any type that implements the Field trait.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<T: Field> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
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
    pub fn new(rows: usize, cols: usize, data: Vec<T>) -> Self {
        assert_eq!(rows * cols, data.len());
        Matrix { rows, cols, data }
    }
    /// Creates a new `Matrix` filled with the zero element of type `T`.
    ///
    /// # Arguments
    /// * `rows` - The number of rows.
    /// * `cols` - The number of columns.
    ///
    /// # Returns
    /// A new `Matrix` of the specified dimensions, with all elements initialized to `T::zero()`.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Matrix {
            rows,
            cols,
            data: vec![T::zero(); rows * cols],
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
    pub fn get(&self, row: usize, col: usize) -> &T {
        &self.data[row * self.cols + col]
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
    pub fn get_mut(&mut self, row: usize, col: usize) -> &mut T {
        &mut self.data[row * self.cols + col]
    }
    /// Returns the number of rows in the matrix.
    pub fn rows(&self) -> usize {
        self.rows
    }
    /// Returns the number of columns in the matrix.
    pub fn cols(&self) -> usize {
        self.cols
    }
    /// Returns an immutable reference to the matrix's internal data vector.
    pub fn data(&self) -> &Vec<T> {
        &self.data
    }
    /// Returns a `Vec` of `Vec<T>` where each inner `Vec` represents a column of the matrix.
    ///
    /// This method effectively transposes the matrix data into a column-major representation.
    ///
    /// # Returns
    /// A `Vec<Vec<T>>` where each inner vector is a column.
    pub fn get_cols(&self) -> Vec<Vec<T>> {
        let mut cols_vec = Vec::with_capacity(self.cols);
        for j in 0..self.cols {
            let mut col = Vec::with_capacity(self.rows);
            for i in 0..self.rows {
                col.push(self.get(i, j).clone());
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
    pub fn rref(&mut self) -> Result<usize, String> {
        let mut pivot_row = 0;
        for j in 0..self.cols {
            if pivot_row >= self.rows {
                break;
            }
            let mut i = pivot_row;
            while i < self.rows && !self.get(i, j).is_invertible() {
                i += 1;
            }
            if i < self.rows {
                self.data.swap(i * self.cols, pivot_row * self.cols);
                let pivot_inv = self.get(pivot_row, j).clone().inverse()?;
                for k in j..self.cols {
                    let val = self.get(pivot_row, k).clone();
                    *self.get_mut(pivot_row, k) = val * pivot_inv.clone();
                }
                for i_prime in 0..self.rows {
                    if i_prime != pivot_row {
                        let factor = self.get(i_prime, j).clone();
                        for k in j..self.cols {
                            let pivot_row_val = self.get(pivot_row, k).clone();
                            let term = factor.clone() * pivot_row_val;
                            let current_val = self.get(i_prime, k).clone();
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
        let mut new_data = vec![T::zero(); self.rows * self.cols];
        for i in 0..self.rows {
            for j in 0..self.cols {
                new_data[j * self.rows + i] = self.get(i, j).clone();
            }
        }
        Matrix::new(self.cols, self.rows, new_data)
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
    /// A new `Matrix` representing the product of the two matrices, or an error if multiplication is not possible.
    pub fn mul_strassen(&self, other: &Self) -> Result<Self, String> {
        if self.cols != other.rows {
            return Err(format!(
                "Matrix multiplication not possible: left.cols ({}) != right.rows ({})",
                self.cols, other.rows
            ));
        }

        let n = self.rows.max(self.cols).max(other.rows).max(other.cols);
        let m = n.next_power_of_two();
        let mut a_padded = Matrix::zeros(m, m);
        let mut b_padded = Matrix::zeros(m, m);

        // Fill the padded matrices
        for i in 0..self.rows {
            for j in 0..self.cols {
                *a_padded.get_mut(i, j) = self.get(i, j).clone();
            }
        }
        for i in 0..other.rows {
            for j in 0..other.cols {
                *b_padded.get_mut(i, j) = other.get(i, j).clone();
            }
        }

        let c_padded = strassen_recursive(&a_padded, &b_padded);

        // Check if the result matrix has the expected size
        if c_padded.rows != m || c_padded.cols != m {
            return Err("Internal error: Strassen result has unexpected dimensions".to_string());
        }

        let mut result_data = vec![T::zero(); self.rows * other.cols];
        for i in 0..self.rows {
            for j in 0..other.cols {
                result_data[i * other.cols + j] = c_padded.get(i, j).clone();
            }
        }

        Ok(Matrix::new(self.rows, other.cols, result_data))
    }
    /// Splits a matrix into four sub-matrices of equal size.
    fn split(&self) -> (Self, Self, Self, Self) {
        // Check that the matrix dimensions are even so we can split it equally
        if self.rows % 2 != 0 || self.cols % 2 != 0 {
            // Return appropriately sized zero matrices if splitting isn't possible
            let half_rows = self.rows / 2;
            let half_cols = self.cols / 2;
            return (
                Matrix::zeros(half_rows, half_cols),
                Matrix::zeros(half_rows, half_cols),
                Matrix::zeros(half_rows, half_cols),
                Matrix::zeros(half_rows, half_cols),
            );
        }

        let new_dim = self.rows / 2;
        let new_col_dim = self.cols / 2; // Make sure to use the correct column division

        let mut a11 = Matrix::zeros(new_dim, new_col_dim);
        let mut a12 = Matrix::zeros(new_dim, new_col_dim);
        let mut a21 = Matrix::zeros(new_dim, new_col_dim);
        let mut a22 = Matrix::zeros(new_dim, new_col_dim);

        for i in 0..new_dim {
            for j in 0..new_col_dim {
                *a11.get_mut(i, j) = self.get(i, j).clone();
                *a12.get_mut(i, j) = self.get(i, j + new_col_dim).clone();
                *a21.get_mut(i, j) = self.get(i + new_dim, j).clone();
                *a22.get_mut(i, j) = self.get(i + new_dim, j + new_col_dim).clone();
            }
        }
        (a11, a12, a21, a22)
    }
    /// Joins four sub-matrices into a single larger matrix.
    fn join(a11: &Self, a12: &Self, a21: &Self, a22: &Self) -> Self {
        // All four submatrices should have the same dimensions
        if a11.rows != a12.rows
            || a11.cols != a12.cols
            || a11.rows != a21.rows
            || a11.cols != a21.cols
            || a11.rows != a22.rows
            || a11.cols != a22.cols
        {
            // Return a zero matrix if dimensions don't match
            return Matrix::zeros(0, 0);
        }

        let result_rows = a11.rows * 2;
        let result_cols = a11.cols * 2;
        let mut result = Matrix::zeros(result_rows, result_cols);

        let sub_row_dim = a11.rows;
        let sub_col_dim = a11.cols;

        // Copy a11 to top-left
        for i in 0..sub_row_dim {
            for j in 0..sub_col_dim {
                *result.get_mut(i, j) = a11.get(i, j).clone();
            }
        }

        // Copy a12 to top-right
        for i in 0..sub_row_dim {
            for j in 0..sub_col_dim {
                *result.get_mut(i, j + sub_col_dim) = a12.get(i, j).clone();
            }
        }

        // Copy a21 to bottom-left
        for i in 0..sub_row_dim {
            for j in 0..sub_col_dim {
                *result.get_mut(i + sub_row_dim, j) = a21.get(i, j).clone();
            }
        }

        // Copy a22 to bottom-right
        for i in 0..sub_row_dim {
            for j in 0..sub_col_dim {
                *result.get_mut(i + sub_row_dim, j + sub_col_dim) = a22.get(i, j).clone();
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
    /// The determinant of the matrix, or an error if the matrix is not square.
    pub fn determinant(&self) -> Result<T, String> {
        if self.rows != self.cols {
            return Err("Matrix must be square to compute the determinant.".to_string());
        }
        if self.rows > 64 && self.rows.is_multiple_of(2) {
            return self.determinant_block();
        }
        if self.rows == 0 {
            return Ok(T::one());
        }
        if self.rows == 1 {
            return Ok(self.get(0, 0).clone());
        }
        if self.rows == 2 {
            let a = self.get(0, 0).clone();
            let b = self.get(0, 1).clone();
            let c = self.get(1, 0).clone();
            let d = self.get(1, 1).clone();
            return Ok(a * d - b * c);
        }
        let (lu, pivots) = self.lu_decomposition()?;
        let mut det = T::one();
        for i in 0..self.rows {
            det *= lu.get(i, i).clone();
        }
        if (pivots.len() % 2) != 0 {
            det = -det;
        }
        Ok(det)
    }
    /// Computes the LU decomposition of a square matrix.
    ///
    /// # Returns
    /// A tuple containing the LU matrix and the number of permutations.
    pub fn lu_decomposition(&self) -> Result<(Matrix<T>, Vec<usize>), String> {
        if self.rows != self.cols {
            return Err("Matrix must be square for LU decomposition.".to_string());
        }
        let n = self.rows;
        let mut lu = self.clone();
        let mut pivots = (0..n).collect::<Vec<usize>>();
        for j in 0..n {
            let mut pivot_row = j;
            for i in j..n {
                if self.get(i, j).is_invertible() {
                    pivot_row = i;
                    break;
                }
            }
            if pivot_row != j {
                pivots.swap(j, pivot_row);
                for k in 0..n {
                    let val1 = lu.get(j, k).clone();
                    let val2 = lu.get(pivot_row, k).clone();
                    *lu.get_mut(j, k) = val2;
                    *lu.get_mut(pivot_row, k) = val1;
                }
            }
            let pivot_val = lu.get(j, j).clone();
            if !pivot_val.is_invertible() {
                return Err("Matrix is singular.".to_string());
            }
            for i in (j + 1)..n {
                let factor = lu.get(i, j).clone() * pivot_val.inverse()?;
                *lu.get_mut(i, j) = factor.clone();
                for k in (j + 1)..n {
                    let val = lu.get(j, k).clone() * factor.clone();
                    let current_val = lu.get(i, k).clone();
                    *lu.get_mut(i, k) = current_val - val;
                }
            }
        }
        Ok((lu, pivots))
    }
    /// # Block Matrix Determinant
    ///
    /// Computes the determinant using block matrix decomposition (Schur complement).
    /// This is efficient for large matrices.
    pub fn determinant_block(&self) -> Result<T, String> {
        if self.rows != self.cols {
            return Err("Matrix must be square.".to_string());
        }
        let n = self.rows;
        if n == 0 {
            return Ok(T::one());
        }
        if !n.is_multiple_of(2) {
            return self.determinant_lu();
        }

        let (a, b, c, d) = self.split();

        // Check that the submatrices have consistent dimensions for block operations
        if a.rows != b.rows || a.cols != c.rows || b.cols != d.cols || c.cols != d.rows {
            return Err(
                "Block matrix decomposition failed due to inconsistent submatrix dimensions."
                    .to_string(),
            );
        }

        // Try to compute determinant using Schur complement if the top-left block is invertible
        match a.inverse() {
            Some(a_inv) => {
                // Calculate Schur complement: S = D - C * A^(-1) * B
                let a_inv_b = a_inv * b.clone();

                let schur_complement = d.clone() - c.clone() * a_inv_b;

                match (a.determinant_lu(), schur_complement.determinant_lu()) {
                    (Ok(det_a), Ok(det_s)) => Ok(det_a * det_s),
                    _ => {
                        // If calculation fails, fallback to LU decomposition
                        self.determinant_lu()
                    }
                }
            }
            None => {
                // If A is not invertible, fallback to LU decomposition
                self.determinant_lu()
            }
        }
    }
    /// Computes the determinant using LU decomposition.
    pub fn determinant_lu(&self) -> Result<T, String> {
        let (lu, pivots) = self.lu_decomposition()?;
        let mut det = T::one();
        for i in 0..self.rows {
            det *= lu.get(i, i).clone();
        }
        if (pivots.len() % 2) != 0 {
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
    pub fn inverse(&self) -> Option<Self> {
        if self.rows != self.cols {
            return None;
        }
        let n = self.rows;
        let mut augmented = Matrix::zeros(n, 2 * n);
        for i in 0..n {
            for j in 0..n {
                *augmented.get_mut(i, j) = self.get(i, j).clone();
                if i == j {
                    *augmented.get_mut(i, j + n) = T::one();
                }
            }
        }

        // Check if RREF operation was successful
        if let Ok(rank) = augmented.rref() {
            if rank == n {
                let mut inv_data = vec![T::zero(); n * n];
                for i in 0..n {
                    for j in 0..n {
                        inv_data[i * n + j] = augmented.get(i, j + n).clone();
                    }
                }
                Some(Matrix::new(n, n, inv_data))
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
    pub fn null_space(&self) -> Result<Matrix<T>, String> {
        let mut rref_matrix = self.clone();
        let rank = rref_matrix.rref()?;
        let mut pivot_cols = Vec::new();
        let mut lead = 0;
        for r in 0..rank {
            if lead >= self.cols {
                break;
            }
            let mut i = lead;
            while !rref_matrix.get(r, i).is_invertible() {
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
        let free_cols: Vec<usize> = (0..self.cols).filter(|c| !pivot_cols.contains(c)).collect();
        let num_free = free_cols.len();
        let mut basis_vectors = Vec::with_capacity(num_free);
        for free_col in free_cols {
            let mut vec = vec![T::zero(); self.cols];
            vec[free_col] = T::one();
            for (i, &pivot_col) in pivot_cols.iter().enumerate() {
                vec[pivot_col] = -rref_matrix.get(i, free_col).clone();
            }
            basis_vectors.push(vec);
        }
        let mut null_space_data = vec![T::zero(); self.cols * num_free];
        for (j, basis_vec) in basis_vectors.iter().enumerate() {
            for (i, val) in basis_vec.iter().enumerate() {
                null_space_data[i * num_free + j] = val.clone();
            }
        }
        Ok(Matrix::new(self.cols, num_free, null_space_data))
    }
}
/// Recursive helper for Strassen's algorithm.
fn strassen_recursive<T: Field>(a: &Matrix<T>, b: &Matrix<T>) -> Matrix<T> {
    // Check that matrices can be multiplied
    if a.cols != b.rows {
        // Return a zero matrix as a fallback (though this shouldn't happen in practice with our padding)
        return Matrix::zeros(a.rows, b.cols);
    }

    let n = a.rows;
    if n <= 64 {
        // For small matrices, use standard multiplication to avoid overhead
        return a.clone() * b.clone();
    }

    // For the Strassen algorithm to work properly, matrix dimensions must be even
    if n % 2 != 0 {
        // Fall back to standard multiplication if dimension is odd
        return a.clone() * b.clone();
    }

    let (a11, a12, a21, a22) = a.split();
    let (b11, b12, b21, b22) = b.split();
    let p1 = strassen_recursive(&(a11.clone() + a22.clone()), &(b11.clone() + b22.clone()));
    let p2 = strassen_recursive(&(a21.clone() + a22.clone()), &b11);
    let p3 = strassen_recursive(&a11, &(b12.clone() - b22.clone()));
    let p4 = strassen_recursive(&a22, &(b21.clone() - b11.clone()));
    let p5 = strassen_recursive(&(a11.clone() + a12.clone()), &b22);
    let p6 = strassen_recursive(&(a21.clone() - a11.clone()), &(b11.clone() + b12.clone()));
    let p7 = strassen_recursive(&(a12.clone() - a22.clone()), &(b21.clone() + b22.clone()));
    let c11 = p1.clone() + p4.clone() - p5.clone() + p7.clone();
    let c12 = p3.clone() + p5.clone();
    let c21 = p2.clone() + p4.clone();
    let c22 = p1.clone() - p2.clone() + p3.clone() + p6.clone();
    Matrix::join(&c11, &c12, &c21, &c22)
}
impl Matrix<f64> {
    /// Creates an identity matrix of a given size.
    ///
    /// An identity matrix is a square matrix with ones on the main diagonal
    /// and zeros elsewhere. It acts as the multiplicative identity in matrix multiplication.
    ///
    /// # Arguments
    /// * `size` - The dimension of the square identity matrix.
    ///
    /// # Returns
    /// A new `Matrix<f64>` representing the identity matrix.
    pub fn identity(size: usize) -> Self {
        let mut data = vec![0.0; size * size];
        for i in 0..size {
            data[i * size + i] = 1.0;
        }
        Matrix::new(size, size, data)
    }
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
    /// and `eigenvectors` is a `Matrix<f64>` (columns are eigenvectors), or an error string if the matrix is not square.
    pub fn jacobi_eigen_decomposition(
        &self,
        max_sweeps: usize,
        tolerance: f64,
    ) -> Result<(Vec<f64>, Matrix<f64>), String> {
        if self.rows != self.cols {
            return Err("Matrix must be square.".to_string());
        }
        let mut a = self.clone();
        let n = self.rows;
        let mut eigenvectors = Matrix::identity(n);
        for _ in 0..max_sweeps {
            let mut off_diagonal_sum = 0.0;
            for p in 0..n {
                for q in (p + 1)..n {
                    off_diagonal_sum += a.get(p, q).abs().powi(2);
                }
            }
            if off_diagonal_sum.sqrt() < tolerance {
                break;
            }
            for p in 0..n {
                for q in (p + 1)..n {
                    let apq = *a.get(p, q);
                    if apq.abs() < tolerance / (n as f64) {
                        continue;
                    }
                    let app = *a.get(p, p);
                    let aqq = *a.get(q, q);
                    let tau = (aqq - app) / (2.0 * apq);
                    let t = if tau >= 0.0 {
                        1.0 / (tau + (1.0 + tau * tau).sqrt())
                    } else {
                        -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let new_app = app - t * apq;
                    let new_aqq = aqq + t * apq;
                    *a.get_mut(p, p) = new_app;
                    *a.get_mut(q, q) = new_aqq;
                    *a.get_mut(p, q) = 0.0;
                    *a.get_mut(q, p) = 0.0;
                    for i in 0..n {
                        if i != p && i != q {
                            let aip = *a.get(i, p);
                            let aiq = *a.get(i, q);
                            *a.get_mut(i, p) = c * aip - s * aiq;
                            *a.get_mut(i, q) = s * aip + c * aiq;
                            *a.get_mut(p, i) = *a.get(i, p);
                            *a.get_mut(q, i) = *a.get(i, q);
                        }
                    }
                    for i in 0..n {
                        let vip = *eigenvectors.get(i, p);
                        let viq = *eigenvectors.get(i, q);
                        *eigenvectors.get_mut(i, p) = c * vip - s * viq;
                        *eigenvectors.get_mut(i, q) = s * vip + c * viq;
                    }
                }
            }
        }
        let eigenvalues = (0..n).map(|i| *a.get(i, i)).collect();
        Ok((eigenvalues, eigenvectors))
    }
}
impl<T: Field> Add for Matrix<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);
        let data = self
            .data
            .into_iter()
            .zip(rhs.data)
            .map(|(a, b)| a + b)
            .collect();
        Matrix::new(self.rows, self.cols, data)
    }
}
impl<T: Field> Sub for Matrix<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);
        let data = self
            .data
            .into_iter()
            .zip(rhs.data)
            .map(|(a, b)| a - b)
            .collect();
        Matrix::new(self.rows, self.cols, data)
    }
}
impl<T: Field> Mul for Matrix<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.cols, rhs.rows);
        let mut data = vec![T::zero(); self.rows * rhs.cols];
        for i in 0..self.rows {
            for j in 0..rhs.cols {
                for k in 0..self.cols {
                    data[i * rhs.cols + j] += self.get(i, k).clone() * rhs.get(k, j).clone();
                }
            }
        }
        Matrix::new(self.rows, rhs.cols, data)
    }
}
/// Scalar multiplication: Matrix * scalar
impl Mul<f64> for Matrix<f64> {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let new_data = self.data.into_iter().map(|x| x * rhs).collect();
        Matrix::new(self.rows, self.cols, new_data)
    }
}
/// Scalar multiplication: &Matrix * scalar
impl Mul<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, rhs: f64) -> Matrix<f64> {
        let new_data = self.data.iter().map(|x| *x * rhs).collect();
        Matrix::new(self.rows, self.cols, new_data)
    }
}
