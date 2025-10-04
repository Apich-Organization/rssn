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
    fn inverse(&self) -> Self;
}

impl Field for f64 {
    fn is_invertible(&self) -> bool {
        *self != 0.0
    }
    fn inverse(&self) -> Self {
        1.0 / self
    }
}

impl Field for PrimeFieldElement {
    fn is_invertible(&self) -> bool {
        !self.value.is_zero()
    }
    fn inverse(&self) -> Self {
        self.inverse().unwrap()
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
    pub fn rref(&mut self) -> usize {
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

                let pivot_inv = self.get(pivot_row, j).clone().inverse();
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
        pivot_row
    }

    /// Computes the transpose of the matrix.
    ///
    /// The transpose of a matrix `A` (denoted `A^T`) is obtained by flipping the matrix
    /// over its diagonal; that is, it switches the row and column indices of the matrix.
    ///
    /// # Returns
    /// A new `Matrix` representing the transpose.
    pub fn transpose(&self) -> Self {
        let mut new_data = vec![T::zero(); self.rows * self.cols];
        for i in 0..self.rows {
            for j in 0..self.cols {
                new_data[j * self.rows + i] = self.get(i, j).clone();
            }
        }
        Matrix::new(self.cols, self.rows, new_data)
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

        if augmented.rref() != n {
            return None;
        } // Not invertible

        let mut inv_data = vec![T::zero(); n * n];
        for i in 0..n {
            for j in 0..n {
                inv_data[i * n + j] = augmented.get(i, j + n).clone();
            }
        }
        Some(Matrix::new(n, n, inv_data))
    }

    /// Computes a basis for the null space (kernel) of the matrix.
    ///
    /// The null space of a matrix `A` is the set of all vectors `x` such that `Ax = 0`.
    /// This method finds the null space by first computing the RREF of the matrix,
    /// identifying pivot and free variables, and then constructing basis vectors.
    ///
    /// # Returns
    /// A `Matrix` whose columns form a basis for the null space.
    pub fn null_space(&self) -> Matrix<T> {
        let mut rref_matrix = self.clone();
        let rank = rref_matrix.rref();

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

        // Transpose the basis vectors to form the columns of the null space matrix
        let mut null_space_data = vec![T::zero(); self.cols * num_free];
        for (j, basis_vec) in basis_vectors.iter().enumerate() {
            for (i, val) in basis_vec.iter().enumerate() {
                null_space_data[i * num_free + j] = val.clone();
            }
        }
        Matrix::new(self.cols, num_free, null_space_data)
    }
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
                break; // Converged
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

                    // Update matrix A = J^T * A * J
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

                    // Update eigenvector matrix V = V * J
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
impl<'a> Mul<f64> for &'a Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, rhs: f64) -> Matrix<f64> {
        let new_data = self.data.iter().map(|x| x.clone() * rhs).collect();
        Matrix::new(self.rows, self.cols, new_data)
    }
}
