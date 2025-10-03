//! # Numerical Input/Output Utilities
//!
//! This module provides input/output utilities for numerical data, primarily focusing
//! on reading and writing `ndarray` arrays to/from `.npy` files. It also includes
//! functions to convert between `Expr::Matrix` and `ndarray::Array2<f64>` for seamless
//! integration with symbolic and numerical computations.

use crate::prelude::Expr;
use ndarray::Array2;
use ndarray_npy::{read_npy, write_npy};
use std::path::Path;

/// Writes a 2D `ndarray::Array` to a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
/// * `arr` - The array to write.
///
/// # Panics
/// Panics if the write fails.
pub fn write_npy_file<P: AsRef<Path>>(filename: P, arr: &Array2<f64>) {
    write_npy(filename, arr).unwrap();
}

/// Reads a 2D `ndarray::Array` from a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
///
/// # Returns
/// The read array as an `ndarray::Array2<f64>`.
///
/// # Panics
/// Panics if the read fails.
pub fn read_npy_file<P: AsRef<Path>>(filename: P) -> Array2<f64> {
    read_npy(filename).unwrap()
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use std::fs;

    #[test]
    pub(crate) fn test_write_read_npy() {
        let arr = array![[1.0, 2.0], [3.0, 4.0]];
        let filename = "test_array.npy";

        write_npy_file(filename, &arr);

        let read_arr = read_npy_file(filename);

        assert_eq!(arr, read_arr);

        // Clean up the created file
        fs::remove_file(filename).unwrap();
    }
}
*/

/// Converts an `Expr::Matrix` to an `ndarray::Array2<f64>` and saves it as a `.npy` file.
///
/// This function acts as a bridge to the existing `ndarray-npy` functionality.
///
/// # Arguments
/// * `path` - The path to the `.npy` file.
/// * `matrix_expr` - The `Expr::Matrix` to save.
///
/// # Returns
/// A `Result` indicating success or an error string if the input is not a matrix
/// or contains non-numerical elements.
pub fn save_expr_as_npy<P: AsRef<Path>>(path: P, matrix_expr: &Expr) -> Result<(), String> {
    if let Expr::Matrix(rows) = matrix_expr {
        if rows.is_empty() {
            let arr: Array2<f64> = Array2::zeros((0, 0));
            write_npy_file(path, &arr);
            return Ok(());
        }
        let num_rows = rows.len();
        let num_cols = rows[0].len();
        let mut arr = Array2::zeros((num_rows, num_cols));

        for (i, row) in rows.iter().enumerate() {
            if row.len() != num_cols {
                return Err("All rows must have the same number of columns".to_string());
            }
            for (j, elem) in row.iter().enumerate() {
                let val = elem
                    .to_f64()
                    .ok_or_else(|| format!("Matrix element at ({},{}) is not a number", i, j))?;
                arr[[i, j]] = val;
            }
        }
        write_npy_file(path, &arr);
        Ok(())
    } else {
        Err("Input expression is not a matrix".to_string())
    }
}

/// Reads a `.npy` file into an `ndarray::Array2<f64>` and converts it to an `Expr::Matrix`.
///
/// # Arguments
/// * `path` - The path to the `.npy` file.
///
/// # Returns
/// A `Result` containing the `Expr::Matrix` representation of the loaded array,
/// or an error string if the read fails.
pub fn load_npy_as_expr<P: AsRef<Path>>(path: P) -> Result<Expr, String> {
    let arr = read_npy_file(path);
    let mut rows = Vec::new();
    for row in arr.outer_iter() {
        let mut expr_row = Vec::new();
        for val in row {
            expr_row.push(Expr::Constant(*val));
        }
        rows.push(expr_row);
    }
    Ok(Expr::Matrix(rows))
}
