//! # Numerical Input/Output Utilities
//!
//! This module provides input/output utilities for numerical data, primarily focusing
//! on reading and writing `ndarray` arrays to/from `.npy` files. It also includes
//! functions to convert between `Expr::Matrix` and `ndarray::Array2<f64>` for seamless
//! integration with symbolic and numerical computations.

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

use ndarray::Array2;
#[cfg(feature = "output")]
use ndarray_npy::read_npy;
#[cfg(feature = "output")]
use ndarray_npy::write_npy;
use serde_json;

use crate::prelude::Expr;

/// Writes a 2D `ndarray::Array` to a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
/// * `arr` - The array to write.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.
///
/// # Panics
/// Panics if the write fails.

#[cfg(feature = "output")]

pub fn write_npy_file<
    P: AsRef<Path>,
>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {

    write_npy(filename, arr)
        .map_err(|e| e.to_string())
}

#[cfg(not(feature = "output"))]

pub fn write_npy_file<
    P: AsRef<Path>,
>(
    _filename: P,
    _arr: &Array2<f64>,
) -> Result<(), String> {

    Err(
        "Feature 'output' is required \
         for .npy support"
            .to_string(),
    )
}

/// Reads a 2D `ndarray::Array` from a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
///
/// # Returns
/// The read array as an `ndarray::Array2<f64>`.
///
/// # Errors
///
/// This function will return an error if the file cannot be read.
///
/// # Panics
/// Panics if the read fails.

#[cfg(feature = "output")]

pub fn read_npy_file<P: AsRef<Path>>(
    filename: P
) -> Result<Array2<f64>, String> {

    read_npy(filename)
        .map_err(|e| e.to_string())
}

#[cfg(not(feature = "output"))]

pub fn read_npy_file<P: AsRef<Path>>(
    _filename: P
) -> Result<Array2<f64>, String> {

    Err(
        "Feature 'output' is required \
         for .npy support"
            .to_string(),
    )
}

/// Writes a 2D `ndarray::Array` to a CSV file.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.

pub fn write_csv_file<
    P: AsRef<Path>,
>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {

    let mut file =
        File::create(filename)
            .map_err(|e| {

                e.to_string()
            })?;

    for row in arr.outer_iter() {

        let line = row
            .iter()
            .map(|&v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");

        writeln!(file, "{line}")
            .map_err(|e| {

                e.to_string()
            })?;
    }

    Ok(())
}

/// Reads a 2D `ndarray::Array` from a CSV file.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.

pub fn read_csv_file<P: AsRef<Path>>(
    filename: P
) -> Result<Array2<f64>, String> {

    let file = File::open(filename)
        .map_err(|e| e.to_string())?;

    let reader = BufReader::new(file);

    let mut data = Vec::new();

    let mut rows = 0;

    let mut cols = 0;

    for line_res in reader.lines() {

        let line =
            line_res.map_err(|e| {

                e.to_string()
            })?;

        let line = line.trim();

        if line.is_empty() {

            continue;
        }

        let row_data: Vec<f64> =
            line.split(',')
                .map(|s| {

                    s.trim()
                        .parse::<f64>()
                        .map_err(|e| {

                            e.to_string(
                            )
                        })
                })
                .collect::<Result<
                    Vec<_>,
                    String,
                >>()?;

        if rows == 0 {

            cols = row_data.len();
        } else if row_data.len() != cols
        {

            return Err(
                "Inconsistent column \
                 count in CSV"
                    .to_string(),
            );
        }

        data.extend(row_data);

        rows += 1;
    }

    if rows == 0 {

        return Ok(Array2::zeros((
            0, 0,
        )));
    }

    Array2::from_shape_vec(
        (rows, cols),
        data,
    )
    .map_err(|e| e.to_string())
}

/// Writes a 2D `ndarray::Array` to a JSON file.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.

pub fn write_json_file<
    P: AsRef<Path>,
>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {

    let file = File::create(filename)
        .map_err(|e| {

        e.to_string()
    })?;

    serde_json::to_writer_pretty(
        file, arr,
    )
    .map_err(|e| e.to_string())
}

/// Reads a 2D `ndarray::Array` from a JSON file.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.

pub fn read_json_file<
    P: AsRef<Path>,
>(
    filename: P
) -> Result<Array2<f64>, String> {

    let file = File::open(filename)
        .map_err(|e| e.to_string())?;

    let reader = BufReader::new(file);

    serde_json::from_reader(reader)
        .map_err(|e| e.to_string())
}

#[cfg(test)]

mod tests {

    use std::fs;

    use ndarray::array;
    use proptest::prelude::*;

    use super::*;

    #[test]

    fn test_write_read_npy() {

        #[cfg(feature = "output")]
        {

            let arr = array![
                [1.0, 2.0],
                [3.0, 4.0]
            ];

            let filename =
                "test_array.npy";

            let _ = write_npy_file(
                filename,
                &arr,
            )
            .unwrap();

            let read_arr =
                read_npy_file(filename)
                    .unwrap();

            assert_eq!(arr, read_arr);

            let _ = fs::remove_file(
                filename,
            );
        }
    }

    #[test]

    fn test_write_read_csv() {

        let arr = array![
            [1.1, 2.2],
            [3.3, 4.4]
        ];

        let filename = "test_array.csv";

        write_csv_file(filename, &arr)
            .unwrap();

        let read_arr =
            read_csv_file(filename)
                .unwrap();

        assert_eq!(
            arr.shape(),
            read_arr.shape()
        );

        for (a, b) in arr
            .iter()
            .zip(read_arr.iter())
        {

            assert!(
                (a - b).abs() < 1e-10
            );
        }

        let _ =
            fs::remove_file(filename);
    }

    #[test]

    fn test_write_read_json() {

        let arr = array![
            [1.0, 2.0],
            [3.0, 4.0]
        ];

        let filename =
            "test_array.json";

        write_json_file(filename, &arr)
            .unwrap();

        let read_arr =
            read_json_file(filename)
                .unwrap();

        assert_eq!(
            arr.shape(),
            read_arr.shape()
        );

        for (a, b) in arr
            .iter()
            .zip(read_arr.iter())
        {

            assert!(
                (a - b).abs() < 1e-10
            );
        }

        let _ =
            fs::remove_file(filename);
    }

    proptest! {
        #[test]
        fn prop_csv_roundtrip(
            rows in 1..20usize,
            cols in 1..20usize,
            data in prop::collection::vec(0.1..100.0f64, 1..400)
        ) {
            prop_assume!(data.len() >= rows * cols);
            let actual_data = data.iter().take(rows * cols).cloned().collect::<Vec<_>>();
            let arr = Array2::from_shape_vec((rows, cols), actual_data).unwrap();
            let filename = "test_prop_csv.csv";
            write_csv_file(filename, &arr).unwrap();
            let read_arr = read_csv_file(filename).unwrap();

            assert_eq!(arr.shape(), read_arr.shape());
            for (a, b) in arr.iter().zip(read_arr.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
            let _ = fs::remove_file(filename);
        }

        #[test]
        fn prop_json_roundtrip(
            rows in 1..20usize,
            cols in 1..20usize,
            data in prop::collection::vec(0.1..100.0f64, 1..400)
        ) {
            prop_assume!(data.len() >= rows * cols);
            let actual_data = data.iter().take(rows * cols).cloned().collect::<Vec<_>>();
            let arr = Array2::from_shape_vec((rows, cols), actual_data).unwrap();
            let filename = "test_prop_json.json";
            write_json_file(filename, &arr).unwrap();
            let read_arr = read_json_file(filename).unwrap();

            assert_eq!(arr.shape(), read_arr.shape());
            for (a, b) in arr.iter().zip(read_arr.iter()) {
                assert!((a - b).abs() < 1e-10);
            }
            let _ = fs::remove_file(filename);
        }
    }

    #[test]

    fn test_dispatchers() {

        let arr = array![
            [1.0, 2.0],
            [3.0, 4.0]
        ];

        let expr = Expr::Matrix(vec![
            vec![
                Expr::Constant(1.0),
                Expr::Constant(2.0),
            ],
            vec![
                Expr::Constant(3.0),
                Expr::Constant(4.0),
            ],
        ]);

        // CSV
        let csv_file =
            "test_dispatch.csv";

        save_expr(csv_file, &expr)
            .unwrap();

        let loaded_csv =
            load_expr(csv_file)
                .unwrap();

        if let Expr::Matrix(rows) =
            loaded_csv
        {

            assert_eq!(rows.len(), 2);

            assert_eq!(
                rows[0][0]
                    .to_f64()
                    .unwrap(),
                1.0
            );
        }

        let _ =
            fs::remove_file(csv_file);

        // JSON
        let json_file =
            "test_dispatch.json";

        save_expr(json_file, &expr)
            .unwrap();

        let loaded_json =
            load_expr(json_file)
                .unwrap();

        assert_eq!(expr, loaded_json);

        let _ =
            fs::remove_file(json_file);

        #[cfg(feature = "output")]
        {

            // NPY
            let npy_file =
                "test_dispatch.npy";

            save_expr(npy_file, &expr)
                .unwrap();

            let loaded_npy =
                load_expr(npy_file)
                    .unwrap();

            assert_eq!(
                expr,
                loaded_npy
            );

            let _ = fs::remove_file(
                npy_file,
            );
        }
    }
}

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
///
/// # Errors
///
/// This function will return an error if the input expression is not a matrix,
/// contains non-numerical elements, or if the file write operation fails.

pub fn save_expr_as_npy<
    P: AsRef<Path>,
>(
    path: P,
    matrix_expr: &Expr,
) -> Result<(), String> {

    if let Expr::Matrix(rows) =
        matrix_expr
    {

        if rows.is_empty() {

            let arr: Array2<f64> =
                Array2::zeros((0, 0));

            write_npy_file(path, &arr)?;

            return Ok(());
        }

        let num_rows = rows.len();

        let num_cols = rows[0].len();

        let mut arr = Array2::zeros((
            num_rows,
            num_cols,
        ));

        for (i, row) in rows
            .iter()
            .enumerate()
        {

            if row.len() != num_cols {

                return Err(
                    "All rows must \
                     have the same \
                     number of columns"
                        .to_string(),
                );
            }

            for (j, elem) in row
                .iter()
                .enumerate()
            {

                let val = elem
                    .to_f64()
                    .ok_or_else(|| {

                        format!(
                            "Matrix element at ({i},{j}) is not a number"
                        )
                    })?;

                arr[[i, j]] = val;
            }
        }

        write_npy_file(path, &arr)?;

        Ok(())
    } else {

        Err(
            "Input expression is not \
             a matrix"
                .to_string(),
        )
    }
}

/// Converts an `Expr::Matrix` to an `ndarray::Array2<f64>` and saves it as a CSV file.
///
/// # Errors
///
/// This function will return an error if the expression is not a matrix or
/// if the file cannot be written to.

pub fn save_expr_as_csv<
    P: AsRef<Path>,
>(
    path: P,
    matrix_expr: &Expr,
) -> Result<(), String> {

    if let Expr::Matrix(rows) =
        matrix_expr
    {

        if rows.is_empty() {

            let arr: Array2<f64> =
                Array2::zeros((0, 0));

            return write_csv_file(
                path, &arr,
            );
        }

        let num_rows = rows.len();

        let num_cols = rows[0].len();

        let mut arr = Array2::zeros((
            num_rows,
            num_cols,
        ));

        for (i, row) in rows
            .iter()
            .enumerate()
        {

            if row.len() != num_cols {

                return Err(
                    "Inconsistent \
                     columns"
                        .to_string(),
                );
            }

            for (j, elem) in row
                .iter()
                .enumerate()
            {

                arr[[i, j]] = elem
                    .to_f64()
                    .ok_or(
                        "Not a number",
                    )?;
            }
        }

        write_csv_file(path, &arr)
    } else {

        Err("Not a matrix".to_string())
    }
}

/// Reads a CSV file into an `ndarray::Array2<f64>` and converts it to an `Expr::Matrix`.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.

pub fn load_csv_as_expr<
    P: AsRef<Path>,
>(
    path: P
) -> Result<Expr, String> {

    let arr = read_csv_file(path)?;

    let mut rows = Vec::new();

    for row in arr.outer_iter() {

        rows.push(
            row.iter()
                .map(|&v| {

                    Expr::Constant(v)
                })
                .collect(),
        );
    }

    Ok(Expr::Matrix(rows))
}

/// Converts an `Expr::Matrix` to an `ndarray::Array2<f64>` and saves it as a JSON file.
///
/// # Errors
///
/// This function will return an error if the expression is not a matrix or
/// if the file cannot be written to.

pub fn save_expr_as_json<
    P: AsRef<Path>,
>(
    path: P,
    matrix_expr: &Expr,
) -> Result<(), String> {

    if let Expr::Matrix(rows) =
        matrix_expr
    {

        if rows.is_empty() {

            let arr: Array2<f64> =
                Array2::zeros((0, 0));

            return write_json_file(
                path, &arr,
            );
        }

        let num_rows = rows.len();

        let num_cols = rows[0].len();

        let mut arr = Array2::zeros((
            num_rows,
            num_cols,
        ));

        for (i, row) in rows
            .iter()
            .enumerate()
        {

            if row.len() != num_cols {

                return Err(
                    "Inconsistent \
                     columns"
                        .to_string(),
                );
            }

            for (j, elem) in row
                .iter()
                .enumerate()
            {

                arr[[i, j]] = elem
                    .to_f64()
                    .ok_or(
                        "Not a number",
                    )?;
            }
        }

        write_json_file(path, &arr)
    } else {

        Err("Not a matrix".to_string())
    }
}

/// Reads a JSON file into an `ndarray::Array2<f64>` and converts it to an `Expr::Matrix`.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.

pub fn load_json_as_expr<
    P: AsRef<Path>,
>(
    path: P
) -> Result<Expr, String> {

    let arr = read_json_file(path)?;

    let mut rows = Vec::new();

    for row in arr.outer_iter() {

        rows.push(
            row.iter()
                .map(|&v| {

                    Expr::Constant(v)
                })
                .collect(),
        );
    }

    Ok(Expr::Matrix(rows))
}

/// Reads a `.npy` file into an `ndarray::Array2<f64>` and converts it to an `Expr::Matrix`.
///
/// # Arguments
/// * `path` - The path to the `.npy` file.
///
/// # Returns
/// A `Result` containing the `Expr::Matrix` representation of the loaded array,
/// or an error string if the read fails.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or if the
/// data cannot be converted to an `Expr::Matrix`.

pub fn load_npy_as_expr<
    P: AsRef<Path>,
>(
    path: P
) -> Result<Expr, String> {

    let arr = read_npy_file(path)?;

    let mut rows = Vec::new();

    for row in arr.outer_iter() {

        let mut expr_row = Vec::new();

        for val in row {

            expr_row.push(
                Expr::Constant(*val),
            );
        }

        rows.push(expr_row);
    }

    Ok(Expr::Matrix(rows))
}

/// Automatically dispatches to the correct saver based on the file extension.
///
/// # Errors
///
/// This function will return an error if the file extension is not supported
/// or if the save operation fails.

pub fn save_expr<P: AsRef<Path>>(
    path: P,
    expr: &Expr,
) -> Result<(), String> {

    let path_ref = path.as_ref();

    let extension = path_ref
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_lowercase)
        .unwrap_or_default();

    match extension.as_str() {
        | "npy" => {
            save_expr_as_npy(path, expr)
        },
        | "csv" => {
            save_expr_as_csv(path, expr)
        },
        | "json" => {
            save_expr_as_json(
                path, expr,
            )
        },
        | _ => {
            Err(format!(
                "Unsupported file \
                 format: .{extension}"
            ))
        },
    }
}

/// Automatically dispatches to the correct loader based on the file extension.
///
/// # Errors
///
/// This function will return an error if the file extension is not supported
/// or if the load operation fails.

pub fn load_expr<P: AsRef<Path>>(
    path: P
) -> Result<Expr, String> {

    let path_ref = path.as_ref();

    let extension = path_ref
        .extension()
        .and_then(|s| s.to_str())
        .map(str::to_lowercase)
        .unwrap_or_default();

    match extension.as_str() {
        | "npy" => {
            load_npy_as_expr(path)
        },
        | "csv" => {
            load_csv_as_expr(path)
        },
        | "json" => {
            load_json_as_expr(path)
        },
        | _ => {
            Err(format!(
                "Unsupported file \
                 format: .{extension}"
            ))
        },
    }
}
