//! Sequential reference algebra and the stable first-product summation order.
use super::DenseError;
use super::Mat;

fn matrix(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<Mat, DenseError> {
    if rows.checked_mul(cols) != Some(data.len()) {
        return Err(DenseError::Shape);
    }
    if data.iter().any(|value| !value.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    Ok(Mat {
        rows,
        cols,
        data: data.to_vec(),
    })
}

/// Square determinant using the independent sequential LU reference kernel.
///
/// # Errors
/// Invalid shapes or non-finite entries.
pub fn det(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<f64, DenseError> {
    let matrix = matrix(rows, cols, data)?;
    match super::lu_factor(&matrix) {
        | Ok(lu) => Ok(lu.det()),
        | Err(DenseError::Singular) => Ok(0.0),
        | Err(error) => Err(error),
    }
}

/// Square inverse using the sequential LU reference kernel.
///
/// # Errors
/// Invalid shapes, non-finite entries, or a singular matrix.
pub fn inv(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<Vec<f64>, DenseError> {
    Ok(super::lu_factor(&matrix(rows, cols, data)?)?.inverse().data)
}

/// Square solve using the sequential LU reference kernel.
///
/// # Errors
/// Invalid shapes, non-finite data, or a singular matrix.
pub fn solve(
    rows: usize,
    cols: usize,
    data: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, DenseError> {
    if rhs.iter().any(|value| !value.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    super::solve(&matrix(rows, cols, data)?, rhs)
}

/// Row-major permutation, lower and upper factors, including zero pivots.
///
/// # Errors
/// Invalid shapes or non-finite entries.
#[allow(clippy::type_complexity)] // Public P, L, U factor layout.
pub fn lu(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), DenseError> {
    super::optimized::lu(rows, cols, data)
}

/// Householder QR with the historical thin Q zero-padded to square storage.
///
/// # Errors
/// Invalid shapes or non-finite entries.
pub fn qr(
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), DenseError> {
    rows.checked_mul(rows).ok_or(DenseError::Shape)?;
    let (mut q, r) = super::qr(&matrix(rows, cols, data)?);
    for row in 0..rows {
        for column in cols.min(rows)..rows {
            q.data[row * rows + column] = 0.0;
        }
    }
    Ok((q.data, r.data))
}

/// Matrix product with the first product as accumulator, preserving IEEE bits.
///
/// Inputs follow the historical unchecked reference API: lengths must match
/// `rows*inner` and `inner*cols`; a zero inner dimension produces positive zeros.
#[must_use]
pub fn matmul(
    rows: usize,
    inner: usize,
    cols: usize,
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    let mut data = vec![0.0; rows * cols];
    if data.is_empty() || inner == 0 {
        return data;
    }
    let left_rows = a[..rows * inner].chunks_exact(inner);
    let right = &b[..inner * cols];
    if cols <= 4 {
        // With only a few columns, keep each dot in a register instead of
        // updating the output row once per factor.
        for (left, output) in left_rows.zip(data.chunks_exact_mut(cols)) {
            for (column, slot) in output.iter_mut().enumerate() {
                let mut value = left[0] * right[column];
                for (k, factor) in left.iter().enumerate().skip(1) {
                    value += factor * right[k * cols + column];
                }
                *slot = value;
            }
        }
        return data;
    }
    let right_rows = right.chunks_exact(cols);
    for (left, output) in left_rows.zip(data.chunks_exact_mut(cols)) {
        for (value, &right) in output.iter_mut().zip(&right[..cols]) {
            *value = left[0] * right;
        }
        for (&factor, right) in left.iter().skip(1).zip(right_rows.clone().skip(1)) {
            // Visit adjacent columns while preserving the dot-product order
            // and the first product's signed zero for each cell.
            for (value, &element) in output.iter_mut().zip(right) {
                *value += factor * element;
            }
        }
    }
    data
}

/// The lower triangle of a matrix, the rest zeroed.
#[must_use]
pub fn tril(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; a.len()];
    for row in 0..rows {
        for column in 0..cols.min(row + 1) {
            out[row * cols + column] = a[row * cols + column];
        }
    }
    out
}

/// The upper triangle of a matrix, the rest zeroed.
#[must_use]
pub fn triu(
    rows: usize,
    cols: usize,
    a: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; a.len()];
    for row in 0..rows {
        for column in row..cols {
            out[row * cols + column] = a[row * cols + column];
        }
    }
    out
}
