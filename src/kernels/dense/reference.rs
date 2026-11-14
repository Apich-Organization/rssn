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
    for row in 0..rows {
        for column in 0..cols {
            if inner == 0 {
                break;
            }
            let mut total = a[row * inner] * b[column];
            for k in 1..inner {
                total += a[row * inner + k] * b[k * cols + column];
            }
            data[row * cols + column] = total;
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
