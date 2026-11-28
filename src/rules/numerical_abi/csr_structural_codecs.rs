// Structural CSR operators share the iterative solver's validated wire format.

use crate::kernels::krylov;
use crate::kernels::numeric_common::scale_product;

record!(krylov::Csr {
    rows: usize,
    cols: usize,
    indptr: Vec<usize>,
    indices: Vec<usize>,
    values: Vec<f64>,
});

impl Arg for krylov::Csr {
    fn ty() -> TypeExpr {
        super::iterative::csr_type()
    }

    fn decode(
        value: &Value,
        _: usize,
    ) -> Result<Self, OpError> {
        super::iterative::decode_rectangular_csr(value)
    }
}

binding!(
    csr_transpose, "num_csr_transpose", "numerical_matrix",
    "Transposes a checked rectangular CSR record, preserving entry units and strictly ordered columns.",
    "crate::kernels::krylov::Csr::transpose", DimensionPolicy::CsrTranspose,
    (matrix: super::krylov::Csr) -> super::krylov::Csr,
    |matrix: super::krylov::Csr| { matrix.transpose() }
);

binding!(
    csr_dense_product, "num_csr_matmul_dense", "numerical_matrix",
    "Multiplies a rectangular CSR matrix by a dense rank-two tensor using SIMD row updates. The result keeps the product of both entry units.",
    "crate::kernels::krylov::Csr::spmm", DimensionPolicy::CsrDenseProduct,
    (matrix: super::krylov::Csr, rhs: dense::Mat) -> Result<dense::Mat, OpError>,
    |matrix: super::krylov::Csr, rhs: dense::Mat| {
        super::csr_dense_multiply(&matrix, &rhs)
    }
);

binding!(
    csr_sparse_product, "num_csr_matmul_sparse", "numerical_matrix",
    "Multiplies two rectangular CSR records. Output columns are strictly ordered and the value budget is checked before allocation. Entry dimensions multiply.",
    "crate::kernels::krylov::Csr::spgemm", DimensionPolicy::CsrSparseProduct,
    (left: super::krylov::Csr, right: super::krylov::Csr) -> Result<super::krylov::Csr, OpError>,
    |left: super::krylov::Csr, right: super::krylov::Csr| {
        super::csr_sparse_multiply(&left, &right)
    }
);

binding!(
    csr_lsqr, "num_csr_lsqr", "numerical_matrix",
    "Solves rectangular sparse least squares using cached working buffers and the original normal equations for convergence. Damping uses the matrix entry dimension. The reported residual is the unaugmented relative residual.",
    "crate::kernels::krylov::lsqr", DimensionPolicy::CsrLeastSquares,
    (matrix: super::krylov::Csr, rhs: Vec<f64>, damping: f64, tolerance: f64, max_iterations: usize) -> Result<super::krylov::IterResult, OpError>,
    |matrix: super::krylov::Csr, rhs: Vec<f64>, damping: f64, tolerance: f64, max_iterations: usize| {
        super::csr_least_squares(matrix, rhs, damping, tolerance, max_iterations)
    }
);

const CSR_MAX_WORK: usize = 50_000_000;
const CSR_VALUE_LIMIT: usize = super::native_object_api::VALUE_LIMIT;

fn csr_dense_multiply(
    matrix: &krylov::Csr,
    rhs: &dense::Mat,
) -> Result<dense::Mat, OpError> {
    if matrix.rows.saturating_mul(rhs.cols) >= CSR_VALUE_LIMIT
        || matrix.values.len().saturating_mul(rhs.cols) > CSR_MAX_WORK
    {
        return Err(OpError::Domain(
            "CSR dense product exceeds the result or work budget".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    matrix.spmm(rhs).map_err(OpError::Domain)
}

fn csr_sparse_multiply(
    left: &krylov::Csr,
    right: &krylov::Csr,
) -> Result<krylov::Csr, OpError> {
    if left.cols != right.rows {
        return Err(OpError::Domain(
            "CSR product inner dimensions differ".into(),
        ));
    }
    let mut work = 0usize;
    let mut upper = 0usize;
    for range in left.indptr.windows(2) {
        let row_work = left.indices[range[0]..range[1]]
            .iter()
            .fold(0usize, |total, &row| {
                total.saturating_add(right.indptr[row + 1] - right.indptr[row])
            });
        work = work.saturating_add(row_work);
        upper = upper.saturating_add(row_work.min(right.cols));
    }
    // The structural result is a five-field record and includes both index
    // and value tensors. Account for their elements and the CSR row pointers.
    if work > CSR_MAX_WORK
        || upper.saturating_mul(2).saturating_add(left.rows + 7) >= CSR_VALUE_LIMIT
    {
        return Err(OpError::Domain(
            "CSR sparse product exceeds the result or work budget".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    left.spgemm(right).map_err(OpError::Domain)
}

fn csr_least_squares(
    mut matrix: krylov::Csr,
    mut rhs: Vec<f64>,
    damping: f64,
    tolerance: f64,
    max_iterations: usize,
) -> Result<krylov::IterResult, OpError> {
    if matrix.rows == 0 || matrix.cols == 0 || matrix.rows != rhs.len() {
        return Err(OpError::Domain(
            "LSQR requires a nonempty rectangular CSR matrix and matching RHS".into(),
        ));
    }
    if damping < 0.0 || tolerance < 0.0 || max_iterations == 0 {
        return Err(OpError::Domain(
            "LSQR requires nonnegative damping and tolerance and a positive iteration limit".into(),
        ));
    }
    let per_iteration = matrix
        .values
        .len()
        .saturating_add(matrix.rows)
        .saturating_add(matrix.cols);
    if per_iteration.saturating_mul(max_iterations) > CSR_MAX_WORK {
        return Err(OpError::Domain(
            "LSQR exceeds the iteration work budget".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    let max_abs = |values: &[f64]| values.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let matrix_scale = max_abs(&matrix.values).max(damping);
    let matrix_scale = if matrix_scale > 0.0 {
        matrix_scale
    } else {
        1.0
    };
    let rhs_scale = max_abs(&rhs);
    let rhs_scale = if rhs_scale > 0.0 {
        rhs_scale
    } else {
        1.0
    };
    for entry in &mut matrix.values {
        *entry /= matrix_scale;
    }
    for entry in &mut rhs {
        *entry /= rhs_scale;
    }
    let transpose = matrix.transpose();
    let mut result = krylov::lsqr(
        |input, output| matrix.matvec(input, output),
        |input, output| transpose.matvec(input, output),
        matrix.rows,
        matrix.cols,
        &rhs,
        damping / matrix_scale,
        tolerance,
        max_iterations,
    );
    for coordinate in &mut result.x {
        *coordinate = scale_product(*coordinate, &[rhs_scale], &[matrix_scale])
            .map_err(|error| OpError::Domain(error.to_string()))?;
    }
    if !result.residual.is_finite() {
        return Err(OpError::Domain("LSQR returned a nonfinite residual".into()));
    }
    crate::cancellation::checkpoint()?;
    Ok(result)
}
