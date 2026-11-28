//! `f64` kernels for compressed sparse matrices (§8.4.1).
//!
//! Everything here takes the CSR arrays of a square-or-rectangular matrix
//! (`row_ptr`, `col_idx`, `values`) and dense right-hand sides. The direct
//! solvers convert CSR to CSC for `faer`'s sparse factorisations. The
//! coefficient matrix stays sparse, and a batch of right-hand sides shares
//! one factorisation. Iterative solvers access the matrix through products.

use crate::kernels::krylov;
use crate::kernels::numeric_common::NumericError;
use faer::linalg::solvers::Solve;
use faer::sparse::FaerError;
use faer::sparse::SparseColMat;
use faer::sparse::SparseRowMatRef;
use faer::sparse::SymbolicSparseRowMatRef;

/// A borrowed CSR matrix.
#[derive(Clone, Copy, Debug)]
pub struct Csr<'a> {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// `rows + 1` offsets into `col_idx`/`values`.
    pub row_ptr: &'a [usize],
    /// Column of each stored entry.
    pub col_idx: &'a [usize],
    /// Value of each stored entry.
    pub values: &'a [f64],
}

fn validate(a: &Csr<'_>) -> Result<(), NumericError> {
    if a.rows.checked_add(1) != Some(a.row_ptr.len())
        || a.col_idx.len() != a.values.len()
        || a.row_ptr.first() != Some(&0)
        || a.row_ptr.last() != Some(&a.values.len())
        || a.row_ptr.windows(2).any(|pair| pair[0] > pair[1])
        || a.col_idx.iter().any(|column| *column >= a.cols)
    {
        return Err(NumericError::ShapeMismatch("invalid CSR arrays".into()));
    }
    if a.values.iter().any(|value| !value.is_finite()) {
        return Err(NumericError::Domain(
            "sparse matrix contains non-finite entries".into(),
        ));
    }
    Ok(())
}

impl Csr<'_> {
    /// `y = A x`.
    ///
    /// # Errors
    /// [`NumericError::ShapeMismatch`] when `x` does not have `cols` entries.
    pub fn matvec(
        &self,
        x: &[f64],
    ) -> Result<Vec<f64>, NumericError> {
        validate(self)?;
        if x.len() != self.cols {
            return Err(NumericError::ShapeMismatch(format!(
                "matrix has {} columns, vector has {} entries",
                self.cols,
                x.len()
            )));
        }
        let mut y = vec![0.0; self.rows];
        self.matvec_validated(x, &mut y);
        Ok(y)
    }

    fn matvec_validated(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) {
        krylov::csr_matvec(self.rows, self.row_ptr, self.col_idx, self.values, x, y);
    }
}

/// Solves `A x = b` by sparse LU.
///
/// # Errors
/// [`NumericError::NotSquare`], [`NumericError::ShapeMismatch`] for a bad
/// right-hand side, and [`NumericError::Singular`] when the factorisation
/// fails or produces a non-finite solution. [`NumericError::Domain`] covers
/// non-finite matrix entries and failed allocation of the dense solution.
pub fn lu_solve(
    a: &Csr<'_>,
    b: &[f64],
) -> Result<Vec<f64>, NumericError> {
    lu_solve_many(a, b, 1)
}

/// Solves `A X = B` with one sparse LU factorisation for all columns of `B`.
///
/// The input and output are row-major arrays of shape `(A.rows, rhs_cols)`.
/// A zero-row matrix or zero-column RHS returns an empty result after input
/// validation.
///
/// # Errors
/// As [`lu_solve`], including an overflow or mismatch in the RHS shape.
pub fn lu_solve_many(
    a: &Csr<'_>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<Vec<f64>, NumericError> {
    validate_direct_rhs(a, b, rhs_cols)?;
    if b.is_empty() {
        return Ok(Vec::new());
    }
    let matrix = to_csc(a).map_err(|_| NumericError::Singular)?;
    solve_lu(&matrix, b, rhs_cols)
}

fn validate_direct_rhs(
    a: &Csr<'_>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<(), NumericError> {
    validate(a)?;
    if a.rows != a.cols {
        return Err(NumericError::NotSquare);
    }
    if a.rows.checked_mul(rhs_cols) != Some(b.len()) {
        return Err(NumericError::ShapeMismatch(if rhs_cols == 1 {
            format!(
                "matrix has {} rows, right-hand side has {} entries",
                a.rows,
                b.len()
            )
        } else {
            format!(
                "right-hand side has {} entries for {} rows and {rhs_cols} columns",
                b.len(),
                a.rows
            )
        }));
    }
    Ok(())
}

fn to_csc(a: &Csr<'_>) -> Result<SparseColMat<usize, f64>, FaerError> {
    if a.cols > isize::MAX as usize || a.rows > isize::MAX as usize {
        return Err(FaerError::IndexOverflow);
    }
    // Faer transposes the CSR storage in linear time and sums duplicates.
    let symbolic =
        SymbolicSparseRowMatRef::new_unsorted_checked(a.rows, a.cols, a.row_ptr, None, a.col_idx);
    SparseRowMatRef::new(symbolic, a.values).to_col_major()
}

fn solve_lu(
    matrix: &SparseColMat<usize, f64>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<Vec<f64>, NumericError> {
    if matrix.val().iter().any(|value| !value.is_finite()) {
        return Err(NumericError::Singular);
    }
    let rows = matrix.nrows();
    if matrix.val().len() == rows
        && matrix.col_ptr().iter().copied().eq(0..=rows)
        && matrix.row_idx().iter().copied().eq(0..rows)
    {
        // Dividing each RHS by the pivot also handles subnormal pivots
        // whose reciprocals would overflow during a triangular solve.
        let mut solution = copy_rhs(b)?;
        for (row, &pivot) in solution.chunks_exact_mut(rhs_cols).zip(matrix.val()) {
            for value in row {
                *value /= pivot;
            }
        }
        return finite_solution(solution, NumericError::Singular);
    }
    let lu = matrix.sp_lu().map_err(|_| NumericError::Singular)?;
    solve_factored(&lu, rows, b, rhs_cols, NumericError::Singular)
}

fn copy_rhs(b: &[f64]) -> Result<Vec<f64>, NumericError> {
    let mut solution = Vec::new();
    solution
        .try_reserve_exact(b.len())
        .map_err(|_| NumericError::Domain("sparse solution allocation exceeds capacity".into()))?;
    solution.extend_from_slice(b);
    Ok(solution)
}

fn finite_solution(
    solution: Vec<f64>,
    failure: NumericError,
) -> Result<Vec<f64>, NumericError> {
    if solution.iter().any(|value| !value.is_finite()) {
        return Err(failure);
    }
    Ok(solution)
}

fn solve_factored(
    factor: &impl Solve<f64>,
    rows: usize,
    b: &[f64],
    rhs_cols: usize,
    failure: NumericError,
) -> Result<Vec<f64>, NumericError> {
    let mut solution = copy_rhs(b)?;
    factor.solve_in_place(faer::MatMut::from_row_major_slice_mut(
        &mut solution,
        rows,
        rhs_cols,
    ));
    finite_solution(solution, failure)
}

/// The result of [`cg`].
#[derive(Clone, Debug)]
pub struct CgResult {
    /// The final iterate.
    pub x: Vec<f64>,
    /// Iterations performed.
    pub iterations: usize,
    /// The final relative residual `‖b − A x‖ / ‖b‖` (`‖r‖` when `b = 0`).
    pub residual: f64,
    /// Whether the relative residual reached the tolerance.
    pub converged: bool,
}

/// Unpreconditioned conjugate gradients for a symmetric positive-definite
/// matrix, starting from zero.
///
/// # Errors
/// Returns shape and input-domain errors, or [`NumericError::NotPositiveDefinite`]
/// when a search direction has nonpositive curvature or the iteration breaks down.
pub fn cg(
    a: &Csr<'_>,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<CgResult, NumericError> {
    validate_iterative_rhs(a, b, tolerance, "CG")?;
    let result = krylov::cg_positive(
        |x, y| a.matvec_validated(x, y),
        b,
        tolerance,
        max_iterations,
    )
    .map_err(|_| NumericError::NotPositiveDefinite)?;
    Ok(CgResult {
        x: result.x,
        iterations: result.iterations,
        residual: result.residual,
        converged: result.converged,
    })
}

fn validate_iterative_rhs(
    a: &Csr<'_>,
    b: &[f64],
    tolerance: f64,
    method: &str,
) -> Result<(), NumericError> {
    validate(a)?;
    if a.rows != a.cols {
        return Err(NumericError::NotSquare);
    }
    if b.len() != a.rows {
        return Err(NumericError::ShapeMismatch(format!(
            "matrix has {} rows, right-hand side has {} entries",
            a.rows,
            b.len()
        )));
    }
    if !tolerance.is_finite() || tolerance < 0.0 || b.iter().any(|value| !value.is_finite()) {
        return Err(NumericError::Domain(format!(
            "{method} needs a finite right-hand side and finite nonnegative tolerance"
        )));
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod cg_contract_tests {
    use super::Csr;
    use super::NumericError;
    use super::PrecondKind;
    use super::cg;
    use super::cg_precond;

    #[test]
    fn native_cg_checks_the_original_residual_at_extreme_rhs_scales() {
        let identity = Csr {
            rows: 1,
            cols: 1,
            row_ptr: &[0, 1],
            col_idx: &[0],
            values: &[1.0],
        };
        for rhs in [1e-200, 1e200] {
            let initial = cg(&identity, &[rhs], 1e-8, 0).unwrap();
            assert!(!initial.converged);
            assert!((initial.residual - 1.0).abs() <= f64::EPSILON);
            let solved = cg(&identity, &[rhs], 1e-8, 4).unwrap();
            assert!(solved.converged);
            assert_eq!(solved.iterations, 1);
            assert!((solved.x[0] / rhs - 1.0).abs() <= f64::EPSILON);
            assert!(solved.residual <= f64::EPSILON);
        }
    }

    #[test]
    fn native_cg_still_rejects_negative_and_zero_curvature() {
        let indefinite = Csr {
            rows: 2,
            cols: 2,
            row_ptr: &[0, 2, 4],
            col_idx: &[0, 1, 0, 1],
            values: &[1.0, 2.0, 2.0, 1.0],
        };
        assert!(matches!(
            cg(&indefinite, &[1.0, -1.0], 1e-8, 4),
            Err(NumericError::NotPositiveDefinite)
        ));
        let singular = Csr {
            values: &[1.0, 1.0, 1.0, 1.0],
            ..indefinite
        };
        assert!(matches!(
            cg(&singular, &[1.0, -1.0], 1e-8, 4),
            Err(NumericError::NotPositiveDefinite)
        ));
    }

    #[test]
    fn automatic_cg_applies_the_incomplete_cholesky_it_selected() {
        let matrix = Csr {
            rows: 3,
            cols: 3,
            row_ptr: &[0, 2, 5, 7],
            col_idx: &[0, 1, 0, 1, 2, 1, 2],
            values: &[4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0],
        };
        let rhs = [3.0, 2.0, 3.0];
        let auto = cg_precond(&matrix, &rhs, 1e-12, 8, PrecondKind::Auto, None).unwrap();
        let ic = cg_precond(&matrix, &rhs, 1e-12, 8, PrecondKind::Ic0, None).unwrap();
        assert!(auto.converged);
        assert_eq!(auto.precond, "Ic0");
        assert_eq!(auto.iterations, ic.iterations);
        assert_eq!(auto.residual.to_bits(), ic.residual.to_bits());
        assert_eq!(auto.x, ic.x);
    }
}

/// Whether `a` is exactly symmetric (pattern and values) with a positive
/// diagonal: the cheap necessary conditions for a Cholesky attempt.
fn looks_spd(matrix: &SparseColMat<usize, f64>) -> bool {
    (0..matrix.nrows()).all(|i| matrix.get(i, i).is_some_and(|v| v.is_finite() && *v > 0.0))
        && matrix.triplet_iter().all(|entry| {
            entry.row == entry.col
                || (entry.val.is_finite()
                    && matrix
                        .get(entry.col, entry.row)
                        .is_some_and(|v| crate::kernels::exact_eq(*v, *entry.val)))
        })
}

#[cfg(test)]
#[path = "sparse_native/direct_tests.rs"]
mod direct_tests;

/// Solves `A x = b` by sparse Cholesky (`A = L Lᵀ`, lower triangle read).
///
/// # Errors
/// [`NumericError::NotPositiveDefinite`] when `A` is not symmetric positive
/// definite (or the factorisation breaks down), plus the shape errors of
/// [`lu_solve`].
pub fn cholesky_solve(
    a: &Csr<'_>,
    b: &[f64],
) -> Result<Vec<f64>, NumericError> {
    cholesky_solve_many(a, b, 1)
}

/// Solves `A X = B` with one sparse Cholesky factorisation for all columns.
///
/// The RHS and result have row-major shape `(A.rows, rhs_cols)`.
/// The stored matrix is checked for symmetry after adding duplicate entries.
/// For a nonempty matrix, zero columns return an empty result after CSR and
/// RHS validation.
///
/// # Errors
/// As [`cholesky_solve`], including an overflow or mismatch in the RHS shape.
pub fn cholesky_solve_many(
    a: &Csr<'_>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<Vec<f64>, NumericError> {
    validate_direct_rhs(a, b, rhs_cols)?;
    if a.rows == 0 {
        return Err(NumericError::NotPositiveDefinite);
    }
    if b.is_empty() {
        return Ok(Vec::new());
    }
    let matrix = to_csc(a).map_err(|_| NumericError::NotPositiveDefinite)?;
    solve_cholesky(&matrix, b, rhs_cols)
}

fn solve_cholesky(
    matrix: &SparseColMat<usize, f64>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<Vec<f64>, NumericError> {
    if !looks_spd(matrix) {
        return Err(NumericError::NotPositiveDefinite);
    }
    let llt = matrix
        .sp_cholesky(faer::Side::Lower)
        .map_err(|_| NumericError::NotPositiveDefinite)?;
    solve_factored(
        &llt,
        matrix.nrows(),
        b,
        rhs_cols,
        NumericError::NotPositiveDefinite,
    )
}

/// Solves `A x = b`: sparse Cholesky when `A` is symmetric positive
/// definite, otherwise (or when the factorisation fails) sparse LU. The
/// second component names the factorisation used (`"Cholesky"` or `"LU"`).
///
/// # Errors
/// As [`lu_solve`].
pub fn solve_auto(
    a: &Csr<'_>,
    b: &[f64],
) -> Result<(Vec<f64>, &'static str), NumericError> {
    solve_auto_many(a, b, 1)
}

/// Solves a batch of right-hand sides with sparse Cholesky or sparse LU.
///
/// The RHS and result use row-major shape `(A.rows, rhs_cols)`. The returned
/// method name is `"Cholesky"` or `"LU"` as in [`solve_auto`].
/// An empty batch is tagged `"LU"` after input validation.
///
/// # Errors
/// As [`lu_solve_many`].
pub fn solve_auto_many(
    a: &Csr<'_>,
    b: &[f64],
    rhs_cols: usize,
) -> Result<(Vec<f64>, &'static str), NumericError> {
    validate_direct_rhs(a, b, rhs_cols)?;
    if b.is_empty() {
        return Ok((Vec::new(), "LU"));
    }
    let matrix = to_csc(a).map_err(|_| NumericError::Singular)?;
    if let Ok(x) = solve_cholesky(&matrix, b, rhs_cols) {
        return Ok((x, "Cholesky"));
    }
    solve_lu(&matrix, b, rhs_cols).map(|x| (x, "LU"))
}

/// The preconditioner a caller asks [`cg_precond`] for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrecondKind {
    /// No preconditioner.
    Identity,
    /// Diagonal scaling.
    Jacobi,
    /// Inverted diagonal blocks (rssn `Bsr::block_jacobi`).
    BlockJacobi,
    /// Incomplete LU with the sparsity pattern (block ILU(0) for a block
    /// operand).
    Ilu0,
    /// Incomplete Cholesky with the sparsity pattern.
    Ic0,
    /// The strongest factorisation that succeeds.
    Auto,
}

/// The result of [`cg_precond`].
#[derive(Clone, Debug)]
pub struct PrecondCgResult {
    /// The final iterate.
    pub x: Vec<f64>,
    /// Iterations performed.
    pub iterations: usize,
    /// The relative residual `‖b − A x‖ / ‖b‖`.
    pub residual: f64,
    /// Whether the tolerance was reached.
    pub converged: bool,
    /// The preconditioner actually used: `"Identity"`, `"Jacobi"`,
    /// `"BlockJacobi"`, `"BlockIlu0"`, `"Ilu0"` or `"Ic0"`.
    pub precond: &'static str,
}

enum Built {
    Identity,
    Jacobi(Vec<f64>),
    Block(crate::kernels::bsr::BlockJacobi),
    BlockIlu(crate::kernels::bsr::Bilu0),
    Ilu(crate::kernels::krylov::Ilu0),
    Ic(crate::kernels::krylov::Ic0),
    Auto(crate::kernels::krylov::AutoPreconditioner),
}

impl Built {
    fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        match self {
            | Self::Identity => z.copy_from_slice(r),
            | Self::Jacobi(inv) => {
                for ((zi, ri), di) in z.iter_mut().zip(r).zip(inv) {
                    *zi = ri * di;
                }
            },
            | Self::Block(p) => p.apply(r, z),
            | Self::BlockIlu(p) => p.apply(r, z),
            | Self::Ilu(p) => p.apply(r, z),
            | Self::Ic(p) => p.apply(r, z),
            | Self::Auto(p) => p.apply(r, z),
        }
    }

    const fn name(&self) -> &'static str {
        match self {
            | Self::Identity => "Identity",
            | Self::Jacobi(_) => "Jacobi",
            | Self::Block(_) => "BlockJacobi",
            | Self::BlockIlu(_) => "BlockIlu0",
            | Self::Ilu(_) => "Ilu0",
            | Self::Ic(_) => "Ic0",
            | Self::Auto(p) => {
                match p.kind() {
                    | krylov::AutoPreconditionerKind::Ic0 => "Ic0",
                    | krylov::AutoPreconditionerKind::Ilu0 => "Ilu0",
                    | krylov::AutoPreconditionerKind::Jacobi => "Jacobi",
                }
            },
        }
    }
}

/// The square block shape the block preconditioners use: the operand's own
/// block when it is square, else the automatic choice.
fn square_block(
    csr: &crate::kernels::krylov::Csr,
    hint: Option<(usize, usize)>,
) -> (usize, usize) {
    match hint {
        | Some((r, c)) if r > 0 && r == c && csr.rows.is_multiple_of(r) => (r, c),
        | _ => crate::kernels::bsr::choose_block_size(csr, &[(4, 4), (3, 3), (2, 2), (1, 1)]),
    }
}

fn build_preconditioner(
    a: &Csr<'_>,
    kind: PrecondKind,
    block: Option<(usize, usize)>,
) -> Built {
    use crate::kernels::bsr::Bsr;
    let jacobi = || {
        Built::Jacobi(krylov::csr_jacobi_inverse(
            a.rows, a.row_ptr, a.col_idx, a.values,
        ))
    };
    let own_pattern = || {
        krylov::Csr {
            rows: a.rows,
            cols: a.cols,
            indptr: a.row_ptr.to_vec(),
            indices: a.col_idx.to_vec(),
            values: a.values.to_vec(),
        }
    };
    let blocked = |csr: &krylov::Csr| {
        let (br, bc) = square_block(csr, block);
        Bsr::from_csr(csr, br, bc).ok()
    };
    match kind {
        | PrecondKind::Identity => Built::Identity,
        | PrecondKind::Jacobi => jacobi(),
        | PrecondKind::BlockJacobi => {
            blocked(&own_pattern())
                .and_then(|b| b.block_jacobi().ok())
                .map_or_else(jacobi, Built::Block)
        },
        | PrecondKind::Ilu0 => {
            let csr = own_pattern();
            if block.is_some() {
                blocked(&csr)
                    .and_then(|b| b.bilu0().ok())
                    .map(Built::BlockIlu)
                    .or_else(|| csr.ilu0().map(Built::Ilu))
                    .unwrap_or_else(jacobi)
            } else {
                csr.ilu0().map_or_else(jacobi, Built::Ilu)
            }
        },
        | PrecondKind::Ic0 => own_pattern().ic0().map_or_else(jacobi, Built::Ic),
        | PrecondKind::Auto => {
            let csr = own_pattern();
            if block.is_some()
                && let Some(b) = blocked(&csr)
            {
                if let Ok(p) = b.bilu0() {
                    return Built::BlockIlu(p);
                }
                if let Ok(p) = b.block_jacobi() {
                    return Built::Block(p);
                }
            }
            Built::Auto(krylov::auto_preconditioner(&csr))
        },
    }
}

/// Preconditioned conjugate gradients from a zero start.
///
/// `block` is the block shape of a block-CSR operand (`None` otherwise); it
/// selects rssn's block Jacobi / block ILU(0) for those kinds. A
/// factorisation that breaks down falls back to Jacobi, and the kind that was
/// really applied is reported in the result.
///
/// # Errors
/// [`NumericError::NotSquare`], [`NumericError::ShapeMismatch`], or
/// [`NumericError::Domain`] for non-finite inputs or invalid tolerance.
pub fn cg_precond(
    a: &Csr<'_>,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
    kind: PrecondKind,
    block: Option<(usize, usize)>,
) -> Result<PrecondCgResult, NumericError> {
    validate_iterative_rhs(a, b, tolerance, "CG")?;
    let built = build_preconditioner(a, kind, block);
    let result = crate::kernels::krylov::cg(
        |x, y| a.matvec_validated(x, y),
        |r, z| built.apply(r, z),
        b,
        None,
        tolerance,
        max_iterations,
    );
    Ok(PrecondCgResult {
        x: result.x,
        iterations: result.iterations,
        residual: result.residual,
        converged: result.converged,
        precond: built.name(),
    })
}

/// Which Krylov method [`krylov_precond`] runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Krylov {
    /// Restarted GMRES(`restart`) with right preconditioning: any square
    /// nonsingular matrix; the residual norm never increases.
    Gmres {
        /// The restart length (the Krylov basis size before restarting).
        restart: usize,
    },
    /// `BiCGSTAB`: any square nonsingular matrix, short recurrences (constant
    /// memory), with a possibly irregular residual history.
    Bicgstab,
}

/// Preconditioned GMRES or `BiCGSTAB` from a zero start, for a **general**
/// (possibly nonsymmetric) square matrix.
///
/// The preconditioners are those of [`cg_precond`] except that `Ic0` and the
/// symmetric `Auto` choice need a symmetric matrix: here `Ic0` is replaced by
/// `Ilu0` and `Auto` tries `Ilu0` then `Jacobi`. A factorisation that breaks
/// down falls back to Jacobi, and the kind actually applied is reported. The
/// Krylov recurrence is sequential. The sparse products can run in parallel
/// with the same summation order per row, so the result is independent of the
/// thread count.
///
/// # Errors
/// [`NumericError::NotSquare`], [`NumericError::ShapeMismatch`], or
/// [`NumericError::Domain`] for non-finite inputs or invalid tolerance.
#[allow(clippy::too_many_arguments)] // Explicit solver, preconditioner and block contracts.
pub fn krylov_precond(
    a: &Csr<'_>,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
    method: Krylov,
    kind: PrecondKind,
    block: Option<(usize, usize)>,
) -> Result<PrecondCgResult, NumericError> {
    validate_iterative_rhs(a, b, tolerance, "Krylov")?;
    let kind = match kind {
        | PrecondKind::Ic0 => PrecondKind::Ilu0,
        | PrecondKind::Auto if block.is_none() => PrecondKind::Ilu0,
        | other => other,
    };
    let built = build_preconditioner(a, kind, block);
    let apply = |x: &[f64], y: &mut [f64]| a.matvec_validated(x, y);
    let precondition = |r: &[f64], z: &mut [f64]| built.apply(r, z);
    let result = match method {
        | Krylov::Gmres { restart } => {
            crate::kernels::krylov::gmres(
                apply,
                precondition,
                b,
                None,
                restart.max(1),
                tolerance,
                max_iterations,
            )
        },
        | Krylov::Bicgstab => {
            crate::kernels::krylov::bicgstab(
                apply,
                precondition,
                b,
                None,
                tolerance,
                max_iterations,
            )
        },
    };
    Ok(PrecondCgResult {
        x: result.x,
        iterations: result.iterations,
        residual: result.residual,
        converged: result.converged,
        precond: built.name(),
    })
}
