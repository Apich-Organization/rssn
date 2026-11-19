//! `f64` kernels for compressed sparse matrices (§8.4.1).
//!
//! Everything here takes the CSR arrays of a square-or-rectangular matrix
//! (`row_ptr`, `col_idx`, `values`) and a dense right-hand side. Nothing is
//! ever densified: the direct solver hands the entries to `faer`'s sparse LU
//! (fill-reducing COLAMD ordering, left-looking numeric factorisation), and
//! the iterative solver only needs the matrix-vector product.

use crate::kernels::krylov;
use crate::kernels::numeric_common::NumericError;
use faer::linalg::solvers::Solve;
use faer::sparse::SparseColMat;
use faer::sparse::Triplet;

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
/// fails or produces a non-finite solution.
pub fn lu_solve(
    a: &Csr<'_>,
    b: &[f64],
) -> Result<Vec<f64>, NumericError> {
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
    let n = a.rows;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut triplets = Vec::with_capacity(a.values.len());
    for row in 0..n {
        for k in a.row_ptr[row]..a.row_ptr[row + 1] {
            triplets.push(Triplet::new(row, a.col_idx[k], a.values[k]));
        }
    }
    let matrix = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
        .map_err(|_| NumericError::Singular)?;
    let lu = matrix.sp_lu().map_err(|_| NumericError::Singular)?;
    let rhs = faer::Col::<f64>::from_fn(n, |i| b[i]);
    let solution = lu.solve(&rhs);
    let out: Vec<f64> = (0..n).map(|i| solution[i]).collect();
    if out.iter().any(|value| !value.is_finite()) {
        return Err(NumericError::Singular);
    }
    Ok(out)
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
        return Err(NumericError::Domain(
            "CG needs a finite right-hand side and finite nonnegative tolerance".into(),
        ));
    }
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod cg_contract_tests {
    use super::Csr;
    use super::NumericError;
    use super::cg;

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
}

/// Whether `a` is exactly symmetric (pattern and values) with a positive
/// diagonal: the cheap necessary conditions for a Cholesky attempt.
fn looks_spd(a: &Csr<'_>) -> bool {
    if a.rows != a.cols || a.rows == 0 {
        return false;
    }
    let entry = |i: usize, j: usize| -> Option<f64> {
        let row = &a.col_idx[a.row_ptr[i]..a.row_ptr[i + 1]];
        row.binary_search(&j)
            .ok()
            .map(|p| a.values[a.row_ptr[i] + p])
    };
    for i in 0..a.rows {
        match entry(i, i) {
            | Some(d) if d > 0.0 => {},
            | _ => return false,
        }
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_idx[k];
            if j != i
                && entry(j, i).is_none_or(|value| !crate::kernels::exact_eq(value, a.values[k]))
            {
                return false;
            }
        }
    }
    true
}

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
    if !looks_spd(a) {
        return Err(NumericError::NotPositiveDefinite);
    }
    let n = a.rows;
    let mut triplets = Vec::with_capacity(a.values.len());
    for row in 0..n {
        for k in a.row_ptr[row]..a.row_ptr[row + 1] {
            triplets.push(Triplet::new(row, a.col_idx[k], a.values[k]));
        }
    }
    let matrix = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
        .map_err(|_| NumericError::NotPositiveDefinite)?;
    let llt = matrix
        .sp_cholesky(faer::Side::Lower)
        .map_err(|_| NumericError::NotPositiveDefinite)?;
    let rhs = faer::Col::<f64>::from_fn(n, |i| b[i]);
    let solution = llt.solve(&rhs);
    let out: Vec<f64> = (0..n).map(|i| solution[i]).collect();
    if out.iter().any(|value| !value.is_finite()) {
        return Err(NumericError::NotPositiveDefinite);
    }
    Ok(out)
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
    if let Ok(x) = cholesky_solve(a, b) {
        return Ok((x, "Cholesky"));
    }
    lu_solve(a, b).map(|x| (x, "LU"))
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
    csr: &crate::kernels::krylov::Csr,
    kind: PrecondKind,
    block: Option<(usize, usize)>,
) -> Built {
    use crate::kernels::bsr::Bsr;
    let jacobi = || Built::Jacobi(csr.jacobi_inverse());
    let blocked = || {
        let (br, bc) = square_block(csr, block);
        Bsr::from_csr(csr, br, bc).ok()
    };
    match kind {
        | PrecondKind::Identity => Built::Identity,
        | PrecondKind::Jacobi => jacobi(),
        | PrecondKind::BlockJacobi => {
            blocked()
                .and_then(|b| b.block_jacobi().ok())
                .map_or_else(jacobi, Built::Block)
        },
        | PrecondKind::Ilu0 => {
            if block.is_some() {
                blocked()
                    .and_then(|b| b.bilu0().ok())
                    .map(Built::BlockIlu)
                    .or_else(|| csr.ilu0().map(Built::Ilu))
                    .unwrap_or_else(jacobi)
            } else {
                csr.ilu0().map_or_else(jacobi, Built::Ilu)
            }
        },
        | PrecondKind::Ic0 => csr.ic0().map_or_else(jacobi, Built::Ic),
        | PrecondKind::Auto => {
            if block.is_some()
                && let Some(b) = blocked()
            {
                if let Ok(p) = b.bilu0() {
                    return Built::BlockIlu(p);
                }
                if let Ok(p) = b.block_jacobi() {
                    return Built::Block(p);
                }
            }
            let auto = crate::kernels::krylov::auto_preconditioner(csr);
            match auto.kind() {
                | crate::kernels::krylov::AutoPreconditionerKind::Ic0 => {
                    csr.ic0().map_or_else(jacobi, Built::Ic)
                },
                | crate::kernels::krylov::AutoPreconditionerKind::Ilu0 => {
                    csr.ilu0().map_or_else(jacobi, Built::Ilu)
                },
                | crate::kernels::krylov::AutoPreconditionerKind::Jacobi => jacobi(),
            }
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
/// [`NumericError::NotSquare`] and [`NumericError::ShapeMismatch`].
pub fn cg_precond(
    a: &Csr<'_>,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
    kind: PrecondKind,
    block: Option<(usize, usize)>,
) -> Result<PrecondCgResult, NumericError> {
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
    let csr = crate::kernels::krylov::Csr {
        rows: a.rows,
        cols: a.cols,
        indptr: a.row_ptr.to_vec(),
        indices: a.col_idx.to_vec(),
        values: a.values.to_vec(),
    };
    let built = build_preconditioner(&csr, kind, block);
    let result = crate::kernels::krylov::cg(
        |x, y| csr.matvec(x, y),
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
/// iteration is sequential (rssn's Krylov kernels), so the result does not
/// depend on the thread count.
///
/// # Errors
/// [`NumericError::NotSquare`] and [`NumericError::ShapeMismatch`].
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
    let csr = crate::kernels::krylov::Csr {
        rows: a.rows,
        cols: a.cols,
        indptr: a.row_ptr.to_vec(),
        indices: a.col_idx.to_vec(),
        values: a.values.to_vec(),
    };
    let kind = match kind {
        | PrecondKind::Ic0 => PrecondKind::Ilu0,
        | PrecondKind::Auto if block.is_none() => PrecondKind::Ilu0,
        | other => other,
    };
    let built = build_preconditioner(&csr, kind, block);
    let apply = |x: &[f64], y: &mut [f64]| csr.matvec(x, y);
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
