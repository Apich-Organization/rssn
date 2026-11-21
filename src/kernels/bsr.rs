//! Block sparse row (BSR) matrices: block matvec, transpose, product, block Jacobi and block ILU(0).
//!
//! # Block sparse row (BSR) matrices
//!
//! A [`Bsr`] stores a matrix as a grid of dense `br x bc` blocks, in
//! compressed-sparse-row layout *over the blocks*. It pays off for systems
//! with a nodal structure (several unknowns per mesh node: elasticity,
//! Navier-Stokes, coupled reaction-diffusion), where the scalar CSR
//! `indices` array repeats the same pattern `br * bc` times and the matvec
//! is a chain of dependent scalar gathers instead of small dense products.
//!
//! # Representation
//!
//! * `nrows`, `ncols` are the *scalar* extents and must be divisible by the
//!   block shape `(br, bc)`; there are `nrows / br` block rows and
//!   `ncols / bc` block columns.
//! * `indptr` has `nrows / br + 1` entries and counts stored blocks per block
//!   row; `indices` holds the block-column index of every stored block,
//!   strictly increasing within a block row.
//! * `values` holds `br * bc` scalars per stored block, row-major inside the
//!   block, blocks contiguous in the order of `indices`.
//!
//! A block is stored iff at least one of its entries is nonzero
//! ([`Bsr::from_csr`]); entries of a stored block that are zero are padding.
//! [`Bsr::fill_ratio`] is the stored-scalars / true-nonzeros ratio and
//! [`choose_block_size`] uses it to pick a block shape.
//!
//! # Operations
//!
//! * [`Bsr::matvec`] is bit-identical to [`Csr::matvec`] on the same matrix
//!   (every output entry is the correctly rounded dot product of its row,
//!   and padding contributes exact zeros); block rows are processed in
//!   parallel from [`BSR_PAR_NNZ`] stored scalars on.
//! * [`Bsr::transpose`], sparse [`Bsr::matmul`], dense batch
//!   [`Bsr::matmul_dense`], and [`Bsr::diagonal_blocks`].
//! * [`Bsr::block_jacobi`] ([`BlockJacobi`]) and [`Bsr::bilu0`] ([`Bilu0`])
//!   are block preconditioners for the matrix-free Krylov solvers in
//!   [`crate::kernels::krylov`]: wrap `apply` in a closure.
//!
//! # Complexity
//!
//! With `B = br * bc`, `nb` stored blocks: storage `nb * B + nb + nbr` words;
//! matvec `O(nb * B)`; transpose `O(nb * B + nbc)`; block Jacobi set-up
//! `O(nbr * b^3)`, apply `O(nbr * b^2)`; BILU(0) set-up
//! `O(sum over rows of (blocks in row)^2 * b^3)`, apply `O(nb * b^2)`.

use super::krylov::Csr;
use super::krylov::dot_exact;
use std::collections::HashMap;

#[path = "bsr/product.rs"]
mod product;

/// Stored scalars from which [`Bsr::matvec`] runs in parallel.
pub const BSR_PAR_NNZ: usize = 60_000;
/// Block rows per parallel chunk of [`Bsr::matvec`].
const BSR_CHUNK_ROWS: usize = 64;
/// The largest fill ratio [`choose_block_size`] accepts.
pub const MAX_FILL_RATIO: f64 = 1.5;

/// Errors of the block sparse routines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BsrError {
    /// A block extent is zero.
    ZeroBlock,
    /// A scalar extent is not a multiple of the block extent.
    NotDivisible {
        /// The scalar extent.
        extent: usize,
        /// The block extent.
        block: usize,
    },
    /// The operand shapes (or block shapes) do not fit the operation.
    ShapeMismatch(String),
    /// The operation needs a square matrix with square blocks.
    NotSquare,
    /// A diagonal block is missing from the pattern (BILU(0)).
    MissingDiagonal {
        /// The block row.
        block_row: usize,
    },
    /// A diagonal or pivot block is numerically singular.
    SingularBlock {
        /// The block row.
        block_row: usize,
    },
    /// A storage shape, allocation, or representation invariant is invalid.
    Invalid(String),
}

impl std::fmt::Display for BsrError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::ZeroBlock => write!(f, "block extents must be positive"),
            | Self::NotDivisible { extent, block } => {
                write!(
                    f,
                    "extent {extent} is not a multiple of the block extent {block}"
                )
            },
            | Self::ShapeMismatch(what) => write!(f, "shape mismatch: {what}"),
            | Self::NotSquare => write!(f, "a square matrix with square blocks is required"),
            | Self::MissingDiagonal { block_row } => {
                write!(f, "diagonal block of block row {block_row} is not stored")
            },
            | Self::SingularBlock { block_row } => {
                write!(f, "diagonal block of block row {block_row} is singular")
            },
            | Self::Invalid(what) => write!(f, "invalid block sparse matrix: {what}"),
        }
    }
}

impl std::error::Error for BsrError {}

/// Block sparse row matrix; see the [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct Bsr {
    /// Number of scalar rows (a multiple of `br`).
    pub nrows: usize,
    /// Number of scalar columns (a multiple of `bc`).
    pub ncols: usize,
    /// Block height.
    pub br: usize,
    /// Block width.
    pub bc: usize,
    /// Block-row pointers (`nrows / br + 1` entries).
    pub indptr: Vec<usize>,
    /// Block-column indices, strictly increasing within each block row.
    pub indices: Vec<usize>,
    /// Dense row-major blocks, `br * bc` scalars each.
    pub values: Vec<f64>,
}

fn check_block(
    nrows: usize,
    ncols: usize,
    br: usize,
    bc: usize,
) -> Result<(), BsrError> {
    if br == 0 || bc == 0 {
        return Err(BsrError::ZeroBlock);
    }
    if !nrows.is_multiple_of(br) {
        return Err(BsrError::NotDivisible {
            extent: nrows,
            block: br,
        });
    }
    if !ncols.is_multiple_of(bc) {
        return Err(BsrError::NotDivisible {
            extent: ncols,
            block: bc,
        });
    }
    br.checked_mul(bc).ok_or_else(size_error)?;
    Ok(())
}

fn size_error() -> BsrError {
    BsrError::Invalid("sparse storage exceeds representable or available allocation size".into())
}

fn filled<T: Clone>(
    len: usize,
    value: T,
) -> Result<Vec<T>, BsrError> {
    let mut data = Vec::new();
    data.try_reserve_exact(len).map_err(|_| size_error())?;
    data.resize(len, value);
    Ok(data)
}

fn extend_zeros(
    values: &mut Vec<f64>,
    len: usize,
) -> Result<(), BsrError> {
    values
        .try_reserve(len.saturating_sub(values.len()))
        .map_err(|_| size_error())?;
    values.resize(len, 0.0);
    Ok(())
}

fn validate_csr(csr: &Csr) -> Result<(), BsrError> {
    if csr.rows.checked_add(1) != Some(csr.indptr.len())
        || csr.indices.len() != csr.values.len()
        || csr.indptr.first() != Some(&0)
        || csr.indptr.last() != Some(&csr.indices.len())
        || csr.indptr.windows(2).any(|w| w[0] > w[1])
        || csr.indices.iter().any(|&col| col >= csr.cols)
    {
        return Err(BsrError::Invalid("malformed CSR arrays".into()));
    }
    Ok(())
}

/// A dense address table for modest column counts, with sparse storage for
/// shapes whose empty columns would otherwise dominate memory and runtime.
enum BlockSlots {
    Dense(Vec<usize>),
    Sparse(HashMap<usize, usize>),
}

impl BlockSlots {
    fn new(
        cols: usize,
        possible_entries: usize,
    ) -> Result<Self, BsrError> {
        let dense_limit = possible_entries.saturating_mul(4).clamp(1024, 262_144);
        if cols <= dense_limit {
            Ok(Self::Dense(filled(cols, usize::MAX)?))
        } else {
            Ok(Self::Sparse(HashMap::new()))
        }
    }

    fn get(
        &self,
        column: usize,
    ) -> Option<usize> {
        match self {
            | Self::Dense(slots) => (slots[column] != usize::MAX).then_some(slots[column]),
            | Self::Sparse(slots) => slots.get(&column).copied(),
        }
    }

    fn insert(
        &mut self,
        column: usize,
        index: usize,
    ) -> Result<(), BsrError> {
        match self {
            | Self::Dense(slots) => slots[column] = index,
            | Self::Sparse(slots) => {
                slots.try_reserve(1).map_err(|_| size_error())?;
                slots.insert(column, index);
            },
        }
        Ok(())
    }

    fn clear(
        &mut self,
        columns: impl IntoIterator<Item = usize>,
    ) {
        match self {
            | Self::Dense(slots) => {
                for column in columns {
                    slots[column] = usize::MAX;
                }
            },
            | Self::Sparse(slots) => slots.clear(),
        }
    }
}

impl Bsr {
    /// An all-zero matrix with no stored blocks.
    ///
    /// # Errors
    /// A zero or incompatible block extent, or an unrepresentable allocation.
    pub fn zeros(
        nrows: usize,
        ncols: usize,
        br: usize,
        bc: usize,
    ) -> Result<Self, BsrError> {
        check_block(nrows, ncols, br, bc)?;
        let ptr_len = (nrows / br).checked_add(1).ok_or_else(size_error)?;
        Ok(Self {
            nrows,
            ncols,
            br,
            bc,
            indptr: filled(ptr_len, 0)?,
            indices: Vec::new(),
            values: Vec::new(),
        })
    }

    /// Number of block rows.
    #[must_use]
    pub const fn block_rows(&self) -> usize {
        self.nrows / self.br
    }

    /// Number of block columns.
    #[must_use]
    pub const fn block_cols(&self) -> usize {
        self.ncols / self.bc
    }

    /// Number of stored blocks.
    #[must_use]
    pub const fn block_count(&self) -> usize {
        self.indices.len()
    }

    /// Number of stored scalars (`block_count * br * bc`, padding included).
    #[must_use]
    pub const fn stored_len(&self) -> usize {
        self.values.len()
    }

    /// Number of nonzero scalars.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.values.iter().filter(|v| **v != 0.0).count()
    }

    /// Stored scalars divided by true nonzeros (`1.0` for an empty matrix).
    ///
    /// `1.0` means the blocks are completely dense; the larger the value the
    /// more padding zeros are stored.
    #[must_use]
    pub fn fill_ratio(&self) -> f64 {
        let nnz = self.nnz();
        if nnz == 0 {
            1.0
        } else {
            self.values.len() as f64 / nnz as f64
        }
    }

    /// Checks the representation invariants.
    ///
    /// # Errors
    /// A description of the first violated invariant.
    pub fn validate(&self) -> Result<(), String> {
        check_block(self.nrows, self.ncols, self.br, self.bc).map_err(|e| e.to_string())?;
        let nbr = self.block_rows();
        let nbc = self.block_cols();
        let nb = self.indices.len();
        if nbr.checked_add(1) != Some(self.indptr.len())
            || self.indptr.first() != Some(&0)
            || self.indptr.last() != Some(&nb)
            || self.indptr.windows(2).any(|w| w[0] > w[1])
        {
            return Err("block-row pointers do not match the indices".into());
        }
        if nb.checked_mul(self.br * self.bc) != Some(self.values.len()) {
            return Err("values do not hold br * bc scalars per block".into());
        }
        for row in self.indptr.windows(2) {
            let columns = &self.indices[row[0]..row[1]];
            if columns.iter().any(|&col| col >= nbc) {
                return Err("block-column index out of range".into());
            }
            if columns.windows(2).any(|w| w[0] >= w[1]) {
                return Err("block-column indices are not strictly increasing".into());
            }
        }
        Ok(())
    }

    /// Converts a CSR matrix with one entry per coordinate. A block is stored
    /// when it contains a nonzero, and zeros inside the block are padding.
    ///
    /// # Errors
    /// An incompatible block shape, malformed CSR arrays, or an
    /// unrepresentable allocation.
    pub fn from_csr(
        csr: &Csr,
        br: usize,
        bc: usize,
    ) -> Result<Self, BsrError> {
        check_block(csr.rows, csr.cols, br, bc)?;
        validate_csr(csr)?;
        let mut out = Self::zeros(csr.rows, csr.cols, br, bc)?;
        if csr.indices.is_empty() {
            return Ok(out);
        }
        let bs = br * bc;
        let mut slot = BlockSlots::new(out.block_cols(), csr.indices.len())?;
        for ib in 0..out.block_rows() {
            let first = out.indices.len();
            for r in 0..br {
                let row = ib * br + r;
                for k in csr.indptr[row]..csr.indptr[row + 1] {
                    if csr.values[k] != 0.0 {
                        let jb = csr.indices[k] / bc;
                        if slot.get(jb).is_none() {
                            slot.insert(jb, 0)?;
                            out.indices.try_reserve(1).map_err(|_| size_error())?;
                            out.indices.push(jb);
                        }
                    }
                }
            }
            let columns = &mut out.indices[first..];
            columns.sort_unstable();
            extend_zeros(
                &mut out.values,
                out.indices.len().checked_mul(bs).ok_or_else(size_error)?,
            )?;
            for (offset, &jb) in out.indices[first..].iter().enumerate() {
                slot.insert(jb, first + offset)?;
            }
            for r in 0..br {
                let row = ib * br + r;
                for k in csr.indptr[row]..csr.indptr[row + 1] {
                    if csr.values[k] != 0.0 {
                        let col = csr.indices[k];
                        let at = slot.get(col / bc).ok_or_else(|| {
                            BsrError::Invalid("converted CSR block is missing".into())
                        })?;
                        out.values[at * bs + r * bc + col % bc] = csr.values[k];
                    }
                }
            }
            slot.clear(out.indices[first..].iter().copied());
            out.indptr[ib + 1] = out.indices.len();
        }
        Ok(out)
    }

    /// Builds a matrix from `(row, col, value)` triplets. Out-of-range
    /// triplets are ignored, duplicate coordinates are summed in input order,
    /// and canceling blocks are omitted.
    ///
    /// # Errors
    /// An incompatible block shape or an unrepresentable allocation.
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        br: usize,
        bc: usize,
        trip: &[(usize, usize, f64)],
    ) -> Result<Self, BsrError> {
        let mut out = Self::zeros(nrows, ncols, br, bc)?;
        let bs = br * bc;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(trip.len())
            .map_err(|_| size_error())?;
        entries.extend(
            trip.iter()
                .filter(|&&(row, col, _)| row < nrows && col < ncols)
                .map(|&(row, col, value)| (row / br, col / bc, (row % br) * bc + col % bc, value)),
        );
        entries.sort_by_key(|&(row, col, entry, _)| (row, col, entry));
        for block in entries.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
            let (ib, jb, _, _) = block[0];
            let offset = out.values.len();
            for duplicates in block.chunk_by(|a, b| a.2 == b.2) {
                let (_, _, at, first) = duplicates[0];
                let value = duplicates[1..]
                    .iter()
                    .fold(first, |total, &(_, _, _, value)| total + value);
                if value == 0.0 {
                    continue;
                }
                if out.values.len() == offset {
                    extend_zeros(
                        &mut out.values,
                        offset.checked_add(bs).ok_or_else(size_error)?,
                    )?;
                    out.indices.try_reserve(1).map_err(|_| size_error())?;
                    out.indices.push(jb);
                    out.indptr[ib + 1] += 1;
                }
                out.values[offset + at] = value;
            }
        }
        for ib in 0..out.block_rows() {
            out.indptr[ib + 1] += out.indptr[ib];
        }
        Ok(out)
    }

    /// Converts to CSR, dropping the zeros inside blocks.
    ///
    /// # Panics
    /// A malformed matrix or a scalar row count too large for the CSR
    /// row-pointer array can panic.
    #[must_use]
    pub fn to_csr(&self) -> Csr {
        let nbr = self.block_rows();
        let bs = self.br * self.bc;
        let mut indptr = vec![0usize; self.nrows + 1];
        let nnz = self.nnz();
        let mut indices = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        for ib in 0..nbr {
            for r in 0..self.br {
                for k in self.indptr[ib]..self.indptr[ib + 1] {
                    for c in 0..self.bc {
                        let v = self.values[k * bs + r * self.bc + c];
                        if v != 0.0 {
                            indices.push(self.indices[k] * self.bc + c);
                            values.push(v);
                        }
                    }
                }
                indptr[ib * self.br + r + 1] = indices.len();
            }
        }
        Csr {
            rows: self.nrows,
            cols: self.ncols,
            indptr,
            indices,
            values,
        }
    }

    /// Scalar entries of a row in ascending column order, skipping padding.
    fn row_entries(
        &self,
        ib: usize,
        r: usize,
    ) -> impl Iterator<Item = (usize, f64)> + Clone + '_ {
        let (bs, bc) = (self.br * self.bc, self.bc);
        let row_start = r * bc;
        (self.indptr[ib]..self.indptr[ib + 1]).flat_map(move |k| {
            let offset = k * bs + row_start;
            let column = self.indices[k] * bc;
            self.values[offset..offset + bc]
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, value)| value != 0.0)
                .map(move |(c, value)| (column + c, value))
        })
    }

    fn write_matvec_rows(
        &self,
        first: usize,
        x: &[f64],
        y: &mut [f64],
    ) {
        for (ib, block) in y.chunks_mut(self.br).enumerate() {
            for (r, value) in block.iter_mut().enumerate() {
                *value = dot_exact(self.row_entries(first + ib, r).map(|(c, a)| (a, x[c])));
            }
        }
    }

    /// Matrix-vector product `y = A x`.
    ///
    /// Every entry of `y` is the correctly rounded dot product of its row, so
    /// the result equals [`Csr::matvec`] of [`Bsr::to_csr`] bit for bit,
    /// whatever the block shape; block rows run in parallel for large
    /// matrices without changing a bit. Only the first `min(y.len(), nrows)`
    /// output entries are written.
    ///
    /// # Panics
    /// A malformed matrix or an input missing a referenced entry can panic.
    pub fn matvec(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) {
        let rows = self.nrows.min(y.len());
        if self.values.len() >= BSR_PAR_NNZ {
            let chunk = BSR_CHUNK_ROWS.saturating_mul(self.br);
            super::par::for_chunks_mut(&mut y[..rows], chunk, 2, |b, block| {
                self.write_matvec_rows(b * BSR_CHUNK_ROWS, x, block);
            });
        } else {
            self.write_matvec_rows(0, x, &mut y[..rows]);
        }
    }

    /// Multiplies a row-major dense matrix of shape `(ncols, rhs_cols)` by
    /// this BSR matrix on the left. Returns row-major `(nrows, rhs_cols)`.
    /// Each output is rounded once from the exact dot product, so running
    /// [`Self::matvec`] on each RHS column gives the same bits.
    ///
    /// # Errors
    /// A malformed BSR matrix, a dense operand that disagrees with its given
    /// width, or an unrepresentable or unavailable output allocation.
    pub fn matmul_dense(
        &self,
        x: &[f64],
        rhs_cols: usize,
    ) -> Result<Vec<f64>, BsrError> {
        self.validate().map_err(BsrError::Invalid)?;
        if self.ncols.checked_mul(rhs_cols) != Some(x.len()) {
            return Err(BsrError::ShapeMismatch(format!(
                "dense right-hand side has {} entries for {} rows and {rhs_cols} columns",
                x.len(),
                self.ncols
            )));
        }
        let len = self.nrows.checked_mul(rhs_cols).ok_or_else(size_error)?;
        let mut output = filled(len, 0.0)?;
        if self.indices.is_empty() || len == 0 {
            return Ok(output);
        }
        let write = |first, output: &mut [f64]| {
            for (offset, cell) in output.iter_mut().enumerate() {
                let index = first + offset;
                let row = index / rhs_cols;
                let rhs = index % rhs_cols;
                let entries = self.row_entries(row / self.br, row % self.br);
                *cell = dot_exact(entries.map(|(c, a)| (a, x[c * rhs_cols + rhs])));
            }
        };
        if self.stored_len().saturating_mul(rhs_cols) >= BSR_PAR_NNZ {
            let chunk = 512;
            super::par::for_chunks_mut(&mut output, chunk, 2, |b, out| write(b * chunk, out));
        } else {
            write(0, &mut output);
        }
        Ok(output)
    }

    /// Transposed product `y = A^T x` (`y` has `ncols` entries).
    ///
    /// Accumulates block by block in plain `f64`, which can round differently
    /// from the exact dot products in `transpose().matvec`.
    ///
    /// # Panics
    /// A malformed matrix or slices shorter than the corresponding extents
    /// can panic.
    pub fn matvec_transpose(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) {
        let bs = self.br * self.bc;
        y[..self.ncols].fill(0.0);
        for ib in 0..self.block_rows() {
            for k in self.indptr[ib]..self.indptr[ib + 1] {
                let jb = self.indices[k];
                for r in 0..self.br {
                    let xr = x[ib * self.br + r];
                    for c in 0..self.bc {
                        let a = self.values[k * bs + r * self.bc + c];
                        if a != 0.0 {
                            y[jb * self.bc + c] += a * xr;
                        }
                    }
                }
            }
        }
    }

    /// The transpose, with block shape `(bc, br)`.
    ///
    /// # Panics
    /// A malformed matrix or a block-column count too large for the
    /// transposed row-pointer array can panic.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let (nbr, nbc) = (self.block_rows(), self.block_cols());
        let bs = self.br * self.bc;
        let mut indptr = vec![0usize; nbc + 1];
        for jb in &self.indices {
            indptr[jb + 1] += 1;
        }
        for j in 0..nbc {
            indptr[j + 1] += indptr[j];
        }
        let mut next = indptr.clone();
        let mut indices = vec![0usize; self.indices.len()];
        let mut values = vec![0.0; self.values.len()];
        for ib in 0..nbr {
            for k in self.indptr[ib]..self.indptr[ib + 1] {
                let jb = self.indices[k];
                let at = next[jb];
                next[jb] += 1;
                indices[at] = ib;
                for r in 0..self.br {
                    for c in 0..self.bc {
                        values[at * bs + c * self.br + r] = self.values[k * bs + r * self.bc + c];
                    }
                }
            }
        }
        Self {
            nrows: self.ncols,
            ncols: self.nrows,
            br: self.bc,
            bc: self.br,
            indptr,
            indices,
            values,
        }
    }

    /// The product `A B`; block shape `(self.br, other.bc)`.
    ///
    /// Blocks that cancel to exactly zero are omitted. Large products share
    /// independent ranges of block rows across workers while keeping the
    /// accumulation order within each output block fixed.
    ///
    /// # Errors
    /// Different inner extents or block sizes, malformed sparse arrays, or
    /// an unrepresentable or unavailable sparse allocation.
    pub fn matmul(
        &self,
        other: &Self,
    ) -> Result<Self, BsrError> {
        if self.ncols != other.nrows || self.bc != other.br {
            return Err(BsrError::ShapeMismatch(format!(
                "{}x{} (block {}x{}) times {}x{} (block {}x{})",
                self.nrows,
                self.ncols,
                self.br,
                self.bc,
                other.nrows,
                other.ncols,
                other.br,
                other.bc
            )));
        }
        self.validate().map_err(BsrError::Invalid)?;
        other.validate().map_err(BsrError::Invalid)?;
        product::matmul(self, other)
    }

    /// The square diagonal blocks, one dense row-major `b x b` block per block
    /// row (all zero when the block is not stored).
    ///
    /// # Errors
    /// Non-square or malformed input, or an unrepresentable allocation.
    pub fn diagonal_blocks(&self) -> Result<Vec<Vec<f64>>, BsrError> {
        let bs = self.square_block_area()?;
        let mut diagonal = Vec::new();
        diagonal
            .try_reserve_exact(self.block_rows())
            .map_err(|_| size_error())?;
        for block_row in 0..self.block_rows() {
            let mut block = filled(bs, 0.0)?;
            if let Some(k) = self.diagonal_index(block_row) {
                block.copy_from_slice(&self.values[k * bs..(k + 1) * bs]);
            }
            diagonal.push(block);
        }
        Ok(diagonal)
    }

    fn square_block_area(&self) -> Result<usize, BsrError> {
        if self.nrows != self.ncols || self.br != self.bc {
            return Err(BsrError::NotSquare);
        }
        self.validate().map_err(BsrError::Invalid)?;
        Ok(self.br * self.bc)
    }

    fn diagonal_index(
        &self,
        block_row: usize,
    ) -> Option<usize> {
        let start = self.indptr[block_row];
        let end = self.indptr[block_row + 1];
        let offset = self.indices[start..end].binary_search(&block_row).ok()?;
        Some(start + offset)
    }

    /// Block-Jacobi preconditioner: the inverse of every diagonal block.
    ///
    /// # Errors
    /// Non-square or malformed input, an unrepresentable allocation, or
    /// [`BsrError::SingularBlock`] for the first absent or singular diagonal.
    pub fn block_jacobi(&self) -> Result<BlockJacobi, BsrError> {
        let b = self.br;
        let bs = self.square_block_area()?;
        let rows = self.block_rows();
        let mut inverse = filled(rows.checked_mul(bs).ok_or_else(size_error)?, 0.0)?;
        if rows != 0 {
            let mut workspace = BlockInverse::new(b, bs)?;
            for (block_row, output) in inverse.chunks_exact_mut(bs).enumerate() {
                let missing = BsrError::SingularBlock { block_row };
                let k = self
                    .diagonal_index(block_row)
                    .ok_or(BsrError::SingularBlock { block_row })?;
                if !workspace.invert(&self.values[k * bs..(k + 1) * bs], output) {
                    return Err(missing);
                }
            }
        }
        Ok(BlockJacobi { b, inverse })
    }

    /// Block ILU(0): incomplete block LU keeping the block sparsity pattern.
    ///
    /// # Errors
    /// Non-square or malformed input, an unrepresentable allocation,
    /// [`BsrError::MissingDiagonal`] for an absent block, or
    /// [`BsrError::SingularBlock`] for a singular pivot.
    pub fn bilu0(&self) -> Result<Bilu0, BsrError> {
        let b = self.br;
        let bs = self.square_block_area()?;
        let nbr = self.block_rows();
        let mut lu = self.clone();
        let mut diag = filled(nbr, 0usize)?;
        let mut inv_diag = filled(nbr.checked_mul(bs).ok_or_else(size_error)?, 0.0)?;
        if nbr == 0 {
            return Ok(Bilu0 { lu, diag, inv_diag });
        }
        // Position of block column `j` inside the current block row.
        let mut pos = filled(nbr, usize::MAX)?;
        let mut tmp = filled(bs, 0.0)?;
        let mut workspace = BlockInverse::new(b, bs)?;
        for i in 0..nbr {
            let (start, end) = (lu.indptr[i], lu.indptr[i + 1]);
            for k in start..end {
                pos[lu.indices[k]] = k;
            }
            let d = pos[i];
            if d == usize::MAX {
                return Err(BsrError::MissingDiagonal { block_row: i });
            }
            diag[i] = d;
            // The U row used in an update was factored earlier. Splitting
            // here borrows its factors without copying a block per update.
            let (previous, remaining) = lu.values.split_at_mut(start * bs);
            for ik in start..d {
                let k = lu.indices[ik];
                // L_ik = A_ik * inv(U_kk)
                let offset = (ik - start) * bs;
                tmp.fill(0.0);
                gemm_acc(
                    &mut tmp,
                    &remaining[offset..offset + bs],
                    &inv_diag[k * bs..(k + 1) * bs],
                    b,
                    b,
                    1.0,
                );
                remaining[offset..offset + bs].copy_from_slice(&tmp);
                // A_ij -= L_ik * U_kj for j > k inside the pattern of row i.
                for kj in (diag[k] + 1)..lu.indptr[k + 1] {
                    let j = lu.indices[kj];
                    if pos[j] != usize::MAX {
                        let u = &previous[kj * bs..(kj + 1) * bs];
                        let offset = (pos[j] - start) * bs;
                        gemm_acc(&mut remaining[offset..offset + bs], &tmp, u, b, b, -1.0);
                    }
                }
            }
            let offset = (d - start) * bs;
            let pivot = &remaining[offset..offset + bs];
            if !workspace.invert(pivot, &mut inv_diag[i * bs..(i + 1) * bs]) {
                return Err(BsrError::SingularBlock { block_row: i });
            }
            for k in start..end {
                pos[lu.indices[k]] = usize::MAX;
            }
        }
        Ok(Bilu0 { lu, diag, inv_diag })
    }
}

/// `c += alpha * a * b` for row-major `a` (`m x k`), `b` (`k x n`), `c` (`m x n`).
fn gemm_acc(
    c: &mut [f64],
    a: &[f64],
    b: &[f64],
    k: usize,
    n: usize,
    alpha: f64,
) {
    for (c_row, a_row) in c.chunks_exact_mut(n).zip(a.chunks_exact(k)) {
        for (&aip, b_row) in a_row.iter().zip(b.chunks_exact(n)) {
            let aip = alpha * aip;
            if aip != 0.0 {
                for (value, &bpj) in c_row.iter_mut().zip(b_row) {
                    *value += aip * bpj;
                }
            }
        }
    }
}

/// In-place LU with partial pivoting. Returns the scale applied before
/// factorization or `None` if it breaks down.
#[allow(
    clippy::single_call_fn,
    reason = "factorization and inverse solves have distinct invariants"
)]
fn lu_factor(
    lu: &mut [f64],
    perm: &mut [usize],
    n: usize,
) -> Option<f64> {
    let scale = lu.iter().try_fold(0.0_f64, |largest, value| {
        value.is_finite().then_some(largest.max(value.abs()))
    })?;
    if scale == 0.0 {
        return None;
    }
    let (pivot_scale, normalization) = if (1e-150..=1e150).contains(&scale) {
        (scale, 1.0)
    } else {
        for value in &mut *lu {
            *value /= scale;
        }
        (1.0, scale)
    };
    let min_pivot = f64::EPSILON * pivot_scale;
    for (i, at) in perm.iter_mut().enumerate() {
        *at = i;
    }
    for col in 0..n {
        let mut best = col;
        for r in (col + 1)..n {
            if lu[r * n + col].abs() > lu[best * n + col].abs() {
                best = r;
            }
        }
        let pivot = lu[best * n + col];
        if !pivot.is_finite() || pivot.abs() <= min_pivot {
            return None;
        }
        if best != col {
            for c in 0..n {
                lu.swap(col * n + c, best * n + c);
            }
            perm.swap(col, best);
        }
        for r in (col + 1)..n {
            let factor = lu[r * n + col] / lu[col * n + col];
            lu[r * n + col] = factor;
            for c in (col + 1)..n {
                lu[r * n + c] -= factor * lu[col * n + c];
            }
        }
    }
    Some(normalization)
}

/// Work buffers shared by all block inversions in one factorization.
struct BlockInverse {
    n: usize,
    lu: Vec<f64>,
    perm: Vec<usize>,
    col: Vec<f64>,
}

impl BlockInverse {
    fn new(
        n: usize,
        area: usize,
    ) -> Result<Self, BsrError> {
        Ok(Self {
            n,
            lu: filled(area, 0.0)?,
            perm: filled(n, 0)?,
            col: filled(n, 0.0)?,
        })
    }

    fn invert(
        &mut self,
        block: &[f64],
        output: &mut [f64],
    ) -> bool {
        let n = self.n;
        self.lu.copy_from_slice(block);
        let Some(normalization) = lu_factor(&mut self.lu, &mut self.perm, n) else {
            return false;
        };
        let (lu, perm, col) = (&self.lu, &self.perm, &mut self.col);
        for e in 0..n {
            for r in 0..n {
                col[r] = if perm[r] == e { 1.0 } else { 0.0 };
            }
            for r in 0..n {
                for c in 0..r {
                    col[r] -= lu[r * n + c] * col[c];
                }
            }
            for r in (0..n).rev() {
                for c in (r + 1)..n {
                    col[r] -= lu[r * n + c] * col[c];
                }
                col[r] /= lu[r * n + r];
            }
            for r in 0..n {
                output[r * n + e] = col[r];
            }
        }
        if normalization.to_bits() != 1.0_f64.to_bits() {
            for value in &mut *output {
                *value /= normalization;
            }
        }
        output.iter().all(|v| v.is_finite())
    }
}

/// `y = M x` for a row-major `n x n` block.
fn block_apply(
    m: &[f64],
    n: usize,
    x: &[f64],
    y: &mut [f64],
) {
    for r in 0..n {
        y[r] = (0..n).map(|c| m[r * n + c] * x[c]).sum();
    }
}

/// Block-Jacobi preconditioner: the inverted diagonal blocks.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockJacobi {
    b: usize,
    inverse: Vec<f64>,
}

impl BlockJacobi {
    /// The block size.
    #[must_use]
    pub const fn block_size(&self) -> usize {
        self.b
    }

    /// `z = M^-1 r`, block by block.
    ///
    /// # Panics
    /// A slice shorter than the matrix extent can panic.
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let (b, bs) = (self.b, self.b * self.b);
        for (i, inv) in self.inverse.chunks_exact(bs).enumerate() {
            block_apply(inv, b, &r[i * b..(i + 1) * b], &mut z[i * b..(i + 1) * b]);
        }
    }
}

/// Block ILU(0) factors; see [`Bsr::bilu0`].
///
/// The strictly lower blocks hold `L_ik = A_ik U_kk^-1` (unit block diagonal),
/// the diagonal and upper blocks hold `U`.
#[derive(Debug, Clone, PartialEq)]
pub struct Bilu0 {
    lu: Bsr,
    diag: Vec<usize>,
    inv_diag: Vec<f64>,
}

impl Bilu0 {
    /// `z = (L U)^-1 r`: block forward then backward substitution.
    ///
    /// # Panics
    /// A slice shorter than the matrix extent can panic.
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let b = self.lu.br;
        let bs = b * b;
        let nbr = self.lu.block_rows();
        if nbr == 0 {
            return;
        }
        let mut acc = super::reuse::f64_scratch(b);
        acc.resize(b, 0.0);
        // L y = r (y stored in z)
        for i in 0..nbr {
            acc.copy_from_slice(&r[i * b..(i + 1) * b]);
            for k in self.lu.indptr[i]..self.diag[i] {
                let j = self.lu.indices[k];
                let blk = &self.lu.values[k * bs..(k + 1) * bs];
                for p in 0..b {
                    for q in 0..b {
                        acc[p] -= blk[p * b + q] * z[j * b + q];
                    }
                }
            }
            z[i * b..(i + 1) * b].copy_from_slice(&acc);
        }
        // U z = y
        for i in (0..nbr).rev() {
            acc.copy_from_slice(&z[i * b..(i + 1) * b]);
            for k in (self.diag[i] + 1)..self.lu.indptr[i + 1] {
                let j = self.lu.indices[k];
                let blk = &self.lu.values[k * bs..(k + 1) * bs];
                for p in 0..b {
                    for q in 0..b {
                        acc[p] -= blk[p * b + q] * z[j * b + q];
                    }
                }
            }
            block_apply(
                &self.inv_diag[i * bs..(i + 1) * bs],
                b,
                &acc,
                &mut z[i * b..(i + 1) * b],
            );
        }
    }
}

/// Number of distinct `br x bc` blocks holding a nonzero.
#[allow(
    clippy::single_call_fn,
    reason = "keep the candidate traversal separate from ranking"
)]
fn count_blocks(
    csr: &Csr,
    br: usize,
    bc: usize,
) -> Result<usize, BsrError> {
    let nbr = csr.rows / br;
    let mut mark = BlockSlots::new(csr.cols / bc, csr.values.len())?;
    let mut columns = Vec::new();
    let mut blocks = 0;
    for ib in 0..nbr {
        for r in 0..br {
            let row = ib * br + r;
            for k in csr.indptr[row]..csr.indptr[row + 1] {
                if csr.values[k] != 0.0 {
                    let jb = csr.indices[k] / bc;
                    if mark.get(jb).is_none() {
                        mark.insert(jb, 0)?;
                        columns.try_reserve(1).map_err(|_| size_error())?;
                        columns.push(jb);
                    }
                }
            }
        }
        blocks += columns.len();
        mark.clear(columns.iter().copied());
        columns.clear();
    }
    Ok(blocks)
}

/// Picks the block shape for `csr` among `candidates`.
///
/// A candidate qualifies when it divides the matrix shape and its fill ratio
/// (stored scalars per true nonzero, see [`Bsr::fill_ratio`]) is at most
/// [`MAX_FILL_RATIO`]. The largest qualifying block wins, ties going to the
/// lower fill ratio. The result is `(1, 1)` for malformed CSR arrays or
/// when no candidate qualifies.
#[must_use]
pub fn choose_block_size(
    csr: &Csr,
    candidates: &[(usize, usize)],
) -> (usize, usize) {
    let mut best = (1usize, 1usize);
    let mut best_key = (1usize, f64::INFINITY);
    if validate_csr(csr).is_err() {
        return best;
    }
    let nnz = csr.values.iter().filter(|&&value| value != 0.0).count();
    if nnz == 0 {
        return best;
    }
    for &(br, bc) in candidates {
        if check_block(csr.rows, csr.cols, br, bc).is_err() {
            continue;
        }
        let area = br * bc;
        if area < best_key.0 {
            continue;
        }
        let Ok(blocks) = count_blocks(csr, br, bc) else {
            continue;
        };
        let ratio = (blocks as f64 * area as f64) / nnz as f64;
        if ratio <= MAX_FILL_RATIO
            && (area > best_key.0 || (area == best_key.0 && ratio < best_key.1))
        {
            best = (br, bc);
            best_key = (area, ratio);
        }
    }
    best
}

#[cfg(test)]
#[path = "bsr/tests.rs"]
mod tests;
