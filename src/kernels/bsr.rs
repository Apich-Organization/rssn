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
//! * [`Bsr::transpose`], [`Bsr::matmul`], [`Bsr::diagonal_blocks`].
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
    /// A representation invariant is violated.
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

const fn check_block(
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
    Ok(())
}

impl Bsr {
    /// An all-zero matrix with no stored blocks.
    ///
    /// # Errors
    /// [`BsrError::ZeroBlock`] or [`BsrError::NotDivisible`].
    pub fn zeros(
        nrows: usize,
        ncols: usize,
        br: usize,
        bc: usize,
    ) -> Result<Self, BsrError> {
        check_block(nrows, ncols, br, bc)?;
        Ok(Self {
            nrows,
            ncols,
            br,
            bc,
            indptr: vec![0; nrows / br + 1],
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
        if self.indptr.len() != nbr + 1 || self.indptr[0] != 0 || self.indptr[nbr] != nb {
            return Err("block-row pointers do not match the indices".into());
        }
        if self.values.len() != nb * self.br * self.bc {
            return Err("values do not hold br * bc scalars per block".into());
        }
        for i in 0..nbr {
            if self.indptr[i] > self.indptr[i + 1] {
                return Err("block-row pointers are not monotone".into());
            }
            for k in self.indptr[i]..self.indptr[i + 1] {
                if self.indices[k] >= nbc {
                    return Err("block-column index out of range".into());
                }
                if k > self.indptr[i] && self.indices[k - 1] >= self.indices[k] {
                    return Err("block-column indices are not strictly increasing".into());
                }
            }
        }
        Ok(())
    }

    /// Converts a CSR matrix; a block is stored iff any of its entries is
    /// nonzero, and the zeros inside a stored block are padding.
    ///
    /// # Errors
    /// [`BsrError::ZeroBlock`] or [`BsrError::NotDivisible`].
    pub fn from_csr(
        csr: &Csr,
        br: usize,
        bc: usize,
    ) -> Result<Self, BsrError> {
        check_block(csr.rows, csr.cols, br, bc)?;
        let nbr = csr.rows / br;
        let nbc = csr.cols / bc;
        let bs = br * bc;
        let mut indptr = vec![0usize; nbr + 1];
        let mut indices: Vec<usize> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        // `slot[jb]` is the position of block column `jb` in the current block
        // row, or `usize::MAX`.
        let mut slot = vec![usize::MAX; nbc];
        for ib in 0..nbr {
            let first = indices.len();
            for r in 0..br {
                let row = ib * br + r;
                for k in csr.indptr[row]..csr.indptr[row + 1] {
                    if csr.values[k] != 0.0 {
                        let jb = csr.indices[k] / bc;
                        if slot[jb] == usize::MAX {
                            slot[jb] = 0;
                            indices.push(jb);
                        }
                    }
                }
            }
            indices[first..].sort_unstable();
            values.resize(indices.len() * bs, 0.0);
            for (offset, jb) in indices[first..].iter().enumerate() {
                slot[*jb] = first + offset;
            }
            for r in 0..br {
                let row = ib * br + r;
                for k in csr.indptr[row]..csr.indptr[row + 1] {
                    if csr.values[k] != 0.0 {
                        let col = csr.indices[k];
                        let at = slot[col / bc];
                        values[at * bs + r * bc + col % bc] = csr.values[k];
                    }
                }
            }
            for jb in &indices[first..] {
                slot[*jb] = usize::MAX;
            }
            indptr[ib + 1] = indices.len();
        }
        Ok(Self {
            nrows: csr.rows,
            ncols: csr.cols,
            br,
            bc,
            indptr,
            indices,
            values,
        })
    }

    /// Builds a matrix from `(row, col, value)` triplets (duplicates summed).
    ///
    /// # Errors
    /// [`BsrError::ZeroBlock`] or [`BsrError::NotDivisible`].
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        br: usize,
        bc: usize,
        trip: &[(usize, usize, f64)],
    ) -> Result<Self, BsrError> {
        check_block(nrows, ncols, br, bc)?;
        Self::from_csr(&Csr::from_triplets(nrows, ncols, trip), br, bc)
    }

    /// Converts to CSR, dropping the zeros inside blocks.
    #[must_use]
    pub fn to_csr(&self) -> Csr {
        let nbr = self.block_rows();
        let bs = self.br * self.bc;
        let mut indptr = vec![0usize; self.nrows + 1];
        let mut indices = Vec::with_capacity(self.nnz());
        let mut values = Vec::with_capacity(self.nnz());
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

    /// One output row of the product: the correctly rounded dot product.
    fn row_dot(
        &self,
        ib: usize,
        r: usize,
        x: &[f64],
    ) -> f64 {
        let (bs, bc) = (self.br * self.bc, self.bc);
        let (values, indices) = (&self.values, &self.indices);
        dot_exact((self.indptr[ib]..self.indptr[ib + 1]).flat_map(move |k| {
            (0..bc).map(move |c| (values[k * bs + r * bc + c], x[indices[k] * bc + c]))
        }))
    }

    /// Matrix-vector product `y = A x`.
    ///
    /// Every entry of `y` is the correctly rounded dot product of its row, so
    /// the result equals [`Csr::matvec`] of [`Bsr::to_csr`] bit for bit,
    /// whatever the block shape; block rows run in parallel for large
    /// matrices without changing a bit.
    pub fn matvec(
        &self,
        x: &[f64],
        y: &mut [f64],
    ) {
        let rows = self.nrows.min(y.len());
        if self.values.len() >= BSR_PAR_NNZ {
            let chunk = BSR_CHUNK_ROWS * self.br;
            super::par::for_chunks_mut(&mut y[..rows], chunk, 2, |b, block| {
                for (r, slot) in block.iter_mut().enumerate() {
                    let row = b * chunk + r;
                    *slot = self.row_dot(row / self.br, row % self.br, x);
                }
            });
        } else {
            for row in 0..rows {
                y[row] = self.row_dot(row / self.br, row % self.br, x);
            }
        }
    }

    /// Transposed product `y = A^T x` (`y` has `ncols` entries).
    ///
    /// Accumulates block by block in plain `f64`, so it is *not* required to
    /// match `transpose().matvec` bit for bit.
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
                        y[jb * self.bc + c] += self.values[k * bs + r * self.bc + c] * xr;
                    }
                }
            }
        }
    }

    /// The transpose, with block shape `(bc, br)`.
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
    /// Blocks that cancel to exactly zero are not stored.
    ///
    /// # Errors
    /// [`BsrError::ShapeMismatch`] when the inner extents or block sizes differ.
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
        let (m, k, n) = (self.br, self.bc, other.bc);
        let nbr = self.block_rows();
        let nbc = other.block_cols();
        let out_bs = m * n;
        let mut indptr = vec![0usize; nbr + 1];
        let mut indices: Vec<usize> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        let mut slot = vec![usize::MAX; nbc];
        let mut acc: Vec<f64> = Vec::new();
        let mut cols: Vec<usize> = Vec::new();
        for ib in 0..nbr {
            acc.clear();
            cols.clear();
            for ka in self.indptr[ib]..self.indptr[ib + 1] {
                let mid = self.indices[ka];
                let a = &self.values[ka * m * k..(ka + 1) * m * k];
                for kb in other.indptr[mid]..other.indptr[mid + 1] {
                    let jb = other.indices[kb];
                    if slot[jb] == usize::MAX {
                        slot[jb] = cols.len();
                        cols.push(jb);
                        acc.resize(cols.len() * out_bs, 0.0);
                    }
                    let b = &other.values[kb * k * n..(kb + 1) * k * n];
                    let c = &mut acc[slot[jb] * out_bs..(slot[jb] + 1) * out_bs];
                    gemm_acc(c, a, b, m, k, n, 1.0);
                }
            }
            let mut order: Vec<usize> = (0..cols.len()).collect();
            order.sort_unstable_by_key(|p| cols[*p]);
            for p in order {
                let block = &acc[p * out_bs..(p + 1) * out_bs];
                if block.iter().any(|v| *v != 0.0) {
                    indices.push(cols[p]);
                    values.extend_from_slice(block);
                }
            }
            for jb in &cols {
                slot[*jb] = usize::MAX;
            }
            indptr[ib + 1] = indices.len();
        }
        Ok(Self {
            nrows: self.nrows,
            ncols: other.ncols,
            br: m,
            bc: n,
            indptr,
            indices,
            values,
        })
    }

    /// The square diagonal blocks, one dense row-major `b x b` block per block
    /// row (all zero when the block is not stored).
    ///
    /// # Errors
    /// [`BsrError::NotSquare`] unless the matrix and its blocks are square.
    pub fn diagonal_blocks(&self) -> Result<Vec<Vec<f64>>, BsrError> {
        if self.nrows != self.ncols || self.br != self.bc {
            return Err(BsrError::NotSquare);
        }
        let bs = self.br * self.bc;
        Ok((0..self.block_rows())
            .map(|ib| {
                let row = &self.indices[self.indptr[ib]..self.indptr[ib + 1]];
                row.binary_search(&ib).map_or_else(
                    |_| vec![0.0; bs],
                    |p| {
                        let k = self.indptr[ib] + p;
                        self.values[k * bs..(k + 1) * bs].to_vec()
                    },
                )
            })
            .collect())
    }

    /// Block-Jacobi preconditioner: the inverse of every diagonal block.
    ///
    /// # Errors
    /// [`BsrError::NotSquare`], or [`BsrError::SingularBlock`] for the first
    /// block row whose diagonal block is singular (or not stored).
    pub fn block_jacobi(&self) -> Result<BlockJacobi, BsrError> {
        let b = self.br;
        let bs = b * b;
        let mut inverse = Vec::with_capacity(self.block_rows() * bs);
        for (ib, block) in self.diagonal_blocks()?.iter().enumerate() {
            let inv = invert(block, b).ok_or(BsrError::SingularBlock { block_row: ib })?;
            inverse.extend_from_slice(&inv);
        }
        Ok(BlockJacobi { b, inverse })
    }

    /// Block ILU(0): incomplete block LU keeping the block sparsity pattern.
    ///
    /// # Errors
    /// [`BsrError::NotSquare`], [`BsrError::MissingDiagonal`] when a diagonal
    /// block is not stored, [`BsrError::SingularBlock`] for a singular pivot.
    pub fn bilu0(&self) -> Result<Bilu0, BsrError> {
        if self.nrows != self.ncols || self.br != self.bc {
            return Err(BsrError::NotSquare);
        }
        let b = self.br;
        let bs = b * b;
        let nbr = self.block_rows();
        let mut lu = self.clone();
        let mut diag = vec![0usize; nbr];
        let mut inv_diag = vec![0.0; nbr * bs];
        // Position of block column `j` inside the current block row.
        let mut pos = vec![usize::MAX; nbr];
        let mut tmp = vec![0.0; bs];
        for i in 0..nbr {
            let (start, end) = (lu.indptr[i], lu.indptr[i + 1]);
            for k in start..end {
                pos[lu.indices[k]] = k;
            }
            let Some(d) = (start..end).find(|k| lu.indices[*k] == i) else {
                return Err(BsrError::MissingDiagonal { block_row: i });
            };
            diag[i] = d;
            for ik in start..d {
                let k = lu.indices[ik];
                // L_ik = A_ik * inv(U_kk)
                tmp.fill(0.0);
                gemm_acc(
                    &mut tmp,
                    &lu.values[ik * bs..(ik + 1) * bs],
                    &inv_diag[k * bs..(k + 1) * bs],
                    b,
                    b,
                    b,
                    1.0,
                );
                lu.values[ik * bs..(ik + 1) * bs].copy_from_slice(&tmp);
                // A_ij -= L_ik * U_kj for j > k inside the pattern of row i.
                for kj in (diag[k] + 1)..lu.indptr[k + 1] {
                    let j = lu.indices[kj];
                    if pos[j] != usize::MAX {
                        let u = lu.values[kj * bs..(kj + 1) * bs].to_vec();
                        let at = pos[j];
                        gemm_acc(
                            &mut lu.values[at * bs..(at + 1) * bs],
                            &tmp,
                            &u,
                            b,
                            b,
                            b,
                            -1.0,
                        );
                    }
                }
            }
            let pivot = &lu.values[d * bs..(d + 1) * bs];
            let inv = invert(pivot, b).ok_or(BsrError::SingularBlock { block_row: i })?;
            inv_diag[i * bs..(i + 1) * bs].copy_from_slice(&inv);
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
    m: usize,
    k: usize,
    n: usize,
    alpha: f64,
) {
    for i in 0..m {
        for p in 0..k {
            let aip = alpha * a[i * k + p];
            if aip != 0.0 {
                for j in 0..n {
                    c[i * n + j] += aip * b[p * n + j];
                }
            }
        }
    }
}

/// LU with partial pivoting of the `n x n` row-major matrix `a`; returns the
/// packed factors and the row permutation, or `None` when singular.
fn lu_factor(
    a: &[f64],
    n: usize,
) -> Option<(Vec<f64>, Vec<usize>)> {
    let mut lu = a.to_vec();
    let scale = lu.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return None;
    }
    let mut perm: Vec<usize> = (0..n).collect();
    for col in 0..n {
        let mut best = col;
        for r in (col + 1)..n {
            if lu[r * n + col].abs() > lu[best * n + col].abs() {
                best = r;
            }
        }
        if lu[best * n + col].abs() <= f64::EPSILON * scale {
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
    Some((lu, perm))
}

/// The inverse of an `n x n` row-major matrix by LU, or `None` if singular.
fn invert(
    a: &[f64],
    n: usize,
) -> Option<Vec<f64>> {
    let (lu, perm) = lu_factor(a, n)?;
    let mut inverse = vec![0.0; n * n];
    let mut col = vec![0.0; n];
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
            inverse[r * n + e] = col[r];
        }
    }
    inverse.iter().all(|v| v.is_finite()).then_some(inverse)
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
    pub fn apply(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let b = self.lu.br;
        let bs = b * b;
        let nbr = self.lu.block_rows();
        let mut tmp = vec![0.0; b];
        let mut acc = vec![0.0; b];
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
            block_apply(&self.inv_diag[i * bs..(i + 1) * bs], b, &acc, &mut tmp);
            z[i * b..(i + 1) * b].copy_from_slice(&tmp);
        }
    }
}

/// Number of distinct `br x bc` blocks holding a nonzero, and the number of
/// nonzeros of `csr`.
fn count_blocks(
    csr: &Csr,
    br: usize,
    bc: usize,
) -> (usize, usize) {
    let nbr = csr.rows / br;
    let mut mark = vec![usize::MAX; csr.cols / bc];
    let (mut blocks, mut nnz) = (0usize, 0usize);
    for ib in 0..nbr {
        for r in 0..br {
            let row = ib * br + r;
            for k in csr.indptr[row]..csr.indptr[row + 1] {
                if csr.values[k] != 0.0 {
                    nnz += 1;
                    let jb = csr.indices[k] / bc;
                    if mark[jb] != ib {
                        mark[jb] = ib;
                        blocks += 1;
                    }
                }
            }
        }
    }
    (blocks, nnz)
}

/// Picks the block shape for `csr` among `candidates`.
///
/// A candidate qualifies when it divides the matrix shape and its fill ratio
/// (stored scalars per true nonzero, see [`Bsr::fill_ratio`]) is at most
/// [`MAX_FILL_RATIO`]. The largest qualifying block wins, ties going to the
/// lower fill ratio; with no qualifying candidate (or no nonzeros) the answer
/// is `(1, 1)`.
#[must_use]
pub fn choose_block_size(
    csr: &Csr,
    candidates: &[(usize, usize)],
) -> (usize, usize) {
    let mut best = (1usize, 1usize);
    let mut best_key = (1usize, f64::INFINITY);
    for &(br, bc) in candidates {
        if br == 0 || bc == 0 || !csr.rows.is_multiple_of(br) || !csr.cols.is_multiple_of(bc) {
            continue;
        }
        let (blocks, nnz) = count_blocks(csr, br, bc);
        if nnz == 0 {
            continue;
        }
        let ratio = (blocks * br * bc) as f64 / nnz as f64;
        let area = br * bc;
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
mod tests {
    use super::*;
    use crate::kernels::krylov::cg;
    use crate::kernels::random::Rng;

    fn random_csr(
        rng: &mut Rng,
        rows: usize,
        cols: usize,
        density: f64,
    ) -> Csr {
        let mut trip = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                if rng.uniform() < density {
                    trip.push((r, c, rng.normal() + 0.1));
                }
            }
        }
        Csr::from_triplets(rows, cols, &trip)
    }

    fn dense(csr: &Csr) -> Vec<f64> {
        let mut out = vec![0.0; csr.rows * csr.cols];
        for r in 0..csr.rows {
            for k in csr.indptr[r]..csr.indptr[r + 1] {
                out[r * csr.cols + csr.indices[k]] = csr.values[k];
            }
        }
        out
    }

    #[test]
    fn roundtrip_csr_bsr_csr() {
        let mut rng = Rng::new(1);
        for (br, bc) in [(1, 1), (2, 3), (3, 2), (4, 3), (12, 18)] {
            let csr = random_csr(&mut rng, 12, 18, 0.2);
            let bsr = Bsr::from_csr(&csr, br, bc).unwrap();
            assert_eq!(bsr.validate(), Ok(()));
            assert_eq!(bsr.to_csr(), csr);
            assert_eq!(bsr.nnz(), csr.values.len());
        }
    }

    #[test]
    fn matvec_equals_csr_bit_for_bit() {
        let mut rng = Rng::new(2);
        let csr = random_csr(&mut rng, 24, 30, 0.15);
        let x: Vec<f64> = (0..30).map(|_| rng.normal()).collect();
        let mut want = vec![0.0; 24];
        csr.matvec(&x, &mut want);
        for (br, bc) in [(1, 1), (2, 3), (4, 5), (3, 6)] {
            let bsr = Bsr::from_csr(&csr, br, bc).unwrap();
            let mut got = vec![9.0; 24];
            bsr.matvec(&x, &mut got);
            assert_eq!(got, want, "block {br}x{bc}");
        }
    }

    #[test]
    fn matvec_parallel_path_is_identical() {
        let mut rng = Rng::new(3);
        let n = 600;
        let csr = random_csr(&mut rng, n, n, 0.2);
        let bsr = Bsr::from_csr(&csr, 4, 4).unwrap();
        assert!(bsr.stored_len() >= BSR_PAR_NNZ);
        let x: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        let (mut a, mut b) = (vec![0.0; n], vec![0.0; n]);
        csr.matvec(&x, &mut a);
        bsr.matvec(&x, &mut b);
        assert_eq!(a, b);
    }

    #[test]
    fn transpose_and_matvec_transpose() {
        let mut rng = Rng::new(4);
        let csr = random_csr(&mut rng, 12, 18, 0.25);
        let bsr = Bsr::from_csr(&csr, 3, 2).unwrap();
        let t = bsr.transpose();
        assert_eq!((t.br, t.bc, t.nrows, t.ncols), (2, 3, 18, 12));
        assert_eq!(t.validate(), Ok(()));
        let d = dense(&csr);
        let dt = dense(&t.to_csr());
        for r in 0..12 {
            for c in 0..18 {
                assert_eq!(d[r * 18 + c], dt[c * 12 + r]);
            }
        }
        assert_eq!(t.transpose(), bsr);
        let x: Vec<f64> = (0..12).map(|_| rng.normal()).collect();
        let (mut y1, mut y2) = (vec![0.0; 18], vec![0.0; 18]);
        bsr.matvec_transpose(&x, &mut y1);
        t.matvec(&x, &mut y2);
        for (a, b) in y1.iter().zip(&y2) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn matmul_matches_dense() {
        let mut rng = Rng::new(5);
        let a = random_csr(&mut rng, 8, 12, 0.3);
        let b = random_csr(&mut rng, 12, 6, 0.3);
        let c = Bsr::from_csr(&a, 2, 3)
            .unwrap()
            .matmul(&Bsr::from_csr(&b, 3, 2).unwrap())
            .unwrap();
        assert_eq!(c.validate(), Ok(()));
        assert_eq!((c.br, c.bc), (2, 2));
        let (da, db, dc) = (dense(&a), dense(&b), dense(&c.to_csr()));
        for i in 0..8 {
            for j in 0..6 {
                let want: f64 = (0..12).map(|k| da[i * 12 + k] * db[k * 6 + j]).sum();
                assert!((dc[i * 6 + j] - want).abs() < 1e-12);
            }
        }
        let bad = Bsr::from_csr(&b, 1, 2).unwrap();
        assert!(matches!(
            Bsr::from_csr(&a, 2, 3).unwrap().matmul(&bad),
            Err(BsrError::ShapeMismatch(_))
        ));
    }

    /// A block-tridiagonal matrix of `n` block rows with random `b x b` blocks.
    fn block_tridiagonal(
        rng: &mut Rng,
        n: usize,
        b: usize,
    ) -> Csr {
        let mut trip = Vec::new();
        for i in 0..n {
            for j in [i.wrapping_sub(1), i, i + 1] {
                if j >= n {
                    continue;
                }
                for r in 0..b {
                    for c in 0..b {
                        let mut v = rng.normal();
                        if i == j && r == c {
                            v += 6.0;
                        }
                        trip.push((i * b + r, j * b + c, v));
                    }
                }
            }
        }
        Csr::from_triplets(n * b, n * b, &trip)
    }

    #[test]
    fn bilu0_is_exact_without_fill() {
        let mut rng = Rng::new(6);
        let (n, b) = (7, 3);
        let csr = block_tridiagonal(&mut rng, n, b);
        let bsr = Bsr::from_csr(&csr, b, b).unwrap();
        let ilu = bsr.bilu0().unwrap();
        let x: Vec<f64> = (0..n * b).map(|_| rng.normal()).collect();
        let mut rhs = vec![0.0; n * b];
        bsr.matvec(&x, &mut rhs);
        let mut z = vec![0.0; n * b];
        ilu.apply(&rhs, &mut z);
        for (a, e) in z.iter().zip(&x) {
            assert!((a - e).abs() < 1e-10, "{a} vs {e}");
        }
    }

    #[test]
    fn block_jacobi_inverts_block_diagonal_matrix() {
        let mut rng = Rng::new(7);
        let csr = block_tridiagonal(&mut rng, 1, 4);
        let bsr = Bsr::from_csr(&csr, 4, 4).unwrap();
        let pre = bsr.block_jacobi().unwrap();
        let x: Vec<f64> = (0..4).map(|_| rng.normal()).collect();
        let mut rhs = vec![0.0; 4];
        bsr.matvec(&x, &mut rhs);
        let mut z = vec![0.0; 4];
        pre.apply(&rhs, &mut z);
        for (a, e) in z.iter().zip(&x) {
            assert!((a - e).abs() < 1e-10);
        }
        let singular = Bsr::from_triplets(
            2,
            2,
            2,
            2,
            &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 4.0)],
        )
        .unwrap();
        assert_eq!(
            singular.block_jacobi(),
            Err(BsrError::SingularBlock { block_row: 0 })
        );
        let missing = Bsr::from_triplets(
            4,
            4,
            2,
            2,
            &[(0, 0, 1.0), (1, 1, 1.0), (0, 2, 1.0), (2, 0, 1.0)],
        )
        .unwrap();
        assert_eq!(
            missing.bilu0(),
            Err(BsrError::MissingDiagonal { block_row: 1 })
        );
    }

    #[test]
    fn shape_validation() {
        let csr = Csr::from_triplets(5, 4, &[(0, 0, 1.0)]);
        assert_eq!(
            Bsr::from_csr(&csr, 2, 2),
            Err(BsrError::NotDivisible { extent: 5, block: 2 })
        );
        assert_eq!(Bsr::from_csr(&csr, 0, 2), Err(BsrError::ZeroBlock));
        let mut bsr = Bsr::from_csr(&Csr::from_triplets(4, 4, &[(0, 0, 1.0)]), 2, 2).unwrap();
        assert_eq!(bsr.validate(), Ok(()));
        bsr.indices[0] = 9;
        assert!(bsr.validate().is_err());
    }

    #[test]
    fn empty_matrix() {
        let empty = Csr::from_triplets(0, 0, &[]);
        let bsr = Bsr::from_csr(&empty, 2, 2).unwrap();
        assert_eq!(bsr.block_count(), 0);
        assert_eq!(bsr.validate(), Ok(()));
        assert_eq!(bsr.to_csr(), empty);
        let zero = Bsr::from_csr(&Csr::from_triplets(4, 6, &[]), 2, 3).unwrap();
        assert_eq!(zero.block_count(), 0);
        assert!((zero.fill_ratio() - 1.0).abs() < f64::EPSILON);
        let mut y = vec![1.0; 4];
        zero.matvec(&[0.0; 6], &mut y);
        assert_eq!(y, vec![0.0; 4]);
        assert_eq!(zero.transpose().validate(), Ok(()));
    }

    #[test]
    fn unit_blocks_are_csr() {
        let mut rng = Rng::new(8);
        let csr = random_csr(&mut rng, 10, 10, 0.3);
        let bsr = Bsr::from_csr(&csr, 1, 1).unwrap();
        assert_eq!(bsr.indptr, csr.indptr);
        assert_eq!(bsr.indices, csr.indices);
        assert_eq!(bsr.values, csr.values);
        assert!((bsr.fill_ratio() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn block_size_selection() {
        // 2x2 nodal blocks on a chain: dense blocks, ratio 1 (4x4 would store 1.86x).
        let mut rng = Rng::new(9);
        let csr = block_tridiagonal(&mut rng, 10, 2);
        let cands = [(4, 4), (3, 3), (2, 2), (1, 1)];
        assert_eq!(choose_block_size(&csr, &cands), (2, 2));
        let bsr = Bsr::from_csr(&csr, 2, 2).unwrap();
        assert!((bsr.fill_ratio() - 1.0).abs() < 1e-12);
        // A scalar tridiagonal matrix gains nothing from blocks.
        let trip: Vec<_> = (0..12)
            .flat_map(|i| [(i, i, 2.0), (i, (i + 1) % 12, -1.0)])
            .collect();
        let scalar = Csr::from_triplets(12, 12, &trip);
        assert_eq!(choose_block_size(&scalar, &cands), (1, 1));
        assert_eq!(
            choose_block_size(&Csr::from_triplets(4, 4, &[]), &cands),
            (1, 1)
        );
    }

    /// `P (x) S` with the 2-D five-point Laplacian `P` on an `m x m` grid and a
    /// dense SPD `s x s` nodal coupling `S`: scalar Jacobi cannot see `S`.
    fn nodal_system(
        m: usize,
        s: usize,
        coupling: &[f64],
    ) -> Csr {
        let mut trip = Vec::new();
        for gy in 0..m {
            for gx in 0..m {
                let node = gy * m + gx;
                let mut nbrs = vec![(node, 4.0)];
                if gx > 0 {
                    nbrs.push((node - 1, -1.0));
                }
                if gx + 1 < m {
                    nbrs.push((node + 1, -1.0));
                }
                if gy > 0 {
                    nbrs.push((node - m, -1.0));
                }
                if gy + 1 < m {
                    nbrs.push((node + m, -1.0));
                }
                for (other, w) in nbrs {
                    for a in 0..s {
                        for b in 0..s {
                            trip.push((node * s + a, other * s + b, w * coupling[a * s + b]));
                        }
                    }
                }
            }
        }
        Csr::from_triplets(m * m * s, m * m * s, &trip)
    }

    fn cg_iterations(
        csr: &Csr,
        pre: &dyn Fn(&[f64], &mut [f64]),
    ) -> usize {
        let n = csr.rows;
        let b: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64).collect();
        let result = cg(|x, y| csr.matvec(x, y), pre, &b, None, 1e-8, 5000);
        assert!(result.converged);
        result.iterations
    }

    #[test]
    fn block_jacobi_beats_scalar_jacobi_in_cg() {
        let couplings: [(usize, Vec<f64>); 2] = [
            (2, vec![1.0, 0.95, 0.95, 1.0]),
            (3, vec![1.0, 0.9, 0.9, 0.9, 1.0, 0.9, 0.9, 0.9, 1.0]),
        ];
        for (s, coupling) in couplings {
            let csr = nodal_system(12, s, &coupling);
            let inv = csr.jacobi_inverse();
            let scalar = cg_iterations(&csr, &|r, z| {
                for i in 0..r.len() {
                    z[i] = r[i] * inv[i];
                }
            });
            let bsr = Bsr::from_csr(&csr, s, s).unwrap();
            assert!((bsr.fill_ratio() - 1.0).abs() < 1e-12);
            let pre = bsr.block_jacobi().unwrap();
            let block = cg_iterations(&csr, &|r, z| pre.apply(r, z));
            assert!(block < scalar, "s={s}: block {block} vs scalar {scalar}");
            let ilu = bsr.bilu0().unwrap();
            let bilu = cg_iterations(&csr, &|r, z| ilu.apply(r, z));
            assert!(bilu < block, "s={s}: bilu {bilu} vs block {block}");
        }
    }
}
