//! Sparse product by independent, ordered groups of block rows.

use super::BlockSlots;
use super::Bsr;
use super::BsrError;
use super::check_block;
use super::extend_zeros;
use super::filled;
use super::gemm_acc;
use super::size_error;
use crate::kernels::par;
use std::ops::Range;

const MIN_BLOCK_ROWS: usize = 64;
const PAR_WORK: usize = 1_000_000;

#[allow(
    clippy::single_call_fn,
    reason = "keep scheduling and row assembly together in this module"
)]
pub(super) fn matmul(
    left: &Bsr,
    right: &Bsr,
) -> Result<Bsr, BsrError> {
    check_block(left.nrows, right.ncols, left.br, right.bc)?;
    if left.indices.is_empty() || right.indices.is_empty() {
        return Bsr::zeros(left.nrows, right.ncols, left.br, right.bc);
    }

    let rows = left.block_rows();
    let enough_work = || {
        let block_work = left.br.saturating_mul(left.bc).saturating_mul(right.bc);
        let needed = PAR_WORK.div_ceil(block_work);
        left.indices
            .iter()
            .try_fold(0usize, |seen, &middle| {
                let next = seen.saturating_add(right.indptr[middle + 1] - right.indptr[middle]);
                (next < needed).then_some(next)
            })
            .is_none()
    };

    let output = if rows >= MIN_BLOCK_ROWS * 2
        && enough_work()
        && crate::execution::should_split(rows, MIN_BLOCK_ROWS * 2)
    {
        let groups = crate::execution::cpu_threads().saturating_mul(2);
        let chunk = rows.div_ceil(groups).max(MIN_BLOCK_ROWS);
        let mut pieces = par::map_indexed(rows.div_ceil(chunk), 2, |part| {
            let first = part * chunk;
            Rows::multiply(left, right, first..first.saturating_add(chunk).min(rows))
        })
        .into_iter();
        let mut result = pieces.next().ok_or_else(size_error)??;
        for piece in pieces {
            result.append(piece?)?;
        }
        result
    } else {
        Rows::multiply(left, right, 0..rows)?
    };

    Ok(Bsr {
        nrows: left.nrows,
        ncols: right.ncols,
        br: left.br,
        bc: right.bc,
        indptr: output.indptr,
        indices: output.indices,
        values: output.values,
    })
}

struct Rows {
    indptr: Vec<usize>,
    indices: Vec<usize>,
    values: Vec<f64>,
}

impl Rows {
    fn multiply(
        left: &Bsr,
        right: &Bsr,
        rows: Range<usize>,
    ) -> Result<Self, BsrError> {
        let (m, k, n) = (left.br, left.bc, right.bc);
        let (left_bs, right_bs, out_bs) = (m * k, k * n, m * n);
        let mut output = Self {
            indptr: filled(rows.len().checked_add(1).ok_or_else(size_error)?, 0)?,
            indices: Vec::new(),
            values: Vec::new(),
        };
        let mut slots = BlockSlots::new(right.block_cols(), right.block_count())?;
        let mut acc = Vec::new();
        let mut columns = Vec::new();

        for (local, row) in rows.enumerate() {
            acc.clear();
            columns.clear();
            for a_idx in left.indptr[row]..left.indptr[row + 1] {
                let middle = left.indices[a_idx];
                let a = &left.values[a_idx * left_bs..(a_idx + 1) * left_bs];
                for b_idx in right.indptr[middle]..right.indptr[middle + 1] {
                    let col = right.indices[b_idx];
                    let slot = if let Some(slot) = slots.get(col) {
                        slot
                    } else {
                        let slot = columns.len();
                        columns.try_reserve(1).map_err(|_| size_error())?;
                        columns.push((col, slot));
                        slots.insert(col, slot)?;
                        extend_zeros(
                            &mut acc,
                            columns.len().checked_mul(out_bs).ok_or_else(size_error)?,
                        )?;
                        slot
                    };
                    let b = &right.values[b_idx * right_bs..(b_idx + 1) * right_bs];
                    let c = &mut acc[slot * out_bs..(slot + 1) * out_bs];
                    gemm_acc(c, a, b, k, n, 1.0);
                }
            }
            columns.sort_unstable();
            for &(col, slot) in &columns {
                let block = &acc[slot * out_bs..(slot + 1) * out_bs];
                if block.iter().any(|&value| value != 0.0) {
                    output.indices.try_reserve(1).map_err(|_| size_error())?;
                    output
                        .values
                        .try_reserve(out_bs)
                        .map_err(|_| size_error())?;
                    output.indices.push(col);
                    output.values.extend_from_slice(block);
                }
            }
            slots.clear(columns.iter().map(|&(col, _)| col));
            output.indptr[local + 1] = output.indices.len();
        }
        Ok(output)
    }

    fn append(
        &mut self,
        mut other: Self,
    ) -> Result<(), BsrError> {
        let first = self.indices.len();
        self.indptr
            .try_reserve(other.indptr.len() - 1)
            .map_err(|_| size_error())?;
        self.indices
            .try_reserve(other.indices.len())
            .map_err(|_| size_error())?;
        self.values
            .try_reserve(other.values.len())
            .map_err(|_| size_error())?;
        self.indptr
            .extend(other.indptr.into_iter().skip(1).map(|ptr| ptr + first));
        self.indices.append(&mut other.indices);
        self.values.append(&mut other.values);
        Ok(())
    }
}
