//! # Deterministic data-parallel helpers
//!
//! Thin wrappers over `rayon` that keep the crate's reproducibility rules in
//! one place:
//!
//! * results are **independent of the thread count**: maps collect in index
//!   order, and floating-point reductions fold fixed-size blocks in order,
//!   so `RAYON_NUM_THREADS=1` and 16 threads
//!   give the same bits;
//! * every helper has a **sequential path** below a grain-size threshold
//!   (small problems must not get slower);
//! * only the outermost loop of a kernel should call these; a helper called
//!   from inside a Rayon task runs on the same pool. Work items inherit the
//!   active cancellation scope.

use std::sync::Mutex;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rayon::prelude::*;

fn indexed_map<T, F>(
    n: usize,
    f: F,
) -> impl IndexedParallelIterator<Item = T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    let cancellation = crate::cancellation::current();
    (0..n).into_par_iter().map_init(
        move || crate::cancellation::install_optional(cancellation.clone()),
        move |_, index| f(index),
    )
}

/// Maps `0..n` through `f` and collects in index order.
///
/// Runs sequentially when `n < grain` (or `n < 2`).
pub fn map_indexed<T, F>(
    n: usize,
    grain: usize,
    f: F,
) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    if crate::execution::should_split(n, grain) {
        indexed_map(n, f).collect()
    } else {
        (0..n).map(f).collect()
    }
}

/// Fallibly maps `0..n` in index order, reporting the error at the lowest
/// failing index. Workers skip work after a known failure when possible.
pub(crate) fn try_map_indexed<T, E, F>(
    n: usize,
    grain: usize,
    f: F,
) -> Result<Vec<T>, E>
where
    T: Send,
    E: Send,
    F: Fn(usize) -> Result<T, E> + Send + Sync,
{
    if !crate::execution::should_split(n, grain) {
        return (0..n).map(f).collect();
    }
    let first_error = AtomicUsize::new(n);
    let results: Vec<_> = indexed_map(n, |index| {
        if index > first_error.load(Ordering::Relaxed) {
            return None;
        }
        let result = f(index);
        if result.is_err() {
            first_error.fetch_min(index, Ordering::Relaxed);
        }
        Some(result)
    })
    .collect();
    // Only indices after an error can be absent. Collect stops at the first
    // error, so a partial result is never returned as a success.
    results.into_iter().flatten().collect()
}

/// Replaces `output` in index order, retaining its allocation when possible.
pub(crate) fn map_indexed_into<T, F>(
    n: usize,
    grain: usize,
    output: &mut Vec<T>,
    f: F,
) where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    if crate::execution::should_split(n, grain) {
        indexed_map(n, f).collect_into_vec(output);
    } else {
        output.clear();
        output.extend((0..n).map(f));
    }
}

/// Fills separate outputs without materializing an intermediate vector of pairs.
pub(crate) fn unzip_indexed_into<A, B, F>(
    n: usize,
    grain: usize,
    first: &mut Vec<A>,
    second: &mut Vec<B>,
    f: F,
) where
    A: Send,
    B: Send,
    F: Fn(usize) -> (A, B) + Sync + Send,
{
    if crate::execution::should_split(n, grain) {
        indexed_map(n, f).unzip_into_vecs(first, second);
    } else {
        first.clear();
        second.clear();
        first.reserve(n);
        second.reserve(n);
        for (a, b) in (0..n).map(f) {
            first.push(a);
            second.push(b);
        }
    }
}

/// Maps the elements of `items` through `f` and collects in order.
pub fn map_slice<A, T, F>(
    items: &[A],
    grain: usize,
    f: F,
) -> Vec<T>
where
    A: Sync,
    T: Send,
    F: Fn(&A) -> T + Sync + Send,
{
    if crate::execution::should_split(items.len(), grain) {
        let cancellation = crate::cancellation::current();
        items
            .par_iter()
            .map_init(
                || crate::cancellation::install_optional(cancellation.clone()),
                |_, item| f(item),
            )
            .collect()
    } else {
        items.iter().map(f).collect()
    }
}

/// Runs `f(i, row)` over mutable chunks of `data` of `chunk` elements each (`i` is the chunk
/// index); sequential when there are fewer than `grain` chunks.
///
/// The chunks are disjoint, so the result cannot depend on the scheduling.
pub fn for_chunks_mut<T, F>(
    data: &mut [T],
    chunk: usize,
    grain: usize,
    f: F,
) where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync + Send,
{
    let chunk = chunk.max(1);
    if crate::execution::should_split(data.len().div_ceil(chunk), grain) {
        let cancellation = crate::cancellation::current();
        data.par_chunks_mut(chunk).enumerate().for_each_init(
            || crate::cancellation::install_optional(cancellation.clone()),
            |_, (i, c)| f(i, c),
        );
    } else {
        data.chunks_mut(chunk)
            .enumerate()
            .for_each(|(i, c)| f(i, c));
    }
}

fn map_reduce_blocks<T, M, R, F>(
    nb: usize,
    minimum: usize,
    init: R,
    map_block: M,
    fold: F,
) -> R
where
    T: Send,
    M: Fn(usize) -> T + Sync + Send,
    F: FnMut(R, T) -> R,
{
    if crate::execution::should_split(nb, minimum) {
        indexed_map(nb, map_block)
            .collect::<Vec<_>>()
            .into_iter()
            .fold(init, fold)
    } else {
        (0..nb).map(map_block).fold(init, fold)
    }
}

const fn neumaier_step(
    sum: f64,
    comp: f64,
    value: f64,
) -> (f64, f64) {
    let next = sum + value;
    let delta = if next.is_finite() {
        if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        }
    } else {
        0.0
    };
    (next, comp + delta)
}

#[allow(clippy::single_call_fn)]
fn dot_block(
    a: &[f64],
    b: &[f64],
) -> f64 {
    crate::kernels::simd::dot_f64(a, b)
}

#[allow(clippy::single_call_fn)]
fn norm2_scaled_fixed(
    a: &[f64],
    block: usize,
    grain: usize,
) -> f64 {
    let n = a.len();
    let nb = n.div_ceil(block);
    let minimum = grain.div_ceil(block).max(2);
    let peak = map_reduce_blocks(
        nb,
        minimum,
        0.0_f64,
        |b_idx| {
            let lo = b_idx.saturating_mul(block);
            let hi = lo.saturating_add(block).min(n);
            a.get(lo..hi)
                .map_or(0.0, crate::kernels::simd::linf_norm_f64)
        },
        |acc, value| {
            if acc.is_nan() || value.is_nan() {
                f64::NAN
            } else {
                acc.max(value)
            }
        },
    );
    if peak == 0.0 || !peak.is_finite() {
        return peak;
    }
    let scaled_sq = map_reduce_blocks(
        nb,
        minimum,
        0.0_f64,
        |b_idx| {
            let lo = b_idx.saturating_mul(block);
            let hi = lo.saturating_add(block).min(n);
            a.get(lo..hi).map_or(0.0, |slice| {
                let mut acc0 = 0.0_f64;
                let mut acc1 = 0.0_f64;
                let mut acc2 = 0.0_f64;
                let mut acc3 = 0.0_f64;
                let (chunks, chunk_tail) = slice.as_chunks::<4>();
                for &[x0, x1, x2, x3] in chunks {
                    let s0 = x0 / peak;
                    let s1 = x1 / peak;
                    let s2 = x2 / peak;
                    let s3 = x3 / peak;
                    acc0 += s0 * s0;
                    acc1 += s1 * s1;
                    acc2 += s2 * s2;
                    acc3 += s3 * s3;
                }
                let mut tail = 0.0_f64;
                for &x in chunk_tail {
                    let s = x / peak;
                    tail += s * s;
                }
                ((acc0 + acc1) + (acc2 + acc3)) + tail
            })
        },
        |total, value| total + value,
    );
    peak * scaled_sq.sqrt()
}

/// Sum of `f(i)` over `0..n` with a fixed association order: the terms are added sequentially
/// inside blocks of `block` indices and the block sums are then added sequentially in block order.
///
/// The result is bit-identical for every thread count (and for the sequential path).
pub fn sum_fixed<F>(
    n: usize,
    block: usize,
    grain: usize,
    f: F,
) -> f64
where
    F: Fn(usize) -> f64 + Sync + Send,
{
    let block = block.max(1);
    let nb = n.div_ceil(block);
    let minimum = grain.div_ceil(block).max(2);
    map_reduce_blocks(
        nb,
        minimum,
        0.0,
        |b| {
            let lo = b.saturating_mul(block);
            let hi = lo.saturating_add(block).min(n);
            let mut s = 0.0;
            for i in lo..hi {
                s += f(i);
            }
            s
        },
        |total, value| total + value,
    )
}

/// Neumaier-compensated sum of `f(i)` over `0..n` using fixed blocks of size `block`.
///
/// Each block accumulates a `(sum, compensation)` pair deterministically and block pairs are
/// folded sequentially in block order with Neumaier compensation, so serial execution
/// (`grain = usize::MAX`) and parallel execution produce bit-identical `f64` values.
#[must_use]
pub fn sum_kahan_fixed<F>(
    n: usize,
    block: usize,
    grain: usize,
    f: F,
) -> f64
where
    F: Fn(usize) -> f64 + Sync + Send,
{
    let block = block.max(1);
    let nb = n.div_ceil(block);
    let minimum = grain.div_ceil(block).max(2);
    let (sum, comp) = map_reduce_blocks(
        nb,
        minimum,
        (0.0_f64, 0.0_f64),
        |b| {
            let lo = b.saturating_mul(block);
            let hi = lo.saturating_add(block).min(n);
            let mut s = 0.0_f64;
            let mut c = 0.0_f64;
            for i in lo..hi {
                (s, c) = neumaier_step(s, c, f(i));
            }
            (s, c)
        },
        |(acc_s, acc_c), (bs, bc)| {
            let (s1, c1) = neumaier_step(acc_s, acc_c, bs);
            neumaier_step(s1, c1, bc)
        },
    );
    sum + comp
}

/// Deterministic dot product of `a` and `b` over their common prefix `a.len().min(b.len())`.
///
/// Each fixed block of `block` elements accumulates with a 4-way unrolled sum, and block sums
/// are folded sequentially in block order so serial execution (`grain = usize::MAX`) and
/// parallel execution produce bit-identical results.
#[must_use]
pub fn dot_fixed(
    a: &[f64],
    b: &[f64],
    block: usize,
    grain: usize,
) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let block = block.max(1);
    let (a_prefix, _) = a.split_at(n);
    let (b_prefix, _) = b.split_at(n);
    let nb = n.div_ceil(block);
    let minimum = grain.div_ceil(block).max(2);
    map_reduce_blocks(
        nb,
        minimum,
        0.0_f64,
        |b_idx| {
            let lo = b_idx.saturating_mul(block);
            let hi = lo.saturating_add(block).min(n);
            match (a_prefix.get(lo..hi), b_prefix.get(lo..hi)) {
                | (Some(sa), Some(sb)) => dot_block(sa, sb),
                | _ => 0.0,
            }
        },
        |total, value| total + value,
    )
}

/// Deterministic Euclidean (`L2`) norm of `a` using fixed blocks of `block` elements.
///
/// Uses [`dot_fixed`] on the ordinary path and falls back to scale-normalized block
/// accumulation when intermediate squares overflow or underflow, keeping serial
/// (`grain = usize::MAX`) and parallel results bit-identical.
#[must_use]
pub fn norm2_fixed(
    a: &[f64],
    block: usize,
    grain: usize,
) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let block = block.max(1);
    let sum_sq = dot_fixed(a, a, block, grain);
    if sum_sq.is_normal() || sum_sq.is_nan() {
        return sum_sq.sqrt();
    }
    norm2_scaled_fixed(a, block, grain)
}

/// Runs `f(i, &mut items[i])` for every element, in parallel when there are
/// at least `grain` elements. Each element is touched by exactly one call, so
/// the outcome is independent of the scheduling.
pub fn for_each_mut<T, F>(
    items: &mut [T],
    grain: usize,
    f: F,
) where
    T: Send,
    F: Fn(usize, &mut T) + Sync + Send,
{
    if crate::execution::should_split(items.len(), grain) {
        let cancellation = crate::cancellation::current();
        items.par_iter_mut().enumerate().for_each_init(
            || crate::cancellation::install_optional(cancellation.clone()),
            |_, (i, x)| f(i, x),
        );
    } else {
        items.iter_mut().enumerate().for_each(|(i, x)| f(i, x));
    }
}

/// Runs `f(row_start, strip)` over disjoint strips of a row-major buffer.
///
/// The buffer has `rows × cols` elements. Each strip holds up to `block_rows` rows and
/// `row_start` is the zero-based row index of its first row.
///
/// Processes `usable_rows = (data.len() / cols).min(rows)` complete rows sequentially when
/// there are fewer than `grain` strips (or fewer than 2).
pub fn for_blocks_2d_mut<T, F>(
    data: &mut [T],
    rows: usize,
    cols: usize,
    block_rows: usize,
    grain: usize,
    f: F,
) where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync + Send,
{
    if rows == 0 || cols == 0 || data.is_empty() {
        return;
    }
    let Some(available_rows) = data.len().checked_div(cols) else {
        return;
    };
    let usable_rows = available_rows.min(rows);
    if usable_rows == 0 {
        return;
    }
    let used_len = usable_rows.saturating_mul(cols);
    let (active, _) = data.split_at_mut(used_len);
    let block_rows = block_rows.min(usable_rows).max(1);
    let strip_len = block_rows.saturating_mul(cols).max(1);
    let num_strips = usable_rows.div_ceil(block_rows);
    if crate::execution::should_split(num_strips, grain) {
        let cancellation = crate::cancellation::current();
        active.par_chunks_mut(strip_len).enumerate().for_each_init(
            || crate::cancellation::install_optional(cancellation.clone()),
            |_, (strip_idx, strip)| {
                f(strip_idx.saturating_mul(block_rows), strip);
            },
        );
    } else {
        active
            .chunks_mut(strip_len)
            .enumerate()
            .for_each(|(strip_idx, strip)| {
                f(strip_idx.saturating_mul(block_rows), strip);
            });
    }
}

/// Runs `f(z, slab)` over disjoint `nx x ny` slabs along `0..nz` of a row-major
/// `nz x ny x nx` buffer (`data` ordered as `z * (nx * ny) + y * nx + x`).
///
/// Processes `usable_slabs = (data.len() / (nx * ny)).min(nz)` complete slabs sequentially when
/// there are fewer than `grain` slabs (or fewer than 2).
pub fn for_blocks_3d_mut<T, F>(
    data: &mut [T],
    nx: usize,
    ny: usize,
    nz: usize,
    grain: usize,
    f: F,
) where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync + Send,
{
    if nx == 0 || ny == 0 || nz == 0 || data.is_empty() {
        return;
    }
    let Some(slab_len) = nx.checked_mul(ny) else {
        return;
    };
    let Some(available_slabs) = data.len().checked_div(slab_len) else {
        return;
    };
    let usable_slabs = available_slabs.min(nz);
    if usable_slabs == 0 {
        return;
    }
    let used_len = usable_slabs.saturating_mul(slab_len);
    let (active, _) = data.split_at_mut(used_len);
    for_chunks_mut(active, slab_len, grain, f);
}

/// Fallibly runs `f(i, chunk_slice)` over mutable chunks of `data` of `chunk` elements each
/// (`i` is the chunk index), returning the error from the lowest failing chunk index.
///
/// Sequential when there are fewer than `grain` chunks (or fewer than 2). On parallel runs,
/// workers skip chunks beyond an already observed failing chunk index atomically.
///
/// # Errors
/// Returns the error produced by the lowest failing chunk index.
pub fn try_for_chunks_mut<T, E, F>(
    data: &mut [T],
    chunk: usize,
    grain: usize,
    f: F,
) -> Result<(), E>
where
    T: Send,
    E: Send,
    F: Fn(usize, &mut [T]) -> Result<(), E> + Sync + Send,
{
    let chunk = chunk.max(1);
    let num_chunks = data.len().div_ceil(chunk);
    if num_chunks == 0 {
        return Ok(());
    }
    if !crate::execution::should_split(num_chunks, grain) {
        for (index, slice) in data.chunks_mut(chunk).enumerate() {
            f(index, slice)?;
        }
        return Ok(());
    }
    let first_error = AtomicUsize::new(num_chunks);
    let best_error: Mutex<Option<(usize, E)>> = Mutex::new(None);
    let cancellation = crate::cancellation::current();
    data.par_chunks_mut(chunk).enumerate().for_each_init(
        || crate::cancellation::install_optional(cancellation.clone()),
        |_, (index, slice)| {
            if index > first_error.load(Ordering::Relaxed) {
                return;
            }
            let Err(err) = f(index, slice) else {
                return;
            };
            let prev = first_error.fetch_min(index, Ordering::Relaxed);
            if index >= prev {
                return;
            }
            if let Ok(mut guard) = best_error.lock() {
                let replace = guard.as_ref().is_none_or(|(best_idx, _)| index < *best_idx);
                if replace {
                    *guard = Some((index, err));
                }
            }
        },
    );
    match best_error.into_inner() {
        | Ok(Some((_, err))) => Err(err),
        | Err(poisoned) => {
            match poisoned.into_inner() {
                | Some((_, err)) => Err(err),
                | None => Ok(()),
            }
        },
        | Ok(None) => Ok(()),
    }
}

/// Deterministic inclusive prefix sum of `data` using a two-pass fixed-block scan.
///
/// Pass 1 scans each contiguous block of `block` elements independently; a sequential scan
/// computes exclusive block offsets left-to-right; pass 2 adds each block's offset to its
/// elements. Serial (`grain = usize::MAX`) and parallel executions produce bit-identical
/// outputs.
#[must_use]
pub fn prefix_sum_f64_fixed(
    data: &[f64],
    block: usize,
    grain: usize,
) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let block = block.max(1);
    let nb = n.div_ceil(block);
    let minimum = grain.div_ceil(block).max(2);
    let mut out = vec![0.0_f64; n];
    let split = crate::execution::should_split(nb, minimum);
    let scan_block = |(dst, src): (&mut [f64], &[f64])| {
        let mut acc = 0.0_f64;
        for (out_slot, &value) in dst.iter_mut().zip(src) {
            acc += value;
            *out_slot = acc;
        }
        acc
    };
    if nb == 1 {
        scan_block((&mut out, data));
        return out;
    }
    let mut block_sums = super::reuse::f64_scratch(nb);
    if split {
        let cancellation = crate::cancellation::current();
        out.par_chunks_mut(block)
            .zip(data.par_chunks(block))
            .map_init(
                || crate::cancellation::install_optional(cancellation.clone()),
                |_, pair| scan_block(pair),
            )
            .collect_into_vec(&mut block_sums);
    } else {
        block_sums.extend(
            out.chunks_mut(block)
                .zip(data.chunks(block))
                .map(scan_block),
        );
    }
    let mut prefix = 0.0_f64;
    for slot in block_sums.as_mut_slice() {
        let block_total = *slot;
        *slot = prefix;
        prefix += block_total;
    }
    let apply_offset = |(dst, &offset): (&mut [f64], &f64)| {
        if offset != 0.0 {
            for out_slot in dst {
                *out_slot += offset;
            }
        }
    };
    if split {
        let cancellation = crate::cancellation::current();
        out.par_chunks_mut(block)
            .zip(block_sums.par_iter())
            .for_each_init(
                || crate::cancellation::install_optional(cancellation.clone()),
                |_, pair| apply_offset(pair),
            );
    } else {
        out.chunks_mut(block)
            .zip(block_sums.iter())
            .for_each(apply_offset);
    }
    out
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    use super::*;
    use crate::cancellation;

    #[test]
    fn fallible_indexed_map_reports_the_lowest_failing_index() {
        for workers in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                let success = try_map_indexed(256, 2, Ok::<_, usize>);
                assert_eq!(success, Ok((0..256).collect()));
                let failure = try_map_indexed(256, 2, |index| {
                    match index {
                        | 3 | 71 | 233 => Err(index),
                        | _ => Ok(index),
                    }
                });
                assert_eq!(failure, Err(3));
            });
        }
    }

    #[test]
    fn indexed_outputs_reuse_storage_and_keep_scope_and_order() {
        for workers in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                let mut indices = Vec::with_capacity(256);
                let mut flags = Vec::with_capacity(256);
                indices.push(99);
                flags.push(false);
                let indices_buffer = indices.as_ptr();
                let flags_buffer = flags.as_ptr();
                {
                    let _scope = cancellation::install(Arc::new(AtomicBool::new(true)));
                    unzip_indexed_into(113, 2, &mut indices, &mut flags, |i| {
                        (i, cancellation::is_cancelled())
                    });
                    assert_eq!(indices, (0..113).collect::<Vec<_>>());
                    assert_eq!(flags, vec![true; 113]);
                }
                map_indexed_into(71, 2, &mut flags, |_| cancellation::is_cancelled());
                assert_eq!(flags, vec![false; 71]);
                assert_eq!(indices.as_ptr(), indices_buffer);
                assert_eq!(flags.as_ptr(), flags_buffer);
                unzip_indexed_into(0, 0, &mut indices, &mut flags, |_| unreachable!());
                assert!(indices.is_empty() && flags.is_empty());
            });
        }
    }

    #[test]
    fn reduction_blocks_and_ragged_chunk_indices_are_worker_independent() {
        let inputs: Vec<_> = [1e16, 1.0, -1e16, -0.0, 0.125]
            .into_iter()
            .cycle()
            .take(123)
            .collect();
        let expected = inputs
            .chunks(7)
            .map(|chunk| chunk.iter().fold(0.0_f64, |sum, value| sum + value))
            .fold(0.0, |total, sum| total + sum);
        for workers in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                let sum = sum_fixed(inputs.len(), 7, 0, |i| inputs[i]);
                assert_eq!(sum.to_bits(), expected.to_bits());
                let mut chunks = vec![0; 23];
                for_chunks_mut(&mut chunks, 4, 2, |i, chunk| chunk.fill(i));
                assert_eq!(chunks, (0..23).map(|i| i / 4).collect::<Vec<_>>());
                for_chunks_mut(&mut chunks, 0, 0, |i, chunk| chunk.fill(i));
                assert_eq!(chunks, (0..23).collect::<Vec<_>>());
            });
        }
    }

    #[test]
    fn rayon_jobs_restore_and_inherit_the_callers_scope() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
        pool.install(|| {
            {
                let _scope = cancellation::install(Arc::new(AtomicBool::new(true)));
                let flags = map_indexed(100, 2, |_| cancellation::is_cancelled());
                assert!(flags.into_iter().all(|flag| flag));
                assert!(cancellation::is_cancelled());
            }
            let mut flags = vec![true; 100];
            for_each_mut(&mut flags, 2, |_, flag| {
                *flag = cancellation::is_cancelled();
            });
            assert!(flags.into_iter().all(|flag| !flag));
        });
    }

    #[test]
    fn sum_kahan_fixed_is_bit_identical_and_compensated() {
        let edge = [1e16, 1.0, -1e16, 0.25, -0.25, 1e-12];
        assert_eq!(
            sum_kahan_fixed(0, 4, 0, |_| 1.0).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            sum_kahan_fixed(1, 4, 0, |_| 2.5).to_bits(),
            2.5_f64.to_bits()
        );
        assert_eq!(
            sum_kahan_fixed(1, 4, 0, |_| -0.0).to_bits(),
            0.0_f64.to_bits()
        );

        for block in [0usize, 1, 2, 4, 16] {
            let serial = sum_kahan_fixed(edge.len(), block, usize::MAX, |i| {
                edge.get(i).copied().unwrap_or(0.0)
            });
            assert!((serial - (1.0 + 1e-12)).abs() <= 1e-15);
            for workers in [1usize, 2, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
                pool.install(|| {
                    let parallel = sum_kahan_fixed(edge.len(), block, 0, |i| {
                        edge.get(i).copied().unwrap_or(0.0)
                    });
                    assert_eq!(parallel.to_bits(), serial.to_bits());
                });
            }
        }

        let long: Vec<f64> = [1e16, 0.5, -1e16, 1e-9, -0.5, 2e-9]
            .into_iter()
            .cycle()
            .take(257)
            .collect();
        let long_serial = sum_kahan_fixed(long.len(), 7, usize::MAX, |i| {
            long.get(i).copied().unwrap_or(0.0)
        });
        for workers in [1usize, 3, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
            pool.install(|| {
                let par_sum =
                    sum_kahan_fixed(long.len(), 7, 0, |i| long.get(i).copied().unwrap_or(0.0));
                assert_eq!(par_sum.to_bits(), long_serial.to_bits());
            });
        }

        let inf_vals = [1.0, f64::INFINITY, 2.0, 3.0];
        let inf_serial = sum_kahan_fixed(inf_vals.len(), 2, usize::MAX, |i| {
            inf_vals.get(i).copied().unwrap_or(0.0)
        });
        assert!(inf_serial.is_infinite() && inf_serial.is_sign_positive());
        let nan_vals = [f64::INFINITY, -f64::INFINITY, 1.0];
        let nan_serial = sum_kahan_fixed(nan_vals.len(), 2, usize::MAX, |i| {
            nan_vals.get(i).copied().unwrap_or(0.0)
        });
        assert!(nan_serial.is_nan());
    }

    #[test]
    fn dot_fixed_and_norm2_fixed_are_bit_identical_and_scale_safe() {
        assert_eq!(dot_fixed(&[], &[], 8, 0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(norm2_fixed(&[], 8, 0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            dot_fixed(&[2.5], &[4.0], 8, 0).to_bits(),
            10.0_f64.to_bits()
        );
        assert_eq!(norm2_fixed(&[-3.0], 8, 0).to_bits(), 3.0_f64.to_bits());
        assert_eq!(norm2_fixed(&[3.0, 4.0], 8, 0).to_bits(), 5.0_f64.to_bits());
        assert_eq!(
            dot_fixed(&[1.0, 2.0, 3.0], &[4.0, 5.0], 2, 0).to_bits(),
            14.0_f64.to_bits()
        );

        let a: Vec<f64> = [1e8, 1.0, -1e8, 0.125, -0.375, 2.25, 0.5]
            .into_iter()
            .cycle()
            .take(131)
            .collect();
        let b: Vec<f64> = [1e8, 0.5, 1e8, -2.0, 1.5, 0.25, -1.0]
            .into_iter()
            .cycle()
            .take(131)
            .collect();

        for block in [0usize, 1, 3, 5, 16, 64] {
            let dot_serial = dot_fixed(&a, &b, block, usize::MAX);
            let norm_serial = norm2_fixed(&a, block, usize::MAX);
            for workers in [1usize, 2, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
                pool.install(|| {
                    assert_eq!(dot_fixed(&a, &b, block, 0).to_bits(), dot_serial.to_bits());
                    assert_eq!(norm2_fixed(&a, block, 0).to_bits(), norm_serial.to_bits());
                });
            }
        }

        let huge = vec![1e200_f64; 37];
        let huge_serial = norm2_fixed(&huge, 5, usize::MAX);
        assert!(huge_serial.is_finite());
        assert!((huge_serial - 1e200 * (37.0_f64).sqrt()).abs() <= 1e186);
        let tiny = vec![1e-200_f64; 37];
        let tiny_serial = norm2_fixed(&tiny, 5, usize::MAX);
        assert!(tiny_serial > 0.0 && tiny_serial.is_finite());
        assert!((tiny_serial - 1e-200 * (37.0_f64).sqrt()).abs() <= 1e-214);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
        pool.install(|| {
            assert_eq!(norm2_fixed(&huge, 5, 0).to_bits(), huge_serial.to_bits());
            assert_eq!(norm2_fixed(&tiny, 5, 0).to_bits(), tiny_serial.to_bits());
        });

        assert_eq!(
            norm2_fixed(&[0.0, -0.0, 0.0], 2, 0).to_bits(),
            0.0_f64.to_bits()
        );
        assert!(norm2_fixed(&[f64::INFINITY, 1.0], 2, 0).is_infinite());
        assert!(norm2_fixed(&[f64::NAN, 1.0], 2, 0).is_nan());
    }

    #[test]
    fn for_blocks_2d_mut_and_3d_mut_cover_every_element_and_match_serial() {
        let mut empty = Vec::<usize>::new();
        for_blocks_2d_mut(&mut empty, 0, 4, 2, 0, |_, _| unreachable!());
        for_blocks_2d_mut(&mut empty, 4, 0, 2, 0, |_, _| unreachable!());
        for_blocks_3d_mut(&mut empty, 0, 2, 2, 0, |_, _| unreachable!());
        for_blocks_3d_mut(&mut empty, 2, 0, 2, 0, |_, _| unreachable!());
        for_blocks_3d_mut(&mut empty, 2, 2, 0, 0, |_, _| unreachable!());

        let mut one = vec![0usize; 1];
        for_blocks_2d_mut(&mut one, 1, 1, 0, 0, |row_start, strip| {
            assert_eq!(row_start, 0);
            assert_eq!(strip.len(), 1);
            strip.fill(7);
        });
        assert_eq!(one, vec![7]);
        for_blocks_3d_mut(&mut one, 1, 1, 1, 0, |z, slab| {
            assert_eq!(z, 0);
            assert_eq!(slab.len(), 1);
            slab.fill(9);
        });
        assert_eq!(one, vec![9]);

        let (rows, cols) = (17usize, 5usize);
        let fill_2d = |row_start: usize, strip: &mut [usize]| {
            for (local_r, row) in strip.chunks_exact_mut(cols).enumerate() {
                let r = row_start.saturating_add(local_r);
                for (c, cell) in row.iter_mut().enumerate() {
                    *cell = r.saturating_mul(100).saturating_add(c);
                }
            }
        };
        for block_rows in [0usize, 1, 3, 17, 100] {
            let mut serial_2d = vec![usize::MAX; rows.saturating_mul(cols)];
            for_blocks_2d_mut(&mut serial_2d, rows, cols, block_rows, usize::MAX, fill_2d);
            for workers in [1usize, 2, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
                pool.install(|| {
                    let mut par_2d = vec![usize::MAX; rows.saturating_mul(cols)];
                    for_blocks_2d_mut(&mut par_2d, rows, cols, block_rows, 0, fill_2d);
                    assert_eq!(par_2d, serial_2d);
                });
            }
        }

        let (nx, ny, nz) = (3usize, 4usize, 11usize);
        let fill_3d = |z: usize, slab: &mut [usize]| {
            assert_eq!(slab.len(), nx.saturating_mul(ny));
            for (offset, cell) in slab.iter_mut().enumerate() {
                *cell = z.saturating_mul(1_000).saturating_add(offset);
            }
        };
        let mut serial_3d = vec![usize::MAX; nx.saturating_mul(ny).saturating_mul(nz)];
        for_blocks_3d_mut(&mut serial_3d, nx, ny, nz, usize::MAX, fill_3d);
        for workers in [1usize, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
            pool.install(|| {
                let mut par_3d = vec![usize::MAX; nx.saturating_mul(ny).saturating_mul(nz)];
                for_blocks_3d_mut(&mut par_3d, nx, ny, nz, 0, fill_3d);
                assert_eq!(par_3d, serial_3d);
            });
        }
    }

    #[test]
    fn try_for_chunks_mut_reports_first_error_deterministically_and_skips_later() {
        let mut empty = Vec::<usize>::new();
        assert_eq!(try_for_chunks_mut(&mut empty, 4, 0, |_, _| Err(())), Ok(()));

        let mut single = vec![0usize];
        assert_eq!(
            try_for_chunks_mut(&mut single, 0, 0, |i, slice| {
                slice.fill(i.saturating_add(1));
                Ok::<_, usize>(())
            }),
            Ok(())
        );
        assert_eq!(single, vec![1]);
        assert_eq!(
            try_for_chunks_mut(&mut single, 4, 0, |i, _| Err::<(), usize>(i)),
            Err(0)
        );

        let ran_serial = AtomicUsize::new(0);
        let mut seq_buf = vec![0usize; 256];
        let seq_res = try_for_chunks_mut(&mut seq_buf, 4, usize::MAX, |index, slice| {
            ran_serial.fetch_add(1, Ordering::Relaxed);
            if matches!(index, 7 | 29 | 51) {
                return Err(index);
            }
            slice.fill(index);
            Ok(())
        });
        assert_eq!(seq_res, Err(7));
        assert_eq!(ran_serial.load(Ordering::Relaxed), 8);

        for workers in [1usize, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
            pool.install(|| {
                let mut ok_buf = vec![0usize; 255];
                let ok_res = try_for_chunks_mut(&mut ok_buf, 4, 2, |index, slice| {
                    slice.fill(index);
                    Ok::<_, usize>(())
                });
                assert_eq!(ok_res, Ok(()));
                assert_eq!(
                    ok_buf,
                    (0_usize..255)
                        .map(|p| p.checked_div(4).unwrap_or(0))
                        .collect::<Vec<_>>()
                );

                let mut err_buf = vec![0usize; 512];
                let err_res = try_for_chunks_mut(&mut err_buf, 4, 2, |index, slice| {
                    if matches!(index, 3 | 71 | 117) {
                        return Err(index);
                    }
                    slice.fill(index);
                    Ok(())
                });
                assert_eq!(err_res, Err(3));
            });
        }
    }

    #[test]
    fn prefix_sum_f64_fixed_is_bit_identical_and_matches_block_scan() {
        assert_eq!(prefix_sum_f64_fixed(&[], 4, 0), [] as [f64; 0]);
        assert_eq!(prefix_sum_f64_fixed(&[3.5], 4, 0), vec![3.5]);
        assert_eq!(
            prefix_sum_f64_fixed(&[-0.0], 4, 0)
                .first()
                .map(|x| x.to_bits()),
            Some(0.0_f64.to_bits())
        );
        assert_eq!(
            prefix_sum_f64_fixed(&[1.0, 2.0, 3.0, 4.0, 5.0], 2, 0),
            vec![1.0, 3.0, 6.0, 10.0, 15.0]
        );

        let inputs: Vec<f64> = [1e16, 1.0, -1e16, -0.0, 0.125]
            .into_iter()
            .cycle()
            .take(123)
            .collect();
        let mut naive_seq = Vec::with_capacity(inputs.len());
        let mut acc = 0.0_f64;
        for &val in &inputs {
            acc += val;
            naive_seq.push(acc);
        }
        let b1 = prefix_sum_f64_fixed(&inputs, 1, usize::MAX);
        assert_eq!(
            b1.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            naive_seq.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );

        for block in [0usize, 1, 7, 64, 1_000] {
            let serial = prefix_sum_f64_fixed(&inputs, block, usize::MAX);
            let serial_bits: Vec<u64> = serial.iter().map(|v| v.to_bits()).collect();
            for workers in [1usize, 2, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
                pool.install(|| {
                    let par_out = prefix_sum_f64_fixed(&inputs, block, 0);
                    let par_bits: Vec<u64> = par_out.iter().map(|v| v.to_bits()).collect();
                    assert_eq!(par_bits, serial_bits);
                });
            }
        }
    }

    #[test]
    fn new_parallel_primitives_inherit_and_restore_cancellation_scope() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
        pool.install(|| {
            {
                let _scope = cancellation::install(Arc::new(AtomicBool::new(true)));
                let kahan_sum = sum_kahan_fixed(64, 4, 2, |_| {
                    if cancellation::is_cancelled() {
                        1.0
                    } else {
                        0.0
                    }
                });
                assert_eq!(kahan_sum.to_bits(), 64.0_f64.to_bits());

                let mut grid2d = vec![false; 32];
                for_blocks_2d_mut(&mut grid2d, 16, 2, 2, 2, |_, strip| {
                    strip.fill(cancellation::is_cancelled());
                });
                assert!(grid2d.into_iter().all(|flag| flag));

                let mut grid3d = vec![false; 48];
                for_blocks_3d_mut(&mut grid3d, 2, 2, 12, 2, |_, slab| {
                    slab.fill(cancellation::is_cancelled());
                });
                assert!(grid3d.into_iter().all(|flag| flag));

                let mut fallible = vec![0usize; 64];
                let res =
                    try_for_chunks_mut(&mut fallible, 4, 2, |_, _| cancellation::checkpoint());
                assert_eq!(res, Err(cancellation::Cancelled));
            }

            let mut clean = vec![true; 32];
            for_blocks_2d_mut(&mut clean, 16, 2, 2, 2, |_, strip| {
                strip.fill(cancellation::is_cancelled());
            });
            assert!(clean.into_iter().all(|flag| !flag));
        });

        crate::execution::with_scoped_pool(Some(1), || {
            let scoped_sum = sum_kahan_fixed(256, 8, 0, |i| i as f64);
            assert_eq!(scoped_sum.to_bits(), (32_640.0_f64).to_bits());
        });
    }
}
