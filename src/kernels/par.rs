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
    use std::sync::atomic::AtomicUsize;
    use std::sync::atomic::Ordering;

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
    let part = |b: usize| {
        let lo = b.saturating_mul(block);
        let hi = lo.saturating_add(block).min(n);
        let mut s = 0.0;
        for i in lo..hi {
            s += f(i);
        }
        s
    };
    let minimum = grain.div_ceil(block).max(2);
    if crate::execution::should_split(nb, minimum) {
        indexed_map(nb, part)
            .collect::<Vec<_>>()
            .into_iter()
            .fold(0.0, |total, value| total + value)
    } else {
        (0..nb).map(part).fold(0.0, |total, value| total + value)
    }
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
}
