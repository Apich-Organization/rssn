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
        let cancellation = crate::cancellation::current();
        (0..n)
            .into_par_iter()
            .map_init(
                || crate::cancellation::install_optional(cancellation.clone()),
                |_, index| f(index),
            )
            .collect()
    } else {
        (0..n).map(f).collect()
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
    let parts = map_indexed(nb, grain.div_ceil(block).max(2), part);
    let mut s = 0.0;
    for p in parts {
        s += p;
    }
    s
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
                *flag = cancellation::is_cancelled()
            });
            assert!(flags.into_iter().all(|flag| !flag));
        });
    }
}
