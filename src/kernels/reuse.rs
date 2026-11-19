//! Bounded reuse of scratch allocations and immutable numerical plans.
//!
//! A scratch lease exclusively owns its buffer while a plan can be shared by
//! concurrent callers. Retention limits apply to idle buffers and cache
//! entries. Active leases and caller-held `Arc`s can outlive cache entries.

use std::collections::VecDeque;
use std::mem;
use std::ops::Deref;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;

use num_complex::Complex64;
use rustfft::Fft;
use rustfft::FftDirection;
use rustfft::FftPlanner;

/// Pool of empty vectors that retain their allocation between computations.
pub struct ScratchPool<T> {
    idle: Mutex<Vec<Vec<T>>>,
    max_idle: usize,
    max_buffer_bytes: usize,
}

impl<T> ScratchPool<T> {
    /// Retains up to `max_idle` buffers, each at most `max_buffer_bytes` bytes.
    /// A zero limit disables retention.
    #[must_use]
    pub const fn new(
        max_idle: usize,
        max_buffer_bytes: usize,
    ) -> Self {
        Self {
            idle: Mutex::new(Vec::new()),
            max_idle,
            max_buffer_bytes,
        }
    }

    /// Takes an empty vector. Resize or reserve it for the current job.
    #[must_use]
    pub fn checkout(&self) -> Scratch<'_, T> {
        Scratch {
            pool: self,
            buffer: self
                .idle
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .pop()
                .unwrap_or_default(),
        }
    }
}

/// Exclusive scratch buffer returned to its pool when dropped.
pub struct Scratch<'a, T> {
    pool: &'a ScratchPool<T>,
    buffer: Vec<T>,
}

impl<T> Deref for Scratch<'_, T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.buffer
    }
}

impl<T> DerefMut for Scratch<'_, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.buffer
    }
}

impl<T> Drop for Scratch<'_, T> {
    fn drop(&mut self) {
        if self.pool.max_buffer_bytes == 0
            || self.buffer.capacity() == 0
            || self.buffer.capacity().saturating_mul(mem::size_of::<T>())
                > self.pool.max_buffer_bytes
        {
            return;
        }
        // Drop elements outside the lock. Only the allocation is reused.
        self.buffer.clear();
        let mut idle = self
            .pool
            .idle
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if idle.len() < self.pool.max_idle {
            idle.push(mem::take(&mut self.buffer));
        }
    }
}

/// Small bounded cache of immutable plans, evicted in insertion order.
pub struct PlanCache<K, V: ?Sized> {
    entries: Mutex<VecDeque<(K, Arc<V>)>>,
    capacity: usize,
}

impl<K: Eq, V: ?Sized> PlanCache<K, V> {
    /// Creates an empty cache. A capacity of zero builds each plan on demand.
    #[must_use]
    pub const fn new(capacity: usize) -> Self {
        Self {
            entries: Mutex::new(VecDeque::new()),
            capacity,
        }
    }

    /// Reuses a plan for `key`, or builds one outside the cache lock.
    /// Concurrent misses can build the same plan. The first plan published
    /// for the key is shared while its entry is present. The key must include
    /// every input that changes the plan.
    pub fn get_or_insert_with(
        &self,
        key: K,
        build: impl FnOnce() -> Arc<V>,
    ) -> Arc<V> {
        if self.capacity == 0 {
            return build();
        }
        let cached = self
            .entries
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .iter()
            .find(|(candidate, _)| *candidate == key)
            .map(|(_, value)| Arc::clone(value));
        if let Some(hit) = cached {
            return hit;
        }
        let plan = build();
        let mut entries = self
            .entries
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, hit)) = entries.iter().find(|(candidate, _)| *candidate == key) {
            return Arc::clone(hit);
        }
        let evicted = if entries.len() == self.capacity {
            entries.pop_front()
        } else {
            None
        };
        entries.push_back((key, Arc::clone(&plan)));
        drop(entries);
        drop(evicted);
        plan
    }
}

/// Plans an `f64` complex FFT using `RustFFT`'s runtime CPU feature detection.
/// Moderately sized FFTs share their plans across algorithms and threads.
#[must_use]
pub fn fft_plan(
    length: usize,
    direction: FftDirection,
) -> Arc<dyn Fft<f64>> {
    static PLANS: PlanCache<(usize, FftDirection), dyn Fft<f64>> = PlanCache::new(16);
    let build = || FftPlanner::new().plan_fft(length, direction);
    if length <= 65_536 {
        PLANS.get_or_insert_with((length, direction), build)
    } else {
        build()
    }
}

/// Checks out a resizable scratch vector for complex `f64` transforms.
#[must_use]
pub fn fft_scratch() -> Scratch<'static, Complex64> {
    static SCRATCH: ScratchPool<Complex64> = ScratchPool::new(8, 4 * 1024 * 1024);
    SCRATCH.checkout()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scratch_is_exclusive_and_returned_empty_under_its_limits() {
        let pool = ScratchPool::new(1, 32);
        let mut first = pool.checkout();
        first.resize(4, 42_u64);
        let ptr = first.as_ptr();
        let mut second = pool.checkout();
        second.resize(16, 17);
        assert_ne!(second.as_ptr(), ptr);
        drop(first);
        drop(second);
        let reused = pool.checkout();
        assert!(reused.is_empty());
        assert_eq!(reused.as_ptr(), ptr);
        assert_eq!(pool.checkout().capacity(), 0);
    }

    #[test]
    fn plans_are_shared_and_eviction_does_not_invalidate_handles() {
        let plans = PlanCache::new(1);
        let a = plans.get_or_insert_with(1, || Arc::new(10));
        let cached = plans.get_or_insert_with(1, || panic!("a cache hit must not build"));
        assert!(Arc::ptr_eq(&a, &cached));
        assert_eq!(*plans.get_or_insert_with(2, || Arc::new(20)), 20);
        assert_eq!(*plans.get_or_insert_with(1, || Arc::new(30)), 30);
        assert_eq!(*a, 10);
    }

    #[test]
    fn fft_plans_are_keyed_by_length_and_direction() {
        let forward = fft_plan(17, FftDirection::Forward);
        let inverse = fft_plan(17, FftDirection::Inverse);
        assert_eq!(forward.len(), 17);
        assert_eq!(forward.fft_direction(), FftDirection::Forward);
        assert_eq!(inverse.fft_direction(), FftDirection::Inverse);
        assert_eq!(fft_plan(18, FftDirection::Forward).len(), 18);
    }
}
