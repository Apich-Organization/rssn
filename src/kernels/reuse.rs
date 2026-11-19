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
use std::sync::MutexGuard;
use std::sync::OnceLock;
use std::sync::TryLockError;

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
    /// Temporary contention is handled by allocating on demand.
    #[must_use]
    pub fn checkout(&self) -> Scratch<'_, T> {
        Scratch {
            pool: self,
            buffer: self
                .available_idle()
                .and_then(|mut idle| idle.pop())
                .unwrap_or_default(),
        }
    }

    /// Takes an empty vector with room for at least `minimum` elements.
    /// Prefers the smallest retained allocation that fits. Requests larger
    /// than the retention limit leave the idle allocations available to other jobs.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn checkout_for(
        &self,
        minimum: usize,
    ) -> Scratch<'_, T> {
        let buffer = if minimum == 0
            || self.max_idle == 0
            || mem::size_of::<T>() == 0
            || minimum.saturating_mul(mem::size_of::<T>()) > self.max_buffer_bytes
        {
            Vec::with_capacity(minimum)
        } else {
            let retained = self.available_idle().and_then(|mut idle| {
                let best = idle
                    .iter()
                    .enumerate()
                    .filter(|(_, buffer)| buffer.capacity() >= minimum)
                    .min_by_key(|(_, buffer)| buffer.capacity())
                    .map(|(index, _)| index);
                best.map(|index| idle.swap_remove(index))
            });
            retained.unwrap_or_else(|| Vec::with_capacity(minimum))
        };
        Scratch { pool: self, buffer }
    }

    fn available_idle(&self) -> Option<MutexGuard<'_, Vec<Vec<T>>>> {
        match self.idle.try_lock() {
            | Ok(idle) => Some(idle),
            | Err(TryLockError::Poisoned(error)) => Some(error.into_inner()),
            | Err(TryLockError::WouldBlock) => None,
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
            || self.pool.max_idle == 0
            || self.buffer.capacity() == 0
            || mem::size_of::<T>() == 0
            || self.buffer.capacity().saturating_mul(mem::size_of::<T>())
                > self.pool.max_buffer_bytes
        {
            return;
        }
        // Drop elements outside the lock. Only the allocation is reused.
        self.buffer.clear();
        let Some(mut idle) = self.pool.available_idle() else {
            return;
        };
        if idle.len() < self.pool.max_idle {
            idle.push(mem::take(&mut self.buffer));
        } else if let Some(smallest) = idle.iter_mut().min_by_key(|buffer| buffer.capacity())
            && smallest.capacity() < self.buffer.capacity()
        {
            // Exchange allocations and let the smaller one drop outside the lock.
            mem::swap(smallest, &mut self.buffer);
        }
    }
}

type PlanSlot<V> = Arc<OnceLock<Arc<V>>>;

/// Small bounded cache of immutable plans, evicted in insertion order.
pub struct PlanCache<K, V: ?Sized> {
    entries: Mutex<VecDeque<(K, PlanSlot<V>)>>,
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
    /// Callers requesting a pending key share its construction. If every
    /// slot holds a pending plan, a miss for a different key builds without
    /// retaining the result. A panicking builder leaves its slot retryable.
    /// The key must include every input that changes the plan. Builders must
    /// avoid cyclic dependencies on other keys in the same cache.
    pub fn get_or_insert_with(
        &self,
        key: K,
        build: impl FnOnce() -> Arc<V>,
    ) -> Arc<V> {
        if self.capacity == 0 {
            return build();
        }
        let mut entries = self
            .entries
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, slot)) = entries.iter().find(|(candidate, _)| *candidate == key) {
            let slot = Arc::clone(slot);
            drop(entries);
            return Arc::clone(slot.get_or_init(build));
        }
        let evicted = if entries.len() == self.capacity {
            // The cache is the sole owner of a failed slot once its callers
            // unwind. It can also evict an initialized slot with active readers.
            let Some(oldest) = entries
                .iter()
                .position(|(_, slot)| slot.get().is_some() || Arc::strong_count(slot) == 1)
            else {
                drop(entries);
                return build();
            };
            entries.remove(oldest)
        } else {
            None
        };
        let slot = Arc::new(OnceLock::new());
        entries.push_back((key, Arc::clone(&slot)));
        drop(entries);
        let plan = Arc::clone(slot.get_or_init(build));
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
    COMPLEX64_SCRATCH.checkout()
}

static COMPLEX64_SCRATCH: ScratchPool<Complex64> = ScratchPool::new(8, 4 * 1024 * 1024);

/// Checks out an empty complex `f64` buffer with the requested capacity.
///
/// # Panics
/// Panics if `minimum` exceeds the maximum capacity of a vector.
#[must_use]
pub fn fft_scratch_for(minimum: usize) -> Scratch<'static, Complex64> {
    COMPLEX64_SCRATCH.checkout_for(minimum)
}

/// Checks out an empty `f64` buffer with the requested capacity. Kernels
/// share this pool for intermediate values and return only their own results.
///
/// # Panics
/// Panics if `minimum` exceeds the maximum capacity of a vector.
#[must_use]
pub fn f64_scratch(minimum: usize) -> Scratch<'static, f64> {
    static SCRATCH: ScratchPool<f64> = ScratchPool::new(8, 4 * 1024 * 1024);
    SCRATCH.checkout_for(minimum)
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
    fn sized_scratch_keeps_allocation_for_matching_requests() {
        let pool = ScratchPool::<u64>::new(2, 128);
        let small = pool.checkout_for(4);
        let large = pool.checkout_for(16);
        let (small_ptr, large_ptr) = (small.as_ptr(), large.as_ptr());
        drop(large);
        drop(small);
        let oversize = pool.checkout_for(32);
        assert!(oversize.capacity() >= 32);
        let medium = pool.checkout_for(8);
        assert_eq!(medium.as_ptr(), large_ptr);
        assert_eq!(pool.checkout_for(4).as_ptr(), small_ptr);
    }

    #[test]
    fn concurrent_requests_share_one_build() {
        use std::sync::Barrier;
        use std::sync::atomic::AtomicUsize;
        use std::sync::atomic::Ordering;

        let cache = PlanCache::new(1);
        let ready = Barrier::new(8);
        let builds = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    ready.wait();
                    let plan = cache.get_or_insert_with(1, || {
                        builds.fetch_add(1, Ordering::Relaxed);
                        Arc::new(41)
                    });
                    assert_eq!(*plan, 41);
                });
            }
        });
        assert_eq!(builds.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn a_panicking_plan_builder_is_retryable() {
        let cache = PlanCache::new(1);
        let failed = std::panic::catch_unwind(|| {
            cache.get_or_insert_with(1, || panic!("plan builder failed"))
        });
        assert!(failed.is_err());
        assert_eq!(*cache.get_or_insert_with(1, || Arc::new(42)), 42);
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
