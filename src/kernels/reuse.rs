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
use std::sync::TryLockError;

use num_complex::Complex64;
use once_cell::sync::OnceCell;
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

type PlanSlot<V> = Arc<OnceCell<Arc<V>>>;

enum PlanAccess<V: ?Sized> {
    Ready(Arc<V>),
    Pending(PlanSlot<V>),
    Uncached,
}

struct PlanEntry<K, V: ?Sized> {
    key: K,
    slot: PlanSlot<V>,
    last_access: u64,
}

struct PlanEntries<K, V: ?Sized> {
    queue: VecDeque<PlanEntry<K, V>>,
    clock: u64,
}

/// Small bounded cache of immutable plans.
pub struct PlanCache<K, V: ?Sized> {
    entries: Mutex<PlanEntries<K, V>>,
    capacity: usize,
    promote_hits: bool,
}

impl<K: Eq, V: ?Sized> PlanCache<K, V> {
    /// Creates an empty cache that evicts in insertion order.
    /// A capacity of zero builds each plan on demand.
    #[must_use]
    pub const fn new(capacity: usize) -> Self {
        Self::with_policy(capacity, false)
    }

    const fn with_lru(capacity: usize) -> Self {
        Self::with_policy(capacity, true)
    }

    const fn with_policy(
        capacity: usize,
        promote_hits: bool,
    ) -> Self {
        Self {
            entries: Mutex::new(PlanEntries {
                queue: VecDeque::new(),
                clock: 0,
            }),
            capacity,
            promote_hits,
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
        match self.lookup(key) {
            | PlanAccess::Ready(plan) => plan,
            | PlanAccess::Pending(slot) => Arc::clone(slot.get_or_init(build)),
            | PlanAccess::Uncached => build(),
        }
    }

    /// Shares construction of a fallible plan. A failed build leaves the slot
    /// empty so callers can retry, including after a transient allocation error.
    pub(crate) fn try_get_or_insert_with<E>(
        &self,
        key: K,
        build: impl FnOnce() -> Result<Arc<V>, E>,
    ) -> Result<Arc<V>, E> {
        match self.lookup(key) {
            | PlanAccess::Ready(plan) => Ok(plan),
            | PlanAccess::Pending(slot) => slot.get_or_try_init(build).map(Arc::clone),
            | PlanAccess::Uncached => build(),
        }
    }

    fn lookup(
        &self,
        key: K,
    ) -> PlanAccess<V> {
        if self.capacity == 0 {
            return PlanAccess::Uncached;
        }
        let mut entries = self
            .entries
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if self.promote_hits && entries.clock == u64::MAX {
            // Preserve order when shrinking the timestamps back into range.
            let queue = entries.queue.make_contiguous();
            queue.sort_unstable_by_key(|entry| entry.last_access);
            for (index, entry) in queue.iter_mut().enumerate() {
                entry.last_access = index as u64;
            }
            entries.clock = queue.len() as u64;
        }
        let stamp = entries.clock;
        entries.clock += u64::from(self.promote_hits);
        if let Some(entry) = entries.queue.iter_mut().find(|entry| entry.key == key) {
            entry.last_access = stamp;
            return entry.slot.get().map_or_else(
                || PlanAccess::Pending(Arc::clone(&entry.slot)),
                |plan| PlanAccess::Ready(Arc::clone(plan)),
            );
        }
        let evicted = if entries.queue.len() == self.capacity {
            // A client may still own a completed plan. A pending build stays
            // visible to concurrent requests until it completes or unwinds.
            let can_evict = |entry: &PlanEntry<K, V>| {
                entry.slot.get().is_some() || Arc::strong_count(&entry.slot) == 1
            };
            let oldest = if self.promote_hits {
                entries
                    .queue
                    .iter()
                    .enumerate()
                    .filter(|(_, entry)| can_evict(entry))
                    .min_by_key(|(_, entry)| entry.last_access)
                    .map(|(index, _)| index)
            } else {
                entries.queue.iter().position(can_evict)
            };
            let Some(oldest) = oldest else {
                return PlanAccess::Uncached;
            };
            entries.queue.remove(oldest)
        } else {
            None
        };
        let slot = Arc::new(OnceCell::new());
        entries.queue.push_back(PlanEntry {
            key,
            slot: Arc::clone(&slot),
            last_access: stamp,
        });
        drop(entries);
        drop(evicted);
        PlanAccess::Pending(slot)
    }
}

/// Plans an `f64` complex FFT using `RustFFT`'s runtime CPU feature detection.
///
/// Small FFTs share their plans across algorithms and threads. A separate
/// two-entry cache admits larger transforms with bounded input and scratch
/// lengths so repeated inverses and convolutions can reuse their plans.
/// The FFT caches evict the least recently used available entry.
#[must_use]
pub fn fft_plan(
    length: usize,
    direction: FftDirection,
) -> Arc<dyn Fft<f64>> {
    static PLANS: PlanCache<(usize, FftDirection), dyn Fft<f64>> = PlanCache::with_lru(16);
    static LARGE_PLANS: PlanCache<(usize, FftDirection), dyn Fft<f64>> = PlanCache::with_lru(2);
    const LARGE_PLAN_MAX_POINTS: usize = 524_288;
    const LARGE_PLAN_MAX_INPUT_AND_SCRATCH: usize = 1_048_576;

    let build = || FftPlanner::new().plan_fft(length, direction);
    if length <= 65_536 {
        PLANS.get_or_insert_with((length, direction), build)
    } else if length <= LARGE_PLAN_MAX_POINTS {
        LARGE_PLANS
            .try_get_or_insert_with((length, direction), || {
                let plan = build();
                if plan.get_inplace_scratch_len().saturating_add(length)
                    <= LARGE_PLAN_MAX_INPUT_AND_SCRATCH
                {
                    Ok(plan)
                } else {
                    Err(plan)
                }
            })
            .unwrap_or_else(std::convert::identity)
    } else {
        build()
    }
}

/// Checks out a resizable scratch vector for complex `f64` transforms.
#[must_use]
pub fn fft_scratch() -> Scratch<'static, Complex64> {
    COMPLEX64_SCRATCH.checkout()
}

const COMPLEX64_SCRATCH_BYTES: usize = 4 * 1024 * 1024;
static COMPLEX64_SCRATCH: ScratchPool<Complex64> = ScratchPool::new(8, COMPLEX64_SCRATCH_BYTES);
static LARGE_COMPLEX64_SCRATCH: ScratchPool<Complex64> = ScratchPool::new(2, 16 * 1024 * 1024);

/// Checks out an empty complex `f64` buffer with the requested capacity.
///
/// Buffers over four MiB use their own retention pool, which retains up to
/// two idle buffers of at most sixteen MiB each.
///
/// # Panics
/// Panics if `minimum` exceeds the maximum capacity of a vector.
#[must_use]
pub fn fft_scratch_for(minimum: usize) -> Scratch<'static, Complex64> {
    if minimum > COMPLEX64_SCRATCH_BYTES / mem::size_of::<Complex64>() {
        LARGE_COMPLEX64_SCRATCH.checkout_for(minimum)
    } else {
        COMPLEX64_SCRATCH.checkout_for(minimum)
    }
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
    fn fft_recency_keeps_an_interleaved_hot_key() {
        let plans = PlanCache::with_lru(2);
        let hot = plans.get_or_insert_with(0, || Arc::new(41));
        for cold in 1..=20 {
            assert_eq!(*plans.get_or_insert_with(cold, || Arc::new(cold)), cold);
            let cached = plans.get_or_insert_with(0, || panic!("hot FFT plan was evicted"));
            assert!(Arc::ptr_eq(&hot, &cached));
        }
    }

    #[test]
    fn fft_recency_preserves_eviction_order_when_the_clock_wraps() {
        let plans = PlanCache::with_lru(2);
        let hot = plans.get_or_insert_with(0, || Arc::new(41));
        let _cold = plans.get_or_insert_with(1, || Arc::new(1));
        plans.entries.lock().unwrap().clock = u64::MAX - 1;
        let _hot = plans.get_or_insert_with(0, || panic!("hot plan was absent before wrap"));
        let _new = plans.get_or_insert_with(2, || Arc::new(2));
        let cached = plans.get_or_insert_with(0, || panic!("hot plan was evicted at wrap"));
        assert!(Arc::ptr_eq(&hot, &cached));
    }

    #[test]
    fn insertion_order_does_not_change_after_a_hit() {
        let plans = PlanCache::new(2);
        let _oldest = plans.get_or_insert_with(0, || Arc::new(40));
        let _second = plans.get_or_insert_with(1, || Arc::new(1));
        let _cached = plans.get_or_insert_with(0, || panic!("first plan was absent"));
        let _new = plans.get_or_insert_with(2, || Arc::new(2));
        assert_eq!(*plans.get_or_insert_with(0, || Arc::new(42)), 42);
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
    fn a_fallible_plan_retries_and_shares_the_successful_value() {
        let cache = PlanCache::new(1);
        let failed = cache.try_get_or_insert_with(1, || Err("no work buffer"));
        assert_eq!(failed, Err("no work buffer"));
        let recovered = cache
            .try_get_or_insert_with(1, || Ok::<_, &str>(Arc::new(7)))
            .unwrap();
        let infallible = cache.get_or_insert_with(1, || panic!("already initialized"));
        assert!(Arc::ptr_eq(&recovered, &infallible));
        let fallible = cache
            .try_get_or_insert_with(1, || -> Result<Arc<_>, &str> {
                panic!("already initialized")
            })
            .unwrap();
        assert!(Arc::ptr_eq(&recovered, &fallible));
    }

    #[test]
    fn an_unrelated_miss_does_not_wait_for_the_only_pending_slot() {
        use std::sync::mpsc::channel;
        use std::time::Duration;

        let cache = PlanCache::new(1);
        let (started, wait_for_start) = channel();
        let (release, wait_for_release) = channel();
        let (done, wait_for_miss) = channel();
        std::thread::scope(|scope| {
            let first_cache = &cache;
            let pending = scope.spawn(move || {
                first_cache.get_or_insert_with(1, || {
                    started.send(()).unwrap();
                    wait_for_release
                        .recv_timeout(Duration::from_secs(5))
                        .unwrap();
                    Arc::new(11)
                })
            });
            wait_for_start.recv_timeout(Duration::from_secs(5)).unwrap();
            let second_cache = &cache;
            scope.spawn(move || {
                let plan = second_cache.get_or_insert_with(2, || Arc::new(22));
                done.send(plan).unwrap();
            });
            assert_eq!(
                *wait_for_miss.recv_timeout(Duration::from_secs(2)).unwrap(),
                22
            );
            release.send(()).unwrap();
            let first = pending.join().unwrap();
            assert!(Arc::ptr_eq(
                &first,
                &cache.get_or_insert_with(1, || panic!("pending slot evicted"))
            ));
            assert_eq!(*cache.get_or_insert_with(2, || Arc::new(33)), 33);
        });
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

    #[test]
    fn large_fft_reuses_plans_and_keeps_large_scratch_separate() {
        let length = 131_072;
        let forward = fft_plan(length, FftDirection::Forward);
        let inverse = fft_plan(length, FftDirection::Inverse);
        assert!(Arc::ptr_eq(
            &forward,
            &fft_plan(length, FftDirection::Forward)
        ));
        assert!(Arc::ptr_eq(
            &inverse,
            &fft_plan(length, FftDirection::Inverse)
        ));

        let points = COMPLEX64_SCRATCH_BYTES / mem::size_of::<Complex64>() + 1;
        let scratch = fft_scratch_for(points);
        assert!(scratch.capacity() >= points);
        assert!(std::ptr::eq(scratch.pool, &LARGE_COMPLEX64_SCRATCH));
        let small = fft_scratch_for(17);
        assert!(std::ptr::eq(small.pool, &COMPLEX64_SCRATCH));
    }
}
