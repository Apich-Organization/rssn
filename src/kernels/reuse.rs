//! Unified compute reuse pool for scratch allocations, immutable plans, and
//! reusable numerical tables.
//!
//! A [`Scratch`] lease exclusively owns its buffer while cached plans and
//! tables are shared by concurrent callers via [`Arc`]. Retention limits apply
//! to idle buffers and cache entries. Active leases and caller-held `Arc`s can
//! outlive cache entries, evictions, and [`ComputeReusePool::clear_all`].
//!
//! # Architecture
//!
//! * [`ScratchPool`] shards idle buffers across cache-line-aligned buckets
//!   indexed by thread shard id, eliminating single-mutex contention under
//!   parallel Rayon workloads while strictly enforcing `max_idle` and
//!   `max_buffer_bytes`.
//! * [`ComputeReusePool`] unifies typed scratch pools (`f64`, `f32`,
//!   `Complex64`, `usize`, `i64`, `u8`) with bounded LRU [`PlanCache`] stores
//!   for FFT plans, quadrature rules, Chebyshev nodes, trig/twiddle tables,
//!   combinatorial/number-theoretic tables, sparse fill-reducing orderings, and
//!   an extensible [`PlanRegistry`] of named caches.
//! * [`global_pool`] provides the process-wide singleton augmented with
//!   fast thread-local scratch slot caches for zero-lock steady-state reuse.

use std::any::Any;
use std::any::TypeId;
use std::cell::RefCell;
use std::collections::VecDeque;
use std::f64::consts::PI;
use std::fmt;
use std::mem;
use std::ops::Deref;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::MutexGuard;
use std::sync::PoisonError;
use std::sync::RwLock;
use std::sync::TryLockError;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use blake3::Hasher;
use num_complex::Complex64;
use once_cell::sync::OnceCell;
use rustfft::Fft;
use rustfft::FftDirection;
use rustfft::FftPlanner;

const SCRATCH_SHARDS: usize = 16;
const SCRATCH_SHARD_MASK: usize = SCRATCH_SHARDS - 1;
const TLS_SLOT_COUNT: usize = 4;
const COMPLEX64_SCRATCH_BYTES: usize = 4 * 1024 * 1024;
const LARGE_COMPLEX64_SCRATCH_BYTES: usize = 16 * 1024 * 1024;
const DEFAULT_SCRATCH_BYTES: usize = 4 * 1024 * 1024;
const DEFAULT_POOL_MAX_IDLE: usize = 32;

/// Snapshot of hit, miss, eviction, and occupancy counters for a cache or pool.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CacheStats {
    /// Lookups satisfied from idle or cached entries.
    pub hits: u64,
    /// Lookups that missed and built or allocated a fresh value.
    pub misses: u64,
    /// Entries or buffers removed due to capacity limits, replacement, or trimming.
    pub evictions: u64,
    /// Number of entries or idle buffers currently stored.
    pub entries: usize,
    /// Estimated retained capacity bytes across idle buffers (zero for plan caches).
    pub retained_bytes: usize,
}

impl CacheStats {
    fn accumulate(
        &mut self,
        other: Self,
    ) {
        self.hits = self.hits.saturating_add(other.hits);
        self.misses = self.misses.saturating_add(other.misses);
        self.evictions = self.evictions.saturating_add(other.evictions);
        self.entries = self.entries.saturating_add(other.entries);
        self.retained_bytes = self.retained_bytes.saturating_add(other.retained_bytes);
    }
}

/// Aggregate telemetry for a [`ComputeReusePool`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReusePoolStats {
    /// Combined hits across all scratch pools and plan/table caches.
    pub hits: u64,
    /// Combined misses across all scratch pools and plan/table caches.
    pub misses: u64,
    /// Combined evictions across all scratch pools and plan/table caches.
    pub evictions: u64,
    /// Total bytes retained in idle scratch buffers (shared shards plus active global TLS).
    pub retained_scratch_bytes: usize,
    /// Total number of idle scratch buffers retained in shared shards.
    pub retained_scratch_buffers: usize,
    /// Total number of entries retained across built-in and named plan caches.
    pub cached_plans: usize,
    /// Breakdown of scratch pool counters.
    pub scratch: CacheStats,
    /// Breakdown of built-in and named plan/table cache counters.
    pub plans: CacheStats,
}

/// Quadrature family key for [`ComputeReusePool::get_or_init_quadrature`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QuadratureRuleKind {
    /// Gauss-Legendre nodes and weights on `[-1, 1]`.
    GaussLegendre,
    /// Gauss-Kronrod nodes and weights on `[-1, 1]`.
    GaussKronrod,
    /// Clenshaw-Curtis nodes and weights on `[-1, 1]`.
    ClenshawCurtis,
    /// Caller-defined quadrature rule family identified by a numeric tag.
    Custom(u64),
}

/// Ordering algorithm family key for [`ComputeReusePool::get_or_build_sparse_ordering`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SparseOrderingKind {
    /// Symmetric Reverse Cuthill-McKee bandwidth-reducing ordering.
    ReverseCuthillMckee,
    /// Approximate minimum degree ordering.
    ApproximateMinimumDegree,
    /// Nested dissection ordering.
    NestedDissection,
    /// Caller-defined ordering family identified by a numeric tag.
    Custom(u64),
}

/// Cached fill-reducing permutation and inverse permutation for a sparse pattern.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseOrderingPlan {
    /// Matrix dimension `n`.
    pub dim: usize,
    /// Forward permutation mapping `perm[new] = old`.
    pub perm: Vec<usize>,
    /// Inverse permutation mapping `inv_perm[old] = new`.
    pub inv_perm: Vec<usize>,
}

#[derive(Clone, PartialEq, Eq)]
struct SparseOrderingKey {
    kind: SparseOrderingKind,
    dim: usize,
    nnz: usize,
    digest: [u8; 32],
}

#[repr(align(64))]
struct ScratchShard<T> {
    idle: Mutex<Vec<Vec<T>>>,
}

impl<T> ScratchShard<T> {
    const fn new() -> Self {
        Self {
            idle: Mutex::new(Vec::new()),
        }
    }
}

struct TlsOps<T> {
    try_pop_any: fn() -> Option<Vec<T>>,
    try_pop_for: fn(usize) -> Option<Vec<T>>,
    try_push: fn(Vec<T>) -> Option<Vec<T>>,
}

impl<T> Clone for TlsOps<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for TlsOps<T> {}

fn lock_shard<T>(shard: &ScratchShard<T>) -> MutexGuard<'_, Vec<Vec<T>>> {
    shard.idle.lock().unwrap_or_else(PoisonError::into_inner)
}

fn try_lock_shard<T>(shard: &ScratchShard<T>) -> Result<MutexGuard<'_, Vec<Vec<T>>>, ()> {
    match shard.idle.try_lock() {
        | Ok(guard) => Ok(guard),
        | Err(TryLockError::Poisoned(error)) => Ok(error.into_inner()),
        | Err(TryLockError::WouldBlock) => Err(()),
    }
}

/// Sharded pool of empty vectors that retain their heap allocations between
/// computations.
///
/// Concurrent threads hash to cache-line-padded shards with neighbor probing,
/// eliminating single-mutex contention while enforcing global `max_idle` and
/// `max_buffer_bytes` retention caps.
pub struct ScratchPool<T> {
    shards: [ScratchShard<T>; SCRATCH_SHARDS],
    max_idle: usize,
    max_buffer_bytes: usize,
    idle_count: AtomicUsize,
    idle_bytes: AtomicUsize,
    hits: AtomicU64,
    misses: AtomicU64,
    evictions: AtomicU64,
    tls: Option<TlsOps<T>>,
}

impl<T> fmt::Debug for ScratchPool<T> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("ScratchPool")
            .field("max_idle", &self.max_idle)
            .field("max_buffer_bytes", &self.max_buffer_bytes)
            .field("idle_count", &self.len())
            .field("idle_bytes", &self.retained_bytes())
            .finish_non_exhaustive()
    }
}

impl<T> ScratchPool<T> {
    /// Retains up to `max_idle` buffers across all shards, each at most
    /// `max_buffer_bytes` bytes. A zero limit disables retention.
    #[must_use]
    pub const fn new(
        max_idle: usize,
        max_buffer_bytes: usize,
    ) -> Self {
        Self::with_optional_tls(max_idle, max_buffer_bytes, None)
    }

    const fn with_optional_tls(
        max_idle: usize,
        max_buffer_bytes: usize,
        tls: Option<TlsOps<T>>,
    ) -> Self {
        Self {
            shards: [
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
                ScratchShard::new(),
            ],
            max_idle,
            max_buffer_bytes,
            idle_count: AtomicUsize::new(0),
            idle_bytes: AtomicUsize::new(0),
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
            evictions: AtomicU64::new(0),
            tls,
        }
    }

    const fn with_tls(
        max_idle: usize,
        max_buffer_bytes: usize,
        tls: TlsOps<T>,
    ) -> Self {
        Self::with_optional_tls(max_idle, max_buffer_bytes, Some(tls))
    }

    /// Returns the maximum number of idle buffers this pool retains.
    #[must_use]
    pub const fn max_idle(&self) -> usize {
        self.max_idle
    }

    /// Returns the maximum allocation size in bytes retained per buffer.
    #[must_use]
    pub const fn max_buffer_bytes(&self) -> usize {
        self.max_buffer_bytes
    }

    /// Returns the current number of idle buffers stored in shared shards.
    #[must_use]
    pub fn len(&self) -> usize {
        self.idle_count.load(Ordering::Relaxed)
    }

    /// Returns `true` when no idle buffers are stored in shared shards.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the total capacity bytes currently retained in shared shards.
    #[must_use]
    pub fn retained_bytes(&self) -> usize {
        self.idle_bytes.load(Ordering::Relaxed)
    }

    /// Returns a snapshot of hit, miss, eviction, and occupancy statistics.
    #[must_use]
    pub fn stats(&self) -> CacheStats {
        CacheStats {
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
            evictions: self.evictions.load(Ordering::Relaxed),
            entries: self.len(),
            retained_bytes: self.retained_bytes(),
        }
    }

    /// Resets the cumulative `hits`, `misses`, and `evictions` counters to zero.
    pub fn reset_stats(&self) {
        self.hits.store(0, Ordering::Relaxed);
        self.misses.store(0, Ordering::Relaxed);
        self.evictions.store(0, Ordering::Relaxed);
    }

    /// Drops all currently retained idle buffers from every shard without
    /// incrementing eviction counters.
    pub fn clear(&self) {
        let _ = self.drain_all_shards(false);
    }

    /// Evicts all currently retained idle buffers from every shard and records
    /// them as evictions.
    pub fn trim_idle(&self) {
        let _ = self.drain_all_shards(true);
    }

    /// Evicts the largest idle buffers until [`Self::retained_bytes`] is at most
    /// `max_bytes`.
    pub fn trim_to_bytes(
        &self,
        max_bytes: usize,
    ) {
        while self.retained_bytes() > max_bytes {
            if self.evict_largest_one() == 0 {
                break;
            }
        }
    }

    /// Takes an empty vector. Resize or reserve it for the current job.
    /// Contention on a shard probes neighbor shards before allocating on demand.
    #[must_use]
    pub fn checkout(&self) -> Scratch<'_, T> {
        if self.max_idle == 0 || self.max_buffer_bytes == 0 || mem::size_of::<T>() == 0 {
            self.misses.fetch_add(1, Ordering::Relaxed);
            return Scratch {
                pool: self,
                buffer: Vec::new(),
            };
        }
        if let Some(tls) = self.tls
            && let Some(buffer) = (tls.try_pop_any)()
        {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Scratch { pool: self, buffer };
        }
        if let Some(buffer) = self.pop_from_shards(None) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Scratch { pool: self, buffer };
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        Scratch {
            pool: self,
            buffer: Vec::new(),
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
        if minimum == 0 {
            return Scratch {
                pool: self,
                buffer: Vec::new(),
            };
        }
        if self.max_idle == 0
            || mem::size_of::<T>() == 0
            || minimum.saturating_mul(mem::size_of::<T>()) > self.max_buffer_bytes
        {
            self.misses.fetch_add(1, Ordering::Relaxed);
            return Scratch {
                pool: self,
                buffer: Vec::with_capacity(minimum),
            };
        }
        if let Some(tls) = self.tls
            && let Some(buffer) = (tls.try_pop_for)(minimum)
        {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Scratch { pool: self, buffer };
        }
        if let Some(buffer) = self.pop_from_shards(Some(minimum)) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Scratch { pool: self, buffer };
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        Scratch {
            pool: self,
            buffer: Vec::with_capacity(minimum),
        }
    }

    fn home_shard(&self) -> usize {
        if self.max_idle <= 2 {
            0
        } else {
            thread_local! {
                static SHARD_INDEX: usize = {
                    static NEXT_SHARD: AtomicUsize = AtomicUsize::new(0);
                    NEXT_SHARD.fetch_add(1, Ordering::Relaxed) & SCRATCH_SHARD_MASK
                };
            }
            SHARD_INDEX.with(|&index| index)
        }
    }

    fn take_from_guard(
        &self,
        idle: &mut Vec<Vec<T>>,
        minimum: Option<usize>,
    ) -> Option<Vec<T>> {
        let taken = match minimum {
            | None => idle.pop(),
            | Some(min_cap) => {
                let best = idle
                    .iter()
                    .enumerate()
                    .filter(|(_, buffer)| buffer.capacity() >= min_cap)
                    .min_by_key(|(_, buffer)| buffer.capacity())
                    .map(|(index, _)| index);
                best.map(|index| idle.swap_remove(index))
            },
        };
        if let Some(ref buf) = taken {
            let bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
            self.idle_count.fetch_sub(1, Ordering::Release);
            self.idle_bytes.fetch_sub(bytes, Ordering::Relaxed);
        }
        taken
    }

    fn pop_from_shards(
        &self,
        minimum: Option<usize>,
    ) -> Option<Vec<T>> {
        if self.idle_count.load(Ordering::Acquire) == 0 {
            return None;
        }
        let home = self.home_shard();
        let home_shard = self.shards.get(home)?;
        let mut home_contended = false;
        match try_lock_shard(home_shard) {
            | Ok(mut guard) => {
                if let Some(buf) = self.take_from_guard(&mut guard, minimum) {
                    return Some(buf);
                }
            },
            | Err(()) => {
                home_contended = true;
            },
        }
        if self.max_idle > 2 {
            let pool_full = self.idle_count.load(Ordering::Acquire) >= self.max_idle;
            for offset in 1..SCRATCH_SHARDS {
                if self.idle_count.load(Ordering::Acquire) == 0 {
                    return None;
                }
                let idx = home.wrapping_add(offset) & SCRATCH_SHARD_MASK;
                if let Some(shard) = self.shards.get(idx)
                    && let Ok(mut guard) = try_lock_shard(shard)
                    && (pool_full || guard.len() > 1)
                    && let Some(buf) = self.take_from_guard(&mut guard, minimum)
                {
                    return Some(buf);
                }
            }
        }
        if home_contended && self.idle_count.load(Ordering::Acquire) > 0 {
            let mut guard = lock_shard(home_shard);
            return self.take_from_guard(&mut guard, minimum);
        }
        None
    }

    #[allow(clippy::single_call_fn)]
    // Keeps Scratch::drop focused on retention guards and TLS dispatch.
    fn return_to_shards(
        &self,
        mut buffer: Vec<T>,
    ) {
        let buf_bytes = buffer.capacity().saturating_mul(mem::size_of::<T>());
        let reserved = self
            .idle_count
            .try_update(Ordering::AcqRel, Ordering::Acquire, |count| {
                (count < self.max_idle).then(|| count.saturating_add(1))
            })
            .is_ok();
        let home = self.home_shard();
        if reserved {
            if let Some(home_shard) = self.shards.get(home)
                && let Ok(mut guard) = try_lock_shard(home_shard)
            {
                guard.push(buffer);
                self.idle_bytes.fetch_add(buf_bytes, Ordering::Relaxed);
                return;
            }
            if self.max_idle > 2 {
                for offset in 1..=4 {
                    let idx = home.wrapping_add(offset) & SCRATCH_SHARD_MASK;
                    if let Some(shard) = self.shards.get(idx)
                        && let Ok(mut guard) = try_lock_shard(shard)
                    {
                        guard.push(buffer);
                        self.idle_bytes.fetch_add(buf_bytes, Ordering::Relaxed);
                        return;
                    }
                }
            }
            if let Some(home_shard) = self.shards.get(home) {
                let mut guard = lock_shard(home_shard);
                guard.push(buffer);
                self.idle_bytes.fetch_add(buf_bytes, Ordering::Relaxed);
                return;
            }
            self.idle_count.fetch_sub(1, Ordering::Release);
            return;
        }
        // Pool is at max_idle: swap with a strictly smaller retained buffer if present.
        let shard_probe_count = if self.max_idle <= 2 { 1 } else { 4 };
        for offset in 0..shard_probe_count {
            let idx = home.wrapping_add(offset) & SCRATCH_SHARD_MASK;
            let Some(shard) = self.shards.get(idx) else {
                continue;
            };
            let Ok(mut guard) = try_lock_shard(shard) else {
                continue;
            };
            if let Some(smallest) = guard.iter_mut().min_by_key(|item| item.capacity())
                && smallest.capacity() < buffer.capacity()
            {
                let old_bytes = smallest.capacity().saturating_mul(mem::size_of::<T>());
                mem::swap(smallest, &mut buffer);
                self.idle_bytes
                    .fetch_add(buf_bytes.saturating_sub(old_bytes), Ordering::Relaxed);
                drop(guard);
                drop(buffer);
                self.evictions.fetch_add(1, Ordering::Relaxed);
                return;
            }
        }
        drop(buffer);
        self.evictions.fetch_add(1, Ordering::Relaxed);
    }

    fn drain_all_shards(
        &self,
        record_evictions: bool,
    ) -> usize {
        let mut total_freed = 0_usize;
        let mut total_count = 0_usize;
        for shard in &self.shards {
            let drained = {
                let mut guard = lock_shard(shard);
                if guard.is_empty() {
                    Vec::new()
                } else {
                    let taken = mem::take(&mut *guard);
                    let count = taken.len();
                    let bytes: usize = taken
                        .iter()
                        .map(|buf| buf.capacity().saturating_mul(mem::size_of::<T>()))
                        .sum();
                    self.idle_count.fetch_sub(count, Ordering::Release);
                    self.idle_bytes.fetch_sub(bytes, Ordering::Relaxed);
                    total_freed = total_freed.saturating_add(bytes);
                    total_count = total_count.saturating_add(count);
                    taken
                }
            };
            drop(drained);
        }
        if record_evictions && total_count > 0 {
            self.evictions
                .fetch_add(total_count as u64, Ordering::Relaxed);
        }
        total_freed
    }

    fn evict_largest_one(&self) -> usize {
        let mut best_shard: Option<usize> = None;
        let mut best_cap: usize = 0;
        for (shard_idx, shard) in self.shards.iter().enumerate() {
            let guard = lock_shard(shard);
            if let Some(max_buf) = guard.iter().max_by_key(|buf| buf.capacity())
                && max_buf.capacity() > best_cap
            {
                best_cap = max_buf.capacity();
                best_shard = Some(shard_idx);
            }
        }
        let Some(shard_idx) = best_shard else {
            return 0;
        };
        let Some(shard) = self.shards.get(shard_idx) else {
            return 0;
        };
        let removed = {
            let mut guard = lock_shard(shard);
            let pos = guard
                .iter()
                .enumerate()
                .max_by_key(|(_, buf)| buf.capacity())
                .map(|(idx, _)| idx);
            let taken = pos.map(|idx| guard.swap_remove(idx));
            if let Some(ref buf) = taken {
                let bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
                self.idle_count.fetch_sub(1, Ordering::Release);
                self.idle_bytes.fetch_sub(bytes, Ordering::Relaxed);
            }
            taken
        };
        removed.map_or(0, |buf| {
            let bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
            self.evictions.fetch_add(1, Ordering::Relaxed);
            drop(buf);
            bytes
        })
    }
}

/// Exclusive scratch buffer returned to its pool (or thread-local cache) when dropped.
pub struct Scratch<'a, T> {
    pool: &'a ScratchPool<T>,
    buffer: Vec<T>,
}

impl<T> fmt::Debug for Scratch<'_, T> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("Scratch")
            .field("len", &self.buffer.len())
            .field("capacity", &self.buffer.capacity())
            .finish()
    }
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
        // Drop elements outside any lock. Only the allocation is reused.
        self.buffer.clear();
        let mut buf = mem::take(&mut self.buffer);
        if let Some(tls) = self.pool.tls {
            let Some(overflow) = (tls.try_push)(buf) else {
                return;
            };
            buf = overflow;
        }
        self.pool.return_to_shards(buf);
    }
}

struct TlsBucket<T> {
    epoch: u64,
    slots: [Vec<T>; TLS_SLOT_COUNT],
    len: usize,
}

impl<T> TlsBucket<T> {
    const fn new() -> Self {
        Self {
            epoch: 0,
            slots: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            len: 0,
        }
    }

    fn sync_epoch(&mut self) {
        let current = GLOBAL_POOL.tls_epoch.load(Ordering::Acquire);
        if self.epoch != current {
            self.clear_all();
            self.epoch = current;
        }
    }

    fn clear_all(&mut self) {
        let mut freed_bytes = 0_usize;
        for slot in self.slots.iter_mut().take(self.len) {
            let cap_bytes = slot.capacity().saturating_mul(mem::size_of::<T>());
            freed_bytes = freed_bytes.saturating_add(cap_bytes);
            *slot = Vec::new();
        }
        self.len = 0;
        if freed_bytes > 0 {
            GLOBAL_POOL
                .tls_bytes
                .fetch_sub(freed_bytes, Ordering::Relaxed);
        }
    }

    fn pop_any(&mut self) -> Option<Vec<T>> {
        self.sync_epoch();
        if self.len == 0 {
            return None;
        }
        let idx = self.len.saturating_sub(1);
        self.len = idx;
        let buf = self.slots.get_mut(idx).map(mem::take)?;
        let bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
        GLOBAL_POOL.tls_bytes.fetch_sub(bytes, Ordering::Relaxed);
        Some(buf)
    }

    fn pop_for(
        &mut self,
        minimum: usize,
    ) -> Option<Vec<T>> {
        self.sync_epoch();
        let best = self
            .slots
            .iter()
            .take(self.len)
            .enumerate()
            .filter(|(_, buf)| buf.capacity() >= minimum)
            .min_by_key(|(_, buf)| buf.capacity())
            .map(|(idx, _)| idx)?;
        let last = self.len.saturating_sub(1);
        self.slots.swap(best, last);
        self.len = last;
        let buf = self.slots.get_mut(last).map(mem::take)?;
        let bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
        GLOBAL_POOL.tls_bytes.fetch_sub(bytes, Ordering::Relaxed);
        Some(buf)
    }

    fn push(
        &mut self,
        mut buf: Vec<T>,
    ) -> Option<Vec<T>> {
        self.sync_epoch();
        let buf_bytes = buf.capacity().saturating_mul(mem::size_of::<T>());
        if self.len < TLS_SLOT_COUNT {
            if let Some(slot) = self.slots.get_mut(self.len) {
                *slot = buf;
                self.len = self.len.saturating_add(1);
                GLOBAL_POOL
                    .tls_bytes
                    .fetch_add(buf_bytes, Ordering::Relaxed);
                return None;
            }
            return Some(buf);
        }
        if let Some(smallest) = self
            .slots
            .iter_mut()
            .min_by_key(|candidate| candidate.capacity())
            && smallest.capacity() < buf.capacity()
        {
            let old_bytes = smallest.capacity().saturating_mul(mem::size_of::<T>());
            mem::swap(smallest, &mut buf);
            GLOBAL_POOL
                .tls_bytes
                .fetch_add(buf_bytes.saturating_sub(old_bytes), Ordering::Relaxed);
        }
        Some(buf)
    }
}

impl<T> Drop for TlsBucket<T> {
    fn drop(&mut self) {
        self.clear_all();
    }
}

macro_rules! define_tls_slot {
    ($mod_name:ident, $ty:ty, $pop_any:ident, $pop_for:ident, $push:ident, $clear:ident) => {
        thread_local! {
            static $mod_name: RefCell<TlsBucket<$ty>> = const { RefCell::new(TlsBucket::new()) };
        }

        fn $pop_any() -> Option<Vec<$ty>> {
            $mod_name
                .try_with(|cell| {
                    cell.try_borrow_mut()
                        .ok()
                        .and_then(|mut bucket| bucket.pop_any())
                })
                .ok()
                .flatten()
        }

        fn $pop_for(minimum: usize) -> Option<Vec<$ty>> {
            $mod_name
                .try_with(|cell| {
                    cell.try_borrow_mut()
                        .ok()
                        .and_then(|mut bucket| bucket.pop_for(minimum))
                })
                .ok()
                .flatten()
        }

        fn $push(buf: Vec<$ty>) -> Option<Vec<$ty>> {
            match $mod_name
                .try_with(|cell| cell.try_borrow_mut().map(|mut bucket| bucket.push(buf)))
            {
                | Ok(Ok(res)) => res,
                | Ok(Err(_)) => None,
                | Err(_) => None,
            }
        }

        fn $clear() {
            let _ = $mod_name.try_with(|cell| {
                if let Ok(mut bucket) = cell.try_borrow_mut() {
                    bucket.clear_all();
                    bucket.epoch = GLOBAL_POOL.tls_epoch.load(Ordering::Acquire);
                }
            });
        }
    };
}

define_tls_slot!(
    TLS_F64,
    f64,
    tls_pop_any_f64,
    tls_pop_for_f64,
    tls_push_f64,
    tls_clear_f64
);
define_tls_slot!(
    TLS_F32,
    f32,
    tls_pop_any_f32,
    tls_pop_for_f32,
    tls_push_f32,
    tls_clear_f32
);
define_tls_slot!(
    TLS_COMPLEX64,
    Complex64,
    tls_pop_any_complex64,
    tls_pop_for_complex64,
    tls_push_complex64,
    tls_clear_complex64
);
define_tls_slot!(
    TLS_USIZE,
    usize,
    tls_pop_any_usize,
    tls_pop_for_usize,
    tls_push_usize,
    tls_clear_usize
);
define_tls_slot!(
    TLS_I64,
    i64,
    tls_pop_any_i64,
    tls_pop_for_i64,
    tls_push_i64,
    tls_clear_i64
);
define_tls_slot!(
    TLS_BYTES,
    u8,
    tls_pop_any_bytes,
    tls_pop_for_bytes,
    tls_push_bytes,
    tls_clear_bytes
);

fn clear_current_thread_tls() {
    tls_clear_f64();
    tls_clear_f32();
    tls_clear_complex64();
    tls_clear_usize();
    tls_clear_i64();
    tls_clear_bytes();
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

/// Small bounded cache of immutable plans supporting FIFO or LRU eviction and
/// coalesced concurrent construction.
pub struct PlanCache<K, V: ?Sized> {
    entries: Mutex<PlanEntries<K, V>>,
    capacity: usize,
    promote_hits: bool,
    hits: AtomicU64,
    misses: AtomicU64,
    evictions: AtomicU64,
}

impl<K, V: ?Sized> fmt::Debug for PlanCache<K, V> {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("PlanCache")
            .field("capacity", &self.capacity)
            .field("promote_hits", &self.promote_hits)
            .field("hits", &self.hits.load(Ordering::Relaxed))
            .field("misses", &self.misses.load(Ordering::Relaxed))
            .field("evictions", &self.evictions.load(Ordering::Relaxed))
            .finish_non_exhaustive()
    }
}

impl<K: Eq, V: ?Sized> PlanCache<K, V> {
    /// Creates an empty cache that evicts in insertion order (FIFO).
    /// A capacity of zero builds each plan on demand.
    #[must_use]
    pub const fn new(capacity: usize) -> Self {
        Self::with_policy(capacity, false)
    }

    /// Creates an empty bounded cache that evicts the least recently used
    /// available entry. A capacity of zero builds each plan on demand.
    #[must_use]
    pub const fn new_lru(capacity: usize) -> Self {
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
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
            evictions: AtomicU64::new(0),
        }
    }

    /// Returns the maximum number of entries retained by this cache.
    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.capacity
    }

    /// Returns the number of entries currently held in this cache.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .queue
            .len()
    }

    /// Returns `true` when the cache holds no entries.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns `true` if `key` currently has an initialized plan in the cache.
    /// Does not modify LRU recency or telemetry counters.
    #[must_use]
    pub fn contains(
        &self,
        key: &K,
    ) -> bool {
        if self.capacity == 0 {
            return false;
        }
        self.entries
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .queue
            .iter()
            .any(|entry| entry.key == *key && entry.slot.get().is_some())
    }

    /// Looks up a previously initialized plan for `key`, promoting LRU recency
    /// on hit.
    #[must_use]
    pub fn get(
        &self,
        key: &K,
    ) -> Option<Arc<V>> {
        if self.capacity == 0 {
            self.misses.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        let mut entries = self.entries.lock().unwrap_or_else(PoisonError::into_inner);
        Self::normalize_clock_if_needed(&mut entries, self.promote_hits);
        let stamp = entries.clock;
        if let Some(entry) = entries.queue.iter_mut().find(|entry| entry.key == *key)
            && let Some(plan) = entry.slot.get()
        {
            entry.last_access = stamp;
            let cloned = Arc::clone(plan);
            entries.clock = entries.clock.saturating_add(u64::from(self.promote_hits));
            drop(entries);
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Some(cloned);
        }
        drop(entries);
        self.misses.fetch_add(1, Ordering::Relaxed);
        None
    }

    /// Removes all entries from the cache without invalidating any `Arc`
    /// already held by callers.
    pub fn clear(&self) {
        let drained = {
            let mut entries = self.entries.lock().unwrap_or_else(PoisonError::into_inner);
            mem::take(&mut entries.queue)
        };
        drop(drained);
    }

    /// Evicts cached entries that have no external `Arc` holders (and abandons
    /// unreferenced empty retry slots), keeping plans currently referenced by
    /// callers.
    pub fn trim_idle(&self) {
        let mut removed = Vec::new();
        {
            let mut entries = self.entries.lock().unwrap_or_else(PoisonError::into_inner);
            let mut index = 0_usize;
            while index < entries.queue.len() {
                let is_idle = entries.queue.get(index).is_some_and(|entry| {
                    let slot_unshared = Arc::strong_count(&entry.slot) == 1;
                    match entry.slot.get() {
                        | Some(plan) => slot_unshared && Arc::strong_count(plan) == 1,
                        | None => slot_unshared,
                    }
                });
                if is_idle {
                    if let Some(evicted) = entries.queue.remove(index) {
                        removed.push(evicted);
                    }
                } else {
                    index = index.saturating_add(1);
                }
            }
        }
        if !removed.is_empty() {
            self.evictions
                .fetch_add(removed.len() as u64, Ordering::Relaxed);
        }
        drop(removed);
    }

    /// Returns a snapshot of hit, miss, eviction, and occupancy statistics.
    #[must_use]
    pub fn stats(&self) -> CacheStats {
        CacheStats {
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
            evictions: self.evictions.load(Ordering::Relaxed),
            entries: self.len(),
            retained_bytes: 0,
        }
    }

    /// Resets the cumulative `hits`, `misses`, and `evictions` counters to zero.
    pub fn reset_stats(&self) {
        self.hits.store(0, Ordering::Relaxed);
        self.misses.store(0, Ordering::Relaxed);
        self.evictions.store(0, Ordering::Relaxed);
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
    ///
    /// # Errors
    /// Propagates the error returned by `build` when construction fails.
    pub fn try_get_or_insert_with<E>(
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

    fn normalize_clock_if_needed(
        entries: &mut PlanEntries<K, V>,
        promote_hits: bool,
    ) {
        if promote_hits && entries.clock == u64::MAX {
            // Preserve order when shrinking the timestamps back into range.
            let queue = entries.queue.make_contiguous();
            queue.sort_unstable_by_key(|entry| entry.last_access);
            for (index, entry) in queue.iter_mut().enumerate() {
                entry.last_access = index as u64;
            }
            entries.clock = queue.len() as u64;
        }
    }

    fn lookup(
        &self,
        key: K,
    ) -> PlanAccess<V> {
        if self.capacity == 0 {
            self.misses.fetch_add(1, Ordering::Relaxed);
            return PlanAccess::Uncached;
        }
        let mut entries = self.entries.lock().unwrap_or_else(PoisonError::into_inner);
        Self::normalize_clock_if_needed(&mut entries, self.promote_hits);
        let stamp = entries.clock;
        entries.clock = entries.clock.saturating_add(u64::from(self.promote_hits));
        if let Some(entry) = entries.queue.iter_mut().find(|entry| entry.key == key) {
            entry.last_access = stamp;
            let access = entry.slot.get().map_or_else(
                || PlanAccess::Pending(Arc::clone(&entry.slot)),
                |plan| PlanAccess::Ready(Arc::clone(plan)),
            );
            drop(entries);
            self.hits.fetch_add(1, Ordering::Relaxed);
            return access;
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
                drop(entries);
                self.misses.fetch_add(1, Ordering::Relaxed);
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
        if evicted.is_some() {
            self.evictions.fetch_add(1, Ordering::Relaxed);
        }
        drop(evicted);
        self.misses.fetch_add(1, Ordering::Relaxed);
        PlanAccess::Pending(slot)
    }
}

trait ManagedPlanCache: Send + Sync {
    fn managed_stats(&self) -> CacheStats;
    fn managed_len(&self) -> usize;
    fn managed_clear(&self);
    fn managed_trim_idle(&self);
    fn managed_reset_stats(&self);
}

impl<K: Eq + Send + Sync, V: Send + Sync + ?Sized> ManagedPlanCache for PlanCache<K, V> {
    fn managed_stats(&self) -> CacheStats {
        self.stats()
    }

    fn managed_len(&self) -> usize {
        self.len()
    }

    fn managed_clear(&self) {
        self.clear();
    }

    fn managed_trim_idle(&self) {
        self.trim_idle();
    }

    fn managed_reset_stats(&self) {
        self.reset_stats();
    }
}

struct RegistryEntry {
    name: String,
    key_type: TypeId,
    val_type: TypeId,
    any_cache: Arc<dyn Any + Send + Sync>,
    managed: Arc<dyn ManagedPlanCache>,
}

/// Thread-safe registry of named, bounded LRU [`PlanCache`] instances.
///
/// Each distinct `(name, TypeId::<K>(), TypeId::<V>())` triple resolves to an
/// isolated [`PlanCache<K, V>`] whose lifecycle and telemetry are aggregated by
/// the registry.
pub struct PlanRegistry {
    caches: RwLock<Vec<RegistryEntry>>,
}

impl Default for PlanRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for PlanRegistry {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("PlanRegistry")
            .field("cache_count", &self.cache_count())
            .field("entries", &self.len())
            .finish_non_exhaustive()
    }
}

impl PlanRegistry {
    /// Creates an empty plan registry.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            caches: RwLock::new(Vec::new()),
        }
    }

    /// Returns the named bounded LRU [`PlanCache<K, V>`], creating one with
    /// `capacity` if it does not yet exist.
    #[must_use]
    pub fn cache<K, V>(
        &self,
        name: &str,
        capacity: usize,
    ) -> Arc<PlanCache<K, V>>
    where
        K: Eq + Send + Sync + 'static,
        V: Send + Sync + ?Sized + 'static,
    {
        let key_type = TypeId::of::<K>();
        let val_type = TypeId::of::<V>();
        {
            let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
            if let Some(found) = Self::find_cache::<K, V>(&guard, name, key_type, val_type) {
                return found;
            }
        }
        let mut guard = self.caches.write().unwrap_or_else(PoisonError::into_inner);
        if let Some(found) = Self::find_cache::<K, V>(&guard, name, key_type, val_type) {
            return found;
        }
        let typed = Arc::new(PlanCache::<K, V>::new_lru(capacity));
        let any_copy = Arc::clone(&typed);
        let managed_copy = Arc::clone(&typed);
        let any_cache: Arc<dyn Any + Send + Sync> = any_copy;
        let managed: Arc<dyn ManagedPlanCache> = managed_copy;
        guard.push(RegistryEntry {
            name: name.to_owned(),
            key_type,
            val_type,
            any_cache,
            managed,
        });
        drop(guard);
        typed
    }

    /// Looks up or inserts a value into the named LRU cache for `(name, K, V)`.
    pub fn get_or_insert_with<K, V>(
        &self,
        name: &str,
        capacity: usize,
        key: K,
        build: impl FnOnce() -> Arc<V>,
    ) -> Arc<V>
    where
        K: Eq + Send + Sync + 'static,
        V: Send + Sync + ?Sized + 'static,
    {
        let cache = self.cache::<K, V>(name, capacity);
        cache.get_or_insert_with(key, build)
    }

    /// Returns the number of registered named caches.
    #[must_use]
    pub fn cache_count(&self) -> usize {
        self.caches
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .len()
    }

    /// Returns the total number of cached plan entries across all registered
    /// named caches.
    #[must_use]
    pub fn len(&self) -> usize {
        let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
        guard.iter().map(|entry| entry.managed.managed_len()).sum()
    }

    /// Returns `true` when no entries are stored in any registered cache.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Clears all entries in every registered named cache.
    pub fn clear(&self) {
        let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
        for entry in guard.iter() {
            entry.managed.managed_clear();
        }
    }

    /// Evicts idle entries (with no external `Arc` holders) in every registered
    /// named cache.
    pub fn trim_idle(&self) {
        let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
        for entry in guard.iter() {
            entry.managed.managed_trim_idle();
        }
    }

    /// Aggregates statistics across all registered named caches.
    #[must_use]
    pub fn stats(&self) -> CacheStats {
        let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
        let mut total = CacheStats::default();
        for entry in guard.iter() {
            total.accumulate(entry.managed.managed_stats());
        }
        total
    }

    /// Resets statistics counters on all registered named caches.
    pub fn reset_stats(&self) {
        let guard = self.caches.read().unwrap_or_else(PoisonError::into_inner);
        for entry in guard.iter() {
            entry.managed.managed_reset_stats();
        }
    }

    fn find_cache<K, V>(
        entries: &[RegistryEntry],
        name: &str,
        key_type: TypeId,
        val_type: TypeId,
    ) -> Option<Arc<PlanCache<K, V>>>
    where
        K: Eq + Send + Sync + 'static,
        V: Send + Sync + ?Sized + 'static,
    {
        entries
            .iter()
            .find(|entry| {
                entry.name == name && entry.key_type == key_type && entry.val_type == val_type
            })
            .and_then(|entry| {
                Arc::clone(&entry.any_cache)
                    .downcast::<PlanCache<K, V>>()
                    .ok()
            })
    }
}

/// Unified compute reuse pool coordinating sharded scratch arenas, bounded LRU
/// plan/table caches, and global thread-local fast-path slots.
pub struct ComputeReusePool {
    f64_pool: ScratchPool<f64>,
    f32_pool: ScratchPool<f32>,
    complex64_pool: ScratchPool<Complex64>,
    complex64_large_pool: ScratchPool<Complex64>,
    usize_pool: ScratchPool<usize>,
    i64_pool: ScratchPool<i64>,
    bytes_pool: ScratchPool<u8>,
    fft_small_plans: PlanCache<(usize, FftDirection), dyn Fft<f64>>,
    fft_large_plans: PlanCache<(usize, FftDirection), dyn Fft<f64>>,
    quadrature_rules: PlanCache<(QuadratureRuleKind, usize), (Vec<f64>, Vec<f64>)>,
    chebyshev_nodes_cache: PlanCache<(bool, usize), Vec<f64>>,
    real_twiddles_cache: PlanCache<usize, [Complex64]>,
    unit_roots_cache: PlanCache<usize, [Complex64]>,
    factorial_tables: PlanCache<usize, Vec<f64>>,
    log_factorial_tables: PlanCache<usize, Vec<f64>>,
    binomial_rows: PlanCache<usize, Vec<f64>>,
    stirling_tables: PlanCache<usize, Vec<f64>>,
    harmonic_tables: PlanCache<usize, Vec<f64>>,
    prime_sieves: PlanCache<usize, Vec<usize>>,
    sparse_orderings: PlanCache<SparseOrderingKey, SparseOrderingPlan>,
    registry: PlanRegistry,
    tls_epoch: AtomicU64,
    tls_bytes: AtomicUsize,
}

static GLOBAL_POOL: ComputeReusePool = ComputeReusePool::new_global();

/// Returns the process-wide shared [`ComputeReusePool`] with thread-local
/// fast-path scratch slots enabled.
#[must_use]
pub fn global_pool() -> &'static ComputeReusePool {
    &GLOBAL_POOL
}

impl Default for ComputeReusePool {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for ComputeReusePool {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let stats = self.stats();
        f.debug_struct("ComputeReusePool")
            .field("retained_scratch_bytes", &stats.retained_scratch_bytes)
            .field("cached_plans", &stats.cached_plans)
            .field("hits", &stats.hits)
            .field("misses", &stats.misses)
            .field("evictions", &stats.evictions)
            .finish_non_exhaustive()
    }
}

impl ComputeReusePool {
    /// Creates an isolated [`ComputeReusePool`] with its own sharded scratch
    /// arenas and bounded LRU plan/table caches.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            f64_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, DEFAULT_SCRATCH_BYTES),
            f32_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, DEFAULT_SCRATCH_BYTES),
            complex64_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, COMPLEX64_SCRATCH_BYTES),
            complex64_large_pool: ScratchPool::new(4, LARGE_COMPLEX64_SCRATCH_BYTES),
            usize_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, DEFAULT_SCRATCH_BYTES),
            i64_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, DEFAULT_SCRATCH_BYTES),
            bytes_pool: ScratchPool::new(DEFAULT_POOL_MAX_IDLE, DEFAULT_SCRATCH_BYTES),
            fft_small_plans: PlanCache::new_lru(16),
            fft_large_plans: PlanCache::new_lru(2),
            quadrature_rules: PlanCache::new_lru(64),
            chebyshev_nodes_cache: PlanCache::new_lru(64),
            real_twiddles_cache: PlanCache::new_lru(32),
            unit_roots_cache: PlanCache::new_lru(32),
            factorial_tables: PlanCache::new_lru(16),
            log_factorial_tables: PlanCache::new_lru(16),
            binomial_rows: PlanCache::new_lru(32),
            stirling_tables: PlanCache::new_lru(16),
            harmonic_tables: PlanCache::new_lru(16),
            prime_sieves: PlanCache::new_lru(16),
            sparse_orderings: PlanCache::new_lru(32),
            registry: PlanRegistry::new(),
            tls_epoch: AtomicU64::new(0),
            tls_bytes: AtomicUsize::new(0),
        }
    }

    #[allow(clippy::single_call_fn)] // initializes the static GLOBAL_POOL with TLS fast-path hooks
    const fn new_global() -> Self {
        Self {
            f64_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                DEFAULT_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_f64,
                    try_pop_for: tls_pop_for_f64,
                    try_push: tls_push_f64,
                },
            ),
            f32_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                DEFAULT_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_f32,
                    try_pop_for: tls_pop_for_f32,
                    try_push: tls_push_f32,
                },
            ),
            complex64_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                COMPLEX64_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_complex64,
                    try_pop_for: tls_pop_for_complex64,
                    try_push: tls_push_complex64,
                },
            ),
            complex64_large_pool: ScratchPool::new(4, LARGE_COMPLEX64_SCRATCH_BYTES),
            usize_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                DEFAULT_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_usize,
                    try_pop_for: tls_pop_for_usize,
                    try_push: tls_push_usize,
                },
            ),
            i64_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                DEFAULT_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_i64,
                    try_pop_for: tls_pop_for_i64,
                    try_push: tls_push_i64,
                },
            ),
            bytes_pool: ScratchPool::with_tls(
                DEFAULT_POOL_MAX_IDLE,
                DEFAULT_SCRATCH_BYTES,
                TlsOps {
                    try_pop_any: tls_pop_any_bytes,
                    try_pop_for: tls_pop_for_bytes,
                    try_push: tls_push_bytes,
                },
            ),
            fft_small_plans: PlanCache::new_lru(16),
            fft_large_plans: PlanCache::new_lru(2),
            quadrature_rules: PlanCache::new_lru(64),
            chebyshev_nodes_cache: PlanCache::new_lru(64),
            real_twiddles_cache: PlanCache::new_lru(32),
            unit_roots_cache: PlanCache::new_lru(32),
            factorial_tables: PlanCache::new_lru(16),
            log_factorial_tables: PlanCache::new_lru(16),
            binomial_rows: PlanCache::new_lru(32),
            stirling_tables: PlanCache::new_lru(16),
            harmonic_tables: PlanCache::new_lru(16),
            prime_sieves: PlanCache::new_lru(16),
            sparse_orderings: PlanCache::new_lru(32),
            registry: PlanRegistry::new(),
            tls_epoch: AtomicU64::new(0),
            tls_bytes: AtomicUsize::new(0),
        }
    }

    fn is_global(&self) -> bool {
        std::ptr::eq(self, global_pool())
    }

    /// Checks out an empty `f64` scratch buffer with room for at least `minimum` elements.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_f64(
        &self,
        minimum: usize,
    ) -> Scratch<'_, f64> {
        self.f64_pool.checkout_for(minimum)
    }

    /// Checks out an empty `f32` scratch buffer with room for at least `minimum` elements.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_f32(
        &self,
        minimum: usize,
    ) -> Scratch<'_, f32> {
        self.f32_pool.checkout_for(minimum)
    }

    /// Checks out an empty `Complex64` scratch buffer with room for at least
    /// `minimum` elements. Requests over 4 MiB use a dedicated large-buffer pool.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_complex64(
        &self,
        minimum: usize,
    ) -> Scratch<'_, Complex64> {
        const MAX_SMALL_POINTS: usize = COMPLEX64_SCRATCH_BYTES / mem::size_of::<Complex64>();
        if minimum > MAX_SMALL_POINTS {
            self.complex64_large_pool.checkout_for(minimum)
        } else {
            self.complex64_pool.checkout_for(minimum)
        }
    }

    /// Checks out an empty `usize` scratch buffer with room for at least `minimum` elements.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_usize(
        &self,
        minimum: usize,
    ) -> Scratch<'_, usize> {
        self.usize_pool.checkout_for(minimum)
    }

    /// Checks out an empty `i64` scratch buffer with room for at least `minimum` elements.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_i64(
        &self,
        minimum: usize,
    ) -> Scratch<'_, i64> {
        self.i64_pool.checkout_for(minimum)
    }

    /// Checks out an empty byte arena (`Vec<u8>`) with room for at least `minimum` bytes.
    ///
    /// # Panics
    /// Panics if `minimum` exceeds the maximum capacity of a vector.
    #[must_use]
    pub fn scratch_bytes(
        &self,
        minimum: usize,
    ) -> Scratch<'_, u8> {
        self.bytes_pool.checkout_for(minimum)
    }

    /// Plans and caches an `f64` complex FFT using `RustFFT`'s runtime CPU feature detection.
    #[must_use]
    pub fn fft_plan(
        &self,
        length: usize,
        direction: FftDirection,
    ) -> Arc<dyn Fft<f64>> {
        const LARGE_PLAN_MAX_POINTS: usize = 524_288;
        const LARGE_PLAN_MAX_INPUT_AND_SCRATCH: usize = 1_048_576;

        let build = || FftPlanner::new().plan_fft(length, direction);
        if length <= 65_536 {
            self.fft_small_plans
                .get_or_insert_with((length, direction), build)
        } else if length <= LARGE_PLAN_MAX_POINTS {
            self.fft_large_plans
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

    /// Returns cached `(nodes, weights)` for any [`QuadratureRuleKind`] and `order`,
    /// building with `build` on a cache miss.
    pub fn get_or_init_quadrature(
        &self,
        kind: QuadratureRuleKind,
        order: usize,
        build: impl FnOnce() -> Arc<(Vec<f64>, Vec<f64>)>,
    ) -> Arc<(Vec<f64>, Vec<f64>)> {
        self.quadrature_rules
            .get_or_insert_with((kind, order), build)
    }

    /// Computes and caches `n`-point Gauss-Legendre quadrature nodes and weights
    /// on `[-1, 1]`.
    #[must_use]
    pub fn gauss_legendre_rule(
        &self,
        n: usize,
    ) -> Arc<(Vec<f64>, Vec<f64>)> {
        if n > 4096 {
            return Arc::new(build_gauss_legendre_rule(n));
        }
        self.get_or_init_quadrature(QuadratureRuleKind::GaussLegendre, n, || {
            Arc::new(build_gauss_legendre_rule(n))
        })
    }

    /// Computes and caches Gauss-Kronrod quadrature nodes and weights on `[-1, 1]`.
    ///
    /// Uses at least `n` points: `0` returns empty vectors, `1` returns `([0.0], [2.0])`,
    /// odd `n = 2m + 1` returns the `(2m + 1)`-point rule (`n = 15` is G7K15 and
    /// `n = 21` is G10K21), and positive even `n` rounds up to `n + 1` points.
    #[must_use]
    pub fn gauss_kronrod_rule(
        &self,
        n: usize,
    ) -> Arc<(Vec<f64>, Vec<f64>)> {
        let points = match n {
            | 0 => 0,
            | 1 => 1,
            | val if val.is_multiple_of(2) => val.saturating_add(1),
            | val => val,
        };
        if points > 1025 {
            return Arc::new(build_gauss_kronrod_rule(points));
        }
        self.get_or_init_quadrature(QuadratureRuleKind::GaussKronrod, points, || {
            Arc::new(build_gauss_kronrod_rule(points))
        })
    }

    /// Computes and caches Clenshaw-Curtis quadrature nodes and weights on `[-1, 1]`
    /// with `n + 1` points (`n` is rounded up to an even integer, at least 2).
    #[must_use]
    pub fn clenshaw_curtis_rule(
        &self,
        n: usize,
    ) -> Arc<(Vec<f64>, Vec<f64>)> {
        let normalized = n.max(2).div_ceil(2).saturating_mul(2);
        if normalized > 4096 {
            return Arc::new(build_clenshaw_curtis_rule(normalized));
        }
        self.get_or_init_quadrature(QuadratureRuleKind::ClenshawCurtis, normalized, || {
            Arc::new(build_clenshaw_curtis_rule(normalized))
        })
    }

    /// Computes and caches `n` Chebyshev-Gauss nodes `cos((2k + 1) * pi / (2n))`
    /// on `[-1, 1]` for `k = 0..n`.
    #[must_use]
    pub fn chebyshev_nodes(
        &self,
        n: usize,
    ) -> Arc<Vec<f64>> {
        if n > 65_536 {
            return Arc::new(build_chebyshev_gauss_nodes(n));
        }
        self.chebyshev_nodes_cache
            .get_or_insert_with((false, n), || Arc::new(build_chebyshev_gauss_nodes(n)))
    }

    /// Computes and caches `n + 1` Chebyshev-Lobatto extrema `cos(k * pi / n)`
    /// on `[-1, 1]` for `k = 0..=n` (matching [`crate::kernels::approx::chebyshev_points`]).
    #[must_use]
    pub fn chebyshev_lobatto_nodes(
        &self,
        n: usize,
    ) -> Arc<Vec<f64>> {
        if n > 65_536 {
            return Arc::new(build_chebyshev_lobatto_nodes(n));
        }
        self.chebyshev_nodes_cache
            .get_or_insert_with((true, n), || Arc::new(build_chebyshev_lobatto_nodes(n)))
    }

    /// Computes and caches forward half-spectrum twiddle factors
    /// `exp(-i * pi * k / half_len)` for `k = 1..half_len`.
    #[must_use]
    pub fn real_fft_twiddles(
        &self,
        half_len: usize,
    ) -> Arc<[Complex64]> {
        if half_len > 32_768 {
            return build_real_twiddles(half_len);
        }
        self.real_twiddles_cache
            .get_or_insert_with(half_len, || build_real_twiddles(half_len))
    }

    /// Computes and caches `n` complex roots of unity `exp(-2 * i * pi * k / n)`
    /// for `k = 0..n`.
    #[must_use]
    pub fn unit_roots(
        &self,
        n: usize,
    ) -> Arc<[Complex64]> {
        if n > 65_536 {
            return build_unit_roots(n);
        }
        self.unit_roots_cache
            .get_or_insert_with(n, || build_unit_roots(n))
    }

    /// Computes and caches factorial values `[0!, 1!, ..., max_n!]` as `f64`
    /// (entries above `170!` saturate to `f64::INFINITY`).
    #[must_use]
    pub fn factorial_table(
        &self,
        max_n: usize,
    ) -> Arc<Vec<f64>> {
        if max_n > 4096 {
            return Arc::new(build_factorial_table(max_n));
        }
        self.factorial_tables
            .get_or_insert_with(max_n, || Arc::new(build_factorial_table(max_n)))
    }

    /// Computes and caches natural log-factorial values `[ln(0!), ..., ln(max_n!)]`.
    #[must_use]
    pub fn log_factorial_table(
        &self,
        max_n: usize,
    ) -> Arc<Vec<f64>> {
        if max_n > 16_384 {
            return Arc::new(build_log_factorial_table(max_n));
        }
        self.log_factorial_tables
            .get_or_insert_with(max_n, || Arc::new(build_log_factorial_table(max_n)))
    }

    /// Computes and caches row `n` of Pascal's triangle `[C(n, 0), ..., C(n, n)]`
    /// as `f64`.
    #[must_use]
    pub fn binomial_row(
        &self,
        n: usize,
    ) -> Arc<Vec<f64>> {
        if n > 4096 {
            return Arc::new(build_binomial_row(n));
        }
        self.binomial_rows
            .get_or_insert_with(n, || Arc::new(build_binomial_row(n)))
    }

    /// Computes and caches a packed lower-triangular table of Stirling numbers
    /// of the second kind `S(n, k)` for `0 <= k <= n <= max_n`.
    ///
    /// Entry `(n, k)` is located at `n * (n + 1) / 2 + k`.
    #[must_use]
    pub fn stirling_second_table(
        &self,
        max_n: usize,
    ) -> Arc<Vec<f64>> {
        if max_n > 256 {
            return Arc::new(build_stirling_second_table(max_n));
        }
        self.stirling_tables
            .get_or_insert_with(max_n, || Arc::new(build_stirling_second_table(max_n)))
    }

    /// Computes and caches harmonic numbers `[H_0 = 0.0, H_1 = 1.0, ..., H_{max_n}]`
    /// using Neumaier compensated summation.
    #[must_use]
    pub fn harmonic_table(
        &self,
        max_n: usize,
    ) -> Arc<Vec<f64>> {
        if max_n > 65_536 {
            return Arc::new(build_harmonic_table(max_n));
        }
        self.harmonic_tables
            .get_or_insert_with(max_n, || Arc::new(build_harmonic_table(max_n)))
    }

    /// Computes and caches all prime numbers `<= limit` via the Sieve of Eratosthenes.
    ///
    /// Returns an empty vector if `limit < 2` or if the sieve allocation cannot
    /// be reserved.
    #[must_use]
    pub fn prime_sieve_cache(
        &self,
        limit: usize,
    ) -> Arc<Vec<usize>> {
        if limit > 1_000_000 {
            return Arc::new(build_prime_sieve(limit));
        }
        self.prime_sieves
            .get_or_insert_with(limit, || Arc::new(build_prime_sieve(limit)))
    }

    /// Computes and caches a symmetric Reverse Cuthill-McKee ordering for the
    /// CSR sparsity pattern `(row_ptrs, col_indices)`.
    #[must_use]
    pub fn sparse_ordering_plan(
        &self,
        row_ptrs: &[usize],
        col_indices: &[usize],
    ) -> Arc<SparseOrderingPlan> {
        self.get_or_build_sparse_ordering(
            SparseOrderingKind::ReverseCuthillMckee,
            row_ptrs,
            col_indices,
            || build_reverse_cuthill_mckee(row_ptrs, col_indices),
        )
    }

    /// Looks up or builds a cached [`SparseOrderingPlan`] for `(kind, row_ptrs, col_indices)`.
    ///
    /// `build` returns the forward permutation `perm[new] = old` over `0..dim`.
    pub fn get_or_build_sparse_ordering(
        &self,
        kind: SparseOrderingKind,
        row_ptrs: &[usize],
        col_indices: &[usize],
        build: impl FnOnce() -> Vec<usize>,
    ) -> Arc<SparseOrderingPlan> {
        let key = sparse_pattern_digest(kind, row_ptrs, col_indices);
        let dim = key.dim;
        self.sparse_orderings.get_or_insert_with(key, || {
            let perm = build();
            let inv_perm = invert_permutation(&perm, dim);
            Arc::new(SparseOrderingPlan { dim, perm, inv_perm })
        })
    }

    /// Returns a named, bounded LRU [`PlanCache<K, V>`] hosted inside this pool's
    /// [`PlanRegistry`].
    #[must_use]
    pub fn named_plan_cache<K, V>(
        &self,
        name: &str,
        capacity: usize,
    ) -> Arc<PlanCache<K, V>>
    where
        K: Eq + Send + Sync + 'static,
        V: Send + Sync + ?Sized + 'static,
    {
        self.registry.cache(name, capacity)
    }

    /// Returns a reference to the extensible [`PlanRegistry`] inside this pool.
    #[must_use]
    pub const fn registry(&self) -> &PlanRegistry {
        &self.registry
    }

    /// Returns aggregated telemetry across all scratch pools, built-in plan
    /// caches, and registered named caches.
    #[must_use]
    pub fn stats(&self) -> ReusePoolStats {
        let mut scratch = CacheStats::default();
        scratch.accumulate(self.f64_pool.stats());
        scratch.accumulate(self.f32_pool.stats());
        scratch.accumulate(self.complex64_pool.stats());
        scratch.accumulate(self.complex64_large_pool.stats());
        scratch.accumulate(self.usize_pool.stats());
        scratch.accumulate(self.i64_pool.stats());
        scratch.accumulate(self.bytes_pool.stats());
        if self.is_global() {
            scratch.retained_bytes = scratch
                .retained_bytes
                .saturating_add(self.tls_bytes.load(Ordering::Relaxed));
        }

        let mut plans = CacheStats::default();
        for cache in self.builtin_managed_caches() {
            plans.accumulate(cache.managed_stats());
        }
        plans.accumulate(self.registry.stats());

        ReusePoolStats {
            hits: scratch.hits.saturating_add(plans.hits),
            misses: scratch.misses.saturating_add(plans.misses),
            evictions: scratch.evictions.saturating_add(plans.evictions),
            retained_scratch_bytes: scratch.retained_bytes,
            retained_scratch_buffers: scratch.entries,
            cached_plans: plans.entries,
            scratch,
            plans,
        }
    }

    /// Resets cumulative hit, miss, and eviction counters across all pools and caches.
    pub fn reset_stats(&self) {
        self.f64_pool.reset_stats();
        self.f32_pool.reset_stats();
        self.complex64_pool.reset_stats();
        self.complex64_large_pool.reset_stats();
        self.usize_pool.reset_stats();
        self.i64_pool.reset_stats();
        self.bytes_pool.reset_stats();
        for cache in self.builtin_managed_caches() {
            cache.managed_reset_stats();
        }
        self.registry.reset_stats();
    }

    /// Clears all idle scratch buffers (including calling-thread global TLS
    /// slots and bumping the global TLS epoch) and all cached plans and tables.
    pub fn clear_all(&self) {
        self.f64_pool.clear();
        self.f32_pool.clear();
        self.complex64_pool.clear();
        self.complex64_large_pool.clear();
        self.usize_pool.clear();
        self.i64_pool.clear();
        self.bytes_pool.clear();
        if self.is_global() {
            self.tls_epoch.fetch_add(1, Ordering::Release);
            clear_current_thread_tls();
        }
        for cache in self.builtin_managed_caches() {
            cache.managed_clear();
        }
        self.registry.clear();
    }

    /// Evicts idle scratch buffers and unreferenced plan entries across all
    /// pools and caches.
    pub fn trim_idle(&self) {
        self.f64_pool.trim_idle();
        self.f32_pool.trim_idle();
        self.complex64_pool.trim_idle();
        self.complex64_large_pool.trim_idle();
        self.usize_pool.trim_idle();
        self.i64_pool.trim_idle();
        self.bytes_pool.trim_idle();
        if self.is_global() {
            self.tls_epoch.fetch_add(1, Ordering::Release);
            clear_current_thread_tls();
        }
        for cache in self.builtin_managed_caches() {
            cache.managed_trim_idle();
        }
        self.registry.trim_idle();
    }

    /// Trims idle scratch buffers (evicting largest buffers first) until total
    /// retained scratch bytes is at most `max_bytes` (or no evictable shared
    /// buffers remain). When `max_bytes == 0`, also trims idle plan entries.
    pub fn trim_to_budget(
        &self,
        max_bytes: usize,
    ) {
        if self.is_global() && self.total_retained_scratch_bytes() > max_bytes {
            clear_current_thread_tls();
            if max_bytes == 0 {
                self.tls_epoch.fetch_add(1, Ordering::Release);
            }
        }
        while self.total_retained_scratch_bytes() > max_bytes {
            let candidates = [
                (0_usize, self.complex64_large_pool.retained_bytes()),
                (1_usize, self.complex64_pool.retained_bytes()),
                (2_usize, self.f64_pool.retained_bytes()),
                (3_usize, self.f32_pool.retained_bytes()),
                (4_usize, self.usize_pool.retained_bytes()),
                (5_usize, self.i64_pool.retained_bytes()),
                (6_usize, self.bytes_pool.retained_bytes()),
            ];
            let Some((pool_idx, _)) = candidates
                .into_iter()
                .filter(|(_, bytes)| *bytes > 0)
                .max_by_key(|(_, bytes)| *bytes)
            else {
                break;
            };
            let freed = match pool_idx {
                | 0 => self.complex64_large_pool.evict_largest_one(),
                | 1 => self.complex64_pool.evict_largest_one(),
                | 2 => self.f64_pool.evict_largest_one(),
                | 3 => self.f32_pool.evict_largest_one(),
                | 4 => self.usize_pool.evict_largest_one(),
                | 5 => self.i64_pool.evict_largest_one(),
                | _ => self.bytes_pool.evict_largest_one(),
            };
            if freed == 0 {
                break;
            }
        }
        if max_bytes == 0 {
            for cache in self.builtin_managed_caches() {
                cache.managed_trim_idle();
            }
            self.registry.trim_idle();
        }
    }

    fn total_retained_scratch_bytes(&self) -> usize {
        let mut total = self
            .f64_pool
            .retained_bytes()
            .saturating_add(self.f32_pool.retained_bytes())
            .saturating_add(self.complex64_pool.retained_bytes())
            .saturating_add(self.complex64_large_pool.retained_bytes())
            .saturating_add(self.usize_pool.retained_bytes())
            .saturating_add(self.i64_pool.retained_bytes())
            .saturating_add(self.bytes_pool.retained_bytes());
        if self.is_global() {
            total = total.saturating_add(self.tls_bytes.load(Ordering::Relaxed));
        }
        total
    }

    fn builtin_managed_caches(&self) -> [&dyn ManagedPlanCache; 13] {
        [
            &self.fft_small_plans,
            &self.fft_large_plans,
            &self.quadrature_rules,
            &self.chebyshev_nodes_cache,
            &self.real_twiddles_cache,
            &self.unit_roots_cache,
            &self.factorial_tables,
            &self.log_factorial_tables,
            &self.binomial_rows,
            &self.stirling_tables,
            &self.harmonic_tables,
            &self.prime_sieves,
            &self.sparse_orderings,
        ]
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
    global_pool().fft_plan(length, direction)
}

/// Checks out a resizable scratch vector for complex `f64` transforms.
#[must_use]
pub fn fft_scratch() -> Scratch<'static, Complex64> {
    global_pool().complex64_pool.checkout()
}

/// Checks out an empty complex `f64` buffer with the requested capacity.
///
/// Buffers over four MiB use their own retention pool, which retains up to
/// four idle buffers of at most sixteen MiB each.
///
/// # Panics
/// Panics if `minimum` exceeds the maximum capacity of a vector.
#[must_use]
pub fn fft_scratch_for(minimum: usize) -> Scratch<'static, Complex64> {
    global_pool().scratch_complex64(minimum)
}

/// Checks out an empty `f64` buffer with the requested capacity. Kernels
/// share this pool for intermediate values and return only their own results.
///
/// # Panics
/// Panics if `minimum` exceeds the maximum capacity of a vector.
#[must_use]
pub fn f64_scratch(minimum: usize) -> Scratch<'static, f64> {
    global_pool().scratch_f64(minimum)
}

pub(crate) fn build_gauss_legendre_rule(n: usize) -> (Vec<f64>, Vec<f64>) {
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let nf = n as f64;
    let half = n.div_ceil(2);
    let pairs = crate::kernels::par::map_indexed(half, 64, |i| {
        let mut z = if !n.is_multiple_of(2) && i == n / 2 {
            0.0
        } else {
            (PI * (i as f64 + 0.75) / (nf + 0.5)).cos()
        };
        if z != 0.0 {
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, z);
                for k in 2..=n {
                    let kf = k as f64;
                    let p2 = ((2.0 * kf - 1.0) * z * p1 - (kf - 1.0) * p0) / kf;
                    p0 = p1;
                    p1 = p2;
                }
                let dp = nf * (z * p1 - p0) / (z * z - 1.0);
                let dz = p1 / dp;
                z -= dz;
                if dz.abs() <= 1e-15 * z.abs().max(1.0) {
                    break;
                }
            }
        }
        let (mut p0, mut p1) = (1.0, z);
        for k in 2..=n {
            let kf = k as f64;
            let p2 = ((2.0 * kf - 1.0) * z * p1 - (kf - 1.0) * p0) / kf;
            p0 = p1;
            p1 = p2;
        }
        let dp = nf * (z * p1 - p0) / (z * z - 1.0);
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        (z, wi)
    });
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for (i, (z, wi)) in pairs.into_iter().enumerate() {
        let j = n.saturating_sub(1).saturating_sub(i);
        if let Some(slot) = x.get_mut(i) {
            *slot = if z == 0.0 { 0.0 } else { -z };
        }
        if let Some(slot) = x.get_mut(j) {
            *slot = z;
        }
        if let Some(slot) = w.get_mut(i) {
            *slot = wi;
        }
        if let Some(slot) = w.get_mut(j) {
            *slot = wi;
        }
    }
    (x, w)
}

#[allow(clippy::single_call_fn)] // isolated Neumann-Adams triple Legendre integral helper
fn adams_triple_legendre(
    a: usize,
    b: usize,
    c: usize,
) -> f64 {
    let sum = a.saturating_add(b).saturating_add(c);
    if !sum.is_multiple_of(2) {
        return 0.0;
    }
    let s = sum / 2;
    if s < a || s < b || s < c {
        return 0.0;
    }
    let p = s.saturating_sub(a);
    let q = s.saturating_sub(b);
    let r = s.saturating_sub(c);
    let mut val = 1.0_f64;
    let mut j = 1_usize;
    for bound in [p, q, r] {
        for i in 1..=bound {
            let iflt = i as f64;
            let jflt = j as f64;
            val *= ((2.0 * iflt - 1.0) * jflt) / (iflt * (2.0 * jflt - 1.0));
            j = j.saturating_add(1);
        }
    }
    (2.0 / (2.0 * (s as f64) + 1.0)) * val
}

fn eval_legendre_and_stieltjes(
    m: usize,
    coeffs: &[f64],
    x: f64,
    poly_buf: &mut [f64],
    dpoly_buf: &mut [f64],
) -> (f64, f64, f64, f64) {
    if let Some(p0) = poly_buf.get_mut(0) {
        *p0 = 1.0;
    }
    if let Some(dp0) = dpoly_buf.get_mut(0) {
        *dp0 = 0.0;
    }
    if let Some(p1) = poly_buf.get_mut(1) {
        *p1 = x;
    }
    if let Some(dp1) = dpoly_buf.get_mut(1) {
        *dp1 = 1.0;
    }
    for k in 1..=m {
        let kf = k as f64;
        let pk = poly_buf.get(k).copied().unwrap_or(0.0);
        let pk_prev = poly_buf.get(k.saturating_sub(1)).copied().unwrap_or(0.0);
        let dpk = dpoly_buf.get(k).copied().unwrap_or(0.0);
        let dpk_prev = dpoly_buf.get(k.saturating_sub(1)).copied().unwrap_or(0.0);
        let next_p = ((2.0 * kf + 1.0) * x * pk - kf * pk_prev) / (kf + 1.0);
        let next_dp = ((2.0 * kf + 1.0) * (pk + x * dpk) - kf * dpk_prev) / (kf + 1.0);
        let next_idx = k.saturating_add(1);
        if let Some(slot) = poly_buf.get_mut(next_idx) {
            *slot = next_p;
        }
        if let Some(slot) = dpoly_buf.get_mut(next_idx) {
            *slot = next_dp;
        }
    }
    let mut pi_val = 0.0;
    let mut dpi_val = 0.0;
    for (j, &cj) in coeffs.iter().enumerate() {
        let idx = m.saturating_add(1).saturating_sub(j.saturating_mul(2));
        pi_val += cj * poly_buf.get(idx).copied().unwrap_or(0.0);
        dpi_val += cj * dpoly_buf.get(idx).copied().unwrap_or(0.0);
    }
    (
        poly_buf.get(m).copied().unwrap_or(0.0),
        dpoly_buf.get(m).copied().unwrap_or(0.0),
        pi_val,
        dpi_val,
    )
}

fn build_gauss_kronrod_rule(points: usize) -> (Vec<f64>, Vec<f64>) {
    if points == 0 {
        return (Vec::new(), Vec::new());
    }
    if points == 1 {
        return (vec![0.0], vec![2.0]);
    }
    let m = points / 2;
    let total = m.saturating_mul(2).saturating_add(1);
    let m_max = m.saturating_add(1) / 2;
    let mut coeffs = vec![0.0; m_max.saturating_add(1)];
    if let Some(c0) = coeffs.get_mut(0) {
        *c0 = 1.0;
    }
    for k in 1..=m_max {
        let l = k.saturating_mul(2).saturating_sub(1);
        let mut rhs = 0.0;
        for (j, &cj) in coeffs.iter().take(k).enumerate() {
            let deg = m.saturating_add(1).saturating_sub(j.saturating_mul(2));
            rhs += cj * adams_triple_legendre(m, deg, l);
        }
        let diag_deg = m.saturating_add(1).saturating_sub(k.saturating_mul(2));
        let diag = adams_triple_legendre(m, diag_deg, l);
        if let Some(ck) = coeffs.get_mut(k) {
            *ck = if diag != 0.0 {
                -rhs / diag
            } else {
                0.0
            };
        }
    }

    let (gx, gw) = build_gauss_legendre_rule(m);
    let mut brackets = Vec::with_capacity(m.saturating_add(2));
    brackets.push(-1.0);
    brackets.extend_from_slice(&gx);
    brackets.push(1.0);

    let mut nodes = vec![0.0; total];
    let mut weights = vec![0.0; total];
    let mut poly_buf = vec![0.0; m.saturating_add(2)];
    let mut dpoly_buf = vec![0.0; m.saturating_add(2)];
    let mf_plus_1 = (m as f64) + 1.0;

    for idx in 0..=(m / 2) {
        let mut z = if m.is_multiple_of(2) && idx == m / 2 {
            0.0
        } else {
            let mut lo = brackets.get(idx).copied().unwrap_or(-1.0);
            let mut hi = brackets.get(idx.saturating_add(1)).copied().unwrap_or(1.0);
            let mut flo =
                eval_legendre_and_stieltjes(m, &coeffs, lo, &mut poly_buf, &mut dpoly_buf).2;
            let mut root = f64::midpoint(lo, hi);
            for _ in 0..80 {
                let (_, _, val, dval) =
                    eval_legendre_and_stieltjes(m, &coeffs, root, &mut poly_buf, &mut dpoly_buf);
                if (val > 0.0) == (flo > 0.0) {
                    lo = root;
                    flo = val;
                } else {
                    hi = root;
                }
                if val.abs() <= 1e-16 {
                    break;
                }
                if dval != 0.0 {
                    let cand = root - val / dval;
                    if cand > lo && cand < hi {
                        let step = (cand - root).abs();
                        root = cand;
                        if step <= 1e-15 * root.abs().max(1.0) {
                            break;
                        }
                        continue;
                    }
                }
                root = f64::midpoint(lo, hi);
                if hi - lo <= 1e-15 * root.abs().max(1.0) {
                    break;
                }
            }
            root
        };
        if z == 0.0 {
            z = 0.0;
        }
        let (pn, _, _, dpi) =
            eval_legendre_and_stieltjes(m, &coeffs, z, &mut poly_buf, &mut dpoly_buf);
        let wi = 2.0 / (mf_plus_1 * pn * dpi);
        let left_pos = idx.saturating_mul(2);
        let right_pos = total.saturating_sub(1).saturating_sub(left_pos);
        if let Some(slot) = nodes.get_mut(left_pos) {
            *slot = z;
        }
        if let Some(slot) = weights.get_mut(left_pos) {
            *slot = wi;
        }
        if let Some(slot) = nodes.get_mut(right_pos) {
            *slot = if z == 0.0 { 0.0 } else { -z };
        }
        if let Some(slot) = weights.get_mut(right_pos) {
            *slot = wi;
        }
    }

    for idx in 0..m.div_ceil(2) {
        let y = gx.get(idx).copied().unwrap_or(0.0);
        let gwi = gw.get(idx).copied().unwrap_or(0.0);
        let (_, dpn, pi_val, _) =
            eval_legendre_and_stieltjes(m, &coeffs, y, &mut poly_buf, &mut dpoly_buf);
        let wi = gwi + 2.0 / (mf_plus_1 * dpn * pi_val);
        let left_pos = idx.saturating_mul(2).saturating_add(1);
        let right_pos = total.saturating_sub(1).saturating_sub(left_pos);
        if let Some(slot) = nodes.get_mut(left_pos) {
            *slot = if y == 0.0 { 0.0 } else { y };
        }
        if let Some(slot) = weights.get_mut(left_pos) {
            *slot = wi;
        }
        if let Some(slot) = nodes.get_mut(right_pos) {
            *slot = if y == 0.0 { 0.0 } else { -y };
        }
        if let Some(slot) = weights.get_mut(right_pos) {
            *slot = wi;
        }
    }

    (nodes, weights)
}

pub(crate) fn build_clenshaw_curtis_rule(n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n.max(2).div_ceil(2).saturating_mul(2);
    let len = n.saturating_add(1);
    let nf = n as f64;
    let half_n = n / 2;
    let mut x = vec![0.0; len];
    let mut w = vec![0.0; len];
    if n >= 64 {
        // Waldvogel's O(n log n) FFT algorithm for Clenshaw-Curtis weights.
        for k in 0..=half_n {
            let xk = if k == 0 {
                1.0
            } else if k == half_n {
                0.0
            } else {
                (PI * (k as f64) / nf).cos()
            };
            let mirror = n.saturating_sub(k);
            if let Some(slot) = x.get_mut(k) {
                *slot = xk;
            }
            if let Some(slot) = x.get_mut(mirror) {
                *slot = if xk == 0.0 { 0.0 } else { -xk };
            }
        }
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(n);
        let mut buf = vec![Complex64::new(0.0, 0.0); n];
        if let Some(slot) = buf.get_mut(0) {
            *slot = Complex64::new(1.0, 0.0);
        }
        for j in 1..half_n {
            let jf = j as f64;
            let val = -1.0 / (4.0 * jf * jf - 1.0);
            if let Some(slot) = buf.get_mut(j) {
                *slot = Complex64::new(val, 0.0);
            }
            if let Some(slot) = buf.get_mut(n.saturating_sub(j)) {
                *slot = Complex64::new(val, 0.0);
            }
        }
        if let Some(slot) = buf.get_mut(half_n) {
            *slot = Complex64::new(-1.0 / (nf * nf - 1.0), 0.0);
        }
        let mut scratch = vec![Complex64::new(0.0, 0.0); fft.get_inplace_scratch_len()];
        fft.process_with_scratch(&mut buf, &mut scratch);
        for k in 0..=half_n {
            let c = if k == 0 || k == n { 1.0 } else { 2.0 };
            let wk = c * buf.get(k).map_or(0.0, |v| v.re) / nf;
            let mirror = n.saturating_sub(k);
            if let Some(slot) = w.get_mut(k) {
                *slot = wk;
            }
            if let Some(slot) = w.get_mut(mirror) {
                *slot = wk;
            }
        }
        return (x, w);
    }
    for k in 0..=half_n {
        let th = PI * (k as f64) / nf;
        let xk = if k == 0 {
            1.0
        } else if k == half_n {
            0.0
        } else {
            th.cos()
        };
        let mut s = 0.0;
        for j in 1..=half_n {
            let b = if j == half_n { 1.0 } else { 2.0 };
            let jf = j as f64;
            s += b * (2.0 * jf * th).cos() / (4.0 * jf * jf - 1.0);
        }
        let c = if k == 0 || k == n { 1.0 } else { 2.0 };
        let wk = c / nf * (1.0 - s);
        let mirror = n.saturating_sub(k);
        if let Some(slot) = x.get_mut(k) {
            *slot = xk;
        }
        if let Some(slot) = x.get_mut(mirror) {
            *slot = if xk == 0.0 { 0.0 } else { -xk };
        }
        if let Some(slot) = w.get_mut(k) {
            *slot = wk;
        }
        if let Some(slot) = w.get_mut(mirror) {
            *slot = wk;
        }
    }
    (x, w)
}

fn build_chebyshev_gauss_nodes(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let denom = 2.0 * (n as f64);
    let odd_n = !n.is_multiple_of(2);
    let mid = n / 2;
    (0..n)
        .map(|k| {
            if odd_n && k == mid {
                0.0
            } else {
                let num = 2.0 * (k as f64) + 1.0;
                (PI * num / denom).cos()
            }
        })
        .collect()
}

fn build_chebyshev_lobatto_nodes(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![0.0];
    }
    let nf = n as f64;
    let even_n = n.is_multiple_of(2);
    let mid = n / 2;
    (0..=n)
        .map(|k| {
            if k == 0 {
                1.0
            } else if k == n {
                -1.0
            } else if even_n && k == mid {
                0.0
            } else {
                (PI * (k as f64) / nf).cos()
            }
        })
        .collect()
}

fn build_real_twiddles(half_len: usize) -> Arc<[Complex64]> {
    if half_len <= 1 {
        return Arc::from([]);
    }
    let scale = -PI / (half_len as f64);
    (1..half_len)
        .map(|bin| Complex64::from_polar(1.0, scale * (bin as f64)))
        .collect()
}

fn build_unit_roots(n: usize) -> Arc<[Complex64]> {
    if n == 0 {
        return Arc::from([]);
    }
    let scale = -2.0 * PI / (n as f64);
    (0..n)
        .map(|bin| Complex64::from_polar(1.0, scale * (bin as f64)))
        .collect()
}

fn build_factorial_table(max_n: usize) -> Vec<f64> {
    let mut table = Vec::with_capacity(max_n.saturating_add(1).min(4097));
    table.push(1.0);
    let mut acc = 1.0_f64;
    for k in 1..=max_n {
        acc = if acc.is_finite() {
            acc * (k as f64)
        } else {
            f64::INFINITY
        };
        table.push(acc);
    }
    table
}

fn build_log_factorial_table(max_n: usize) -> Vec<f64> {
    let mut table = Vec::with_capacity(max_n.saturating_add(1).min(16_385));
    table.push(0.0);
    let mut sum = 0.0_f64;
    let mut comp = 0.0_f64;
    for k in 1..=max_n {
        let term = (k as f64).ln();
        let temp = sum + term;
        if sum.abs() >= term.abs() {
            comp += (sum - temp) + term;
        } else {
            comp += (term - temp) + sum;
        }
        sum = temp;
        table.push(sum + comp);
    }
    table
}

fn build_binomial_row(n: usize) -> Vec<f64> {
    let len = n.saturating_add(1);
    let mut row = vec![0.0; len];
    if let Some(first) = row.get_mut(0) {
        *first = 1.0;
    }
    if let Some(last) = row.get_mut(n) {
        *last = 1.0;
    }
    let mut acc = 1.0_f64;
    for k in 1..=(n / 2) {
        let numer = n.saturating_sub(k).saturating_add(1) as f64;
        let denom = k as f64;
        acc = if acc.is_finite() {
            (acc * numer) / denom
        } else {
            f64::INFINITY
        };
        let mirror = n.saturating_sub(k);
        if let Some(slot) = row.get_mut(k) {
            *slot = acc;
        }
        if let Some(slot) = row.get_mut(mirror) {
            *slot = acc;
        }
    }
    row
}

fn build_stirling_second_table(max_n: usize) -> Vec<f64> {
    let rows = max_n.saturating_add(1);
    let total = rows.saturating_mul(rows.saturating_add(1)) / 2;
    let mut table = vec![0.0; total];
    if let Some(s00) = table.get_mut(0) {
        *s00 = 1.0;
    }
    for n in 1..=max_n {
        let row_start = n.saturating_mul(n.saturating_add(1)) / 2;
        let prev_start = n.saturating_sub(1).saturating_mul(n) / 2;
        if let Some(diag) = table.get_mut(row_start.saturating_add(n)) {
            *diag = 1.0;
        }
        for k in 1..n {
            let same = table
                .get(prev_start.saturating_add(k))
                .copied()
                .unwrap_or(0.0);
            let left = table
                .get(prev_start.saturating_add(k.saturating_sub(1)))
                .copied()
                .unwrap_or(0.0);
            if let Some(slot) = table.get_mut(row_start.saturating_add(k)) {
                *slot = (k as f64) * same + left;
            }
        }
    }
    table
}

fn build_harmonic_table(max_n: usize) -> Vec<f64> {
    let mut table = Vec::with_capacity(max_n.saturating_add(1).min(65_537));
    table.push(0.0);
    let mut sum = 0.0_f64;
    let mut comp = 0.0_f64;
    for k in 1..=max_n {
        let term = 1.0 / (k as f64);
        let temp = sum + term;
        if sum.abs() >= term.abs() {
            comp += (sum - temp) + term;
        } else {
            comp += (term - temp) + sum;
        }
        sum = temp;
        table.push(sum + comp);
    }
    table
}

fn build_prime_sieve(limit: usize) -> Vec<usize> {
    if limit < 2 {
        return Vec::new();
    }
    let Some(len) = limit.checked_add(1) else {
        return Vec::new();
    };
    let mut is_prime: Vec<bool> = Vec::new();
    if is_prime.try_reserve_exact(len).is_err() {
        return Vec::new();
    }
    is_prime.resize(len, true);
    if let Some(s0) = is_prime.get_mut(0) {
        *s0 = false;
    }
    if let Some(s1) = is_prime.get_mut(1) {
        *s1 = false;
    }
    let mut p = 2_usize;
    while p.saturating_mul(p) <= limit {
        if is_prime.get(p).copied().unwrap_or(false) {
            let start = p.saturating_mul(p);
            for mult in (start..=limit).step_by(p) {
                if let Some(flag) = is_prime.get_mut(mult) {
                    *flag = false;
                }
            }
        }
        p = p.saturating_add(1);
    }
    is_prime
        .into_iter()
        .enumerate()
        .filter(|&(_, prime)| prime)
        .map(|(value, _)| value)
        .collect()
}

#[allow(clippy::single_call_fn)] // isolates blake3 pattern hashing for sparse ordering cache keys
fn sparse_pattern_digest(
    kind: SparseOrderingKind,
    row_ptrs: &[usize],
    col_indices: &[usize],
) -> SparseOrderingKey {
    let mut hasher = Hasher::new();
    let kind_tag: u8 = match kind {
        | SparseOrderingKind::ReverseCuthillMckee => 1,
        | SparseOrderingKind::ApproximateMinimumDegree => 2,
        | SparseOrderingKind::NestedDissection => 3,
        | SparseOrderingKind::Custom(_) => 4,
    };
    hasher.update(&[kind_tag]);
    if let SparseOrderingKind::Custom(id) = kind {
        hasher.update(&id.to_le_bytes());
    }
    hasher.update(&(row_ptrs.len() as u64).to_le_bytes());
    for &ptr in row_ptrs {
        hasher.update(&(ptr as u64).to_le_bytes());
    }
    hasher.update(&(col_indices.len() as u64).to_le_bytes());
    for &col in col_indices {
        hasher.update(&(col as u64).to_le_bytes());
    }
    SparseOrderingKey {
        kind,
        dim: row_ptrs.len().saturating_sub(1),
        nnz: col_indices.len(),
        digest: *hasher.finalize().as_bytes(),
    }
}

#[allow(clippy::single_call_fn)] // dedicated CSR Reverse Cuthill-McKee ordering builder
fn build_reverse_cuthill_mckee(
    row_ptrs: &[usize],
    col_indices: &[usize],
) -> Vec<usize> {
    let dim = row_ptrs.len().saturating_sub(1);
    if dim == 0 {
        return Vec::new();
    }
    let mut adj = vec![Vec::new(); dim];
    for u in 0..dim {
        let start = row_ptrs.get(u).copied().unwrap_or(col_indices.len());
        let end = row_ptrs
            .get(u.saturating_add(1))
            .copied()
            .unwrap_or(col_indices.len())
            .min(col_indices.len());
        if start >= end {
            continue;
        }
        if let Some(cols) = col_indices.get(start..end) {
            for &v in cols {
                if v < dim && v != u {
                    if let Some(row_u) = adj.get_mut(u) {
                        row_u.push(v);
                    }
                    if let Some(row_v) = adj.get_mut(v) {
                        row_v.push(u);
                    }
                }
            }
        }
    }
    for neighbors in &mut adj {
        neighbors.sort_unstable();
        neighbors.dedup();
    }
    let degrees: Vec<usize> = adj.iter().map(Vec::len).collect();
    let mut visited = vec![false; dim];
    let mut order = Vec::with_capacity(dim);
    let mut queue = VecDeque::new();
    let mut component = Vec::new();
    let mut level_buf = Vec::new();

    for start in 0..dim {
        if visited.get(start).copied().unwrap_or(true) {
            continue;
        }
        // Pick a minimum-degree vertex in the remaining unvisited vertices as seed,
        // then run one BFS sweep to find a pseudo-peripheral vertex.
        let seed = (start..dim)
            .filter(|&v| !visited.get(v).copied().unwrap_or(true))
            .min_by_key(|&v| degrees.get(v).copied().unwrap_or(0))
            .unwrap_or(start);

        queue.clear();
        level_buf.clear();
        if let Some(flag) = visited.get_mut(seed) {
            *flag = true;
        }
        queue.push_back(seed);
        while let Some(u) = queue.pop_front() {
            level_buf.push(u);
            if let Some(neighbors) = adj.get(u) {
                for &v in neighbors {
                    if !visited.get(v).copied().unwrap_or(true) {
                        if let Some(flag) = visited.get_mut(v) {
                            *flag = true;
                        }
                        queue.push_back(v);
                    }
                }
            }
        }
        let peripheral = level_buf
            .iter()
            .rev()
            .copied()
            .min_by_key(|&v| degrees.get(v).copied().unwrap_or(0))
            .unwrap_or(seed);
        for &v in &level_buf {
            if let Some(flag) = visited.get_mut(v) {
                *flag = false;
            }
        }

        component.clear();
        queue.clear();
        if let Some(flag) = visited.get_mut(peripheral) {
            *flag = true;
        }
        queue.push_back(peripheral);
        while let Some(u) = queue.pop_front() {
            component.push(u);
            if let Some(neighbors) = adj.get(u) {
                let mut unvisited: Vec<usize> = neighbors
                    .iter()
                    .copied()
                    .filter(|&v| !visited.get(v).copied().unwrap_or(true))
                    .collect();
                unvisited.sort_unstable_by_key(|&v| degrees.get(v).copied().unwrap_or(0));
                for v in unvisited {
                    if !visited.get(v).copied().unwrap_or(true) {
                        if let Some(flag) = visited.get_mut(v) {
                            *flag = true;
                        }
                        queue.push_back(v);
                    }
                }
            }
        }
        component.reverse();
        order.extend_from_slice(&component);
    }
    order
}

#[allow(clippy::single_call_fn)] // paired inverse permutation builder for SparseOrderingPlan
fn invert_permutation(
    perm: &[usize],
    dim: usize,
) -> Vec<usize> {
    let mut inv = vec![0; dim];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        if let Some(slot) = inv.get_mut(old_idx) {
            *slot = new_idx;
        }
    }
    inv
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(dead_code)]
    const fn assert_send_sync<T: Send + Sync>() {}

    #[test]
    fn compile_time_send_sync_bounds() {
        assert_send_sync::<ComputeReusePool>();
        assert_send_sync::<PlanRegistry>();
        assert_send_sync::<ScratchPool<f64>>();
        assert_send_sync::<PlanCache<u64, Vec<f64>>>();
        assert_send_sync::<Scratch<'_, f64>>();
        assert_send_sync::<CacheStats>();
        assert_send_sync::<ReusePoolStats>();
        assert_send_sync::<SparseOrderingPlan>();
    }

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
        let plans = PlanCache::new_lru(2);
        let hot = plans.get_or_insert_with(0, || Arc::new(41));
        for cold in 1..=20 {
            assert_eq!(*plans.get_or_insert_with(cold, || Arc::new(cold)), cold);
            let cached = plans.get_or_insert_with(0, || panic!("hot FFT plan was evicted"));
            assert!(Arc::ptr_eq(&hot, &cached));
        }
    }

    #[test]
    fn fft_recency_preserves_eviction_order_when_the_clock_wraps() {
        let plans = PlanCache::new_lru(2);
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
        assert!(std::ptr::eq(
            scratch.pool,
            &global_pool().complex64_large_pool
        ));
        let small = fft_scratch_for(17);
        assert!(std::ptr::eq(small.pool, &global_pool().complex64_pool));
        assert!(!std::ptr::eq(scratch.pool, small.pool));
    }

    #[test]
    fn sharded_scratch_pool_scales_across_8_and_16_threads_with_steady_state_reuse() {
        use std::sync::Barrier;

        for threads in [8_usize, 16_usize] {
            let pool = ScratchPool::<f64>::new(32, 64 * 1024);
            let barrier = Barrier::new(threads);
            let iterations = 60_usize;

            std::thread::scope(|scope| {
                for _ in 0..threads {
                    scope.spawn(|| {
                        barrier.wait();
                        for step in 0..iterations {
                            let mut buf = pool.checkout_for(128);
                            buf.resize(128, step as f64);
                        }
                    });
                }
            });

            let stats = pool.stats();
            let total_ops = (threads.saturating_mul(iterations)) as u64;
            assert!(pool.len() <= pool.max_idle());
            // Shard indices come from one process-wide counter, so threads
            // outside this test can share a shard with one of ours; the buffer
            // a thread gets back is therefore not always the one it returned.
            // What the pool guarantees is that checkouts are served from
            // retained buffers, which is what the counters measure.
            assert!(
                stats.hits >= total_ops.saturating_mul(9) / 10,
                "expected steady-state reuse under {threads} threads, got {stats:?}"
            );
            assert!(pool.retained_bytes() > 0);
            pool.clear();
            assert_eq!(pool.len(), 0);
            assert_eq!(pool.retained_bytes(), 0);
        }
    }

    #[test]
    fn global_pool_tls_fast_path_reuses_multiple_simultaneous_slots() {
        use std::sync::Barrier;

        let barrier = Barrier::new(8);
        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    barrier.wait();
                    let a = global_pool().scratch_f64(64);
                    let b = global_pool().scratch_f64(128);
                    let c = global_pool().scratch_f64(256);
                    let (p0, p1, p2) = (a.as_ptr(), b.as_ptr(), c.as_ptr());
                    drop(a);
                    drop(b);
                    drop(c);

                    let a2 = global_pool().scratch_f64(64);
                    let b2 = global_pool().scratch_f64(128);
                    let c2 = global_pool().scratch_f64(256);
                    assert_eq!(a2.as_ptr(), p0);
                    assert_eq!(b2.as_ptr(), p1);
                    assert_eq!(c2.as_ptr(), p2);
                });
            }
        });
    }

    #[test]
    fn isolated_compute_reuse_pool_accounting_budget_and_clear() {
        const MIN_EXPECTED_BYTES: usize = 64 * 8 + 32 * 4 + 16 * 16 + 24 * 8 + 12 * 8 + 256;
        let p1 = ComputeReusePool::new();
        let p2 = ComputeReusePool::new();

        {
            let mut f64_s = p1.scratch_f64(64);
            let mut f32_s = p1.scratch_f32(32);
            let mut c64_s = p1.scratch_complex64(16);
            let mut u_s = p1.scratch_usize(24);
            let mut i_s = p1.scratch_i64(12);
            let mut b_s = p1.scratch_bytes(256);
            f64_s.push(1.0);
            f32_s.push(2.0);
            c64_s.push(Complex64::new(1.0, -1.0));
            u_s.push(7);
            i_s.push(-3);
            b_s.push(0xAB);
        }
        {
            let mut p2_buf = p2.scratch_f64(128);
            p2_buf.push(3.5);
        }

        let _ = p1.gauss_legendre_rule(8);
        let _ = p1.factorial_table(12);
        let p1_stats = p1.stats();
        let p2_stats_before = p2.stats();

        assert!(p1_stats.retained_scratch_bytes >= MIN_EXPECTED_BYTES);
        assert_eq!(p1_stats.retained_scratch_buffers, 6);
        assert_eq!(p1_stats.cached_plans, 2);
        assert_eq!(p2_stats_before.retained_scratch_buffers, 1);
        assert_eq!(p2_stats_before.cached_plans, 0);

        let half_budget = p1_stats.retained_scratch_bytes / 2;
        p1.trim_to_budget(half_budget);
        let trimmed_stats = p1.stats();
        assert!(trimmed_stats.retained_scratch_bytes <= half_budget);
        assert!(trimmed_stats.evictions >= 1);

        p1.clear_all();
        let cleared_stats = p1.stats();
        assert_eq!(cleared_stats.retained_scratch_bytes, 0);
        assert_eq!(cleared_stats.retained_scratch_buffers, 0);
        assert_eq!(cleared_stats.cached_plans, 0);

        // Pool p2 must remain completely isolated from p1.
        let p2_stats_after = p2.stats();
        assert_eq!(p2_stats_after.retained_scratch_buffers, 1);
        assert_eq!(
            p2_stats_after.retained_scratch_bytes,
            p2_stats_before.retained_scratch_bytes
        );
    }

    #[test]
    fn plan_cache_lru_stats_get_and_trim_idle() {
        let cache = PlanCache::new_lru(2);
        assert!(cache.is_empty());
        assert_eq!(cache.capacity(), 2);
        assert!(cache.get(&1).is_none());
        assert!(!cache.contains(&1));

        let held = cache.get_or_insert_with(1, || Arc::new(100));
        assert!(cache.contains(&1));
        let got = cache.get(&1);
        assert!(got.as_ref().is_some_and(|v| Arc::ptr_eq(v, &held)));
        drop(got);

        let _ = cache.get_or_insert_with(2, || Arc::new(200));
        assert_eq!(cache.len(), 2);

        // Entry 2 has no external Arc holders; trim_idle must evict only entry 2.
        cache.trim_idle();
        assert_eq!(cache.len(), 1);
        assert!(cache.contains(&1));
        assert!(!cache.contains(&2));

        // Now drop `held`; entry 1 becomes idle and is trimmed too.
        drop(held);
        cache.trim_idle();
        assert_eq!(cache.len(), 0);

        let stats = cache.stats();
        assert_eq!(stats.hits, 1);
        assert_eq!(stats.misses, 3);
        assert_eq!(stats.evictions, 2);

        let _ = cache.get_or_insert_with(3, || Arc::new(300));
        cache.clear();
        assert!(cache.is_empty());
        cache.reset_stats();
        assert_eq!(cache.stats(), CacheStats::default());
    }

    #[test]
    fn quadrature_rules_are_accurate_and_cached() {
        let pool = ComputeReusePool::new();

        // Gauss-Legendre edge cases & polynomial exactness up to degree 2n - 1.
        let gl0 = pool.gauss_legendre_rule(0);
        assert!(gl0.0.is_empty() && gl0.1.is_empty());
        let gl1 = pool.gauss_legendre_rule(1);
        assert_eq!(gl1.0, vec![0.0]);
        assert_eq!(gl1.1, vec![2.0]);

        let gl8 = pool.gauss_legendre_rule(8);
        assert!(Arc::ptr_eq(&gl8, &pool.gauss_legendre_rule(8)));
        assert_eq!(gl8.0.len(), 8);
        let weight_sum: f64 = gl8.1.iter().sum();
        assert!((weight_sum - 2.0).abs() < 1e-14);
        // Exact for degree 2*8 - 1 = 15 (check x^14).
        let x14: f64 = gl8
            .0
            .iter()
            .zip(&gl8.1)
            .map(|(&x, &w)| w * x.powi(14))
            .sum();
        assert!((x14 - 2.0 / 15.0).abs() < 1e-14);

        // Gauss-Kronrod edge cases, G7K15, G10K21, and higher order G15K31.
        let gk0 = pool.gauss_kronrod_rule(0);
        assert!(gk0.0.is_empty());
        let gk1 = pool.gauss_kronrod_rule(1);
        assert_eq!(gk1.0, vec![0.0]);
        assert_eq!(gk1.1, vec![2.0]);

        let gk15 = pool.gauss_kronrod_rule(15);
        assert!(Arc::ptr_eq(&gk15, &pool.gauss_kronrod_rule(14)));
        assert_eq!(gk15.0.len(), 15);
        assert!(gk15.1.iter().all(|&w| w > 0.0));
        let gk15_x22: f64 = gk15
            .0
            .iter()
            .zip(&gk15.1)
            .map(|(&x, &w)| w * x.powi(22))
            .sum();
        assert!((gk15_x22 - 2.0 / 23.0).abs() < 1e-14);

        let gk21 = pool.gauss_kronrod_rule(21);
        assert_eq!(gk21.0.len(), 21);
        let gk21_x30: f64 = gk21
            .0
            .iter()
            .zip(&gk21.1)
            .map(|(&x, &w)| w * x.powi(30))
            .sum();
        assert!((gk21_x30 - 2.0 / 31.0).abs() < 1e-14);

        let gk31 = pool.gauss_kronrod_rule(31);
        assert_eq!(gk31.0.len(), 31);
        let gk31_x46: f64 = gk31
            .0
            .iter()
            .zip(&gk31.1)
            .map(|(&x, &w)| w * x.powi(46))
            .sum();
        assert!((gk31_x46 - 2.0 / 47.0).abs() < 1e-13);

        // Clenshaw-Curtis nodes & exactness up to degree n.
        let cc8 = pool.clenshaw_curtis_rule(7);
        assert!(Arc::ptr_eq(&cc8, &pool.clenshaw_curtis_rule(8)));
        assert_eq!(cc8.0.len(), 9);
        assert!((cc8.0.first().copied().unwrap_or(0.0) - 1.0).abs() < 1e-15);
        assert!((cc8.0.last().copied().unwrap_or(0.0) + 1.0).abs() < 1e-15);
        let cc8_x8: f64 = cc8.0.iter().zip(&cc8.1).map(|(&x, &w)| w * x.powi(8)).sum();
        assert!((cc8_x8 - 2.0 / 9.0).abs() < 1e-14);

        // Custom quadrature panic safety + retryability.
        let failed = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.get_or_init_quadrature(QuadratureRuleKind::Custom(99), 3, || {
                panic!("quadrature build failed")
            })
        }));
        assert!(failed.is_err());
        let recovered = pool.get_or_init_quadrature(QuadratureRuleKind::Custom(99), 3, || {
            Arc::new((vec![-0.5, 0.5], vec![1.0, 1.0]))
        });
        assert_eq!(recovered.0, vec![-0.5, 0.5]);

        // Chebyshev nodes.
        let cg5 = pool.chebyshev_nodes(5);
        assert!(Arc::ptr_eq(&cg5, &pool.chebyshev_nodes(5)));
        assert_eq!(cg5.get(2).copied(), Some(0.0));
        let cl4 = pool.chebyshev_lobatto_nodes(4);
        assert_eq!(cl4.len(), 5);
        assert_eq!(cl4.first().copied(), Some(1.0));
        assert_eq!(cl4.get(2).copied(), Some(0.0));
        assert_eq!(cl4.last().copied(), Some(-1.0));
    }

    #[test]
    fn numeric_tables_twiddles_and_sieve_are_accurate_and_cached() {
        let pool = ComputeReusePool::new();

        let fact = pool.factorial_table(171);
        assert!(Arc::ptr_eq(&fact, &pool.factorial_table(171)));
        assert_eq!(fact.first().copied(), Some(1.0));
        assert_eq!(fact.get(5).copied(), Some(120.0));
        assert!(fact.get(170).copied().is_some_and(f64::is_finite));
        assert_eq!(fact.get(171).copied(), Some(f64::INFINITY));

        let log_fact = pool.log_factorial_table(6);
        assert!((log_fact.get(5).copied().unwrap_or(0.0) - 120.0_f64.ln()).abs() < 1e-15);

        let binom = pool.binomial_row(10);
        assert!(Arc::ptr_eq(&binom, &pool.binomial_row(10)));
        assert_eq!(binom.get(0).copied(), Some(1.0));
        assert_eq!(binom.get(5).copied(), Some(252.0));
        assert_eq!(binom.get(10).copied(), Some(1.0));

        let stirling = pool.stirling_second_table(5);
        assert!(Arc::ptr_eq(&stirling, &pool.stirling_second_table(5)));
        // S(5, 3) at index 5*6/2 + 3 = 18 -> 25.0
        assert_eq!(stirling.get(18).copied(), Some(25.0));

        let harmonic = pool.harmonic_table(4);
        assert!(Arc::ptr_eq(&harmonic, &pool.harmonic_table(4)));
        let expected_h4 = 1.0 + 0.5 + 1.0 / 3.0 + 0.25;
        assert!((harmonic.get(4).copied().unwrap_or(0.0) - expected_h4).abs() < 1e-15);

        let primes = pool.prime_sieve_cache(100);
        assert!(Arc::ptr_eq(&primes, &pool.prime_sieve_cache(100)));
        assert_eq!(primes.len(), 25);
        assert_eq!(primes.first().copied(), Some(2));
        assert_eq!(primes.last().copied(), Some(97));

        let twiddles = pool.real_fft_twiddles(8);
        assert!(Arc::ptr_eq(&twiddles, &pool.real_fft_twiddles(8)));
        assert_eq!(twiddles.len(), 7);
        assert!(twiddles.iter().all(|z| (z.norm() - 1.0).abs() < 1e-14));

        let roots = pool.unit_roots(8);
        assert!(Arc::ptr_eq(&roots, &pool.unit_roots(8)));
        assert_eq!(roots.len(), 8);
        let root0 = roots.first().copied().unwrap_or_default();
        assert!((root0 - Complex64::new(1.0, 0.0)).norm() < 1e-15);
    }

    #[test]
    fn sparse_ordering_plan_and_named_plan_registry_work() {
        let pool = ComputeReusePool::new();

        // 5-vertex path graph 0-1-2-3-4 in CSR.
        let row_ptrs = [0_usize, 1, 3, 5, 7, 8];
        let col_indices = [1_usize, 0, 2, 1, 3, 2, 4, 3];
        let plan = pool.sparse_ordering_plan(&row_ptrs, &col_indices);
        assert!(Arc::ptr_eq(
            &plan,
            &pool.sparse_ordering_plan(&row_ptrs, &col_indices)
        ));
        assert_eq!(plan.dim, 5);
        assert_eq!(plan.perm.len(), 5);
        assert_eq!(plan.inv_perm.len(), 5);
        for (new_idx, &old_idx) in plan.perm.iter().enumerate() {
            assert_eq!(plan.inv_perm.get(old_idx).copied(), Some(new_idx));
        }

        let custom = pool.get_or_build_sparse_ordering(
            SparseOrderingKind::Custom(7),
            &row_ptrs,
            &col_indices,
            || vec![4, 3, 2, 1, 0],
        );
        assert_eq!(custom.perm, vec![4, 3, 2, 1, 0]);
        assert_eq!(custom.inv_perm, vec![4, 3, 2, 1, 0]);

        // Named registry isolates distinct (name, K, V) caches and participates in clear_all().
        let c_vec = pool.named_plan_cache::<u32, Vec<f64>>("interp", 4);
        let c_str = pool.named_plan_cache::<u32, String>("interp", 4);
        let v1 = c_vec.get_or_insert_with(1, || Arc::new(vec![1.0, 2.0]));
        let v2 = pool
            .registry()
            .get_or_insert_with::<u32, Vec<f64>>("interp", 4, 1, || {
                panic!("same named cache must hit")
            });
        assert!(Arc::ptr_eq(&v1, &v2));
        let s1 = c_str.get_or_insert_with(1, || Arc::new("ok".to_owned()));
        assert_eq!(s1.as_str(), "ok");
        assert_eq!(pool.registry().cache_count(), 2);
        assert_eq!(pool.registry().len(), 2);

        drop(v1);
        drop(v2);
        drop(s1);
        pool.trim_idle();
        assert_eq!(pool.registry().len(), 0);

        let _ = c_vec.get_or_insert_with(2, || Arc::new(vec![3.0]));
        pool.clear_all();
        assert_eq!(c_vec.len(), 0);
        assert_eq!(pool.stats().cached_plans, 0);
    }
}
