//! CPU parallelism, SIMD capability probing, and execution policy for the kernel
//! and embedding hosts.

use std::cell::Cell;
use std::num::NonZeroUsize;
use std::sync::OnceLock;

const DEFAULT_PARALLEL_MIN_WORK: usize = 16_384;
const DEFAULT_L1_TILE_F64: usize = 1_024;
const DEFAULT_L2_TILE_F64: usize = 16_384;
const TARGET_CHUNKS_PER_THREAD: usize = 4;

thread_local! {
    static SCOPED_THREAD_LIMIT: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Runtime-detected CPU SIMD instruction-set capabilities.
#[expect(
    clippy::struct_excessive_bools,
    reason = "Host diagnostics report each independent SIMD capability"
)]
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct SimdSupport {
    /// x86/`x86_64` AVX-512 Foundation (`avx512f`).
    pub avx512f: bool,
    /// x86/`x86_64` AVX2 (`avx2`).
    pub avx2: bool,
    /// x86/`x86_64` FMA (`fma`).
    pub fma: bool,
    /// x86/`x86_64` SSE4.2 (`sse4.2`).
    pub sse42: bool,
    /// `AArch64` NEON (`neon`).
    pub neon: bool,
}

impl SimdSupport {
    /// Probes runtime CPU SIMD capabilities without spawning worker threads.
    #[must_use]
    pub fn detect() -> Self {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            Self {
                avx512f: std::arch::is_x86_feature_detected!("avx512f"),
                avx2: std::arch::is_x86_feature_detected!("avx2"),
                fma: std::arch::is_x86_feature_detected!("fma"),
                sse42: std::arch::is_x86_feature_detected!("sse4.2"),
                neon: false,
            }
        }
        #[cfg(target_arch = "aarch64")]
        {
            Self {
                avx512f: false,
                avx2: false,
                fma: false,
                sse42: false,
                neon: std::arch::is_aarch64_feature_detected!("neon"),
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
        {
            Self::default()
        }
    }

    /// Preferred `f64` SIMD lane width (`8`, `4`, `2`, or `1`).
    #[must_use]
    pub const fn f64_lanes(self) -> usize {
        if self.avx512f {
            8
        } else if self.avx2 {
            4
        } else if self.sse42 || self.neon {
            2
        } else {
            1
        }
    }
}

/// Execution policy for parallelizable numeric and symbolic workloads.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum ExecutionMode {
    /// Automatically uses the active Rayon pool or detected host CPUs.
    #[default]
    Auto,
    /// Forces single-threaded sequential execution.
    Serial,
    /// Requests parallel execution capped at the specified worker count.
    Parallel(usize),
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum ParallelMinWork {
    Default,
    Threshold(usize),
    Multiplier(f64),
}

/// Cached host hardware and runtime execution capabilities.
///
/// Probing runs once via [`host_environment`] without spawning worker threads.
#[derive(Clone, Debug, PartialEq)]
pub struct HostEnvironment {
    cpu_threads: usize,
    simd: SimdSupport,
    f64_lanes: usize,
    l1_tile_f64: usize,
    l2_tile_f64: usize,
    execution_mode: ExecutionMode,
    min_work_override: ParallelMinWork,
    parallel_min_work: usize,
}

impl HostEnvironment {
    /// Logical CPU threads detected on the host.
    #[must_use]
    pub const fn cpu_threads(&self) -> usize {
        self.cpu_threads
    }

    /// Runtime-detected SIMD capabilities.
    #[must_use]
    pub const fn simd(&self) -> SimdSupport {
        self.simd
    }

    /// Preferred `f64` SIMD lane width (`8`, `4`, `2`, or `1`).
    #[must_use]
    pub const fn f64_lanes(&self) -> usize {
        self.f64_lanes
    }

    /// Cache-friendly L1 tile size in `f64` elements.
    #[must_use]
    pub const fn l1_tile_f64(&self) -> usize {
        self.l1_tile_f64
    }

    /// Cache-friendly L2 tile size in `f64` elements.
    #[must_use]
    pub const fn l2_tile_f64(&self) -> usize {
        self.l2_tile_f64
    }

    /// Resolved default execution policy from environment configuration.
    #[must_use]
    pub const fn execution_mode(&self) -> ExecutionMode {
        self.execution_mode
    }

    /// Resolved minimum work threshold from `RSSN_PARALLEL_MIN_WORK`.
    #[must_use]
    pub const fn parallel_min_work(&self) -> usize {
        self.parallel_min_work
    }

    #[allow(clippy::single_call_fn)]
    fn adjusted_min_items(
        &self,
        minimum: usize,
    ) -> usize {
        if minimum == usize::MAX {
            return usize::MAX;
        }
        let base = minimum.max(2);
        match self.min_work_override {
            | ParallelMinWork::Multiplier(factor) => scale_count(base, factor, 2),
            | ParallelMinWork::Default | ParallelMinWork::Threshold(_) => base,
        }
    }

    #[allow(clippy::single_call_fn)]
    fn from_env_inputs(
        rssn_mode: Option<&str>,
        rayon_num_threads: Option<&str>,
        rssn_min_work: Option<&str>,
        detected_threads: usize,
        max_threads: usize,
        simd: SimdSupport,
    ) -> Self {
        let cpu_threads = detected_threads.max(1);
        let f64_lanes = simd.f64_lanes().max(1);
        let l1_tile_f64 = DEFAULT_L1_TILE_F64
            .checked_div(f64_lanes)
            .map_or(DEFAULT_L1_TILE_F64, |blocks| {
                blocks.max(1).saturating_mul(f64_lanes)
            });
        let l2_tile_f64 = DEFAULT_L2_TILE_F64
            .checked_div(l1_tile_f64)
            .map_or(l1_tile_f64, |blocks| {
                blocks.max(1).saturating_mul(l1_tile_f64)
            });
        let execution_mode =
            parse_execution_mode(rssn_mode, rayon_num_threads, cpu_threads, max_threads);
        let min_work_override = parse_parallel_min_work(rssn_min_work);
        let parallel_min_work = match min_work_override {
            | ParallelMinWork::Default => DEFAULT_PARALLEL_MIN_WORK,
            | ParallelMinWork::Threshold(value) => value.max(1),
            | ParallelMinWork::Multiplier(factor) => {
                scale_count(DEFAULT_PARALLEL_MIN_WORK, factor, 1)
            },
        };
        Self {
            cpu_threads,
            simd,
            f64_lanes,
            l1_tile_f64,
            l2_tile_f64,
            execution_mode,
            min_work_override,
            parallel_min_work,
        }
    }

    #[allow(clippy::single_call_fn)]
    fn probe() -> Self {
        let rssn_mode = std::env::var("RSSN_EXECUTION_MODE").ok();
        let rayon_num_threads = std::env::var("RAYON_NUM_THREADS").ok();
        let rssn_min_work = std::env::var("RSSN_PARALLEL_MIN_WORK").ok();
        Self::from_env_inputs(
            rssn_mode.as_deref(),
            rayon_num_threads.as_deref(),
            rssn_min_work.as_deref(),
            detected_cpu_threads(),
            rayon::max_num_threads(),
            SimdSupport::detect(),
        )
    }
}

struct ScopedLimitGuard {
    previous: Option<usize>,
}

impl ScopedLimitGuard {
    fn install(limit: Option<usize>) -> Self {
        let previous = SCOPED_THREAD_LIMIT.with(|cell| cell.replace(limit));
        Self { previous }
    }
}

impl Drop for ScopedLimitGuard {
    fn drop(&mut self) {
        SCOPED_THREAD_LIMIT.with(|cell| cell.set(self.previous));
    }
}

fn parse_positive_usize(raw: Option<&str>) -> Option<usize> {
    raw.map(str::trim)
        .filter(|text| !text.is_empty())
        .and_then(|text| text.parse::<NonZeroUsize>().ok())
        .map(NonZeroUsize::get)
}

#[allow(clippy::single_call_fn)]
fn parse_execution_mode(
    rssn_mode: Option<&str>,
    rayon_num_threads: Option<&str>,
    detected_threads: usize,
    max_threads: usize,
) -> ExecutionMode {
    let cap = max_threads.max(1);
    let rayon_threads = parse_positive_usize(rayon_num_threads).map(|n| n.min(cap));
    let fallback_parallel = rayon_threads.unwrap_or_else(|| detected_threads.max(1).min(cap));

    if let Some(raw) = rssn_mode.map(str::trim).filter(|text| !text.is_empty()) {
        if raw.eq_ignore_ascii_case("auto") {
            return ExecutionMode::Auto;
        }
        if raw.eq_ignore_ascii_case("serial") {
            return ExecutionMode::Serial;
        }
        if raw.eq_ignore_ascii_case("parallel") {
            return ExecutionMode::Parallel(fallback_parallel);
        }
        if let Some((prefix, count)) = raw.split_once(':') {
            let parsed = prefix
                .trim()
                .eq_ignore_ascii_case("parallel")
                .then(|| parse_positive_usize(Some(count)))
                .flatten();
            return parsed.map_or(ExecutionMode::Auto, |n| ExecutionMode::Parallel(n.min(cap)));
        }
        return ExecutionMode::Auto;
    }

    match rayon_threads {
        | Some(1) => ExecutionMode::Serial,
        | Some(n) => ExecutionMode::Parallel(n),
        | None => ExecutionMode::Auto,
    }
}

fn parse_positive_f64(raw: &str) -> Option<f64> {
    raw.trim()
        .parse::<f64>()
        .ok()
        .filter(|factor| factor.is_finite() && *factor > 0.0)
}

#[allow(clippy::single_call_fn)]
fn parse_parallel_min_work(raw: Option<&str>) -> ParallelMinWork {
    let Some(text) = raw.map(str::trim).filter(|value| !value.is_empty()) else {
        return ParallelMinWork::Default;
    };
    if let Some(multiplier) = text
        .strip_suffix(['x', 'X'])
        .or_else(|| text.strip_prefix('*'))
        .and_then(parse_positive_f64)
    {
        return ParallelMinWork::Multiplier(multiplier);
    }
    if text.ends_with(['x', 'X']) || text.starts_with('*') {
        return ParallelMinWork::Default;
    }
    if let Some(threshold) = parse_positive_usize(Some(text)) {
        return ParallelMinWork::Threshold(threshold);
    }
    if text.contains('.') {
        return parse_positive_f64(text)
            .map_or(ParallelMinWork::Default, ParallelMinWork::Multiplier);
    }
    ParallelMinWork::Default
}

fn scale_count(
    base: usize,
    factor: f64,
    floor: usize,
) -> usize {
    let scaled = (base as f64) * factor;
    if !scaled.is_finite() || scaled <= 0.0 {
        return floor;
    }
    let rounded = scaled.round();
    if rounded >= usize::MAX as f64 {
        usize::MAX
    } else {
        (rounded as usize).max(floor)
    }
}

#[allow(clippy::single_call_fn)]
fn should_split_work_with(
    items: usize,
    work_per_item: usize,
    min_work: usize,
    threads: usize,
) -> bool {
    items >= 2 && threads > 1 && items.saturating_mul(work_per_item) >= min_work.max(1)
}

#[allow(clippy::single_call_fn)]
fn choose_chunk_size_for(
    items: usize,
    min_chunk: usize,
    threads: usize,
) -> usize {
    let min_chunk = min_chunk.max(1);
    if threads <= 1 {
        return items.max(min_chunk);
    }
    let target_chunks = threads.saturating_mul(TARGET_CHUNKS_PER_THREAD).max(1);
    items.div_ceil(target_chunks).max(min_chunk)
}

/// Returns the cached [`HostEnvironment`] probed once without starting worker threads.
#[must_use]
pub fn host_environment() -> &'static HostEnvironment {
    static HOST_ENV: OnceLock<HostEnvironment> = OnceLock::new();
    HOST_ENV.get_or_init(HostEnvironment::probe)
}

/// Returns the resolved default [`ExecutionMode`] from environment configuration.
#[must_use]
pub fn default_execution_mode() -> ExecutionMode {
    host_environment().execution_mode()
}

/// CPU parallelism reported by the host when first queried, without
/// starting worker threads. Embedders can use this to size their own pools.
#[must_use]
pub fn detected_cpu_threads() -> usize {
    static DETECTED: OnceLock<usize> = OnceLock::new();
    *DETECTED.get_or_init(|| std::thread::available_parallelism().map_or(1, usize::from))
}

/// Default worker count for a host creating its own CPU pool. Uses a positive
/// `RAYON_NUM_THREADS` override or the CPUs available to the current process.
/// This query does not start worker threads.
#[must_use]
pub fn default_cpu_threads() -> usize {
    let env_threads = std::env::var("RAYON_NUM_THREADS").ok();
    parse_positive_usize(env_threads.as_deref())
        .unwrap_or_else(detected_cpu_threads)
        .min(rayon::max_num_threads())
}

/// Number of threads in the active Rayon pool.
///
/// Outside a Rayon worker this is the global pool. It detects available CPUs and respects
/// `RAYON_NUM_THREADS`. Nested kernel work uses the enclosing pool.
#[must_use]
pub fn cpu_threads() -> usize {
    rayon::current_num_threads()
}

/// Effective worker thread count respecting [`ExecutionMode`], scoped thread limits,
/// and the active Rayon pool.
#[must_use]
pub fn effective_threads() -> usize {
    let scoped = SCOPED_THREAD_LIMIT.with(Cell::get);
    if scoped.is_some_and(|limit| limit <= 1) {
        return 1;
    }
    let mode_threads = match default_execution_mode() {
        | ExecutionMode::Serial => return 1,
        | ExecutionMode::Parallel(limit) if limit <= 1 => return 1,
        | ExecutionMode::Auto => cpu_threads().max(1),
        | ExecutionMode::Parallel(limit) => cpu_threads().min(limit).max(1),
    };
    scoped.map_or(mode_threads, |limit| mode_threads.min(limit).max(1))
}

/// Whether an indexed workload has enough work and threads to split.
///
/// Honors [`ExecutionMode::Serial`] (returning `false` without spawning Rayon
/// threads) and defaults to [`ExecutionMode::Auto`] using the active Rayon pool.
#[must_use]
pub fn should_split(
    items: usize,
    minimum: usize,
) -> bool {
    let env = host_environment();
    match env.execution_mode {
        | ExecutionMode::Serial => return false,
        | ExecutionMode::Parallel(limit) if limit <= 1 => return false,
        | ExecutionMode::Auto | ExecutionMode::Parallel(_) => {},
    }
    if SCOPED_THREAD_LIMIT
        .with(Cell::get)
        .is_some_and(|limit| limit <= 1)
    {
        return false;
    }
    let min_items = env.adjusted_min_items(minimum);
    items >= min_items && effective_threads() > 1
}

/// Whether an indexed workload with `work_per_item` estimated cost per element
/// exceeds `RSSN_PARALLEL_MIN_WORK` and has parallel workers available.
#[must_use]
pub fn should_split_work(
    items: usize,
    work_per_item: usize,
) -> bool {
    let env = host_environment();
    match env.execution_mode {
        | ExecutionMode::Serial => return false,
        | ExecutionMode::Parallel(limit) if limit <= 1 => return false,
        | ExecutionMode::Auto | ExecutionMode::Parallel(_) => {},
    }
    if SCOPED_THREAD_LIMIT
        .with(Cell::get)
        .is_some_and(|limit| limit <= 1)
    {
        return false;
    }
    if items < 2 || items.saturating_mul(work_per_item) < env.parallel_min_work() {
        return false;
    }
    should_split_work_with(
        items,
        work_per_item,
        env.parallel_min_work(),
        effective_threads(),
    )
}

/// Selects a positive chunk size for `items` elements that respects `min_chunk`
/// and balances work across [`effective_threads`].
#[must_use]
pub fn choose_chunk_size(
    items: usize,
    min_chunk: usize,
) -> usize {
    choose_chunk_size_for(items, min_chunk, effective_threads())
}

/// Returns a thread budget up to `requested` that never oversubscribes the
/// active execution scope or enclosing Rayon pool.
///
/// Inside an active Rayon worker thread this returns `1` so nested callers do
/// not spawn additional worker pools; outside a pool it caps `requested` by
/// [`effective_threads`].
#[must_use]
pub fn nested_thread_limit(requested: usize) -> usize {
    if requested <= 1 || rayon::current_thread_index().is_some() {
        1
    } else {
        effective_threads().min(requested).max(1)
    }
}

/// Runs `f` within a thread budget capped at `max_threads` without oversubscribing.
///
/// * In [`ExecutionMode::Serial`] or when `max_threads` is `Some(0)` or `Some(1)`,
///   `f` runs directly on the calling thread with nested splitting disabled.
/// * When already inside a Rayon worker thread, `f` runs on the enclosing pool
///   under a tightened thread-local budget rather than spawning a nested pool.
/// * Otherwise, when a smaller positive thread cap is requested, a scoped Rayon
///   pool of that size runs `f` (falling back to serial execution on the calling
///   thread if thread-pool creation fails).
pub fn with_scoped_pool<R, F>(
    max_threads: Option<usize>,
    f: F,
) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    let env = host_environment();
    let base_cap = match env.execution_mode {
        | ExecutionMode::Serial => 1,
        | ExecutionMode::Parallel(limit) => limit.max(1).min(rayon::max_num_threads()),
        | ExecutionMode::Auto => default_cpu_threads().max(1),
    };
    let current_scoped = SCOPED_THREAD_LIMIT.with(Cell::get);
    let requested_cap = max_threads.map_or(base_cap, |limit| limit.max(1).min(base_cap));
    let target = current_scoped.map_or(requested_cap, |prev| prev.max(1).min(requested_cap));

    if target <= 1 {
        let _guard = ScopedLimitGuard::install(Some(1));
        return f();
    }

    if rayon::current_thread_index().is_some() {
        let enclosing = rayon::current_num_threads().max(1);
        let bounded = target.min(enclosing);
        let _guard = ScopedLimitGuard::install(Some(bounded));
        return f();
    }

    let pool_cap = default_cpu_threads().max(1);
    let needs_custom_pool = target < pool_cap
        && match env.execution_mode {
            | ExecutionMode::Parallel(_) => true,
            | ExecutionMode::Auto => max_threads.is_some(),
            | ExecutionMode::Serial => false,
        };

    if needs_custom_pool {
        return match rayon::ThreadPoolBuilder::new().num_threads(target).build() {
            | Ok(pool) => {
                pool.install(|| {
                    let _guard = ScopedLimitGuard::install(Some(target));
                    f()
                })
            },
            | Err(_) => {
                let _guard = ScopedLimitGuard::install(Some(1));
                f()
            },
        };
    }

    let _guard = ScopedLimitGuard::install(Some(target));
    f()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn host_environment_probes_consistently_and_caches_singleton() {
        let first = host_environment();
        let second = host_environment();
        assert!(std::ptr::eq(first, second));
        assert_eq!(first.cpu_threads(), detected_cpu_threads());
        assert!(first.cpu_threads() >= 1);
        assert!(matches!(first.f64_lanes(), 1 | 2 | 4 | 8));
        assert_eq!(first.f64_lanes(), first.simd().f64_lanes());
        assert!(first.l1_tile_f64() > 0);
        assert!(first.l2_tile_f64() >= first.l1_tile_f64());
        assert_eq!(first.l1_tile_f64().checked_rem(first.f64_lanes()), Some(0));
        assert_eq!(
            first.l2_tile_f64().checked_rem(first.l1_tile_f64()),
            Some(0)
        );
        assert!(first.parallel_min_work() >= 1);
        assert_eq!(default_execution_mode(), first.execution_mode());
    }

    #[test]
    fn simd_support_f64_lanes_covers_all_widths() {
        assert_eq!(SimdSupport::default().f64_lanes(), 1);
        assert_eq!(
            SimdSupport {
                neon: true,
                ..SimdSupport::default()
            }
            .f64_lanes(),
            2
        );
        assert_eq!(
            SimdSupport {
                sse42: true,
                ..SimdSupport::default()
            }
            .f64_lanes(),
            2
        );
        assert_eq!(
            SimdSupport {
                avx2: true,
                sse42: true,
                ..SimdSupport::default()
            }
            .f64_lanes(),
            4
        );
        assert_eq!(
            SimdSupport {
                avx512f: true,
                avx2: true,
                fma: true,
                sse42: true,
                neon: false,
            }
            .f64_lanes(),
            8
        );

        let probed = SimdSupport::detect();
        if probed.avx512f {
            assert_eq!(probed.f64_lanes(), 8);
        } else if probed.avx2 {
            assert_eq!(probed.f64_lanes(), 4);
        } else if probed.sse42 || probed.neon {
            assert_eq!(probed.f64_lanes(), 2);
        } else {
            assert_eq!(probed.f64_lanes(), 1);
        }
    }

    #[test]
    fn execution_mode_parsing_handles_all_variants_and_fallbacks() {
        assert_eq!(parse_execution_mode(None, None, 8, 64), ExecutionMode::Auto);
        assert_eq!(
            parse_execution_mode(Some(""), Some("4"), 8, 64),
            ExecutionMode::Parallel(4)
        );
        assert_eq!(
            parse_execution_mode(Some("auto"), Some("4"), 8, 64),
            ExecutionMode::Auto
        );
        assert_eq!(
            parse_execution_mode(Some("  SERIAL  "), Some("8"), 8, 64),
            ExecutionMode::Serial
        );
        assert_eq!(
            parse_execution_mode(Some("parallel"), None, 8, 64),
            ExecutionMode::Parallel(8)
        );
        assert_eq!(
            parse_execution_mode(Some("parallel"), Some("4"), 8, 64),
            ExecutionMode::Parallel(4)
        );
        assert_eq!(
            parse_execution_mode(Some("parallel:3"), Some("8"), 8, 64),
            ExecutionMode::Parallel(3)
        );
        assert_eq!(
            parse_execution_mode(Some("PARALLEL : 12"), None, 8, 64),
            ExecutionMode::Parallel(12)
        );
        assert_eq!(
            parse_execution_mode(Some("parallel:999"), None, 8, 16),
            ExecutionMode::Parallel(16)
        );
        assert_eq!(
            parse_execution_mode(Some("parallel:0"), None, 8, 64),
            ExecutionMode::Auto
        );
        assert_eq!(
            parse_execution_mode(Some("parallel:abc"), None, 8, 64),
            ExecutionMode::Auto
        );
        assert_eq!(
            parse_execution_mode(Some("bogus"), None, 8, 64),
            ExecutionMode::Auto
        );
        assert_eq!(
            parse_execution_mode(None, Some("1"), 8, 64),
            ExecutionMode::Serial
        );
        assert_eq!(
            parse_execution_mode(None, Some("6"), 8, 64),
            ExecutionMode::Parallel(6)
        );
        assert_eq!(
            parse_execution_mode(None, Some("0"), 8, 64),
            ExecutionMode::Auto
        );
    }

    #[test]
    fn parallel_min_work_parsing_handles_thresholds_and_multipliers() {
        assert_eq!(parse_parallel_min_work(None), ParallelMinWork::Default);
        assert_eq!(
            parse_parallel_min_work(Some("4096")),
            ParallelMinWork::Threshold(4_096)
        );
        assert_eq!(
            parse_parallel_min_work(Some(" 2x ")),
            ParallelMinWork::Multiplier(2.0)
        );
        assert_eq!(
            parse_parallel_min_work(Some("0.5X")),
            ParallelMinWork::Multiplier(0.5)
        );
        assert_eq!(
            parse_parallel_min_work(Some("*4")),
            ParallelMinWork::Multiplier(4.0)
        );
        assert_eq!(
            parse_parallel_min_work(Some("1.5")),
            ParallelMinWork::Multiplier(1.5)
        );
        for invalid in ["", "0", "-4", "0x", "NaNx", "infx", "bogus"] {
            assert_eq!(
                parse_parallel_min_work(Some(invalid)),
                ParallelMinWork::Default
            );
        }

        let env_mult = HostEnvironment::from_env_inputs(
            Some("auto"),
            None,
            Some("2x"),
            8,
            64,
            SimdSupport::default(),
        );
        assert_eq!(
            env_mult.parallel_min_work(),
            DEFAULT_PARALLEL_MIN_WORK.saturating_mul(2)
        );
        assert_eq!(env_mult.adjusted_min_items(10), 20);
        assert_eq!(env_mult.adjusted_min_items(0), 4);
        assert_eq!(env_mult.adjusted_min_items(usize::MAX), usize::MAX);

        let env_thresh = HostEnvironment::from_env_inputs(
            Some("auto"),
            None,
            Some("512"),
            8,
            64,
            SimdSupport::default(),
        );
        assert_eq!(env_thresh.parallel_min_work(), 512);
        assert_eq!(env_thresh.adjusted_min_items(2), 2);
        assert_eq!(env_thresh.adjusted_min_items(0), 2);

        let env_tiny = HostEnvironment::from_env_inputs(
            Some("auto"),
            None,
            Some("0.25x"),
            8,
            64,
            SimdSupport::default(),
        );
        assert_eq!(env_tiny.adjusted_min_items(usize::MAX), usize::MAX);
        assert_eq!(env_tiny.adjusted_min_items(100), 25);
        assert_eq!(env_tiny.adjusted_min_items(0), 2);
    }

    #[test]
    fn should_split_work_and_chunk_size_invariants() {
        assert!(!should_split_work_with(0, 100_000, 16_384, 8));
        assert!(!should_split_work_with(1, 100_000, 16_384, 8));
        assert!(!should_split_work_with(1_000, 0, 16_384, 8));
        assert!(!should_split_work_with(1_000, 16, 16_384, 8));
        assert!(should_split_work_with(1_024, 16, 16_384, 8));
        assert!(!should_split_work_with(1_024, 16, 16_384, 1));
        assert!(should_split_work_with(usize::MAX, usize::MAX, 16_384, 4));

        assert!(!should_split_work(0, 100_000));
        assert!(!should_split_work(1, 100_000));

        for items in [0usize, 1, 5, 64, 1_000, 100_000] {
            for min_chunk in [0usize, 1, 16, 64] {
                for threads in [1usize, 2, 8] {
                    let chunk = choose_chunk_size_for(items, min_chunk, threads);
                    assert!(chunk >= 1);
                    assert!(chunk >= min_chunk.max(1));
                }
                let live_chunk = choose_chunk_size(items, min_chunk);
                assert!(live_chunk >= min_chunk.max(1));
            }
        }
        assert_eq!(choose_chunk_size_for(1_000, 16, 1), 1_000);
        assert_eq!(choose_chunk_size_for(1_000, 16, 4), 63);
    }

    #[test]
    fn scoped_pool_limits_threads_and_prevents_nested_oversubscription() {
        with_scoped_pool(Some(1), || {
            assert_eq!(effective_threads(), 1);
            assert!(!should_split(10_000, 2));
            assert!(!should_split_work(10_000, 1_000));
            assert_eq!(nested_thread_limit(8), 1);
        });

        with_scoped_pool(Some(2), || {
            assert!(effective_threads() <= 2);
            assert!(effective_threads() >= 1);
            with_scoped_pool(Some(8), || {
                assert!(effective_threads() <= 2);
            });
            with_scoped_pool(Some(0), || {
                assert_eq!(effective_threads(), 1);
                assert!(!should_split(10_000, 2));
            });
            assert!(effective_threads() <= 2);
        });

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap_or_else(|error| panic!("cannot start test pool: {error}"));
        pool.install(|| {
            assert_eq!(nested_thread_limit(8), 1);
            with_scoped_pool(Some(4), || {
                assert!(effective_threads() <= 2);
                assert_eq!(rayon::current_num_threads(), 2);
            });
            with_scoped_pool(Some(1), || {
                assert_eq!(effective_threads(), 1);
                assert!(!should_split(10_000, 2));
            });
        });

        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap_or_else(|error| panic!("cannot start single-thread pool: {error}"));
        one.install(|| {
            assert_eq!(effective_threads(), 1);
            assert!(!should_split(10_000, 2));
            assert!(!should_split_work(10_000, 1_000));
        });
    }

    #[test]
    fn backwards_compatible_thread_queries_are_consistent() {
        assert!(detected_cpu_threads() >= 1);
        assert!(default_cpu_threads() >= 1);
        assert!(cpu_threads() >= 1);
        assert!(!should_split(0, 0));
        assert!(!should_split(1, 0));
        assert!(!should_split(100, usize::MAX));
    }
}
