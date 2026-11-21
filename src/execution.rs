//! CPU parallelism for the kernel and embedding hosts.

/// CPU parallelism reported by the host when first queried, without
/// starting worker threads. Embedders can use this to size their own pools.
#[must_use]
pub fn detected_cpu_threads() -> usize {
    static DETECTED: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *DETECTED.get_or_init(|| std::thread::available_parallelism().map_or(1, usize::from))
}

/// Default worker count for a host creating its own CPU pool. Uses a positive
/// `RAYON_NUM_THREADS` override or the CPUs available to the current process.
/// This query does not start worker threads.
#[must_use]
pub fn default_cpu_threads() -> usize {
    std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|text| text.parse::<std::num::NonZeroUsize>().ok())
        .map_or_else(detected_cpu_threads, std::num::NonZeroUsize::get)
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

/// Whether an indexed workload has enough work and threads to split.
#[must_use]
pub fn should_split(
    items: usize,
    minimum: usize,
) -> bool {
    items >= minimum.max(2) && cpu_threads() > 1
}
