//! The quadratic sieve and its self-initialising variant (MPQS).
//!
//! `qs_factor(n, limits)` finds a non-trivial divisor of an odd composite
//! `n` that has no small prime factor: a factor base of primes with
//! `(n/p) = 1`, sieving of `Q(x) = (ax + b)^2 - n` over a symmetric
//! interval, a single-large-prime variation, and a GF(2) dependency solve
//! by sparse elimination. The rule set's `factorint` uses it for
//! composites that survive trial division, perfect-power detection, BPSW
//! and ECM.

/// Limits for one quadratic-sieve attempt.
#[derive(Copy, Clone, Debug)]
pub struct QsLimits {
    /// Factor-base bound.
    pub bound: u64,
    /// Half-width of the sieving interval.
    pub interval: u64,
    /// Maximum number of polynomials tried before giving up.
    pub polynomials: usize,
}
