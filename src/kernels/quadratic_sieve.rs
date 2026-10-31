//! The quadratic sieve and its self-initialising variant (MPQS).
//!
//! `qs_factor(n, limits)` finds a non-trivial divisor of an odd composite
//! `n` that has no small prime factor: a factor base of primes with
//! `(n/p) = 1`, sieving of `Q(x) = (ax + b)^2 - n` over a symmetric
//! interval, a large-prime variation, and a structured GF(2) dependency
//! solve. The rule set's `factorint` uses it for composites that survive
//! trial division, perfect-power detection and BPSW; below 40 digits the
//! elliptic curves get the first attempt.
//!
//! # Algorithm
//!
//! This is a multiple-polynomial quadratic sieve (MPQS) of the simple
//! `a = q^2` form. Every polynomial picks a prime `q` near
//! `sqrt(sqrt(2 n) / M)`, computes `b` with `b^2 = n (mod q^2)` by
//! Tonelli-Shanks followed by one Hensel lift, and sets
//! `a = q^2`, `c = (b^2 - n) / a`. Then
//!
//! ```text
//! Q(x) = (a x + b)^2 - n = a * g(x),
//! g(x) = a x^2 + 2 b x + c
//! ```
//!
//! is sieved in `x in [-M, M]` with a logarithmic approximation (the
//! factor `a = q^2` is a square and only enters a relation through `q`);
//! the sign of `g` is the `-1` exponent of the relation. Both parities of
//! `x` are sieved: `g(x)` is even for one parity and odd for the other.
//! Candidates whose accumulated `log p` is close to `log |g(x)|` are
//! factored exactly over the factor base; a cofactor below
//! `large_prime_factor * bound` becomes a partial relation and is combined
//! with the next relation that has the same cofactor. Every prime factor
//! of the cofactor is larger than the factor-base bound, so a cofactor
//! below `bound^2` is prime; with the default multiplier (< bound) a
//! partial relation therefore carries exactly one large prime.
//!
//! Full relations are reduced over GF(2) by a structured elimination.
//! Columns are taken lightest-first: a singleton deletes its row, and a
//! heavier column xors its shortest row into the other rows containing it
//! (tracking which original relations each row is the xor of) until fill-in
//! stops paying off; the small remaining core is eliminated densely with
//! the same tracking. The elimination stops after 64 dependencies, which is
//! far more than needed because each one succeeds with probability about
//! one half. For every dependency the two square roots
//! `X = prod (a x + b)` and `Y = prod q * prod p^(e/2)` are formed and
//! `gcd(X - Y, n)`, `gcd(X + Y, n)` are tried.
//!
//! The sieve is not trial division and not Pollard's rho: the work is
//! dominated by the logarithmic sieve over the interval, as a real QS.
//!
//! # Practical range and work bound
//!
//! [`QsLimits::for_digits`] follows the standard quadratic-sieve
//! complexity heuristic. The factor-base bound is `L_n[1/2, 1/2]` (the
//! optimum for the relation-collection cost), rounded to two significant
//! digits; the sieve interval is of the same order; the large-prime bound
//! is a digit-dependent multiple of the factor base; and the relation
//! target is the factor-base width plus a digit-dependent safety margin.
//! The polynomial budget and the explicit `max_sieve_positions` budget
//! come from `L_n[1/2, 1]`, capped so that a hopeless input exhausts a
//! bounded amount of work instead of running forever.
//!
//! In the test profile (`opt-level = 1`) the design point is roughly
//! 40-55 digits; 55 digits typically factor in seconds, 60 digits in tens
//! of seconds on a current multi-core machine, and 65 digits is
//! best-effort (factor within budget or `None`). The parameters are
//! calibrated through 80 digits so that callers can ask for a larger
//! budget, but the factor base is capped at
//! [`MAX_FACTOR_BASE_COLUMNS`] columns: an 80-digit input builds its
//! parameter set, then returns `None` before allocating the sieve because
//! the relation matrix would not fit the documented memory bound.
//! [`QsLimits::MAX_DIGITS`] is the largest calibrated size.
//!
//! Every step is deterministic: the polynomial primes are taken in
//! increasing order, relations are visited in a fixed order inside each
//! fixed-size batch (the batches may run in parallel, but they are
//! collected in index order), and the elimination is deterministic. The
//! divisor found therefore does not depend on the thread count, and no
//! random source is used.
//!
//! [`MAX_FACTOR_BASE_COLUMNS`]: QsLimits::MAX_FACTOR_BASE_COLUMNS

// The sieve hot loops index fixed-size arrays in step with their formulas
// and do index arithmetic that is bounded by the array lengths; the code
// reads closest to the mathematical description this way.
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)]
// The helpers below are small, named pieces of the algorithm and are used
// from exactly one place; splitting them further would obscure the flow.
#![allow(clippy::single_call_fn)]
// Prime and exponent counts become float limits for the logarithmic sieve.
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::kernels::number_theory::primes_sieve;

/// The logarithmic sieve stores `log2(p)` scaled by this factor.
const LOG_SCALE: f64 = 256.0;
/// Above this many decimal digits a request returns `None` immediately.
///
/// The parameters are calibrated through this size; the practical limit is
/// lower because the factor base may not exceed
/// [`QsLimits::MAX_FACTOR_BASE_COLUMNS`] columns.
const MAX_DIGITS: u32 = 80;
/// Largest number of factor-base columns (including the sign column) that
/// the relation matrix is allowed to reach.
///
/// The structured solve stores one bitset over the relations per row while
/// reducing the matrix, so this bounds both memory (roughly
/// `columns^2 / 8` bytes at the peak) and the elimination time. A larger
/// factor base is not attempted: the sieve returns `None` instead of
/// allocating an unbounded matrix.
pub const MAX_FACTOR_BASE_COLUMNS: usize = 40_000;
/// Trial division by the primes below this before the sieve starts.
const SMALL_TRIAL: u64 = 1_000;
/// Relation-count target over the factor-base size (dependencies needed).
const EXTRA_DEPENDENCIES: usize = 8;
/// Polynomials sieved per deterministic parallel batch.
///
/// The batch size is fixed, so the set of polynomials processed before a
/// solve, and therefore the result, does not depend on the thread count.
const POLYNOMIAL_BATCH: usize = 32;
/// Partials kept per factor-base column (a memory cap, not a math cap).
const PARTIALS_PER_COLUMN: usize = 2;
/// Heaviest column pivoted in the sparse phase of the GF(2) solve.
///
/// Singleton and doubleton columns are always worth eliminating; continuing
/// up to this weight shrinks the dense core by roughly two thirds on QS
/// matrices without letting fill-in take over.
const SPARSE_MAX_WEIGHT: usize = 64;
/// Longest pivot row used in the sparse phase (guards against fill-in).
const SPARSE_MAX_ROW: usize = 1_024;
/// Stop the sparse phase when at most this many rows remain.
///
/// The dense phase is cubic in the core width, so trading cheap sparse row
/// operations for a smaller core is a large win. The remaining nullity means
/// the dense column count is a little below this target.
const DENSE_CORE_TARGET: usize = 4_000;
/// Number of dependencies that is plenty to find a non-trivial factor.
///
/// Each dependency fails with probability about one half (odd `-1` or a
/// trivial gcd), so 64 are already overwhelming. Stopping the elimination
/// early avoids most of the dense work.
const MAX_DEPENDENCIES: usize = 64;

/// Limits for one quadratic-sieve attempt.
///
/// The defaults come from [`QsLimits::for_digits`] for 30-digit input.
/// The fields scale with the input size: a larger factor base needs more
/// relations (a larger GF(2) matrix), a larger interval yields more
/// candidates per polynomial, the large-prime factor widens the partial
/// relations, and [`QsLimits::polynomials`] /
/// [`QsLimits::max_sieve_positions`] are the hard work bounds.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct QsLimits {
    /// Factor-base bound: primes `p <= bound` with `(n/p) = 1` are used.
    pub bound: u64,
    /// Half-width `M` of the sieving interval `[-M, M]`.
    pub interval: u64,
    /// Maximum number of polynomials tried before giving up.
    pub polynomials: usize,
    /// Multiplier of the factor-base bound that caps the cofactor of a
    /// partial relation. Kept below `bound` so a partial has exactly one
    /// large prime.
    pub large_prime_factor: u64,
    /// Extra relations collected beyond the factor-base width, on top of
    /// the eight needed to have dependencies to try.
    pub extra_relations: usize,
    /// Hard bound on the number of sieved positions (polynomials times
    /// interval length) for the whole attempt.
    pub max_sieve_positions: u64,
}

impl QsLimits {
    /// Largest calibrated input size in decimal digits.
    ///
    /// Requests for larger numbers are refused before any sieving work;
    /// the caller can still use ECM there. The practical limit is lower
    /// (see [`MAX_FACTOR_BASE_COLUMNS`]).
    pub const MAX_DIGITS: u32 = MAX_DIGITS;
    /// Largest number of factor-base columns this implementation will try.
    pub const MAX_FACTOR_BASE_COLUMNS: usize = MAX_FACTOR_BASE_COLUMNS;

    /// Parameter set for a `digits`-digit number.
    ///
    /// The factor-base bound is the standard QS optimum
    /// `L_n[1/2, 1/2] = exp(0.5 sqrt(ln n ln ln n))`, rounded up to two
    /// significant digits; the sieve interval is of the same order. The
    /// large-prime multiplier grows with the size, the relation margin is
    /// a small fraction of the estimated factor-base width, and the work
    /// bound is `L_n[1/2, 1] = exp(sqrt(ln n ln ln n))` capped at
    /// `2 * 10^10` positions so that even a hopeless input has a bounded
    /// cost.
    ///
    /// For 80 digits this yields a factor base around five million and a
    /// sieve interval of two million, which is beyond the column limit:
    /// the attempt returns `None` before allocating the sieve. At 65 digits
    /// the matrix is already large (tens of thousands of columns and
    /// hundreds of megabytes at the peak), but the work bound still
    /// terminates; beyond about 66 digits the factor base exceeds the
    /// column limit, so the attempt returns `None` after construction.
    #[must_use]
    pub fn for_digits(digits: u32) -> Self {
        let digits = digits.clamp(4, MAX_DIGITS);
        let ln_n = f64::from(digits) * std::f64::consts::LN_10;
        let ln_ln_n = ln_n.ln();
        let optimum = (0.5 * (ln_n * ln_ln_n).sqrt()).exp();
        let bound = round_significant(optimum, 2).clamp(120.0, 6_000_000.0) as u64;
        let interval = round_significant(optimum, 2).clamp(3_000.0, 2_000_000.0) as u64;
        let large_prime_factor = match digits {
            | 0..=32 => 16,
            | 33..=52 => 24,
            | 53..=64 => 32,
            | _ => 48,
        };
        // Estimated factor-base width: about half the primes up to `bound`
        // have `(n/p) = 1`; use `bound / (2 ln bound)`.
        let estimated_columns = bound as f64 / (2.0 * (bound as f64).ln());
        let extra_relations = ((estimated_columns / 128.0).ceil() as usize).max(EXTRA_DEPENDENCIES);
        let heuristic = (ln_n * ln_ln_n).sqrt().exp();
        let max_sieve_positions = heuristic.min(2.0e10) as u64;
        let polynomials = usize::try_from((max_sieve_positions / interval).max(64))
            .unwrap_or(64)
            .min(20_000);
        Self {
            bound,
            interval,
            polynomials,
            large_prime_factor,
            extra_relations,
            max_sieve_positions,
        }
    }
}

/// Rounds `value` to `significant` decimal significant digits (upwards,
/// to a two-digit value such as `4_700_000`).
fn round_significant(
    value: f64,
    significant: i32,
) -> f64 {
    if !value.is_finite() || value <= 0.0 {
        return 0.0;
    }
    let exponent = value.log10().floor() as i32 - significant + 1;
    let scale = 10_f64.powi(exponent);
    (value / scale).ceil() * scale
}

impl Default for QsLimits {
    fn default() -> Self {
        Self::for_digits(30)
    }
}

/// Configuration of one sieve attempt.
///
/// `large_primes` turns the large-prime variation on or off; the
/// comparison is useful when measuring how much the variation helps. The
/// large-prime bound itself scales with the input and lives in
/// [`QsLimits::large_prime_factor`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct QsConfig {
    /// Factor-base, interval and polynomial budget.
    pub limits: QsLimits,
    /// Whether partial relations with one large prime are collected.
    pub large_primes: bool,
}

impl Default for QsConfig {
    fn default() -> Self {
        Self {
            limits: QsLimits::default(),
            large_primes: true,
        }
    }
}

/// Counters describing what one sieve attempt did.
///
/// They are informational only; the factorisation itself depends on the
/// input and the configuration, never on how the counters are read.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct QsStats {
    /// Polynomials actually sieved.
    pub polynomials: usize,
    /// Sieve positions that passed the logarithmic threshold.
    pub candidates: u64,
    /// Full (smooth) relations collected.
    pub full_relations: usize,
    /// Partial relations with one bounded large cofactor.
    pub partial_relations: usize,
    /// Partial relations that found a partner and became full relations.
    pub combined_relations: usize,
    /// Sieve positions actually processed (the work-budget unit).
    pub sieve_positions: u64,
    /// Rows in the matrix handed to the last GF(2) solve.
    pub matrix_rows: usize,
    /// Columns in the matrix handed to the last GF(2) solve.
    pub matrix_columns: usize,
    /// Dependencies found by the GF(2) elimination in the last solve.
    pub dependencies: usize,
}

/// Finds a non-trivial divisor of the odd composite `n`.
///
/// Returns `Some(d)` with `1 < d < n` and `d | n`, or `None` for invalid
/// input, for a prime, when `n` has more than [`QsLimits::MAX_DIGITS`]
/// decimal digits, when the relation matrix would exceed
/// [`QsLimits::MAX_FACTOR_BASE_COLUMNS`], and when the polynomial or
/// [`QsLimits::max_sieve_positions`] budget is exhausted. An even `n > 2`
/// returns `2` directly.
///
/// The input may contain small factors: trial division by the primes
/// below 1000 runs first and returns a factor directly. The digit limit is
/// checked before that trial division, so an oversized request never does
/// any sieving work. Perfect powers, prime input and even input are
/// otherwise handled before the sieve. The limits are not clamped upwards
/// or downwards to hide a bad choice: they are used as given, apart from a
/// safety clamp on `bound` and `interval`.
#[must_use]
pub fn qs_factor(
    n: &BigInt,
    limits: QsLimits,
) -> Option<BigInt> {
    qs_factor_with(
        n,
        QsConfig {
            limits,
            ..QsConfig::default()
        },
        &mut QsStats::default(),
    )
}

/// [`qs_factor`] with an explicit configuration and run statistics.
///
/// The same guarantees hold: the returned divisor is always in
/// `(1, n)` and divides `n`; `None` means "no factor found", not "n is
/// composite" or "n is prime".
#[must_use]
pub fn qs_factor_with(
    n: &BigInt,
    config: QsConfig,
    stats: &mut QsStats,
) -> Option<BigInt> {
    *stats = QsStats::default();
    let divisor = qs_core(n, config, stats)?;
    (divisor > BigInt::one() && divisor < *n).then_some(divisor)
}

/// Immutable parameters shared by every polynomial of one attempt.
#[derive(Copy, Clone, Debug)]
struct SieveParams {
    /// Half-width of the interval.
    interval: usize,
    /// Number of sieve positions, `2 * interval + 1`.
    length: usize,
    /// Largest cofactor accepted as a partial relation.
    large_bound: u64,
    /// Whether partial relations are collected.
    large_primes: bool,
    /// Log-scale slack for the candidate test.
    tolerance: f64,
}

/// Everything one sieved polynomial contributes.
#[derive(Default)]
struct SieveOutput {
    /// The polynomial was built and sieved.
    polynomial: bool,
    /// Per-polynomial counters.
    stats: QsStats,
    /// Full (smooth) relations, in increasing `x`.
    full: Vec<Relation>,
    /// `(large prime, relation)` pairs, in increasing `x`.
    partial: Vec<(u64, Relation)>,
    /// A factor found without a dependency (for example `q | n`).
    factor: Option<BigInt>,
}

/// The sieve proper; its result is validated by [`qs_factor_with`].
fn qs_core(
    n: &BigInt,
    config: QsConfig,
    stats: &mut QsStats,
) -> Option<BigInt> {
    if !n.is_positive() || *n < BigInt::from(2) {
        return None;
    }
    if decimal_digits(n) > MAX_DIGITS {
        return None;
    }
    if !n.bit(0) {
        return (*n > BigInt::from(2)).then_some(BigInt::from(2));
    }
    if let Some(p) = small_prime_divisor(n) {
        return Some(BigInt::from(p));
    }
    if is_probable_prime(n) {
        return None;
    }
    if let Some(root) = perfect_power_root(n) {
        return Some(root);
    }
    let limits = config.limits;
    let bound = limits.bound.clamp(30, 6_000_000);
    let interval = usize::try_from(limits.interval.clamp(100, 2_000_000)).ok()?;
    if limits.polynomials == 0 || limits.max_sieve_positions == 0 {
        return None;
    }
    let factor_base = match build_factor_base(n, bound) {
        | Ok(base) => base,
        | Err(divisor) => return Some(divisor),
    };
    if factor_base.len() < 3 {
        return None;
    }
    let columns = factor_base.len() + 1;
    if columns > MAX_FACTOR_BASE_COLUMNS {
        return None;
    }
    let divisors: Vec<BigInt> = factor_base.iter().map(|p| BigInt::from(p.prime)).collect();
    let large_prime_factor = limits
        .large_prime_factor
        .clamp(1, bound.saturating_sub(1).max(1));
    let large_bound = if config.large_primes {
        bound.saturating_mul(large_prime_factor)
    } else {
        1
    };
    // The sieve counts each factor-base prime once; the slack covers the
    // unsieved large cofactor and the rounding of the logarithmic sums.
    let tolerance = (log2_approx(large_bound as f64) + 8.0) * LOG_SCALE;
    let params = SieveParams {
        interval,
        length: 2 * interval + 1,
        large_bound,
        large_primes: config.large_primes,
        tolerance,
    };
    let length_u64 = u64::try_from(params.length).ok()?;
    let mut q_cursor = polynomial_prime_start(n, interval);
    let mut relations: Vec<Relation> = Vec::new();
    let mut partials: HashMap<u64, Relation> = HashMap::new();
    let extra = limits.extra_relations.max(EXTRA_DEPENDENCIES);
    let mut target = columns + extra;
    let max_relations = 3 * columns + 256;
    let max_partials = columns.saturating_mul(PARTIALS_PER_COLUMN);
    let position_budget = limits.max_sieve_positions;
    while stats.polynomials < limits.polynomials && stats.sieve_positions < position_budget {
        let work_left = position_budget - stats.sieve_positions;
        let by_work = usize::try_from(work_left / length_u64).unwrap_or(usize::MAX);
        if by_work == 0 {
            break;
        }
        let remaining = limits.polynomials - stats.polynomials;
        let batch = POLYNOMIAL_BATCH.min(remaining).min(by_work);
        let mut primes = Vec::with_capacity(batch);
        for _ in 0..batch {
            let q = match next_polynomial_prime(n, q_cursor) {
                | Ok(prime) => prime,
                | Err(factor) => return Some(factor),
            };
            q_cursor = q.saturating_add(2);
            primes.push(q);
        }
        // Batches are independent; `map_slice` returns them in index order,
        // so the relation stream (and the result) does not depend on the
        // number of worker threads.
        let outputs = crate::kernels::par::map_slice(&primes, 4, |&q| {
            sieve_polynomial(n, q, &factor_base, &divisors, params)
        });
        for output in outputs {
            if let Some(factor) = output.factor {
                return Some(factor);
            }
            if !output.polynomial {
                continue;
            }
            stats.polynomials += output.stats.polynomials;
            stats.candidates += output.stats.candidates;
            stats.full_relations += output.stats.full_relations;
            stats.partial_relations += output.stats.partial_relations;
            stats.sieve_positions += length_u64;
            for relation in output.full {
                relations.push(relation);
            }
            for (cofactor, relation) in output.partial {
                if let Some(previous) = partials.remove(&cofactor) {
                    stats.combined_relations += 1;
                    stats.full_relations += 1;
                    relations.push(combine(&previous, &relation, cofactor, n));
                } else if partials.len() < max_partials {
                    partials.insert(cofactor, relation);
                }
            }
            if relations.len() >= max_relations {
                break;
            }
        }
        if relations.len() >= target {
            if let Some(divisor) = solve_dependencies(&relations, &factor_base, n, stats) {
                return Some(divisor);
            }
            target = relations.len().saturating_add((columns / 16).max(extra));
            if target > max_relations {
                return None;
            }
        }
        if relations.len() >= max_relations {
            break;
        }
    }
    if relations.len() > columns {
        solve_dependencies(&relations, &factor_base, n, stats)
    } else {
        None
    }
}

/// Sieves one polynomial and extracts its relations.
fn sieve_polynomial(
    n: &BigInt,
    q: u64,
    factor_base: &[FbPrime],
    divisors: &[BigInt],
    params: SieveParams,
) -> SieveOutput {
    let mut output = SieveOutput::default();
    let Some(polynomial) = build_polynomial(n, q) else {
        return output;
    };
    let length = params.length;
    let interval = params.interval;
    let skip = factor_base.iter().position(|entry| entry.prime == q);
    let (roots1, roots2) = sieve_roots(&polynomial, factor_base, skip, divisors);
    // g(x) is even iff x and c have the same parity; convert that to an
    // index parity through x = index - interval.
    let x_parity = usize::from(polynomial.c.is_odd());
    let index_parity = (x_parity + interval % 2) % 2;
    let mut sieve = vec![0_u32; length];
    let two_log = scaled_log(2);
    for index in (index_parity..length).step_by(2) {
        sieve[index] += two_log;
    }
    for (index, entry) in factor_base.iter().enumerate().skip(1) {
        if Some(index) == skip {
            continue;
        }
        let prime = usize::try_from(entry.prime).unwrap_or(0);
        if prime == 0 {
            break;
        }
        let log = entry.log;
        let first = (usize::try_from(roots1[index]).unwrap_or(0) + interval) % prime;
        let mut position = first;
        while position < length {
            sieve[position] += log;
            position += prime;
        }
        let second = roots2[index];
        if second != roots1[index] {
            let mut position = (usize::try_from(second).unwrap_or(0) + interval) % prime;
            while position < length {
                sieve[position] += log;
                position += prime;
            }
        }
    }
    // Every `x` is considered. `g(x)` is even for exactly one parity (the
    // special prime 2 above); the other parity gives odd, still smoothable
    // values.
    let mut candidate = vec![false; length];
    let mut square_factor = None;
    for (index, log) in sieve.iter().enumerate() {
        let x = index as f64 - interval as f64;
        let value = (polynomial.a_f64 * x + polynomial.two_b_f64) * x + polynomial.c_f64;
        if value == 0.0 {
            // (a x + b)^2 = n exactly; this is a factor with no sieve.
            let ax_b = &polynomial.a * BigInt::from(index as i64 - interval as i64) + &polynomial.b;
            square_factor = Some(Integer::gcd(&ax_b, n));
            continue;
        }
        candidate[index] =
            f64::from(*log) + params.tolerance >= log2_approx(value.abs()) * LOG_SCALE;
    }
    if let Some(factor) = square_factor
        && factor > BigInt::one()
        && factor < *n
    {
        output.factor = Some(factor);
        return output;
    }
    if !candidate.iter().any(|&hit| hit) {
        output.polynomial = true;
        output.stats.polynomials = 1;
        return output;
    }
    let mut hits: Vec<(u32, u32)> = Vec::new();
    for (index, hit) in candidate.iter().enumerate() {
        if *hit && index % 2 == index_parity {
            hits.push((index as u32, 0));
        }
    }
    for (index, entry) in factor_base.iter().enumerate().skip(1) {
        if Some(index) == skip {
            continue;
        }
        let prime = usize::try_from(entry.prime).unwrap_or(0);
        if prime == 0 {
            break;
        }
        let mut positions = [roots1[index], roots2[index]];
        if positions[0] == positions[1] {
            positions[1] = u32::MAX;
        }
        for root in positions {
            if root == u32::MAX {
                continue;
            }
            let mut position = (usize::try_from(root).unwrap_or(0) + interval) % prime;
            while position < length {
                if candidate[position] {
                    hits.push((position as u32, index as u32));
                }
                position += prime;
            }
        }
    }
    hits.sort_unstable_by_key(|&(position, prime)| (position, prime));
    let known = &polynomial.q % n;
    let mut next_hit = 0_usize;
    for index in 0..length {
        if !candidate[index] {
            continue;
        }
        output.stats.candidates += 1;
        let x = index as i64 - interval as i64;
        let value = evaluate_g(&polynomial, x);
        let negative = value.is_negative();
        let mut magnitude = value.abs();
        let mut exponents: Vec<(u32, u32)> = Vec::new();
        while next_hit < hits.len()
            && usize::try_from(hits[next_hit].0).unwrap_or(usize::MAX) == index
        {
            let prime_index = usize::try_from(hits[next_hit].1).unwrap_or(usize::MAX);
            if prime_index >= divisors.len() {
                next_hit += 1;
                continue;
            }
            let mut exponent = 0_u32;
            loop {
                let (quotient, remainder) = magnitude.div_rem(&divisors[prime_index]);
                if !remainder.is_zero() {
                    break;
                }
                magnitude = quotient;
                exponent += 1;
                if magnitude.is_one() {
                    break;
                }
            }
            if exponent > 0 {
                exponents.push((hits[next_hit].1, exponent));
            }
            next_hit += 1;
        }
        let relation = Relation {
            x: mod_pos(&(&polynomial.a * BigInt::from(x) + &polynomial.b), n),
            known: known.clone(),
            negative,
            exponents,
        };
        if magnitude.is_one() {
            output.stats.full_relations += 1;
            output.full.push(relation);
        } else if let Some(cofactor) = magnitude
            .to_u64()
            .filter(|&value| value <= params.large_bound)
        {
            output.stats.partial_relations += 1;
            if params.large_primes {
                output.partial.push((cofactor, relation));
            }
        }
        if output.full.len() + output.partial.len() >= 4 * divisors.len() + 256 {
            break;
        }
    }
    output.polynomial = true;
    output.stats.polynomials = 1;
    output
}

/// Returns the smallest prime divisor `p <= 1000` of the odd `n > 1`.
fn small_prime_divisor(n: &BigInt) -> Option<u64> {
    let primes = primes_sieve(usize::try_from(SMALL_TRIAL).ok()?);
    for prime in primes.into_iter().skip(1) {
        let divisor = BigInt::from(prime);
        if (n % &divisor).is_zero() {
            return (divisor != *n).then_some(prime as u64);
        }
    }
    None
}

/// Decimal digit count of a positive integer.
fn decimal_digits(n: &BigInt) -> u32 {
    u32::try_from(n.to_str_radix(10).len()).unwrap_or(u32::MAX)
}

/// A Baillie-PSW probable-prime test.
///
/// Deterministic Miller-Rabin with the first twelve prime bases for
/// machine words (which is a proof for all of `u64`), and a base-2 strong
/// probable-prime test plus a strong Lucas test with Selfridge's
/// parameters for larger integers. No composite is known to pass the
/// latter, in particular the smallest strong pseudoprime to those twelve
/// Miller-Rabin bases (318665857834031151167461) does not.
#[must_use]
pub fn is_probable_prime(n: &BigInt) -> bool {
    if let Some(value) = n.to_u64() {
        return is_probable_prime_u64(value);
    }
    if !n.is_positive() || !n.bit(0) {
        return false;
    }
    for prime in [
        3_u64, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
        89, 97,
    ] {
        let prime = BigInt::from(prime);
        if (n % &prime).is_zero() {
            return false;
        }
    }
    strong_probable_prime(n, &BigInt::from(2)) && strong_lucas(n)
}

/// Strong probable-prime test of the odd `n > 2` to the base `a`.
fn strong_probable_prime(
    n: &BigInt,
    a: &BigInt,
) -> bool {
    let one = BigInt::one();
    let minus_one = n - &one;
    let shift = minus_one.trailing_zeros().unwrap_or(0);
    let odd = &minus_one >> shift;
    let mut power = a.modpow(&odd, n);
    if power == one || power == minus_one {
        return true;
    }
    for _ in 1..shift {
        power = &power * &power % n;
        if power == minus_one {
            return true;
        }
    }
    false
}

/// Strong Lucas probable-prime test of the odd `n > 2` with Selfridge's
/// parameters (`D` the first of `5, -7, 9, -11, ...` with Jacobi symbol
/// `-1`, `P = 1`, `Q = (1 - D)/4`): `U_d = 0` or `V_{d 2^r} = 0` for some
/// `r < s`, where `n + 1 = d 2^s` with `d` odd.
fn strong_lucas(n: &BigInt) -> bool {
    let root = n.sqrt();
    if &root * &root == *n {
        return false;
    }
    let mut d: i64 = 5;
    loop {
        match jacobi_big(&BigInt::from(d), n) {
            | -1 => break,
            // A common factor with D: composite unless n is D itself.
            | 0 if BigInt::from(d.abs()) != *n => return false,
            | _ => {},
        }
        d = if d > 0 { -(d + 2) } else { -d + 2 };
    }
    let (d_big, q) = (BigInt::from(d), BigInt::from((1 - d) / 4));
    let n1: BigInt = n + 1;
    let shift = n1.trailing_zeros().unwrap_or(0);
    let odd = &n1 >> shift;
    // Halving modulo the odd n.
    let half = |x: BigInt| -> BigInt {
        let x = mod_pos(&x, n);
        if x.bit(0) {
            (x + n) >> 1
        } else {
            x >> 1
        }
    };
    let (mut u, mut v, mut qk) = (BigInt::one(), BigInt::one(), mod_pos(&q, n));
    for i in (0..odd.bits().saturating_sub(1)).rev() {
        u = &u * &v % n;
        v = mod_pos(&(&v * &v - 2 * &qk), n);
        qk = &qk * &qk % n;
        if odd.bit(i) {
            let next_u = half(&u + &v);
            let next_v = half(&d_big * &u + &v);
            u = next_u;
            v = next_v;
            qk = mod_pos(&(&qk * &q), n);
        }
    }
    if u.is_zero() || v.is_zero() {
        return true;
    }
    for _ in 1..shift {
        v = mod_pos(&(&v * &v - 2 * &qk), n);
        qk = &qk * &qk % n;
        if v.is_zero() {
            return true;
        }
    }
    false
}

/// The Jacobi symbol `(a/n)` for odd positive `n`, on arbitrary-precision
/// integers (used only by the strong Lucas test above).
fn jacobi_big(
    a: &BigInt,
    n: &BigInt,
) -> i32 {
    if !n.is_positive() || !n.bit(0) {
        return 0;
    }
    let mut a = mod_pos(a, n);
    let mut n = n.clone();
    let mut result = 1_i32;
    while !a.is_zero() {
        while !a.bit(0) {
            a >>= 1;
            let residue = (&n % 8_u32).to_u32().unwrap_or(0);
            if residue == 3 || residue == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        let a_mod4 = (&a % 4_u32).to_u32().unwrap_or(0);
        let n_mod4 = (&n % 4_u32).to_u32().unwrap_or(0);
        if a_mod4 == 3 && n_mod4 == 3 {
            result = -result;
        }
        a = mod_pos(&a, &n);
    }
    if n.is_one() { result } else { 0 }
}

/// Deterministic Miller-Rabin for a machine word.
fn is_probable_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for prime in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == prime {
            return true;
        }
        if n.is_multiple_of(prime) {
            return false;
        }
    }
    let shift = (n - 1).trailing_zeros();
    let odd = (n - 1) >> shift;
    'witness: for base in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut power = powmod(base, odd, n);
        if power == 1 || power == n - 1 {
            continue;
        }
        for _ in 1..shift {
            power = mulmod(power, power, n);
            if power == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// A root `r` of `r^k = n` for some `k >= 2`, if `n` is a perfect power.
///
/// Only prime exponents below 64 need checking: a perfect power also has
/// a prime exponent, and with at most 80 digits a larger exponent forces
/// a base below 21, which trial division has already removed.
fn perfect_power_root(n: &BigInt) -> Option<BigInt> {
    const EXPONENTS: [u32; 18] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    ];
    for exponent in EXPONENTS {
        if u64::from(exponent) > n.bits() {
            break;
        }
        let root = n.nth_root(exponent);
        if root > BigInt::one() && root.pow(exponent) == *n {
            return Some(root);
        }
    }
    None
}

/// The factor base entry for one sieving prime.
#[derive(Copy, Clone, Debug)]
struct FbPrime {
    /// The prime.
    prime: u64,
    /// A square root of `n` modulo `prime` (unused for `prime = 2`).
    sqrt_n: u64,
    /// `log2(prime)` scaled by [`LOG_SCALE`].
    log: u32,
}

/// Primes `p <= bound` that can divide `g(x)`, plus the special prime 2.
///
/// Returns `Err(p)` if `p` divides `n`, which is a factor of `n`.
fn build_factor_base(
    n: &BigInt,
    bound: u64,
) -> Result<Vec<FbPrime>, BigInt> {
    let limit = usize::try_from(bound).map_err(|_| BigInt::one())?;
    let mut base = Vec::new();
    base.push(FbPrime {
        prime: 2,
        sqrt_n: 1,
        log: scaled_log(2),
    });
    for prime in primes_sieve(limit).into_iter().skip(1) {
        let prime = u64::try_from(prime).map_err(|_| BigInt::one())?;
        let big = BigInt::from(prime);
        let residue = mod_pos(n, &big);
        let residue_u64 = residue.to_u64().ok_or_else(|| big.clone())?;
        if residue_u64 == 0 {
            return Err(big);
        }
        if jacobi(residue_u64, prime) != 1 {
            continue;
        }
        let Some(root) = sqrt_mod_prime(residue_u64, prime) else {
            continue;
        };
        base.push(FbPrime {
            prime,
            sqrt_n: root,
            log: scaled_log(prime),
        });
    }
    Ok(base)
}

/// `log2(p)` scaled and rounded to the sieve unit.
fn scaled_log(p: u64) -> u32 {
    (log2_approx(p as f64) * LOG_SCALE).round() as u32
}

/// A polynomial `a x^2 + 2 b x + c = ((a x + b)^2 - n) / a`, `a = q^2`.
struct Polynomial {
    /// `a = q^2`.
    a: BigInt,
    /// `b` with `b^2 = n (mod a)` and `|b| <= a/2`.
    b: BigInt,
    /// `2 b`, the linear coefficient.
    two_b: BigInt,
    /// `c = (b^2 - n) / a`.
    c: BigInt,
    /// The polynomial prime `q`.
    q: BigInt,
    /// `a` as a float for the logarithmic target.
    a_f64: f64,
    /// `2 b` as a float.
    two_b_f64: f64,
    /// `c` as a float.
    c_f64: f64,
}

/// Builds the polynomial for the prime `q` (Tonelli-Shanks plus one
/// Hensel lift).
fn build_polynomial(
    n: &BigInt,
    q: u64,
) -> Option<Polynomial> {
    let prime = BigInt::from(q);
    let residue = mod_pos(n, &prime).to_u64()?;
    let root = sqrt_mod_prime(residue, q)?;
    if root == 0 {
        return None;
    }
    let a = &prime * &prime;
    let root_big = BigInt::from(root);
    let difference = (n - &root_big * &root_big) / &prime;
    let difference_mod = mod_pos(&difference, &prime).to_u64()?;
    let inverse = inverse(2 * root % q, q)?;
    let lift = mulmod(difference_mod, inverse, q);
    let unsigned = &root_big + &prime * BigInt::from(lift);
    let half = &a >> 1;
    let b = if unsigned > half {
        &unsigned - &a
    } else {
        unsigned
    };
    let c = (&b * &b - n) / &a;
    let a_f64 = a.to_f64()?;
    let two_b: BigInt = &b * 2;
    let two_b_f64 = two_b.to_f64()?;
    let c_f64 = c.to_f64()?;
    Some(Polynomial {
        a,
        b,
        two_b,
        c,
        q: prime,
        a_f64,
        two_b_f64,
        c_f64,
    })
}

/// The first prime `q` to try: the integer square root of
/// `sqrt(2 n) / M`.
fn polynomial_prime_start(
    n: &BigInt,
    interval: usize,
) -> u64 {
    let doubled: BigInt = n * 2;
    let target = doubled.sqrt() / BigInt::from(interval);
    let start = target.sqrt().to_u64().unwrap_or(u64::MAX);
    start.max(3)
}

/// The next prime `q >= start` with `(n/q) = 1` and `q` not dividing `n`.
///
/// `Err(p)` means `p` divides `n` and is a factor.
fn next_polynomial_prime(
    n: &BigInt,
    start: u64,
) -> Result<u64, BigInt> {
    let mut candidate = start.max(3) | 1;
    for _ in 0..1_000_000 {
        let big = BigInt::from(candidate);
        let residue = (n % &big).to_u64().unwrap_or(0);
        if residue == 0 {
            return Err(big);
        }
        if is_probable_prime_u64(candidate) && jacobi(residue, candidate) == 1 {
            return Ok(candidate);
        }
        candidate = candidate.checked_add(2).ok_or_else(BigInt::one)?;
    }
    Err(BigInt::one())
}

/// The two roots of `g(x) = 0 (mod p)` for every factor-base prime.
///
/// `roots2` equals `roots1` for a double root. The prime `q` that defines
/// `a` is skipped: it divides `a`, is accounted for separately through
/// the square `a = q^2`, and never takes part in the smooth part.
fn sieve_roots(
    polynomial: &Polynomial,
    base: &[FbPrime],
    skip: Option<usize>,
    divisors: &[BigInt],
) -> (Vec<u32>, Vec<u32>) {
    let mut roots1 = vec![0_u32; base.len()];
    let mut roots2 = vec![0_u32; base.len()];
    let q = polynomial.q.to_u64().unwrap_or(0);
    for (index, entry) in base.iter().enumerate().skip(1) {
        if Some(index) == skip {
            roots2[index] = roots1[index];
            continue;
        }
        let prime = entry.prime;
        let inverse_q = inverse(q % prime, prime).unwrap_or(0);
        let inverse_a = mulmod(inverse_q, inverse_q, prime);
        let b_mod = mod_pos(&polynomial.b, &divisors[index])
            .to_u64()
            .unwrap_or(0)
            % prime;
        let t_mod = entry.sqrt_n % prime;
        let first = mulmod(mod_sub(t_mod, b_mod, prime), inverse_a, prime);
        let second = mulmod(mod_sub(prime - t_mod, b_mod, prime), inverse_a, prime);
        roots1[index] = u32::try_from(first).unwrap_or(0);
        roots2[index] = u32::try_from(second).unwrap_or(0);
    }
    (roots1, roots2)
}

/// `a - b mod m` for `0 <= a, b < m`.
const fn mod_sub(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    (a + m - b) % m
}

/// `g(x)` exactly.
fn evaluate_g(
    polynomial: &Polynomial,
    x: i64,
) -> BigInt {
    let x = BigInt::from(x);
    &polynomial.a * &x * &x + &polynomial.two_b * &x + &polynomial.c
}

/// One relation: `X = prod (a x + b)` and the factorisation of `g`.
#[derive(Clone, Debug)]
struct Relation {
    /// `(a x + b) mod n`, the term whose square is congruent to `Q(x)`.
    x: BigInt,
    /// The product of the polynomial primes `q`, the known square root
    /// part of `Q` outside the factor base.
    known: BigInt,
    /// Whether `g(x)` is negative (the `-1` exponent is odd).
    negative: bool,
    /// `(factor-base index, exponent)` pairs, sorted by index.
    exponents: Vec<(u32, u32)>,
}

/// Combines two partial relations with the same cofactor `cofactor`.
fn combine(
    first: &Relation,
    second: &Relation,
    cofactor: u64,
    n: &BigInt,
) -> Relation {
    let mut exponents: Vec<(u32, u32)> = first
        .exponents
        .iter()
        .chain(second.exponents.iter())
        .copied()
        .collect();
    exponents.sort_unstable_by_key(|&(index, _)| index);
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(exponents.len());
    for (index, exponent) in exponents {
        match merged.last_mut() {
            | Some((last, total)) if *last == index => *total += exponent,
            | _ => merged.push((index, exponent)),
        }
    }
    Relation {
        x: &first.x * &second.x % n,
        known: &first.known * &second.known % n * BigInt::from(cofactor) % n,
        negative: first.negative ^ second.negative,
        exponents: merged,
    }
}

/// One row of the GF(2) matrix: the parity of `g`'s factorisation.
fn parity_row(
    relation: &Relation,
    columns: usize,
) -> Vec<u64> {
    let mut row = vec![0_u64; columns.div_ceil(64)];
    if relation.negative {
        row[0] |= 1;
    }
    for &(index, exponent) in &relation.exponents {
        if exponent % 2 == 1 {
            let column = index as usize + 1;
            row[column / 64] |= 1_u64 << (column % 64);
        }
    }
    row
}

/// Highest set bit of a bitset, if any.
fn highest_bit(bits: &[u64]) -> Option<usize> {
    for (index, word) in bits.iter().enumerate().rev() {
        if *word != 0 {
            return Some(index * 64 + 63 - word.leading_zeros() as usize);
        }
    }
    None
}

/// Extracts the set bits of a dense row as sorted column indices.
fn sparse_columns(
    row: &[u64],
    columns: usize,
) -> Vec<u32> {
    let mut out = Vec::new();
    for (word_index, word) in row.iter().enumerate() {
        let mut bits = *word;
        while bits != 0 {
            let bit = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            let column = word_index * 64 + bit;
            if column < columns {
                out.push(u32::try_from(column).unwrap_or(u32::MAX));
            }
        }
    }
    out
}

/// Xor of two sorted column lists.
fn xor_sorted(
    first: &[u32],
    second: &[u32],
) -> Vec<u32> {
    let mut out = Vec::with_capacity(first.len() + second.len());
    let (mut i, mut j) = (0_usize, 0_usize);
    while i < first.len() && j < second.len() {
        match first[i].cmp(&second[j]) {
            | std::cmp::Ordering::Less => {
                out.push(first[i]);
                i += 1;
            },
            | std::cmp::Ordering::Greater => {
                out.push(second[j]);
                j += 1;
            },
            | std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            },
        }
    }
    out.extend_from_slice(&first[i..]);
    out.extend_from_slice(&second[j..]);
    out
}

/// Xors `source` into `target` (same length).
fn xor_bits(
    target: &mut [u64],
    source: &[u64],
) {
    for (left, right) in target.iter_mut().zip(source) {
        *left ^= right;
    }
}

/// The combination bitset of a row.
///
/// Original rows are kept implicit (only their own bit is set) until they
/// take part in a merge; this keeps the memory proportional to the number
/// of merged rows rather than to `rows * rows / 64`.
fn materialize_combo(
    row: usize,
    combos: &mut [Option<Vec<u64>>],
    combo_words: usize,
) -> Vec<u64> {
    match combos.get_mut(row).and_then(Option::take) {
        | Some(bits) => bits,
        | None => {
            let mut bits = vec![0_u64; combo_words];
            bits[row / 64] |= 1_u64 << (row % 64);
            bits
        },
    }
}

/// Removes a row from all of its columns and marks it dead.
fn detach_row(
    row: usize,
    row_cols: &[Vec<u32>],
    col_rows: &mut [Vec<u32>],
    weights: &mut [usize],
    heap: &mut BinaryHeap<Reverse<(usize, usize)>>,
) {
    let row_id = u32::try_from(row).unwrap_or(u32::MAX);
    for &column in &row_cols[row] {
        let column = column as usize;
        if let Some(position) = col_rows[column].iter().position(|&value| value == row_id) {
            col_rows[column].swap_remove(position);
        }
        weights[column] = weights[column].saturating_sub(1);
        heap.push(Reverse((weights[column], column)));
    }
}

/// Adds a row to all of its columns.
fn attach_row(
    row: usize,
    columns: &[u32],
    col_rows: &mut [Vec<u32>],
    weights: &mut [usize],
    heap: &mut BinaryHeap<Reverse<(usize, usize)>>,
) {
    let row_id = u32::try_from(row).unwrap_or(u32::MAX);
    for &column in columns {
        let column = column as usize;
        col_rows[column].push(row_id);
        weights[column] += 1;
        heap.push(Reverse((weights[column], column)));
    }
}

/// Rows whose xor is zero, as bitsets over the row indices.
///
/// Structured elimination processes columns from the lightest to the
/// heaviest. A singleton column deletes its row (that row cannot be part of
/// a dependency); a heavier column xors its shortest row into every other
/// row containing it, which removes the column and discards the pivot row.
/// The row operations are reversible and each row carries the xor of the
/// original relations it is a combination of, so a row that becomes empty
/// yields a dependency directly. Once the remaining core is small enough
/// (or fill-in makes the sparse phase unprofitable) the core is eliminated
/// densely with the same combination tracking.
fn find_dependencies(
    rows: &[Vec<u64>],
    columns: usize,
) -> Vec<Vec<u64>> {
    let row_count = rows.len();
    let combo_words = row_count.div_ceil(64);
    let mut row_cols: Vec<Vec<u32>> = rows
        .iter()
        .map(|row| sparse_columns(row, columns))
        .collect();
    let mut alive = vec![true; row_count];
    let mut combos: Vec<Option<Vec<u64>>> = vec![None; row_count];
    let mut col_rows: Vec<Vec<u32>> = vec![Vec::new(); columns];
    for (row, cols) in row_cols.iter().enumerate() {
        let row_id = u32::try_from(row).unwrap_or(u32::MAX);
        for &column in cols {
            col_rows[column as usize].push(row_id);
        }
    }
    let mut weights: Vec<usize> = col_rows.iter().map(Vec::len).collect();
    let mut heap: BinaryHeap<Reverse<(usize, usize)>> = weights
        .iter()
        .enumerate()
        .map(|(column, &weight)| Reverse((weight, column)))
        .collect();
    let mut dependencies: Vec<Vec<u64>> = Vec::new();
    let mut alive_count = row_count;
    // A relation that is already a square is a dependency on its own.
    for row in 0..row_count {
        if row_cols[row].is_empty() {
            let combo = materialize_combo(row, &mut combos, combo_words);
            dependencies.push(combo);
            alive[row] = false;
            alive_count -= 1;
        }
    }
    if dependencies.len() >= MAX_DEPENDENCIES {
        return dependencies;
    }
    let mut steps = 0_usize;
    let step_budget = row_count.saturating_mul(8).saturating_add(columns);
    while let Some(Reverse((weight, column))) = heap.pop() {
        if dependencies.len() >= MAX_DEPENDENCIES {
            return dependencies;
        }
        steps += 1;
        if steps > step_budget || alive_count == 0 || alive_count <= DENSE_CORE_TARGET {
            break;
        }
        if weights[column] != weight {
            continue;
        }
        if weight == 0 {
            continue;
        }
        if weight > SPARSE_MAX_WEIGHT {
            // The heap is ordered, so every remaining column is heavier.
            break;
        }
        let active_rows: Vec<usize> = col_rows[column]
            .iter()
            .filter_map(|&row| {
                let row = row as usize;
                (alive.get(row).copied().unwrap_or(false) && !row_cols[row].is_empty())
                    .then_some(row)
            })
            .collect();
        if active_rows.len() != weight {
            // Stale adjacency; the weight will be refreshed by a pending
            // heap entry.
            continue;
        }
        if weight == 1 {
            let row = active_rows[0];
            alive[row] = false;
            alive_count -= 1;
            detach_row(row, &row_cols, &mut col_rows, &mut weights, &mut heap);
            row_cols[row].clear();
            continue;
        }
        // Markowitz-style pivot: the shortest row keeps fill-in small.
        let pivot_row = active_rows
            .iter()
            .copied()
            .min_by_key(|&row| (row_cols[row].len(), row))
            .unwrap_or(active_rows[0]);
        if row_cols[pivot_row].len() > SPARSE_MAX_ROW {
            break;
        }
        let pivot_columns = row_cols[pivot_row].clone();
        let pivot_combo = materialize_combo(pivot_row, &mut combos, combo_words);
        for &row in &active_rows {
            if row == pivot_row {
                continue;
            }
            let mut combined = materialize_combo(row, &mut combos, combo_words);
            xor_bits(&mut combined, &pivot_combo);
            let new_cols = xor_sorted(&row_cols[row], &pivot_columns);
            detach_row(row, &row_cols, &mut col_rows, &mut weights, &mut heap);
            if new_cols.is_empty() {
                alive[row] = false;
                alive_count -= 1;
                row_cols[row].clear();
                dependencies.push(combined);
                continue;
            }
            row_cols[row] = new_cols;
            combos[row] = Some(combined);
            attach_row(row, &row_cols[row], &mut col_rows, &mut weights, &mut heap);
        }
        // The pivot column appears only in the pivot row now; that row
        // cannot belong to a dependency and is discarded.
        detach_row(pivot_row, &row_cols, &mut col_rows, &mut weights, &mut heap);
        alive[pivot_row] = false;
        alive_count -= 1;
        row_cols[pivot_row].clear();
    }
    if dependencies.len() >= MAX_DEPENDENCIES {
        return dependencies;
    }
    // Dense elimination of the remaining core.
    let mut column_map = vec![usize::MAX; columns];
    let mut core_columns = 0_usize;
    for row in 0..row_count {
        if !alive[row] {
            continue;
        }
        for &column in &row_cols[row] {
            let column = column as usize;
            if column_map[column] == usize::MAX {
                column_map[column] = core_columns;
                core_columns += 1;
            }
        }
    }
    let words = core_columns.div_ceil(64);
    let mut core: Vec<(Vec<u64>, Vec<u64>)> = Vec::new();
    for row in 0..row_count {
        if !alive[row] {
            continue;
        }
        let mut dense = vec![0_u64; words];
        for &column in &row_cols[row] {
            let mapped = column_map[column as usize];
            if mapped != usize::MAX {
                dense[mapped / 64] |= 1_u64 << (mapped % 64);
            }
        }
        let combo = materialize_combo(row, &mut combos, combo_words);
        core.push((dense, combo));
    }
    let mut basis: Vec<Option<(Vec<u64>, Vec<u64>)>> = vec![None; core_columns];
    for (mut reduced, mut combo) in core {
        if dependencies.len() >= MAX_DEPENDENCIES {
            break;
        }
        loop {
            match highest_bit(&reduced) {
                | Some(pivot) => {
                    match basis[pivot].take() {
                        | Some((row_bits, row_combo)) => {
                            xor_bits(&mut reduced, &row_bits);
                            xor_bits(&mut combo, &row_combo);
                            basis[pivot] = Some((row_bits, row_combo));
                        },
                        | None => {
                            basis[pivot] = Some((reduced, combo));
                            break;
                        },
                    }
                },
                | None => {
                    dependencies.push(combo);
                    break;
                },
            }
        }
    }
    dependencies
}

/// Extracts a factor from one dependency, if the congruence is non-trivial.
fn dependency_factor(
    combo: &[u64],
    relations: &[Relation],
    base: &[FbPrime],
    n: &BigInt,
) -> Option<BigInt> {
    let mut x = BigInt::one();
    let mut known = BigInt::one();
    let mut negative = false;
    let mut exponents = vec![0_u32; base.len()];
    for (word_index, word) in combo.iter().enumerate() {
        let mut bits = *word;
        while bits != 0 {
            let bit = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            let row = word_index * 64 + bit;
            let relation = relations.get(row)?;
            x = x * &relation.x % n;
            known = known * &relation.known % n;
            negative ^= relation.negative;
            for &(index, exponent) in &relation.exponents {
                exponents[index as usize] += exponent;
            }
        }
    }
    if negative {
        return None;
    }
    let mut y = known;
    for (index, exponent) in exponents.iter().enumerate() {
        if exponent % 2 == 1 {
            return None;
        }
        if *exponent > 0 {
            let prime = BigInt::from(base[index].prime);
            y = y * prime.modpow(&BigInt::from(*exponent / 2), n) % n;
        }
    }
    let one = BigInt::one();
    [Integer::gcd(&(&x - &y), n), Integer::gcd(&(&x + &y), n)]
        .into_iter()
        .find(|candidate| candidate > &one && candidate < n)
}

/// Solves the linear system and tries the dependencies in order.
fn solve_dependencies(
    relations: &[Relation],
    base: &[FbPrime],
    n: &BigInt,
    stats: &mut QsStats,
) -> Option<BigInt> {
    let columns = base.len() + 1;
    let rows: Vec<Vec<u64>> = relations
        .iter()
        .map(|relation| parity_row(relation, columns))
        .collect();
    stats.matrix_rows = rows.len();
    stats.matrix_columns = columns;
    let dependencies = find_dependencies(&rows, columns);
    stats.dependencies = dependencies.len();
    for combo in dependencies {
        if let Some(divisor) = dependency_factor(&combo, relations, base, n) {
            return Some(divisor);
        }
    }
    None
}

/// `x mod m` in `[0, m)`.
fn mod_pos(
    x: &BigInt,
    m: &BigInt,
) -> BigInt {
    let remainder = x % m;
    if remainder.is_negative() {
        remainder + m
    } else {
        remainder
    }
}

/// The Jacobi symbol `(a/n)` for odd `n`.
fn jacobi(
    mut a: u64,
    mut n: u64,
) -> i32 {
    if n == 0 {
        return i32::from(a == 1);
    }
    a %= n;
    let mut result = 1_i32;
    while a != 0 {
        while a.is_multiple_of(2) {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a %= n;
    }
    if n == 1 { result } else { 0 }
}

/// `a * b mod m` without overflow (`m < 2^64`).
fn mulmod(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(m)) as u64
}

/// `base^exponent mod m`.
fn powmod(
    mut base: u64,
    mut exponent: u64,
    m: u64,
) -> u64 {
    let mut accumulator = 1 % m;
    base %= m;
    while exponent > 0 {
        if exponent & 1 == 1 {
            accumulator = mulmod(accumulator, base, m);
        }
        base = mulmod(base, base, m);
        exponent >>= 1;
    }
    accumulator
}

/// The inverse of `a` modulo the odd prime `m`.
fn inverse(
    a: u64,
    m: u64,
) -> Option<u64> {
    let a = a % m;
    if a == 0 {
        return None;
    }
    let modulus = i128::from(m);
    let (mut old_r, mut r) = (modulus, i128::from(a));
    let (mut old_t, mut t) = (0_i128, 1_i128);
    while r != 0 {
        let quotient = old_r / r;
        let next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;
        let next_t = old_t - quotient * t;
        old_t = t;
        t = next_t;
    }
    if old_r != 1 {
        return None;
    }
    Some((((old_t % modulus) + modulus) % modulus) as u64)
}

/// Tonelli-Shanks: a square root of `a` modulo the odd prime `p`.
fn sqrt_mod_prime(
    a: u64,
    p: u64,
) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if p == 2 {
        return Some(a);
    }
    if powmod(a, (p - 1) / 2, p) != 1 {
        return None;
    }
    let mut odd = p - 1;
    let mut shift = 0_u64;
    while odd.is_multiple_of(2) {
        odd /= 2;
        shift += 1;
    }
    let mut non_residue = 2_u64;
    while jacobi(non_residue, p) != -1 {
        non_residue += 1;
        if non_residue >= p {
            return None;
        }
    }
    let mut m = shift;
    let mut c = powmod(non_residue, odd, p);
    let mut t = powmod(a, odd, p);
    let mut root = powmod(a, odd.div_ceil(2), p);
    while t != 1 {
        let mut inner = 0_u64;
        let mut probe = t;
        while probe != 1 {
            probe = mulmod(probe, probe, p);
            inner += 1;
            if inner >= m {
                return None;
            }
        }
        let exponent = 1_u64.checked_shl(u32::try_from(m - inner - 1).ok()?)?;
        let factor = powmod(c, exponent, p);
        root = mulmod(root, factor, p);
        c = mulmod(factor, factor, p);
        t = mulmod(t, c, p);
        m = inner;
    }
    Some(root)
}

/// `log2(x)` to about `1e-7`, using the exponent bits and an `atanh`
/// series on the mantissa (no call into the platform logarithm).
fn log2_approx(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return 0.0;
    }
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i64 - 1023;
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let z = (mantissa - 1.0) / (mantissa + 1.0);
    let z2 = z * z;
    let series = z2
        * (1.0 / 3.0
            + z2 * (1.0 / 5.0
                + z2 * (1.0 / 7.0 + z2 * (1.0 / 9.0 + z2 * (1.0 / 11.0 + z2 / 13.0)))));
    let logarithm = 2.0 * z * (1.0 + series);
    exponent as f64 + logarithm * std::f64::consts::LOG2_E
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(s: &str) -> BigInt {
        s.parse().unwrap_or_else(|_| BigInt::zero())
    }

    #[test]
    fn jacobi_matches_a_direct_definition() {
        for p in [3_u64, 5, 7, 11, 13, 17, 19, 97] {
            for a in 1..p {
                let direct = (1..p).any(|x| (x * x) % p == a);
                let symbol = jacobi(a, p);
                assert_eq!(symbol == 1, direct, "({a}/{p})");
                assert_eq!(symbol == 0, a % p == 0, "({a}/{p})");
            }
        }
    }

    #[test]
    fn square_roots_are_roots() {
        for p in [3_u64, 5, 7, 11, 13, 17, 101, 65_537] {
            for a in 1..p.min(30) {
                if jacobi(a, p) != 1 {
                    continue;
                }
                let root = sqrt_mod_prime(a, p).expect("a residue has a root");
                assert_eq!(mulmod(root, root, p), a % p, "sqrt({a} mod {p})");
            }
        }
    }

    #[test]
    fn inverse_round_trips() {
        for p in [3_u64, 5, 7, 101, 65_537] {
            for a in 1..p.min(50) {
                let inverse = inverse(a, p).expect("coprime");
                assert_eq!(mulmod(a, inverse, p), 1, "1/{a} mod {p}");
            }
        }
    }

    #[test]
    fn baillie_psw_rejects_the_twelve_base_pseudoprime() {
        // The smallest strong pseudoprime to the first twelve prime bases;
        // the old twelve-base Miller-Rabin test accepted it as prime.
        let pseudoprime = n("318665857834031151167461");
        assert!(!is_probable_prime(&pseudoprime));
        assert!(is_probable_prime(&n(
            "170141183460469231731687303715884105727"
        )));
        assert!(is_probable_prime(&n("1000003")));
        assert!(!is_probable_prime(&(n("1000003") * n("1000033"))));
        assert!(!is_probable_prime(&(n("1000003") * n("1000003"))));
    }

    #[test]
    fn polynomial_and_roots_satisfy_the_defining_congruences() {
        let n = n("1000003") * n("1000033");
        let base = build_factor_base(&n, 500).expect("no factor in the base");
        let divisors: Vec<BigInt> = base.iter().map(|p| BigInt::from(p.prime)).collect();
        let q = next_polynomial_prime(&n, polynomial_prime_start(&n, 1_000)).expect("a prime");
        let polynomial = build_polynomial(&n, q).expect("a polynomial");
        // a = q^2, b^2 = n (mod a), c = (b^2 - n)/a.
        assert_eq!(polynomial.a, BigInt::from(q) * BigInt::from(q));
        let b_squared = &polynomial.b * &polynomial.b;
        assert!((&b_squared - &n).is_multiple_of(&polynomial.a));
        assert_eq!(&polynomial.c * &polynomial.a, &b_squared - &n);
        assert_eq!(polynomial.two_b, &polynomial.b * 2);
        // (a x + b)^2 - n = a g(x) for the sieved x values.
        for x in -200_i64..=200 {
            let ax_b = &polynomial.a * BigInt::from(x) + &polynomial.b;
            assert_eq!(
                &ax_b * &ax_b - &n,
                &polynomial.a * evaluate_g(&polynomial, x)
            );
        }
        // Every stored root is a root of g modulo its factor-base prime.
        let skip = base.iter().position(|entry| entry.prime == q);
        let (roots1, roots2) = sieve_roots(&polynomial, &base, skip, &divisors);
        for (index, entry) in base.iter().enumerate().skip(1) {
            if Some(index) == skip {
                continue;
            }
            let modulus = BigInt::from(entry.prime);
            for root in [roots1[index], roots2[index]] {
                let value = evaluate_g(&polynomial, i64::from(root));
                assert!(
                    mod_pos(&value, &modulus).is_zero(),
                    "p = {}, root = {root}",
                    entry.prime
                );
            }
        }
    }

    #[test]
    fn combining_partials_squares_the_shared_cofactor() {
        let n = BigInt::from(1_000_003_u64);
        let first = Relation {
            x: BigInt::from(11),
            known: BigInt::from(3),
            negative: true,
            exponents: vec![(0, 1), (2, 3)],
        };
        let second = Relation {
            x: BigInt::from(13),
            known: BigInt::from(7),
            negative: true,
            exponents: vec![(2, 1), (5, 2)],
        };
        let combined = combine(&first, &second, 15, &n);
        assert_eq!(combined.x, BigInt::from(11 * 13));
        assert_eq!(combined.known, BigInt::from(3 * 7 * 15));
        assert!(!combined.negative);
        assert_eq!(combined.exponents, vec![(0, 1), (2, 4), (5, 2)]);
    }

    #[test]
    fn log2_approximation_is_accurate() {
        for value in [0.5_f64, 1.0, 1.5, 2.0, 3.0, 10.0, 1234.5, 1e18, 1e30] {
            let error = (log2_approx(value) - value.log2()).abs();
            assert!(error < 1e-6, "{value}: {error}");
        }
    }

    #[test]
    fn dependencies_are_valid_over_gf_two() {
        // A small deterministic generator.
        let mut state = 0x1234_5678_9abc_def0_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state
        };
        let columns = 37;
        let words = (columns + 63) / 64;
        let mask = (1_u64 << (columns % 64)) - 1;
        let mut rows: Vec<Vec<u64>> = (0..20)
            .map(|_| (0..words).map(|_| next() & mask).collect())
            .collect();
        // Row 19 is the xor of rows 2, 5, 13, and row 18 of 7 and 11:
        // known dependencies.
        rows.push(
            (0..words)
                .map(|w| rows[2][w] ^ rows[5][w] ^ rows[13][w])
                .collect(),
        );
        rows.push((0..words).map(|w| rows[7][w] ^ rows[11][w]).collect());
        let dependencies = find_dependencies(&rows, columns);
        assert!(!dependencies.is_empty());
        for combo in &dependencies {
            let mut xor = vec![0_u64; words];
            for (word_index, word) in combo.iter().enumerate() {
                let mut bits = *word;
                while bits != 0 {
                    let bit = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    for (target, source) in xor.iter_mut().zip(&rows[word_index * 64 + bit]) {
                        *target ^= source;
                    }
                }
            }
            assert!(xor.iter().all(|&word| word == 0));
        }
    }

    #[test]
    fn limits_scale_up_to_eighty_digits() {
        assert_eq!(QsLimits::MAX_DIGITS, 80);
        let forty = QsLimits::for_digits(40);
        let fifty_five = QsLimits::for_digits(55);
        let sixty = QsLimits::for_digits(60);
        let seventy = QsLimits::for_digits(70);
        let eighty = QsLimits::for_digits(80);
        assert!(forty.bound < fifty_five.bound);
        assert!(fifty_five.bound < sixty.bound);
        assert!(sixty.bound < seventy.bound);
        assert!(seventy.bound < eighty.bound);
        assert!(fifty_five.interval < sixty.interval);
        assert!(sixty.interval < seventy.interval);
        assert!(fifty_five.polynomials > 0);
        assert!(sixty.polynomials > 0);
        assert!(eighty.polynomials > 0);
        assert!(fifty_five.large_prime_factor < fifty_five.bound);
        assert!(sixty.large_prime_factor < sixty.bound);
        assert!(sixty.extra_relations >= EXTRA_DEPENDENCIES);
        assert!(eighty.max_sieve_positions <= 20_000_000_000);
        // The hard ceiling applies before any factorization work.
        let oversized = BigInt::from(10_u8).pow(81);
        assert!(qs_factor(&oversized, eighty).is_none());
    }

    #[test]
    fn small_composites_factor() {
        for text in [
            "1000003", "1000033", "1000037", "1000039", "1000081", "1000099", "1000117",
        ] {
            let prime = n(text);
            if is_probable_prime(&prime) {
                assert!(qs_factor(&prime, QsLimits::for_digits(8)).is_none());
            }
        }
        let n1 = n("1000003");
        let n2 = n("1000033");
        let composite = &n1 * &n2;
        let divisor = qs_factor(&composite, QsLimits::for_digits(14)).expect("semiprime");
        assert!(&composite % &divisor == BigInt::zero());
        assert!(divisor > BigInt::one() && divisor < composite);
    }
}
