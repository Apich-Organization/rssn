//! The quadratic sieve and its self-initialising variant (MPQS).
//!
//! `qs_factor(n, limits)` finds a non-trivial divisor of an odd composite
//! `n` that has no small prime factor: a factor base of primes with
//! `(n/p) = 1`, sieving of `Q(x) = (ax + b)^2 - n` over a symmetric
//! interval, a single-large-prime variation, and a GF(2) dependency solve
//! by sparse elimination. The rule set's `factorint` uses it for
//! composites that survive trial division, perfect-power detection, BPSW
//! and ECM.
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
//! the sign of `g` is the `-1` exponent of the relation. Candidates whose
//! accumulated `log p` is close to `log |g(x)|` are
//! factored exactly over the factor base; a cofactor below
//! `16 * bound` becomes a partial relation and is combined with the next
//! relation that has the same cofactor (treated as one large prime, which
//! is valid even when the cofactor is composite). Full relations are
//! reduced over GF(2) with a bitset elimination that tracks which rows
//! combine into each dependency; for every dependency the two square
//! roots `X = prod (a x + b)` and `Y = prod q * prod p^(e/2)` are formed
//! and `gcd(X - Y, n)`, `gcd(X + Y, n)` are tried.
//!
//! The sieve is not trial division and not Pollard's rho: the work is
//! dominated by the logarithmic sieve over the interval, as a real QS.
//!
//! # Practical range
//!
//! A 30-40 digit semiprime is the design point (well under ten seconds
//! on a current machine, in release mode and in the test profile, which
//! uses `opt-level = 1`). 45-50 digit semiprimes are attempted with the
//! larger default limits; they either factor within the work bound or
//! return [`None`] when the polynomial budget is exhausted. Inputs above
//! [`QsLimits::MAX_DIGITS`] digits return `None` immediately. The three
//! [`QsLimits`] fields bound the work: `bound` fixes the factor base,
//! `interval` the number of sieved `x` values per polynomial and
//! `polynomials` the number of polynomials tried.
//!
//! Every step is deterministic: the polynomial primes are taken in
//! increasing order, the sieve is sequential, the candidates are visited
//! in increasing `x` and the elimination is deterministic, so the result
//! is independent of the thread count. No random source is used.

// The sieve hot loops index fixed-size arrays in step with their formulas
// and do index arithmetic that is bounded by the array lengths; the code
// reads closest to the mathematical description this way.
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)]
// The helpers below are small, named pieces of the algorithm and are used
// from exactly one place; splitting them further would obscure the flow.
#![allow(clippy::single_call_fn)]
// Prime and exponent counts become float limits for the logarithmic sieve.
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

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
const MAX_DIGITS: u32 = 50;
/// Trial division by the primes below this before the sieve starts.
const SMALL_TRIAL: u64 = 1_000;
/// Relation-count target over the factor-base size (dependencies needed).
const EXTRA_DEPENDENCIES: usize = 8;

/// Limits for one quadratic-sieve attempt.
///
/// The defaults come from [`QsLimits::for_digits`] for 30-digit input.
/// Larger values enlarge the factor base (more relations needed) and the
/// sieve interval (more candidates per polynomial) but raise the cost of
/// one polynomial; `polynomials` is the hard work bound.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct QsLimits {
    /// Factor-base bound: primes `p <= bound` with `(n/p) = 1` are used.
    pub bound: u64,
    /// Half-width `M` of the sieving interval `[-M, M]`.
    pub interval: u64,
    /// Maximum number of polynomials tried before giving up.
    pub polynomials: usize,
}

impl QsLimits {
    /// Largest supported input size in decimal digits.
    ///
    /// Requests for larger numbers are refused so that no call can run
    /// for an unbounded time; the caller can still use ECM there.
    pub const MAX_DIGITS: u32 = MAX_DIGITS;

    /// Parameter set for an `digits`-digit number.
    ///
    /// The bound grows roughly like the square root of the
    /// `L_n[1/2, 1]` heuristic, and the interval with it, keeping the
    /// number of polynomials per attempt in the hundreds.
    #[must_use]
    pub const fn for_digits(digits: u32) -> Self {
        let (bound, interval, polynomials) = match digits {
            | 0..=12 => (120, 3_000, 400),
            | 13..=16 => (350, 5_000, 500),
            | 17..=20 => (800, 8_000, 600),
            | 21..=24 => (1_500, 12_000, 700),
            | 25..=28 => (3_000, 18_000, 800),
            | 29..=32 => (6_000, 25_000, 900),
            | 33..=36 => (12_000, 35_000, 1_000),
            | 37..=40 => (25_000, 50_000, 1_200),
            | 41..=44 => (45_000, 70_000, 1_400),
            | 45..=50 => (90_000, 90_000, 1_200),
            | _ => (120_000, 120_000, 800),
        };
        Self {
            bound,
            interval,
            polynomials,
        }
    }
}

impl Default for QsLimits {
    fn default() -> Self {
        Self::for_digits(30)
    }
}

/// Configuration of one sieve attempt.
///
/// `large_primes` turns the single-large-prime variation on or off; the
/// comparison is useful when measuring how much the variation helps.
/// `large_prime_factor` is the multiplier of the factor-base bound that
/// caps the large cofactor: a cofactor of at most
/// `large_prime_factor * bound` is kept, matched with another relation
/// that has the same cofactor, and the two are merged into a full
/// relation. The default of 16 keeps the number of possible cofactors
/// small enough that pairs are found quickly.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct QsConfig {
    /// Factor-base, interval and polynomial budget.
    pub limits: QsLimits,
    /// Whether partial relations with one large prime are collected.
    pub large_primes: bool,
    /// Multiplier of the factor-base bound for the large cofactor.
    pub large_prime_factor: u64,
}

impl Default for QsConfig {
    fn default() -> Self {
        Self {
            limits: QsLimits::default(),
            large_primes: true,
            large_prime_factor: 16,
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
    /// Dependencies found by the GF(2) elimination in the last solve.
    pub dependencies: usize,
}

/// Finds a non-trivial divisor of the odd composite `n`.
///
/// Returns `Some(d)` with `1 < d < n` and `d | n`, or `None` for invalid
/// input, for a prime, when `n` has more than [`QsLimits::MAX_DIGITS`]
/// decimal digits, and when the work bound [`QsLimits::polynomials`] is
/// exhausted. An even `n > 2` returns `2` directly.
///
/// The input may contain small factors: trial division by the primes
/// below 1000 runs first and returns a factor directly. The digit limit is
/// checked before that trial division, so an oversized request never does
/// any sieving work. Perfect powers, prime input and even input are
/// otherwise handled before the sieve.
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
    let bound = config.limits.bound.clamp(30, 500_000);
    let interval = usize::try_from(config.limits.interval.clamp(100, 2_000_000)).ok()?;
    if config.limits.polynomials == 0 {
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
    let divisors: Vec<BigInt> = factor_base.iter().map(|p| BigInt::from(p.prime)).collect();
    let large_bound = if config.large_primes {
        bound.saturating_mul(config.large_prime_factor.max(1))
    } else {
        1
    };
    let tolerance = (log2_approx(large_bound as f64) + log2_approx(bound as f64) + 2.0) * LOG_SCALE;
    let mut q_cursor = polynomial_prime_start(n, interval);
    let mut relations: Vec<Relation> = Vec::new();
    let mut partials: HashMap<u64, Relation> = HashMap::new();
    let mut target = columns + EXTRA_DEPENDENCIES;
    let max_relations = 4 * columns + 128;
    let length = 2 * interval + 1;
    for _ in 0..config.limits.polynomials {
        let q = match next_polynomial_prime(n, q_cursor) {
            | Ok(prime) => prime,
            | Err(factor) => return Some(factor),
        };
        q_cursor = q.saturating_add(2);
        let Some(polynomial) = build_polynomial(n, q) else {
            continue;
        };
        stats.polynomials += 1;
        let skip = factor_base.iter().position(|entry| entry.prime == q);
        let (roots1, roots2) = sieve_roots(&polynomial, &factor_base, skip, &divisors);
        let parity = usize::from(polynomial.c.is_odd());
        let mut sieve = vec![0_u32; length];
        let two_log = scaled_log(2);
        for index in (parity..length).step_by(2) {
            sieve[index] += two_log;
        }
        for (index, entry) in factor_base.iter().enumerate().skip(1) {
            if Some(index) == skip {
                continue;
            }
            let prime = usize::try_from(entry.prime).ok()?;
            let log = entry.log;
            let first = (usize::try_from(roots1[index]).ok()? + interval) % prime;
            let mut position = first;
            while position < length {
                sieve[position] += log;
                position += prime;
            }
            let second = roots2[index];
            if second != roots1[index] {
                let mut position = (usize::try_from(second).ok()? + interval) % prime;
                while position < length {
                    sieve[position] += log;
                    position += prime;
                }
            }
        }
        let mut candidate = vec![false; length];
        for (index, log) in sieve.iter().enumerate() {
            let x = index as f64 - interval as f64;
            let value = (polynomial.a_f64 * x + polynomial.two_b_f64) * x + polynomial.c_f64;
            if value != 0.0 {
                candidate[index] =
                    f64::from(*log) + tolerance >= log2_approx(value.abs()) * LOG_SCALE;
            }
        }
        let mut hits: Vec<(u32, u32)> = Vec::new();
        for index in (parity..length).step_by(2) {
            if candidate[index] {
                hits.push((index as u32, 0));
            }
        }
        for (index, entry) in factor_base.iter().enumerate().skip(1) {
            if Some(index) == skip {
                continue;
            }
            let prime = usize::try_from(entry.prime).ok()?;
            let mut positions = [roots1[index], roots2[index]];
            if positions[0] == positions[1] {
                positions[1] = u32::MAX;
            }
            for root in positions {
                if root == u32::MAX {
                    continue;
                }
                let mut position = (usize::try_from(root).ok()? + interval) % prime;
                while position < length {
                    if candidate[position] {
                        hits.push((position as u32, index as u32));
                    }
                    position += prime;
                }
            }
        }
        hits.sort_unstable_by_key(|&(position, prime)| (position, prime));
        let mut next_hit = 0_usize;
        for index in 0..length {
            if !candidate[index] {
                continue;
            }
            stats.candidates += 1;
            let x = index as i64 - interval as i64;
            let value = evaluate_g(&polynomial, x);
            let negative = value.is_negative();
            let mut magnitude = value.abs();
            let mut exponents: Vec<(u32, u32)> = Vec::new();
            while next_hit < hits.len() && usize::try_from(hits[next_hit].0).ok()? == index {
                let prime_index = usize::try_from(hits[next_hit].1).ok()?;
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
                known: &polynomial.q % n,
                negative,
                exponents,
            };
            if magnitude.is_one() {
                stats.full_relations += 1;
                relations.push(relation);
            } else if let Some(cofactor) = magnitude.to_u64().filter(|&value| value <= large_bound)
            {
                stats.partial_relations += 1;
                if config.large_primes {
                    match partials.remove(&cofactor) {
                        | Some(previous) => {
                            let combined = combine(&previous, &relation, cofactor, n);
                            stats.combined_relations += 1;
                            stats.full_relations += 1;
                            relations.push(combined);
                        },
                        | None => {
                            partials.insert(cofactor, relation);
                        },
                    }
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
            target = relations
                .len()
                .saturating_add((columns / 8).max(EXTRA_DEPENDENCIES));
            if target > max_relations {
                return None;
            }
        }
        if relations.len() >= max_relations {
            return None;
        }
    }
    if relations.len() > columns {
        solve_dependencies(&relations, &factor_base, n, stats)
    } else {
        None
    }
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
fn is_probable_prime(n: &BigInt) -> bool {
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
/// a prime exponent, and with at most 55 digits a larger exponent forces
/// a base below 8, which trial division has already removed.
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

/// Rows whose xor is zero, as bitsets over the row indices.
fn find_dependencies(
    rows: &[Vec<u64>],
    columns: usize,
) -> Vec<Vec<u64>> {
    let mut basis: Vec<Option<(Vec<u64>, Vec<u64>)>> = vec![None; columns];
    let combo_words = rows.len().div_ceil(64);
    let mut dependencies = Vec::new();
    for (row_index, row) in rows.iter().enumerate() {
        let mut reduced = row.clone();
        let mut combo = vec![0_u64; combo_words];
        combo[row_index / 64] |= 1_u64 << (row_index % 64);
        loop {
            match highest_bit(&reduced) {
                | Some(pivot) => match basis[pivot].take() {
                    | Some((row_bits, row_combo)) => {
                        for (target, source) in reduced.iter_mut().zip(&row_bits) {
                            *target ^= source;
                        }
                        for (target, source) in combo.iter_mut().zip(&row_combo) {
                            *target ^= source;
                        }
                        basis[pivot] = Some((row_bits, row_combo));
                    },
                    | None => {
                        basis[pivot] = Some((reduced, combo));
                        break;
                    },
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
