//! Integration tests for the quadratic sieve kernel and its use by
//! `factorint`.
//!
//! The tests are deterministic: primes are generated from fixed starting
//! points with a deterministic Miller-Rabin test, and every sieve call is
//! seeded by its own input.

use std::time::Instant;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use rssn::kernels::quadratic_sieve::QsConfig;
use rssn::kernels::quadratic_sieve::QsLimits;
use rssn::kernels::quadratic_sieve::QsStats;
use rssn::kernels::quadratic_sieve::qs_factor;
use rssn::kernels::quadratic_sieve::qs_factor_with;
use rssn::prelude::*;

/// `a * b mod m` without overflow.
fn mulmod(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(m)) as u64
}

/// `base^exponent mod m` for machine words.
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

/// Deterministic Miller-Rabin for a machine word.
fn is_prime_u64(n: u64) -> bool {
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

/// Deterministic Miller-Rabin for an arbitrary integer.
fn is_prime_big(n: &BigInt) -> bool {
    if let Some(value) = n.to_u64() {
        return is_prime_u64(value);
    }
    if !n.is_positive() || !n.bit(0) {
        return false;
    }
    let one = BigInt::one();
    let minus_one = n - &one;
    let shift = minus_one.trailing_zeros().unwrap_or(0);
    let odd = &minus_one >> shift;
    for base in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let base = BigInt::from(base);
        if base >= *n {
            continue;
        }
        let mut power = base.modpow(&odd, n);
        if power == one || power == minus_one {
            continue;
        }
        let mut composite = true;
        for _ in 1..shift {
            power = &power * &power % n;
            if power == minus_one {
                composite = false;
                break;
            }
        }
        if composite {
            return false;
        }
    }
    true
}

/// Smallest prime `>= n`.
fn next_prime_u64(n: u64) -> u64 {
    let mut candidate = n.max(2);
    if candidate > 2 {
        candidate |= 1;
    }
    loop {
        if is_prime_u64(candidate) {
            return candidate;
        }
        candidate = candidate.saturating_add(if candidate == 2 { 1 } else { 2 });
    }
}

/// Smallest prime `>= n` for arbitrary integers.
fn next_prime_big(n: &BigInt) -> BigInt {
    let mut candidate = if *n < BigInt::from(2) {
        BigInt::from(2)
    } else {
        n.clone()
    };
    if candidate > BigInt::from(2) && !candidate.bit(0) {
        candidate += 1;
    }
    loop {
        if is_prime_big(&candidate) {
            return candidate;
        }
        candidate += if candidate == BigInt::from(2) {
            1
        } else {
            2
        };
    }
}

/// A minimal deterministic xorshift generator for the statistical test.
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        let mut value = self.0;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        self.0 = value;
        value
    }
}

/// `10^exponent` as a machine word.
fn ten_pow(exponent: u32) -> u64 {
    (0..exponent).fold(1_u64, |accumulator, _| accumulator * 10)
}

/// Decimal digit count.
fn digits(n: &BigInt) -> u32 {
    u32::try_from(n.to_str_radix(10).len()).unwrap_or(u32::MAX)
}

/// Asserts that `divisor` is a proper factor of `n` and returns it.
fn check_divisor(
    n: &BigInt,
    divisor: Option<BigInt>,
) -> BigInt {
    let divisor = divisor.expect("a divisor");
    assert!(divisor > BigInt::one(), "divisor is not > 1");
    assert!(divisor < *n, "divisor is not < n");
    assert!((n % &divisor).is_zero(), "divisor does not divide n");
    divisor
}

#[test]
fn semiprimes_across_thirty_digits() {
    // Balanced and unbalanced products of known primes over 10^12..10^30.
    let cases: Vec<(u64, u64)> = vec![
        (next_prime_u64(1_000_000), next_prime_u64(1_000_500)),
        (next_prime_u64(1_000_003), next_prime_u64(1_000_000_007)),
        (next_prime_u64(10_000_003), next_prime_u64(10_000_019)),
        (next_prime_u64(1_000_000_007), next_prime_u64(1_000_000_009)),
        (next_prime_u64(16_777_219), next_prime_u64(1_000_000_000)),
        (next_prime_u64(1_000_000_000_039), next_prime_u64(1_000_003)),
        (
            next_prime_u64(1_000_000_000_000),
            next_prime_u64(1_000_000_000_061),
        ),
    ];
    for (p, q) in cases {
        let n = BigInt::from(p) * BigInt::from(q);
        let digits = digits(&n);
        assert!((12..=30).contains(&digits), "{n} has {digits} digits");
        let divisor = check_divisor(&n, qs_factor(&n, QsLimits::for_digits(digits)));
        assert!(
            divisor == BigInt::from(p) || divisor == BigInt::from(q),
            "{divisor} is not a factor of {n}"
        );
        assert!(is_prime_u64(p) && is_prime_u64(q));
    }
}

#[test]
fn squares_three_factors_small_factors_and_degenerate_input() {
    let p = next_prime_u64(1_000_000_000_000);
    let q = next_prime_u64(1_000_000_000_061);
    let r = next_prime_u64(1_000_003);

    // A prime square is handled as a perfect power.
    let square = BigInt::from(p) * BigInt::from(p);
    let divisor = check_divisor(&square, qs_factor(&square, QsLimits::for_digits(25)));
    assert_eq!(divisor, BigInt::from(p));

    // A cube.
    let cube = BigInt::from(r) * BigInt::from(r) * BigInt::from(r);
    let divisor = check_divisor(&cube, qs_factor(&cube, QsLimits::for_digits(18)));
    assert_eq!(divisor, BigInt::from(r));

    // Three prime factors.
    let triple = BigInt::from(p) * BigInt::from(q) * BigInt::from(r);
    let divisor = check_divisor(&triple, qs_factor(&triple, QsLimits::for_digits(30)));
    assert!(is_prime_u64(p) && is_prime_u64(q) && is_prime_u64(r));
    let _ = divisor;

    // A small factor is found before the sieve runs.
    let mixed = BigInt::from(3) * BigInt::from(p) * BigInt::from(q);
    assert_eq!(
        check_divisor(&mixed, qs_factor(&mixed, QsLimits::for_digits(28))),
        BigInt::from(3)
    );

    // Even, prime, negative, zero and one.
    let even = BigInt::from(2) * BigInt::from(p) * BigInt::from(q);
    assert_eq!(
        check_divisor(&even, qs_factor(&even, QsLimits::for_digits(30))),
        BigInt::from(2)
    );
    let prime = BigInt::from(next_prime_u64(1_000_000_000_000));
    assert!(qs_factor(&prime, QsLimits::for_digits(13)).is_none());
    assert!(qs_factor(&BigInt::zero(), QsLimits::default()).is_none());
    assert!(qs_factor(&BigInt::one(), QsLimits::default()).is_none());
    assert!(qs_factor(&BigInt::from(-12), QsLimits::default()).is_none());
    assert!(qs_factor(&BigInt::from(97), QsLimits::default()).is_none());
    // An even composite with a large prime factor is still fine.
    let twice = BigInt::from(2) * BigInt::from(next_prime_u64(1_000_000_003));
    assert_eq!(
        check_divisor(&twice, qs_factor(&twice, QsLimits::default())),
        BigInt::from(2)
    );
}

#[test]
fn one_forty_digit_semiprime() {
    // Two primes on either side of 10^19 and 10^20 make a 40-digit number.
    let p = next_prime_big(&BigInt::from(10_u64).pow(19));
    let q = next_prime_big(&BigInt::from(10_u64).pow(20));
    let n = &p * &q;
    assert_eq!(digits(&n), 40);
    let start = Instant::now();
    let mut stats = QsStats::default();
    let divisor = check_divisor(
        &n,
        qs_factor_with(
            &n,
            QsConfig {
                limits: QsLimits::for_digits(40),
                large_primes: true,
                ..QsConfig::default()
            },
            &mut stats,
        ),
    );
    let elapsed = start.elapsed();
    assert!(divisor == p || divisor == q);
    assert!(is_prime_big(&divisor));
    assert!(stats.dependencies > 0, "no GF(2) dependency was extracted");
    println!("40-digit semiprime factored in {elapsed:?}: {stats:?}");
    // Generous bound: the documented design point is well under ten
    // seconds; debug-instrumented runs on a slow machine may take longer.
    assert!(elapsed.as_secs() < 120, "QS took {elapsed:?}");
}

#[test]
fn statistical_fifty_semiprimes() {
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    let start = Instant::now();
    let mut factors = Vec::new();
    for _ in 0..50 {
        let total = 20 + random.next() % 11;
        let high_digits = u32::try_from(total / 2).unwrap_or(10);
        let low_digits = u32::try_from(total - u64::from(high_digits)).unwrap_or(10);
        // Both factors start with a 9, so the product keeps `total`
        // decimal digits instead of dropping to `total - 1`.
        let p =
            next_prime_u64(9 * ten_pow(high_digits - 1) + random.next() % ten_pow(high_digits - 1));
        let q =
            next_prime_u64(9 * ten_pow(low_digits - 1) + random.next() % ten_pow(low_digits - 1));
        let n = BigInt::from(p) * BigInt::from(q);
        let digits = digits(&n);
        assert!((20..=30).contains(&digits), "{n} has {digits} digits");
        let divisor = check_divisor(&n, qs_factor(&n, QsLimits::for_digits(digits)));
        assert!(
            divisor == BigInt::from(p) || divisor == BigInt::from(q),
            "{divisor} does not split {p} * {q}"
        );
        // The product of the two parts reconstructs the input.
        assert_eq!(&divisor * (&n / &divisor), n);
        factors.push((p, q));
    }
    let elapsed = start.elapsed();
    assert_eq!(factors.len(), 50);
    println!("50 random 20-30 digit semiprimes factored in {elapsed:?}");
    assert!(elapsed.as_secs() < 120, "statistical run took {elapsed:?}");
}

#[test]
fn large_prime_variation_adds_relations() {
    let p = next_prime_u64(100_000_000_000_003);
    let q = next_prime_u64(1_000_000_000_000_003);
    let n = BigInt::from(p) * BigInt::from(q);
    let digits = digits(&n);
    let limits = QsLimits::for_digits(digits);
    let mut with_stats = QsStats::default();
    let with = qs_factor_with(
        &n,
        QsConfig {
            limits,
            large_primes: true,
            ..QsConfig::default()
        },
        &mut with_stats,
    );
    let mut without_stats = QsStats::default();
    let without = qs_factor_with(
        &n,
        QsConfig {
            limits,
            large_primes: false,
            ..QsConfig::default()
        },
        &mut without_stats,
    );
    assert!(check_divisor(&n, with).is_positive());
    assert!(
        with_stats.partial_relations > 0,
        "no large-prime candidates"
    );
    assert!(
        with_stats.combined_relations > 0,
        "no partial relation paired"
    );
    assert!(
        with_stats.polynomials <= without_stats.polynomials,
        "large primes needed {} polynomials, without {}",
        with_stats.polynomials,
        without_stats.polynomials
    );
    // The variation can only help: the run without it either needs more
    // polynomials or does not finish inside the budget at all.
    println!(
        "with large primes: {:?} relations in {} polynomials; without: {:?} relations in {}",
        with_stats.full_relations,
        with_stats.polynomials,
        without_stats.full_relations,
        without_stats.polynomials
    );
    if let Some(divisor) = without {
        assert!(divisor > BigInt::one() && divisor < n && (&n % &divisor).is_zero());
    }
}

#[test]
fn deterministic_across_runs_and_thread_counts() {
    let p = next_prime_u64(10_000_000_000_000_000);
    let q = next_prime_u64(1_000_000_000_000_000_003);
    let n = BigInt::from(p) * BigInt::from(q);
    let limits = QsLimits::for_digits(digits(&n));
    let first = check_divisor(&n, qs_factor(&n, limits));
    for _ in 0..2 {
        assert_eq!(qs_factor(&n, limits), Some(first.clone()));
    }
    let pool_one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("thread pool");
    let one = pool_one.install(|| qs_factor(&n, limits));
    let pool_four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("thread pool");
    let four = pool_four.install(|| qs_factor(&n, limits));
    assert_eq!(one, Some(first.clone()));
    assert_eq!(four, Some(first));
}

/// The `factorint` text for a request.
fn factorint(
    session: &Session,
    request: &str,
) -> String {
    let term = session
        .parse(&format!("factorint({request})"))
        .expect("parse");
    session
        .compute(term, &Config::new())
        .expect("compute")
        .term
        .to_string()
}

#[test]
fn factorint_integration() {
    let session = Session::new();
    // The Mersenne and Fermat cases that the existing number-theory tests
    // cover must keep working through the same schedule.
    assert_eq!(
        factorint(&session, "2^101 - 1"),
        "list(list(7432339208719, 1), list(341117531003194129, 1))"
    );
    assert_eq!(
        factorint(&session, "2^128 + 1"),
        "list(list(59649589127497217, 1), list(5704689200685129054721, 1))"
    );

    // A 40-digit semiprime goes through the quadratic sieve.
    let p = next_prime_big(&BigInt::from(10_u64).pow(19));
    let q = next_prime_big(&BigInt::from(10_u64).pow(20));
    let n = &p * &q;
    let expected = if p < q {
        format!("list(list({p}, 1), list({q}, 1))")
    } else {
        format!("list(list({q}, 1), list({p}, 1))")
    };
    let start = Instant::now();
    assert_eq!(factorint(&session, &n.to_string()), expected);
    println!(
        "factorint of a 40-digit semiprime took {:?}",
        start.elapsed()
    );
}
