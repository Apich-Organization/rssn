//! Edge cases and bounds of the machine-word number theory kernels:
//! no panics on any input, deterministic primality over all of `u64`,
//! fast factorisation of semiprimes near `2^64`, modular square roots and
//! discrete logarithms.

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::number_theory as nt;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0xA11CE),
        failure_persistence: None,
        cases: 256,
        ..ProptestConfig::default()
    }
}

#[test]
fn mod_pow_degenerate_moduli() {
    assert_eq!(
        nt::mod_pow(2, 10, 0),
        0,
        "modulus 0 returns 0, no division by zero"
    );
    assert_eq!(nt::mod_pow(7, 0, 1), 0);
    assert_eq!(nt::mod_pow(0, 0, 5), 1);
    // The full u64 range: u128 intermediates must not overflow.
    let m = u64::MAX - 58;
    let r = nt::mod_pow(u128::MAX, u64::MAX, m);
    assert!(r < m);
    assert_eq!(
        nt::mod_pow(2, m - 1, m),
        1,
        "Fermat: m = 2^64 - 59 is prime"
    );
}

#[test]
fn mod_inverse_range_and_signs() {
    assert_eq!(nt::mod_inverse(3, 0), None);
    assert_eq!(nt::mod_inverse(3, 7), Some(5));
    assert_eq!(nt::mod_inverse(-3, 7), Some(2));
    assert_eq!(nt::mod_inverse(3, -7), Some(5), "negative modulus uses |m|");
    assert_eq!(nt::mod_inverse(5, 1), Some(0));
    assert_eq!(nt::mod_inverse(0, 7), None);
    assert_eq!(nt::mod_inverse(4, 8), None);
    // Extreme operands used to overflow `i64` arithmetic.
    let m = i64::MAX;
    let inv = nt::mod_inverse(3, m).unwrap();
    assert_eq!(i128::from(inv) * 3 % i128::from(m), 1);
    assert_eq!(
        nt::mod_inverse(i64::MIN, 7).map(|x| (x as i128 * i64::MIN as i128).rem_euclid(7)),
        Some(1)
    );
    assert_eq!(
        nt::mod_inverse(1, i64::MIN),
        Some(1),
        "|i64::MIN| = 2^63 is handled in i128"
    );
    assert_eq!(nt::mod_inverse(2, i64::MIN), None);
}

#[test]
fn gcd_lcm_extremes() {
    assert_eq!(nt::gcd(u64::MAX, u64::MAX - 1), 1);
    assert_eq!(nt::gcd(0, 0), 0);
    assert_eq!(
        nt::lcm(u64::MAX, u64::MAX - 1),
        u64::MAX,
        "saturates instead of overflowing"
    );
    assert_eq!(nt::checked_lcm(u64::MAX, u64::MAX - 1), None);
    assert_eq!(nt::checked_lcm(0, u64::MAX), Some(0));
    assert_eq!(nt::checked_lcm(4, 6), Some(12));
}

#[test]
fn miller_rabin_is_deterministic_for_all_u64() {
    // Largest primes below 2^64 and strong pseudoprimes to many bases.
    for p in [
        18_446_744_073_709_551_557_u64,
        18_446_744_073_709_551_533,
        4_294_967_291,
        2,
        3,
        97,
    ] {
        assert!(nt::is_prime_miller_rabin(p), "{p}");
    }
    for c in [
        0_u64,
        1,
        u64::MAX,
        // Strong pseudoprime to the twelve bases 2..=37 below 3.3e24 exists,
        // but none below 2^64; this one fools bases 2..=31 only.
        3_825_123_056_546_413_051,
        // Carmichael numbers and strong pseudoprimes to base 2.
        561,
        41_041,
        825_265,
        2_047,
        3_215_031_751,
        // Product of two 32-bit primes.
        4_294_967_291 * 4_294_967_279,
        // 2^64 - 1 is divisible by 3, so is 2^64 - 3 * 7; evens.
        18_446_744_073_709_551_614,
    ] {
        assert!(!nt::is_prime_miller_rabin(c), "{c} is composite");
    }
    // Cross-check the whole range below 100000 with a sieve.
    let sieve = nt::primes_sieve(100_000);
    let set: std::collections::HashSet<usize> = sieve.iter().copied().collect();
    for n in 0..=100_000_usize {
        assert_eq!(nt::is_prime_miller_rabin(n as u64), set.contains(&n), "{n}");
    }
}

#[test]
fn factorize_large_semiprimes_quickly() {
    let t = std::time::Instant::now();
    // Largest 32-bit primes: the product is a 64-bit semiprime that trial
    // division could not touch in reasonable time.
    let (p, q) = (4_294_967_291_u64, 4_294_967_279_u64);
    assert_eq!(nt::factorize(p * q), vec![q, p]);
    assert_eq!(
        nt::factorize(18_446_744_073_709_551_557),
        vec![18_446_744_073_709_551_557]
    );
    assert_eq!(
        nt::factorize(u64::MAX),
        vec![3, 5, 17, 257, 641, 65_537, 6_700_417]
    );
    assert_eq!(nt::factorize(1 << 63), vec![2; 63]);
    // A prime square exercises the perfect-square path of rho.
    assert_eq!(
        nt::factorize(4_294_967_291 * 4_294_967_291),
        vec![4_294_967_291; 2]
    );
    assert_eq!(nt::factorize(0), Vec::<u64>::new());
    assert_eq!(nt::factorize(1), Vec::<u64>::new());
    assert!(t.elapsed().as_secs_f64() < 10.0, "took {:?}", t.elapsed());
}

#[test]
fn phi_no_overflow_and_known_values() {
    assert_eq!(nt::phi(0), 0);
    assert_eq!(nt::phi(1), 1);
    assert_eq!(nt::phi(36), 12);
    assert_eq!(nt::phi(97), 96);
    // u64::MAX = 3 * 5 * 17 * 257 * 641 * 65537 * 6700417
    assert_eq!(
        nt::phi(u64::MAX),
        2 * 4 * 16 * 256 * 640 * 65_536 * 6_700_416
    );
    let p = 4_294_967_291_u64;
    assert_eq!(nt::phi(p * p), p * (p - 1));
}

#[test]
fn primes_sieve_degenerate() {
    assert!(nt::primes_sieve(0).is_empty());
    assert!(nt::primes_sieve(1).is_empty());
    assert_eq!(nt::primes_sieve(2), vec![2]);
    assert_eq!(
        nt::primes_sieve(30),
        vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    );
    assert!(
        nt::primes_sieve(usize::MAX).is_empty(),
        "overflow of limit + 1 is handled"
    );
    assert_eq!(nt::primes_sieve(1_000_000).len(), 78_498);
}

#[test]
fn tonelli_shanks_roots() {
    // p = 1 mod 8 forces the full Tonelli-Shanks loop.
    for &p in &[
        17_u64,
        41,
        73,
        97,
        998_244_353,
        1_000_000_007,
        4_294_967_291,
        18_446_744_073_709_551_557,
    ] {
        for n in [0_u64, 1, 2, 3, 4, 5, 10, 123_456_789] {
            match nt::sqrt_mod_prime(n, p) {
                | Some(r) => {
                    let sq = (u128::from(r) * u128::from(r)) % u128::from(p);
                    assert_eq!(sq, u128::from(n % p), "sqrt({n}) mod {p}");
                    assert!(r <= p - r || r == 0, "smaller root returned");
                },
                | None => {
                    // Euler's criterion says non-residue.
                    assert_eq!(
                        nt::mod_pow(u128::from(n), (p - 1) / 2, p),
                        p - 1,
                        "{n} mod {p}"
                    );
                },
            }
        }
    }
    assert_eq!(nt::sqrt_mod_prime(4, 2), Some(0));
    assert_eq!(nt::sqrt_mod_prime(3, 2), Some(1));
    assert_eq!(nt::sqrt_mod_prime(2, 15), None, "composite modulus");
    assert_eq!(nt::sqrt_mod_prime(2, 0), None);
    assert_eq!(nt::sqrt_mod_prime(2, 1), None);
}

#[test]
fn discrete_log_baby_giant() {
    assert_eq!(nt::discrete_log(2, 1, 101), Some(0));
    let (g, p) = (3_u64, 1_000_000_007_u64);
    for x in [0_u64, 1, 2, 12345, 999_999, 500_000_003] {
        let h = nt::mod_pow(u128::from(g), x, p);
        let got = nt::discrete_log(g, h, p).unwrap();
        assert_eq!(nt::mod_pow(u128::from(g), got, p), h);
        assert!(got <= x, "least exponent");
    }
    // 4 is not in <2> mod 7? 2^2 = 4: yes; 3 is not a power of 2 mod 7.
    assert_eq!(nt::discrete_log(2, 4, 7), Some(2));
    assert_eq!(nt::discrete_log(2, 3, 7), None);
    assert_eq!(nt::discrete_log(2, 3, 0), None);
    assert_eq!(nt::discrete_log(2, 3, 1), None);
    assert_eq!(
        nt::discrete_log(2, 3, u64::MAX),
        None,
        "table would be too large"
    );
    // Non-invertible base.
    assert_eq!(nt::discrete_log(0, 0, 7), Some(1));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn nothing_panics(a in any::<u64>(), b in any::<u64>(), c in any::<u64>(), i in any::<i64>(), j in any::<i64>()) {
        let _ = nt::gcd(a, b);
        let _ = nt::lcm(a, b);
        let _ = nt::mod_pow(u128::from(a) * u128::from(b), c, a);
        let _ = nt::mod_inverse(i, j);
        let _ = nt::is_prime_miller_rabin(a);
        let _ = nt::sqrt_mod_prime(a, b);
        let _ = nt::discrete_log(a, b, c % (1 << 20));
    }

    #[test]
    fn factorize_multiplies_back_to_n(n in 2_u64..u64::MAX) {
        let f = nt::factorize(n);
        prop_assert!(f.windows(2).all(|w| w[0] <= w[1]), "sorted");
        let product = f.iter().fold(1_u128, |acc, &p| acc * u128::from(p));
        prop_assert_eq!(product, u128::from(n));
        for p in f {
            prop_assert!(nt::is_prime_miller_rabin(p), "{} is not prime", p);
        }
    }

    #[test]
    fn mod_inverse_is_an_inverse(a in any::<i64>(), m in 2_i64..i64::MAX) {
        match nt::mod_inverse(a, m) {
            Some(x) => {
                prop_assert!((0..m).contains(&x));
                prop_assert_eq!((i128::from(a) * i128::from(x)).rem_euclid(i128::from(m)), 1);
            }
            None => prop_assert!(nt::gcd(a.unsigned_abs(), m as u64) != 1),
        }
    }

    #[test]
    fn sqrt_mod_prime_squares_back(n in any::<u64>(), pi in 0_usize..2000) {
        let primes = nt::primes_sieve(20_000);
        let p = primes[pi % primes.len()] as u64;
        if let Some(r) = nt::sqrt_mod_prime(n, p) {
            prop_assert_eq!(u128::from(r) * u128::from(r) % u128::from(p), u128::from(n % p));
        }
    }
}
