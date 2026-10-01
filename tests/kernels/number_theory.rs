//! Integer number theory (ported from `numerical_number_theory_test.rs`).

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::number_theory::{
    factorize, gcd, is_prime_miller_rabin, lcm, mod_inverse, mod_pow, phi, primes_sieve,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn gcd_lcm_modpow_modinverse() {
    assert_eq!(gcd(48, 18), 6);
    assert_eq!(lcm(48, 18), 144);
    assert_eq!(mod_pow(2, 10, 1000), 24); // 1024 % 1000
    assert_eq!(mod_inverse(3, 11), Some(4)); // 3 * 4 = 12 = 1 (mod 11)
}

#[test]
fn gcd_lcm_edge_cases() {
    assert_eq!(gcd(0, 7), 7);
    assert_eq!(gcd(7, 0), 7);
    assert_eq!(gcd(17, 13), 1);
    assert_eq!(lcm(0, 5), 0);
    assert_eq!(lcm(4, 6), 12);
}

#[test]
fn modular_inverse_fails_when_not_coprime() {
    assert_eq!(mod_inverse(4, 8), None);
    assert_eq!(mod_inverse(6, 9), None);
    assert_eq!(mod_inverse(7, 26), Some(15)); // 7 * 15 = 105 = 4 * 26 + 1
}

#[test]
fn mod_pow_reference_values() {
    assert_eq!(mod_pow(3, 0, 7), 1);
    assert_eq!(mod_pow(5, 117, 19), 1); // ord(5) = 9 mod 19 and 9 | 117
    let direct = (0..117).fold(1u128, |acc, _| acc * 5 % 19) as u64;
    assert_eq!(mod_pow(5, 117, 19), direct);
    // Fermat's little theorem with a 61-bit Mersenne prime.
    let p = (1u64 << 61) - 1;
    assert_eq!(mod_pow(123_456_789, p - 1, p), 1);
}

#[test]
fn mod_pow_modulus_one_is_zero() {
    assert_eq!(mod_pow(5, 0, 1), 0);
}

#[test]
fn primality_small() {
    for p in [2u64, 3, 5, 7, 17, 104_729] {
        assert!(is_prime_miller_rabin(p), "{p} is prime");
    }
    for c in [0u64, 1, 4, 9, 100, 104_730] {
        assert!(!is_prime_miller_rabin(c), "{c} is composite or unit");
    }
}

#[test]
fn primality_against_sieve() {
    let sieve = primes_sieve(5000);
    let set: std::collections::HashSet<u64> = sieve.iter().map(|&p| p as u64).collect();
    for n in 0..=5000u64 {
        assert_eq!(is_prime_miller_rabin(n), set.contains(&n), "n = {n}");
    }
}

#[test]
fn primality_hard_cases() {
    // Carmichael numbers and strong pseudoprimes to many bases.
    for c in [
        561u64,
        1105,
        1729,
        2465,
        2821,
        6601,
        3_215_031_751,
        3_825_123_056_546_413_051,
    ] {
        assert!(!is_prime_miller_rabin(c), "{c} is composite");
    }
    // Mersenne prime 2^61 - 1 and the largest 64-bit prime.
    assert!(is_prime_miller_rabin((1u64 << 61) - 1));
    assert!(is_prime_miller_rabin(18_446_744_073_709_551_557));
    assert!(!is_prime_miller_rabin(18_446_744_073_709_551_557 - 2));
}

#[test]
fn totient() {
    assert_eq!(phi(0), 0);
    assert_eq!(phi(1), 1);
    assert_eq!(phi(10), 4); // 1, 3, 7, 9
    assert_eq!(phi(11), 10);
    assert_eq!(phi(36), 12);
    assert_eq!(phi(1024), 512);
}

#[test]
fn factorization() {
    assert_eq!(factorize(12), vec![2, 2, 3]);
    assert_eq!(factorize(60), vec![2, 2, 3, 5]);
    assert_eq!(factorize(17), vec![17]);
    assert_eq!(factorize(1), Vec::<u64>::new());
    assert_eq!(factorize(0), Vec::<u64>::new());
    assert_eq!(factorize(1 << 20), vec![2; 20]);
    assert_eq!(factorize(600_851_475_143), vec![71, 839, 1471, 6857]);
}

#[test]
fn sieve() {
    assert_eq!(primes_sieve(20), vec![2, 3, 5, 7, 11, 13, 17, 19]);
    assert_eq!(primes_sieve(2), vec![2]);
    assert!(primes_sieve(1).is_empty());
    assert!(primes_sieve(0).is_empty());
    assert_eq!(primes_sieve(1000).len(), 168);
    assert_eq!(primes_sieve(100_000).len(), 9592);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_gcd_times_lcm_is_product(a in 1..10_000u64, b in 1..10_000u64) {
        prop_assert_eq!(a * b, gcd(a, b) * lcm(a, b));
    }

    #[test]
    fn prop_gcd_divides_both(a in 1..1_000_000u64, b in 1..1_000_000u64) {
        let g = gcd(a, b);
        prop_assert!(a % g == 0 && b % g == 0);
        prop_assert_eq!(gcd(a / g, b / g), 1);
    }

    #[test]
    fn prop_mod_inverse_is_an_inverse_when_it_exists(a in 1..1000i64, m in 2..1000i64) {
        match mod_inverse(a, m) {
            Some(inv) => prop_assert_eq!(a % m * inv % m, 1 % m),
            None => prop_assert!(gcd(a as u64, m as u64) != 1, "coprime {a}, {m} must have an inverse"),
        }
    }

    #[test]
    fn prop_mod_pow_matches_naive(base in 0..50u128, exp in 0..40u64, modulus in 2..10_000u64) {
        let naive = (0..exp).fold(1u128, |acc, _| acc * base % u128::from(modulus)) as u64;
        prop_assert_eq!(mod_pow(base, exp, modulus), naive);
    }

    #[test]
    fn prop_factorize_multiplies_back_to_n_with_prime_factors(n in 2..1_000_000u64) {
        let f = factorize(n);
        prop_assert_eq!(f.iter().product::<u64>(), n);
        prop_assert!(f.windows(2).all(|w| w[0] <= w[1]));
        prop_assert!(f.iter().all(|&p| is_prime_miller_rabin(p)));
    }

    #[test]
    fn prop_phi_is_multiplicative_on_coprimes(a in 1..300u64, b in 1..300u64) {
        prop_assume!(gcd(a, b) == 1);
        prop_assert_eq!(phi(a * b), phi(a) * phi(b));
    }

    #[test]
    fn prop_euler_theorem(a in 2..500u64, n in 3..500u64) {
        prop_assume!(gcd(a, n) == 1);
        prop_assert_eq!(mod_pow(u128::from(a), phi(n), n), 1);
    }
}

#[test]
fn mod_pow_zero_exponent_and_small_cases() {
    // x^0 = 1 for any modulus > 1, including x = 0 by convention.
    assert_eq!(mod_pow(0, 0, 7), 1);
    assert_eq!(mod_pow(9, 0, 2), 1);
    // Everything is 0 modulo 1.
    assert_eq!(mod_pow(123, 45, 1), 0);
    // 2^10 = 1024 = 4 (mod 10).
    assert_eq!(mod_pow(2, 10, 10), 4);
}
