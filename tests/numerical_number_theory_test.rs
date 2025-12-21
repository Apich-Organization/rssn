use rssn::prelude::*;
use rssn::prelude::numerical::*;
use proptest::prelude::*;

#[test]
fn test_nt_basic() {
    assert_eq!(numerical_gcd(48, 18), 6);
    assert_eq!(numerical_lcm(48, 18), 144);
    assert_eq!(numerical_mod_pow(2, 10, 1000), 24); // 1024 % 1000
    assert_eq!(numerical_mod_inverse(3, 11), Some(4)); // 3*4 = 12 ≡ 1 mod 11
}

#[test]
fn test_nt_primality() {
    assert!(numerical_is_prime(2));
    assert!(numerical_is_prime(3));
    assert!(numerical_is_prime(17));
    assert!(numerical_is_prime(104729)); // 10000th prime
    assert!(!numerical_is_prime(1));
    assert!(!numerical_is_prime(4));
    assert!(!numerical_is_prime(100));
}

#[test]
fn test_nt_phi() {
    assert_eq!(numerical_phi(1), 1);
    assert_eq!(numerical_phi(10), 4); // 1, 3, 7, 9
    assert_eq!(numerical_phi(11), 10);
}

#[test]
fn test_nt_factorize() {
    assert_eq!(numerical_factorize(12), vec![2, 2, 3]);
    assert_eq!(numerical_factorize(60), vec![2, 2, 3, 5]);
    assert_eq!(numerical_factorize(17), vec![17]);
}

#[test]
fn test_nt_sieve() {
    let primes = numerical_primes_sieve(20);
    assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
}

proptest! {
    #[test]
    fn prop_gcd_lcm(a in 1..10000u64, b in 1..10000u64) {
        let g = numerical_gcd(a, b);
        let l = numerical_lcm(a, b);
        assert_eq!(a * b, g * l);
    }

    #[test]
    fn prop_mod_inverse(a in 1..1000i64, m in 2..1000i64) {
        if let Some(inv) = numerical_mod_inverse(a, m) {
            let res = (a % m * inv) % m;
            let expected = 1 % m;
            assert_eq!(res, expected);
        }
    }
}
