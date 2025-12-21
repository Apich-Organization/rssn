//! # Numerical Number Theory
//!
//! This module provides numerical implementations for number theory algorithms.
//! It includes functions for computing the greatest common divisor (GCD),
//! modular exponentiation, modular inverse, and primality testing using the Miller-Rabin algorithm.

/// Computes the greatest common divisor (GCD) of two numbers using the Euclidean algorithm.
///
/// # Arguments
/// * `a` - The first number.
/// * `b` - The second number.
///
/// # Returns
/// The greatest common divisor of `a` and `b`.
pub fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Computes `(base^exp) % modulus` efficiently using modular exponentiation (binary exponentiation).
///
/// # Arguments
/// * `base` - The base.
/// * `exp` - The exponent.
/// * `modulus` - The modulus.
///
/// # Returns
/// The result of `(base^exp) % modulus`.
pub fn mod_pow(mut base: u128, mut exp: u64, modulus: u64) -> u64 {
    let mut res = 1;
    base %= u128::from(modulus);
    while exp > 0 {
        if exp % 2 == 1 {
            res = (res * base) % u128::from(modulus);
        }
        base = (base * base) % u128::from(modulus);
        exp /= 2;
    }
    res as u64
}

/// Finds the modular multiplicative inverse of a number.
///
/// Solves for `x` in `ax ≡ 1 (mod m)` using the Extended Euclidean Algorithm.
///
/// # Arguments
/// * `a` - The number for which to find the inverse.
/// * `m` - The modulus.
///
/// # Returns
/// An `Option<i64>` containing the modular inverse if it exists, otherwise `None`.
pub fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(a, m);
    if g == 1 {
        Some((x % m + m) % m)
    } else {
        None
    }
}

/// Extended Euclidean algorithm for i64.
pub(crate) fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

/// Performs a Miller-Rabin primality test.
///
/// The Miller-Rabin test is a probabilistic primality test. This implementation
/// uses a set of bases sufficient to deterministically test all `u64` numbers.
///
/// # Arguments
/// * `n` - The number to test for primality.
///
/// # Returns
/// `true` if `n` is prime, `false` otherwise.
pub fn is_prime_miller_rabin(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &a in &bases {
        if n == a {
            return true;
        }
        let mut x = mod_pow(u128::from(a), d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut composite = true;
        for _ in 0..s - 1 {
            x = mod_pow(u128::from(x), 2, n);
            if x == n - 1 {
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

/// Computes the least common multiple (LCM) of two numbers.
pub fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 {
        return 0;
    }
    (a / gcd(a, b)) * b
}

/// Computes Euler's totient function φ(n).
pub fn phi(mut n: u64) -> u64 {
    if n == 0 { return 0; }
    let mut result = n;
    let mut i = 2;
    while i * i <= n {
        if n % i == 0 {
            while n % i == 0 {
                n /= i;
            }
            result -= result / i;
        }
        i += 1;
    }
    if n > 1 {
        result -= result / n;
    }
    result
}

/// Returns a list of prime factors of n.
pub fn factorize(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    if n < 2 { return factors; }
    let mut i = 2;
    while i * i <= n {
        while n % i == 0 {
            factors.push(i);
            n /= i;
        }
        i += 1;
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

/// Generates primes up to limit using Sieve of Eratosthenes.
pub fn primes_sieve(limit: usize) -> Vec<usize> {
    if limit < 2 { return Vec::new(); }
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;
    for p in 2..=(limit as f64).sqrt() as usize {
        if is_prime[p] {
            let mut i = p * p;
            while i <= limit {
                is_prime[i] = false;
                i += p;
            }
        }
    }
    is_prime.into_iter().enumerate()
        .filter(|&(_, p)| p)
        .map(|(i, _)| i)
        .collect()
}
