//! Exact integer and signed-word number theory for language frontends.
use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

/// Prime factorization of the magnitude of `n`, grouped in increasing order.
#[must_use]
pub fn factor(n: i64) -> Vec<(i64, u32)> {
    let factors = super::number_theory::factorize(n.unsigned_abs());
    let mut result: Vec<(i64, u32)> = Vec::new();
    for prime in factors {
        // The only magnitude above i64::MAX is abs(i64::MIN), whose factors
        // are all two. Every prime factor is therefore representable.
        let prime = prime as i64;
        match result.last_mut() {
            | Some((last, count)) if *last == prime => *count += 1,
            | _ => result.push((prime, 1)),
        }
    }
    result
}
/// Deterministic primality for every unsigned machine word.
#[must_use]
pub const fn is_prime(n: u64) -> bool {
    super::number_theory::is_prime_miller_rabin(n)
}
/// Prime numbers strictly below `limit`; negative limits yield an empty list.
#[must_use]
pub fn primes(limit: i64) -> Vec<i64> {
    if limit <= 2 {
        return Vec::new();
    }
    let Ok(bound) = usize::try_from(limit - 1) else {
        return Vec::new();
    };
    super::number_theory::primes_sieve(bound)
        .into_iter()
        .map(|n| n as i64)
        .collect()
}
/// Exact binomial coefficient; impossible or negative selections yield zero.
#[must_use]
pub fn binomial(
    n: i64,
    k: i64,
) -> BigInt {
    if k < 0 || n < 0 || k > n {
        return BigInt::zero();
    }
    let (n, k) = (n as u64, k as u64);
    let k = k.min(n - k);
    if k > 64 {
        return word_product(n - k + 1, n) / word_product(2, k);
    }
    let mut result = BigInt::one();
    for step in 0..k {
        result *= n - step;
        result /= step + 1;
    }
    result
}
/// Exact factorial, conventionally zero for negative inputs.
#[must_use]
pub fn permutations(n: i64) -> BigInt {
    if n < 0 {
        return BigInt::zero();
    }
    word_product(2, n as u64)
}

/// Exact number of ordered selections of `k` distinct members from `n`.
/// An empty selection has one ordering. Negative or impossible inputs return zero.
#[must_use]
pub fn falling_factorial(
    n: i64,
    k: i64,
) -> BigInt {
    if k < 0 || n < 0 || k > n {
        return BigInt::zero();
    }
    word_product((n - k) as u64 + 1, n as u64)
}

fn word_product(
    first: u64,
    last: u64,
) -> BigInt {
    if last.saturating_sub(first) < 32 {
        let mut result = BigInt::one();
        let mut block = 1_u128;
        for word in first..=last {
            let word = u128::from(word);
            if let Some(product) = block.checked_mul(word) {
                block = product;
            } else {
                result *= block;
                block = word;
            }
        }
        return result * block;
    }
    let middle = first.midpoint(last);
    word_product(first, middle) * word_product(middle + 1, last)
}
/// Greatest common divisor of signed magnitudes.
///
/// Preserves the signed-word ABI: the sole unrepresentable result, 2^63,
/// uses its `i64::MIN` bit representation. Use arbitrary precision for that case.
#[must_use]
pub const fn gcd(
    a: i64,
    b: i64,
) -> i64 {
    super::number_theory::gcd(a.unsigned_abs(), b.unsigned_abs()) as i64
}
/// Nonnegative least common multiple, saturated at `i64::MAX` on overflow.
#[must_use]
pub fn lcm(
    a: i64,
    b: i64,
) -> i64 {
    if a == 0 || b == 0 {
        return 0;
    }
    let divisor = super::number_theory::gcd(a.unsigned_abs(), b.unsigned_abs());
    let magnitude = u128::from(a.unsigned_abs() / divisor) * u128::from(b.unsigned_abs());
    i64::try_from(magnitude).unwrap_or(i64::MAX)
}
/// Whether an exact integer is negative.
#[must_use]
pub fn is_negative(value: &BigInt) -> bool {
    value.is_negative()
}
