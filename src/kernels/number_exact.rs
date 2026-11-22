//! Exact integer and signed-word number theory for language frontends.
#![allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn
)]

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::number_theory::MontU64;
use super::number_theory::gcd_u128;
use super::number_theory::is_prime_miller_rabin;
use super::number_theory::isqrt;
use super::number_theory::mod_inverse;
use super::number_theory::mul_mod;
use super::number_theory::primes_sieve;
use crate::kernels::reuse::global_pool;

const PRIME_POWER_BINOMIAL_MIN_K: u64 = 256;
const PRIME_POWER_BINOMIAL_MAX_N: u64 = 2_000_000;

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
    is_prime_miller_rabin(n)
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
    primes_sieve(bound).into_iter().map(|n| n as i64).collect()
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
    if k == 0 {
        return BigInt::one();
    }
    if k == 1 {
        return BigInt::from(n);
    }
    if k > PRIME_POWER_BINOMIAL_MIN_K && n <= PRIME_POWER_BINOMIAL_MAX_N {
        return binomial_prime_power(n as usize, k as usize);
    }
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

fn word_product_leaf(
    first: u64,
    last: u64,
) -> BigUint {
    let mut result = BigUint::one();
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
    result * block
}

fn word_product_uint(
    first: u64,
    last: u64,
) -> BigUint {
    if first > last {
        return BigUint::one();
    }
    if first == 0 {
        return BigUint::zero();
    }
    let span = last - first;
    if span == 0 {
        return BigUint::from(first);
    }
    if span == 1 {
        return BigUint::from(u128::from(first) * u128::from(last));
    }
    if span < 32 {
        return word_product_leaf(first, last);
    }
    let middle = first.midpoint(last);
    let work_items = usize::try_from(span).unwrap_or(usize::MAX);
    if crate::execution::should_split_work(work_items, 8) {
        let token = crate::cancellation::current();
        let left_token = token.clone();
        let (left, right) = rayon::join(
            move || {
                let _guard = crate::cancellation::install_optional(left_token);
                word_product_uint(first, middle)
            },
            move || {
                let _guard = crate::cancellation::install_optional(token);
                word_product_uint(middle + 1, last)
            },
        );
        left * right
    } else {
        word_product_uint(first, middle) * word_product_uint(middle + 1, last)
    }
}

fn word_product(
    first: u64,
    last: u64,
) -> BigInt {
    BigInt::from(word_product_uint(first, last))
}

fn product_tree_u128(leaves: &[u128]) -> BigUint {
    if leaves.is_empty() {
        return BigUint::one();
    }
    if leaves.len() <= 16 {
        let mut acc = BigUint::from(leaves[0]);
        for &leaf in &leaves[1..] {
            acc *= leaf;
        }
        return acc;
    }
    let mid = leaves.len() / 2;
    let (left_slice, right_slice) = leaves.split_at(mid);
    if crate::execution::should_split_work(leaves.len(), 64) {
        let token = crate::cancellation::current();
        let left_token = token.clone();
        let (left, right) = rayon::join(
            move || {
                let _guard = crate::cancellation::install_optional(left_token);
                product_tree_u128(left_slice)
            },
            move || {
                let _guard = crate::cancellation::install_optional(token);
                product_tree_u128(right_slice)
            },
        );
        left * right
    } else {
        product_tree_u128(left_slice) * product_tree_u128(right_slice)
    }
}

fn binomial_prime_power(
    n: usize,
    k: usize,
) -> BigInt {
    let primes = primes_sieve(n);
    let m = n - k;
    let root = isqrt(n as u64) as usize;
    let half_n = n / 2;
    let mut leaves = Vec::with_capacity(primes.len() / 4 + 8);
    let mut block = 1_u128;

    let mut push_factor = |factor: u64| {
        let factor = u128::from(factor);
        if let Some(prod) = block.checked_mul(factor) {
            block = prod;
        } else {
            leaves.push(block);
            block = factor;
        }
    };

    for &p in &primes {
        if p <= root {
            let mut exp = 0_u32;
            let mut nn = n;
            let mut kk = k;
            let mut mm = m;
            while nn >= p {
                nn /= p;
                kk /= p;
                mm /= p;
                exp += (nn - kk - mm) as u32;
            }
            if exp > 0 {
                push_factor((p as u64).pow(exp));
            }
        } else if p <= half_n {
            let exp = (n / p) - (k / p) - (m / p);
            if exp > 0 {
                push_factor(p as u64);
            }
        } else if p > m {
            // For p > n / 2 >= k, only primes strictly above n - k divide C(n, k), with exponent 1.
            push_factor(p as u64);
        }
    }
    if block > 1 || leaves.is_empty() {
        leaves.push(block);
    }
    BigInt::from(product_tree_u128(&leaves))
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

/// Nonnegative greatest common divisor of two arbitrary-precision integers.
#[must_use]
pub fn bigint_gcd(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    if let (Some(ua), Some(ub)) = (a.magnitude().to_u128(), b.magnitude().to_u128()) {
        return BigInt::from(gcd_u128(ua, ub));
    }
    BigInt::from(a.magnitude().gcd(b.magnitude()))
}

/// Extended greatest common divisor `(g, x, y)` of arbitrary-precision integers
/// such that `g = gcd(a, b) >= 0` and `a * x + b * y == g`.
#[must_use]
pub fn bigint_xgcd(
    a: &BigInt,
    b: &BigInt,
) -> (BigInt, BigInt, BigInt) {
    if let (Some(ax), Some(by)) = (a.to_i64(), b.to_i64()) {
        let (g, x, y) = super::number_ext::xgcd(ax, by);
        return (BigInt::from(g), BigInt::from(x), BigInt::from(y));
    }
    let (mut r0, mut r1) = (a.abs(), b.abs());
    let (mut s0, mut s1) = (BigInt::one(), BigInt::zero());
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    while !r1.is_zero() {
        let (q, r2) = r0.div_rem(&r1);
        let s2 = &s0 - &q * &s1;
        let t2 = &t0 - &q * &t1;
        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    if a.is_negative() {
        s0 = -s0;
    }
    if b.is_negative() {
        t0 = -t0;
    }
    (r0, s0, t0)
}

/// Modular multiplicative inverse of `a` modulo `modulus` in `[0, modulus)`.
///
/// Returns `None` if `modulus <= 0` or `gcd(a, modulus) != 1` (except `modulus == 1`,
/// which returns `Some(0)`).
#[must_use]
pub fn bigint_mod_inverse(
    a: &BigInt,
    modulus: &BigInt,
) -> Option<BigInt> {
    if !modulus.is_positive() {
        return None;
    }
    if modulus.is_one() {
        return Some(BigInt::zero());
    }
    if let (Some(av), Some(mv)) = (a.to_i64(), modulus.to_i64()) {
        return mod_inverse(av, mv).map(BigInt::from);
    }
    let mut r1 = a % modulus;
    if r1.is_negative() {
        r1 += modulus;
    }
    if r1.is_zero() {
        return None;
    }
    let mut r0 = modulus.clone();
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    while !r1.is_zero() {
        let (q, r2) = r0.div_rem(&r1);
        let t2 = &t0 - &q * &t1;
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    if !r0.is_one() {
        return None;
    }
    if t0.is_negative() {
        t0 += modulus;
    }
    Some(t0)
}

/// Floor integer square root `floor(sqrt(n))` for nonnegative `n`, or `None` for `n < 0`.
#[must_use]
pub fn bigint_isqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if let Some(v) = n.to_u64() {
        return Some(BigInt::from(isqrt(v)));
    }
    Some(n.sqrt())
}

/// Truncated integer `k`-th root of `n` (`Some(r)` with `r^k` closest to `n`
/// toward zero without exceeding `|n|` in magnitude).
///
/// Returns `None` when `k == 0` or when `k` is even and `n < 0`.
#[must_use]
pub fn bigint_nth_root(
    n: &BigInt,
    k: u32,
) -> Option<BigInt> {
    if k == 0 {
        return None;
    }
    if k == 1 {
        return Some(n.clone());
    }
    if n.is_negative() && k.is_multiple_of(2) {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    if u64::from(k) >= n.magnitude().bits() {
        return Some(if n.is_negative() {
            BigInt::from(-1_i8)
        } else {
            BigInt::one()
        });
    }
    Some(n.nth_root(k))
}

/// `a^(-1) mod m`, from [`crate::kernels::modular::mod_inv_u64`], with the
/// historical convention that `m == 1` has the unique residue `0`.
const fn inv_mod_u64(
    a: u64,
    m: u64,
) -> Option<u64> {
    if m == 1 {
        return Some(0);
    }
    crate::kernels::modular::mod_inv_u64(a, m)
}

fn garner_reconstruct_tree(
    coeffs: &[usize],
    primes: &[u64],
) -> (BigUint, BigUint) {
    let len = coeffs.len();
    if len <= 32 {
        let mut acc = BigUint::from(coeffs[len - 1] as u64);
        let mut prod = BigUint::from(primes[len - 1]);
        for i in (0..len - 1).rev() {
            let p_i = primes[i];
            acc *= p_i;
            acc += coeffs[i] as u64;
            prod *= p_i;
        }
        return (acc, prod);
    }
    let mid = len / 2;
    let (c_left, c_right) = coeffs.split_at(mid);
    let (p_left, p_right) = primes.split_at(mid);
    let ((val_l, mod_l), (val_r, mod_r)) = if crate::execution::should_split_work(len, 64) {
        let token = crate::cancellation::current();
        let left_token = token.clone();
        rayon::join(
            move || {
                let _guard = crate::cancellation::install_optional(left_token);
                garner_reconstruct_tree(c_left, p_left)
            },
            move || {
                let _guard = crate::cancellation::install_optional(token);
                garner_reconstruct_tree(c_right, p_right)
            },
        )
    } else {
        (
            garner_reconstruct_tree(c_left, p_left),
            garner_reconstruct_tree(c_right, p_right),
        )
    };
    (val_l + val_r * &mod_l, mod_l * mod_r)
}

fn garner_reconstruct_u64(
    coeffs: &[u64],
    primes: &[u64],
) -> BigInt {
    let mut acc = BigUint::from(coeffs[coeffs.len() - 1]);
    for i in (0..coeffs.len() - 1).rev() {
        acc *= primes[i];
        acc += coeffs[i];
    }
    BigInt::from(acc)
}

fn garner_step_with(
    i: usize,
    rem: u64,
    p_i: u64,
    primes: &[u64],
    coeff_at: impl Fn(usize) -> u64,
) -> Option<u64> {
    if p_i == 0 {
        return None;
    }
    let r_i = rem % p_i;
    if i == 0 {
        return Some(r_i);
    }
    if p_i == 1 {
        return Some(0);
    }
    if let Some(mont) = MontU64::new(p_i) {
        let mut u_m = mont.to_mont(coeff_at(i - 1));
        let mut prod_m = mont.to_mont(primes[i - 1]);
        for j in (0..i - 1).rev() {
            let pj_m = mont.to_mont(primes[j]);
            u_m = mont.add(mont.mul_mont(u_m, pj_m), mont.to_mont(coeff_at(j)));
            prod_m = mont.mul_mont(prod_m, pj_m);
        }
        if prod_m == 0 {
            return None;
        }
        let inv_m = if is_prime_miller_rabin(p_i) {
            mont.pow_mont(prod_m, p_i - 2)
        } else {
            let prod = mont.from_mont(prod_m);
            let inv = inv_mod_u64(prod, p_i)?;
            mont.to_mont(inv)
        };
        let diff_m = mont.sub(mont.to_mont(r_i), u_m);
        Some(mont.from_mont(mont.mul_mont(diff_m, inv_m)))
    } else {
        let mut u = coeff_at(i - 1) % p_i;
        let mut prod = primes[i - 1] % p_i;
        for j in (0..i - 1).rev() {
            let pj = primes[j] % p_i;
            let term = mul_mod(u, pj, p_i);
            let cj = coeff_at(j) % p_i;
            u = ((u128::from(term) + u128::from(cj)) % u128::from(p_i)) as u64;
            prod = mul_mod(prod, pj, p_i);
        }
        let inv = inv_mod_u64(prod, p_i)?;
        let diff = if r_i >= u {
            r_i - u
        } else {
            p_i - (u - r_i)
        };
        Some(mul_mod(diff, inv, p_i))
    }
}

/// Solves `x ≡ remainders[i] (mod primes[i])` via Garner's mixed-radix algorithm,
/// returning the unique canonical solution `x` in `[0, prod(primes))`.
///
/// Uses pooled scratch storage from [`global_pool`] and Montgomery modular arithmetic
/// (`MontU64`) for odd moduli. Returns `None` if the slices are empty, have mismatched
/// lengths, contain a zero modulus, or contain non-coprime moduli.
#[must_use]
pub fn bigint_crt_garner(
    remainders: &[u64],
    primes: &[u64],
) -> Option<BigInt> {
    if remainders.is_empty() || remainders.len() != primes.len() {
        return None;
    }
    let k = primes.len();
    let fits_usize = primes.iter().all(|&p| usize::try_from(p).is_ok());
    if !fits_usize {
        let mut coeffs = Vec::with_capacity(k);
        for (i, (&rem, &p_i)) in remainders.iter().zip(primes.iter()).enumerate() {
            let v_i = garner_step_with(i, rem, p_i, primes, |j| coeffs[j])?;
            coeffs.push(v_i);
        }
        return Some(garner_reconstruct_u64(&coeffs, primes));
    }

    let mut coeffs = global_pool().scratch_usize(k);
    for (i, (&rem, &p_i)) in remainders.iter().zip(primes.iter()).enumerate() {
        let v_i = garner_step_with(i, rem, p_i, primes, |j| coeffs[j] as u64)?;
        coeffs.push(v_i as usize);
    }

    if k <= 64 {
        let mut acc = BigUint::from(coeffs[k - 1] as u64);
        for i in (0..k - 1).rev() {
            acc *= primes[i];
            acc += coeffs[i] as u64;
        }
        Some(BigInt::from(acc))
    } else {
        let (acc, _) = garner_reconstruct_tree(&coeffs, primes);
        Some(BigInt::from(acc))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    #[test]
    fn word_product_and_large_binomial_match_exact_oracles() {
        let mut expected = BigInt::one();
        for n in 1..=3_000_i64 {
            expected *= n;
            if matches!(n, 1 | 31 | 32 | 64 | 500 | 2_048 | 2_500 | 3_000) {
                assert_eq!(permutations(n), expected);
            }
        }
        assert_eq!(permutations(0), BigInt::one());
        assert_eq!(permutations(-5), BigInt::zero());

        // Exercises the Kummer/Legendre prime-power binomial branch (k > 256 && n <= 2_000_000).
        for (n, k) in [(600_i64, 300_i64), (1_000, 350), (2_500, 1_250)] {
            let direct = word_product((n - k + 1) as u64, n as u64) / word_product(2, k as u64);
            assert_eq!(binomial(n, k), direct);
            assert_eq!(binomial(n, n - k), direct);
        }
        // Large binomial beyond 2_000_000 falls back to parallel word_product ratio cleanly.
        assert_eq!(
            binomial(3_000_000, 300),
            binomial(3_000_000, 3_000_000 - 300)
        );
        assert!(binomial(100_000, 50_000) > BigInt::zero());
    }

    #[test]
    fn parallel_word_product_is_deterministic_and_inherits_cancellation() {
        use std::sync::Arc;
        use std::sync::atomic::AtomicBool;

        let expected = word_product(2, 5_000);
        for threads in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                assert_eq!(word_product(2, 5_000), expected);
                let token = Arc::new(AtomicBool::new(true));
                let _guard = crate::cancellation::install(token);
                assert_eq!(word_product(2, 5_000), expected);
                assert!(crate::cancellation::is_cancelled());
            });
        }
    }

    #[test]
    fn bigint_gcd_xgcd_and_mod_inverse_satisfy_identities() {
        let a = BigInt::from(10_u8).pow(40) * 17;
        let b = BigInt::from(10_u8).pow(35) * 34;
        let (g, x, y) = bigint_xgcd(&a, &b);
        assert_eq!(g, bigint_gcd(&a, &b));
        assert_eq!(&a * &x + &b * &y, g);

        for sa in [-1, 1] {
            for sb in [-1, 1] {
                let aa = &a * sa;
                let bb = &b * sb;
                let (gg, xx, yy) = bigint_xgcd(&aa, &bb);
                assert_eq!(gg, bigint_gcd(&aa, &bb));
                assert!(gg >= BigInt::zero());
                assert_eq!(&aa * &xx + &bb * &yy, gg);
            }
        }
        assert_eq!(bigint_gcd(&BigInt::zero(), &(-a.clone())), a);

        let p = (BigInt::one() << 127) - 1;
        let val = (BigInt::one() << 90) + 123_456_789;
        let inv = bigint_mod_inverse(&val, &p).unwrap();
        assert!(inv >= BigInt::zero() && inv < p);
        assert_eq!((&val * &inv) % &p, BigInt::one());

        let neg_val = -(&val);
        let neg_inv = bigint_mod_inverse(&neg_val, &p).unwrap();
        assert_eq!((&neg_val * &neg_inv).mod_floor(&p), BigInt::one());

        assert_eq!(bigint_mod_inverse(&BigInt::from(6), &BigInt::from(9)), None);
        assert_eq!(
            bigint_mod_inverse(&BigInt::from(5), &BigInt::one()),
            Some(BigInt::zero())
        );
        assert_eq!(bigint_mod_inverse(&BigInt::from(5), &BigInt::zero()), None);
        assert_eq!(
            bigint_mod_inverse(&BigInt::from(5), &BigInt::from(-7)),
            None
        );
    }

    #[test]
    fn bigint_roots_and_crt_garner_match_exact_invariants() {
        assert_eq!(bigint_isqrt(&BigInt::from(-1)), None);
        assert_eq!(bigint_isqrt(&BigInt::zero()), Some(BigInt::zero()));
        let big: BigInt = (BigInt::one() << 200) + 987_654_321;
        let sq = &big * &big + 42;
        assert_eq!(bigint_isqrt(&sq), Some(big.clone()));

        assert_eq!(bigint_nth_root(&BigInt::from(8), 0), None);
        assert_eq!(bigint_nth_root(&BigInt::from(-16), 2), None);
        assert_eq!(
            bigint_nth_root(&BigInt::from(-27), 3),
            Some(BigInt::from(-3))
        );
        assert_eq!(
            bigint_nth_root(&BigInt::from(-28), 3),
            Some(BigInt::from(-3))
        );
        assert_eq!(
            bigint_nth_root(&BigInt::from(-1), 3),
            Some(BigInt::from(-1))
        );
        assert_eq!(
            bigint_nth_root(&BigInt::from(-2), 3),
            Some(BigInt::from(-1))
        );
        let fifth = big.pow(5) + 17;
        assert_eq!(bigint_nth_root(&fifth, 5), Some(big));

        // Small Garner CRT including prime 2, composite coprime, and rejection cases.
        assert_eq!(bigint_crt_garner(&[], &[]), None);
        assert_eq!(bigint_crt_garner(&[1], &[2, 3]), None);
        assert_eq!(bigint_crt_garner(&[1, 2], &[4, 6]), None);
        assert_eq!(bigint_crt_garner(&[1], &[0]), None);
        assert_eq!(
            bigint_crt_garner(&[1, 2, 4], &[2, 3, 5]),
            Some(BigInt::from(29))
        );
        assert_eq!(bigint_crt_garner(&[3, 5], &[4, 9]), Some(BigInt::from(23)));

        // Large Garner CRT exercising > 64 primes and balanced tree reconstruction.
        let small_primes: Vec<u64> = primes_sieve(600).into_iter().map(|p| p as u64).collect();
        assert!(small_primes.len() > 80);
        let secret: BigInt = (BigInt::one() << 250) + BigInt::from(987_654_321_123_456_789_u64);
        let rems: Vec<u64> = small_primes
            .iter()
            .map(|&p| (&secret % BigInt::from(p)).to_u64().unwrap())
            .collect();
        let solved = bigint_crt_garner(&rems, &small_primes).unwrap();
        assert_eq!(solved, secret);
    }
}
