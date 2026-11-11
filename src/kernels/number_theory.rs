//! # Numerical Number Theory
//!
//! Machine-word integer number theory: greatest common divisors, modular
//! arithmetic, deterministic primality testing, factorisation (trial
//! division plus Pollard–Brent rho), square roots modulo a prime
//! (Tonelli–Shanks) and discrete logarithms (baby-step giant-step).
//!
//! None of the functions panic on any input: undefined results are `None`
//! (or `0` where the signature predates `Option` and says so), and every
//! intermediate product is carried out in `u128`.

/// Computes the greatest common divisor (GCD) of two numbers using the Euclidean algorithm.
///
/// # Arguments
/// * `a` - The first number.
/// * `b` - The second number.
///
/// # Returns
/// The greatest common divisor of `a` and `b`.
#[must_use]
pub const fn gcd(
    mut a: u64,
    mut b: u64,
) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// `a * b mod m` without overflow. `m` must be non-zero.
#[inline]
const fn mul_mod(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Computes `(base^exp) % modulus` efficiently using modular exponentiation (binary exponentiation).
///
/// `modulus == 0` has no meaning; the function returns `0` for it instead of
/// dividing by zero. `x^0 mod 1` is `0`.
///
/// # Arguments
/// * `base` - The base.
/// * `exp` - The exponent.
/// * `modulus` - The modulus.
///
/// # Returns
/// The result of `(base^exp) % modulus`.
#[must_use]
pub const fn mod_pow(
    base: u128,
    mut exp: u64,
    modulus: u64,
) -> u64 {
    if modulus == 0 {
        return 0;
    }
    let mut res = 1 % modulus;
    let mut base = (base % modulus as u128) as u64;
    while exp > 0 {
        if exp & 1 == 1 {
            res = mul_mod(res, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exp >>= 1;
    }
    res
}

/// Finds the modular multiplicative inverse of a number.
///
/// Solves for `x` in `ax ≡ 1 (mod m)` using the extended Euclidean
/// algorithm in `i128`. The modulus may be negative (its absolute value is
/// used); `m == 0` has no inverse. The result lies in `0..|m|`, and
/// `|m| == 1` gives `Some(0)`.
///
/// # Returns
/// An `Option<i64>` containing the modular inverse if it exists, otherwise `None`.
#[must_use]
pub fn mod_inverse(
    a: i64,
    m: i64,
) -> Option<i64> {
    if m == 0 {
        return None;
    }
    let m = i128::from(m).abs();
    // Iterative extended Euclid on (a mod m, m); coefficients stay below m.
    let (mut r0, mut r1) = (m, i128::from(a).rem_euclid(m));
    let (mut t0, mut t1) = (0_i128, 1_i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    if r0 != 1 {
        // gcd(a, m) != 1; but m == 1 has the unique residue 0.
        return if m == 1 { Some(0) } else { None };
    }
    i64::try_from(t0.rem_euclid(m)).ok()
}

/// Performs a Miller-Rabin primality test.
///
/// The twelve prime bases 2..=37 are a deterministic witness set for every
/// `n < 3.3e24`, hence for all of `u64`. Small primes are recognised by
/// trial division first.
///
/// # Arguments
/// * `n` - The number to test for primality.
///
/// # Returns
/// `true` if `n` is prime, `false` otherwise.
#[must_use]
pub const fn is_prime_miller_rabin(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    let mut i = 0;
    while i < BASES.len() {
        if n == BASES[i] {
            return true;
        }
        if n % BASES[i] == 0 {
            return false;
        }
        i += 1;
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    let mut i = 0;
    while i < BASES.len() {
        let mut x = mod_pow(BASES[i] as u128, d, n);
        if x != 1 && x != n - 1 {
            let mut composite = true;
            let mut r = 1;
            while r < s {
                x = mul_mod(x, x, n);
                if x == n - 1 {
                    composite = false;
                    break;
                }
                r += 1;
            }
            if composite {
                return false;
            }
        }
        i += 1;
    }
    true
}

/// Computes the least common multiple (LCM) of two numbers.
///
/// Returns `0` when either argument is `0`, and `u64::MAX` when the true
/// value does not fit (use [`checked_lcm`] to tell the two apart).
#[must_use]
pub fn lcm(
    a: u64,
    b: u64,
) -> u64 {
    checked_lcm(a, b).unwrap_or(u64::MAX)
}

/// The least common multiple, or `None` when it overflows `u64`.
#[must_use]
pub fn checked_lcm(
    a: u64,
    b: u64,
) -> Option<u64> {
    if a == 0 || b == 0 {
        return Some(0);
    }
    (a / gcd(a, b)).checked_mul(b)
}

/// Computes Euler's totient function φ(n) (`φ(0) = 0`).
///
/// Uses [`factorize`], so it is fast for every `u64`.
#[must_use]
pub fn phi(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut result = n;
    let mut last = 0;
    for p in factorize(n) {
        if p != last {
            result -= result / p;
            last = p;
        }
    }
    result
}

/// Pollard–Brent rho: a non-trivial factor of the odd composite `n`.
fn pollard_brent(n: u64) -> Option<u64> {
    for c in 1..64_u64 {
        let f = |x: u64| (mul_mod(x, x, n) + c) % n;
        let (mut x, mut y, mut g, mut q) = (2_u64, 2_u64, 1_u64, 1_u64);
        let (mut r, m) = (1_u64, 128_u64);
        let mut ys = y;
        'outer: while g == 1 {
            x = y;
            for _ in 0..r {
                y = f(y);
            }
            let mut k = 0;
            while k < r && g == 1 {
                ys = y;
                for _ in 0..m.min(r - k) {
                    y = f(y);
                    q = mul_mod(q, x.abs_diff(y), n);
                }
                g = gcd(q, n);
                k += m;
            }
            r *= 2;
            if r > 1 << 40 {
                break 'outer;
            }
        }
        if g == n {
            // The batch overshot: replay it one step at a time.
            loop {
                ys = f(ys);
                g = gcd(x.abs_diff(ys), n);
                if g > 1 {
                    break;
                }
            }
        }
        if g != n && g != 1 {
            return Some(g);
        }
    }
    None
}

/// Returns the prime factors of `n` with multiplicity, in ascending order.
///
/// `0` and `1` have no factors. Small primes are removed by trial division;
/// the cofactor is split with Miller–Rabin and Pollard–Brent rho, so every
/// `u64` factors in milliseconds.
#[must_use]
pub fn factorize(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    if n < 2 {
        return factors;
    }
    for p in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] {
        while n % p == 0 {
            factors.push(p);
            n /= p;
        }
    }
    let mut stack = if n > 1 { vec![n] } else { Vec::new() };
    while let Some(m) = stack.pop() {
        if m == 1 {
            continue;
        }
        if is_prime_miller_rabin(m) {
            factors.push(m);
            continue;
        }
        // Perfect squares defeat rho with some polynomials; peel them off.
        let r = isqrt(m);
        if r * r == m {
            stack.push(r);
            stack.push(r);
            continue;
        }
        match pollard_brent(m) {
            | Some(d) => {
                stack.push(d);
                stack.push(m / d);
            },
            // Unreachable in practice; fall back to trial division.
            | None => {
                let mut i = 53_u64;
                let mut rest = m;
                while i <= rest / i {
                    while rest % i == 0 {
                        factors.push(i);
                        rest /= i;
                    }
                    i += 2;
                }
                if rest > 1 {
                    factors.push(rest);
                }
            },
        }
    }
    factors.sort_unstable();
    factors
}

/// Floor of the square root of `n`, exact for all `u64`.
#[must_use]
pub fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r.checked_mul(r).is_none_or(|sq| sq > n) {
        r -= 1;
    }
    while (r + 1).checked_mul(r + 1).is_some_and(|sq| sq <= n) {
        r += 1;
    }
    r
}

/// Generates primes up to limit using Sieve of Eratosthenes.
///
/// Returns an empty vector when the sieve table cannot be allocated
/// (`limit` near `usize::MAX`) instead of aborting.
#[must_use]
pub fn primes_sieve(limit: usize) -> Vec<usize> {
    if limit < 2 {
        return Vec::new();
    }
    let Some(len) = limit.checked_add(1) else {
        return Vec::new();
    };
    let mut is_prime: Vec<bool> = Vec::new();
    if is_prime.try_reserve_exact(len).is_err() {
        return Vec::new();
    }
    is_prime.resize(len, true);
    is_prime[0] = false;
    is_prime[1] = false;
    let mut p = 2_usize;
    while p <= limit / p {
        if is_prime[p] {
            let mut i = p * p;
            while i <= limit {
                is_prime[i] = false;
                i += p;
            }
        }
        p += 1;
    }
    is_prime
        .into_iter()
        .enumerate()
        .filter(|&(_, p)| p)
        .map(|(i, _)| i)
        .collect()
}

/// A square root of `n` modulo the prime `p` (Tonelli–Shanks): some `x` with
/// `x² ≡ n (mod p)`, the smaller of the two roots.
///
/// `None` when `p` is not an odd prime or `n` is a non-residue.
#[must_use]
pub fn sqrt_mod_prime(
    n: u64,
    p: u64,
) -> Option<u64> {
    if p == 2 {
        return Some(n & 1);
    }
    if p < 3 || !is_prime_miller_rabin(p) {
        return None;
    }
    let n = n % p;
    if n == 0 {
        return Some(0);
    }
    if mod_pow(u128::from(n), (p - 1) / 2, p) != 1 {
        return None;
    }
    let s = (p - 1).trailing_zeros();
    let q = (p - 1) >> s;
    // Any non-residue.
    let mut z = 2_u64;
    while mod_pow(u128::from(z), (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = mod_pow(u128::from(z), q, p);
    let mut t = mod_pow(u128::from(n), q, p);
    let mut r = mod_pow(u128::from(n), q.div_ceil(2), p);
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mul_mod(tt, tt, p);
            i += 1;
        }
        let mut b = c;
        for _ in 0..m - i - 1 {
            b = mul_mod(b, b, p);
        }
        m = i;
        c = mul_mod(b, b, p);
        t = mul_mod(t, c, p);
        r = mul_mod(r, b, p);
    }
    Some(r.min(p - r))
}

/// The discrete logarithm: the least `x >= 0` with `g^x ≡ h (mod p)`
/// (baby-step giant-step, `O(sqrt p)` time and memory), or `None` when none
/// exists or `p` is too large for the table (`p > 2^48`).
#[must_use]
pub fn discrete_log(
    g: u64,
    h: u64,
    p: u64,
) -> Option<u64> {
    if p < 2 || p > 1 << 48 {
        return None;
    }
    let (g, h) = (g % p, h % p);
    if h == 1 % p {
        return Some(0);
    }
    let m = isqrt(p) + 1;
    let mut table = std::collections::HashMap::with_capacity(m as usize);
    let mut e = 1 % p;
    for j in 0..m {
        table.entry(e).or_insert(j);
        e = mul_mod(e, g, p);
    }
    // g^{-m}, if g is invertible; otherwise only exponents below m can match.
    let factor = mod_inverse(i64::try_from(g).ok()?, i64::try_from(p).ok()?)
        .map(|inv| mod_pow(u128::from(inv as u64), m, p));
    let mut gamma = h;
    for i in 0..m {
        if let Some(&j) = table.get(&gamma) {
            return Some(i * m + j);
        }
        gamma = mul_mod(gamma, factor?, p);
    }
    None
}
