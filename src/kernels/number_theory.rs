//! # Numerical Number Theory
//!
//! Machine-word integer number theory: greatest common divisors, division-free
//! Montgomery modular arithmetic, deterministic primality testing (`u64` and
//! `u128` Miller–Rabin and Baillie–PSW), factorisation (trial division plus
//! Montgomery Pollard–Brent rho), Jacobi/Legendre symbols, arithmetic
//! functions ($\phi$, $\lambda$, $\tau$, $\sigma_k$, divisors, primitive
//! roots), square roots modulo a prime (Tonelli–Shanks) and discrete
//! logarithms for prime or composite moduli.
//!
//! None of the functions panic on any input: undefined results are `None`
//! (or `0` where the signature predates `Option` and says so), and every
//! intermediate product is carried out in `u128` or wider precision.

#![allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn
)]

/// Division-free Montgomery arithmetic for an odd modulus.
///
/// The implementation is [`crate::kernels::modular::Montgomery`]; the name is
/// kept because this module's algorithms (Miller–Rabin, Pollard–Brent, Tonelli–
/// Shanks) were written against it.
pub use crate::kernels::modular::Montgomery as MontU64;

/// Computes the greatest common divisor (GCD) of two numbers using Stein's binary GCD algorithm
/// after one initial Euclidean reduction step.
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
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    a %= b;
    if a == 0 {
        return b;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    while b != 0 {
        b >>= b.trailing_zeros();
        if a > b {
            let t = a;
            a = b;
            b = t;
        }
        b -= a;
    }
    a << shift
}

/// Computes the greatest common divisor (GCD) of two 128-bit unsigned integers
/// using Stein's binary GCD algorithm.
#[must_use]
pub const fn gcd_u128(
    mut a: u128,
    mut b: u128,
) -> u128 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    if a <= u64::MAX as u128 && b <= u64::MAX as u128 {
        return gcd(a as u64, b as u64) as u128;
    }
    a %= b;
    if a == 0 {
        return b;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    while b != 0 {
        b >>= b.trailing_zeros();
        if a <= u64::MAX as u128 && b <= u64::MAX as u128 {
            return (gcd(a as u64, b as u64) as u128) << shift;
        }
        if a > b {
            let t = a;
            a = b;
            b = t;
        }
        b -= a;
    }
    a << shift
}

/// `a * b mod m` without overflow. `m` must be non-zero.
#[must_use]
#[inline]
pub const fn mul_mod(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    let product = a as u128 * b as u128;
    if product <= u64::MAX as u128 {
        product as u64 % m
    } else {
        (product % m as u128) as u64
    }
}

/// Computes `(a * b) % m` for 128-bit unsigned integers without overflow.
///
/// Returns `0` when `m <= 1`.
#[must_use]
pub const fn mul_mod_u128(
    mut a: u128,
    mut b: u128,
    m: u128,
) -> u128 {
    if m <= 1 {
        return 0;
    }
    if a >= m {
        a %= m;
    }
    if b >= m {
        b %= m;
    }
    if a == 0 || b == 0 {
        return 0;
    }
    if let Some(prod) = a.checked_mul(b) {
        return prod % m;
    }
    if a < b {
        let t = a;
        a = b;
        b = t;
    }
    let mut res = 0_u128;
    while b != 0 {
        if (b & 1) == 1 {
            res = if res >= m - a {
                res - (m - a)
            } else {
                res + a
            };
        }
        b >>= 1;
        if b != 0 {
            a = if a >= m - a {
                a - (m - a)
            } else {
                a << 1
            };
        }
    }
    res
}

/// Computes `(base^exp) % modulus` for 128-bit unsigned integers without overflow.
///
/// Returns `0` when `modulus <= 1`.
#[must_use]
pub const fn mod_pow_u128(
    mut base: u128,
    mut exp: u128,
    modulus: u128,
) -> u128 {
    if modulus <= 1 {
        return 0;
    }
    if exp == 0 {
        return 1;
    }
    if modulus <= u64::MAX as u128 && exp <= u64::MAX as u128 {
        return mod_pow(base, exp as u64, modulus as u64) as u128;
    }
    base %= modulus;
    let mut res = 1_u128;
    while exp > 1 {
        if (exp & 1) == 1 {
            res = mul_mod_u128(res, base, modulus);
        }
        base = mul_mod_u128(base, base, modulus);
        exp >>= 1;
    }
    if exp == 1 {
        mul_mod_u128(res, base, modulus)
    } else {
        res
    }
}

/// Computes `(base^exp) % modulus` efficiently using modular exponentiation
/// (binary exponentiation).
///
/// Uses division-free Montgomery multiplication (`MontU64`) when `modulus > 1` is odd.
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
    if modulus <= 1 {
        return 0;
    }
    if exp == 0 {
        return 1;
    }
    if let Some(mont) = MontU64::new(modulus) {
        let base_m = mont.to_mont_u128(base);
        return mont.from_mont(mont.pow_mont(base_m, exp));
    }
    let mut res = 1;
    let mut base = if base < modulus as u128 {
        base as u64
    } else {
        (base % modulus as u128) as u64
    };
    while exp > 1 {
        if exp & 1 == 1 {
            res = mul_mod(res, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exp >>= 1;
    }
    if exp == 1 {
        mul_mod(res, base, modulus)
    } else {
        res
    }
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

/// Performs a deterministic Miller-Rabin primality test over `u64` using
/// division-free Montgomery multiplication (`MontU64`).
///
/// The first twelve prime bases are a deterministic witness set over the
/// entire `u64` range (and smaller proven subsets are used for smaller `n`).
/// Divisibility by those primes is checked first.
///
/// # Arguments
/// * `n` - The number to test for primality.
///
/// # Returns
/// `true` if `n` is prime, `false` otherwise.
#[must_use]
pub const fn is_prime_miller_rabin(n: u64) -> bool {
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if n < 2 {
        return false;
    }
    let mut i = 0;
    while i < BASES.len() {
        if n == BASES[i] {
            return true;
        }
        if n.is_multiple_of(BASES[i]) {
            return false;
        }
        i += 1;
    }
    if n < 41 * 41 {
        return true;
    }
    let Some(mont) = MontU64::new(n) else {
        return false;
    };
    let one_m = mont.one();
    let minus_one_m = n - one_m;
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    let witness_count = if n < 3_215_031_751 {
        4
    } else if n < 3_825_123_056_546_413_051 {
        9
    } else {
        12
    };
    let mut i = 0;
    while i < witness_count {
        let base_m = mont.to_mont(BASES[i]);
        let mut x = mont.pow_mont(base_m, d);
        if x != one_m && x != minus_one_m {
            let mut composite = true;
            let mut r = 1;
            while r < s {
                x = mont.mul_mont(x, x);
                if x == minus_one_m {
                    composite = false;
                    break;
                }
                if x == one_m {
                    return false;
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

/// Computes the Jacobi symbol $\left(\frac{a}{n}\right)$ for an odd positive integer `n`.
///
/// Returns `0` when `n == 0`, `n` is even, or $\gcd(a, n) > 1$. For `n == 1`, returns `1`.
#[must_use]
pub const fn jacobi_symbol(
    a: i128,
    n: u64,
) -> i32 {
    if n == 0 || (n & 1) == 0 {
        return 0;
    }
    if n == 1 {
        return 1;
    }
    let mut a_u = a.rem_euclid(n as i128) as u64;
    let mut n_u = n;
    let mut result = 1_i32;
    while a_u != 0 {
        let zeros = a_u.trailing_zeros();
        a_u >>= zeros;
        if (zeros & 1) == 1 {
            let r = n_u & 7;
            if r == 3 || r == 5 {
                result = -result;
            }
        }
        let t = a_u;
        a_u = n_u;
        n_u = t;
        if (a_u & 3) == 3 && (n_u & 3) == 3 {
            result = -result;
        }
        a_u %= n_u;
    }
    if n_u == 1 { result } else { 0 }
}

/// Computes the Legendre symbol $\left(\frac{a}{p}\right)$ for an odd prime `p`.
///
/// Returns `0` if `p < 3`, `p` is composite, or `p` divides `a`; `1` if `a` is
/// a quadratic residue modulo `p`; and `-1` if `a` is a quadratic non-residue.
#[must_use]
pub const fn legendre_symbol(
    a: u64,
    p: u64,
) -> i32 {
    if p < 3 || !is_prime_miller_rabin(p) {
        return 0;
    }
    jacobi_symbol(a as i128, p)
}

/// Deterministic Baillie–PSW primality test over `u64` (strong base-2 Miller–Rabin
/// test combined with a strong Lucas probable-prime test using Selfridge's parameters,
/// executed in the Montgomery domain).
#[must_use]
pub const fn is_prime_bpsw(n: u64) -> bool {
    const SMALL_PRIMES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if n < 2 {
        return false;
    }
    let mut i = 0;
    while i < SMALL_PRIMES.len() {
        let p = SMALL_PRIMES[i];
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
        i += 1;
    }
    if n < 41 * 41 {
        return true;
    }
    let Some(mont) = MontU64::new(n) else {
        return false;
    };
    let one_m = mont.one();
    let minus_one_m = n - one_m;
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;

    // 1. Strong base-2 test in Montgomery domain.
    let two_m = mont.add(one_m, one_m);
    let mut x = mont.pow_mont(two_m, d);
    if x != one_m && x != minus_one_m {
        let mut passed = false;
        let mut r = 1;
        while r < s {
            x = mont.mul_mont(x, x);
            if x == minus_one_m {
                passed = true;
                break;
            }
            if x == one_m {
                return false;
            }
            r += 1;
        }
        if !passed {
            return false;
        }
    }

    // 2. Reject perfect squares before Selfridge parameter search.
    let root = isqrt(n);
    if root * root == n {
        return false;
    }

    // 3. Selfridge parameters: D in 5, -7, 9, -11, ... with (D/n) == -1.
    let mut d_param: i64 = 5;
    loop {
        let j = jacobi_symbol(d_param as i128, n);
        if j == -1 {
            break;
        }
        if j == 0 && d_param.unsigned_abs() < n {
            return false;
        }
        d_param = if d_param > 0 {
            -(d_param + 2)
        } else {
            -d_param + 2
        };
    }

    // 4. Strong Lucas test with P = 1, Q = (1 - D) / 4 in Montgomery domain.
    let d_mod = (d_param as i128).rem_euclid(n as i128) as u64;
    let q_param = (1 - d_param) / 4;
    let q_mod = (q_param as i128).rem_euclid(n as i128) as u64;
    let d_m = mont.to_mont(d_mod);
    let q_m = mont.to_mont(q_mod);

    let n_plus_1 = (n as u128) + 1;
    let s_lucas = n_plus_1.trailing_zeros();
    let odd_lucas = (n_plus_1 >> s_lucas) as u64;

    let mut u_m = one_m;
    let mut v_m = one_m;
    let mut qk_m = q_m;
    let top_bit = 63 - odd_lucas.leading_zeros();
    let mut bit = top_bit;
    while bit > 0 {
        bit -= 1;
        // Doubling: U_{2k} = U_k * V_k, V_{2k} = V_k^2 - 2 Q^k, Q^{2k} = (Q^k)^2.
        u_m = mont.mul_mont(u_m, v_m);
        let v_sq = mont.mul_mont(v_m, v_m);
        let two_qk = mont.add(qk_m, qk_m);
        v_m = mont.sub(v_sq, two_qk);
        qk_m = mont.mul_mont(qk_m, qk_m);

        if ((odd_lucas >> bit) & 1) == 1 {
            let next_u = mont.half(mont.add(u_m, v_m));
            let du = mont.mul_mont(d_m, u_m);
            let next_v = mont.half(mont.add(du, v_m));
            u_m = next_u;
            v_m = next_v;
            qk_m = mont.mul_mont(qk_m, q_m);
        }
    }

    if u_m == 0 || v_m == 0 {
        return true;
    }
    let mut r = 1;
    while r < s_lucas {
        let v_sq = mont.mul_mont(v_m, v_m);
        let two_qk = mont.add(qk_m, qk_m);
        v_m = mont.sub(v_sq, two_qk);
        if v_m == 0 {
            return true;
        }
        qk_m = mont.mul_mont(qk_m, qk_m);
        r += 1;
    }
    false
}

/// Primality test for 128-bit unsigned integers combining deterministic `u64`
/// Miller–Rabin with 13-base strong Miller–Rabin and a strong Selfridge Lucas
/// test (Baillie–PSW) for values above `u64::MAX`.
#[must_use]
pub fn is_prime_u128(n: u128) -> bool {
    if n <= u128::from(u64::MAX) {
        return is_prime_miller_rabin(n as u64);
    }
    if (n & 1) == 0 {
        return false;
    }
    const SMALL_PRIMES: [u128; 30] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127,
    ];
    for p in SMALL_PRIMES {
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    const MR_BASES: [u128; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
    for base in MR_BASES {
        let mut x = mod_pow_u128(base, d, n);
        if x != 1 && x != n - 1 {
            let mut passed = false;
            for _ in 1..s {
                x = mul_mod_u128(x, x, n);
                if x == n - 1 {
                    passed = true;
                    break;
                }
                if x == 1 {
                    return false;
                }
            }
            if !passed {
                return false;
            }
        }
    }
    strong_lucas_u128(n)
}

const fn isqrt_u128(n: u128) -> u128 {
    if n <= u64::MAX as u128 {
        return isqrt(n as u64) as u128;
    }
    let shift = (128 - n.leading_zeros()).div_ceil(2);
    let mut x = 1_u128 << shift;
    loop {
        let next = (x + n / x) >> 1;
        if next >= x {
            return x;
        }
        x = next;
    }
}

fn jacobi_u128(
    a: i128,
    n: u128,
) -> i32 {
    if n == 0 || (n & 1) == 0 {
        return 0;
    }
    if n <= u128::from(u64::MAX) {
        return jacobi_symbol(a, n as u64);
    }
    let mut a_u = if a >= 0 {
        (a as u128) % n
    } else {
        let mag = a.unsigned_abs() % n;
        if mag == 0 { 0 } else { n - mag }
    };
    let mut n_u = n;
    let mut result = 1_i32;
    while a_u != 0 {
        let zeros = a_u.trailing_zeros();
        a_u >>= zeros;
        if (zeros & 1) == 1 {
            let r = n_u & 7;
            if r == 3 || r == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a_u, &mut n_u);
        if (a_u & 3) == 3 && (n_u & 3) == 3 {
            result = -result;
        }
        a_u %= n_u;
    }
    if n_u == 1 { result } else { 0 }
}

const fn add_mod_u128(
    a: u128,
    b: u128,
    n: u128,
) -> u128 {
    if a >= n - b {
        a - (n - b)
    } else {
        a + b
    }
}

const fn sub_mod_u128(
    a: u128,
    b: u128,
    n: u128,
) -> u128 {
    if a >= b { a - b } else { n - (b - a) }
}

const fn half_mod_u128(
    x: u128,
    n: u128,
) -> u128 {
    if (x & 1) == 0 {
        x >> 1
    } else {
        (x >> 1) + (n >> 1) + 1
    }
}

fn strong_lucas_u128(n: u128) -> bool {
    let r = isqrt_u128(n);
    if r * r == n {
        return false;
    }
    let mut d_param: i128 = 5;
    loop {
        let j = jacobi_u128(d_param, n);
        if j == -1 {
            break;
        }
        if j == 0 && d_param.unsigned_abs() < n {
            return false;
        }
        d_param = if d_param > 0 {
            -(d_param + 2)
        } else {
            -d_param + 2
        };
    }
    let d_mod = if d_param >= 0 {
        (d_param as u128) % n
    } else {
        n - (d_param.unsigned_abs() % n)
    };
    let q_param = (1 - d_param) / 4;
    let q_mod = if q_param >= 0 {
        (q_param as u128) % n
    } else {
        n - (q_param.unsigned_abs() % n)
    };
    let (n1, overflow) = n.overflowing_add(1);
    let (s_lucas, odd_lucas) = if overflow {
        (128_u32, 1_u128)
    } else {
        let s = n1.trailing_zeros();
        (s, n1 >> s)
    };

    let mut u = 1_u128;
    let mut v = 1_u128;
    let mut qk = q_mod;
    let top_bit = 127 - odd_lucas.leading_zeros();
    for bit in (0..top_bit).rev() {
        u = mul_mod_u128(u, v, n);
        let v_sq = mul_mod_u128(v, v, n);
        let two_qk = add_mod_u128(qk, qk, n);
        v = sub_mod_u128(v_sq, two_qk, n);
        qk = mul_mod_u128(qk, qk, n);
        if ((odd_lucas >> bit) & 1) == 1 {
            let next_u = half_mod_u128(add_mod_u128(u, v, n), n);
            let du = mul_mod_u128(d_mod, u, n);
            let next_v = half_mod_u128(add_mod_u128(du, v, n), n);
            u = next_u;
            v = next_v;
            qk = mul_mod_u128(qk, q_mod, n);
        }
    }
    if u == 0 || v == 0 {
        return true;
    }
    for _ in 1..s_lucas {
        let v_sq = mul_mod_u128(v, v, n);
        let two_qk = add_mod_u128(qk, qk, n);
        v = sub_mod_u128(v_sq, two_qk, n);
        if v == 0 {
            return true;
        }
        qk = mul_mod_u128(qk, qk, n);
    }
    false
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
pub const fn checked_lcm(
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

/// Computes Carmichael's reduced totient function $\lambda(n)$ (`carmichael_lambda(0) = 0`).
///
/// $\lambda(n)$ is the exponent of the multiplicative group $(\mathbb{Z}/n\mathbb{Z})^\times$,
/// i.e. the smallest positive integer $m$ such that $a^m \equiv 1 \pmod n$ for all
/// $\gcd(a, n) = 1$.
#[must_use]
pub fn carmichael_lambda(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    if n == 1 {
        return 1;
    }
    let factors = factorize(n);
    let mut result = 1_u64;
    let mut idx = 0_usize;
    while idx < factors.len() {
        let p = factors[idx];
        let mut exp = 0_u32;
        let mut pk = 1_u64;
        while idx < factors.len() && factors[idx] == p {
            exp += 1;
            pk *= p;
            idx += 1;
        }
        let lambda_pk = if p == 2 {
            if exp == 1 {
                1
            } else if exp == 2 {
                2
            } else {
                pk >> 2
            }
        } else {
            (pk / p) * (p - 1)
        };
        result = lcm(result, lambda_pk);
    }
    result
}

/// Returns all positive divisors of `n` in ascending order (`divisors(0)` is empty).
#[must_use]
pub fn divisors(n: u64) -> Vec<u64> {
    if n == 0 {
        return Vec::new();
    }
    let factors = factorize(n);
    let mut divs = vec![1_u64];
    let mut idx = 0_usize;
    while idx < factors.len() {
        let p = factors[idx];
        let mut count = 0_usize;
        while idx < factors.len() && factors[idx] == p {
            count += 1;
            idx += 1;
        }
        let prev_len = divs.len();
        let mut p_pow = 1_u64;
        for _ in 0..count {
            p_pow *= p;
            for i in 0..prev_len {
                let d = divs[i] * p_pow;
                divs.push(d);
            }
        }
    }
    divs.sort_unstable();
    divs
}

/// Returns the number of positive divisors $\tau(n) = \sigma_0(n)$ (`divisor_count(0) = 0`).
#[must_use]
pub fn divisor_count(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let factors = factorize(n);
    let mut total = 1_u64;
    let mut idx = 0_usize;
    while idx < factors.len() {
        let p = factors[idx];
        let mut count = 0_u64;
        while idx < factors.len() && factors[idx] == p {
            count += 1;
            idx += 1;
        }
        total = total.saturating_mul(count + 1);
    }
    total
}

/// Returns the sum of the $k$-th powers of the positive divisors of `n`,
/// $\sigma_k(n) = \sum_{d \mid n} d^k$, saturating at `u128::MAX` (`divisor_sigma(0, _) = 0`).
#[must_use]
pub fn divisor_sigma(
    n: u64,
    k: u32,
) -> u128 {
    if n == 0 {
        return 0;
    }
    if k == 0 {
        return u128::from(divisor_count(n));
    }
    let factors = factorize(n);
    let mut total = 1_u128;
    let mut idx = 0_usize;
    while idx < factors.len() {
        let p = u128::from(factors[idx]);
        let mut count = 0_u32;
        while idx < factors.len() && u128::from(factors[idx]) == p {
            count += 1;
            idx += 1;
        }
        let pk = p.checked_pow(k).unwrap_or(u128::MAX);
        let mut term_sum = 1_u128;
        let mut cur = 1_u128;
        for _ in 0..count {
            cur = cur.saturating_mul(pk);
            term_sum = term_sum.saturating_add(cur);
        }
        total = total.saturating_mul(term_sum);
    }
    total
}

/// Returns the smallest positive primitive root modulo `n` if the multiplicative
/// group $(\mathbb{Z}/n\mathbb{Z})^\times$ is cyclic ($n \in \{2, 4, p^k, 2p^k\}$
/// for an odd prime $p$), or `None` otherwise.
#[must_use]
pub fn primitive_root(n: u64) -> Option<u64> {
    if n < 2 {
        return None;
    }
    if n == 2 {
        return Some(1);
    }
    if n == 4 {
        return Some(3);
    }
    let odd_part = if (n & 1) == 0 { n >> 1 } else { n };
    if (odd_part & 1) == 0 {
        return None;
    }
    let factors = factorize(odd_part);
    let p = *factors.first()?;
    if factors.last() != Some(&p) {
        return None;
    }
    let k = factors.len();
    let phi_n = (odd_part / p) * (p - 1);
    let mut prime_divs = factorize(p - 1);
    if k > 1 {
        prime_divs.push(p);
    }
    prime_divs.sort_unstable();
    prime_divs.dedup();

    let mont = MontU64::new(n);
    for g in 2..n {
        if g.is_multiple_of(p) || ((n & 1) == 0 && (g & 1) == 0) {
            continue;
        }
        let is_generator = if let Some(m) = mont {
            let g_m = m.to_mont(g);
            let one_m = m.one();
            prime_divs
                .iter()
                .all(|&q| m.pow_mont(g_m, phi_n / q) != one_m)
        } else {
            prime_divs
                .iter()
                .all(|&q| mod_pow(u128::from(g), phi_n / q, n) != 1)
        };
        if is_generator {
            return Some(g);
        }
    }
    None
}

/// Pollard–Brent rho: a non-trivial factor of the odd composite `n` using
/// division-free Montgomery arithmetic (`MontU64`).
fn pollard_brent(n: u64) -> Option<u64> {
    let mont = MontU64::new(n)?;
    let two_m = mont.to_mont(2);
    let one_m = mont.one();
    for c in 1..64_u64 {
        let c_m = mont.to_mont(c);
        let f = |x_m: u64| mont.add(mont.mul_mont(x_m, x_m), c_m);
        let (mut x, mut y, mut g, mut q) = (two_m, two_m, 1_u64, one_m);
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
                    q = mont.mul_mont(q, mont.sub(x, y));
                }
                // Since gcd(R, n) == 1 for odd n, gcd(q_mont, n) == gcd(q, n).
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
                g = gcd(mont.sub(x, ys), n);
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
    let twos = n.trailing_zeros();
    if twos > 0 {
        factors.resize(twos as usize, 2);
        n >>= twos;
    }
    for p in [3_u64, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] {
        while n.is_multiple_of(p) {
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
                    while rest.is_multiple_of(i) {
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
pub const fn isqrt(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let shift = (64 - n.leading_zeros()).div_ceil(2);
    let mut r = 1_u64 << shift;
    loop {
        let next = (r + n / r) >> 1;
        if next >= r {
            return r;
        }
        r = next;
    }
}

/// Generates primes up to `limit` using the Sieve of Eratosthenes.
///
/// For `limit <= 1_000_000`, reuses the cached prime table in
/// [`crate::kernels::reuse::global_pool`]. For larger limits, executes an
/// odd-only bit-packed L1-cache-resident segmented sieve leasing scratch
/// storage from [`crate::kernels::reuse::global_pool`].
///
/// Returns an empty vector when the sieve table cannot be allocated
/// (`limit` near `usize::MAX`) instead of aborting.
#[must_use]
pub fn primes_sieve(limit: usize) -> Vec<usize> {
    if limit < 2 {
        return Vec::new();
    }
    if limit <= 1_000_000 {
        return crate::kernels::reuse::global_pool()
            .prime_sieve_cache(limit)
            .as_ref()
            .clone();
    }
    // Verify that the result and sieve extent are realistically addressable.
    // Match the traditional byte-per-entry capacity guard so astronomical
    // limits (e.g. usize::MAX / 2) return Vec::new() instead of OOM-aborting.
    if limit > (usize::MAX >> 2) {
        return Vec::new();
    }
    let approx_primes = (limit as f64 / ((limit as f64).ln() - 1.1).max(1.0)) as usize + 64;
    let mut primes = Vec::new();
    if primes.try_reserve(approx_primes).is_err() {
        return Vec::new();
    }
    primes.push(2);

    let root = isqrt(limit as u64) as usize;
    let base_primes = crate::kernels::reuse::global_pool().prime_sieve_cache(root);

    // 32 KiB L1-resident odd-only bit-packed segment (covers 524,288 integers per segment).
    const SEGMENT_BYTES: usize = 32 * 1024;
    const SEGMENT_BITS: usize = SEGMENT_BYTES * 8;
    const SEGMENT_SPAN: usize = SEGMENT_BITS * 2;

    let max_odd = limit - ((limit & 1) ^ 1);
    if max_odd < 3 {
        return primes;
    }
    let num_segments = (max_odd - 3) / SEGMENT_SPAN + 1;

    let sieve_segment = |seg_idx: usize, out: &mut Vec<usize>| {
        let seg_low = 3_usize + seg_idx * SEGMENT_SPAN;
        let seg_high = seg_low.saturating_add(SEGMENT_SPAN - 2).min(max_odd);
        if seg_high < seg_low {
            return;
        }
        let odd_count = (seg_high - seg_low) / 2 + 1;
        let used_bytes = odd_count.div_ceil(8);
        let mut scratch = crate::kernels::reuse::global_pool().scratch_bytes(SEGMENT_BYTES);
        scratch.resize(SEGMENT_BYTES, 0);
        scratch[..used_bytes].fill(0);

        for &p in base_primes.iter().skip(1) {
            let p_sq = p.saturating_mul(p);
            if p_sq > seg_high {
                break;
            }
            let mut start = if p_sq >= seg_low {
                p_sq
            } else {
                let rem = seg_low % p;
                let first = if rem == 0 {
                    seg_low
                } else {
                    seg_low + (p - rem)
                };
                if (first & 1) == 0 {
                    first + p
                } else {
                    first
                }
            };
            if start < seg_low {
                start = seg_low;
            }
            let mut bit_idx = (start - seg_low) >> 1;
            while bit_idx < odd_count {
                scratch[bit_idx >> 3] |= 1_u8 << (bit_idx & 7);
                bit_idx += p;
            }
        }

        let full_bytes = odd_count >> 3;
        for (byte_idx, &byte) in scratch[..full_bytes].iter().enumerate() {
            let mut free = !byte;
            let base_num = seg_low + (byte_idx << 4);
            while free != 0 {
                let bit = free.trailing_zeros() as usize;
                free &= free - 1;
                out.push(base_num + (bit << 1));
            }
        }
        let rem_bits = odd_count & 7;
        if rem_bits > 0 {
            let mask = (1_u8 << rem_bits) - 1;
            let mut free = (!scratch[full_bytes]) & mask;
            let base_num = seg_low + (full_bytes << 4);
            while free != 0 {
                let bit = free.trailing_zeros() as usize;
                free &= free - 1;
                out.push(base_num + (bit << 1));
            }
        }
    };

    if crate::execution::should_split(num_segments, 4) {
        let chunks = crate::kernels::par::map_indexed(num_segments, 4, |seg_idx| {
            let mut local = Vec::new();
            sieve_segment(seg_idx, &mut local);
            local
        });
        for chunk in chunks {
            primes.extend(chunk);
        }
    } else {
        for seg_idx in 0..num_segments {
            sieve_segment(seg_idx, &mut primes);
        }
    }

    primes
}

/// A square root of `n` modulo the prime `p`, choosing the smaller root.
///
/// Accepts `p = 2`. Returns `None` for a composite or invalid modulus or
/// when `n` is a non-residue. Use [`PrimeSqrtPlan`] for repeated queries
/// with the same prime modulus.
#[must_use]
pub fn sqrt_mod_prime(
    n: u64,
    p: u64,
) -> Option<u64> {
    if p == 2 {
        return Some(n & 1);
    }
    let method = PrimeSqrtMethod::for_modulus(p)?;
    method.root(n, p, || method.correction(p))
}

/// Prepared modular square roots for a fixed machine-word prime.
///
/// The plan checks primality and prepares the root of unity used by its
/// square-root method. Repeated queries reuse these values.
#[derive(Clone, Copy, Debug)]
pub struct PrimeSqrtPlan {
    modulus: u64,
    method: PrimeSqrtMethod,
    correction: u64,
}

impl PrimeSqrtPlan {
    /// Prepares an odd prime modulus or `2`. Returns `None` for other moduli.
    #[must_use]
    pub fn new(modulus: u64) -> Option<Self> {
        let method = PrimeSqrtMethod::for_modulus(modulus)?;
        Some(Self {
            modulus,
            method,
            correction: method.correction(modulus)?,
        })
    }

    /// The prime modulus used by this plan.
    #[must_use]
    pub const fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Returns the smaller root of `n` modulo the prepared prime.
    ///
    /// The argument is reduced modulo the prime. Non-residues give `None`.
    #[must_use]
    #[inline]
    pub fn sqrt(
        &self,
        n: u64,
    ) -> Option<u64> {
        if self.modulus == 2 {
            return Some(n & 1);
        }
        self.method.root(n, self.modulus, || Some(self.correction))
    }
}

#[derive(Clone, Copy, Debug)]
enum PrimeSqrtMethod {
    Binary,
    ThreeModFour,
    FiveModEight,
    Tonelli { odd: u64, twos: u32 },
}

impl PrimeSqrtMethod {
    const fn for_modulus(p: u64) -> Option<Self> {
        if p == 2 {
            return Some(Self::Binary);
        }
        if p < 3 || !is_prime_miller_rabin(p) {
            return None;
        }
        match p % 8 {
            | 3 | 7 => Some(Self::ThreeModFour),
            | 5 => Some(Self::FiveModEight),
            | _ => {
                let twos = (p - 1).trailing_zeros();
                Some(Self::Tonelli {
                    odd: (p - 1) >> twos,
                    twos,
                })
            },
        }
    }

    fn correction(
        self,
        p: u64,
    ) -> Option<u64> {
        match self {
            | Self::Binary | Self::ThreeModFour => Some(1),
            | Self::FiveModEight => Some(mod_pow(2, p / 4, p)),
            | Self::Tonelli { odd, twos } => {
                let mont = MontU64::new(p)?;
                let minus_one_m = p - mont.one();
                (2..p).find_map(|z| {
                    let root_m = mont.pow_mont(mont.to_mont(z), odd);
                    (repeated_square_mont(mont, root_m, twos - 1) == minus_one_m)
                        .then(|| mont.from_mont(root_m))
                })
            },
        }
    }

    fn root(
        self,
        n: u64,
        p: u64,
        correction: impl FnOnce() -> Option<u64>,
    ) -> Option<u64> {
        let a = n % p;
        if a == 0 {
            return Some(0);
        }
        let root = match self {
            | Self::Binary => return Some(a),
            | Self::ThreeModFour => {
                let mont = MontU64::new(p)?;
                let a_m = mont.to_mont(a);
                let root_m = mont.pow_mont(a_m, p / 4 + 1);
                if mont.mul_mont(root_m, root_m) != a_m {
                    return None;
                }
                mont.from_mont(root_m)
            },
            | Self::FiveModEight => {
                let mont = MontU64::new(p)?;
                let a_m = mont.to_mont(a);
                let root_m = mont.pow_mont(a_m, p / 8 + 1);
                let squared_m = mont.mul_mont(root_m, root_m);
                if squared_m == a_m {
                    mont.from_mont(root_m)
                } else if squared_m == mont.sub(0, a_m) {
                    let corr_m = mont.to_mont(correction()?);
                    mont.from_mont(mont.mul_mont(root_m, corr_m))
                } else {
                    return None;
                }
            },
            | Self::Tonelli { odd, twos } => tonelli_root(a, p, odd, twos, correction)?,
        };
        Some(root.min(p - root))
    }
}

const fn repeated_square_mont(
    mont: MontU64,
    mut value_m: u64,
    times: u32,
) -> u64 {
    let mut i = 0;
    while i < times {
        value_m = mont.mul_mont(value_m, value_m);
        i += 1;
    }
    value_m
}

// Find the least i in 1..limit with value^(2^i) = 1 in Montgomery domain.
fn tonelli_order_mont(
    mont: MontU64,
    mut value_m: u64,
    limit: u32,
) -> Option<u32> {
    let one_m = mont.one();
    (1..limit).find(|_| {
        value_m = mont.mul_mont(value_m, value_m);
        value_m == one_m
    })
}

fn tonelli_root(
    a: u64,
    p: u64,
    odd: u64,
    mut limit: u32,
    correction: impl FnOnce() -> Option<u64>,
) -> Option<u64> {
    let mont = MontU64::new(p)?;
    let one_m = mont.one();
    let a_m = mont.to_mont(a);
    // One exponentiation supplies both a^((odd+1)/2) and a^odd.
    let power_m = mont.pow_mont(a_m, odd / 2);
    let mut root_m = mont.mul_mont(a_m, power_m);
    let mut residue_m = mont.mul_mont(root_m, power_m);
    if residue_m == one_m {
        return Some(mont.from_mont(root_m));
    }
    let mut order = tonelli_order_mont(mont, residue_m, limit)?;
    let mut unity_m = mont.to_mont(correction()?);
    loop {
        let factor_m = repeated_square_mont(mont, unity_m, limit - order - 1);
        unity_m = mont.mul_mont(factor_m, factor_m);
        residue_m = mont.mul_mont(residue_m, unity_m);
        root_m = mont.mul_mont(root_m, factor_m);
        if residue_m == one_m {
            return Some(mont.from_mont(root_m));
        }
        limit = order;
        order = tonelli_order_mont(mont, residue_m, limit)?;
    }
}

/// The least `x >= 0` with `g^x ≡ h (mod p)`.
///
/// Accepts prime and composite moduli, including bases that share factors with
/// the modulus. Returns `None` if there is no solution or `p` is outside
/// `2..=2^48`. After removing common factors, baby-step giant-step takes
/// `O(sqrt p)` time and memory in the worst case.
#[must_use]
pub fn discrete_log(
    g: u64,
    h: u64,
    p: u64,
) -> Option<u64> {
    if !(2..=1 << 48).contains(&p) {
        return None;
    }
    let g = g % p;
    let (mut target, mut modulus) = (h % p, p);
    let (mut removed, mut prefix) = (0, 1);
    loop {
        if target == prefix {
            return Some(removed);
        }
        let divisor = gcd(g, modulus);
        if divisor == 1 {
            break;
        }
        // For each further factor of g, divide the congruence by its gcd
        // with the current modulus. The prefix accounts for those factors.
        if !target.is_multiple_of(divisor) {
            return None;
        }
        target /= divisor;
        modulus /= divisor;
        prefix = mul_mod(prefix, g / divisor, modulus);
        removed += 1;
    }
    if gcd(target, modulus) != 1 {
        return None;
    }

    let inverse = mod_inverse(i64::try_from(prefix).ok()?, i64::try_from(modulus).ok()?)?;
    let target = mul_mod(target, u64::try_from(inverse).ok()?, modulus);
    discrete_log_coprime(g, target, modulus).map(|x| x + removed)
}

fn discrete_log_coprime(
    g: u64,
    h: u64,
    p: u64,
) -> Option<u64> {
    let m = isqrt(p) + 1;
    let mut table = std::collections::HashMap::new();
    let mut e = 1;
    for j in 0..m {
        if e == h {
            return Some(j);
        }
        if j > 0 && e == 1 {
            // The cycle containing 1 did not reach h.
            return None;
        }
        table.insert(e, j);
        e = mul_mod(e, g, p);
    }
    // e = g^m; stepping backwards by m scans exponents in ascending blocks.
    let inverse = mod_inverse(i64::try_from(e).ok()?, i64::try_from(p).ok()?)?;
    let factor = u64::try_from(inverse).ok()?;
    let mut gamma = h;
    for i in 1..=m {
        gamma = mul_mod(gamma, factor, p);
        if let Some(&j) = table.get(&gamma) {
            return Some(i * m + j);
        }
    }
    None
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    #[test]
    fn montgomery_arithmetic_round_trips_and_matches_u128_oracle() {
        assert!(MontU64::new(0).is_none());
        assert!(MontU64::new(1).is_none());
        assert!(MontU64::new(2).is_none());
        assert!(MontU64::new(100).is_none());

        for modulus in [
            3_u64,
            5,
            17,
            97,
            65_537,
            1_000_000_007,
            (1_u64 << 61) - 1,
            u64::MAX - 58,
            u64::MAX,
        ] {
            let mont = MontU64::new(modulus).expect("odd modulus > 1");
            assert_eq!(mont.modulus(), modulus);
            assert_eq!(modulus.wrapping_mul(mont.n_prime()), u64::MAX);
            assert_eq!(mont.from_mont(mont.one()), 1);
            assert_eq!(mont.redc(u128::from(mont.r2_mod())), mont.one());
            let wide = (1_u128 << 100) + 123_457;
            assert_eq!(
                mont.from_mont(mont.to_mont_u128(wide)),
                (wide % u128::from(modulus)) as u64
            );

            for a in [0_u64, 1, 2, 7, 42, modulus / 2, modulus - 2, modulus - 1] {
                let a_m = mont.to_mont(a);
                assert_eq!(mont.from_mont(a_m), a % modulus);
                let half_a = mont.half(a);
                assert_eq!(mul_mod(half_a, 2, modulus), a % modulus);
                for b in [0_u64, 1, 3, 19, modulus / 3, modulus - 1] {
                    let b_m = mont.to_mont(b);
                    let prod = mont.from_mont(mont.mul_mont(a_m, b_m));
                    assert_eq!(prod, mul_mod(a, b, modulus));
                    // `add`/`sub` take residues in `[0, N)`, so the loop
                    // values are reduced first; the oracle is the u128 one.
                    let (ra, rb) = (a % modulus, b % modulus);
                    assert_eq!(
                        mont.add(ra, rb),
                        ((u128::from(ra) + u128::from(rb)) % u128::from(modulus)) as u64
                    );
                    assert_eq!(
                        mont.sub(ra, rb),
                        ((u128::from(ra) + u128::from(modulus) - u128::from(rb))
                            % u128::from(modulus)) as u64
                    );
                }
                for exp in [0_u64, 1, 2, 13, 64, 1_000_003] {
                    let expected =
                        mod_pow_u128(u128::from(a), u128::from(exp), u128::from(modulus)) as u64;
                    assert_eq!(mont.pow(a, exp), expected);
                    assert_eq!(mod_pow(u128::from(a), exp, modulus), expected);
                }
            }
        }
    }

    #[test]
    fn binary_gcd_and_u128_primitives_match_invariants() {
        for &(a, b, expected) in &[
            (0_u64, 0_u64, 0_u64),
            (0, 17, 17),
            (17, 0, 17),
            (48, 18, 6),
            (270, 192, 6),
            (u64::MAX, u64::MAX - 1, 1),
            (1_u64 << 63, 1_u64 << 40, 1_u64 << 40),
        ] {
            assert_eq!(gcd(a, b), expected);
            assert_eq!(gcd_u128(u128::from(a), u128::from(b)), u128::from(expected));
        }
        let big_g = (1_u128 << 100) + 15;
        let a = big_g * 97;
        let b = big_g * 89;
        assert_eq!(gcd_u128(a, b), big_g);

        assert_eq!(mul_mod_u128(5, 7, 0), 0);
        assert_eq!(mul_mod_u128(5, 7, 1), 0);
        assert_eq!(mul_mod_u128(0, 7, 11), 0);
        assert_eq!(mod_pow_u128(5, 0, 1), 0);
        assert_eq!(mod_pow_u128(0, 0, 7), 1);

        let m = u128::MAX - 158;
        let x = m - 7;
        let y = m - 11;
        assert_eq!(mul_mod_u128(x, y, m), 77);
        // Fermat test on Mersenne prime 2^127 - 1 with > 64-bit exponent.
        let m127 = (1_u128 << 127) - 1;
        assert_eq!(mod_pow_u128(3, m127 - 1, m127), 1);
    }

    #[test]
    fn primality_tests_miller_rabin_bpsw_and_u128_agree() {
        for n in 0..2_000_u64 {
            let mr = is_prime_miller_rabin(n);
            let bpsw = is_prime_bpsw(n);
            let p128 = is_prime_u128(u128::from(n));
            assert_eq!(mr, bpsw, "mismatch at n = {n}");
            assert_eq!(mr, p128, "u128 mismatch at n = {n}");
        }
        // Carmichael numbers, prime squares, and strong pseudoprimes.
        for composite in [
            561_u64,
            1_105,
            1_681,
            1_729,
            2_047,
            1_373_653,
            25_326_001,
            3_215_031_751,
            2_152_302_898_747,
            3_474_749_660_383,
            341_550_071_728_321,
            3_825_123_056_546_413_051,
            u64::MAX,
        ] {
            assert!(!is_prime_miller_rabin(composite), "{composite}");
            assert!(!is_prime_bpsw(composite), "{composite}");
            assert!(!is_prime_u128(u128::from(composite)), "{composite}");
        }
        for prime in [
            2_u64,
            3,
            65_537,
            2_147_483_647,
            (1_u64 << 61) - 1,
            u64::MAX - 58,
        ] {
            assert!(is_prime_miller_rabin(prime), "{prime}");
            assert!(is_prime_bpsw(prime), "{prime}");
            assert!(is_prime_u128(u128::from(prime)), "{prime}");
        }
        // 12-base Miller-Rabin strong pseudoprime 318665857834031151167461 (~80 bits),
        // u128::MAX, and Mersenne prime 2^127 - 1.
        let sp12: u128 = 318_665_857_834_031_151_167_461;
        assert!(!is_prime_u128(sp12));
        assert!(!is_prime_u128(u128::MAX));
        let mersenne_127: u128 = (1_u128 << 127) - 1;
        assert!(is_prime_u128(mersenne_127));
    }

    #[test]
    fn jacobi_and_legendre_symbols_match_quadratic_residues() {
        for p in [3_u64, 5, 7, 11, 13, 17, 97] {
            for a in 0..p {
                let is_res = (1..p).any(|x| (x * x) % p == a);
                let expected = if a == 0 {
                    0
                } else if is_res {
                    1
                } else {
                    -1
                };
                assert_eq!(legendre_symbol(a, p), expected);
                assert_eq!(jacobi_symbol(i128::from(a), p), expected);
                assert_eq!(jacobi_symbol(i128::from(a) - i128::from(p), p), expected);
            }
        }
        assert_eq!(legendre_symbol(2, 15), 0);
        assert_eq!(jacobi_symbol(2, 15), 1);
        assert_eq!(jacobi_symbol(5, 0), 0);
        assert_eq!(jacobi_symbol(5, 4), 0);
        assert_eq!(jacobi_symbol(-1, 1), 1);
    }

    #[test]
    fn arithmetic_functions_and_primitive_root_hold() {
        assert_eq!(carmichael_lambda(0), 0);
        assert_eq!(carmichael_lambda(1), 1);
        assert_eq!(carmichael_lambda(2), 1);
        assert_eq!(carmichael_lambda(4), 2);
        assert_eq!(carmichael_lambda(8), 2);
        assert_eq!(carmichael_lambda(15), 4);
        assert_eq!(carmichael_lambda(16), 4);
        assert_eq!(carmichael_lambda(72), 6);
        assert_eq!(carmichael_lambda(561), 80);
        for n in 2..=80_u64 {
            let lam = carmichael_lambda(n);
            for a in 1..n {
                if gcd(a, n) == 1 {
                    assert_eq!(mod_pow(u128::from(a), lam, n), 1);
                }
            }
        }

        assert_eq!(divisors(0), Vec::<u64>::new());
        assert_eq!(divisors(1), vec![1]);
        assert_eq!(divisors(28), vec![1, 2, 4, 7, 14, 28]);
        assert_eq!(divisor_count(0), 0);
        assert_eq!(divisor_count(28), 6);
        assert_eq!(divisor_sigma(0, 1), 0);
        assert_eq!(divisor_sigma(1, 3), 1);
        assert_eq!(divisor_sigma(28, 0), 6);
        assert_eq!(divisor_sigma(28, 1), 56);
        assert_eq!(divisor_sigma(6, 2), 1 + 4 + 9 + 36);
        // sigma_4(p) = 1 + p + p^2 + p^3 + p^4, exact in u128 for p = 1e9+7.
        assert_eq!(
            divisor_sigma(1_000_000_007, 4),
            1_000_000_028_000_000_294_000_001_372_000_002_402
        );

        assert_eq!(primitive_root(0), None);
        assert_eq!(primitive_root(1), None);
        assert_eq!(primitive_root(2), Some(1));
        assert_eq!(primitive_root(3), Some(2));
        assert_eq!(primitive_root(4), Some(3));
        assert_eq!(primitive_root(5), Some(2));
        assert_eq!(primitive_root(7), Some(3));
        assert_eq!(primitive_root(8), None);
        assert_eq!(primitive_root(9), Some(2));
        assert_eq!(primitive_root(12), None);
        assert_eq!(primitive_root(14), Some(3));
        assert_eq!(primitive_root(15), None);
        assert_eq!(primitive_root(18), Some(5));
        assert_eq!(primitive_root(25), Some(2));
        assert_eq!(primitive_root(31), Some(3));
        assert_eq!(primitive_root(50), Some(3));
        assert_eq!(primitive_root(65_537), Some(3));
    }

    #[test]
    fn primes_sieve_uses_cache_and_segmented_sieve_above_one_million() {
        assert!(primes_sieve(0).is_empty());
        assert!(primes_sieve(1).is_empty());
        assert_eq!(primes_sieve(2), vec![2]);
        let small = primes_sieve(100);
        assert_eq!(small.len(), 25);
        assert_eq!(small.first(), Some(&2));
        assert_eq!(small.last(), Some(&97));

        let one_m = primes_sieve(1_000_000);
        assert_eq!(one_m.len(), 78_498);
        assert_eq!(one_m.last(), Some(&999_983));

        // Prime immediately above the cache cutoff tests exact boundary & tail bits.
        let just_above = primes_sieve(1_000_003);
        assert_eq!(just_above.len(), 78_499);
        assert_eq!(just_above.last(), Some(&1_000_003));

        // π(1_500_000) == 114_155
        let large = primes_sieve(1_500_000);
        assert_eq!(large.len(), 114_155);
        assert_eq!(large.first(), Some(&2));
        assert_eq!(large.last(), Some(&1_499_977));
    }
}
