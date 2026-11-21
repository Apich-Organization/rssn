//! # Numerical Number Theory
//!
//! Machine-word integer number theory: greatest common divisors, modular
//! arithmetic, deterministic primality testing, factorisation (trial
//! division plus Pollard–Brent rho), square roots modulo a prime
//! (Tonelli–Shanks) and discrete logarithms for prime or composite moduli.
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
pub(crate) const fn mul_mod(
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
    if modulus <= 1 {
        return 0;
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

/// Performs a Miller-Rabin primality test.
///
/// The first twelve prime bases are a deterministic witness set over the
/// entire `u64` range. Divisibility by those primes is checked first.
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
                (2..p).find_map(|z| {
                    let root = mod_pow(u128::from(z), odd, p);
                    (repeated_square(root, twos - 1, p) == p - 1).then_some(root)
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
                let root = mod_pow(u128::from(a), p / 4 + 1, p);
                if mul_mod(root, root, p) != a {
                    return None;
                }
                root
            },
            | Self::FiveModEight => {
                let root = mod_pow(u128::from(a), p / 8 + 1, p);
                let squared = mul_mod(root, root, p);
                if squared == a {
                    root
                } else if squared == p - a {
                    mul_mod(root, correction()?, p)
                } else {
                    return None;
                }
            },
            | Self::Tonelli { odd, twos } => tonelli_root(a, p, odd, twos, correction)?,
        };
        Some(root.min(p - root))
    }
}

fn repeated_square(
    mut value: u64,
    times: u32,
    p: u64,
) -> u64 {
    for _ in 0..times {
        value = mul_mod(value, value, p);
    }
    value
}

// Find the least i in 1..limit with value^(2^i) = 1. In the
// first Tonelli step a non-residue has order 2^limit.
fn tonelli_order(
    mut value: u64,
    limit: u32,
    p: u64,
) -> Option<u32> {
    (1..limit).find(|_| {
        value = mul_mod(value, value, p);
        value == 1
    })
}

fn tonelli_root(
    a: u64,
    p: u64,
    odd: u64,
    mut limit: u32,
    correction: impl FnOnce() -> Option<u64>,
) -> Option<u64> {
    // One exponentiation supplies both a^((odd+1)/2) and a^odd.
    let power = mod_pow(u128::from(a), odd / 2, p);
    let mut root = mul_mod(a, power, p);
    let mut residue = mul_mod(root, power, p);
    if residue == 1 {
        return Some(root);
    }
    let mut order = tonelli_order(residue, limit, p)?;
    let mut unity = correction()?;
    loop {
        let factor = repeated_square(unity, limit - order - 1, p);
        unity = mul_mod(factor, factor, p);
        residue = mul_mod(residue, unity, p);
        root = mul_mod(root, factor, p);
        if residue == 1 {
            return Some(root);
        }
        limit = order;
        order = tonelli_order(residue, limit, p)?;
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
