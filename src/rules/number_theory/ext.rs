//! Further number theory: quadratic residues and modular square roots,
//! orders and primitive roots, discrete logarithms, arithmetic functions,
//! sums of two squares, linear congruences and Farey sequences.
//!
//! | operator | value |
//! |---|---|
//! | `legendre(a, p)` | the Legendre symbol `(a/p)` for an odd prime `p` (`-1`, `0` or `1`) |
//! | `kronecker(a, n)` | the Kronecker symbol, extending the Jacobi symbol to even and negative `n` |
//! | `sqrtmod(a, p)` | the smaller square root of `a` modulo the prime `p` (Tonelli–Shanks); stays unreduced when `a` is a non-residue |
//! | `quadratic_residues(n)` | the sorted non-zero squares modulo `n` (`n <= 100000`) |
//! | `mult_order(a, n)` | the multiplicative order of `a` modulo `n` (`gcd(a, n) = 1`) |
//! | `carmichael(n)` | the Carmichael function `λ(n)`, the exponent of the group of units |
//! | `primitive_root(n)` | the smallest primitive root modulo `n` (only for `n = 2, 4, p^k, 2 p^k`) |
//! | `dlog(g, h, n)` | the least `x >= 0` with `g^x ≡ h (mod n)`, by baby-step giant-step; stays unreduced when none exists |
//! | `mobius(n)`, `liouville(n)` | the Möbius and Liouville functions |
//! | `omega(n)`, `bigomega(n)` | the number of distinct prime factors, and with multiplicity |
//! | `radical(n)` | the product of the distinct prime factors |
//! | `sigma(n, k)` | the sum of the `k`-th powers of the divisors of `n` |
//! | `divisors(n)` | the sorted positive divisors |
//! | `is_squarefree(n)`, `is_square(n)` | truth values |
//! | `isqrt(n)`, `iroot(n, k)` | the floor of the square and `k`-th root |
//! | `perfect_power(n)` | `list(b, e)` with `n = b^e` and `e` maximal (`list(n, 1)` when `n` is not a perfect power) |
//! | `prevprime(n)`, `prime_pi(n)`, `nth_prime(k)` | the largest prime below `n`, the number of primes up to `n` (`n <= 100000`), the `k`-th prime |
//! | `two_squares(n)` | `list(a, b)`, `0 <= a <= b`, with `a² + b² = n` (one representation; stays unreduced when none exists) |
//! | `lincong(a, b, m)` | the solutions of `a x ≡ b (mod m)` as `list(x0, step)`, or `list()` when there is none |
//! | `binomial_mod(n, k, p)` | `C(n, k) mod p` for a prime `p`, by Lucas' theorem |
//! | `farey(n)` | the Farey sequence of order `n` as `list(list(num, den), ...)` |

use num_bigint::BigInt;
use num_integer::Roots;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::egcd;
use super::exact;
use super::factor;
use super::gcd;
use super::invmod;
use super::is_prime;
use super::jacobi;
use super::modulo;
use super::Value;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;

const SET: &str = "number_theory";
const MAX_TABLE: u64 = 100_000;
const MAX_BSGS: u64 = 1 << 22;

fn big(v: u64) -> BigInt {
    BigInt::from(v)
}

fn small(n: &BigInt) -> Option<u64> {
    n.to_u64().filter(|&v| v >= 1)
}

fn int(v: impl Into<BigInt>) -> Option<Value> {
    Some(Value::Int(v.into()))
}

fn count(n: usize) -> Option<Value> {
    u64::try_from(n).ok().and_then(int)
}

/// `a^e mod m` on machine words.
fn pow_mod(
    a: u64,
    e: u64,
    m: u64,
) -> u64 {
    big(a).modpow(&big(e), &big(m)).to_u64().unwrap_or(0)
}

fn mulm(
    a: u64,
    b: u64,
    m: u64,
) -> u64 {
    u64::try_from(u128::from(a) * u128::from(b) % u128::from(m)).unwrap_or(0)
}

/// A square root of `a` modulo the prime `p` by Tonelli–Shanks.
fn tonelli_shanks(
    a: u64,
    p: u64,
) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if p == 2 {
        return Some(a);
    }
    if pow_mod(a, (p - 1) / 2, p) != 1 {
        return None;
    }
    let (mut q, mut s) = (p - 1, 0_u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while pow_mod(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = pow_mod(z, q, p);
    let mut t = pow_mod(a, q, p);
    let mut r = pow_mod(a, q.div_ceil(2), p);
    while t != 1 {
        let mut i = 0_u32;
        let mut probe = t;
        while probe != 1 {
            probe = pow_mod(probe, 2, p);
            i += 1;
            if i >= m {
                return None;
            }
        }
        let b = pow_mod(c, 1 << (m - i - 1), p);
        m = i;
        c = pow_mod(b, 2, p);
        t = mulm(t, c, p);
        r = mulm(r, b, p);
    }
    Some(r.min(p - r))
}

/// Carmichael's `λ(n)`.
fn carmichael(n: &BigInt) -> Option<BigInt> {
    let mut lambda = BigInt::one();
    for (p, e) in factor(n)? {
        let part = if p == 2 {
            match e {
                | 1 => BigInt::one(),
                | 2 => big(2),
                | _ => big(1) << (e - 2),
            }
        } else {
            big(p - 1) * big(p).pow(e - 1)
        };
        lambda = &lambda / gcd(&lambda, &part) * part;
    }
    Some(lambda)
}

fn mult_order(
    a: &BigInt,
    n: &BigInt,
) -> Option<BigInt> {
    if !n.is_positive() {
        return None;
    }
    if n.is_one() {
        return Some(BigInt::one());
    }
    let a = modulo(a, n)?;
    if !gcd(&a, n).is_one() {
        return None;
    }
    let mut order = carmichael(n)?;
    let primes: Vec<u64> = factor(&order)?.into_iter().map(|(p, _)| p).collect();
    for p in primes {
        let p = big(p);
        while (&order % &p).is_zero() && a.modpow(&(&order / &p), n).is_one() {
            order /= &p;
        }
    }
    Some(order)
}

fn primitive_root(n: &BigInt) -> Option<BigInt> {
    let modulus = small(n)?;
    match modulus {
        | 1 => return Some(BigInt::zero()),
        | 2 => return Some(BigInt::one()),
        | 4 => return Some(big(3)),
        | _ => {},
    }
    let primes = factor(n)?;
    let cyclic = match primes.as_slice() {
        | [(p, _)] => *p != 2,
        | [(2, 1), _] => true,
        | _ => false,
    };
    if !cyclic {
        return None;
    }
    let phi: u64 = primes.iter().map(|&(p, e)| (p - 1) * p.pow(e - 1)).product();
    let prime_factors: Vec<u64> = factor(&big(phi))?.into_iter().map(|(p, _)| p).collect();
    (2..modulus)
        .find(|&g| {
            gcd(&big(g), n).is_one() && prime_factors.iter().all(|&q| pow_mod(g, phi / q, modulus) != 1)
        })
        .map(BigInt::from)
}

/// Least `x >= 0` with `g^x = h (mod n)`.
fn dlog(
    g: &BigInt,
    h: &BigInt,
    n: &BigInt,
) -> Option<BigInt> {
    let modulus = small(n)?;
    if modulus == 1 {
        return Some(BigInt::zero());
    }
    let (g, h) = (modulo(g, n)?, modulo(h, n)?);
    // Without invertibility, scan (the sequence is eventually periodic).
    if !gcd(&g, n).is_one() {
        let mut value = BigInt::one();
        for x in 0..modulus.min(MAX_BSGS) {
            if value == h {
                return Some(big(x));
            }
            value = &value * &g % n;
        }
        return None;
    }
    let order = mult_order(&g, n)?.to_u64()?;
    let m = order.isqrt() + 1;
    if m > MAX_BSGS {
        return None;
    }
    let mut table = std::collections::HashMap::with_capacity(usize::try_from(m).ok()?);
    let mut value = BigInt::one();
    for j in 0..m {
        table.entry(value.clone()).or_insert(j);
        value = &value * &g % n;
    }
    let giant_inverse = invmod(&g.modpow(&big(m), n), n)?;
    let mut gamma = h;
    for i in 0..=m {
        if let Some(&j) = table.get(&gamma) {
            return Some(big(i * m + j));
        }
        gamma = gamma * &giant_inverse % n;
    }
    None
}

/// Divisors, ascending.
fn divisors(n: &BigInt) -> Option<Vec<BigInt>> {
    let mut out = vec![BigInt::one()];
    for (p, e) in factor(n)? {
        let current = out.clone();
        let mut power = BigInt::one();
        for _ in 1..=e {
            power *= p;
            out.extend(current.iter().map(|d| d * &power));
        }
    }
    out.sort();
    Some(out)
}

/// Gaussian multiplication `(a + bi)(c + di)`.
fn gauss_mul(
    x: &(BigInt, BigInt),
    y: &(BigInt, BigInt),
) -> (BigInt, BigInt) {
    (&x.0 * &y.0 - &x.1 * &y.1, &x.0 * &y.1 + &x.1 * &y.0)
}

fn two_squares(n: &BigInt) -> Option<(BigInt, BigInt)> {
    if n.is_zero() {
        return Some((BigInt::zero(), BigInt::zero()));
    }
    let mut acc = (BigInt::one(), BigInt::zero());
    for (p, e) in factor(n)? {
        match p % 4 {
            | 3 => {
                if e % 2 == 1 {
                    return None;
                }
                let scale = big(p).pow(e / 2);
                acc = (&acc.0 * &scale, &acc.1 * &scale);
            },
            | 2 => {
                for _ in 0..e {
                    acc = gauss_mul(&acc, &(BigInt::one(), BigInt::one()));
                }
            },
            | _ => {
                let root = tonelli_shanks(p - 1, p)?;
                // Euclid on (p, root) until the remainder drops below sqrt p.
                let (mut a, mut b) = (p, root);
                let bound = p.isqrt();
                while b > bound {
                    (a, b) = (b, a % b);
                }
                let rest = p - b * b;
                let c = rest.isqrt();
                if c * c != rest {
                    return None;
                }
                let prime_part = (big(b), big(c));
                for _ in 0..e {
                    acc = gauss_mul(&acc, &prime_part);
                }
            },
        }
    }
    let (x, y) = (acc.0.abs(), acc.1.abs());
    Some(if x <= y { (x, y) } else { (y, x) })
}

/// The Kronecker symbol `(a/n)`.
fn kronecker(
    a: &BigInt,
    n: &BigInt,
) -> Option<i8> {
    if n.is_zero() {
        return Some(i8::from(a.abs().is_one()));
    }
    let mut result = 1_i8;
    let mut n = n.clone();
    if n.is_negative() {
        n = -n;
        if a.is_negative() {
            result = -result;
        }
    }
    let twos = n.trailing_zeros().unwrap_or(0);
    n >>= twos;
    if twos > 0 {
        if !a.bit(0) {
            return Some(0);
        }
        let low = modulo(a, &big(8))?.to_u8()?;
        if twos % 2 == 1 && matches!(low, 3 | 5) {
            result = -result;
        }
    }
    Some(result * jacobi(a, &n)?)
}

fn lucas_binomial(
    mut n: u64,
    mut k: u64,
    p: u64,
) -> u64 {
    let mut result = 1_u64;
    while n > 0 || k > 0 {
        let (ni, ki) = (n % p, k % p);
        if ki > ni {
            return 0;
        }
        // C(ni, ki) mod p with ni < p, as a product of small factors.
        let mut term = 1_u64;
        for j in 0..ki {
            term = mulm(mulm(term, (ni - j) % p, p), pow_mod((j + 1) % p, p - 2, p), p);
        }
        result = mulm(result, term, p);
        n /= p;
        k /= p;
    }
    result
}

fn farey(order: u64) -> Vec<Vec<BigInt>> {
    let (mut a, mut b, mut c, mut d) = (0_u64, 1_u64, 1_u64, order);
    let mut out = vec![vec![big(a), big(b)]];
    while c <= order {
        let k = (order + b) / d;
        (a, b, c, d) = (c, d, k * c - a, k * d - b);
        out.push(vec![big(a), big(b)]);
    }
    out
}

#[allow(clippy::too_many_lines)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));
    let ternary = |name: &str| OpDescriptor::new(name, Arity::Fixed(3));

    exact(i, SET, binary("legendre"), |a| match a {
        | [x, p] if *p > big(2) && is_prime(p) => jacobi(x, p).and_then(int),
        | _ => None,
    })?;
    exact(i, SET, binary("kronecker"), |a| match a {
        | [x, n] => kronecker(x, n).and_then(int),
        | _ => None,
    })?;
    exact(i, SET, binary("sqrtmod"), |a| match a {
        | [x, p] if is_prime(p) => {
            let p = p.to_u64()?;
            let x = modulo(x, &big(p))?.to_u64()?;
            tonelli_shanks(x, p).and_then(int)
        },
        | _ => None,
    })?;
    exact(i, SET, unary("quadratic_residues"), |a| {
        let [n] = a else { return None };
        let n = small(n).filter(|&v| v <= MAX_TABLE)?;
        let mut seen: Vec<u64> = (1..n).map(|x| x * x % n).filter(|&r| r != 0).collect();
        seen.sort_unstable();
        seen.dedup();
        Some(Value::Flat(seen.into_iter().map(BigInt::from).collect()))
    })?;
    exact(i, SET, binary("mult_order"), |a| match a {
        | [x, n] => mult_order(x, n).map(Value::Int),
        | _ => None,
    })?;
    exact(i, SET, unary("carmichael"), |a| match a {
        | [n] => carmichael(n).map(Value::Int),
        | _ => None,
    })?;
    exact(i, SET, unary("primitive_root"), |a| match a {
        | [n] => primitive_root(n).map(Value::Int),
        | _ => None,
    })?;
    exact(i, SET, ternary("dlog"), |a| match a {
        | [g, h, n] => dlog(g, h, n).map(Value::Int),
        | _ => None,
    })?;
    exact(i, SET, unary("mobius"), |a| {
        let [n] = a else { return None };
        let primes = factor(n)?;
        if primes.iter().any(|&(_, e)| e > 1) {
            return int(0);
        }
        int(if primes.len() % 2 == 0 { 1 } else { -1 })
    })?;
    exact(i, SET, unary("liouville"), |a| {
        let [n] = a else { return None };
        let total: u32 = factor(n)?.iter().map(|&(_, e)| e).sum();
        int(if total % 2 == 0 { 1 } else { -1 })
    })?;
    exact(i, SET, unary("omega"), |a| {
        let [n] = a else { return None };
        count(factor(n)?.len())
    })?;
    exact(i, SET, unary("bigomega"), |a| {
        let [n] = a else { return None };
        int(factor(n)?.iter().map(|&(_, e)| u64::from(e)).sum::<u64>())
    })?;
    exact(i, SET, unary("radical"), |a| {
        let [n] = a else { return None };
        int(factor(n)?.iter().map(|&(p, _)| big(p)).product::<BigInt>())
    })?;
    exact(i, SET, binary("sigma"), |a| {
        let [n, k] = a else { return None };
        let k = k.to_u32().filter(|&k| k <= 64)?;
        let mut total = BigInt::one();
        for (p, e) in factor(n)? {
            let q = big(p).pow(k);
            total *= (0..=e).fold(BigInt::zero(), |acc, j| acc + q.pow(j));
        }
        Some(Value::Int(total))
    })?;
    exact(i, SET, unary("divisors"), |a| {
        let [n] = a else { return None };
        divisors(n).map(Value::Flat)
    })?;
    exact(i, SET, unary("is_squarefree"), |a| {
        let [n] = a else { return None };
        Some(Value::Bool(factor(n)?.iter().all(|&(_, e)| e == 1)))
    })?;
    exact(i, SET, unary("is_square"), |a| {
        let [n] = a else { return None };
        if n.is_negative() {
            return Some(Value::Bool(false));
        }
        let r = n.sqrt();
        Some(Value::Bool(&r * &r == *n))
    })?;
    exact(i, SET, unary("isqrt"), |a| match a {
        | [n] if !n.is_negative() => int(n.sqrt()),
        | _ => None,
    })?;
    exact(i, SET, binary("iroot"), |a| match a {
        | [n, k] if !n.is_negative() && k.is_positive() => int(n.nth_root(k.to_u32().filter(|&k| k <= 4096)?)),
        | _ => None,
    })?;
    exact(i, SET, unary("perfect_power"), |a| {
        let [n] = a else { return None };
        if *n < big(2) {
            return None;
        }
        let bits = u32::try_from(n.bits()).ok()?;
        for e in (2..=bits).rev() {
            let r = n.nth_root(e);
            if r.pow(e) == *n {
                return Some(Value::Flat(vec![r, BigInt::from(e)]));
            }
        }
        Some(Value::Flat(vec![n.clone(), BigInt::one()]))
    })?;
    exact(i, SET, unary("prevprime"), |a| {
        let [n] = a else { return None };
        let mut candidate = n - 1;
        while candidate >= big(2) {
            if is_prime(&candidate) {
                return Some(Value::Int(candidate));
            }
            candidate -= 1;
        }
        None
    })?;
    exact(i, SET, unary("prime_pi"), |a| {
        let [n] = a else { return None };
        let limit = usize::try_from(n.to_u64().filter(|&v| v <= MAX_TABLE)?).ok()?;
        count(crate::kernels::number_theory::primes_sieve(limit).len())
    })?;
    exact(i, SET, unary("nth_prime"), |a| {
        let [k] = a else { return None };
        let k = usize::try_from(k.to_u64().filter(|&v| v >= 1)?).ok()?;
        let mut limit = 100;
        while limit <= 50_000_000 {
            let primes = crate::kernels::number_theory::primes_sieve(limit);
            if let Some(&p) = primes.get(k - 1) {
                return count(p);
            }
            limit *= 4;
        }
        None
    })?;
    exact(i, SET, unary("two_squares"), |a| {
        let [n] = a else { return None };
        if n.is_negative() {
            return None;
        }
        two_squares(n).map(|(x, y)| Value::Flat(vec![x, y]))
    })?;
    exact(i, SET, ternary("lincong"), |a| {
        let [x, b, m] = a else { return None };
        if !m.is_positive() {
            return None;
        }
        let (g, s, _) = egcd(&modulo(x, m)?, m);
        let b = modulo(b, m)?;
        if !(&b % &g).is_zero() {
            return Some(Value::Flat(Vec::new()));
        }
        let step = m / &g;
        let x0 = modulo(&(s * (b / &g)), &step)?;
        Some(Value::Flat(vec![x0, step]))
    })?;
    exact(i, SET, ternary("binomial_mod"), |a| {
        let [n, k, p] = a else { return None };
        if !is_prime(p) || n.is_negative() || k.is_negative() {
            return None;
        }
        int(lucas_binomial(n.to_u64()?, k.to_u64()?, p.to_u64()?))
    })?;
    exact(i, SET, unary("farey"), |a| {
        let [n] = a else { return None };
        let n = n.to_u64().filter(|&v| (1..=200).contains(&v))?;
        Some(Value::Nested(farey(n)))
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::rules::number_theory::number_theory;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[number_theory()], src)
    }

    #[test]
    fn legendre_kronecker_and_sqrtmod() {
        assert_eq!(s("legendre(2, 7)"), "1");
        assert_eq!(s("legendre(3, 7)"), "-1");
        assert_eq!(s("legendre(14, 7)"), "0");
        assert_eq!(s("kronecker(2, 15)"), "1");
        assert_eq!(s("kronecker(3, 8)"), "-1");
        assert_eq!(s("kronecker(-1, 4)"), "1");
        assert_eq!(s("sqrtmod(10, 13)"), "6");
        assert_eq!(s("sqrtmod(13, 17)"), "8");
        assert_eq!(s("mod(sqrtmod(2, 1000000007)^2, 1000000007)"), "2");
        assert_eq!(s("sqrtmod(3, 7)"), "sqrtmod(3, 7)");
        assert_eq!(s("quadratic_residues(11)"), "list(1, 3, 4, 5, 9)");
    }

    #[test]
    fn orders_and_primitive_roots() {
        assert_eq!(s("mult_order(2, 7)"), "3");
        assert_eq!(s("mult_order(3, 7)"), "6");
        assert_eq!(s("mult_order(2, 15)"), "4");
        assert_eq!(s("carmichael(561)"), "80");
        assert_eq!(s("carmichael(8)"), "2");
        assert_eq!(s("primitive_root(7)"), "3");
        assert_eq!(s("primitive_root(25)"), "2");
        assert_eq!(s("primitive_root(14)"), "3");
        assert_eq!(s("primitive_root(8)"), "primitive_root(8)");
    }

    #[test]
    fn discrete_logarithm() {
        assert_eq!(s("dlog(3, 13, 17)"), "4");
        assert_eq!(s("powmod(2, dlog(2, 1234, 10007), 10007)"), "1234");
        assert_eq!(s("dlog(2, 3, 7)"), "dlog(2, 3, 7)");
        assert_eq!(s("dlog(4, 8, 12)"), "dlog(4, 8, 12)");
    }

    #[test]
    fn arithmetic_functions() {
        assert_eq!(s("mobius(30)"), "-1");
        assert_eq!(s("mobius(12)"), "0");
        assert_eq!(s("liouville(12)"), "-1");
        assert_eq!(s("omega(360)"), "3");
        assert_eq!(s("bigomega(360)"), "6");
        assert_eq!(s("radical(360)"), "30");
        assert_eq!(s("sigma(12, 2)"), "210");
        assert_eq!(s("divisors(12)"), "list(1, 2, 3, 4, 6, 12)");
        assert_eq!(s("is_squarefree(30)"), "true");
        assert_eq!(s("is_square(49)"), "true");
        assert_eq!(s("isqrt(99)"), "9");
        assert_eq!(s("iroot(1000, 3)"), "10");
        assert_eq!(s("perfect_power(64)"), "list(2, 6)");
        assert_eq!(s("perfect_power(10)"), "list(10, 1)");
        assert_eq!(s("prevprime(100)"), "97");
        assert_eq!(s("prime_pi(100)"), "25");
        assert_eq!(s("nth_prime(10)"), "29");
    }

    #[test]
    fn squares_congruences_and_farey() {
        assert_eq!(s("two_squares(5)"), "list(1, 2)");
        assert_eq!(s("two_squares(10)"), "list(1, 3)");
        assert_eq!(s("two_squares(13)"), "list(2, 3)");
        assert_eq!(s("two_squares(3)"), "two_squares(3)");
        assert_eq!(s("lincong(6, 4, 10)"), "list(4, 5)");
        assert_eq!(s("lincong(2, 1, 4)"), "list()");
        assert_eq!(s("binomial_mod(10, 3, 7)"), "1");
        assert_eq!(s("binomial_mod(100, 50, 7)"), "4");
        assert_eq!(s("farey(3)"), "list(list(0, 1), list(1, 3), list(1, 2), list(2, 3), list(1, 1))");
    }
}
