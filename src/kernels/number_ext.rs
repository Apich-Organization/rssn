//! Extended number theory and exact combinatorics kernels.

use num_bigint::BigInt;
use num_bigint::BigUint;

use crate::kernels::number_theory;

/// Error returned by extended number-theoretic operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NumberExtError {
    /// Input slices were empty.
    Empty,
    /// Input slice lengths did not match.
    ShapeMismatch(String),
    /// Invalid input or a result exceeding the addressable capacity.
    Domain(String),
}

impl std::fmt::Display for NumberExtError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("empty input"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
        }
    }
}

impl std::error::Error for NumberExtError {}

const fn xgcd_i128(
    a: i128,
    b: i128,
) -> (i128, i128, i128) {
    let mut old_r = a;
    let mut r = b;
    let mut old_s = 1i128;
    let mut s = 0i128;
    let mut old_t = 0i128;
    let mut t = 1i128;
    while r != 0 {
        let q = old_r / r;
        let next_r = old_r - q * r;
        old_r = r;
        r = next_r;
        let next_s = old_s - q * s;
        old_s = s;
        s = next_s;
        let next_t = old_t - q * t;
        old_t = t;
        t = next_t;
    }
    if old_r < 0 {
        (-old_r, -old_s, -old_t)
    } else {
        (old_r, old_s, old_t)
    }
}

/// Extended greatest common divisor `(g, x, y)` such that `a * x + b * y == g` and `g >= 0`.
/// Wide outputs represent `gcd(i64::MIN, 0)` without overflow.
#[must_use]
pub const fn xgcd(
    a: i64,
    b: i64,
) -> (i128, i128, i128) {
    xgcd_i128(a as i128, b as i128)
}

/// Modular exponentiation: `base^exponent mod modulus`, in `[0, modulus)`.
///
/// # Errors
/// The modulus must be positive.
pub fn mod_pow(
    base: i64,
    exponent: u64,
    modulus: i64,
) -> Result<i64, NumberExtError> {
    if modulus <= 0 {
        return Err(NumberExtError::Domain(
            "mod_pow requires a positive modulus".to_string(),
        ));
    }
    if modulus == 1 {
        return Ok(0);
    }
    let base = base.rem_euclid(modulus) as u128;
    Ok(number_theory::mod_pow(base, exponent, modulus as u64) as i64)
}

/// Chinese Remainder Theorem: solves `x ≡ remainders[i] (mod moduli[i])`.
///
/// Returns `(first, period)` with `0 <= first < period`. The period is the
/// least common multiple of the moduli. Non-coprime moduli are supported.
///
/// # Errors
/// The inputs must be nonempty, of the same length, and have positive moduli
/// and consistent congruences. The period must fit in `i64`.
pub fn crt(
    remainders: &[i64],
    moduli: &[i64],
) -> Result<(i64, i64), NumberExtError> {
    if remainders.len() != moduli.len() {
        return Err(NumberExtError::ShapeMismatch(format!(
            "crt needs equal-length slices, found {} and {}",
            remainders.len(),
            moduli.len()
        )));
    }
    if remainders.is_empty() {
        return Err(NumberExtError::Empty);
    }
    if moduli[0] <= 0 {
        return Err(NumberExtError::Domain(
            "crt requires positive moduli".to_string(),
        ));
    }
    let mut m1 = i128::from(moduli[0]);
    let mut r1 = i128::from(remainders[0]).rem_euclid(m1);

    for (&rem, &mod_i) in remainders.iter().zip(moduli.iter()).skip(1) {
        if mod_i <= 0 {
            return Err(NumberExtError::Domain(
                "crt requires positive moduli".to_string(),
            ));
        }
        let m2 = i128::from(mod_i);
        let r2 = i128::from(rem).rem_euclid(m2);
        let (g, s, _) = xgcd_i128(m1, m2);
        let diff = r2 - r1;
        if diff % g != 0 {
            return Err(NumberExtError::Domain(
                "inconsistent congruences in crt".to_string(),
            ));
        }
        let m2_g = m2 / g;
        let lcm = m1
            .checked_mul(m2_g)
            .ok_or_else(|| NumberExtError::Domain("crt modulus overflow".to_string()))?;
        if lcm > i128::from(i64::MAX) {
            return Err(NumberExtError::Domain("crt modulus overflow".to_string()));
        }
        let step = ((diff / g).rem_euclid(m2_g) * s.rem_euclid(m2_g)).rem_euclid(m2_g);
        r1 = (r1 + m1 * step).rem_euclid(lcm);
        m1 = lcm;
    }

    let final_m = i64::try_from(m1)
        .map_err(|_| NumberExtError::Domain("crt modulus overflow".to_string()))?;
    Ok((r1 as i64, final_m))
}

/// Solves `a * x ≡ b (mod modulus)` over the integers.
///
/// A solution `(first, period)` represents exactly the integers
/// `first + k * period` for any integer `k`, with `0 <= first < period`.
/// The period is `modulus / gcd(a, modulus)`. Returns `None` when the
/// congruence has no solution. Modulus one admits all integers.
/// For example, `6 * x ≡ 8 (mod 14)` gives `(6, 7)`.
///
/// # Errors
/// The modulus must be nonzero.
pub fn solve_linear_congruence(
    a: i64,
    b: i64,
    modulus: u64,
) -> Result<Option<(u64, u64)>, NumberExtError> {
    if modulus == 0 {
        return Err(NumberExtError::Domain(
            "linear congruence requires a nonzero modulus".to_owned(),
        ));
    }
    let modulus = i128::from(modulus);
    let (gcd, coefficient, _) = xgcd_i128(i128::from(a), modulus);
    let b = i128::from(b);
    if b % gcd != 0 {
        return Ok(None);
    }
    let period = modulus / gcd;
    // The coefficient is bounded by modulus and |b| by 2^63, so the
    // product fits in i128 even when the modulus uses all 64 bits.
    let first = (coefficient * (b / gcd)).rem_euclid(period);
    Ok(Some((first as u64, period as u64)))
}

/// Möbius function `μ(n)`, extended to return zero for `n <= 0`.
#[must_use]
pub fn mobius(n: i64) -> i64 {
    if n <= 0 {
        return 0;
    }
    if n == 1 {
        return 1;
    }
    let factors = number_theory::factorize(n as u64);
    if factors.windows(2).any(|w| w[0] == w[1]) {
        return 0;
    }
    if factors.len().is_multiple_of(2) {
        1
    } else {
        -1
    }
}

/// Computes `μ(0)` through `μ(limit)` with a linear sieve.
///
/// The value at index zero is zero. Each composite is visited once using
/// its smallest prime factor, so consecutive values need no separate
/// factorizations. The sieve uses `O(limit)` time and storage.
///
/// # Errors
/// The inclusive length must fit the address space and its storage must be
/// allocatable.
pub fn mobius_sieve(limit: usize) -> Result<Vec<i64>, NumberExtError> {
    let capacity_error = || NumberExtError::Domain("Möbius sieve extent is too large".to_owned());
    let len = limit.checked_add(1).ok_or_else(capacity_error)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| capacity_error())?;
    // Two is a marker for integers not yet visited by the sieve.
    values.resize(len, 2);
    values[0] = 0;
    if limit > 0 {
        values[1] = 1;
    }
    let mut primes = Vec::new();
    for number in 2..=limit {
        if values[number] == 2 {
            primes.try_reserve(1).map_err(|_| capacity_error())?;
            primes.push(number);
            values[number] = -1;
        }
        for &prime in &primes {
            let Some(composite) = number.checked_mul(prime).filter(|&x| x <= limit) else {
                break;
            };
            if number.is_multiple_of(prime) {
                values[composite] = 0;
                break;
            }
            values[composite] = -values[number];
        }
    }
    Ok(values)
}

/// Integer partition function `p(n)` via Euler's pentagonal number theorem.
///
/// Negative arguments produce zero. The recurrence retains values through `n`.
///
/// # Panics
/// Panics when the table does not fit the host's address space.
#[must_use]
pub fn partition_count(n: i64) -> BigInt {
    if n < 0 {
        return BigInt::from(0u8);
    }
    let n = usize::try_from(n)
        .unwrap_or_else(|_| panic!("partition count table length overflows usize"));
    let len = n
        .checked_add(1)
        .unwrap_or_else(|| panic!("partition count table length overflows usize"));
    let mut table = vec![BigInt::from(0u8); len];
    table[0] = BigInt::from(1u8);
    for i in 1..=n {
        let mut total = BigInt::from(0u8);
        let mut k = 1usize;
        loop {
            let first = k * (3 * k - 1) / 2;
            if first > i {
                break;
            }
            if k % 2 == 1 {
                total += &table[i - first];
            } else {
                total -= &table[i - first];
            }
            let second = k * (3 * k + 1) / 2;
            if second <= i {
                if k % 2 == 1 {
                    total += &table[i - second];
                } else {
                    total -= &table[i - second];
                }
            }
            k += 1;
        }
        table[i] = total;
    }
    table.pop().unwrap_or_else(|| BigInt::from(0u8))
}

/// Unsigned Stirling numbers of the first kind `c(n, k)`.
///
/// Returns zero when `n < 0` or `k` lies outside `0..=n`.
///
/// # Panics
/// Panics when the recurrence's workspace length cannot fit the address space.
#[must_use]
pub fn stirling_first(
    n: i64,
    k: i64,
) -> BigInt {
    stirling::<true>(n, k)
}

/// Signed Stirling numbers of the first kind `s(n, k) = (-1)^(n - k) c(n, k)`.
///
/// Uses the same domain and workspace as [`stirling_first`].
#[must_use]
pub fn stirling_first_signed(
    n: i64,
    k: i64,
) -> BigInt {
    let unsigned = stirling_first(n, k);
    if (i128::from(n) - i128::from(k)) % 2 != 0 {
        -unsigned
    } else {
        unsigned
    }
}

/// Exact Stirling numbers of the second kind `S(n, k)`.
///
/// Returns zero when `n < 0` or `k` lies outside `0..=n`.
///
/// # Panics
/// Panics when the recurrence's workspace length cannot fit the address space.
#[must_use]
pub fn stirling_second_exact(
    n: i64,
    k: i64,
) -> BigInt {
    stirling::<false>(n, k)
}

fn stirling<const FIRST_KIND: bool>(
    n: i64,
    k: i64,
) -> BigInt {
    if n < 0 || k < 0 || k > n {
        return BigInt::default();
    }
    if n == k || (!FIRST_KIND && k == 1) {
        return BigInt::from(1u8);
    }
    if k == 0 {
        return BigInt::default();
    }
    let k = k as u64;
    let n = n as u64;
    let deficit = n - k;
    if deficit == 1 {
        return BigInt::from(u128::from(n) * u128::from(n - 1) / 2);
    }

    // Index the Stirling recurrence by (k, d = n-k). Each cell combines
    // (k-1, d) and (k, d-1), so only the k by d rectangle is needed. The
    // smaller axis is the rolling row. See DLMF 26.8.18 and 26.8.22.
    let width = usize::try_from(k.min(deficit) + 1)
        .unwrap_or_else(|_| panic!("Stirling number workspace width overflows usize"));
    let mut row;
    if k <= deficit {
        // At d = 0 all entries are diagonal ones. At k = 0 all
        // subsequent rows vanish.
        row = vec![BigUint::from(1u8); width];
        row[0] = BigUint::default();
        for d in 1..=deficit {
            for column in 1..width {
                let k = column as u64;
                let factor = if FIRST_KIND { k + d - 1 } else { k };
                let (before, after) = row.split_at_mut(column);
                after[0] *= factor;
                after[0] += &before[column - 1];
            }
        }
    } else {
        // At k = 0 only d = 0 is nonzero. Reuse the product allocation
        // while combining the left-hand value with the preceding row.
        row = vec![BigUint::default(); width];
        row[0] = BigUint::from(1u8);
        let mut product = BigUint::default();
        for k in 1..=k {
            for column in 1..width {
                let d = column as u64;
                let factor = if FIRST_KIND { k + d - 1 } else { k };
                let (before, after) = row.split_at_mut(column);
                product.clone_from(&before[column - 1]);
                product *= factor;
                after[0] += &product;
            }
        }
    }
    BigInt::from(row.pop().unwrap_or_default())
}

/// Exact `n`-th Catalan number `C_n = C(2n, n) / (n + 1)`.
/// Returns zero for negative arguments.
#[must_use]
pub fn catalan_exact(n: i64) -> BigInt {
    if n < 0 {
        return BigInt::from(0u8);
    }
    let n = n as u64;
    let mut binom = BigUint::from(1u8);
    for i in 0..n {
        binom *= 2 * n - i;
        binom /= i + 1;
    }
    BigInt::from(binom / (n + 1))
}

/// Modular multiplicative inverse in `[0, modulus)`.
///
/// # Errors
/// The modulus must be positive and relatively prime to `a`.
pub fn mod_inv(
    a: i64,
    modulus: i64,
) -> Result<i64, NumberExtError> {
    if modulus <= 0 {
        return Err(NumberExtError::Domain(
            "mod_inv requires a positive modulus".to_string(),
        ));
    }
    number_theory::mod_inverse(a, modulus)
        .ok_or_else(|| NumberExtError::Domain("value has no modular inverse".to_owned()))
}

#[cfg(test)]
mod tests;
