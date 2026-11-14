//! Extended number theory and exact combinatorics kernels.

use num_bigint::BigInt;

/// Error returned by extended number-theoretic operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NumberExtError {
    /// Input slices were empty.
    Empty,
    /// Input slice lengths did not match.
    ShapeMismatch(String),
    /// Invalid modulus or inconsistent congruences.
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

fn modular_power_u64(
    base: u64,
    mut exp: u64,
    modulus: u64,
) -> u64 {
    let mut acc = 1u128;
    let mut b = u128::from(base % modulus);
    let m = u128::from(modulus);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = (acc * b) % m;
        }
        b = (b * b) % m;
        exp >>= 1;
    }
    acc as u64
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
    let m = modulus as u64;
    let b = base.rem_euclid(modulus) as u64;
    Ok(modular_power_u64(b, exponent, m) as i64)
}

/// Chinese Remainder Theorem: solves `x ≡ remainders[i] (mod moduli[i])`.
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

/// Möbius function `μ(n)`.
#[must_use]
pub fn mobius(n: i64) -> i64 {
    if n <= 0 {
        return 0;
    }
    if n == 1 {
        return 1;
    }
    let mut factors = crate::kernels::number_theory::factorize(n as u64);
    factors.sort_unstable();
    if factors.windows(2).any(|w| w[0] == w[1]) {
        return 0;
    }
    if factors.len().is_multiple_of(2) {
        1
    } else {
        -1
    }
}

/// Integer partition function `p(n)` via Euler's pentagonal number theorem.
#[must_use]
pub fn partition_count(n: i64) -> BigInt {
    if n < 0 {
        return BigInt::from(0u8);
    }
    let n = n as usize;
    let mut table = vec![BigInt::from(0u8); n + 1];
    table[0] = BigInt::from(1u8);
    for i in 1..=n {
        let mut total = BigInt::from(0u8);
        let mut k = 1usize;
        loop {
            let g1 = k * (3 * k - 1) / 2;
            if g1 > i {
                break;
            }
            if k % 2 == 1 {
                total += &table[i - g1];
            } else {
                total -= &table[i - g1];
            }
            let g2 = k * (3 * k + 1) / 2;
            if g2 <= i {
                if k % 2 == 1 {
                    total += &table[i - g2];
                } else {
                    total -= &table[i - g2];
                }
            }
            k += 1;
        }
        table[i] = total;
    }
    table.pop().unwrap_or_else(|| BigInt::from(0u8))
}

/// Unsigned Stirling numbers of the first kind `c(n, k)`.
#[must_use]
pub fn stirling_first(
    n: i64,
    k: i64,
) -> BigInt {
    if n < 0 || k < 0 || k > n {
        return BigInt::from(0u8);
    }
    if n == 0 && k == 0 {
        return BigInt::from(1u8);
    }
    if k == 0 {
        return BigInt::from(0u8);
    }
    let n = n as usize;
    let k = k as usize;
    let mut row = vec![BigInt::from(0u8); k + 1];
    row[0] = BigInt::from(1u8);
    for i in 1..=n {
        let factor = BigInt::from(i - 1);
        for j in (1..=k.min(i)).rev() {
            row[j] = &row[j - 1] + &factor * &row[j];
        }
        if i == 1 {
            row[0] = BigInt::from(0u8);
        }
    }
    row[k].clone()
}

/// Signed Stirling numbers of the first kind `s(n, k) = (-1)^(n - k) c(n, k)`.
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
#[must_use]
pub fn stirling_second_exact(
    n: i64,
    k: i64,
) -> BigInt {
    if n < 0 || k < 0 || k > n {
        return BigInt::from(0u8);
    }
    if n == 0 && k == 0 {
        return BigInt::from(1u8);
    }
    if k == 0 {
        return BigInt::from(0u8);
    }
    let n = n as usize;
    let k = k as usize;
    let mut row = vec![BigInt::from(0u8); k + 1];
    row[0] = BigInt::from(1u8);
    for i in 1..=n {
        for j in (1..=k.min(i)).rev() {
            row[j] = &row[j - 1] + BigInt::from(j) * &row[j];
        }
        if i == 1 {
            row[0] = BigInt::from(0u8);
        }
    }
    row[k].clone()
}

/// Exact `n`-th Catalan number `C_n = C(2n, n) / (n + 1)`.
#[must_use]
pub fn catalan_exact(n: i64) -> BigInt {
    if n < 0 {
        return BigInt::from(0u8);
    }
    let n = n as u64;
    let mut binom = BigInt::from(1u8);
    for i in 0..n {
        binom = (binom * BigInt::from(2 * n - i)) / BigInt::from(i + 1);
    }
    binom / BigInt::from(n + 1)
}

/// Modular multiplicative inverse in `[0, modulus)`.
pub fn mod_inv(
    a: i64,
    modulus: i64,
) -> Result<i64, NumberExtError> {
    if modulus <= 0 {
        return Err(NumberExtError::Domain(
            "mod_inv requires a positive modulus".to_string(),
        ));
    }
    let (g, x, _) = xgcd(a, modulus);
    if g != 1 {
        return Err(NumberExtError::Domain(
            "value has no modular inverse".to_string(),
        ));
    }
    Ok(x.rem_euclid(i128::from(modulus)) as i64)
}
