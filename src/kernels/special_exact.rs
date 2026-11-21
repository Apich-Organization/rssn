//! Classical orthogonal polynomials over the rationals and modular Chebyshev values.
//!
//! Three-term recurrences and binary doubling preserve exactness for every
//! nonnegative degree.
// Arbitrary-precision rational arithmetic cannot overflow; recurrence
// denominators are the positive integer order+1.
#![allow(clippy::arithmetic_side_effects)]
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

fn rational(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

/// `P_n(x)` over the rationals.
#[must_use]
pub fn legendre_p_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    let mut previous = BigRational::one();
    let mut current = x.clone();
    for order in 1..n {
        let k = rational(i64::from(order));
        let next = ((rational(2) * k.clone() + BigRational::one()) * x.clone() * current.clone()
            - k.clone() * previous.clone())
            / (k + BigRational::one());
        previous = current;
        current = next;
    }
    current
}

/// `H_n(x)` over the rationals.
#[must_use]
pub fn hermite_h_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    let mut previous = BigRational::one();
    let mut current = rational(2) * x.clone();
    for order in 1..n {
        let next = rational(2) * x.clone() * current.clone()
            - rational(2 * i64::from(order)) * previous.clone();
        previous = current;
        current = next;
    }
    current
}

/// `L_n(x)` over the rationals.
#[must_use]
pub fn laguerre_l_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    let mut previous = BigRational::one();
    let mut current = BigRational::one() - x.clone();
    for order in 1..n {
        let k = rational(i64::from(order));
        let next = ((rational(2) * k.clone() + BigRational::one() - x.clone()) * current.clone()
            - k.clone() * previous.clone())
            / (k + BigRational::one());
        previous = current;
        current = next;
    }
    current
}

/// Chebyshev polynomial of the first kind `T_n(x)` over the rationals.
///
/// Uses logarithmically many integer operations before reduction.
#[must_use]
pub fn chebyshev_t_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    chebyshev_exact(n, x, false)
}

/// Chebyshev polynomial of the second kind `U_n(x)` over the rationals.
///
/// Uses logarithmically many integer operations before reduction.
#[must_use]
pub fn chebyshev_u_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    chebyshev_exact(n, x, true)
}

fn chebyshev_exact(
    n: u32,
    x: &BigRational,
    second_kind: bool,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    if n == 1 {
        return if second_kind {
            x * rational(2)
        } else {
            x.clone()
        };
    }
    let (p, q) = (x.numer(), x.denom());
    if p.is_zero() {
        return rational([1, 0, -1, 0][(n % 4) as usize]);
    }
    let sign = if n % 2 == 1 && x.is_negative() {
        -1
    } else {
        1
    };
    if p.magnitude() == q.magnitude() {
        let magnitude = if second_kind {
            i64::from(n) + 1
        } else {
            1
        };
        return rational(sign * magnitude);
    }
    if (p.magnitude() << 1_usize) == *q.magnitude() {
        // x = cos(pi / 3) gives a period of six. Parity handles -x.
        let i = (n % 6) as usize;
        return if second_kind {
            rational(sign * [1, 1, 0, -1, -1, 0][i])
        } else {
            BigRational::new(BigInt::from(sign * [2, 1, -1, -2, -1, 1][i]), 2.into())
        };
    }

    // With V_k = q^k U_k(p/q), V_(-1) = 0 and V_0 = 1, the doubling
    // identities are V_(2k) = V_k² - q² V_(k-1)² and
    // V_(2k-1) = 2 V_(k-1) (V_k - p V_(k-1)).
    let q_squared = q * q;
    let (mut previous, mut current) = (BigInt::zero(), BigInt::one());
    for bit in (0..n.bit_width()).rev() {
        let even = &current * &current - &q_squared * (&previous * &previous);
        let before = (&previous * (&current - p * &previous)) << 1_usize;
        if (n >> bit) & 1 == 0 {
            (previous, current) = (before, even);
        } else {
            let next = ((p * &even) << 1_usize) - &q_squared * before;
            (previous, current) = (even, next);
        }
    }
    // T_n = U_n - x U_(n-1). Both values have denominator q^n.
    let numerator = if second_kind {
        current
    } else {
        current - p * previous
    };
    BigRational::new(numerator, q.pow(n))
}

/// Returns `(T_n(x), U_n(x))` reduced modulo a positive modulus.
///
/// The residues are in `0..modulus`, including for negative `x` and
/// composite moduli. Runs in `O(log n)` modular operations for any
/// `u64` degree. A zero modulus gives `None`.
#[must_use]
pub fn chebyshev_tu_mod(
    n: u64,
    x: i64,
    modulus: u64,
) -> Option<(u64, u64)> {
    if modulus == 0 {
        return None;
    }
    if modulus == 1 {
        return Some((0, 0));
    }
    let x = i128::from(x).rem_euclid(i128::from(modulus)) as u64;
    let wide = u128::from(modulus);
    let product = |a, b| ((u128::from(a) * u128::from(b)) % wide) as u64;
    let difference = |a, b| {
        if a >= b {
            a - b
        } else {
            modulus - (b - a)
        }
    };
    let two = 2 % modulus;
    let twice_x = product(two, x);
    let (mut previous, mut current) = (0, 1);
    for bit in (0..n.bit_width()).rev() {
        let even = difference(product(current, current), product(previous, previous));
        let before = product(
            two,
            product(previous, difference(current, product(x, previous))),
        );
        if (n >> bit) & 1 == 0 {
            (previous, current) = (before, even);
        } else {
            (previous, current) = (even, difference(product(twice_x, even), before));
        }
    }
    Some((difference(current, product(x, previous)), current))
}
