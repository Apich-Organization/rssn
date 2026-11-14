//! Classical orthogonal polynomials evaluated over exact rational inputs.
//!
//! Three-term recurrences preserve exactness for every nonnegative degree.
// Arbitrary-precision rational arithmetic cannot overflow; recurrence
// denominators are the positive integer order+1.
#![allow(clippy::arithmetic_side_effects)]
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

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

/// `T_n(x)` over the rationals.
#[must_use]
pub fn chebyshev_t_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    let mut previous = BigRational::one();
    let mut current = x.clone();
    for _ in 1..n {
        let next = rational(2) * x.clone() * current.clone() - previous.clone();
        previous = current;
        current = next;
    }
    current
}

/// `U_n(x)` over the rationals.
#[must_use]
pub fn chebyshev_u_exact(
    n: u32,
    x: &BigRational,
) -> BigRational {
    if n == 0 {
        return BigRational::one();
    }
    let mut previous = BigRational::one();
    let mut current = rational(2) * x.clone();
    for _ in 1..n {
        let next = rational(2) * x.clone() * current.clone() - previous.clone();
        previous = current;
        current = next;
    }
    current
}
