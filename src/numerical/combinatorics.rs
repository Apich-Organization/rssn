//! # Numerical Combinatorics
//!
//! This module provides numerical implementations for combinatorial functions.
//! It includes functions for calculating factorials, permutations, and combinations,
//! as well as a numerical solver for linear recurrence relations.

/// Computes the factorial of `n` (`n!`) as an `f64`.
///
/// # Arguments
/// * `n` - The non-negative integer for which to compute the factorial.
///
/// # Returns
/// The factorial of `n` as an `f64`. Returns `f64::INFINITY` if `n` is too large to fit in `f64`.
#[must_use]

pub fn factorial(n: u64) -> f64 {

    if n > 170 {

        return f64::INFINITY;
    }

    (1..=n)
        .map(|i| i as f64)
        .product()
}

/// Computes the number of permutations `P(n, k) = n! / (n-k)!`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of permutations as an `f64`. Returns `0.0` if `k > n`.
#[must_use]

pub fn permutations(
    n: u64,
    k: u64,
) -> f64 {

    if k > n {

        return 0.0;
    }

    (n - k + 1..=n)
        .map(|i| i as f64)
        .product()
}

/// Computes the number of combinations `C(n, k) = n! / (k! * (n-k)!)`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of combinations as an `f64`. Returns `0.0` if `k > n`.
#[must_use]

pub fn combinations(
    n: u64,
    k: u64,
) -> f64 {

    if k > n {

        return 0.0;
    }

    if k == 0 || k == n {

        return 1.0;
    }

    if k > n / 2 {

        return combinations(n, n - k);
    }

    let mut res = 1.0;

    for i in 1..=k {

        res = res * (n - i + 1) as f64
            / i as f64;
    }

    res
}

/// Numerically solves a linear recurrence relation by unrolling it.
///
/// The recurrence relation is assumed to be of the form:
/// `a_n = coeffs[0]*a_{n-1} + coeffs[1]*a_{n-2} + ... + coeffs[order-1]*a_{n-order}`.
///
/// # Arguments
/// * `coeffs` - A slice of `f64` representing the coefficients of the recurrence relation.
/// * `initial_conditions` - A slice of `f64` representing the initial values `a_0, a_1, ..., a_{order-1}`.
/// * `target_n` - The index `n` for which to compute `a_n`.
///
/// # Returns
/// A `Result` containing the numerical value of `a_n`, or an error string if input dimensions mismatch.

pub fn solve_recurrence_numerical(
    coeffs: &[f64],
    initial_conditions: &[f64],
    target_n: usize,
) -> Result<f64, String> {

    let order = coeffs.len();

    if initial_conditions.len() != order
    {

        return Err(
            "Number of initial \
             conditions must match \
             the order of the \
             recurrence."
                .to_string(),
        );
    }

    if target_n < order {

        return Ok(initial_conditions
            [target_n]);
    }

    let mut values =
        initial_conditions.to_vec();

    for n in order..=target_n {

        let mut next_val = 0.0;

        for i in 0..order {

            next_val += coeffs[i]
                * values[n - 1 - i];
        }

        values.push(next_val);
    }

    match values.last() {
        | Some(v) => Ok(*v),
        | None => {
            Err("Failed to compute \
                 the recurrence \
                 relation, values \
                 vector was empty."
                .to_string())
        },
    }
}

/// Computes the Stirling numbers of the second kind S(n, k).
/// This is the number of ways to partition a set of `n` elements into `k` non-empty subsets.
#[must_use]

pub fn stirling_second(
    n: u64,
    k: u64,
) -> f64 {

    if k > n {

        return 0.0;
    }

    if k == 0 {

        return if n == 0 {

            1.0
        } else {

            0.0
        };
    }

    if k == n {

        return 1.0;
    }

    if k == 1 {

        return 1.0;
    }

    // Recurrence: S(n, k) = k*S(n-1, k) + S(n-1, k-1)
    // For numerical implementation, explicit formula is better for large n, k:
    // S(n, k) = (1/k!) * sum_{j=0}^k (-1)^(k-j) * C(k, j) * j^n

    let mut sum = 0.0;

    for j in 0..=k {

        let term = combinations(k, j)
            * (j as f64).powf(n as f64);

        if (k - j) % 2 == 1 {

            sum -= term;
        } else {

            sum += term;
        }
    }

    sum / factorial(k)
}

/// Computes the Bell number B(n), which is the number of partitions of a set of `n` elements.
/// B(n) = sum_{k=0}^n S(n, k)
#[must_use]

pub fn bell(n: u64) -> f64 {

    (0..=n)
        .map(|k| stirling_second(n, k))
        .sum()
}

/// Computes the nth Catalan number C_n.
/// C_n = (1 / (n + 1)) * C(2n, n)
#[must_use]

pub fn catalan(n: u64) -> f64 {

    combinations(2 * n, n)
        / ((n + 1) as f64)
}

/// Computes the rising factorial (Pochhammer symbol) x^(n) = x(x+1)...(x+n-1).
#[must_use]

pub fn rising_factorial(
    x: f64,
    n: u64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    let mut res = 1.0;

    for i in 0..n {

        res *= x + (i as f64);
    }

    res
}

/// Computes the falling factorial (x)_n = x(x-1)...(x-n+1).
#[must_use]

pub fn falling_factorial(
    x: f64,
    n: u64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    let mut res = 1.0;

    for i in 0..n {

        res *= x - (i as f64);
    }

    res
}
