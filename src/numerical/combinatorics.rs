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
pub fn factorial(n: u64) -> f64 {
    if n > 170 {
        return f64::INFINITY;
    } // f64 overflows around 171!
    (1..=n).map(|i| i as f64).product()
}

/// Computes the number of permutations `P(n, k) = n! / (n-k)!`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of permutations as an `f64`. Returns `0.0` if `k > n`.
pub fn permutations(n: u64, k: u64) -> f64 {
    if k > n {
        return 0.0;
    }
    (n - k + 1..=n).map(|i| i as f64).product()
}

/// Computes the number of combinations `C(n, k) = n! / (k! * (n-k)!)`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// The number of combinations as an `f64`. Returns `0.0` if `k > n`.
pub fn combinations(n: u64, k: u64) -> f64 {
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
        res = res * (n - i + 1) as f64 / i as f64;
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
    if initial_conditions.len() != order {
        return Err(
            "Number of initial conditions must match the order of the recurrence.".to_string(),
        );
    }

    if target_n < order {
        return Ok(initial_conditions[target_n]);
    }

    let mut values = initial_conditions.to_vec();
    for n in order..=target_n {
        let mut next_val = 0.0;
        for i in 0..order {
            next_val += coeffs[i] * values[n - 1 - i];
        }
        values.push(next_val);
    }

    Ok(*values.last().unwrap())
}
