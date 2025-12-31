//! # Numerical Convergence Analysis
//!
//! This module provides numerical methods for analyzing and accelerating the convergence
//! of series and sequences. It includes functions for summing series up to a given tolerance
//! and for accelerating sequence convergence using techniques like Aitken's delta-squared process.

use std::collections::HashMap;

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;

/// Numerically sums a series until the term is smaller than a tolerance or `max_terms` is reached.
///
/// # Arguments
/// * `term_expr` - The symbolic expression for the `n`-th term of the series.
/// * `var` - The variable representing the index `n`.
/// * `start_n` - The starting index for the summation.
/// * `max_terms` - The maximum number of terms to sum.
/// * `tolerance` - The threshold for the absolute value of a term, below which summation stops.
///
/// # Returns
/// A `Result` containing the numerical sum.
///
/// # Errors
///
/// Returns an error if the expression evaluation fails for any term.

pub fn sum_series_numerical(
    term_expr: &Expr,
    var: &str,
    start_n: usize,
    max_terms: usize,
    tolerance: f64,
) -> Result<f64, String> {

    let mut sum = 0.0;

    let mut vars = HashMap::new();

    for i in
        start_n .. (start_n + max_terms)
    {

        vars.insert(
            var.to_string(),
            i as f64,
        );

        let term_val = eval_expr(
            term_expr,
            &vars,
        )?;

        if term_val.abs() < tolerance {

            break;
        }

        sum += term_val;
    }

    Ok(sum)
}

/// Accelerates the convergence of a sequence using Aitken's delta-squared process.
///
/// Aitken's method is a sequence acceleration technique that can improve the rate
/// of convergence of a slowly converging sequence. It is particularly effective
/// for linearly converging sequences.
///
/// # Arguments
/// * `sequence` - The original sequence `s_n`.
///
/// # Returns
/// A `Vec<f64>` representing the accelerated sequence.
#[must_use]

pub fn aitken_acceleration(
    sequence: &[f64]
) -> Vec<f64> {

    if sequence.len() < 3 {

        return vec![];
    }

    let mut accelerated_seq =
        Vec::new();

    for i in 0 .. (sequence.len() - 2) {

        let s_n = sequence[i];

        let s_n1 = sequence[i + 1];

        let s_n2 = sequence[i + 2];

        let denominator = 2.0f64
            .mul_add(-s_n1, s_n2)
            + s_n;

        if denominator.abs() > 1e-9 {

            let aitken_s = s_n
                - (s_n1 - s_n).powi(2)
                    / denominator;

            accelerated_seq
                .push(aitken_s);
        }
    }

    accelerated_seq
}

/// Numerically finds the limit of a sequence by generating terms and applying acceleration.
///
/// This function generates terms of a sequence defined by `term_expr` and then repeatedly
/// applies Aitken's delta-squared process to accelerate its convergence until a specified
/// tolerance is met or `max_terms` is reached.
///
/// # Arguments
/// * `term_expr` - The symbolic expression for the `n`-th term of the sequence.
/// * `var` - The variable representing the index `n`.
/// * `max_terms` - The maximum number of terms to generate before giving up.
/// * `tolerance` - The desired precision for the limit.
///
/// # Returns
/// A `Result` containing the numerical limit.
///
/// # Errors
///
/// Returns an error if:
/// - The expression evaluation fails for any generated term.
/// - The sequence becomes empty unexpectedly during acceleration.
/// - Convergence is not found within the specified tolerance and maximum terms.

pub fn find_sequence_limit(
    term_expr: &Expr,
    var: &str,
    max_terms: usize,
    tolerance: f64,
) -> Result<f64, String> {

    let mut sequence = Vec::new();

    let mut vars = HashMap::new();

    for i in 0 .. max_terms {

        vars.insert(
            var.to_string(),
            i as f64,
        );

        sequence.push(eval_expr(
            term_expr,
            &vars,
        )?);
    }

    let mut accelerated =
        aitken_acceleration(&sequence);

    while accelerated.len() > 1 {

        let last =
            match accelerated.last() {
                | Some(l) => l,
                | None => {

                    return Err(
                    "Unexpected empty \
                     sequence in \
                     convergence loop."
                        .to_string(),
                );
                },
            };

        let second_last = accelerated
            [accelerated.len() - 2];

        if (last - second_last).abs()
            < tolerance
        {

            return Ok(*last);
        }

        accelerated =
            aitken_acceleration(
                &accelerated,
            );
    }

    accelerated
        .last()
        .copied()
        .ok_or_else(|| {

            "Convergence not found"
                .to_string()
        })
}

/// Performs Richardson extrapolation on a sequence of approximations.
///
/// This method assumes the input sequence consists of approximations $A(h), A(h/2), A(h/4), \dots$,
/// where the error term is $O(h^k)$. This is commonly used in numerical differentiation and
/// Romberg integration.
///
/// # Arguments
/// * `sequence` - A slice of approximations with halving step sizes.
///
/// # Returns
/// A `Vec<f64>` containing the extrapolated values. The last element is the highest order extrapolation.
#[must_use]

pub fn richardson_extrapolation(
    sequence: &[f64]
) -> Vec<f64> {

    if sequence.is_empty() {

        return vec![];
    }

    let n = sequence.len();

    let mut table =
        vec![vec![0.0; n]; n];

    // Initialize the first column with the input sequence
    for (i, &val) in sequence
        .iter()
        .enumerate()
    {

        table[i][0] = val;
    }

    // Compute the extrapolation table
    for j in 1 .. n {

        for i in j .. n {

            let power_of_4 =
                4.0f64.powi(j as i32);

            table[i][j] =
                power_of_4.mul_add(
                    table[i][j - 1],
                    -table[i - 1]
                        [j - 1],
                ) / (power_of_4 - 1.0);
        }
    }

    // The diagonal elements are the best approximations for each order
    (0 .. n)
        .map(|i| table[i][i])
        .collect()
}

/// Applies Wynn's epsilon algorithm to accelerate the convergence of a sequence.
///
/// This is a recursive scheme to implement the Shanks transformation, which effectively
/// models the sequence as a sum of exponentials/transients.
///
/// # Arguments
/// * `sequence` - The sequence to accelerate.
///
/// # Returns
/// A `Vec<f64>` of accelerated terms.
#[must_use]

pub fn wynn_epsilon(
    sequence: &[f64]
) -> Vec<f64> {

    let n = sequence.len();

    if n < 3 {

        return Vec::from(sequence);
    }

    // Epsilon table
    // eps[k] stores the current column's value for the current row
    // We need to store previous columns to compute the next.
    // Wynn's epsilon algorithm:
    // eps(-1, n) = 0
    // eps(0, n) = S_n (original sequence)
    // eps(k+1, n) = eps(k-1, n+1) + 1 / (eps(k, n+1) - eps(k, n))

    // Just implementing a simple version that returns the diagonal or best equivalents.
    // However, the table is 2D. Let's return the simplified "shanks" equivalent column if possible,
    // or just the best single estimate.
    // For general utility, returning a vector of "improved" estimates is good.
    // Let's perform the algorithm and return the values from the even columns (which correspond to sequence estimates).

    // k=0: Original sequence
    // k=1: 1/(dS) ...

    // We limit k to n.
    // Re-initialization for clarity:
    let mut eps =
        vec![vec![0.0; n]; n + 1];

    eps[0][.. n].copy_from_slice(
        &sequence[.. n],
    );

    for k in 0 .. n - 1 {

        for i in 0 .. n - k - 1 {

            let numerator = 1.0;

            let denominator = eps[k]
                [i + 1]
                - eps[k][i];

            if denominator.abs() < 1e-12
            {

                // If denominator is too small, we might have convergence or numerical instability.
                // We propagate the previous value or stop.
                // For this implementation, let's just use a very large number relative to the prev to avoid NaN,
                // or break.
                eps[k + 1][i] = 1e12; // Placeholder for infinity
            } else {

                let prev_term = if k
                    == 0
                {

                    0.0
                } else {

                    eps[k - 1][i + 1]
                };

                eps[k + 1][i] = prev_term + numerator / denominator;
            }
        }
    }

    // The even columns k=0, 2, 4... represent sequence approximations.
    // The odd columns k=1, 3... are auxiliary quantities.
    // We want to return the best estimates.
    // Usually the logic is to look at the "lower diagonal" or the last computed even column.

    // Let's collect the values from the highest available even k for each index.


    // Note: Aitken is eps[2]. So this gives the same as Aitken.
    // To give more power, we should perhaps return the "diagonal": eps[2k][0].

    let mut diag = Vec::new();

    let mut k = 0;

    loop {

        if k >= n {

            break;
        }

        diag.push(eps[k][0]);

        k += 2;
    }

    if diag.is_empty() {

        Vec::from(sequence)
    } else {

        diag
    }
}
