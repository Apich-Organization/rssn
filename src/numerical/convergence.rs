//! # Numerical Convergence Analysis
//!
//! This module provides numerical methods for analyzing and accelerating the convergence
//! of series and sequences. It includes functions for summing series up to a given tolerance
//! and for accelerating sequence convergence using techniques like Aitken's delta-squared process.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

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
/// A `Result` containing the numerical sum, or an error string if evaluation fails.
pub fn sum_series_numerical(
    term_expr: &Expr,
    var: &str,
    start_n: usize,
    max_terms: usize,
    tolerance: f64,
) -> Result<f64, String> {
    let mut sum = 0.0;
    let mut vars = HashMap::new();
    for i in start_n..(start_n + max_terms) {
        vars.insert(var.to_string(), i as f64);
        let term_val = eval_expr(term_expr, &vars)?;
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
pub fn aitken_acceleration(sequence: &[f64]) -> Vec<f64> {
    if sequence.len() < 3 {
        return vec![];
    }
    let mut accelerated_seq = Vec::new();
    for i in 0..(sequence.len() - 2) {
        let s_n = sequence[i];
        let s_n1 = sequence[i + 1];
        let s_n2 = sequence[i + 2];
        let denominator = s_n2 - 2.0 * s_n1 + s_n;
        if denominator.abs() > 1e-9 {
            let aitken_s = s_n - (s_n1 - s_n).powi(2) / denominator;
            accelerated_seq.push(aitken_s);
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
/// A `Result` containing the numerical limit, or an error string if convergence is not found.
pub fn find_sequence_limit(
    term_expr: &Expr,
    var: &str,
    max_terms: usize,
    tolerance: f64,
) -> Result<f64, String> {
    let mut sequence = Vec::new();
    let mut vars = HashMap::new();
    for i in 0..max_terms {
        vars.insert(var.to_string(), i as f64);
        sequence.push(eval_expr(term_expr, &vars)?);
    }

    let mut accelerated = aitken_acceleration(&sequence);
    while accelerated.len() > 1 {
        let last = accelerated.last().unwrap();
        let second_last = accelerated[accelerated.len() - 2];
        if (last - second_last).abs() < tolerance {
            return Ok(*last);
        }
        accelerated = aitken_acceleration(&accelerated);
    }

    accelerated
        .last()
        .cloned()
        .ok_or_else(|| "Convergence not found".to_string())
}
