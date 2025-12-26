//! # Symbolic Statistical Inference
//!
//! This module provides functions for symbolic statistical inference, focusing on
//! hypothesis testing. It allows for the construction of symbolic formulas for
//! test statistics, degrees of freedom, and p-values.

use crate::symbolic::core::Expr;
use crate::symbolic::stats::{mean, variance};
use std::sync::Arc;

/// Represents a formal hypothesis test.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]

pub struct HypothesisTest {
    pub null_hypothesis: Expr,
    pub alternative_hypothesis: Expr,
    pub test_statistic: Expr,
    pub p_value_formula: Expr,
    pub degrees_of_freedom: Option<Expr>,
}

#[must_use]

pub fn one_sample_t_test_symbolic(sample: &[Expr], target_mean: &Expr) -> HypothesisTest {

    let n = Expr::Constant(sample.len() as f64);

    let mu = mean(sample);

    let var = variance(sample); // Population variance (Sum/n)

    // Unbiased sample variance s^2 = Sum/(n-1)
    // We have var = Sum/n, so s^2 = var * n / (n-1)
    // Standard Error = s / sqrt(n) = sqrt(s^2 / n) = sqrt( (var * n / (n-1)) / n ) = sqrt( var / (n-1) )

    let n_minus_1 = Expr::new_sub(n, Expr::Constant(1.0));

    let standard_error_sq = Expr::new_div(var, n_minus_1.clone());

    let standard_error = Expr::new_sqrt(standard_error_sq);

    let test_statistic = Expr::new_div(Expr::new_sub(mu, target_mean.clone()), standard_error);

    let df = n_minus_1;

    // Two-tailed p-value: 2 * (1 - CDF(|t|))
    let p_value = Expr::new_mul(
        Expr::Constant(2.0),
        Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_apply(
                Expr::Variable("t_dist_cdf".to_string()),
                Expr::Tuple(vec![
                    Expr::new_abs(test_statistic.clone()),
                    df.clone(),
                ]),
            ),
        ),
    );

    HypothesisTest {
        null_hypothesis: Expr::Eq(
            Arc::new(Expr::Variable("mu".to_string())),
            Arc::new(target_mean.clone()),
        ),
        alternative_hypothesis: Expr::new_not(Expr::Eq(
            Arc::new(Expr::Variable("mu".to_string())),
            Arc::new(target_mean.clone()),
        )),
        test_statistic,
        p_value_formula: p_value,
        degrees_of_freedom: Some(df),
    }
}

/// Constructs a symbolic two-sample t-test (Welch's t-test).
#[must_use]

pub fn two_sample_t_test_symbolic(
    sample1: &[Expr],
    sample2: &[Expr],
    mu_diff: &Expr,
) -> HypothesisTest {

    let n1 = Expr::Constant(sample1.len() as f64);

    let n2 = Expr::Constant(sample2.len() as f64);

    let mean1 = mean(sample1);

    let mean2 = mean(sample2);

    let var1 = variance(sample1);

    let var2 = variance(sample2);

    // Terms for SE and DF should use unbiased sample variance variance contribution: s^2 / n
    // s^2 / n = var / (n-1)

    let term1 = Expr::new_div(var1, Expr::new_sub(n1.clone(), Expr::Constant(1.0)));

    let term2 = Expr::new_div(var2, Expr::new_sub(n2.clone(), Expr::Constant(1.0)));

    let test_statistic = Expr::new_div(
        Expr::new_sub(Expr::new_sub(mean1, mean2), mu_diff.clone()),
        Expr::new_sqrt(Expr::new_add(term1.clone(), term2.clone())),
    );

    // Welch-Satterthwaite equation for df
    let df_num = Expr::new_pow(
        Expr::new_add(term1.clone(), term2.clone()),
        Expr::Constant(2.0),
    );

    let df_den1 = Expr::new_div(
        Expr::new_pow(term1, Expr::Constant(2.0)),
        Expr::new_sub(n1, Expr::Constant(1.0)),
    );

    let df_den2 = Expr::new_div(
        Expr::new_pow(term2, Expr::Constant(2.0)),
        Expr::new_sub(n2, Expr::Constant(1.0)),
    );

    let df = Expr::new_div(df_num, Expr::new_add(df_den1, df_den2));

    let p_value_formula = Expr::new_mul(
        Expr::Constant(2.0),
        Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_apply(
                Expr::Variable("t_dist_cdf".to_string()),
                Expr::Tuple(vec![
                    Expr::new_abs(test_statistic.clone()),
                    df.clone(),
                ]),
            ),
        ),
    );

    HypothesisTest {
        null_hypothesis: Expr::Eq(
            Arc::new(Expr::new_sub(
                Expr::Variable("mu1".to_string()),
                Expr::Variable("mu2".to_string()),
            )),
            Arc::new(mu_diff.clone()),
        ),
        alternative_hypothesis: Expr::new_not(Expr::Eq(
            Arc::new(Expr::new_sub(
                Expr::Variable("mu1".to_string()),
                Expr::Variable("mu2".to_string()),
            )),
            Arc::new(mu_diff.clone()),
        )),
        test_statistic,
        p_value_formula,
        degrees_of_freedom: Some(df),
    }
}

/// Constructs a symbolic Z-test.
///
/// Tests whether the mean of a population is equal to `target_mean` when the population
/// standard deviation is known.
#[must_use]

pub fn z_test_symbolic(sample: &[Expr], target_mean: &Expr, pop_std_dev: &Expr) -> HypothesisTest {

    let n = Expr::Constant(sample.len() as f64);

    let mu = mean(sample);

    let standard_error = Expr::new_div(pop_std_dev.clone(), Expr::new_sqrt(n));

    let test_statistic = Expr::new_div(Expr::new_sub(mu, target_mean.clone()), standard_error);

    let p_value = Expr::new_mul(
        Expr::Constant(2.0),
        Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_apply(
                Expr::Variable("normal_cdf".to_string()),
                Expr::Tuple(vec![Expr::new_abs(
                    test_statistic.clone(),
                )]),
            ),
        ),
    );

    HypothesisTest {
        null_hypothesis: Expr::Eq(
            Arc::new(Expr::Variable("mu".to_string())),
            Arc::new(target_mean.clone()),
        ),
        alternative_hypothesis: Expr::new_not(Expr::Eq(
            Arc::new(Expr::Variable("mu".to_string())),
            Arc::new(target_mean.clone()),
        )),
        test_statistic,
        p_value_formula: p_value,
        degrees_of_freedom: None,
    }
}
