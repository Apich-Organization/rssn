//! # Symbolic Statistical Inference
//!
//! This module provides functions for symbolic statistical inference, particularly
//! focusing on hypothesis testing. It allows for the construction of symbolic
//! representations of test statistics and p-value formulas for various tests,
//! such as the two-sample t-test.
use crate::symbolic::core::Expr;
use crate::symbolic::stats::{mean, variance};
use std::sync::Arc;
/// Represents a formal hypothesis test.
#[derive(Debug, Clone)]
pub struct HypothesisTest {
    pub null_hypothesis: Expr,
    pub alternative_hypothesis: Expr,
    pub test_statistic: Expr,
    pub p_value_formula: Expr,
    pub degrees_of_freedom: Option<Expr>,
}
/// Constructs a symbolic two-sample t-test.
///
/// This function generates the symbolic formulas for the test statistic and degrees of freedom
/// for a two-sample t-test, assuming unequal variances (Welch's t-test).
///
/// # Arguments
/// * `sample1` - The first data sample as a slice of expressions.
/// * `sample2` - The second data sample as a slice of expressions.
/// * `mu_diff` - The hypothesized difference in means (often 0).
///
/// # Returns
/// A `HypothesisTest` struct containing the symbolic formulas for the test.
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
    let test_statistic = Expr::new_div(
        Expr::new_sub(Expr::new_sub(mean1.clone(), mean2.clone()), mu_diff.clone()),
        Expr::new_sqrt(Expr::new_add(
            Expr::new_div(var1.clone(), n1.clone()),
            Expr::new_div(var2.clone(), n2.clone()),
        )),
    );
    let term1 = Expr::new_div(var1, n1.clone());
    let term2 = Expr::new_div(var2, n2.clone());
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
    let p_value_formula = Expr::new_apply(
        Expr::Variable("t_dist_cdf".to_string()),
        Expr::Tuple(vec![test_statistic.clone(), df.clone()]),
    );
    HypothesisTest {
        null_hypothesis: Expr::Eq(
            Arc::new(Expr::Variable("mu1".to_string())),
            Arc::new(Expr::Variable("mu2".to_string())),
        ),
        alternative_hypothesis: Expr::new_not(Expr::Eq(
            Arc::new(Expr::Variable("mu1".to_string())),
            Arc::new(Expr::Variable("mu2".to_string())),
        )),
        test_statistic,
        p_value_formula,
        degrees_of_freedom: Some(df),
    }
}
