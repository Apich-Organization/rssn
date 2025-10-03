//! # Symbolic Statistical Inference
//!
//! This module provides functions for symbolic statistical inference, particularly
//! focusing on hypothesis testing. It allows for the construction of symbolic
//! representations of test statistics and p-value formulas for various tests,
//! such as the two-sample t-test.

use crate::symbolic::core::Expr;
use crate::symbolic::stats::{mean, variance};

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
    mu_diff: Expr,
) -> HypothesisTest {
    let n1 = Expr::Constant(sample1.len() as f64);
    let n2 = Expr::Constant(sample2.len() as f64);

    let mean1 = mean(sample1);
    let mean2 = mean(sample2);
    let var1 = variance(sample1);
    let var2 = variance(sample2);

    // t = ( (mean1 - mean2) - mu_diff ) / sqrt(var1/n1 + var2/n2)
    let test_statistic = Expr::Div(
        Box::new(Expr::Sub(
            Box::new(Expr::Sub(Box::new(mean1.clone()), Box::new(mean2.clone()))),
            Box::new(mu_diff.clone()),
        )),
        Box::new(Expr::Sqrt(Box::new(Expr::Add(
            Box::new(Expr::Div(Box::new(var1.clone()), Box::new(n1.clone()))),
            Box::new(Expr::Div(Box::new(var2.clone()), Box::new(n2.clone()))),
        )))),
    );

    // Degrees of freedom for Welch's t-test (Satterthwaite equation)
    let term1 = Expr::Div(Box::new(var1), Box::new(n1.clone()));
    let term2 = Expr::Div(Box::new(var2), Box::new(n2.clone()));
    let df_num = Expr::Power(
        Box::new(Expr::Add(Box::new(term1.clone()), Box::new(term2.clone()))),
        Box::new(Expr::Constant(2.0)),
    );
    let df_den1 = Expr::Div(
        Box::new(Expr::Power(Box::new(term1), Box::new(Expr::Constant(2.0)))),
        Box::new(Expr::Sub(Box::new(n1), Box::new(Expr::Constant(1.0)))),
    );
    let df_den2 = Expr::Div(
        Box::new(Expr::Power(Box::new(term2), Box::new(Expr::Constant(2.0)))),
        Box::new(Expr::Sub(Box::new(n2), Box::new(Expr::Constant(1.0)))),
    );
    let df = Expr::Div(
        Box::new(df_num),
        Box::new(Expr::Add(Box::new(df_den1), Box::new(df_den2))),
    );

    // p-value is the CDF of the t-distribution. We represent this symbolically.
    let p_value_formula = Expr::Apply(
        Box::new(Expr::Variable("t_dist_cdf".to_string())),
        Box::new(Expr::Tuple(vec![test_statistic.clone(), df.clone()])),
    );

    HypothesisTest {
        null_hypothesis: Expr::Eq(
            Box::new(Expr::Variable("mu1".to_string())),
            Box::new(Expr::Variable("mu2".to_string())),
        ),
        alternative_hypothesis: Expr::Not(Box::new(Expr::Eq(
            Box::new(Expr::Variable("mu1".to_string())),
            Box::new(Expr::Variable("mu2".to_string())),
        ))),
        test_statistic,
        p_value_formula,
        degrees_of_freedom: Some(df),
    }
}
