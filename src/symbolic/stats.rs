//! # Symbolic Statistics
//!
//! This module provides functions for symbolic statistical calculations.
//! It includes basic descriptive statistics such as mean, variance, standard deviation,
//! covariance, and correlation, all expressed symbolically.

use num_bigint::BigInt;
use crate::symbolic::core::Expr;
use crate::symbolic::egraph::simplify;

/// Computes the symbolic mean of a set of expressions.
///
/// Encapsulated as an E-Graph Facade evaluated via the saturation pipeline.
#[must_use]
pub fn mean(data: &[Expr]) -> Expr {
    let call = Expr::NaryList("stats_mean".to_string(), data.to_vec());
    crate::symbolic::egraph::simplify(&call)
}

/// Internal solver for symbolic mean.
#[must_use]
pub fn mean_internal(data: &[Expr]) -> Expr {
    let n = data.len();

    if n == 0 {
        return Expr::Constant(0.0);
    }

    let sum = data
        .iter()
        .cloned()
        .reduce(|acc, e| simplify(&Expr::new_add(acc, e)))
        .unwrap_or(Expr::Constant(0.0));

    simplify(&Expr::new_div(sum, Expr::Constant(n as f64)))
}

/// Computes the symbolic variance of a set of expressions.
///
/// Encapsulated as an E-Graph Facade evaluated via the saturation pipeline.
#[must_use]
pub fn variance(data: &[Expr]) -> Expr {
    let call = Expr::NaryList("stats_variance".to_string(), data.to_vec());
    crate::symbolic::egraph::simplify(&call)
}

/// Internal solver for symbolic variance.
#[must_use]
pub fn variance_internal(data: &[Expr]) -> Expr {
    let n = data.len();

    if n == 0 {
        return Expr::Constant(0.0);
    }

    let mu = mean_internal(data);

    let squared_diffs = data
        .iter()
        .map(|x_i| {
            let diff = Expr::new_sub(x_i.clone(), mu.clone());

            Expr::new_pow(diff, Expr::Constant(2.0))
        })
        .reduce(|acc, e| simplify(&Expr::new_add(acc, e)))
        .unwrap_or(Expr::Constant(0.0));

    simplify(&Expr::new_div(squared_diffs, Expr::Constant(n as f64)))
}

/// Computes the symbolic standard deviation of a set of expressions.
///
/// Encapsulated as an E-Graph Facade evaluated via the saturation pipeline.
#[must_use]
pub fn std_dev(data: &[Expr]) -> Expr {
    let call = Expr::NaryList("stats_std_dev".to_string(), data.to_vec());
    crate::symbolic::egraph::simplify(&call)
}

/// Internal solver for symbolic standard deviation.
#[must_use]
pub fn std_dev_internal(data: &[Expr]) -> Expr {
    simplify(&Expr::new_sqrt(variance_internal(data)))
}

/// Computes the symbolic covariance of two sets of expressions.
///
/// Encapsulated as an E-Graph Facade evaluated via the saturation pipeline.
#[must_use]
pub fn covariance(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {
    let mut args = vec![Expr::BigInt(BigInt::from(data1.len()))];
    args.extend(data1.iter().cloned());
    args.extend(data2.iter().cloned());
    let call = Expr::NaryList("stats_covariance".to_string(), args);
    crate::symbolic::egraph::simplify(&call)
}

/// Internal solver for symbolic covariance.
#[must_use]
pub fn covariance_internal(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {
    if data1.len() != data2.len() || data1.is_empty() {
        return Expr::Constant(0.0);
    }

    let n = data1.len();

    let mu_x = mean_internal(data1);

    let mu_y = mean_internal(data2);

    let sum_of_products = data1
        .iter()
        .zip(data2.iter())
        .map(|(x_i, y_i)| {
            let diff_x = Expr::new_sub(x_i.clone(), mu_x.clone());

            let diff_y = Expr::new_sub(y_i.clone(), mu_y.clone());

            Expr::new_mul(diff_x, diff_y)
        })
        .reduce(|acc, e| simplify(&Expr::new_add(acc, e)))
        .unwrap_or(Expr::Constant(0.0));

    simplify(&Expr::new_div(sum_of_products, Expr::Constant(n as f64)))
}

/// Computes the symbolic Pearson correlation coefficient.
///
/// Encapsulated as an E-Graph Facade evaluated via the saturation pipeline.
#[must_use]
pub fn correlation(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {
    let mut args = vec![Expr::BigInt(BigInt::from(data1.len()))];
    args.extend(data1.iter().cloned());
    args.extend(data2.iter().cloned());
    let call = Expr::NaryList("stats_correlation".to_string(), args);
    crate::symbolic::egraph::simplify(&call)
}

/// Internal solver for symbolic Pearson correlation coefficient.
#[must_use]
pub fn correlation_internal(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {
    let cov_xy = covariance_internal(data1, data2);

    let std_dev_x = std_dev_internal(data1);

    let std_dev_y = std_dev_internal(data2);

    simplify(&Expr::new_div(cov_xy, Expr::new_mul(std_dev_x, std_dev_y)))
}
