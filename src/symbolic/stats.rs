//! # Symbolic Statistics
//!
//! This module provides functions for symbolic statistical calculations.
//! It includes basic descriptive statistics such as mean, variance, standard deviation,
//! covariance, and correlation, all expressed symbolically.

use std::sync::Arc;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;

/// Computes the symbolic mean of a set of expressions.
///
/// The mean (average) is a measure of central tendency. For a set of `n` data points `x_i`,
/// it is defined as `(1/n) * Σx_i`.
///
/// # Arguments
/// * `data` - A slice of `Expr` representing the data points.
///
/// # Returns
/// An `Expr` representing the symbolic mean.
pub fn mean(data: &[Expr]) -> Expr {
    let n = data.len();
    if n == 0 {
        return Expr::Constant(0.0);
    }
    let sum = data
        .iter()
        .cloned()
        .reduce(|acc, e| simplify(Expr::Add(Arc::new(acc), Arc::new(e))))
        .unwrap_or(Expr::Constant(0.0));
    simplify(Expr::Div(Arc::new(sum), Arc::new(Expr::Constant(n as f64))))
}

/// Computes the symbolic variance of a set of expressions.
///
/// The variance is a measure of the spread or dispersion of a set of data.
/// For a set of `n` data points `x_i` with mean `μ`, it is defined as `(1/n) * Σ(x_i - μ)²`.
///
/// # Arguments
/// * `data` - A slice of `Expr` representing the data points.
///
/// # Returns
/// An `Expr` representing the symbolic variance.
pub fn variance(data: &[Expr]) -> Expr {
    let n = data.len();
    if n == 0 {
        return Expr::Constant(0.0);
    }
    let mu = mean(data);
    let squared_diffs = data
        .iter()
        .map(|x_i| {
            let diff = Expr::Sub(Arc::new(x_i.clone()), Arc::new(mu.clone()));
            Expr::Power(Arc::new(diff), Arc::new(Expr::Constant(2.0)))
        })
        .reduce(|acc, e| simplify(Expr::Add(Arc::new(acc), Arc::new(e))))
        .unwrap_or(Expr::Constant(0.0));
    simplify(Expr::Div(
        Arc::new(squared_diffs),
        Arc::new(Expr::Constant(n as f64)),
    ))
}

/// Computes the symbolic standard deviation of a set of expressions.
///
/// The standard deviation is the square root of the variance, providing a measure
/// of data dispersion in the same units as the data itself.
///
/// # Arguments
/// * `data` - A slice of `Expr` representing the data points.
///
/// # Returns
/// An `Expr` representing the symbolic standard deviation.
pub fn std_dev(data: &[Expr]) -> Expr {
    simplify(Expr::Sqrt(Arc::new(variance(data))))
}

/// Computes the symbolic covariance of two sets of expressions.
///
/// Covariance measures the joint variability of two random variables. For two sets
/// of `n` data points `x_i` and `y_i` with means `μ_x` and `μ_y`, it is defined as
/// `(1/n) * Σ((x_i - μ_x) * (y_i - μ_y))`.
///
/// # Arguments
/// * `data1` - The first set of data points.
/// * `data2` - The second set of data points.
///
/// # Returns
/// An `Expr` representing the symbolic covariance.
pub fn covariance(data1: &[Expr], data2: &[Expr]) -> Expr {
    if data1.len() != data2.len() || data1.is_empty() {
        return Expr::Constant(0.0); // Or return an error
    }
    let n = data1.len();
    let mu_x = mean(data1);
    let mu_y = mean(data2);
    let sum_of_products = data1
        .iter()
        .zip(data2.iter())
        .map(|(x_i, y_i)| {
            let diff_x = Expr::Sub(Arc::new(x_i.clone()), Arc::new(mu_x.clone()));
            let diff_y = Expr::Sub(Arc::new(y_i.clone()), Arc::new(mu_y.clone()));
            Expr::Mul(Arc::new(diff_x), Arc::new(diff_y))
        })
        .reduce(|acc, e| simplify(Expr::Add(Arc::new(acc), Arc::new(e))))
        .unwrap_or(Expr::Constant(0.0));
    simplify(Expr::Div(
        Arc::new(sum_of_products),
        Arc::new(Expr::Constant(n as f64)),
    ))
}

/// Computes the symbolic Pearson correlation coefficient.
///
/// The Pearson correlation coefficient `ρ` measures the linear correlation between
/// two sets of data. It is defined as `cov(X, Y) / (σ_x * σ_y)`,
/// where `cov` is the covariance and `σ` is the standard deviation.
///
/// # Arguments
/// * `data1` - The first set of data points.
/// * `data2` - The second set of data points.
///
/// # Returns
/// An `Expr` representing the symbolic Pearson correlation coefficient.
pub fn correlation(data1: &[Expr], data2: &[Expr]) -> Expr {
    let cov_xy = covariance(data1, data2);
    let std_dev_x = std_dev(data1);
    let std_dev_y = std_dev(data2);
    simplify(Expr::Div(
        Arc::new(cov_xy),
        Arc::new(Expr::Mul(Arc::new(std_dev_x), Arc::new(std_dev_y))),
    ))
}
