//! # Symbolic Statistics
//!
//! This module provides functions for symbolic statistical calculations.
//! It includes basic descriptive statistics such as mean, variance, standard deviation,
//! covariance, and correlation, all expressed symbolically.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

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
#[must_use]

pub fn mean(data: &[Expr]) -> Expr {

    let n = data.len();

    if n == 0 {

        return Expr::new_constant(0.0);
    }

    let sum = data
        .iter()
        .cloned()
        .reduce(|acc, e| {

            simplify(&Expr::new_add(
                acc, e,
            ))
        })
        .unwrap_or(Expr::new_constant(0.0));

    simplify(&Expr::new_div(
        sum,
        Expr::new_constant(n as f64),
    ))
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
#[must_use]

pub fn variance(data: &[Expr]) -> Expr {

    let n = data.len();

    if n == 0 {

        return Expr::new_constant(0.0);
    }

    let mu = mean(data);

    let squared_diffs = data
        .iter()
        .map(|x_i| {

            let diff = Expr::new_sub(
                x_i.clone(),
                mu.clone(),
            );

            Expr::new_pow(
                diff,
                Expr::new_constant(2.0),
            )
        })
        .reduce(|acc, e| {

            simplify(&Expr::new_add(
                acc, e,
            ))
        })
        .unwrap_or(Expr::new_constant(0.0));

    simplify(&Expr::new_div(
        squared_diffs,
        Expr::new_constant(n as f64),
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
#[must_use]

pub fn std_dev(data: &[Expr]) -> Expr {

    simplify(&Expr::new_sqrt(
        variance(data),
    ))
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
#[must_use]

pub fn covariance(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {

    if data1.len() != data2.len()
        || data1.is_empty()
    {

        return Expr::new_constant(0.0);
    }

    let n = data1.len();

    let mu_x = mean(data1);

    let mu_y = mean(data2);

    let sum_of_products = data1
        .iter()
        .zip(data2.iter())
        .map(|(x_i, y_i)| {

            let diff_x = Expr::new_sub(
                x_i.clone(),
                mu_x.clone(),
            );

            let diff_y = Expr::new_sub(
                y_i.clone(),
                mu_y.clone(),
            );

            Expr::new_mul(
                diff_x,
                diff_y,
            )
        })
        .reduce(|acc, e| {

            simplify(&Expr::new_add(
                acc, e,
            ))
        })
        .unwrap_or(Expr::new_constant(0.0));

    simplify(&Expr::new_div(
        sum_of_products,
        Expr::new_constant(n as f64),
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
#[must_use]

pub fn correlation(
    data1: &[Expr],
    data2: &[Expr],
) -> Expr {

    let cov_xy =
        covariance(data1, data2);

    let std_dev_x = std_dev(data1);

    let std_dev_y = std_dev(data2);

    simplify(&Expr::new_div(
        cov_xy,
        Expr::new_mul(
            std_dev_x,
            std_dev_y,
        ),
    ))
}
