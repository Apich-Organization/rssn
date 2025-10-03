//! # Symbolic Information Theory
//!
//! This module provides functions for symbolic information theory calculations.
//! It includes measures such as Shannon entropy, Kullback-Leibler divergence,
//! cross-entropy, joint entropy, conditional entropy, mutual information, and Gini impurity,
//! all expressed symbolically for discrete probability distributions.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;

/// Computes the symbolic Shannon entropy of a discrete probability distribution.
///
/// Shannon entropy `H(X)` quantifies the average amount of information or uncertainty
/// in a random variable. For a discrete distribution `p(x)`, it is defined as
/// `H(X) = -Σ p(x) * log2(p(x))`.
///
/// # Arguments
/// * `probs` - A slice of `Expr` representing the probabilities `p(x)`.
///
/// # Returns
/// An `Expr` representing the symbolic Shannon entropy.
pub fn shannon_entropy(probs: &[Expr]) -> Expr {
    let log2 = Expr::Log(Box::new(Expr::Constant(2.0)));
    let sum = probs
        .iter()
        .map(|p| {
            let log2_p = Expr::Div(
                Box::new(Expr::Log(Box::new(p.clone()))),
                Box::new(log2.clone()),
            );
            Expr::Mul(Box::new(p.clone()), Box::new(log2_p))
        })
        .reduce(|acc, e| simplify(Expr::Add(Box::new(acc), Box::new(e))))
        .unwrap_or(Expr::Constant(0.0));

    simplify(Expr::Neg(Box::new(sum)))
}

/// Computes the symbolic Kullback-Leibler (KL) divergence between two distributions.
///
/// KL divergence `D_KL(P || Q)` measures how one probability distribution `Q` is different
/// from a second, reference probability distribution `P`. It is defined as
/// `D_KL(P || Q) = Σ p(x) * log2(p(x) / q(x))`.
///
/// # Arguments
/// * `p_dist` - A slice of `Expr` representing the reference probabilities `p(x)`.
/// * `q_dist` - A slice of `Expr` representing the comparison probabilities `q(x)`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the symbolic KL divergence, or an error
/// if the distributions have different lengths.
pub fn kl_divergence(p_dist: &[Expr], q_dist: &[Expr]) -> Result<Expr, String> {
    if p_dist.len() != q_dist.len() {
        return Err("Distributions must have the same length".to_string());
    }
    let log2 = Expr::Log(Box::new(Expr::Constant(2.0)));
    let sum = p_dist
        .iter()
        .zip(q_dist.iter())
        .map(|(p, q)| {
            let ratio = Expr::Div(Box::new(p.clone()), Box::new(q.clone()));
            let log2_ratio =
                Expr::Div(Box::new(Expr::Log(Box::new(ratio))), Box::new(log2.clone()));
            Expr::Mul(Box::new(p.clone()), Box::new(log2_ratio))
        })
        .reduce(|acc, e| simplify(Expr::Add(Box::new(acc), Box::new(e))))
        .unwrap_or(Expr::Constant(0.0));

    Ok(simplify(sum))
}

/// Computes the symbolic cross-entropy between two distributions.
///
/// Cross-entropy `H(P, Q)` measures the average number of bits needed to identify an event
/// drawn from `P` if a coding scheme optimized for `Q` is used. It is defined as
/// `H(P, Q) = -Σ p(x) * log2(q(x))`.
///
/// # Arguments
/// * `p_dist` - A slice of `Expr` representing the true probabilities `p(x)`.
/// * `q_dist` - A slice of `Expr` representing the predicted probabilities `q(x)`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the symbolic cross-entropy, or an error
/// if the distributions have different lengths.
pub fn cross_entropy(p_dist: &[Expr], q_dist: &[Expr]) -> Result<Expr, String> {
    if p_dist.len() != q_dist.len() {
        return Err("Distributions must have the same length".to_string());
    }
    let log2 = Expr::Log(Box::new(Expr::Constant(2.0)));
    let sum = p_dist
        .iter()
        .zip(q_dist.iter())
        .map(|(p, q)| {
            let log2_q = Expr::Div(
                Box::new(Expr::Log(Box::new(q.clone()))),
                Box::new(log2.clone()),
            );
            Expr::Mul(Box::new(p.clone()), Box::new(log2_q))
        })
        .reduce(|acc, e| simplify(Expr::Add(Box::new(acc), Box::new(e))))
        .unwrap_or(Expr::Constant(0.0));

    Ok(simplify(Expr::Neg(Box::new(sum))))
}

/// Computes the symbolic Joint Entropy of two discrete probability distributions.
///
/// Joint entropy `H(X, Y)` quantifies the uncertainty associated with a set of random
/// variables. For a joint distribution `p(x,y)`, it is defined as
/// `H(X, Y) = -Σ_x Σ_y p(x,y) * log2(p(x,y))`.
///
/// # Arguments
/// * `joint_probs` - An `Expr::Matrix` representing the joint probabilities `p(x,y)`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the symbolic joint entropy, or an error
/// if the input is not a matrix.
pub fn joint_entropy(joint_probs: &Expr) -> Result<Expr, String> {
    if let Expr::Matrix(rows) = joint_probs {
        let flat_probs: Vec<Expr> = rows.iter().flatten().cloned().collect();
        Ok(shannon_entropy(&flat_probs))
    } else {
        Err("Input must be a matrix of joint probabilities.".to_string())
    }
}

/// Computes the symbolic Conditional Entropy `H(Y|X)`.
///
/// Conditional entropy `H(Y|X)` measures the uncertainty of a random variable `Y`
/// given that the value of another random variable `X` is known. It is defined as
/// `H(Y|X) = H(X, Y) - H(X)`.
///
/// # Arguments
/// * `joint_probs` - An `Expr::Matrix` representing the joint probabilities `p(x,y)`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the symbolic conditional entropy, or an error
/// if the input is not a matrix.
pub fn conditional_entropy(joint_probs: &Expr) -> Result<Expr, String> {
    if let Expr::Matrix(rows) = joint_probs {
        let p_x: Vec<Expr> = rows
            .iter()
            .map(|row| {
                row.iter()
                    .cloned()
                    .reduce(|a, b| simplify(Expr::Add(Box::new(a), Box::new(b))))
                    .unwrap_or(Expr::Constant(0.0))
            })
            .collect();

        let h_x = shannon_entropy(&p_x);
        let h_xy = joint_entropy(joint_probs)?;

        Ok(simplify(Expr::Sub(Box::new(h_xy), Box::new(h_x))))
    } else {
        Err("Input must be a matrix of joint probabilities.".to_string())
    }
}

/// Computes the symbolic Mutual Information `I(X;Y)`.
///
/// Mutual information `I(X;Y)` measures the amount of information obtained about one
/// random variable by observing another random variable. It is defined as
/// `I(X;Y) = H(X) + H(Y) - H(X,Y)`.
///
/// # Arguments
/// * `joint_probs` - An `Expr::Matrix` representing the joint probabilities `p(x,y)`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the symbolic mutual information, or an error
/// if the input is not a matrix.
pub fn mutual_information(joint_probs: &Expr) -> Result<Expr, String> {
    if let Expr::Matrix(rows) = joint_probs {
        let p_x: Vec<Expr> = rows
            .iter()
            .map(|row| {
                row.iter()
                    .cloned()
                    .reduce(|a, b| simplify(Expr::Add(Box::new(a), Box::new(b))))
                    .unwrap_or(Expr::Constant(0.0))
            })
            .collect();

        let num_cols = rows.get(0).map_or(0, |r| r.len());
        let mut p_y = vec![Expr::Constant(0.0); num_cols];
        for row in rows {
            for (j, p_ij) in row.iter().enumerate() {
                p_y[j] = simplify(Expr::Add(Box::new(p_y[j].clone()), Box::new(p_ij.clone())));
            }
        }

        let h_x = shannon_entropy(&p_x);
        let h_y = shannon_entropy(&p_y);
        let h_xy = joint_entropy(joint_probs)?;

        Ok(simplify(Expr::Sub(
            Box::new(Expr::Add(Box::new(h_x), Box::new(h_y))),
            Box::new(h_xy),
        )))
    } else {
        Err("Input must be a matrix of joint probabilities.".to_string())
    }
}

/// Computes the symbolic Gini Impurity of a discrete probability distribution.
///
/// Gini impurity `G` is a measure of the impurity or disorder of a set of elements.
/// For a discrete distribution `p_i`, it is defined as `G = Σ p_i * (1 - p_i) = 1 - Σ p_i^2`.
/// It is commonly used in decision tree algorithms.
///
/// # Arguments
/// * `probs` - A slice of `Expr` representing the probabilities `p_i`.
///
/// # Returns
/// An `Expr` representing the symbolic Gini impurity.
pub fn gini_impurity(probs: &[Expr]) -> Expr {
    let sum_of_squares = probs
        .iter()
        .map(|p| Expr::Power(Box::new(p.clone()), Box::new(Expr::Constant(2.0))))
        .reduce(|acc, e| simplify(Expr::Add(Box::new(acc), Box::new(e))))
        .unwrap_or(Expr::Constant(0.0));

    simplify(Expr::Sub(
        Box::new(Expr::Constant(1.0)),
        Box::new(sum_of_squares),
    ))
}
