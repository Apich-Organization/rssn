//! # Symbolic Probability Distributions
//!
//! This module provides symbolic representations of common probability distributions,
//! both discrete and continuous. It includes structures for Normal, Uniform, Binomial,
//! Poisson, Bernoulli, Exponential, Gamma, Beta, and Student's t-distributions,
//! along with methods to generate their symbolic PDF/PMF, CDF, expectation, and variance.

use crate::symbolic::combinatorics::combinations;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use std::f64::consts::PI;

/// Represents a Normal (Gaussian) distribution with symbolic parameters.
pub struct Normal {
    pub mean: Expr,
    pub std_dev: Expr,
}

impl Normal {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of a normal distribution is given by: `f(x) = (1 / (σ * sqrt(2π))) * exp(- (x - μ)² / (2σ²))`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, x: &Expr) -> Expr {
        let pi = Expr::Constant(PI);
        let two = Expr::Constant(2.0);
        let one = Expr::Constant(1.0);

        let term1 = Expr::Div(
            Box::new(one.clone()),
            Box::new(Expr::Sqrt(Box::new(Expr::Mul(
                Box::new(two.clone()),
                Box::new(pi),
            )))),
        );
        let term2 = Expr::Div(Box::new(term1), Box::new(self.std_dev.clone()));

        let exp_arg_num = Expr::Neg(Box::new(Expr::Power(
            Box::new(Expr::Sub(Box::new(x.clone()), Box::new(self.mean.clone()))),
            Box::new(two.clone()),
        )));
        let exp_arg_den = Expr::Mul(
            Box::new(two.clone()),
            Box::new(Expr::Power(
                Box::new(self.std_dev.clone()),
                Box::new(two.clone()),
            )),
        );
        let exp_arg = Expr::Div(Box::new(exp_arg_num), Box::new(exp_arg_den));

        simplify(Expr::Mul(
            Box::new(term2),
            Box::new(Expr::Exp(Box::new(exp_arg))),
        ))
    }

    /// Returns the symbolic expression for the cumulative distribution function (CDF).
    ///
    /// The CDF of a normal distribution is given by: `F(x) = 0.5 * (1 + erf((x - μ) / (σ * sqrt(2))))`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the CDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic CDF.
    pub fn cdf(&self, x: &Expr) -> Expr {
        let one = Expr::Constant(1.0);
        let two = Expr::Constant(2.0);
        let arg = Expr::Div(
            Box::new(Expr::Sub(Box::new(x.clone()), Box::new(self.mean.clone()))),
            Box::new(Expr::Mul(
                Box::new(self.std_dev.clone()),
                Box::new(Expr::Sqrt(Box::new(two))),
            )),
        );
        simplify(Expr::Mul(
            Box::new(Expr::Constant(0.5)),
            Box::new(Expr::Add(Box::new(one), Box::new(Expr::Erf(Box::new(arg))))),
        ))
    }

    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Normal distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `μ`.
        self.mean.clone()
    }

    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Normal distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `σ²`.
        simplify(Expr::Power(
            Box::new(self.std_dev.clone()),
            Box::new(Expr::Constant(2.0)),
        ))
    }
}

/// Represents a Uniform distribution with symbolic parameters.
pub struct Uniform {
    pub min: Expr,
    pub max: Expr,
}

impl Uniform {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of a uniform distribution over `[min, max]` is `1 / (max - min)` for `min <= x <= max`,
    /// and `0` otherwise.
    ///
    /// # Arguments
    /// * `_x` - The value at which to evaluate the PDF (ignored for the constant value within range).
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, _x: &Expr) -> Expr {
        // A full implementation would return a piecewise function.
        // For now, we return the value within the range.
        simplify(Expr::Div(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Sub(
                Box::new(self.max.clone()),
                Box::new(self.min.clone()),
            )),
        ))
    }

    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Uniform distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `(min + max) / 2`.
        simplify(Expr::Div(
            Box::new(Expr::Add(
                Box::new(self.max.clone()),
                Box::new(self.min.clone()),
            )),
            Box::new(Expr::Constant(2.0)),
        ))
    }
}

/// Represents a Binomial distribution with symbolic parameters.
pub struct Binomial {
    pub n: Expr, // number of trials
    pub p: Expr, // probability of success
}

impl Binomial {
    /// Returns the symbolic expression for the probability mass function (PMF).
    ///
    /// The PMF of a binomial distribution is given by: `P(X=k) = C(n, k) * p^k * (1-p)^(n-k)`.
    ///
    /// # Arguments
    /// * `k` - The number of successes.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PMF.
    pub fn pmf(&self, k: &Expr) -> Expr {
        let n_choose_k = combinations(self.n.clone(), k.clone());
        let p_k = Expr::Power(Box::new(self.p.clone()), Box::new(k.clone()));
        let one_minus_p = Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(self.p.clone()));
        let n_minus_k = Expr::Sub(Box::new(self.n.clone()), Box::new(k.clone()));
        let one_minus_p_pow = Expr::Power(Box::new(one_minus_p), Box::new(n_minus_k));
        simplify(Expr::Mul(
            Box::new(n_choose_k),
            Box::new(Expr::Mul(Box::new(p_k), Box::new(one_minus_p_pow))),
        ))
    }

    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `n * p`.
        simplify(Expr::Mul(
            Box::new(self.n.clone()),
            Box::new(self.p.clone()),
        ))
    }

    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `n * p * (1 - p)`.
        let one_minus_p = Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(self.p.clone()));
        simplify(Expr::Mul(
            Box::new(self.n.clone()),
            Box::new(Expr::Mul(Box::new(self.p.clone()), Box::new(one_minus_p))),
        ))
    }
}

// =====================================================================================
// region: More Discrete Distributions
// =====================================================================================

/// Represents a Poisson distribution with symbolic rate parameter λ.
pub struct Poisson {
    pub rate: Expr, // lambda
}

impl Poisson {
    /// Returns the symbolic expression for the probability mass function (PMF).
    ///
    /// The PMF of a Poisson distribution is given by: `P(X=k) = (λ^k * e^(-λ)) / k!`.
    ///
    /// # Arguments
    /// * `k` - The number of occurrences.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PMF.
    pub fn pmf(&self, k: &Expr) -> Expr {
        let lambda_k = Expr::Power(Box::new(self.rate.clone()), Box::new(k.clone()));
        let exp_neg_lambda = Expr::Exp(Box::new(Expr::Neg(Box::new(self.rate.clone()))));
        let k_factorial = Expr::Factorial(Box::new(k.clone()));
        simplify(Expr::Div(
            Box::new(Expr::Mul(Box::new(lambda_k), Box::new(exp_neg_lambda))),
            Box::new(k_factorial),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Poisson distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `λ`.
        self.rate.clone()
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Poisson distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `λ`.
        self.rate.clone()
    }
}

/// Represents a Bernoulli distribution with symbolic probability p.
pub struct Bernoulli {
    pub p: Expr,
}

impl Bernoulli {
    /// Returns the symbolic expression for the probability mass function (PMF).
    ///
    /// The PMF of a Bernoulli distribution is `p` for `k=1` (success) and `1-p` for `k=0` (failure).
    ///
    /// # Arguments
    /// * `k` - The outcome (0 or 1).
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PMF.
    pub fn pmf(&self, k: &Expr) -> Expr {
        // p if k=1, 1-p if k=0
        let one_minus_p = Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(self.p.clone()));
        let p_term = Expr::Mul(Box::new(self.p.clone()), Box::new(k.clone()));
        let one_minus_p_term = Expr::Mul(
            Box::new(one_minus_p),
            Box::new(Expr::Sub(
                Box::new(Expr::Constant(1.0)),
                Box::new(k.clone()),
            )),
        );
        simplify(Expr::Add(Box::new(p_term), Box::new(one_minus_p_term))) // This is a trick: k*p + (1-k)*(1-p)
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Bernoulli distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `p`.
        self.p.clone()
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Bernoulli distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `p * (1 - p)`.
        simplify(Expr::Mul(
            Box::new(self.p.clone()),
            Box::new(Expr::Sub(
                Box::new(Expr::Constant(1.0)),
                Box::new(self.p.clone()),
            )),
        ))
    }
}

// =====================================================================================
// region: More Continuous Distributions
// =====================================================================================

/// Represents an Exponential distribution with symbolic rate λ.
pub struct Exponential {
    pub rate: Expr, // lambda
}

impl Exponential {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of an exponential distribution is `f(x) = λ * e^(-λx)` for `x >= 0`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, x: &Expr) -> Expr {
        simplify(Expr::Mul(
            Box::new(self.rate.clone()),
            Box::new(Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
                Box::new(self.rate.clone()),
                Box::new(x.clone()),
            )))))),
        ))
    }
    /// Returns the symbolic expression for the cumulative distribution function (CDF).
    ///
    /// The CDF of an exponential distribution is `F(x) = 1 - e^(-λx)` for `x >= 0`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the CDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic CDF.
    pub fn cdf(&self, x: &Expr) -> Expr {
        simplify(Expr::Sub(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
                Box::new(self.rate.clone()),
                Box::new(x.clone()),
            )))))),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `1 / λ`.
        simplify(Expr::Div(
            Box::new(Expr::Constant(1.0)),
            Box::new(self.rate.clone()),
        ))
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `1 / λ²`.
        simplify(Expr::Div(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Power(
                Box::new(self.rate.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
        ))
    }
}

/// Represents a Gamma distribution with symbolic shape α and rate β.
pub struct Gamma {
    pub shape: Expr, // alpha
    pub rate: Expr,  // beta
}

impl Gamma {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of a Gamma distribution is `f(x; α, β) = (β^α / Γ(α)) * x^(α-1) * e^(-βx)`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, x: &Expr) -> Expr {
        let term1_num = Expr::Power(Box::new(self.rate.clone()), Box::new(self.shape.clone()));
        let term1_den = Expr::Gamma(Box::new(self.shape.clone()));
        let term1 = Expr::Div(Box::new(term1_num), Box::new(term1_den));

        let term2 = Expr::Power(
            Box::new(x.clone()),
            Box::new(Expr::Sub(
                Box::new(self.shape.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
        );
        let term3 = Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Mul(
            Box::new(self.rate.clone()),
            Box::new(x.clone()),
        )))));

        simplify(Expr::Mul(
            Box::new(term1),
            Box::new(Expr::Mul(Box::new(term2), Box::new(term3))),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Gamma distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `α / β`.
        simplify(Expr::Div(
            Box::new(self.shape.clone()),
            Box::new(self.rate.clone()),
        ))
    }
}

/// Represents a Beta distribution with symbolic parameters α and β.
pub struct Beta {
    pub alpha: Expr,
    pub beta: Expr,
}

impl Beta {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of a Beta distribution is `f(x; α, β) = (1 / B(α, β)) * x^(α-1) * (1-x)^(β-1)`.
    ///
    /// # Arguments
    /// * `x` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, x: &Expr) -> Expr {
        let num1 = Expr::Power(
            Box::new(x.clone()),
            Box::new(Expr::Sub(
                Box::new(self.alpha.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
        );
        let one_minus_x = Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(x.clone()));
        let num2 = Expr::Power(
            Box::new(one_minus_x),
            Box::new(Expr::Sub(
                Box::new(self.beta.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
        );
        let den = Expr::Beta(Box::new(self.alpha.clone()), Box::new(self.beta.clone()));
        simplify(Expr::Div(
            Box::new(Expr::Mul(Box::new(num1), Box::new(num2))),
            Box::new(den),
        ))
    }
}

/// Represents a Student's t-distribution with symbolic degrees of freedom ν.
pub struct StudentT {
    pub nu: Expr, // degrees of freedom
}

impl StudentT {
    /// Returns the symbolic expression for the probability density function (PDF).
    ///
    /// The PDF of a Student's t-distribution is `f(t; ν) = (Γ((ν+1)/2) / (sqrt(νπ) * Γ(ν/2))) * (1 + t²/ν)^(-(ν+1)/2)`.
    ///
    /// # Arguments
    /// * `t` - The value at which to evaluate the PDF.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic PDF.
    pub fn pdf(&self, t: &Expr) -> Expr {
        let term1_num = Expr::Gamma(Box::new(Expr::Div(
            Box::new(Expr::Add(
                Box::new(self.nu.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
            Box::new(Expr::Constant(2.0)),
        )));
        let term1_den_sqrt = Expr::Sqrt(Box::new(Expr::Mul(
            Box::new(self.nu.clone()),
            Box::new(Expr::Pi),
        )));
        let term1_den_gamma = Expr::Gamma(Box::new(Expr::Div(
            Box::new(self.nu.clone()),
            Box::new(Expr::Constant(2.0)),
        )));
        let term1 = Expr::Div(
            Box::new(term1_num),
            Box::new(Expr::Mul(
                Box::new(term1_den_sqrt),
                Box::new(term1_den_gamma),
            )),
        );

        let term2_base = Expr::Add(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Div(
                Box::new(Expr::Power(
                    Box::new(t.clone()),
                    Box::new(Expr::Constant(2.0)),
                )),
                Box::new(self.nu.clone()),
            )),
        );
        let term2_exp = Expr::Neg(Box::new(Expr::Div(
            Box::new(Expr::Add(
                Box::new(self.nu.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
            Box::new(Expr::Constant(2.0)),
        )));
        let term2 = Expr::Power(Box::new(term2_base), Box::new(term2_exp));

        simplify(Expr::Mul(Box::new(term1), Box::new(term2)))
    }
}
