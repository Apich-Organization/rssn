//! # Symbolic Probability Distributions
//!
//! This module provides symbolic representations of common probability distributions,
//! both discrete and continuous. It includes structures for Normal, Uniform, Binomial,
//! Poisson, Bernoulli, Exponential, Gamma, Beta, and Student's t-distributions,
//! along with methods to generate their symbolic PDF/PMF, CDF, expectation, and variance.

use std::sync::Arc;

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
            Arc::new(one.clone()),
            Arc::new(Expr::Sqrt(Arc::new(Expr::Mul(
                Arc::new(two.clone()),
                Arc::new(pi),
            )))),
        );
        let term2 = Expr::Div(Arc::new(term1), Arc::new(self.std_dev.clone()));

        let exp_arg_num = Expr::Neg(Arc::new(Expr::Power(
            Arc::new(Expr::Sub(Arc::new(x.clone()), Arc::new(self.mean.clone()))),
            Arc::new(two.clone()),
        )));
        let exp_arg_den = Expr::Mul(
            Arc::new(two.clone()),
            Arc::new(Expr::Power(
                Arc::new(self.std_dev.clone()),
                Arc::new(two.clone()),
            )),
        );
        let exp_arg = Expr::Div(Arc::new(exp_arg_num), Arc::new(exp_arg_den));

        simplify(Expr::Mul(
            Arc::new(term2),
            Arc::new(Expr::Exp(Arc::new(exp_arg))),
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
            Arc::new(Expr::Sub(Arc::new(x.clone()), Arc::new(self.mean.clone()))),
            Arc::new(Expr::Mul(
                Arc::new(self.std_dev.clone()),
                Arc::new(Expr::Sqrt(Arc::new(two))),
            )),
        );
        simplify(Expr::Mul(
            Arc::new(Expr::Constant(0.5)),
            Arc::new(Expr::Add(Arc::new(one), Arc::new(Expr::Erf(Arc::new(arg))))),
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
            Arc::new(self.std_dev.clone()),
            Arc::new(Expr::Constant(2.0)),
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
            Arc::new(Expr::Constant(1.0)),
            Arc::new(Expr::Sub(
                Arc::new(self.max.clone()),
                Arc::new(self.min.clone()),
            )),
        ))
    }

    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Uniform distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `(min + max) / 2`.
        simplify(Expr::Div(
            Arc::new(Expr::Add(
                Arc::new(self.max.clone()),
                Arc::new(self.min.clone()),
            )),
            Arc::new(Expr::Constant(2.0)),
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
        let p_k = Expr::Power(Arc::new(self.p.clone()), Arc::new(k.clone()));
        let one_minus_p = Expr::Sub(Arc::new(Expr::Constant(1.0)), Arc::new(self.p.clone()));
        let n_minus_k = Expr::Sub(Arc::new(self.n.clone()), Arc::new(k.clone()));
        let one_minus_p_pow = Expr::Power(Arc::new(one_minus_p), Arc::new(n_minus_k));
        simplify(Expr::Mul(
            Arc::new(n_choose_k),
            Arc::new(Expr::Mul(Arc::new(p_k), Arc::new(one_minus_p_pow))),
        ))
    }

    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `n * p`.
        simplify(Expr::Mul(
            Arc::new(self.n.clone()),
            Arc::new(self.p.clone()),
        ))
    }

    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `n * p * (1 - p)`.
        let one_minus_p = Expr::Sub(Arc::new(Expr::Constant(1.0)), Arc::new(self.p.clone()));
        simplify(Expr::Mul(
            Arc::new(self.n.clone()),
            Arc::new(Expr::Mul(Arc::new(self.p.clone()), Arc::new(one_minus_p))),
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
        let lambda_k = Expr::Power(Arc::new(self.rate.clone()), Arc::new(k.clone()));
        let exp_neg_lambda = Expr::Exp(Arc::new(Expr::Neg(Arc::new(self.rate.clone()))));
        let k_factorial = Expr::Factorial(Arc::new(k.clone()));
        simplify(Expr::Div(
            Arc::new(Expr::Mul(Arc::new(lambda_k), Arc::new(exp_neg_lambda))),
            Arc::new(k_factorial),
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
        let one_minus_p = Expr::Sub(Arc::new(Expr::Constant(1.0)), Arc::new(self.p.clone()));
        let p_term = Expr::Mul(Arc::new(self.p.clone()), Arc::new(k.clone()));
        let one_minus_p_term = Expr::Mul(
            Arc::new(one_minus_p),
            Arc::new(Expr::Sub(
                Arc::new(Expr::Constant(1.0)),
                Arc::new(k.clone()),
            )),
        );
        simplify(Expr::Add(Arc::new(p_term), Arc::new(one_minus_p_term))) // This is a trick: k*p + (1-k)*(1-p)
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
            Arc::new(self.p.clone()),
            Arc::new(Expr::Sub(
                Arc::new(Expr::Constant(1.0)),
                Arc::new(self.p.clone()),
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
            Arc::new(self.rate.clone()),
            Arc::new(Expr::Exp(Arc::new(Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(self.rate.clone()),
                Arc::new(x.clone()),
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
            Arc::new(Expr::Constant(1.0)),
            Arc::new(Expr::Exp(Arc::new(Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(self.rate.clone()),
                Arc::new(x.clone()),
            )))))),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `1 / λ`.
        simplify(Expr::Div(
            Arc::new(Expr::Constant(1.0)),
            Arc::new(self.rate.clone()),
        ))
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `1 / λ²`.
        simplify(Expr::Div(
            Arc::new(Expr::Constant(1.0)),
            Arc::new(Expr::Power(
                Arc::new(self.rate.clone()),
                Arc::new(Expr::Constant(2.0)),
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
        let term1_num = Expr::Power(Arc::new(self.rate.clone()), Arc::new(self.shape.clone()));
        let term1_den = Expr::Gamma(Arc::new(self.shape.clone()));
        let term1 = Expr::Div(Arc::new(term1_num), Arc::new(term1_den));

        let term2 = Expr::Power(
            Arc::new(x.clone()),
            Arc::new(Expr::Sub(
                Arc::new(self.shape.clone()),
                Arc::new(Expr::Constant(1.0)),
            )),
        );
        let term3 = Expr::Exp(Arc::new(Expr::Neg(Arc::new(Expr::Mul(
            Arc::new(self.rate.clone()),
            Arc::new(x.clone()),
        )))));

        simplify(Expr::Mul(
            Arc::new(term1),
            Arc::new(Expr::Mul(Arc::new(term2), Arc::new(term3))),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Gamma distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `α / β`.
        simplify(Expr::Div(
            Arc::new(self.shape.clone()),
            Arc::new(self.rate.clone()),
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
            Arc::new(x.clone()),
            Arc::new(Expr::Sub(
                Arc::new(self.alpha.clone()),
                Arc::new(Expr::Constant(1.0)),
            )),
        );
        let one_minus_x = Expr::Sub(Arc::new(Expr::Constant(1.0)), Arc::new(x.clone()));
        let num2 = Expr::Power(
            Arc::new(one_minus_x),
            Arc::new(Expr::Sub(
                Arc::new(self.beta.clone()),
                Arc::new(Expr::Constant(1.0)),
            )),
        );
        let den = Expr::Beta(Arc::new(self.alpha.clone()), Arc::new(self.beta.clone()));
        simplify(Expr::Div(
            Arc::new(Expr::Mul(Arc::new(num1), Arc::new(num2))),
            Arc::new(den),
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
        let term1_num = Expr::Gamma(Arc::new(Expr::Div(
            Arc::new(Expr::Add(
                Arc::new(self.nu.clone()),
                Arc::new(Expr::Constant(1.0)),
            )),
            Arc::new(Expr::Constant(2.0)),
        )));
        let term1_den_sqrt = Expr::Sqrt(Arc::new(Expr::Mul(
            Arc::new(self.nu.clone()),
            Arc::new(Expr::Pi),
        )));
        let term1_den_gamma = Expr::Gamma(Arc::new(Expr::Div(
            Arc::new(self.nu.clone()),
            Arc::new(Expr::Constant(2.0)),
        )));
        let term1 = Expr::Div(
            Arc::new(term1_num),
            Arc::new(Expr::Mul(
                Arc::new(term1_den_sqrt),
                Arc::new(term1_den_gamma),
            )),
        );

        let term2_base = Expr::Add(
            Arc::new(Expr::Constant(1.0)),
            Arc::new(Expr::Div(
                Arc::new(Expr::Power(
                    Arc::new(t.clone()),
                    Arc::new(Expr::Constant(2.0)),
                )),
                Arc::new(self.nu.clone()),
            )),
        );
        let term2_exp = Expr::Neg(Arc::new(Expr::Div(
            Arc::new(Expr::Add(
                Arc::new(self.nu.clone()),
                Arc::new(Expr::Constant(1.0)),
            )),
            Arc::new(Expr::Constant(2.0)),
        )));
        let term2 = Expr::Power(Arc::new(term2_base), Arc::new(term2_exp));

        simplify(Expr::Mul(Arc::new(term1), Arc::new(term2)))
    }
}
