//! # Symbolic Probability Distributions
//!
//! This module provides symbolic representations of common probability distributions,
//! both discrete and continuous. It includes structures for Normal, Uniform, Binomial,
//! Poisson, Bernoulli, Exponential, Gamma, Beta, and Student's t-distributions,
//! along with methods to generate their symbolic PDF/PMF, CDF, expectation, and variance.
use crate::symbolic::combinatorics::combinations;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use std::f64::consts::PI;
use std::sync::Arc;
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
        let term1 = Expr::new_div(one, Expr::new_sqrt(Expr::new_mul(two.clone(), pi)));
        let term2 = Expr::new_div(term1, self.std_dev.clone());
        let exp_arg_num = Expr::new_neg(Expr::new_pow(
            Expr::new_sub(x.clone(), self.mean.clone()),
            two.clone(),
        ));
        let exp_arg_den = Expr::new_mul(
            two.clone(),
            Expr::new_pow(self.std_dev.clone(), two),
        );
        let exp_arg = Expr::new_div(exp_arg_num, exp_arg_den);
        simplify(&Expr::new_mul(term2, Expr::new_exp(exp_arg)))
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
        let arg = Expr::new_div(
            Expr::new_sub(x.clone(), self.mean.clone()),
            Expr::new_mul(self.std_dev.clone(), Expr::new_sqrt(two)),
        );
        simplify(&Expr::new_mul(
            Expr::Constant(0.5),
            Expr::new_add(one, Expr::new_erf(arg)),
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
        simplify(&Expr::new_pow(self.std_dev.clone(), Expr::Constant(2.0)))
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
        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            Expr::new_sub(self.max.clone(), self.min.clone()),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Uniform distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `(min + max) / 2`.
        simplify(&Expr::new_div(
            Expr::new_add(self.max.clone(), self.min.clone()),
            Expr::Constant(2.0),
        ))
    }
}
/// Represents a Binomial distribution with symbolic parameters.
pub struct Binomial {
    pub n: Expr,
    pub p: Expr,
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
        let n_choose_k = combinations(&self.n.clone(), k.clone());
        let p_k = Expr::new_pow(self.p.clone(), k.clone());
        let one_minus_p = Expr::new_sub(Expr::Constant(1.0), self.p.clone());
        let n_minus_k = Expr::new_sub(self.n.clone(), k.clone());
        let one_minus_p_pow = Expr::new_pow(one_minus_p, n_minus_k);
        simplify(&Expr::new_mul(
            n_choose_k,
            Expr::new_mul(p_k, one_minus_p_pow),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `n * p`.
        simplify(&Expr::new_mul(self.n.clone(), self.p.clone()))
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Binomial distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `n * p * (1 - p)`.
        let one_minus_p = Expr::new_sub(Expr::Constant(1.0), self.p.clone());
        simplify(&Expr::new_mul(
            self.n.clone(),
            Expr::new_mul(self.p.clone(), one_minus_p),
        ))
    }
}
/// Represents a Poisson distribution with symbolic rate parameter λ.
pub struct Poisson {
    pub rate: Expr,
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
        let lambda_k = Expr::new_pow(self.rate.clone(), k.clone());
        let exp_neg_lambda = Expr::new_exp(Expr::new_neg(self.rate.clone()));
        let k_factorial = Expr::Factorial(Arc::new(k.clone()));
        simplify(&Expr::new_div(
            Expr::new_mul(lambda_k, exp_neg_lambda),
            k_factorial,
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
        let one_minus_p = Expr::new_sub(Expr::Constant(1.0), self.p.clone());
        let p_term = Expr::new_mul(self.p.clone(), k.clone());
        let one_minus_p_term =
            Expr::new_mul(one_minus_p, Expr::new_sub(Expr::Constant(1.0), k.clone()));
        simplify(&Expr::new_add(p_term, one_minus_p_term))
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
        simplify(&Expr::new_mul(
            self.p.clone(),
            Expr::new_sub(Expr::Constant(1.0), self.p.clone()),
        ))
    }
}
/// Represents an Exponential distribution with symbolic rate λ.
pub struct Exponential {
    pub rate: Expr,
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
        simplify(&Expr::new_mul(
            self.rate.clone(),
            Expr::new_exp(Expr::new_neg(Expr::new_mul(self.rate.clone(), x.clone()))),
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
        simplify(&Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_exp(Expr::new_neg(Expr::new_mul(self.rate.clone(), x.clone()))),
        ))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `1 / λ`.
        simplify(&Expr::new_div(Expr::Constant(1.0), self.rate.clone()))
    }
    pub fn variance(&self) -> Expr {
        /// Returns the symbolic variance of the Exponential distribution.
        ///
        /// # Returns
        /// An `Expr` representing the variance `1 / λ²`.
        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            Expr::new_pow(self.rate.clone(), Expr::Constant(2.0)),
        ))
    }
}
/// Represents a Gamma distribution with symbolic shape α and rate β.
pub struct Gamma {
    pub shape: Expr,
    pub rate: Expr,
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
        let term1_num = Expr::new_pow(self.rate.clone(), self.shape.clone());
        let term1_den = Expr::new_gamma(self.shape.clone());
        let term1 = Expr::new_div(term1_num, term1_den);
        let term2 = Expr::new_pow(
            x.clone(),
            Expr::new_sub(self.shape.clone(), Expr::Constant(1.0)),
        );
        let term3 = Expr::new_exp(Expr::new_neg(Expr::new_mul(self.rate.clone(), x.clone())));
        simplify(&Expr::new_mul(term1, Expr::new_mul(term2, term3)))
    }
    pub fn expectation(&self) -> Expr {
        /// Returns the symbolic expectation (mean) of the Gamma distribution.
        ///
        /// # Returns
        /// An `Expr` representing the mean `α / β`.
        simplify(&Expr::new_div(self.shape.clone(), self.rate.clone()))
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
        let num1 = Expr::new_pow(
            x.clone(),
            Expr::new_sub(self.alpha.clone(), Expr::Constant(1.0)),
        );
        let one_minus_x = Expr::new_sub(Expr::Constant(1.0), x.clone());
        let num2 = Expr::new_pow(
            one_minus_x,
            Expr::new_sub(self.beta.clone(), Expr::Constant(1.0)),
        );
        let den = Expr::new_beta(self.alpha.clone(), self.beta.clone());
        simplify(&Expr::new_div(Expr::new_mul(num1, num2), den))
    }
}
/// Represents a Student's t-distribution with symbolic degrees of freedom ν.
pub struct StudentT {
    pub nu: Expr,
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
        let term1_num = Expr::new_gamma(Expr::new_div(
            Expr::new_add(self.nu.clone(), Expr::Constant(1.0)),
            Expr::Constant(2.0),
        ));
        let term1_den_sqrt = Expr::new_sqrt(Expr::new_mul(self.nu.clone(), Expr::Pi));
        let term1_den_gamma = Expr::new_gamma(Expr::new_div(self.nu.clone(), Expr::Constant(2.0)));
        let term1 = Expr::new_div(term1_num, Expr::new_mul(term1_den_sqrt, term1_den_gamma));
        let term2_base = Expr::new_add(
            Expr::Constant(1.0),
            Expr::new_div(
                Expr::new_pow(t.clone(), Expr::Constant(2.0)),
                self.nu.clone(),
            ),
        );
        let term2_exp = Expr::new_neg(Expr::new_div(
            Expr::new_add(self.nu.clone(), Expr::Constant(1.0)),
            Expr::Constant(2.0),
        ));
        let term2 = Expr::new_pow(term2_base, term2_exp);
        simplify(&Expr::new_mul(term1, term2))
    }
}
