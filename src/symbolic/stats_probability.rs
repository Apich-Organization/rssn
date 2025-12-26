//! # Symbolic Probability Distributions
//!
//! This module provides symbolic representations of common probability distributions,
//! both discrete and continuous. It includes structures for Normal, Uniform, Binomial,
//! Poisson, Bernoulli, Exponential, Gamma, Beta, and Student's t-distributions,
//! along with methods to generate their symbolic PDF/PMF, CDF, expectation, and variance.

use crate::symbolic::combinatorics::combinations;
use crate::symbolic::core::{
    Distribution,
    Expr,
};
use crate::symbolic::simplify_dag::simplify;
use std::f64::consts::PI;
use std::fmt::Debug;
use std::sync::Arc;

/// Represents a Normal (Gaussian) distribution with symbolic parameters.
#[derive(Debug, Clone)]

pub struct Normal {
    pub mean: Expr,
    pub std_dev: Expr,
}

impl Distribution for Normal {
    fn pdf(
        &self,
        x: &Expr,
    ) -> Expr {

        let pi = Expr::Constant(PI);

        let two = Expr::Constant(2.0);

        let one = Expr::Constant(1.0);

        let term1 = Expr::new_div(
            one,
            Expr::new_sqrt(Expr::new_mul(
                two.clone(),
                pi,
            )),
        );

        let term2 = Expr::new_div(
            term1,
            self.std_dev.clone(),
        );

        let exp_arg_num = Expr::new_neg(Expr::new_pow(
            Expr::new_sub(
                x.clone(),
                self.mean.clone(),
            ),
            two.clone(),
        ));

        let exp_arg_den = Expr::new_mul(
            two.clone(),
            Expr::new_pow(
                self.std_dev.clone(),
                two,
            ),
        );

        let exp_arg = Expr::new_div(
            exp_arg_num,
            exp_arg_den,
        );

        simplify(&Expr::new_mul(
            term2,
            Expr::new_exp(exp_arg),
        ))
    }

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr {

        let one = Expr::Constant(1.0);

        let two = Expr::Constant(2.0);

        let arg = Expr::new_div(
            Expr::new_sub(
                x.clone(),
                self.mean.clone(),
            ),
            Expr::new_mul(
                self.std_dev.clone(),
                Expr::new_sqrt(two),
            ),
        );

        simplify(&Expr::new_mul(
            Expr::Constant(0.5),
            Expr::new_add(
                one,
                Expr::new_erf(arg),
            ),
        ))
    }

    fn expectation(&self) -> Expr {

        self.mean.clone()
    }

    fn variance(&self) -> Expr {

        simplify(&Expr::new_pow(
            self.std_dev.clone(),
            Expr::Constant(2.0),
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // exp(mu*t + \sigma^2*t^2 / 2)
        let mu_t = Expr::new_mul(
            self.mean.clone(),
            t.clone(),
        );

        let sigma_sq = Expr::new_pow(
            self.std_dev.clone(),
            Expr::Constant(2.0),
        );

        let t_sq = Expr::new_pow(
            t.clone(),
            Expr::Constant(2.0),
        );

        let half_sigma_sq_t_sq = Expr::new_div(
            Expr::new_mul(sigma_sq, t_sq),
            Expr::Constant(2.0),
        );

        simplify(&Expr::new_exp(
            Expr::new_add(
                mu_t,
                half_sigma_sq_t_sq,
            ),
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Uniform distribution with symbolic parameters.
#[derive(Debug, Clone)]

pub struct Uniform {
    pub min: Expr,
    pub max: Expr,
}

impl Distribution for Uniform {
    fn pdf(
        &self,
        _x: &Expr,
    ) -> Expr {

        // Technically this should be piece-wise, but for symbolic simplification we often return the density inside the support
        // Or we could return a Piecewise if supported. For now, sticking to the main density.
        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            Expr::new_sub(
                self.max.clone(),
                self.min.clone(),
            ),
        ))
    }

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr {

        // (x - min) / (max - min)
        let num = Expr::new_sub(
            x.clone(),
            self.min.clone(),
        );

        let den = Expr::new_sub(
            self.max.clone(),
            self.min.clone(),
        );

        simplify(&Expr::new_div(
            num, den,
        ))
    }

    fn expectation(&self) -> Expr {

        simplify(&Expr::new_div(
            Expr::new_add(
                self.max.clone(),
                self.min.clone(),
            ),
            Expr::Constant(2.0),
        ))
    }

    fn variance(&self) -> Expr {

        // (max - min)^2 / 12
        let range = Expr::new_sub(
            self.max.clone(),
            self.min.clone(),
        );

        let num = Expr::new_pow(
            range,
            Expr::Constant(2.0),
        );

        simplify(&Expr::new_div(
            num,
            Expr::Constant(12.0),
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // (e^(tb) - e^(ta)) / (t(b-a))
        let etb = Expr::new_exp(Expr::new_mul(
            t.clone(),
            self.max.clone(),
        ));

        let eta = Expr::new_exp(Expr::new_mul(
            t.clone(),
            self.min.clone(),
        ));

        let num = Expr::new_sub(etb, eta);

        let range = Expr::new_sub(
            self.max.clone(),
            self.min.clone(),
        );

        let den = Expr::new_mul(t.clone(), range);

        simplify(&Expr::new_div(
            num, den,
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Binomial distribution with symbolic parameters.
#[derive(Debug, Clone)]

pub struct Binomial {
    pub n: Expr,
    pub p: Expr,
}

impl Distribution for Binomial {
    fn pdf(
        &self,
        k: &Expr,
    ) -> Expr {

        let n_choose_k = combinations(
            &self.n.clone(),
            k.clone(),
        );

        let p_k = Expr::new_pow(
            self.p.clone(),
            k.clone(),
        );

        let one_minus_p = Expr::new_sub(
            Expr::Constant(1.0),
            self.p.clone(),
        );

        let n_minus_k = Expr::new_sub(
            self.n.clone(),
            k.clone(),
        );

        let one_minus_p_pow = Expr::new_pow(
            one_minus_p,
            n_minus_k,
        );

        simplify(&Expr::new_mul(
            n_choose_k,
            Expr::new_mul(p_k, one_minus_p_pow),
        ))
    }

    fn cdf(
        &self,
        k: &Expr,
    ) -> Expr {

        // CDF is related to Regularized Incomplete Beta Function I_{1-p}(n-k, 1+k)
        // For now, let's represent it abstractly or using Summation
        Expr::Sum {
            body: Arc::new(
                self.pdf(&Expr::new_variable(
                    "i",
                )),
            ), // Use a dummy variable
            var: Arc::new(Expr::new_variable(
                "i",
            )),
            from: Arc::new(Expr::Constant(0.0)),
            to: Arc::new(k.clone()),
        }
    }

    fn expectation(&self) -> Expr {

        simplify(&Expr::new_mul(
            self.n.clone(),
            self.p.clone(),
        ))
    }

    fn variance(&self) -> Expr {

        let one_minus_p = Expr::new_sub(
            Expr::Constant(1.0),
            self.p.clone(),
        );

        simplify(&Expr::new_mul(
            self.n.clone(),
            Expr::new_mul(
                self.p.clone(),
                one_minus_p,
            ),
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // (1 - p + p e^t)^n
        let et = Expr::new_exp(t.clone());

        let one_minus_p = Expr::new_sub(
            Expr::Constant(1.0),
            self.p.clone(),
        );

        let pet = Expr::new_mul(self.p.clone(), et);

        let base = Expr::new_add(one_minus_p, pet);

        simplify(&Expr::new_pow(
            base,
            self.n.clone(),
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Poisson distribution with symbolic rate parameter λ.
#[derive(Debug, Clone)]

pub struct Poisson {
    pub rate: Expr,
}

impl Distribution for Poisson {
    fn pdf(
        &self,
        k: &Expr,
    ) -> Expr {

        let lambda_k = Expr::new_pow(
            self.rate.clone(),
            k.clone(),
        );

        let exp_neg_lambda = Expr::new_exp(Expr::new_neg(
            self.rate.clone(),
        ));

        let k_factorial = Expr::Factorial(Arc::new(k.clone()));

        simplify(&Expr::new_div(
            Expr::new_mul(
                lambda_k,
                exp_neg_lambda,
            ),
            k_factorial,
        ))
    }

    fn cdf(
        &self,
        k: &Expr,
    ) -> Expr {

        // Regularized Gamma Q(floor(k+1), lambda)
        // Or Summation
        Expr::Sum {
            body: Arc::new(
                self.pdf(&Expr::new_variable(
                    "i",
                )),
            ),
            var: Arc::new(Expr::new_variable(
                "i",
            )),
            from: Arc::new(Expr::Constant(0.0)),
            to: Arc::new(Expr::Floor(
                Arc::new(k.clone()),
            )),
        }
    }

    fn expectation(&self) -> Expr {

        self.rate.clone()
    }

    fn variance(&self) -> Expr {

        self.rate.clone()
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // exp(lambda * (e^t - 1))
        let et_minus_1 = Expr::new_sub(
            Expr::new_exp(t.clone()),
            Expr::Constant(1.0),
        );

        let arg = Expr::new_mul(
            self.rate.clone(),
            et_minus_1,
        );

        simplify(&Expr::new_exp(arg))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Bernoulli distribution with symbolic probability p.
#[derive(Debug, Clone)]

pub struct Bernoulli {
    pub p: Expr,
}

impl Distribution for Bernoulli {
    fn pdf(
        &self,
        k: &Expr,
    ) -> Expr {

        // p^k * (1-p)^(1-k) for k in {0, 1}
        let one_minus_p = Expr::new_sub(
            Expr::Constant(1.0),
            self.p.clone(),
        );

        let p_k = Expr::new_pow(
            self.p.clone(),
            k.clone(),
        );

        let one_minus_k = Expr::new_sub(
            Expr::Constant(1.0),
            k.clone(),
        );

        let one_minus_p_pow = Expr::new_pow(
            one_minus_p,
            one_minus_k,
        );

        simplify(&Expr::new_mul(
            p_k,
            one_minus_p_pow,
        ))
    }

    fn cdf(
        &self,
        k: &Expr,
    ) -> Expr {

        // Piecewise: 0 if k<0, 1-p if 0<=k<1, 1 if k>=1
        // Simplified symbolic rep?
        // Let's use logic/piecewise/conditions later if available. For now, just return a Sum.
        // Or specific for Bernoulli
        Expr::Sum {
            body: Arc::new(
                self.pdf(&Expr::new_variable(
                    "i",
                )),
            ),
            var: Arc::new(Expr::new_variable(
                "i",
            )),
            from: Arc::new(Expr::Constant(0.0)),
            to: Arc::new(Expr::Floor(
                Arc::new(k.clone()),
            )),
        }
    }

    fn expectation(&self) -> Expr {

        self.p.clone()
    }

    fn variance(&self) -> Expr {

        simplify(&Expr::new_mul(
            self.p.clone(),
            Expr::new_sub(
                Expr::Constant(1.0),
                self.p.clone(),
            ),
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // 1 - p + p e^t
        let et = Expr::new_exp(t.clone());

        let one_minus_p = Expr::new_sub(
            Expr::Constant(1.0),
            self.p.clone(),
        );

        let pet = Expr::new_mul(self.p.clone(), et);

        simplify(&Expr::new_add(
            one_minus_p,
            pet,
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents an Exponential distribution with symbolic rate λ.
#[derive(Debug, Clone)]

pub struct Exponential {
    pub rate: Expr,
}

impl Distribution for Exponential {
    fn pdf(
        &self,
        x: &Expr,
    ) -> Expr {

        simplify(&Expr::new_mul(
            self.rate.clone(),
            Expr::new_exp(Expr::new_neg(
                Expr::new_mul(
                    self.rate.clone(),
                    x.clone(),
                ),
            )),
        ))
    }

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr {

        simplify(&Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_exp(Expr::new_neg(
                Expr::new_mul(
                    self.rate.clone(),
                    x.clone(),
                ),
            )),
        ))
    }

    fn expectation(&self) -> Expr {

        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            self.rate.clone(),
        ))
    }

    fn variance(&self) -> Expr {

        simplify(&Expr::new_div(
            Expr::Constant(1.0),
            Expr::new_pow(
                self.rate.clone(),
                Expr::Constant(2.0),
            ),
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // lambda / (lambda - t)
        simplify(&Expr::new_div(
            self.rate.clone(),
            Expr::new_sub(
                self.rate.clone(),
                t.clone(),
            ),
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Gamma distribution with symbolic shape α and rate β.
#[derive(Debug, Clone)]

pub struct Gamma {
    pub shape: Expr,
    pub rate: Expr,
}

impl Distribution for Gamma {
    fn pdf(
        &self,
        x: &Expr,
    ) -> Expr {

        let term1_num = Expr::new_pow(
            self.rate.clone(),
            self.shape.clone(),
        );

        let term1_den = Expr::new_gamma(self.shape.clone());

        let term1 = Expr::new_div(term1_num, term1_den);

        let term2 = Expr::new_pow(
            x.clone(),
            Expr::new_sub(
                self.shape.clone(),
                Expr::Constant(1.0),
            ),
        );

        let term3 = Expr::new_exp(Expr::new_neg(
            Expr::new_mul(
                self.rate.clone(),
                x.clone(),
            ),
        ));

        simplify(&Expr::new_mul(
            term1,
            Expr::new_mul(term2, term3),
        ))
    }

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr {

        // Regularized Gamma P(alpha, beta*x)
        // GammaIncLower(alpha, beta*x) / Gamma(alpha)
        // We lack explicit GammaIncLower in basic Expr ops shown so far, maybe extend or use Integral
        Expr::Integral {
            integrand: Arc::new(
                self.pdf(&Expr::new_variable(
                    "t",
                )),
            ),
            var: Arc::new(Expr::new_variable(
                "t",
            )),
            lower_bound: Arc::new(Expr::Constant(0.0)),
            upper_bound: Arc::new(x.clone()),
        }
    }

    fn expectation(&self) -> Expr {

        simplify(&Expr::new_div(
            self.shape.clone(),
            self.rate.clone(),
        ))
    }

    fn variance(&self) -> Expr {

        // alpha / beta^2
        let beta_sq = Expr::new_pow(
            self.rate.clone(),
            Expr::Constant(2.0),
        );

        simplify(&Expr::new_div(
            self.shape.clone(),
            beta_sq,
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // (1 - t/beta)^(-alpha)
        let t_over_beta = Expr::new_div(
            t.clone(),
            self.rate.clone(),
        );

        let base = Expr::new_sub(
            Expr::Constant(1.0),
            t_over_beta,
        );

        let exp = Expr::new_neg(self.shape.clone());

        simplify(&Expr::new_pow(
            base, exp,
        ))
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Beta distribution with symbolic parameters α and β.
#[derive(Debug, Clone)]

pub struct Beta {
    pub alpha: Expr,
    pub beta: Expr,
}

impl Distribution for Beta {
    fn pdf(
        &self,
        x: &Expr,
    ) -> Expr {

        let num1 = Expr::new_pow(
            x.clone(),
            Expr::new_sub(
                self.alpha.clone(),
                Expr::Constant(1.0),
            ),
        );

        let one_minus_x = Expr::new_sub(
            Expr::Constant(1.0),
            x.clone(),
        );

        let num2 = Expr::new_pow(
            one_minus_x,
            Expr::new_sub(
                self.beta.clone(),
                Expr::Constant(1.0),
            ),
        );

        let den = Expr::new_beta(
            self.alpha.clone(),
            self.beta.clone(),
        );

        simplify(&Expr::new_div(
            Expr::new_mul(num1, num2),
            den,
        ))
    }

    fn cdf(
        &self,
        x: &Expr,
    ) -> Expr {

        // Regularized Incomplete Beta
        Expr::Integral {
            integrand: Arc::new(
                self.pdf(&Expr::new_variable(
                    "t",
                )),
            ),
            var: Arc::new(Expr::new_variable(
                "t",
            )),
            lower_bound: Arc::new(Expr::Constant(0.0)),
            upper_bound: Arc::new(x.clone()),
        }
    }

    fn expectation(&self) -> Expr {

        // alpha / (alpha + beta)
        let sum = Expr::new_add(
            self.alpha.clone(),
            self.beta.clone(),
        );

        simplify(&Expr::new_div(
            self.alpha.clone(),
            sum,
        ))
    }

    fn variance(&self) -> Expr {

        // (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
        let sum = Expr::new_add(
            self.alpha.clone(),
            self.beta.clone(),
        );

        let sum_sq = Expr::new_pow(
            sum.clone(),
            Expr::Constant(2.0),
        );

        let sum_plus_1 = Expr::new_add(
            sum,
            Expr::Constant(1.0),
        );

        let den = Expr::new_mul(sum_sq, sum_plus_1);

        let num = Expr::new_mul(
            self.alpha.clone(),
            self.beta.clone(),
        );

        simplify(&Expr::new_div(
            num, den,
        ))
    }

    fn mgf(
        &self,
        t: &Expr,
    ) -> Expr {

        // 1 + \sum_{k=1}^\infty ( \prod_{r=0}^{k-1} \frac{\alpha+r}{\alpha+\beta+r} ) \frac{t^k}{k!}
        // Represent as Confluent Hypergeometric Function 1F1(alpha; alpha+beta; t)
        // If we don't have Hypergeometric, return integral definition or Series
        // For now, let's leave as Integral of e^(tx) * pdf(x)
        Expr::Integral {
            integrand: Arc::new(Expr::new_mul(
                Expr::new_exp(Expr::new_mul(
                    t.clone(),
                    Expr::new_variable("x"),
                )),
                self.pdf(&Expr::new_variable(
                    "x",
                )),
            )),
            var: Arc::new(Expr::new_variable(
                "x",
            )),
            lower_bound: Arc::new(Expr::Constant(0.0)),
            upper_bound: Arc::new(Expr::Constant(1.0)),
        }
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}

/// Represents a Student's t-distribution with symbolic degrees of freedom ν.
#[derive(Debug, Clone)]

pub struct StudentT {
    pub nu: Expr,
}

impl Distribution for StudentT {
    fn pdf(
        &self,
        t: &Expr,
    ) -> Expr {

        let term1_num = Expr::new_gamma(Expr::new_div(
            Expr::new_add(
                self.nu.clone(),
                Expr::Constant(1.0),
            ),
            Expr::Constant(2.0),
        ));

        let term1_den_sqrt = Expr::new_sqrt(Expr::new_mul(
            self.nu.clone(),
            Expr::Pi,
        ));

        let term1_den_gamma = Expr::new_gamma(Expr::new_div(
            self.nu.clone(),
            Expr::Constant(2.0),
        ));

        let term1 = Expr::new_div(
            term1_num,
            Expr::new_mul(
                term1_den_sqrt,
                term1_den_gamma,
            ),
        );

        let term2_base = Expr::new_add(
            Expr::Constant(1.0),
            Expr::new_div(
                Expr::new_pow(
                    t.clone(),
                    Expr::Constant(2.0),
                ),
                self.nu.clone(),
            ),
        );

        let term2_exp = Expr::new_neg(Expr::new_div(
            Expr::new_add(
                self.nu.clone(),
                Expr::Constant(1.0),
            ),
            Expr::Constant(2.0),
        ));

        let term2 = Expr::new_pow(
            term2_base, term2_exp,
        );

        simplify(&Expr::new_mul(
            term1, term2,
        ))
    }

    fn cdf(
        &self,
        t: &Expr,
    ) -> Expr {

        // 1/2 + t * Gamma(...) * Hypergeometric(...)
        // Use Integral form
        Expr::Integral {
            integrand: Arc::new(
                self.pdf(&Expr::new_variable(
                    "x",
                )),
            ),
            var: Arc::new(Expr::new_variable(
                "x",
            )),
            lower_bound: Arc::new(Expr::NegativeInfinity),
            upper_bound: Arc::new(t.clone()),
        }
    }

    fn expectation(&self) -> Expr {

        // 0 if nu > 1, undefined otherwise
        // Assuming nu > 1
        Expr::Constant(0.0)
    }

    fn variance(&self) -> Expr {

        // nu / (nu - 2) for nu > 2
        simplify(&Expr::new_div(
            self.nu.clone(),
            Expr::new_sub(
                self.nu.clone(),
                Expr::Constant(2.0),
            ),
        ))
    }

    fn mgf(
        &self,
        _t: &Expr,
    ) -> Expr {

        // Does not exist for Student's T
        Expr::NoSolution // Or undefined
    }

    fn clone_box(&self) -> Arc<dyn Distribution> {

        Arc::new(self.clone())
    }
}
