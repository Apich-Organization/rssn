//! # Numerical Polynomial Operations
//!
//! This module provides a `Polynomial` struct and associated functions for numerical
//! operations on polynomials with `f64` coefficients. It supports evaluation,
//! differentiation, arithmetic (addition, subtraction, multiplication, division),
//! and finding real roots.

use crate::numerical::real_roots;
use std::ops::{Add, Div, Mul, Sub};

/// Represents a polynomial with f64 coefficients for numerical operations.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<f64>, // Coefficients in descending order of degree
}

impl Polynomial {
    /// Evaluates the polynomial at a given point `x` using Horner's method.
    ///
    /// Horner's method is an efficient algorithm for evaluating polynomials.
    /// For a polynomial `P(x) = c_n*x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0`,
    /// it computes `P(x) = (...((c_n*x + c_{n-1})*x + c_{n-2})*x + ... + c_1)*x + c_0`.
    ///
    /// # Arguments
    /// * `x` - The point at which to evaluate the polynomial.
    ///
    /// # Returns
    /// The value of the polynomial at `x` as an `f64`.
    pub fn eval(&self, x: f64) -> f64 {
        self.coeffs.iter().fold(0.0, |acc, &c| acc * x + c)
    }

    /// Finds the real roots of the polynomial.
    ///
    /// This method combines Sturm's theorem for root isolation with Newton's method
    /// for refining the roots. Sturm's theorem provides disjoint intervals, each
    /// containing exactly one real root, which are then used as starting points
    /// for Newton's method to converge to the root.
    ///
    /// # Returns
    /// A `Result` containing a `Vec<f64>` of the real roots found, or an error string
    /// if root isolation or refinement fails.
    pub fn find_roots(&self) -> Result<Vec<f64>, String> {
        let derivative = self.derivative();

        let isolating_intervals = real_roots::isolate_real_roots(self, 1e-9)?;
        let mut roots = Vec::new();

        for (a, b) in isolating_intervals {
            let mut guess = (a + b) / 2.0;
            // Use Newton's method to refine the root within the isolating interval
            for _ in 0..30 {
                // Max 30 iterations for refinement
                let f_val = self.eval(guess);
                let f_prime_val = derivative.eval(guess);
                if f_prime_val.abs() < 1e-12 {
                    break;
                } // Avoid division by zero
                let next_guess = guess - f_val / f_prime_val;
                if (next_guess - guess).abs() < 1e-12 {
                    guess = next_guess;
                    break;
                }
                guess = next_guess;
            }
            roots.push(guess);
        }
        Ok(roots)
    }

    /// Returns the derivative of the polynomial.
    ///
    /// The derivative is computed by applying the power rule to each term.
    /// For a term `c*x^n`, its derivative is `(c*n)*x^(n-1)`.
    ///
    /// # Returns
    /// A new `Polynomial` representing the derivative.
    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Polynomial { coeffs: vec![0.0] };
        }
        let mut new_coeffs = Vec::with_capacity(self.coeffs.len() - 1);
        let n = (self.coeffs.len() - 1) as f64;
        for (i, &c) in self.coeffs.iter().enumerate().take(self.coeffs.len() - 1) {
            new_coeffs.push(c * (n - i as f64));
        }
        Polynomial { coeffs: new_coeffs }
    }

    /// Performs polynomial long division.
    ///
    /// This method divides the current polynomial (dividend) by another polynomial (divisor).
    ///
    /// # Arguments
    /// * `divisor` - The polynomial to divide by.
    ///
    /// # Returns
    /// A tuple `(quotient, remainder)` as `Polynomial`s.
    pub fn long_division(mut self, divisor: Self) -> (Self, Self) {
        let mut quotient = vec![0.0; self.coeffs.len()];
        let divisor_lead = divisor.coeffs[0];

        while self.coeffs.len() >= divisor.coeffs.len() {
            let lead_coeff = self.coeffs[0];
            let q_coeff = lead_coeff / divisor_lead;
            let deg_diff = self.coeffs.len() - divisor.coeffs.len();
            quotient[deg_diff] = q_coeff;

            for i in 0..divisor.coeffs.len() {
                self.coeffs[i] -= divisor.coeffs[i] * q_coeff;
            }
            self.coeffs.remove(0);
        }
        (Polynomial { coeffs: quotient }, self)
    }
}

impl Add for Polynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut new_coeffs = vec![0.0; max_len];
        let self_pad = max_len - self.coeffs.len();
        let rhs_pad = max_len - rhs.coeffs.len();

        for i in 0..max_len {
            let c1 = if i >= self_pad {
                self.coeffs[i - self_pad]
            } else {
                0.0
            };
            let c2 = if i >= rhs_pad {
                rhs.coeffs[i - rhs_pad]
            } else {
                0.0
            };
            new_coeffs[i] = c1 + c2;
        }
        Polynomial { coeffs: new_coeffs }
    }
}

impl Sub for Polynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut new_coeffs = vec![0.0; max_len];
        let self_pad = max_len - self.coeffs.len();
        let rhs_pad = max_len - rhs.coeffs.len();

        for i in 0..max_len {
            let c1 = if i >= self_pad {
                self.coeffs[i - self_pad]
            } else {
                0.0
            };
            let c2 = if i >= rhs_pad {
                rhs.coeffs[i - rhs_pad]
            } else {
                0.0
            };
            new_coeffs[i] = c1 - c2;
        }
        Polynomial { coeffs: new_coeffs }
    }
}

impl Mul for Polynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Polynomial { coeffs: vec![] };
        }
        let mut new_coeffs = vec![0.0; self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, &c1) in self.coeffs.iter().enumerate() {
            for (j, &c2) in rhs.coeffs.iter().enumerate() {
                new_coeffs[i + j] += c1 * c2;
            }
        }
        Polynomial { coeffs: new_coeffs }
    }
}

impl Div for Polynomial {
    type Output = Self; // Returns the quotient
    fn div(self, rhs: Self) -> Self {
        self.long_division(rhs).0
    }
}

impl Mul<f64> for Polynomial {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let new_coeffs = self.coeffs.iter().map(|&c| c * rhs).collect();
        Polynomial { coeffs: new_coeffs }
    }
}

impl Div<f64> for Polynomial {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        if rhs == 0.0 {
            panic!("Division by zero scalar");
        }
        let new_coeffs = self.coeffs.iter().map(|&c| c / rhs).collect();
        Polynomial { coeffs: new_coeffs }
    }
}
