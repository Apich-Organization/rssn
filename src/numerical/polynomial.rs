//! # Numerical Polynomial Operations
//!
//! This module provides a `Polynomial` struct and associated functions for numerical
//! operations on polynomials with `f64` coefficients. It supports evaluation,
//! differentiation, arithmetic (addition, subtraction, multiplication, division),
//! and finding real roots.

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::real_roots;

/// Represents a polynomial with f64 coefficients for numerical operations.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Polynomial {
    /// Coefficients of the polynomial, from highest degree to lowest.
    pub coeffs: Vec<f64>,
}

impl Polynomial {
    /// Creates a new `Polynomial` from coefficients (highest degree first).

    #[must_use]

    pub const fn new(
        coeffs: Vec<f64>
    ) -> Self {

        Self {
            coeffs,
        }
    }

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
    #[must_use]

    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {

        self.coeffs
            .iter()
            .fold(0.0, |acc, &c| {

                acc.mul_add(x, c)
            })
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
    ///
    /// # Errors
    /// Returns an error if root isolation or refinement fails.

    pub fn find_roots(
        &self
    ) -> Result<Vec<f64>, String> {

        let derivative =
            self.derivative();

        let isolating_intervals = real_roots::isolate_real_roots(self, 1e-9)?;

        let mut roots = Vec::new();

        for (a, b) in
            isolating_intervals
        {

            let mut guess =
                f64::midpoint(a, b);

            for _ in 0 .. 30 {

                let f_val =
                    self.eval(guess);

                let f_prime_val =
                    derivative
                        .eval(guess);

                if f_prime_val.abs()
                    < 1e-12
                {

                    break;
                }

                let next_guess = guess
                    - f_val
                        / f_prime_val;

                if (next_guess - guess)
                    .abs()
                    < 1e-12
                {

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
    #[must_use]

    pub fn derivative(&self) -> Self {

        if self.coeffs.len() <= 1 {

            return Self {
                coeffs: vec![0.0],
            };
        }

        let mut new_coeffs =
            Vec::with_capacity(
                self.coeffs.len() - 1,
            );

        let n = (self.coeffs.len() - 1)
            as f64;

        for (i, &c) in self
            .coeffs
            .iter()
            .enumerate()
            .take(self.coeffs.len() - 1)
        {

            new_coeffs.push(
                c * (n - i as f64),
            );
        }

        Self {
            coeffs: new_coeffs,
        }
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
    #[must_use]

    pub fn long_division(
        mut self,
        divisor: &Self,
    ) -> (Self, Self) {

        let mut quotient = vec![
                0.0;
                self.coeffs.len()
            ];

        let divisor_lead =
            divisor.coeffs[0];

        while self.coeffs.len()
            >= divisor.coeffs.len()
        {

            let lead_coeff =
                self.coeffs[0];

            let q_coeff = lead_coeff
                / divisor_lead;

            let deg_diff = self
                .coeffs
                .len()
                - divisor.coeffs.len();

            quotient[deg_diff] =
                q_coeff;

            for i in 0 .. divisor
                .coeffs
                .len()
            {

                self.coeffs[i] -=
                    divisor.coeffs[i]
                        * q_coeff;
            }

            self.coeffs
                .remove(0);
        }

        (
            Self {
                coeffs: quotient,
            },
            self,
        )
    }

    /// Returns the degree of the polynomial.
    #[must_use]

    pub const fn degree(
        &self
    ) -> usize {

        if self
            .coeffs
            .is_empty()
        {

            return 0;
        }

        self.coeffs.len() - 1
    }

    /// Checks if the polynomial is zero (within epsilon).
    #[must_use]

    pub fn is_zero(
        &self,
        epsilon: f64,
    ) -> bool {

        self.coeffs
            .iter()
            .all(|&c| c.abs() < epsilon)
    }

    /// Returns the indefinite integral of the polynomial (constant of integration = 0).
    #[must_use]

    pub fn integral(&self) -> Self {

        if self
            .coeffs
            .is_empty()
        {

            return Self {
                coeffs: vec![0.0],
            };
        }

        let mut new_coeffs =
            Vec::with_capacity(
                self.coeffs.len() + 1,
            );

        let d = self.degree() as f64;

        for (i, &c) in self
            .coeffs
            .iter()
            .enumerate()
        {

            new_coeffs.push(
                c / (d - i as f64
                    + 1.0),
            );
        }

        new_coeffs.push(0.0);

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Add for Polynomial {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        let max_len = self
            .coeffs
            .len()
            .max(rhs.coeffs.len());

        let mut new_coeffs =
            vec![0.0; max_len];

        let self_pad =
            max_len - self.coeffs.len();

        let rhs_pad =
            max_len - rhs.coeffs.len();

        for (i, var) in new_coeffs
            .iter_mut()
            .enumerate()
            .take(max_len)
        {

            let c1 = if i >= self_pad {

                self.coeffs
                    [i - self_pad]
            } else {

                0.0
            };

            let c2 = if i >= rhs_pad {

                rhs.coeffs[i - rhs_pad]
            } else {

                0.0
            };

            *var = c1 + c2;
        }

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        let max_len = self
            .coeffs
            .len()
            .max(rhs.coeffs.len());

        let mut new_coeffs =
            vec![0.0; max_len];

        let self_pad =
            max_len - self.coeffs.len();

        let rhs_pad =
            max_len - rhs.coeffs.len();

        for (i, var) in new_coeffs
            .iter_mut()
            .enumerate()
            .take(max_len)
        {

            let c1 = if i >= self_pad {

                self.coeffs
                    [i - self_pad]
            } else {

                0.0
            };

            let c2 = if i >= rhs_pad {

                rhs.coeffs[i - rhs_pad]
            } else {

                0.0
            };

            *var = c1 - c2;
        }

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {

        if self
            .coeffs
            .is_empty()
            || rhs
                .coeffs
                .is_empty()
        {

            return Self {
                coeffs: vec![],
            };
        }

        let mut new_coeffs = vec![
                0.0;
                self.coeffs.len()
                    + rhs.coeffs.len()
                    - 1
            ];

        for (i, &c1) in self
            .coeffs
            .iter()
            .enumerate()
        {

            for (j, &c2) in rhs
                .coeffs
                .iter()
                .enumerate()
            {

                new_coeffs[i + j] +=
                    c1 * c2;
            }
        }

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Div for Polynomial {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {

        self.long_division(&rhs)
            .0
    }
}

impl Mul<f64> for Polynomial {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {

        let new_coeffs = self
            .coeffs
            .iter()
            .map(|&c| c * rhs)
            .collect();

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Div<f64> for Polynomial {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {

        let new_coeffs = self
            .coeffs
            .iter()
            .map(|&c| c / rhs)
            .collect();

        Self {
            coeffs: new_coeffs,
        }
    }
}

impl Polynomial {
    /// Divides the polynomial by a scalar.
    ///
    /// # Arguments
    /// * `rhs` - The scalar to divide by.
    ///
    /// # Returns
    /// A `Result` containing the new `Polynomial`, or an error if the scalar is zero.
    ///
    /// # Errors
    /// Returns an error if the divisor scalar is zero.

    pub fn div_scalar(
        self,
        rhs: f64,
    ) -> Result<Self, String> {

        if rhs == 0.0 {

            return Err("Division by \
                        zero scalar"
                .to_string());
        }

        let new_coeffs = self
            .coeffs
            .iter()
            .map(|&c| c / rhs)
            .collect();

        Ok(Self {
            coeffs: new_coeffs,
        })
    }
}
