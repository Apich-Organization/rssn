//! # Integral Equations
//!
//! This module provides structures and methods for solving various types of integral equations.
//! An integral equation is an equation in which an unknown function appears under an integral sign.
//! It includes solvers for Fredholm and Volterra integral equations of the second kind,
//! using methods like successive approximations (Neumann Series) and conversion to ODEs.
//! Singular integral equations are also supported.

use std::sync::Arc;
/*
use crate::symbolic::core::Expr;
use crate::symbolic::calculus::{differentiate, integrate, substitute};
use crate::symbolic::simplify::{is_zero, simplify};
use crate::symbolic::solve::solve_linear_system;
/// Represents a Fredholm integral equation of the second kind:
/// `y(x) = f(x) + lambda * integral_a_b(K(x, t) * y(t) dt)`
#[derive(Debug, Clone)]
*/

use crate::symbolic::calculus::{differentiate, integrate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use crate::symbolic::solve::solve_linear_system;

// =====================================================================================
// region: Fredholm Integral Equations
// =====================================================================================

/// Represents a Fredholm integral equation of the second kind.
///
/// The equation has the form: `y(x) = f(x) + lambda * integral_a_b(K(x, t) * y(t) dt)`,
/// where `y(x)` is the unknown function to be solved for.
#[derive(Debug, Clone)]
pub struct FredholmEquation {
    /// The unknown function `y(x)`.
    pub y_x: Expr,
    /// The known function `f(x)`.
    pub f_x: Expr,
    /// The constant parameter `lambda`.
    pub lambda: Expr,
    /// The kernel of the integral, `K(x, t)`.
    pub kernel: Expr,
    /// The lower bound of integration, `a`.
    pub lower_bound: Expr,
    /// The upper bound of integration, `b`.
    pub upper_bound: Expr,
    /// The main variable of the functions, `x`.
    pub var_x: String,
    /// The integration variable, `t`.
    pub var_t: String,
}

impl FredholmEquation {
    /// Creates a new instance of a Fredholm integral equation of the second kind.
    ///
    /// # Arguments
    /// * `y_x` - The unknown function `y(x)`.
    /// * `f_x` - The known function `f(x)`.
    /// * `lambda` - The constant parameter `lambda`.
    /// * `kernel` - The kernel of the integral, `K(x, t)`.
    /// * `lower_bound` - The lower bound of integration, `a`.
    /// * `upper_bound` - The upper bound of integration, `b`.
    /// * `var_x` - The main variable of the functions, `x`.
    /// * `var_t` - The integration variable, `t`.
    ///
    /// # Returns
    /// A new `FredholmEquation` instance.
    pub fn new(
        y_x: Expr,
        f_x: Expr,
        lambda: Expr,
        kernel: Expr,
        lower_bound: Expr,
        upper_bound: Expr,
        var_x: String,
        var_t: String,
    ) -> Self {
        FredholmEquation {
            y_x,
            f_x,
            lambda,
            kernel,
            lower_bound,
            upper_bound,
            var_x,
            var_t,
        }
    }

    /// Solves the Fredholm integral equation using the method of successive approximations (Neumann Series).
    ///
    /// This iterative method constructs a sequence of functions `y_n(x)` that converges to the solution.
    /// The sequence is defined by:
    /// `y_0(x) = f(x)`
    /// `y_{n+1}(x) = f(x) + lambda * integral(K(x, t) * y_n(t) dt)`
    /// This method is generally applicable when `lambda` is small enough for the series to converge.
    ///
    /// # Arguments
    /// * `iterations` - The number of iterations to perform.
    ///
    /// # Returns
    /// An `Expr` representing the approximate solution `y(x)`.
    pub fn solve_neumann_series(&self, iterations: usize) -> Expr {
        let mut y_n = self.f_x.clone(); // y_0(x) = f(x)

        for _ in 0..iterations {
            let integral_term = Expr::Mul(
                Arc::new(self.kernel.clone()),
                Arc::new(substitute(
                    &y_n,
                    &self.var_x,
                    &Expr::Variable(self.var_t.clone()),
                )),
            );
            let integrated_val = integrate(
                &integral_term,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&self.upper_bound),
            );
            let next_y_n = simplify(Expr::Add(
                Arc::new(self.f_x.clone()),
                Arc::new(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(integrated_val),
                )),
            ));
            y_n = next_y_n;
        }
        y_n
    }

    /// Solves a Fredholm integral equation of the second kind with a separable (or degenerate) kernel.
    ///
    /// A separable kernel can be written as a finite sum of products of functions of a single variable:
    /// `K(x, t) = sum_{i=1 to m} a_i(x) * b_i(t)`.
    /// This method transforms the integral equation into a system of linear algebraic equations for
    /// unknown coefficients `c_i`, and then constructs the solution.
    ///
    /// # Arguments
    /// * `a_funcs` - A vector of `Expr` representing the `a_i(x)` functions.
    /// * `b_funcs` - A vector of `Expr` representing the `b_i(t)` functions.
    ///
    /// # Returns
    /// A `Result<Expr, String>` which is the solution `y(x)` on success, or an error message.
    pub fn solve_separable_kernel(
        &self,
        a_funcs: Vec<Expr>,
        b_funcs: Vec<Expr>,
    ) -> Result<Expr, String> {
        if a_funcs.len() != b_funcs.len() {
            return Err("Number of a_i(x) functions must match b_i(t) functions".to_string());
        }

        let m = a_funcs.len();
        let c_vars: Vec<String> = (0..m).map(|i| format!("c{}", i)).collect();
        let mut system_eqs: Vec<Expr> = Vec::new();

        // The solution y(x) is of the form: y(x) = f(x) + lambda * sum(c_i * a_i(x))
        // where the coefficients c_i are given by: c_i = integral(b_i(t) * y(t) dt)
        // Substituting y(t) into the definition of c_i leads to a system of linear equations.

        for k in 0..m {
            let b_k_t = substitute(
                &b_funcs[k],
                &self.var_x,
                &Expr::Variable(self.var_t.clone()),
            );
            let f_t = substitute(&self.f_x, &self.var_x, &Expr::Variable(self.var_t.clone()));

            // beta_k = integral(b_k(t) * f(t) dt)
            let beta_k_integrand = simplify(Expr::Mul(Arc::new(b_k_t.clone()), Arc::new(f_t)));
            let beta_k = integrate(
                &beta_k_integrand,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&self.upper_bound),
            );

            let mut lhs_sum_terms = Vec::new();
            for i in 0..m {
                let a_i_t = substitute(
                    &a_funcs[i],
                    &self.var_x,
                    &Expr::Variable(self.var_t.clone()),
                );

                // alpha_ki = integral(b_k(t) * a_i(t) dt)
                let alpha_ki_integrand =
                    simplify(Expr::Mul(Arc::new(b_k_t.clone()), Arc::new(a_i_t)));
                let alpha_ki = integrate(
                    &alpha_ki_integrand,
                    &self.var_t,
                    Some(&self.lower_bound),
                    Some(&self.upper_bound),
                );

                let c_i_var = Expr::Variable(c_vars[i].clone());
                let term = simplify(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(Expr::Mul(Arc::new(c_i_var), Arc::new(alpha_ki))),
                ));
                lhs_sum_terms.push(term);
            }

            // Equation: c_k - lambda * sum(c_i * alpha_ki) = beta_k
            let c_k_var = Expr::Variable(c_vars[k].clone());
            let sum_of_terms = lhs_sum_terms
                .into_iter()
                .fold(Expr::Constant(0.0), |acc, x| {
                    simplify(Expr::Add(Arc::new(acc), Arc::new(x)))
                });

            let equation_lhs = simplify(Expr::Sub(Arc::new(c_k_var), Arc::new(sum_of_terms)));
            system_eqs.push(Expr::Eq(Arc::new(equation_lhs), Arc::new(beta_k)));
        }

        let c_solved = solve_linear_system(&Expr::System(system_eqs), &c_vars)?;

        // Reconstruct y(x) = f(x) + lambda * sum(c_i * a_i(x))
        let mut solution_sum_terms = Vec::new();
        for i in 0..m {
            let c_i_val = c_solved[i].clone();
            let a_i_x = a_funcs[i].clone();
            let term = simplify(Expr::Mul(Arc::new(c_i_val), Arc::new(a_i_x)));
            solution_sum_terms.push(term);
        }

        let sum_of_solution_terms = solution_sum_terms
            .into_iter()
            .fold(Expr::Constant(0.0), |acc, x| {
                simplify(Expr::Add(Arc::new(acc), Arc::new(x)))
            });

        let final_solution = simplify(Expr::Add(
            Arc::new(self.f_x.clone()),
            Arc::new(Expr::Mul(
                Arc::new(self.lambda.clone()),
                Arc::new(sum_of_solution_terms),
            )),
        ));

        Ok(final_solution)
    }
}

// =====================================================================================
// endregion: Fredholm Integral Equations
// =====================================================================================

// =====================================================================================
// region: Volterra Integral Equations
// =====================================================================================

/// Represents a Volterra integral equation of the second kind.
///
/// The equation has the form: `y(x) = f(x) + lambda * integral_a_x(K(x, t) * y(t) dt)`.
/// It is similar to the Fredholm equation, but the upper limit of integration is the variable `x`.
#[derive(Debug, Clone)]
pub struct VolterraEquation {
    /// The unknown function `y(x)`.
    pub y_x: Expr,
    /// The known function `f(x)`.
    pub f_x: Expr,
    /// The constant parameter `lambda`.
    pub lambda: Expr,
    /// The kernel of the integral, `K(x, t)`.
    pub kernel: Expr,
    /// The lower bound of integration, `a`.
    pub lower_bound: Expr,
    /// The main variable of the functions, `x`.
    pub var_x: String,
    /// The integration variable, `t`.
    pub var_t: String,
}

impl VolterraEquation {
    /// Creates a new instance of a Volterra integral equation of the second kind.
    ///
    /// # Arguments
    /// * `y_x` - The unknown function `y(x)`.
    /// * `f_x` - The known function `f(x)`.
    /// * `lambda` - The constant parameter `lambda`.
    /// * `kernel` - The kernel of the integral, `K(x, t)`.
    /// * `lower_bound` - The lower bound of integration, `a`.
    /// * `var_x` - The main variable of the functions, `x`.
    /// * `var_t` - The integration variable, `t`.
    ///
    /// # Returns
    /// A new `VolterraEquation` instance.
    pub fn new(
        y_x: Expr,
        f_x: Expr,
        lambda: Expr,
        kernel: Expr,
        lower_bound: Expr,
        var_x: String,
        var_t: String,
    ) -> Self {
        VolterraEquation {
            y_x,
            f_x,
            lambda,
            kernel,
            lower_bound,
            var_x,
            var_t,
        }
    }

    /// Solves the Volterra integral equation using the method of successive approximations.
    ///
    /// This iterative method is analogous to the Neumann series for Fredholm equations.
    /// The sequence is defined by:
    /// `y_0(x) = f(x)`
    /// `y_{n+1}(x) = f(x) + lambda * integral_a_x(K(x, t) * y_n(t) dt)`
    ///
    /// # Arguments
    /// * `iterations` - The number of iterations to perform.
    ///
    /// # Returns
    /// An `Expr` representing the approximate solution `y(x)`.
    pub fn solve_successive_approximations(&self, iterations: usize) -> Expr {
        let mut y_n = self.f_x.clone(); // y_0(x) = f(x)

        for _ in 0..iterations {
            let y_n_t = substitute(&y_n, &self.var_x, &Expr::Variable(self.var_t.clone()));
            let integral_term = Expr::Mul(Arc::new(self.kernel.clone()), Arc::new(y_n_t));

            let integrated_val = integrate(
                &integral_term,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&Expr::Variable(self.var_x.clone())),
            );

            let next_y_n = simplify(Expr::Add(
                Arc::new(self.f_x.clone()),
                Arc::new(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(integrated_val),
                )),
            ));
            y_n = next_y_n;
        }
        y_n
    }

    /// Solves the Volterra equation by converting it into an Ordinary Differential Equation (ODE).
    ///
    /// This method is applicable if the kernel `K(x, t)` is a function of `x` only, or if
    /// differentiating the equation with respect to `x` (using the Leibniz integral rule)
    /// results in a solvable ODE.
    ///
    /// # Returns
    /// A `Result<Expr, String>` which is the solution `y(x)` on success, or an error message.
    pub fn solve_by_differentiation(&self) -> Result<Expr, String> {
        // Differentiate the entire equation y(x) = f(x) + lambda * integral_a_x(K(x,t)y(t)dt) w.r.t. x
        // Using Leibniz rule: d/dx integral_a(x)_b(x) F(x,t) dt = F(x,b(x))*b'(x) - F(x,a(x))*a'(x) + integral_a(x)_b(x) dF/dx dt
        // Here, a(x)=a (const), b(x)=x, F(x,t) = K(x,t)y(t)
        // y'(x) = f'(x) + lambda * [ K(x,x)y(x) + integral_a_x(dK/dx * y(t) dt) ]

        let y_prime = differentiate(&self.y_x, &self.var_x);
        let f_prime = differentiate(&self.f_x, &self.var_x);

        // K(x,x)
        let k_x_x = substitute(
            &self.kernel,
            &self.var_t,
            &Expr::Variable(self.var_x.clone()),
        );
        let term1 = simplify(Expr::Mul(Arc::new(k_x_x), Arc::new(self.y_x.clone())));

        // integral_a_x(dK/dx * y(t) dt)
        let dk_dx = differentiate(&self.kernel, &self.var_x);
        let y_t = substitute(&self.y_x, &self.var_x, &Expr::Variable(self.var_t.clone()));
        let integrand = simplify(Expr::Mul(Arc::new(dk_dx), Arc::new(y_t)));
        let integral_term = integrate(
            &integrand,
            &self.var_t,
            Some(&self.lower_bound),
            Some(&Expr::Variable(self.var_x.clone())),
        );

        let rhs = simplify(Expr::Add(
            Arc::new(f_prime),
            Arc::new(Expr::Mul(
                Arc::new(self.lambda.clone()),
                Arc::new(Expr::Add(Arc::new(term1), Arc::new(integral_term))),
            )),
        ));

        // Now we have an integro-differential equation. If the integral term disappears
        // (i.e., if K(x,t) does not depend on x), we get a simple ODE.
        let ode_expr = Expr::Eq(Arc::new(y_prime), Arc::new(rhs));

        // This is a placeholder for a more advanced solver.
        // A real implementation would need to check if `ode_expr` is a valid ODE
        // (i.e., no integral terms left) and then call the ODE solver.
        // For now, we return the derived ODE expression.
        // In a full system, you would call: solve_ode(&ode_expr, &self.y_x)
        Err(format!("Conversion to ODE resulted in: {}", ode_expr))
    }
}

// =====================================================================================
// endregion: Volterra Integral Equations
// =====================================================================================

// =====================================================================================
// region: Singular Integral Equations
// =====================================================================================

/// Solves the airfoil singular integral equation.
///
/// The equation is a specific type of Cauchy-type singular integral equation given by:
/// `(1/π) * ∫[-1, 1] y(t)/(t-x) dt = f(x)` for `x` in `(-1, 1)`.
/// The integral is taken as the Cauchy Principal Value.
///
/// The solution is known in closed form:
/// `y(x) = (-1 / (π * sqrt(1-x^2))) * ∫[-1, 1] (sqrt(1-t^2)/(t-x)) * f(t) dt + C / sqrt(1-x^2)`
/// where C is an arbitrary constant.
///
/// # Arguments
/// * `f_x` - The known function `f(x)`.
/// * `var_x` - The variable `x`.
/// * `var_t` - The integration variable `t`.
///
/// # Returns
/// An `Expr` representing the solution `y(x)` with a constant of integration `C`.
pub fn solve_airfoil_equation(f_x: &Expr, var_x: &str, var_t: &str) -> Expr {
    let one = Expr::Constant(1.0);
    let neg_one = Expr::Constant(-1.0);
    let pi = Expr::Pi;

    // sqrt(1-t^2)
    let sqrt_1_minus_t2 = Expr::Sqrt(Arc::new(Expr::Sub(
        Arc::new(one.clone()),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(var_t.to_string())),
            Arc::new(Expr::Constant(2.0)),
        )),
    )));

    // t - x
    let t_minus_x = Expr::Sub(
        Arc::new(Expr::Variable(var_t.to_string())),
        Arc::new(Expr::Variable(var_x.to_string())),
    );

    // f(t)
    let f_t = substitute(f_x, var_x, &Expr::Variable(var_t.to_string()));

    // Integrand: (sqrt(1-t^2)/(t-x)) * f(t)
    let integrand = Expr::Mul(
        Arc::new(Expr::Div(Arc::new(sqrt_1_minus_t2), Arc::new(t_minus_x))),
        Arc::new(f_t),
    );

    // The integral part: ∫[-1, 1] ... dt
    let integral_part = integrate(&integrand, var_t, Some(&neg_one), Some(&one));

    // sqrt(1-x^2)
    let sqrt_1_minus_x2 = Expr::Sqrt(Arc::new(Expr::Sub(
        Arc::new(one.clone()),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(var_x.to_string())),
            Arc::new(Expr::Constant(2.0)),
        )),
    )));

    // First term: (-1 / (π * sqrt(1-x^2))) * integral
    let factor1 = Expr::Div(
        Arc::new(Expr::Constant(-1.0)),
        Arc::new(Expr::Mul(
            Arc::new(pi.clone()),
            Arc::new(sqrt_1_minus_x2.clone()),
        )),
    );
    let term1 = Expr::Mul(Arc::new(factor1), Arc::new(integral_part));

    // Second term: C / sqrt(1-x^2)
    let const_c = Expr::Variable("C".to_string());
    let term2 = Expr::Div(Arc::new(const_c), Arc::new(sqrt_1_minus_x2));

    simplify(Expr::Add(Arc::new(term1), Arc::new(term2)))
}

// =====================================================================================
// endregion: Singular Integral Equations
// =====================================================================================

/*
impl FredholmEquation {
    pub fn new(
        y_x: Expr,
        f_x: Expr,
        lambda: Expr,
        kernel: Expr,
        lower_bound: Expr,
        upper_bound: Expr,
        var_x: String,
        var_t: String,
    ) -> Self {
        FredholmEquation {
            y_x,
            f_x,
            lambda,
            kernel,
            lower_bound,
            upper_bound,
            var_x,
            var_t,
        }
    }

    /// Solves the Fredholm integral equation using the method of successive approximations (Neumann Series).
    /// This method is applicable when `lambda` is small enough.
    ///
    /// # Arguments
    /// * `iterations` - The number of iterations for the successive approximation.
    ///
    /// # Returns
    /// An `Expr` representing the approximate solution `y(x)`.
    pub fn solve_neumann_series(&self, iterations: usize) -> Expr {
        let mut y_n = self.f_x.clone(); // y_0(x) = f(x)

        for _ in 0..iterations {
            let integral_term = Expr::Mul(
                Arc::new(self.kernel.clone()),
                Arc::new(substitute(
                    &y_n,
                    &self.var_x,
                    &Expr::Variable(self.var_t.clone()),
                )),
            );
            let integrated_val = integrate(
                &integral_term,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&self.upper_bound),
            );
            let next_y_n = simplify(Expr::Add(
                Arc::new(self.f_x.clone()),
                Arc::new(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(integrated_val),
                )),
            ));
            y_n = next_y_n;
        }
        y_n
    }

    /// Solves a Fredholm integral equation of the second kind with a separable kernel.
    /// A separable kernel can be written as `K(x, t) = sum(a_i(x) * b_i(t))`.
    ///
    /// # Arguments
    /// * `a_funcs` - A vector of `Expr` representing `a_i(x)` functions.
    /// * `b_funcs` - A vector of `Expr` representing `b_i(t)` functions.
    ///
    /// # Returns
    /// An `Expr` representing the solution `y(x)`, or an error if the system is singular.
    pub fn solve_separable_kernel(
        &self,
        a_funcs: Vec<Expr>,
        b_funcs: Vec<Expr>,
    ) -> Result<Expr, String> {
        if a_funcs.len() != b_funcs.len() {
            return Err("Number of a_i(x) functions must match b_i(t) functions".to_string());
        }

        let m = a_funcs.len();
        let mut c_vars: Vec<String> = (0..m).map(|i| format!("c{}", i)).collect();
        let mut system_eqs: Vec<Expr> = Vec::new();

        // The solution y(x) will be of the form y(x) = f(x) + lambda * sum(c_i * a_i(x))
        // where c_i = integral_a_b(b_i(t) * y(t) dt)
        // Substituting y(t) into the c_i definition leads to a system of linear equations for c_i.
        // c_k = integral_a_b(b_k(t) * (f(t) + lambda * sum(c_i * a_i(t))) dt)
        // c_k = integral_a_b(b_k(t) * f(t) dt) + lambda * sum(c_i * integral_a_b(b_k(t) * a_i(t) dt))
        // Let beta_k = integral_a_b(b_k(t) * f(t) dt)
        // Let alpha_ki = integral_a_b(b_k(t) * a_i(t) dt)
        // Then the system is: c_k - lambda * sum(c_i * alpha_ki) = beta_k

        for k in 0..m {
            let b_k_t = substitute(
                &b_funcs[k],
                &self.var_x,
                &Expr::Variable(self.var_t.clone()),
            );
            let f_t = substitute(&self.f_x, &self.var_x, &Expr::Variable(self.var_t.clone()));

            // Calculate beta_k = integral_a_b(b_k(t) * f(t) dt)
            let beta_k_integrand = simplify(Expr::Mul(Arc::new(b_k_t.clone()), Arc::new(f_t)));
            let beta_k = integrate(
                &beta_k_integrand,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&self.upper_bound),
            );

            let mut lhs_sum_terms = Vec::new();
            for i in 0..m {
                let a_i_t = substitute(
                    &a_funcs[i],
                    &self.var_x,
                    &Expr::Variable(self.var_t.clone()),
                );

                // Calculate alpha_ki = integral_a_b(b_k(t) * a_i(t) dt)
                let alpha_ki_integrand =
                    simplify(Expr::Mul(Arc::new(b_k_t.clone()), Arc::new(a_i_t)));
                let alpha_ki = integrate(
                    &alpha_ki_integrand,
                    &self.var_t,
                    Some(&self.lower_bound),
                    Some(&self.upper_bound),
                );

                // Term: lambda * c_i * alpha_ki
                let c_i_var = Expr::Variable(c_vars[i].clone());
                let term = simplify(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(Expr::Mul(Arc::new(c_i_var), Arc::new(alpha_ki))),
                ));
                lhs_sum_terms.push(term);
            }

            // Construct the equation: c_k - sum(lambda * c_i * alpha_ki) = beta_k
            let c_k_var = Expr::Variable(c_vars[k].clone());
            let sum_of_terms = if lhs_sum_terms.is_empty() {
                Expr::Constant(0.0)
            } else {
                lhs_sum_terms
                    .into_iter()
                    .fold(Expr::Constant(0.0), |acc, x| {
                        Expr::Add(Arc::new(acc), Arc::new(x))
                    })
            };

            let equation_lhs = simplify(Expr::Sub(Arc::new(c_k_var), Arc::new(sum_of_terms)));
            system_eqs.push(Expr::Eq(Arc::new(equation_lhs), Arc::new(beta_k)));
        }

        // Solve the system for c_i
        let c_solved = solve_linear_system(&Expr::System(system_eqs), &c_vars)?;

        // Reconstruct y(x) = f(x) + lambda * sum(c_i * a_i(x))
        let mut solution_sum_terms = Vec::new();
        // The `c_solved` is now `Vec<Expr>`, directly containing the solutions for c_i.
        let c_vec = c_solved;

        for i in 0..m {
            let c_i_val = c_vec[i].clone();
            let a_i_x = a_funcs[i].clone();
            let term = simplify(Expr::Mul(Arc::new(c_i_val), Arc::new(a_i_x)));
            solution_sum_terms.push(term);
        }

        let sum_of_solution_terms = if solution_sum_terms.is_empty() {
            Expr::Constant(0.0)
        } else {
            solution_sum_terms
                .into_iter()
                .fold(Expr::Constant(0.0), |acc, x| {
                    Expr::Add(Arc::new(acc), Arc::new(x))
                })
        };

        let final_solution = simplify(Expr::Add(
            Arc::new(self.f_x.clone()),
            Arc::new(Expr::Mul(
                Arc::new(self.lambda.clone()),
                Arc::new(sum_of_solution_terms),
            )),
        ));

        Ok(final_solution)
    }
}

/// Represents a Volterra integral equation of the second kind:
/// `y(x) = f(x) + lambda * integral_a_x(K(x, t) * y(t) dt)`
#[derive(Debug, Clone)]
pub struct VolterraEquation {
    pub y_x: Expr,    // The unknown function y(x)
    pub f_x: Expr,    // The known function f(x)
    pub lambda: Expr, // The constant lambda
    pub kernel: Expr, // The kernel K(x, t)
    pub lower_bound: Expr,
    pub var_x: String, // Variable for y(x) and f(x)
    pub var_t: String, // Integration variable for kernel and y(t)
}

impl VolterraEquation {
    pub fn new(
        y_x: Expr,
        f_x: Expr,
        lambda: Expr,
        kernel: Expr,
        lower_bound: Expr,
        var_x: String,
        var_t: String,
    ) -> Self {
        VolterraEquation {
            y_x,
            f_x,
            lambda,
            kernel,
            lower_bound,
            var_x,
            var_t,
        }
    }

    /// Solves the Volterra integral equation using the method of successive approximations.
    ///
    /// # Arguments
    /// * `iterations` - The number of iterations for the successive approximation.
    ///
    /// # Returns
    /// An `Expr` representing the approximate solution `y(x)`.
    pub fn solve_successive_approximations(&self, iterations: usize) -> Expr {
        let mut y_n = self.f_x.clone(); // y_0(x) = f(x)

        for _ in 0..iterations {
            let integral_term = Expr::Mul(
                Arc::new(self.kernel.clone()),
                Arc::new(substitute(
                    &y_n,
                    &self.var_x,
                    &Expr::Variable(self.var_t.clone()),
                )),
            );
            let integrated_val = integrate(
                &integral_term,
                &self.var_t,
                Some(&self.lower_bound),
                Some(&Expr::Variable(self.var_x.clone())),
            );
            let next_y_n = simplify(Expr::Add(
                Arc::new(self.f_x.clone()),
                Arc::new(Expr::Mul(
                    Arc::new(self.lambda.clone()),
                    Arc::new(integrated_val),
                )),
            ));
            y_n = next_y_n;
        }
        y_n
    }
}
*/
