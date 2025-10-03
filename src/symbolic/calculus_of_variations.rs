//! # Calculus of Variations Module
//!
//! The `calculus_of_variations` module provides tools for solving problems in the calculus of variations.
//! This field of mathematical analysis deals with finding functions that maximize or minimize functionals
//! (integrals of a function and its derivatives).
//!
//! The core functionality of this module is the derivation and solving of the Euler-Lagrange equation,
//! which is a fundamental equation in this field.

use crate::symbolic::calculus::{differentiate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::ode::solve_ode;
use crate::symbolic::simplify::simplify;

/// # Euler-Lagrange Equation
///
/// Computes the Euler-Lagrange equation for a given Lagrangian.
///
/// The Euler-Lagrange equation is a second-order ordinary differential equation that describes
/// the path a system will take to extremize the action functional. The equation is given by:
/// `d/dt (∂L/∂q') - ∂L/∂q = 0`
/// where `L(t, q, q')` is the Lagrangian, `q` is the generalized coordinate, and `q'` is the generalized velocity.
///
/// ## Arguments
/// * `lagrangian` - An `Expr` representing the Lagrangian `L`. It should be an expression
///   depending on the independent variable, the function, and its first derivative.
/// * `func` - A string slice representing the name of the function `q` (e.g., "y" in y(x)).
/// * `var` - A string slice representing the name of the independent variable `t` (e.g., "x" in y(x)).
///
/// ## Returns
/// An `Expr` representing the left-hand side of the Euler-Lagrange equation. Setting this expression
/// to zero gives the equation of motion for the system.
///
/// ## Implementation Notes
/// To compute `∂L/∂q'`, where `q'` is the derivative of `q` with respect to `var`, we introduce a
/// temporary placeholder variable for `q'`. This allows us to treat `L` as a standard multivariate
/// function and perform partial differentiation. After the partial derivative is computed, the
/// placeholder is substituted back with the original derivative expression before the total
/// derivative with respect to `var` is taken.
pub fn euler_lagrange(lagrangian: &Expr, func: &str, var: &str) -> Expr {
    // `q` represents the generalized coordinate function, e.g., y(x). Here, `func` is "y" and `var` is "x".
    let q = Expr::Variable(func.to_string());

    // `q_prime` represents the generalized velocity, the derivative of the function with respect to the variable, e.g., y'(x).
    // We use a unique temporary name to treat q' as an independent variable for differentiation.
    let q_prime_str = format!("{}__prime", func);
    let q_prime_var = Expr::Variable(q_prime_str.clone());

    // Substitute the derivative expression for q' in the Lagrangian with the temporary variable.
    // This allows us to take the partial derivative of L with respect to q'.
    // Note: This relies on the string representation of the derivative, which can be brittle.
    // A more robust system might use a more structured way to identify derivatives.
    let lagrangian_sub = substitute(
        lagrangian,
        &differentiate(&q, var).to_string(),
        &q_prime_var,
    );

    // Calculate ∂L/∂q (the partial derivative of L with respect to q)
    let dl_dq = differentiate(&lagrangian_sub, func);

    // Calculate ∂L/∂q' (the partial derivative of L with respect to q')
    let dl_dq_prime = differentiate(&lagrangian_sub, &q_prime_str);

    // Substitute the temporary variable for q' back to its original derivative expression.
    let q_prime_expr = differentiate(&q, var);
    let dl_dq_prime_full = substitute(&dl_dq_prime, &q_prime_str, &q_prime_expr);

    // Calculate d/dt (∂L/∂q') (the total derivative with respect to the independent variable `var`)
    let d_dt_dl_dq_prime = differentiate(&dl_dq_prime_full, var);

    // The Euler-Lagrange equation is d/dt(dL/dq') - dL/dq = 0.
    // We return the simplified left-hand side.
    simplify(Expr::Sub(Box::new(d_dt_dl_dq_prime), Box::new(dl_dq)))
}

/// # Solve Euler-Lagrange Equation
///
/// Generates and attempts to solve the Euler-Lagrange equation for a given Lagrangian.
///
/// ## Workflow
/// 1. Calls the `euler_lagrange` function to generate the ordinary differential equation (ODE)
///    that describes the system's behavior.
/// 2. Passes the resulting equation (of the form `F(t, q, q', q'') = 0`) to an ODE solver.
///
/// ## Arguments
/// * `lagrangian` - The Lagrangian `L` of the system.
/// * `func` - The name of the generalized coordinate `q`.
/// * `var` - The name of the independent variable `t`.
///
/// ## Returns
/// An `Expr` representing the solution to the ODE. If the solver cannot find a solution,
/// it may return an expression representing the unsolved equation.
pub fn solve_euler_lagrange(lagrangian: &Expr, func: &str, var: &str) -> Expr {
    // 1. Generate the Euler-Lagrange equation, which is an ODE.
    let el_equation = euler_lagrange(lagrangian, func, var);

    // The equation is of the form F(t, q, q', q'') = 0.
    let ode_to_solve = Expr::Eq(Box::new(el_equation), Box::new(Expr::Constant(0.0)));

    // 2. Pass the resulting ODE to the ODE solver.
    // The `solve_ode` function should be capable of handling this form of second-order ODE.
    solve_ode(&ode_to_solve, func, var, None)
}

/// # Hamilton's Principle
///
/// Applies Hamilton's Principle to derive the equations of motion for a system.
///
/// Hamilton's Principle, also known as the principle of least action, states that the true
/// evolutionary path of a physical system is the one that makes the action functional
/// (the integral of the Lagrangian over time) stationary.
///
/// Applying the calculus of variations to this principle directly yields the Euler-Lagrange equations.
/// Therefore, this function is functionally equivalent to `euler_lagrange`.
///
/// ## Arguments
/// * `lagrangian` - The Lagrangian `L(t, q, q')` of the system.
/// * `func` - The name of the generalized coordinate (e.g., "q").
/// * `var` - The name of the independent variable (e.g., "t").
///
/// ## Returns
/// An `Expr` representing the Euler-Lagrange equation(s) derived from the Lagrangian,
/// which are the equations of motion according to Hamilton's Principle.
pub fn hamiltons_principle(lagrangian: &Expr, func: &str, var: &str) -> Expr {
    // Hamilton's Principle leads directly to the Euler-Lagrange equations.
    // We can therefore simply call the existing euler_lagrange function.
    euler_lagrange(lagrangian, func, var)
}
