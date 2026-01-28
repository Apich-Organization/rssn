//! # Calculus of Variations Module
//!
//! The `calculus_of_variations` module provides tools for solving problems in the calculus of variations.
//! This field of mathematical analysis deals with finding functions that maximize or minimize functionals
//! (integrals of a function and its derivatives).
//!
//! The core functionality of this module is the derivation and solving of the Euler-Lagrange equation,
//! which is a fundamental equation in this field.

use std::sync::Arc;

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::ode::solve_ode;
use crate::symbolic::simplify_dag::simplify;

/// # Euler-Lagrange Equation
///
/// Computes the Euler-Lagrange equation for a given Lagrangian functional.
///
/// The Euler-Lagrange equation is the fundamental equation of the calculus of variations.
/// For a functional $S = \int L(t, q, \dot{q}) dt$, the condition for $S$ to be stationary is:
///
/// $$\frac{d}{dt} \left( \frac{\partial L}{\partial \dot{q}} \right) - \frac{\partial L}{\partial q} = 0$$
///
/// where:
/// - $L$ is the Lagrangian (the integrand).
/// - $t$ is the independent variable.
/// - $q(t)$ is the dependent variable (generalized coordinate).
/// - $\dot{q} = dq/dt$ is the generalized velocity.
///
/// ## Arguments
/// * `lagrangian` - An [`Expr`] representing the Lagrangian $L$.
/// * `func` - The name of the dependent function $q$ as a string.
/// * `var` - The name of the independent variable $t$ as a string.
///
/// ## Returns
/// An [`Expr`] representing the left-hand side of the Euler-Lagrange equation.
///
/// ## Example: Free Particle
/// The Lagrangian for a free particle of mass $m$ is $L = \frac{1}{2} m \dot{x}^2$.
/// ```rust
/// 
/// use rssn::symbolic::calculus_of_variations::euler_lagrange;
/// use rssn::symbolic::core::Expr;
/// use std::sync::Arc;
///
/// let m = Expr::new_variable("m");
///
/// let x = Expr::new_variable("x");
///
/// let t = Expr::new_variable("t");
///
/// let x_prime = Expr::new_derivative(x.clone(), "t".to_string());
///
/// // L = 1/2 * m * (x')^2
/// let lagrangian = Expr::new_mul(
///     Expr::new_mul(
///         Expr::new_constant(0.5),
///         m,
///     ),
///     Expr::new_pow(
///         &x_prime,
///         Expr::new_constant(2.0),
///     ),
/// );
///
/// let eq = euler_lagrange(
///     &lagrangian,
///     "x",
///     "t",
/// );
/// // Result should be simplified to: m * d^2x/dt^2
/// ```
#[must_use]

pub fn euler_lagrange(
    lagrangian: &Expr,
    func: &str,
    var: &str,
) -> Expr {

    let q = Expr::Variable(
        func.to_string(),
    );

    let q_prime_str =
        format!("{func}__prime");

    let q_prime_var = Expr::Variable(
        q_prime_str.clone(),
    );

    // We need to substitute q' (which appears as Derivative(q, var) in the expression)
    // with a temporary variable q_prime_var to perform partial differentiation.
    let q_prime_expr = Expr::Derivative(
        Arc::new(q),
        var.to_string(),
    );

    let lagrangian_sub = crate::symbolic::calculus::substitute_expr(
        lagrangian,
        &q_prime_expr,
        &q_prime_var,
    );

    let dl_dq = differentiate(
        &lagrangian_sub,
        func,
    );

    let dl_dq_prime = differentiate(
        &lagrangian_sub,
        &q_prime_str,
    );

    // Substitute q' back into the partial derivative result
    let dl_dq_prime_full = crate::symbolic::calculus::substitute_expr(
        &dl_dq_prime,
        &q_prime_var,
        &q_prime_expr,
    );

    // Now take the total time derivative: d/dt (dl/dq')
    let d_dt_dl_dq_prime =
        differentiate(
            &dl_dq_prime_full,
            var,
        );

    simplify(&Expr::new_sub(
        d_dt_dl_dq_prime,
        dl_dq,
    ))
}

/// # Solve Euler-Lagrange Equation
///
/// Automatically generates and attempts to solve the Euler-Lagrange equation for a system.
///
/// This is a convenience function that computes the Euler-Lagrange equation as an ODE
/// and immediately passes it to the `solve_ode` engine.
///
/// ## Arguments
/// * `lagrangian` - The Lagrangian functional.
/// * `func` - The name of the function to solve for.
/// * `var` - The independent variable.
///
/// ## Returns
/// An [`Expr`] representing the general or particular solution to the system's motion.
#[must_use]

pub fn solve_euler_lagrange(
    lagrangian: &Expr,
    func: &str,
    var: &str,
) -> Expr {

    let el_equation = euler_lagrange(
        lagrangian,
        func,
        var,
    );

    let ode_to_solve = Expr::Eq(
        Arc::new(el_equation),
        Arc::new(Expr::new_constant(0.0)),
    );

    solve_ode(
        &ode_to_solve,
        func,
        var,
        None,
    )
}

/// # Hamilton's Principle (Least Action)
///
/// Derives the equations of motion for a physical system using the principle of stationary action.
///
/// Hamilton's Principle states that for a conservative system, the actual path $q(t)$ taken
/// by the system makes the action $S = \int L dt$ stationary.
///
/// This function is an alias for [`euler_lagrange`], providing the terminology used in physics.
///
/// ## Arguments
/// * `lagrangian` - The Lagrangian $L = T - V$ (Kinetic - Potential energy).
/// * `func` - The generalized coordinate $q$.
/// * `var` - The time variable $t$.
#[must_use]

pub fn hamiltons_principle(
    lagrangian: &Expr,
    func: &str,
    var: &str,
) -> Expr {

    euler_lagrange(
        lagrangian,
        func,
        var,
    )
}
